"""Prepare independently resumable genome files without a merged FASTA copy.

The installed k2 downloader still selects/checks the same reference genomes.
Masking uses the unchanged k2mask algorithm on whole sequences. Only scheduling,
checkpointing and the on-disk library layout change.
"""
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from reference_layout import source_path
from taxonomy import atomic_json

SCHEMA = 'metag-genome-batches-v1'
INPUTS = 'library.inputs.json'
GENOMIC = {'archaea', 'bacteria', 'viral', 'human', 'fungi', 'plant',
           'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'invertebrate'}


def enabled(args):
    return os.environ.get('METAG_BATCHED_REFERENCE') == '1' and args.library in GENOMIC


def identity(path):
    info = Path(path).stat()
    if not Path(path).is_file() or not info.st_size:
        raise RuntimeError(f'Missing/empty reference product: {path}')
    return dict(bytes=info.st_size, mtime_ns=info.st_mtime_ns)


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def frozen_summary(args, selected=None):
    """Reuse the exact first assembly selection after an interrupted batch."""
    if not enabled(args):
        return selected
    folder = Path(args.db).resolve() / 'library' / args.library
    frozen = folder / 'metag.assembly_summary.txt'
    receipt = folder / 'metag.assembly_summary.json'
    if receipt.exists():
        data = json.loads(receipt.read_text())
        if data != dict(schema=SCHEMA, sha256=checksum(frozen)):
            raise RuntimeError(f'Frozen assembly summary changed: {frozen}')
        if selected and checksum(selected) != data['sha256']:
            raise RuntimeError('Frozen batch and species-selection profile disagree')
        return str(frozen)
    if frozen.exists():
        # No completed freeze receipt: this copy is not yet a trusted input.
        return selected
    return selected


def freeze_current_summary(args):
    if not enabled(args):
        return
    source = Path('assembly_summary.txt')
    frozen = Path('metag.assembly_summary.txt')
    receipt = Path('metag.assembly_summary.json')
    digest = checksum(source)
    if receipt.exists():
        if json.loads(receipt.read_text()) != dict(schema=SCHEMA, sha256=digest):
            raise RuntimeError('Assembly selection changed during a batch resume')
        return
    temp = frozen.with_suffix('.tmp')
    shutil.copyfile(source, temp)
    temp.replace(frozen)
    atomic_json(receipt, dict(schema=SCHEMA, sha256=digest))


def prepare_one(filepath, taxid, accessions, url, masker, masker_sha, threads, scan):
    """Publish a receipt only after the masker and both outputs succeed."""
    from k2_compat import write_taxid_blocks
    source = Path(source_path(str(Path(filepath).absolute())))
    if not str(filepath).endswith('.gz'):
        raise RuntimeError('Batched genomic inputs must be compressed source archives')
    output = Path(filepath).absolute().with_suffix('')
    mapping = output.with_name(output.name + '.prelim.tsv')
    receipt = output.with_name(output.name + '.metag.json')
    parameters = dict(schema=SCHEMA, source=str(source), source_identity=identity(source),
                      taxid=taxid, accessions=accessions if not taxid else {},
                      url=url, masker_sha256=masker_sha, masking='k2mask -outfmt fasta -r x')
    if receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved.get('parameters') != parameters:
            raise RuntimeError(f'Genome input/masking parameters changed; use a new database: {source}')
        try:
            valid = (identity(output) == saved['output_identity'] and
                     identity(mapping) == saved['map_identity'])
        except (OSError, RuntimeError):
            valid = False
        if valid:
            return dict(saved, reused=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + '.metag.partial')
    map_temp = mapping.with_name(mapping.name + '.partial')
    process = None
    try:
        with temporary.open('wb') as out_file:
            process = subprocess.Popen([str(masker), '-outfmt', 'fasta', '-threads', str(threads),
                                        '-r', 'x'], stdin=subprocess.PIPE, stdout=out_file)
            with gzip.open(source, 'rt') as in_file:
                sequences, bases, headers = write_taxid_blocks(
                    in_file, process.stdin, filepath, {filepath: taxid}, accessions,
                    {filepath: url}, binary=True)
            process.stdin.close()
            if process.wait() != 0:
                raise RuntimeError(f'k2mask failed for {source}; no completion receipt written')
        if not sequences or not bases or not temporary.stat().st_size:
            raise RuntimeError(f'Genome has no usable sequence: {source}')
        if identity(source) != parameters['source_identity']:
            raise RuntimeError(f'Source genome changed during processing: {source}')
        # k2's scanner iterates these already collected headers, not the FASTA.
        with temporary.open() as fasta, map_temp.open('w') as out_map:
            scan(fasta, out_map, sequence_to_url=headers)
        identity(map_temp)
        temporary.replace(output)
        map_temp.replace(mapping)
        saved = dict(parameters=parameters, output=str(output), map=str(mapping),
                     output_identity=identity(output), map_identity=identity(mapping),
                     sequences=sequences, bases=bases)
        atomic_json(receipt, saved)
        return dict(saved, reused=False)
    finally:
        if process is not None:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            if process.stdin is not None:
                try:
                    process.stdin.close()
                except (BrokenPipeError, OSError):
                    pass
        temporary.unlink(missing_ok=True)
        map_temp.unlink(missing_ok=True)


def prepare_library(args, manifest, accessions, urls, scan, log):
    if args.protein or args.no_masking:
        raise RuntimeError('Batched references require the existing nucleotide/masked settings')
    jobs = int(os.environ.get('METAG_MASK_JOBS', '2'))
    if jobs < 1 or args.masker_threads < 1:
        raise ValueError('Masking concurrency must be positive')
    masker = Path(shutil.which('k2')).resolve().with_name('k2mask')
    masker_sha = checksum(masker)
    files = sorted(manifest)
    if not files:
        raise RuntimeError('Empty genome manifest')
    # Only one bounded window of futures exists; no library-wide header dictionary.
    log.info(f'Batched reference: {len(files)} genomes; mask jobs={jobs}; '
             f'threads/job={args.masker_threads}; merged FASTA disabled\n')
    results, pending = {}, {}
    todo = iter(files)
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        def submit_one():
            path = next(todo, None)
            if path is None:
                return False
            future = pool.submit(prepare_one, path, manifest[path], accessions,
                                 urls[path], masker, masker_sha, args.masker_threads, scan)
            pending[future] = path
            return True
        for _ in range(jobs):
            submit_one()
        while pending:
            finished, _ = wait(pending, timeout=60, return_when=FIRST_COMPLETED)
            if not finished:
                log.info(f'BATCH {len(results)}/{len(files)}: {len(pending)} genomes active; '
                         f'elapsed={(time.monotonic()-started)/60:.1f} min\n')
                continue
            for future in finished:
                path = pending.pop(future)
                result = future.result()
                results[path] = result
                log.info(f'BATCH {len(results)}/{len(files)}: taxid={manifest[path]}; '
                         f'{"reused" if result["reused"] else "completed"} {Path(path).name}; '
                         f'elapsed={(time.monotonic()-started)/60:.1f} min\n')
            for _ in finished:
                submit_one()
    # Mapping is small; concatenate it in the same genome order as upstream.
    map_temp = Path('prelim_map.txt.partial')
    with map_temp.open('wb') as target:
        target.write(f'# prelim_map for {args.library}\n'.encode())
        for path in files:
            with Path(results[path]['map']).open('rb') as source:
                shutil.copyfileobj(source, target, 1 << 20)
    map_temp.replace('prelim_map.txt')
    root = Path.cwd()
    rows = [dict(path=str(Path(results[p]['output']).relative_to(root)),
                 identity=results[p]['output_identity']) for p in files]
    atomic_json(root / INPUTS, dict(schema=SCHEMA, library=args.library,
                                  files=rows, map_identity=identity('prelim_map.txt')))
    log.info(f'Completed {args.library}: {len(files)} genome files ready for direct indexing; '
             'no library.fna copy created\n')
    return None  # The per-genome maps were already published above.


def library_inputs(folder):
    """Reject missing/modified batches before any bytes reach Kraken2."""
    folder = Path(folder).resolve()
    document = json.loads((folder / INPUTS).read_text())
    if document.get('schema') != SCHEMA or document.get('library') != folder.name:
        raise RuntimeError(f'Invalid batched library manifest: {folder}')
    if identity(folder / 'prelim_map.txt') != document['map_identity']:
        raise RuntimeError(f'Batched taxonomy map changed: {folder}')
    paths = []
    for row in document['files']:
        relative = Path(row['path'])
        path = folder / relative
        if relative.is_absolute() or '..' in relative.parts or not path.resolve().is_relative_to(folder):
            raise RuntimeError(f'Genome path escapes its library: {path}')
        if identity(path) != row['identity']:
            raise RuntimeError(f'Genome output changed/incomplete: {path}')
        paths.append(path)
    if not paths or len(set(paths)) != len(paths):
        raise RuntimeError(f'Empty/duplicate genome inputs: {folder}')
    return paths
