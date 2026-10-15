#!/usr/bin/env python3
"""Verify completed direct-UKB runs before exporting disposable VCFs.

Exit codes for check: 0 = complete, 1 = pending/unsupported, 2 = incompatible.
Original input stat records remain strict. Disposable sample paths are replaced
by the independently reconstructed cohort and sample-order hashes, never by
file presence alone. Published archives supply immutable historical evidence.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import sys
import tarfile
import tempfile
import zlib
from functools import lru_cache

sys.dont_write_bytecode = True
import ibdmix
import ukb_prepare as ukb

F = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


@lru_cache(maxsize=128)
def file_digest(path, size, mtime_ns):
    with open(path, 'rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def sha(path):
    path = Path(path)
    info = path.stat()
    return file_digest(str(path), info.st_size, info.st_mtime_ns)


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.resume-', dir=path.parent)
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(value, handle)
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def wanted(name):
    return (name in {'.complete', 'run.meta.tsv', 'cache.meta.tsv', 'inputs/ukb.samples.tsv'}
            or (name.startswith('inputs/ukb-scoring/') and name.endswith('/target.json'))
            or (name.startswith('qc/masks/') and name.endswith('.preparation.json')))


def evidence(run):
    """Freeze small receipts before stage-qc can overwrite historical inputs."""
    record = (run / 'run.meta.tsv').read_text()
    key = digest(record.encode())
    receipt = run / '.ibdmix.resume.json'
    archive = run / 'ibdmix.raw.tar.gz'
    archive_sig = ukb.signature([archive]) if archive.is_file() else []
    completion = (run / '.complete').exists()
    if receipt.is_file():
        saved = json.loads(receipt.read_text())
        if (saved.get('schema') == 1 and saved.get('record_sha256') == key
                and saved.get('archive') == archive_sig
                and (archive_sig or saved.get('completion') == completion)):
            return saved
    texts, outputs = {}, {}
    if archive.is_file():
        observed = {}
        manifest = None
        with tarfile.open(archive, 'r|gz') as bundle:
            for member in bundle:
                name = member.name
                if not (wanted(name) or name.startswith('final/') or name == 'GU-MANIFEST.json'):
                    continue
                require(member.isfile(), f'nonregular archived receipt: {name}')
                with bundle.extractfile(member) as handle:
                    data = handle.read()
                observed[name] = digest(data)
                if wanted(name):
                    texts[name] = data.decode()
                elif name == 'GU-MANIFEST.json':
                    manifest = json.loads(data)
                elif name.endswith('.segments.tsv.gz'):
                    with gzip.GzipFile(fileobj=io.BytesIO(data)) as handle:
                        while handle.read(1024 * 1024):
                            pass
                    outputs[name] = observed[name]
        require(manifest is not None, 'archive lacks GU-MANIFEST.json')
        declared = {row['path']: row for row in manifest['files']}
        for name, checksum in observed.items():
            if name != 'GU-MANIFEST.json':
                require(declared.get(name, {}).get('sha256') == checksum, f'archive checksum mismatch: {name}')
        require(texts.get('run.meta.tsv') == record, 'published and working run contracts differ')
    else:
        for path in run.rglob('*'):
            name = path.relative_to(run).as_posix()
            if path.is_file() and wanted(name):
                texts[name] = path.read_text()
        if completion:
            for path in (run / 'final').glob('*.segments.tsv.gz'):
                with gzip.open(path, 'rb') as handle:
                    while handle.read(1024 * 1024):
                        pass
                outputs[path.relative_to(run).as_posix()] = sha(path)
    saved = dict(schema=1, record_sha256=key, archive=archive_sig, completion=completion, texts=texts, outputs=outputs)
    atomic_json(receipt, saved)
    return saved


def command(path):
    env, argv = dict(os.environ), []
    for line in path.read_text().splitlines():
        words = shlex.split(line)
        if words[:1] == ['export']:
            env.update(word.split('=', 1) for word in words[1:])
        elif words[:1] == ['exec']:
            argv = words[1:]
    require(len(argv) >= 2 and argv[1] in {'ibdmix', 'ibdmix-GPU'}, 'not an IBDmix worker command')
    env['IBDMIX_BACKEND'] = 'cuda' if argv[1] == 'ibdmix-GPU' else 'cpu'
    args = argv[2:]
    if args[:1] == ['run']:
        args = args[1:]
    if '--keep-males' in args:
        args.insert(args.index('--keep-males') + 1, 'TRUE')
    require(len(args) % 2 == 0 and all(x.startswith('--') for x in args[::2]), 'unsupported worker arguments')
    return env, dict(zip(args[::2], args[1::2]))


def records(text):
    return [line.replace('\\t', '\t').split('\t') for line in text.splitlines() if line]


def scoring(saved, chrom):
    return json.loads(saved['texts'][f'inputs/ukb-scoring/chr{chrom}/target.json'])


def disposable(path):
    p = Path(path)
    return ((p.parent.parent == Path('/tmp/gu-psam-chunks') and p.name == 'chunk.keep')
            or (p.parent.parent == Path('/tmp/gu-cohorts') and p.name in {'samples.txt', 'cohort.keep'}))


def stable_scoring(value):
    """Only generated selection files have a content-based replacement below."""
    value = json.loads(json.dumps(value))
    value['preparation']['inputs'] = [x for x in value['preparation']['inputs'] if not disposable(x['path'])]
    return value


def verify_stats(items):
    for item in items:
        require(ukb.signature([item['path']])[0] == item, f'input changed: {item["path"]}')


@lru_cache(maxsize=32)
def psam_rows(path, size, mtime_ns):
    return ukb.read_psam(path)


def samples_at(path):
    path = Path(path).resolve()
    info = path.stat()
    # attach_panel mutates rows; cached PSAM data must stay pristine.
    return [dict(row) for row in psam_rows(str(path), info.st_size, info.st_mtime_ns)]


def verify_cohort(prep, env, opts):
    chrom = opts['--chr']
    root = Path(opts['--target-dir']).parent.resolve()
    require(prep['pgen_root'] == str(root), 'target PGEN directory changed')
    require(Path(opts['--target-dir']).name == 'chr', 'unsupported UKB PGEN prefix')
    options = prep['options']
    male = opts.get('--keep-males') == 'TRUE'
    require(options.get('existing_pgen') is True and options['male_only'] == male, 'PGEN mode or sex selection changed')
    require(options['grch'] == int(opts['--grch']) and options['ukb_source'] == env['GU_UKB_SOURCE'], 'build or UKB source changed')
    require(options['ukb_require_phase'] == ukb.boolean(env.get('UKB_REQUIRE_PHASE', 'false')), 'phase requirement changed')
    require(not env.get('UKB_INFO_FILE') and not env.get('UKB_QC_PASS_VARIANTS'), 'external UKB variant filters require full validation')
    metadata = root / ('chr1.psam' if chrom == 'X' else f'chr{chrom}.psam')
    samples = samples_at(metadata)
    if '--keep-psam' in opts:
        name, size, index = opts['--keep-psam'].rsplit(',', 2)
        path, size, index = Path(name).resolve(), int(size), int(index)
        require(size > 0 and index > 0, 'invalid PSAM chunk')
        info = options.get('keep_psam', {})
        require((str(path), size, index, sha(path)) == (info.get('source'), info.get('chunk_size'), info.get('chunk_index'), info.get('source_sha256')), 'PSAM chunk selection changed')
        keep = [s['sample'] for s in samples_at(path)[(index - 1) * size:index * size]]
    elif opts.get('--keep') or env.get('GU_KEEP'):
        keep = ukb.read_ids(opts.get('--keep') or env['GU_KEEP'])
        require(not options.get('keep_psam'), 'PSAM chunk selection changed')
    else:
        raise ValueError('missing fixed-cohort selection')
    panel = opts.get('--sample-panel')
    samples = ukb.selected_samples(ukb.attach_panel(samples, panel), keep, '1')
    if male:
        samples = [s for s in samples if s['sex'] == 'male']
    ukb.validate_cohort(samples, argparse.Namespace(sample_panel=panel))
    require(ukb.cohort_hash(samples) == prep['cohort_metadata_sha256'], 'sample cohort/sex/population metadata changed')
    _, paths = ukb.pfile(root, chrom)
    full = ukb.attach_panel(samples_at(paths[1]), panel)
    keep = [s['sample'] for s in samples if chrom != 'X' or s['sex'] == 'male']
    if chrom == 'X':
        present = {s['sample'] for s in full}
        keep = [s for s in keep if s in present]
    selected = ukb.selected_samples(full, keep, chrom)
    require(ukb.cohort_hash(selected) == prep['scored_sample_metadata_sha256'] and ukb.sample_hash(selected) == prep['sample_order_sha256'], 'chromosome sample set/order changed')
    recorded = {x['path'] for x in prep['inputs']}
    require(all(str(p.resolve()) in recorded for p in paths), 'PGEN input paths changed')
    fasta = Path(env.get('UKB_REF_FASTA', str(Path(env.get('GU_REF_ROOT', '/mnt/f/gen')) / f'1hgp/{opts["--grch"]}/GRCH{opts["--grch"]}.fasta.gz'))).resolve()
    require(str(fasta) in recorded and str(fasta) + '.fai' in recorded, 'preparation reference FASTA changed')
    verify_stats(x for x in prep['inputs'] if not disposable(x['path']))


def verify_request(saved, env, opts, complete):
    chrom = opts['--chr']
    old = records(saved['texts']['run.meta.tsv'])
    ibdmix.contract(saved['texts']['run.meta.tsv'])
    values = {row[0]: row[1] for row in old if len(row) == 2}
    cuda = env['IBDMIX_BACKEND'] == 'cuda'
    version = re.search(r'^pipeline_version=([^\n]+)', (F / 'ibdmix.sh').read_text(), re.M)[1]
    if cuda:
        version += '.cuda-segment-v1'
    expected = dict(pipeline_version=version, profile=env.get('IBDMIX_PROFILE', 'multi_reference'),
                    refs=env['IBDMIX_REFS'], genome_build='b' + opts['--grch'],
                    background_filter=env.get('IBDMIX_AFR_DENISOVAN_FILTER', '0'),
                    custom_masks=env.get('IBDMIX_MASK_DIR', ''),
                    lod_cut=env.get('IBDMIX_LOD', '4'), len_cut=env.get('IBDMIX_MIN_BP', '50000'),
                    emit_lod_cut=env.get('emit_lod_cut', env.get('IBDMIX_LOD', '4')),
                    minor_allele_count=env.get('minor_allele_count', '1'),
                    archaic_error=env.get('archaic_error', '0.01'),
                    modern_error_max=env.get('modern_error_max', '0.002'),
                    modern_error_proportion=env.get('modern_error_proportion', '2'))
    if chrom == 'X':
        expected.update(x_scoring_profile=env.get('IBDMIX_X_PROFILE', 'nonpar-v2'),
                        x_minor_allele_count=env.get('IBDMIX_X_MINOR_ALLELE_COUNT', '1'))
    for key, value in expected.items():
        require(values.get(key) == value, f'parameter changed: {key} ({values.get(key)!r} -> {value!r})')
    soft = Path(env.get('GU_SOFT', '/mnt/d/software/gu'))
    prepare = Path(env.get('IBDMIX_PREPARE_BIN', soft / 'ibdmix-prepare/prepare'))
    require(sha(prepare) == values['genotype_binary_sha256'], 'genotype preparation binary changed')
    if cuda:
        require(sha(env.get('IBDMIX_GPU_BIN', soft / 'ibdmix-gpu/ibdmix-gpu')) == values['gpu_binary_sha256'], 'GPU caller binary changed')
    else:
        runtime = Path(env.get('IBDMIX_RUNTIME', soft / 'ibdmix'))
        if not runtime.is_dir():
            runtime = soft / 'IBDmix'
        caller = runtime / 'build/src/ibdmix'
        info = caller.stat()
        require(['software', f'{caller}:{info.st_size}:{int(info.st_mtime)}'] in old, 'CPU caller binary changed')
    value = scoring(saved, chrom)
    contract = value['contract']
    verify_cohort(contract['preparation'], env, opts)
    require(contract['references'] == env['IBDMIX_REFS'].split(), 'archaic reference selection changed')
    source_root = Path(opts['--archaic-path']).resolve()
    sources = [ukb._ukb_scoring_source(source_root, ref, chrom).resolve() for ref in contract['references']]
    require([x['path'] for x in contract['inputs'] if Path(x['path']) != F / 'ukb_prepare.py'] == list(map(str, sources)), 'archaic source paths changed')
    verify_stats(contract['inputs'])
    mask_root = str(Path(env.get('IBDMIX_MASK_ROOT', source_root.parent / 'mask')).resolve())
    require(contract['original_mask_root'] == mask_root, 'mask root changed')
    fallback = source_root.parent / 'mask'
    require(contract.get('original_fallback_mask_root') in {None, str(fallback)}, 'fallback mask root changed')
    fasta = env.get('IBDMIX_REFERENCE_FASTA', str(Path(env.get('GU_REF_ROOT', '/mnt/f/gen')) / '1hgp/37/GRCH37.fasta.gz'))
    require(values['reference_fasta'].rsplit(':', 2)[0] == fasta, 'mask reference FASTA changed')
    # The scoring receipt must belong to these calls, including after a failed
    # rerun has overwritten the working inputs directory.
    prep_hash = digest(json.dumps(contract['preparation'], sort_keys=True).encode())
    require(['modern_source', 'chr' + chrom, 'preparation_contract_sha256', prep_hash] in old,
            'scoring receipt does not match the original call contract')
    if complete:
        masks = [json.loads(text) for name, text in saved['texts'].items() if name.startswith('qc/masks/')]
        require(masks, 'completed run lacks mask provenance')
        for mask in masks:
            sig = mask['signature']
            require(sig['code'] == sha(F / 'ibdmix.py'), 'mask implementation changed')
            require(sig['refs'] == contract['references'], 'mask references changed')
            verify_stats(dict(path=p, bytes=n, mtime_ns=t) for p, n, t in sig['inputs'])
            if chrom == 'X':
                enabled = all(ibdmix.local_mask_ready(ibdmix.axt_reference(chrom, species, None)[1]) for species in ['panTro2', 'ponAbe2', 'rheMac2'])
                require(sig['cpg_enabled'] == enabled, 'chrX CpG resources changed')


def check(cmd):
    env, opts = command(cmd)
    if (env.get('GU_UKB_DIRECT') != '1' or '--chr' not in opts
            or opts.get('--replace-ibdmix', 'FALSE') != 'FALSE'):
        return 1
    run = Path(env['IBDMIX_OUT'])
    if not (run / 'run.meta.tsv').is_file():
        return 1
    saved = evidence(run)
    complete = '.complete' in saved['texts']
    verify_request(saved, env, opts, complete)
    if not complete:
        print(f'[GU CMD] RESUME unit=chr{opts["--chr"]} reason=partial_calls_verified', file=sys.stderr)
        return 1
    final = f'final/all_archaic_refs.lod{env.get("IBDMIX_LOD", "4")}.len{env.get("IBDMIX_MIN_BP", "50000")}.segments.tsv.gz'
    require(final in saved['outputs'], 'completed run lacks verified final segments')
    if (run / final).exists():
        require(sha(run / final) == saved['outputs'][final], 'final segments changed or corrupted')
    else:
        require(bool(saved['archive']), 'final segments are missing')
    print(run / final if (run / final).exists() else run / 'ibdmix.raw.tar.gz')
    return 0


def provenance(previous, current, run):
    before, after = previous.read_text(), current.read_text()
    if not ibdmix.differences(before, after):
        return 0
    receipt = run / '.ibdmix.resume.json'
    if receipt.is_file():
        saved = json.loads(receipt.read_text())
        require(saved['texts']['run.meta.tsv'] == (run / 'run.meta.tsv').read_text(), 'resume receipt no longer matches original run')
        old_paths = [n for n in saved['texts'] if n.startswith('inputs/ukb-scoring/') and n.endswith('/target.json')]
        require(len(old_paths) == 1, 'expected one chromosome scoring receipt')
        old = json.loads(saved['texts'][old_paths[0]])['contract']
        new = json.loads((run / old_paths[0]).read_text())['contract']
        require(stable_scoring(old) == stable_scoring(new), 'original sources, cohort, preparation options or scored sites changed')
        for text, value in [(before, old), (after, new)]:
            prep_hash = digest(json.dumps(value['preparation'], sort_keys=True).encode())
            require(['modern_source', 'chr' + value['chromosome'], 'preparation_contract_sha256', prep_hash] in records(text),
                    'scoring receipt does not match the compared call contract')
        # Only exchange-file identities are disposable. Masks, caller binaries,
        # thresholds and all other provenance still use the existing comparator.
        def normalize(text):
            return '\n'.join('\t'.join(row) for row in records(text)
                             if row[0] not in {'sample_file', 'archaic_vcf'}
                             and not (row[0] == 'modern_source' and row[2] == 'preparation_contract_sha256')) + '\n'
        diff = ibdmix.differences(normalize(before), normalize(after))
    else:
        diff = ibdmix.differences(before, after)
    if diff:
        print(diff, end='', file=sys.stderr)
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('check')
    p.add_argument('command', type=Path)
    p = sub.add_parser('provenance')
    p.add_argument('previous', type=Path)
    p.add_argument('current', type=Path)
    p.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    try:
        return check(args.command) if args.action == 'check' else provenance(args.previous, args.current, args.run)
    except (ValueError, OSError, KeyError, IndexError, TypeError, EOFError, tarfile.TarError, zlib.error) as error:
        print(f'ERROR: IBDmix resume check: {error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
