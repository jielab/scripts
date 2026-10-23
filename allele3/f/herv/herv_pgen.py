#!/usr/bin/env python3
"""Prepare the existing HERV sequence models directly from phased PLINK 2 files.

Only a compact, temporary single-sample variant stream bridges to the unchanged
VCF preparation engine. No individual chromosome/genome FASTA is materialized.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import fcntl
import gzip
import importlib.metadata
import io
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time
from types import SimpleNamespace

import herv_assembly as engine


def source_files(root, chrom):
    base = root/chrom
    pvar = Path(str(base)+'.pvar.zst')
    if not pvar.is_file():
        pvar = Path(str(base)+'.pvar')
    files = {'pgen': Path(str(base)+'.pgen'), 'pvar': pvar,
             'psam': Path(str(base)+'.psam')}
    for path in files.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    pgi = Path(str(base)+'.pgen.pgi')
    if pgi.exists():
        files['pgi'] = pgi
    return files


def file_state(path):
    stat = path.stat()
    return dict(path=str(path.resolve()), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def sample_index(psam, sample):
    found, total = [], 0
    with psam.open() as stream:
        header = stream.readline().lstrip('#').split()
        if 'IID' not in header:
            raise ValueError(f'PSAM requires an IID header: {psam}')
        column = header.index('IID')
        for line in stream:
            if not line.strip():
                continue
            row = line.split()
            if len(row) != len(header):
                raise ValueError(f'Malformed PSAM: {psam}')
            if row[column] == sample:
                found.append(total)
            total += 1
    if len(found) != 1:
        raise ValueError(f'Expected exactly one IID={sample} in {psam}; found {len(found)}')
    return found[0], total


@contextmanager
def pvar_text(path):
    if path.suffix == '.zst':
        import zstandard
        with path.open('rb') as raw, zstandard.ZstdDecompressor().stream_reader(raw) as reader:
            with io.TextIOWrapper(reader) as stream:
                yield stream
    else:
        with path.open() as stream:
            yield stream


def compact_info(info):
    # These are all INFO fields consulted by the unchanged reconstruction policy.
    return ';'.join(x for x in info.split(';')
                    if x.partition('=')[0] in {'END', 'IMPRECISE', 'CIPOS', 'CIEND'}) or '.'


def export_selected(files, sample, chrom, length, output):
    import numpy as np
    import pgenlib
    index, total = sample_index(files['psam'], sample)
    counts = Counter()
    # The minimal native PVAR reader supplies multiallelic offsets, while the
    # streaming text pass preserves FILTER and structural-variant metadata.
    with pgenlib.PvarReader(os.fsencode(files['pvar']), omit_chrom=True, omit_pos=True) as variants, \
         pgenlib.PgenReader(os.fsencode(files['pgen']), raw_sample_ct=total,
             sample_subset=np.asarray([index], dtype=np.uint32), pvar=variants) as genotypes, \
         pvar_text(files['pvar']) as source, gzip.open(output, 'wt', compresslevel=1) as dest:
        if not genotypes.hardcall_phase_present():
            raise ValueError(f'PGEN contains no hardcall phase information: {files["pgen"]}; use the phased VCF')
        n = variants.get_variant_ct()
        if n != genotypes.get_variant_ct():
            raise ValueError('PVAR/PGEN variant count mismatch')
        columns = None
        for line in source:
            if line.startswith('#CHROM'):
                columns = line.lstrip('#').rstrip('\n').split('\t')
                break
            if not line.startswith('##'):
                raise ValueError('PVAR requires a #CHROM header')
        required = {'CHROM', 'POS', 'ID', 'REF', 'ALT'}
        if columns is None or not required <= set(columns):
            raise ValueError('Missing required PVAR columns')
        lookup = {name: columns.index(name) for name in columns}
        dest.write('##fileformat=VCFv4.3\n'+f'##contig=<ID={chrom},length={length}>\n')
        dest.write('##FILTER=<ID=PGEN_FILTERED,Description="Non-PASS FILTER retained from PVAR">\n')
        for name, number, kind in [('END', '1', 'Integer'), ('IMPRECISE', '0', 'Flag'),
                                   ('CIPOS', '2', 'Integer'), ('CIEND', '2', 'Integer')]:
            dest.write(f'##INFO=<ID={name},Number={number},Type={kind},Description="Preserved from PVAR">\n')
        dest.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased PGEN hardcalls">\n')
        dest.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')
        batch = 8192
        alleles = np.empty((batch, 2), dtype=np.int32)
        phase = np.empty((batch, 1), dtype=np.uint8)
        previous_pos = 0
        for ordinal, line in enumerate(source):
            if ordinal >= n or not line.strip() or line.startswith('#'):
                raise ValueError('PVAR row count/order differs from the native reader')
            if ordinal % batch == 0:
                genotypes.read_alleles_and_phasepresent_range(ordinal, min(ordinal+batch, n), alleles, phase)
            row = line.rstrip('\n').split('\t')
            if len(row) != len(columns):
                raise ValueError(f'Malformed PVAR row {ordinal+1}')
            position = int(row[lookup['POS']])
            row_chrom = row[lookup['CHROM']]
            if ('chr'+row_chrom.removeprefix('chr')) != chrom or not previous_pos <= position <= length:
                raise ValueError(f'PVAR chromosome/position/order mismatch at row {ordinal+1}')
            previous_pos = position
            counts['total_variants'] += 1
            a, b = map(int, alleles[ordinal % batch])
            if a == b == 0:
                continue
            phased = bool(phase[ordinal % batch, 0])
            if a < 0 or b < 0:
                counts['missing_GT'] += 1
            elif a != b:
                counts['phased_heterozygotes' if phased else 'unphased_heterozygotes'] += 1
            else:
                counts['homozygous_alternate'] += 1
            # Homozygous calls have no haplotype ambiguity even when PGEN does
            # not retain their original VCF separator. Missing calls stay missing.
            separator = '|' if phased or (a == b and a >= 0) else '/'
            gt = separator.join('.' if x < 0 else str(x) for x in (a, b))
            info = compact_info(row[lookup['INFO']]) if 'INFO' in lookup else '.'
            filt = row[lookup['FILTER']] if 'FILTER' in lookup else '.'
            # FILTER names can be arbitrary; the existing engine only needs
            # PASS/dot vs failure, so declare one deterministic failure code.
            if filt not in ('.', 'PASS'):
                filt = 'PGEN_FILTERED'
            dest.write('\t'.join([chrom, str(position), row[lookup['ID']], row[lookup['REF']],
                row[lookup['ALT']], '.', filt, info, 'GT', gt])+'\n')
            counts['selected_records'] += 1
        if counts['total_variants'] != n:
            raise ValueError('PVAR ended before all PGEN variants were read')
    return dict(counts)


def prepare_chrom(args, chrom, length):
    files = source_files(args.pfile_dir, chrom)
    sources = {name: file_state(path) for name, path in files.items()}
    inputs = dict(method=engine.METHOD, input_source='pgen', sample=args.sample,
        chromosome=chrom, length=length, fasta_sha256=args.fasta_sha256,
        uncertain_flank=args.uncertain_flank, sources=sources,
        preparer_sha256=engine.digest(__file__), pgenlib_version=importlib.metadata.version('pgenlib'))
    root = args.out/'variants'/chrom
    receipt = root/'done.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old['inputs'] != inputs or any(engine.digest(root/n) != h for n,h in old['files'].items()):
            raise ValueError(f'Changed PGEN inputs/outputs: {root}')
        return old
    t0 = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='herv-pgen-', dir=args.scratch) as folder:
        temp = Path(folder)
        # PGEN performs many small reads; staging also avoids WSL/NTFS latency.
        # At most --jobs chromosome filesets are staged and each is removed
        # after its compact models are published.
        staged = {}
        for name, source in files.items():
            staged[name] = temp/source.name
            shutil.copyfile(source, staged[name])
        staging_seconds = time.monotonic()-t0
        selected = temp/(chrom+'.vcf.gz')
        export_counts = export_selected(staged, args.sample, chrom, length, selected)
        export_seconds = time.monotonic()-t0
        if sources != {name: file_state(path) for name, path in files.items()}:
            raise ValueError('PGEN input files changed during reading')
        shared = SimpleNamespace(out=temp/'models', vcf_dir=temp, fasta=args.fasta,
            fasta_sha256=args.fasta_sha256, sample=args.sample, uncertain_flank=args.uncertain_flank,
            vcf_scratch=args.scratch)
        result = engine.prepare_chrom(shared, chrom, length)
        result.update(inputs=inputs, pgen_export_counts=export_counts, export_seconds=export_seconds,
                      staging_seconds=staging_seconds,
                      seconds=time.monotonic()-t0, completed=engine.now())
        result['query_command_note'] = 'Temporary compact PGEN hardcall export; removed after preparation.'
        root.mkdir(parents=True, exist_ok=True)
        source_root = shared.out/'variants'/chrom
        for name in [*result['files'], 'query.log']:
            shutil.copy2(source_root/name, root/name)
        # Completion is published only after every model/audit file is copied.
        engine.save(receipt, result)
    print(json.dumps(dict(event='prepared_pgen', chromosome=chrom,
        seconds=round(result['seconds'], 2), export_counts=export_counts)), flush=True)
    return result


def main():
    os.environ['PATH'] = str(Path(sys.executable).parent)+os.pathsep+os.environ.get('PATH', '')
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['prepare'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--fasta', type=Path, required=True)
    p.add_argument('--pfile-dir', type=Path, required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--chromosomes', default=','.join('chr'+str(i) for i in range(1,23)))
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--uncertain-flank', type=int, default=1000)
    p.add_argument('--scratch', type=Path, default=Path('/tmp/metag-herv-pgen'))
    args = p.parse_args()
    args.chromosomes = args.chromosomes.split(',')
    if args.jobs < 1 or args.uncertain_flank < 0 or len(set(args.chromosomes)) != len(args.chromosomes) or \
       any(not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])', c) for c in args.chromosomes):
        p.error('Specify positive jobs, nonnegative flank and unique chromosomes chr1..chr22')
    args.out = args.out.resolve()
    args.fasta = args.fasta.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    args.scratch.mkdir(parents=True, exist_ok=True)
    with (args.out/'prepare.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        lengths = {r.split('\t')[0]: int(r.split('\t')[1])
                   for r in Path(str(args.fasta)+'.fai').read_text().splitlines()}
        args.fasta_sha256 = engine.digest(args.fasta)
        # PVAR parsing is CPU-bound Python work; separate processes let the
        # default two chromosome preparations use both allocated CPUs.
        with ProcessPoolExecutor(args.jobs) as pool:
            futures = [pool.submit(prepare_chrom, args, c, lengths[c]) for c in args.chromosomes]
            for future in futures:
                future.result()
        engine.save(args.out/'prepared.json', dict(method=engine.METHOD, input_source='pgen',
            sample=args.sample, chromosomes=args.chromosomes, fasta=str(args.fasta),
            fasta_sha256=args.fasta_sha256, uncertain_flank=args.uncertain_flank,
            pfile_dir=str(args.pfile_dir.resolve()), preparer_sha256=engine.digest(__file__),
            policy='Same reconstruction/masking policy as herv_assembly.prepare_chrom; '
                   'PGEN hardcalls and phase, PVAR REF/ALT/FILTER/SV INFO. No dosage rounding.',
            completed=engine.now()))


if __name__ == '__main__':
    main()
