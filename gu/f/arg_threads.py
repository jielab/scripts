#!/usr/bin/env python3
"""PGEN -> Threads -> chromosome-coordinate, mutation-bearing TRACE trees.

Runs in an isolated environment. Stage manifests bind checkpoints to inputs,
options, software versions and output stamps; incomplete stages are never reused.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time

import numpy as np
import pgenlib
import tskit
import tszip

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arg import trace_validate, needle_x_panel


def env(key, default=None):
    value = os.environ.get('REFGEN_' + key, default)
    if value is None:
        raise ValueError(f'Missing REFGEN_{key}')
    return value


def stamp(path):
    p = Path(path).resolve(strict=True)
    s = p.stat()
    return [str(p), s.st_size, s.st_mtime_ns]


def atomic_json(path, data):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.next')
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
    tmp.replace(path)


def checkpoint(path, request, files):
    atomic_json(path, dict(request=request, outputs=[stamp(p) for p in files]))


def verified(path, request=None):
    try:
        data = json.loads(Path(path).read_text())
        return ((request is None or data['request'] == request)
                and all(stamp(row[0]) == row for row in data['outputs']))
    except (OSError, ValueError, KeyError, TypeError):
        return False


def run(args, log):
    import shlex
    args = list(map(str, args))
    print('[THREADS] ' + shlex.join(args), flush=True)
    with Path(log).open('a') as handle:
        handle.write('\nCOMMAND ' + shlex.join(args) + '\n'); handle.flush()
        start = time.monotonic()
        proc = subprocess.Popen(args, stdout=handle, stderr=subprocess.STDOUT)
        try:
            while True:
                try:
                    code = proc.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    print(f'PROGRESS elapsed={time.monotonic()-start:.0f}s log={log}', flush=True)
            if code:
                with Path(log).open('rb') as tail:
                    tail.seek(max(0, tail.seek(0, 2) - 5000))
                    print(tail.read().decode(errors='replace'), file=sys.stderr)
                raise subprocess.CalledProcessError(code, args)
        finally:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill(); proc.wait()


def map_for(chrom):
    pattern = env('MAP_PATTERN', '')
    if pattern:
        path = Path(pattern.replace('{chr}', chrom).replace('%CHR%', chrom))
        if path.is_file():
            return path
    root = Path(env('MAP_DIR'))
    build = env('GRCH')
    for name in (f'chr{chrom}.map', f'chr{chrom}.txt', f'chr{chrom}.map.gz',
                 f'chr{chrom}.b{build}.gmap.gz', f'genetic_map_GRCh{build}_chr{chrom}.txt',
                 f'genetic_map_chr{chrom}_combined_b{build}.txt',
                 f'genetic_map_chr{chrom}_combined_b{build}.txt.map'):
        for suffix in ('', '.gz'):
            if (root / (name + suffix)).is_file():
                return root / (name + suffix)
    raise ValueError(f'No genetic map for chr{chrom} in {root}; use --map-pattern')


def read_psam(path):
    with Path(path).open() as h:
        fields = h.readline().lstrip('#').split()
        rows = [dict(zip(fields, line.split())) for line in h if line.strip()]
    names = [r['IID'] for r in rows]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate sample IID')
    return rows


def repack(source, output, chrom, source_map):
    """Stream PGEN; enforce phase and unique positions, pack only real X copies."""
    from build_argneedle_map import read_source
    samples = read_psam(str(source) + '.psam')
    ids, columns = [], []
    for i, row in enumerate(samples):
        if chrom == 'X' and row.get('SEX') not in ('1', '2'):
            raise ValueError('Non-PAR X requires known sex for every selected sample')
        copies = 1 if chrom == 'X' and row['SEX'] == '1' else 2
        for hap in range(copies):
            ids.append((row['IID'], hap + 1)); columns.append(2*i + hap)
    if len(ids) < 4 or len(ids) % 2:
        raise ValueError('Threads needs an even number of >=4 real haplotypes; check X panel policy')
    variants = []
    with Path(str(source) + '.pvar').open() as h:
        for line in h:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                fields = line.lstrip('#').split(); continue
            variants.append(dict(zip(fields, line.split())))
    bp, cm = read_source(source_map)
    selected, positions, seen = [], [], set()
    previous = -1
    for i, row in enumerate(variants):
        pos = int(row['POS'])
        if pos < previous:
            raise ValueError('Unsorted PVAR')
        previous = pos
        if pos not in seen:
            selected.append(i); positions.append(pos); seen.add(pos)
    if len(selected) < 2:
        raise ValueError('Need at least two distinct polymorphic sites')
    cols = np.asarray(columns, dtype=np.int32)
    alleles = np.empty(2*len(samples), dtype=np.int32)
    phase = np.empty(len(samples), dtype=np.uint8)
    with pgenlib.PgenReader(os.fsencode(str(source) + '.pgen')) as reader:
        with pgenlib.PgenWriter(os.fsencode(str(output) + '.pgen'), len(ids)//2,
                               variant_ct=len(selected), hardcall_phase_present=True) as writer:
            for index in selected:
                reader.read_alleles_and_phasepresent(index, alleles, phase)
                if np.any((alleles < 0) | (alleles > 1)):
                    raise ValueError(f'Missing/non-biallelic genotype at variant {index}; use --missing-max 0')
                het = alleles[::2] != alleles[1::2]
                if np.any(het & (phase == 0)):
                    raise ValueError(f'Unphased heterozygote at variant {index}')
                if chrom == 'X' and any(het[i] for i, row in enumerate(samples) if row['SEX'] == '1'):
                    raise ValueError('Heterozygous male non-PAR X genotype')
                writer.append_alleles(np.ascontiguousarray(alleles[cols]), all_phased=True)
    with Path(str(output) + '.pvar').open('w') as h:
        h.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for i in selected:
            r = variants[i]
            h.write(f"{chrom}\t{r['POS']}\t{r['ID']}\t{r['REF']}\t{r['ALT']}\n")
    with Path(str(output) + '.psam').open('w') as h:
        h.write('#IID\tSEX\n')
        for i in range(len(ids)//2):
            h.write(f'pair_{i}\t2\n')
    with Path(str(output) + '.haplotypes.tsv').open('w') as h:
        h.write('sample\thaplotype\n')
        for sample, hap in ids:
            h.write(f'{sample}\t{hap}\n')
    with Path(str(output) + '.map').open('w') as h:
        h.write('pos\tchr\tcM\n')
        for pos, genetic in zip(positions, np.interp(positions, bp, cm)):
            h.write(f'{pos}\t{chrom}\t{genetic:.12g}\n')
    atomic_json(str(output) + '.qc.json', dict(chromosome=chrom, haplotypes=len(ids),
                individuals=len(samples), variants=len(selected),
                duplicate_positions_removed=len(variants)-len(selected), errors=[]))


def prepare(chrom, prep, work, grid, replace, require=False):
    p = prep / ('chr' + chrom)
    source = Path(env('PFILE_DIR')) / ('chr' + chrom)
    var = Path(str(source) + '.pvar')
    if not var.is_file():
        var = Path(str(source) + '.pvar.zst')
    source_map = map_for(chrom)
    inputs = [str(source)+'.pgen', var, str(source)+'.psam', source_map,
              env('ANCESTRY_FILE'), __file__, HERE/'arg.py',
              grid/'make_arg_keep.py', grid/'build_argneedle_map.py', shutil.which('plink2')]
    for option in ('ARG_KEEP', 'EXTRACT'):
        if env(option, ''):
            inputs.append(env(option))
    request = dict(inputs=[stamp(x) for x in inputs], options={k:env(k) for k in
                   ('GRCH','ARG_MAF_MIN','ARG_GENO_MAX','SEED','MAX_INDIVIDUALS',
                    'ANCHORS_PER_POP','ARG_FULL','X_ODD_MALE')})
    marker = Path(str(p)+'.prepared.json')
    if not replace and verified(marker, request):
        print(f'SKIP verified Threads PGEN chr{chrom}', flush=True); return p
    if require:
        raise ValueError(f'Prepared chr{chrom} missing/stale; run prep_gen with the same options')
    marker.unlink(missing_ok=True)
    log = prep/'log'/f'arg.prepare.chr{chrom}.log'
    with tempfile.TemporaryDirectory(prefix=f'prepare.chr{chrom}.', dir=work) as tmp:
        tmp = Path(tmp); keep = tmp/'keep.tsv'; order = tmp/'order.tsv'
        cmd = [sys.executable, grid/'make_arg_keep.py', '--sample', str(source)+'.psam',
               '--ancestry', env('ANCESTRY_FILE'), '--out-keep', keep, '--out-panel', tmp/'panel.tsv',
               '--seed', env('SEED'), '--anchors-per-pop', env('ANCHORS_PER_POP'),
               '--max-individuals', env('MAX_INDIVIDUALS')]
        if env('ARG_FULL') == 'TRUE': cmd += ['--full']
        if env('ARG_KEEP', ''): cmd += ['--custom-keep', env('ARG_KEEP')]
        run(cmd, log)
        if chrom == 'X':
            needle_x_panel(keep, str(source)+'.psam', env('X_ODD_MALE'), tmp/'selection.qc.json')
            shutil.copy2(tmp/'selection.qc.json', str(p)+'.selection.qc.json')
        order.write_text('\n'.join(keep.read_text().splitlines()[1:])+'\n')
        subset = tmp/'subset'
        cmd = ['plink2', '--pfile', source]
        if str(var).endswith('.zst'): cmd += ['vzs']
        cmd += ['--keep', keep, '--indiv-sort', 'f', order, '--snps-only', 'just-acgt',
                '--max-alleles', '2', '--maf', env('ARG_MAF_MIN'), '--mac', '1',
                '--geno', env('ARG_GENO_MAX'), '--sort-vars', '--make-pgen',
                '--threads', env('THREADS'), '--out', subset]
        if chrom == 'X':
            lo, hi = {'37':(2699521,154931043), '38':(2781480,155701382)}[env('GRCH')]
            cmd += ['--chr', 'X', '--from-bp', str(lo), '--to-bp', str(hi)]
        if env('EXTRACT', ''): cmd += ['--extract', env('EXTRACT')]
        run(cmd, log)
        repack(subset, tmp/'ready', chrom, source_map)
        extensions = ('pgen','pvar','psam','haplotypes.tsv','map','qc.json')
        for ext in extensions:
            shutil.copy2(tmp/f'ready.{ext}', str(p)+'.'+ext+'.next')
            Path(str(p)+'.'+ext+'.next').replace(str(p)+'.'+ext)
        shutil.copy2(keep, str(p)+'.keep.tsv')
    files = [str(p)+'.'+ext for ext in (*extensions, 'keep.tsv')]
    if chrom == 'X': files.append(str(p)+'.selection.qc.json')
    checkpoint(marker, request, files)
    print(f'PREPARED Threads chr{chrom}: {p}', flush=True)
    return p


def export_tree(raw, dest, identities, sample_map, pgen_prefix):
    ts = tszip.decompress(raw)
    with Path(identities).open() as h:
        ids = [(r['sample'], int(r['haplotype'])) for r in csv.DictReader(h, delimiter='\t')]
    if ts.num_samples != len(ids):
        raise ValueError('Threads conversion changed sample count')
    tables = ts.dump_tables()
    offset = ts.metadata.get('offset', 0) if isinstance(ts.metadata, dict) else 0
    if offset:
        tables.sequence_length += offset
        tables.edges.left += offset; tables.edges.right += offset
        tables.sites.position += offset
    tables.time_units = 'generations'
    tables.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.metadata = dict(source_method='threads', source_arg_offset=offset,
                           offset=0, coordinate_system='chromosome_bp')
    tables.provenances.add_row(record=json.dumps(dict(software={'name':'gu-threads'},
                                parameters={'source_offset':offset, 'time_units':'generations',
                                            'ancestral_state_source':'tree_parsimony_not_external_AA'})))
    # arg-needle-lib 1.1.3's mutation converter assumes infinite sites; Threads
    # parsimony can have several mutations at one site. Build one site per PVAR
    # row and preserve parent links for recurrent/back mutations instead.
    topology = tables.tree_sequence()
    marginal = topology.first()
    genotypes = np.empty(len(ids), dtype=np.int32)
    with pgenlib.PgenReader(os.fsencode(str(pgen_prefix)+'.pgen')) as reader:
        with Path(str(pgen_prefix)+'.pvar').open() as variants:
            index = 0
            for line in variants:
                if line.startswith('#'): continue
                chrom, position, vid, ref, alt = line.split()[:5]
                position = int(position)
                if not 0 <= position < topology.sequence_length:
                    raise ValueError(f'Variant outside inferred ARG: {position}')
                reader.read_alleles(index, genotypes)
                marginal.seek(position)
                ancestor, mutations = marginal.map_mutations(genotypes, alleles=[ref,alt])
                site = tables.sites.add_row(position=position, ancestral_state=ancestor)
                base = len(tables.mutations)
                for mutation in mutations:
                    parent = tskit.NULL if mutation.parent == tskit.NULL else base + mutation.parent
                    tables.mutations.add_row(site=site,node=mutation.node,
                                             derived_state=mutation.derived_state,parent=parent)
                index += 1
                if index % 100000 == 0:
                    print(f'CONVERT mapped {index} sites: {dest}', flush=True)
    tables.sort()
    result = tables.tree_sequence()
    result.dump(str(dest)+'.next'); Path(str(dest)+'.next').replace(dest)
    Path(sample_map).write_text('tree_node_id\tsample\thaplotype\n' + ''.join(
        f'{int(node)}\t{name}\t{hap}\n' for node,(name,hap) in zip(result.samples(),ids)))
    return result, [(int(node), name, hap) for node,(name,hap) in zip(result.samples(),ids)]


def build(chrom, p, out, work, model, replace, action):
    instruction = out/'threads'/f'chr{chrom}.threads'
    marker = out/'threads'/f'chr{chrom}.infer.json'
    request = dict(prepared=stamp(str(p)+'.prepared.json'), model=model,
                   mode=env('THREADS_MODE'), fit=env('THREADS_FIT'),
                   mutation_rate=env('THREADS_MUTATION_RATE'), threads=env('THREADS'),
                   code=[stamp(__file__),stamp(HERE/'arg_threads_worker.py')], versions={k:importlib.metadata.version(k) for k in
                       ('threads_arg','arg-needle-lib','numpy','pgenlib','tskit')})
    if replace or not verified(marker, request):
        if action == 'convert': raise ValueError(f'No compatible inference checkpoint for chr{chrom}')
        marker.unlink(missing_ok=True)
        # Native Linux copies for repeated PGEN scans by all worker processes.
        local = work/f'chr{chrom}'
        for ext in ('pgen','pvar','psam','map'):
            shutil.copyfile(str(p)+'.'+ext, str(local)+'.'+ext)
        cmd = [sys.executable,HERE/'arg_threads_worker.py','infer','--pgen',str(local)+'.pgen','--map',str(local)+'.map',
               '--demography',model['path'],'--mutation_rate',env('THREADS_MUTATION_RATE'),
               '--mode',env('THREADS_MODE'),'--num_threads',env('THREADS'),
               '--out',str(instruction)+'.next']
        if env('THREADS_FIT') == 'TRUE': cmd += ['--fit_to_data']
        run(cmd, out/'log'/f'arg.infer.chr{chrom}.log')
        Path(str(instruction)+'.next').replace(instruction)
        checkpoint(marker, request, [instruction])
    else:
        print(f'SKIP verified Threads inference chr{chrom}', flush=True)
    if action == 'infer': return
    tree = out/'trees'/f'chr{chrom}.trees'; smap = tree.with_suffix('.sample_map.tsv')
    converted = out/'trees'/f'chr{chrom}.converted.json'
    conversion = dict(inference=stamp(marker), identities=stamp(str(p)+'.haplotypes.tsv'), code=stamp(__file__))
    if replace or not verified(converted, conversion):
        converted.unlink(missing_ok=True)
        raw = work/f'chr{chrom}.tsz'
        run(['threads','convert','--threads',instruction,'--tsz',raw],
            out/'log'/f'arg.convert.chr{chrom}.log')
        ts, rows = export_tree(raw, tree, str(p)+'.haplotypes.tsv', smap, p)
        stats = trace_validate(ts, rows, chrom)
        if ts.num_mutations < 1: raise ValueError('Threads tree has no mutations')
        atomic_json(tree.with_suffix('.arg_qc.json'), stats)
        checkpoint(converted, conversion, [tree,smap,tree.with_suffix('.arg_qc.json')])
        raw.unlink()
    if env('ARG_FORMAT') == 'trace':
        run([sys.executable,HERE/'arg.py','trace','--method','threads','--tree',tree,
             '--sample-map',smap,'--out-dir',out/'trace/threads','--chr',chrom],
            out/'log'/f'arg.trace.chr{chrom}.log')


def check(chrom, out):
    tree = out/'trees'/f'chr{chrom}.trees'; smap = tree.with_suffix('.sample_map.tsv')
    if not verified(out/'threads'/f'chr{chrom}.infer.json') or not verified(out/'trees'/f'chr{chrom}.converted.json'):
        raise ValueError(f'Missing/modified Threads checkpoint chr{chrom}')
    converted = json.loads((out/'trees'/f'chr{chrom}.converted.json').read_text())
    if converted['request']['inference'] != stamp(out/'threads'/f'chr{chrom}.infer.json'):
        raise ValueError(f'Tree does not match current Threads inference chr{chrom}')
    from arg import trace_read_map
    ts = tskit.load(tree)
    stats = trace_validate(ts, trace_read_map(smap, 'threads'), chrom)
    if ts.time_units != 'generations' or ts.num_mutations < 1:
        raise ValueError('Invalid Threads time units or missing mutations')
    if env('ARG_FORMAT') == 'trace':
        # Read-only validation: no logs or manifests are written here.
        subprocess.run([sys.executable,str(HERE/'arg.py'),'trace','--method','threads',
                        '--tree',str(tree),'--sample-map',str(smap),'--out-dir',str(out/'trace/threads'),
                        '--chr',chrom,'--check'], check=True)
    print(f'Threads CHECK PASS chr{chrom}: {stats}', flush=True)


def main():
    def stopped(signum, frame):
        raise SystemExit(128 + signum)
    for sig in (signal.SIGTERM, signal.SIGHUP):
        signal.signal(sig, stopped)
    out = Path(env('ARG_DIR')); prep = Path(env('GEN4ARG_DIR'))
    action = env('ARG_ACTION'); chrs = env('ARG_CHRS').split()
    if action == 'check':
        for chrom in chrs: check(chrom, out)
        return
    grid = Path(os.environ.get('REFGEN_GRID_ROOT', HERE.parents[1]/'grid'))/'f'
    sys.path.insert(0, str(grid))
    if not shutil.which('plink2'): raise ValueError('plink2 is missing from PATH')
    for chrom in chrs:
        source = Path(env('PFILE_DIR')) / ('chr' + chrom)
        for ext in ('pgen', 'psam'):
            stamp(str(source)+'.'+ext)
        if not any(Path(str(source)+'.'+ext).is_file() for ext in ('pvar','pvar.zst')):
            raise ValueError(f'Missing PVAR for chr{chrom}')
        map_for(chrom)
    stamp(env('ANCESTRY_FILE'))
    work = Path(env('ARG_SCRATCH'))/'refgen_threads'
    if str(work).startswith('/mnt/'): raise ValueError('--scratch must use native Linux storage')
    for path in (work, prep/'log',out/'threads',out/'trees',out/'log'):
        path.mkdir(parents=True, exist_ok=True)
    # Ray uses Unix sockets, whose paths must fit the ~108-byte OS limit.
    ray_tmp = tempfile.TemporaryDirectory(prefix='gu-ray-')
    os.environ['RAY_TMPDIR'] = ray_tmp.name
    import fcntl
    lock = (work/'.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    replace = env('REPLACE') == 'TRUE'
    mu = float(env('THREADS_MUTATION_RATE'))
    if not math.isfinite(mu) or mu <= 0: raise ValueError('--mutation-rate must be finite and positive')
    custom = env('THREADS_DEMOGRAPHY', '')
    if custom:
        demo = Path(custom).resolve()
    else:
        demo = work/'constant.haploid20000.demo'
        if not demo.exists(): demo.write_text('0\t20000\n')
    values = np.loadtxt(demo, ndmin=2)
    if (values.shape[1] != 2 or not np.isfinite(values).all() or values[0,0] != 0
        or np.any(np.diff(values[:,0]) <= 0) or np.any(values[:,1] <= 0)):
        raise ValueError('Demography must start at generation 0 with increasing times and positive haploid Ne')
    model = dict(path=str(demo), sha256=hashlib.sha256(demo.read_bytes()).hexdigest(),
                 description='custom haploid Ne' if custom else 'constant haploid Ne=20000 (diploid Ne=10000)')
    print(f'Threads model: {model}; workers={env("THREADS")}; chromosome jobs=1', flush=True)
    for chrom in chrs:
        p = prepare(chrom, prep, work, grid, replace, require=action in ('infer','convert'))
        if action != 'prepare': build(chrom, p, out, work, model, replace, action)
    if action not in ('prepare','infer'):
        atomic_json(out/'ARG_THREADS_BUILD.json',dict(method='threads',chromosomes=chrs,model=model))
        if env('ARG_FORMAT') == 'trace':
            (out/'trace/threads/ARG_TRACE_BUILD.tsv').write_text(
                f'method\tthreads\nformat\ttrace\nbuild\tGRCh{env("GRCH")}\nchromosomes\t{" ".join(chrs)}\n')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr, flush=True)
        raise
