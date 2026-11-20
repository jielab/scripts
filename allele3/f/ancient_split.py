"""Resumable per-chromosome export and conservative, allele-preserving liftover."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import fcntl
import gzip
import itertools
import json
import os
from pathlib import Path

import sys as _path_sys
_path_sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "0f"))
from workflow_paths import configured_path
import re
import shlex
import shutil
import subprocess
import tempfile

SCHEMA = 2


def run(argv, **kw):
    print('+ ' + shlex.join(map(str, argv)), flush=True)
    subprocess.run(list(map(str, argv)), check=True, **kw)


def state(path):
    p = Path(path).absolute()
    s = p.stat()
    return [str(p), s.st_size, s.st_mtime_ns]


def save(path, data):
    p = Path(path)
    temp = p.with_suffix(p.suffix + '.partial')
    temp.write_text(json.dumps(data, indent=2) + '\n')
    temp.replace(p)


def reusable(receipt, identity, outputs, replace):
    if replace or not receipt.is_file() or not all(p.is_file() for p in outputs):
        return False
    prior = json.loads(receipt.read_text())
    return prior.get('identity') == identity and prior.get('outputs') == [state(p) for p in outputs]


def record(receipt, identity, outputs, **kw):
    save(receipt, dict(identity=identity, outputs=[state(p) for p in outputs], **kw))


def canonical(chrom):
    chrom = chrom.removeprefix('chr')
    return {'23':'X', '24':'Y', '25':'XY', '26':'MT', 'M':'MT'}.get(chrom, chrom)


def plink(args, *extra):
    run(['plink2', '--threads', args.threads, '--memory', args.memory_mb, *extra])


def publish(source, target):
    """Stage on the destination filesystem before atomic publication."""
    target = Path(target)
    pending = target.with_name(target.name+'.partial')
    shutil.copyfile(source, pending)
    pending.replace(target)


def link(source, target):
    target = Path(target)
    pending = target.with_name(target.name+'.partial')
    pending.unlink(missing_ok=True)
    pending.symlink_to(os.path.relpath(source, target.parent))
    pending.replace(target)


def rows(path):
    with Path(path).open() as handle:
        for line in handle:
            if not line.startswith('#'): yield line.rstrip('\n').split('\t')


def fasta_aliases(references, chromosomes):
    """Resolve each requested chromosome to exactly one FASTA contig."""
    by_chrom = {}
    for name in references:
        by_chrom.setdefault(canonical(name), []).append(name)
    result = {}
    for chrom in sorted(chromosomes):
        matches = by_chrom.get(canonical(chrom), [])
        if len(matches) > 1:
            raise ValueError('Ambiguous FASTA chromosome aliases for ' + chrom + ': ' + ', '.join(matches))
        if matches:
            result[chrom] = matches[0]
    if not result:
        raise ValueError('No PVAR chromosomes match the FASTA contigs: ' + ', '.join(sorted(chromosomes)))
    return result


def lift_pvar(source, target, chain, fasta, work, report):
    """Keep all PVAR rows/alleles; change POS only, failures become -1."""
    import pysam
    bed, mapped, rejected, ordered = [work/n for n in ('source.bed','mapped.bed','unmapped.bed','ordered.bed')]
    total = 0
    chromosomes = set()
    with bed.open('w') as dst:
        for fields in rows(source):
            total += 1
            chrom, pos, _, ref = fields[:4]
            c = 'chr' + canonical(chrom)
            if c == 'chrMT': c = 'chrM'
            chromosomes.add(c)
            if int(pos) > 0 and re.fullmatch('[ACGTNacgtn]+', ref):
                dst.write(f'{c}\t{int(pos)-1}\t{int(pos)-1+len(ref)}\t{total}\t0\t+\n')
    # Read only the FASTA index here, before invoking liftOver or opening outputs.
    with pysam.FastaFile(str(fasta)) as refgen:
        reference_names = fasta_aliases(refgen.references, chromosomes)
    run(['liftOver', '-minMatch=1.0', bed, chain, mapped, rejected])
    with ordered.open('w') as dst:
        run(['sort', '-T', work, '-k4,4n', mapped], stdout=dst, env=dict(os.environ, LC_ALL='C'))
    counts = Counter(total=total)
    with pysam.FastaFile(str(fasta)) as refgen, ordered.open() as maps, Path(source).open() as src, Path(target).open('w') as dst, gzip.open(report, 'wt') as reject:
        groups = itertools.groupby((l.rstrip('\n').split('\t') for l in maps), key=lambda x:int(x[3]))
        current = next(groups, None)
        row = 0
        reject.write('source_row\tCHROM37\tPOS37\tID\tREF\tALT\treason\n')
        for line in src:
            if line.startswith('##reference='): continue
            if line.startswith('#CHROM'): dst.write('##reference=GRCh38\n')
            if line.startswith('#'): dst.write(line); continue
            row += 1
            fields = line.rstrip('\n').split('\t')
            chrom, pos, vid, ref, alt = fields[:5]
            expected = 'chr' + canonical(chrom)
            if expected == 'chrMT': expected = 'chrM'
            reason = 'unmapped_or_unsupported_REF'
            hits = []
            if current is not None and current[0] == row:
                hits = list(current[1]); current = next(groups, None)
            new_pos = -1
            if len(hits) > 1: reason = 'multiple_mappings'
            elif hits:
                c, start, end, _, _, strand = hits[0]
                start, end = int(start), int(end)
                reference_name = reference_names.get(c)
                if strand != '+': reason = 'reverse_strand'
                elif c != expected: reason = 'different_chromosome'
                elif end-start != len(ref): reason = 'changed_REF_span'
                elif reference_name is None or start < 0 or end > refgen.get_reference_length(reference_name): reason = 'outside_reference'
                elif refgen.fetch(reference_name, start, end).upper() != ref.upper(): reason = 'REF_mismatch'
                else:
                    new_pos = start+1
                    reason = 'accepted'
            counts[reason] += 1
            if new_pos == -1: reject.write(f'{row}\t{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{reason}\n')
            fields[1] = str(new_pos)
            dst.write('\t'.join(fields)+'\n')
    if sum(v for k,v in counts.items() if k != 'total') != total:
        raise RuntimeError('Liftover accounting mismatch')
    return dict(counts)


def chromosome(args, chrom, source_state, scratch):
    d37, d38 = map(Path, (args.out37, args.out38))
    label = 'chr' + ('M' if chrom == 'MT' else chrom)
    prefix37, prefix38 = d37/label, d38/label
    receipt37, receipt38 = d37/'.ancient_split', d38/'.ancient_split'
    work = scratch/label; work.mkdir()
    identity = dict(schema=SCHEMA, source=source_state, chrom=chrom, export='make-pgen', build=37)
    outputs = [Path(str(prefix37)+suffix) for suffix in ('.pgen','.pvar','.psam')]
    if not reusable(receipt37/(label+'.json'), identity, outputs, args.replace):
        staged = work/label
        plink(args, '--pfile', args.input_prefix, '--chr', chrom, '--make-pgen', '--out', staged)
        for dst,suffix in zip(outputs,('.pgen','.pvar','.psam')): publish(str(staged)+suffix,dst)
        publish(str(staged)+'.log',Path(args.runtime)/'logs'/(label+'.split.log'))
        record(receipt37/(label+'.json'), identity, outputs)
        for suffix in ('.pgen','.pvar','.psam'): Path(str(staged)+suffix).unlink()
    identity = dict(schema=SCHEMA, source=[state(p) for p in outputs], chain=state(args.chain), fasta=state(args.fasta), policy='unique_forward_same_chrom_exact_REF_POS_minus_1_fasta_alias_v1')
    rejected = receipt38/(label+'.rejected.tsv.gz')
    lifted = [Path(str(prefix38)+suffix) for suffix in ('.pgen','.pvar','.psam')]
    outputs = lifted+[rejected]
    if not reusable(receipt38/(label+'.json'), identity, outputs, args.replace):
        counts = lift_pvar(str(prefix37)+'.pvar', work/'lifted.pvar', args.chain, args.fasta, work, work/'rejected.tsv.gz')
        publish(work/'lifted.pvar',lifted[1]);publish(work/'rejected.tsv.gz',rejected)
        link(str(prefix37)+'.pgen',lifted[0]);link(str(prefix37)+'.psam',lifted[2])
        record(receipt38/(label+'.json'),identity,outputs,counts=counts)
    shutil.rmtree(work)
    return prefix38


def worker(args):
    if not os.environ.get('LSB_JOBID'): raise RuntimeError('Production worker requires an LSF allocation')
    if int(os.environ.get('LSB_DJOB_NUMPROC', '0')) < args.jobs*args.threads: raise RuntimeError('Insufficient allocated CPU slots')
    for name in ('plink2','liftOver','sort'):
        if not shutil.which(name): raise RuntimeError('Missing executable: '+name)
    for path in (args.out37,args.out38,str(Path(args.out37)/'.ancient_split'),str(Path(args.out38)/'.ancient_split'),str(Path(args.runtime)/'logs')):
        Path(path).mkdir(parents=True,exist_ok=True)
    with (Path(args.out38)/'.ancient_split.lock').open('a') as lock, tempfile.TemporaryDirectory(prefix='ancient-split-'+os.environ['LSB_JOBID']+'-',dir=args.scratch_root) as temp:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        scratch=Path(temp)
        progress=Path(args.out38)/'ancient.prepare.inprogress'
        progress.write_text('LSF job '+os.environ['LSB_JOBID']+'\n')
        source=[state(args.input_prefix+s) for s in ('.pgen','.pvar','.psam')]
        for p in (args.chain,args.fasta,args.fasta+'.fai'):state(p)
        if args.fasta.endswith(('.gz','.bgz','.bgzf')):state(args.fasta+'.gzi')
        available=set(canonical(r[0]) for r in rows(args.input_prefix+'.pvar'))
        order=list(map(str,range(1,23)))+['X','Y','MT']
        chosen=[c for c in order if c in available] if args.chromosomes=='auto' else [canonical(c) for c in args.chromosomes.split(',')]
        if not chosen or len(set(chosen))!=len(chosen) or set(chosen)-available:raise ValueError('Missing or duplicate requested chromosomes')
        if args.chromosomes=='auto' and available-set(order):raise ValueError('Unsupported contigs: '+str(available-set(order)))
        print('Chromosomes: '+','.join(chosen)+'; scratch='+str(scratch),flush=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            prefixes=list(pool.map(lambda c:chromosome(args,c,source,scratch),chosen))
        if source != [state(args.input_prefix+s) for s in ('.pgen','.pvar','.psam')]:raise ValueError('Source pfile changed during preparation')
        # Reconstruct original global row order, reusing the original whole PGEN.
        # This avoids merge/sort changes and preserves -1 sentinel rows exactly.
        if set(chosen)==available:
            identity=dict(schema=SCHEMA,source=source,parts=[state(str(p)+'.pvar') for p in prefixes])
            final=Path(args.out38)/'ancient'
            outputs=[Path(str(final)+s) for s in ('.pgen','.pvar','.psam')]
            receipt=Path(args.out38)/'ancient.prepare.json'
            if not reusable(receipt,identity,outputs,args.replace):
                receipt.unlink(missing_ok=True)
                streams={c:iter(rows(str(prefix)+'.pvar')) for c,prefix in zip(chosen,prefixes)}
                with open(args.input_prefix+'.pvar') as src,(scratch/'ancient.pvar').open('w') as dst:
                    for line in src:
                        if line.startswith('##reference='):continue
                        if line.startswith('#CHROM'):dst.write('##reference=GRCh38\n')
                        if line.startswith('#'):dst.write(line);continue
                        original=line.rstrip('\n').split('\t')
                        mapped=next(streams[canonical(original[0])])
                        if canonical(mapped[0])!=canonical(original[0]) or mapped[2:5]!=original[2:5]:raise ValueError('Split PVAR variant order/alleles changed')
                        original[1]=mapped[1];dst.write('\t'.join(original)+'\n')
                if any(next(it,None) is not None for it in streams.values()):raise ValueError('Extra split rows')
                publish(scratch/'ancient.pvar',outputs[1])
                link(args.input_prefix+'.pgen',outputs[0]);link(args.input_prefix+'.psam',outputs[2])
                record(receipt,identity,outputs,build='GRCh38',complete=True)
            progress.unlink()
            print('COMPLETE '+str(receipt),flush=True)
        else:
            print('Selected chromosomes complete; whole-panel marker retained until all chromosomes are prepared.',flush=True)


def main():
    work = configured_path('DIR_WORK')
    gen = configured_path('DIR_GEN')
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grch', type=int, choices=[38], default=38, help='Target build after liftOver [38]')
    p.add_argument('--input-prefix',default=gen+'/ancient/37/ancient')
    p.add_argument('--out37',default=gen+'/ancient/37')
    p.add_argument('--out38',default=gen+'/ancient/38')
    p.add_argument('--runtime',default=configured_path('ANALYSIS_ROOT')+'/allele3/38/ancient_split')
    p.add_argument('--chain',default=configured_path('DIR_FILES')+'/liftOver/hg19ToHg38.over.chain.gz')
    p.add_argument('--fasta',default=gen+'/1hgp/38/GRCH38.fasta.gz')
    p.add_argument('--chromosomes',default='auto')
    p.add_argument('--scratch-root',default='/tmp')
    for name,default in [('jobs',4),('threads',4),('memory-mb',8192),('job-memory-gb',64),('job-hours',120)]: p.add_argument('--'+name,type=int,default=default)
    p.add_argument('--queue',default='33');p.add_argument('--submit',choices=['hpc'],default='hpc')
    p.add_argument('--replace',choices=['TRUE','FALSE'],default='FALSE')
    p.add_argument('--plan-only',action='store_true');p.add_argument('--worker',action='store_true')
    args=p.parse_args()
    if min(args.jobs,args.threads,args.memory_mb,args.job_memory_gb,args.job_hours)<=0: p.error('Resources must be positive')
    if args.jobs*(args.memory_mb+2048)>args.job_memory_gb*1024: p.error('Total memory must cover concurrent PLINK workers and overhead')
    for name in ('input_prefix','out37','out38','runtime','chain','fasta','scratch_root'):setattr(args,name,str(Path(getattr(args,name)).absolute()))
    if Path(args.out37)==Path(args.out38):p.error('GRCh37 and GRCh38 outputs must differ')
    args.replace=args.replace=='TRUE'
    if args.worker: worker(args);return
    entry=Path(__file__).resolve().parents[1]/'0.ancient_split.sh'
    argv=['bash',str(entry),'--worker']
    for k,v in vars(args).items():
        if k in ('worker','plan_only'):continue
        argv.extend(['--'+k.replace('_','-'),str(v).upper() if isinstance(v,bool) else str(v)])
    command=shlex.join(argv)
    if args.plan_only: print(command);return
    logs=Path(args.runtime)/'jobs';logs.mkdir(parents=True,exist_ok=True)
    script='#!/bin/bash\nset -euo pipefail\nexec '+command+'\n'
    (logs/'worker.sh').write_text(script)
    subprocess.run(['bsub','-q',args.queue,'-n',str(args.jobs*args.threads),'-R',f'span[hosts=1] rusage[mem={args.job_memory_gb}GB/host]','-W',f'{args.job_hours}:00','-J','allele3.ancient_split','-o',str(logs/'%J.out'),'-e',str(logs/'%J.err')],input=script,text=True,check=True)

if __name__=='__main__':main()
