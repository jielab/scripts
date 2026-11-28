#!/usr/bin/env python3
"""Run and compare a reproducible, population-stratified HERV sequence census.

The sequence engine is census.py: phased GRCh38 variants condition a reference
sequence, and RepeatMasker annotates both reconstructed haplotypes.  This is
not a de novo assembly or a search for unknown inserted sequence.  Comparisons
below use homologous reference positions and an explicit callable denominator.
"""
from __future__ import annotations

import argparse
import bisect
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import fcntl
import filecmp
import gzip
import html
import json
import math
import os
from pathlib import Path
from output_layout import data_path
from herv_layout import product_path, ensure_layout
import queue
import random
import re
import shlex
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT.parent/'0f'))
from workflow_paths import configured_path

SCHEMA = 1
FAMILIES = ('ERV1', 'ERVK', 'ERVL', 'ERVL-MaLR')
POPULATIONS = {
    'ACB':'AFR', 'ASW':'AFR', 'ESN':'AFR', 'GWD':'AFR', 'LWK':'AFR', 'MSL':'AFR', 'YRI':'AFR',
    'CLM':'AMR', 'MXL':'AMR', 'PEL':'AMR', 'PUR':'AMR',
    'CDX':'EAS', 'CHB':'EAS', 'CHS':'EAS', 'JPT':'EAS', 'KHV':'EAS',
    'CEU':'EUR', 'FIN':'EUR', 'GBR':'EUR', 'IBS':'EUR', 'TSI':'EUR',
    'BEB':'SAS', 'GIH':'SAS', 'ITU':'SAS', 'PJL':'SAS', 'STU':'SAS',
}
CHROMS = tuple('chr'+str(i) for i in range(1,23))
LABEL = 'GRCh38 variant-conditioned reference reconstruction; not independent de novo assemblies'
HIGHLIGHT_START = '<!-- HERV_ENRICHMENT_START -->'
HIGHLIGHT_END = '<!-- HERV_ENRICHMENT_END -->'


def now():
    return datetime.now(timezone.utc).isoformat()


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name+f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(value, indent=2)+'\n')
    tmp.replace(path)


def read_json(path):
    return json.loads(Path(path).read_text())


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def file_state(path):
    p = Path(path)
    st = p.stat()
    return dict(path=str(p.resolve()), bytes=st.st_size, mtime_ns=st.st_mtime_ns)


@contextmanager
def file_lock(path, blocking=True):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | (0 if blocking else fcntl.LOCK_NB))
        yield


def chrom_key(chrom):
    c = str(chrom).removeprefix('chr')
    return int(c) if c.isdigit() else 99


def union(intervals):
    out=[]
    for a,b in sorted(intervals):
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1],b)
        else:
            out.append([int(a),int(b)])
    return out


def subtract(intervals, masks):
    """Subtract sorted, disjoint masks from sorted disjoint intervals."""
    out=[]
    j=0
    for a,b in intervals:
        while j<len(masks) and masks[j][1]<=a:
            j+=1
        cursor=a
        k=j
        while k<len(masks) and masks[k][0]<b:
            x,y=masks[k]
            if x>cursor:
                out.append([cursor,min(x,b)])
            cursor=max(cursor,y)
            if cursor>=b:
                break
            k+=1
        if cursor<b:
            out.append([cursor,b])
    return out


def overlap_lengths(loci, intervals):
    """Length of intersection for each sorted locus; loci may overlap."""
    import numpy as np
    answer=np.zeros(len(loci),dtype=np.uint32)
    ends=[b for a,b in intervals]
    for k,(a,b) in enumerate(loci):
        j=bisect.bisect_right(ends,a)
        n=0
        while j<len(intervals) and intervals[j][0]<b:
            x,y=intervals[j]
            n+=max(0,min(b,y)-max(a,x))
            j+=1
        answer[k]=n
    return answer


class HomologyMap:
    """Exact mapping through substitutions, insertions, deletions and delins.

    Equal-length replacements retain positional correspondence.  An inserted
    nucleotide has no reference counterpart.  A resolved deletion contributes
    an observed absence, not missing callability.  Unequal-length delins have
    ambiguous internal homology, so their replaced reference span is excluded.
    """
    def __init__(self, length, edits):
        self.length=int(length)
        self.segments=[]
        self.ambiguous=[]
        self.insertions=[]
        self.deletions=[]
        r=h=0
        def add(ha,hb,ra):
            if ha==hb:
                return
            if self.segments and self.segments[-1][1]==ha and self.segments[-1][2]+ha-self.segments[-1][0]==ra:
                old=self.segments[-1]
                self.segments[-1]=(old[0],hb,old[2])
            else:
                self.segments.append((ha,hb,ra))
        for a,b,alt in edits:
            a,b=int(a),int(b)
            if not 0<=r<=a<=b<=length:
                raise ValueError('Non-disjoint/out-of-range prepared sequence edits')
            add(h,h+a-r,r)
            h+=a-r
            width=b-a
            n=len(alt)
            if width==n:
                add(h,h+n,a)
            elif width==0:
                self.insertions.append([h,h+n])
            elif n==0:
                self.deletions.append([a,b])
            else:
                self.ambiguous.append([a,b])
            h+=n
            r=b
        add(h,h+length-r,r)
        self.haplotype_length=h+length-r
        self.ends=[x[1] for x in self.segments]

    def project(self, intervals):
        result=[]
        for a,b in intervals:
            if a<0 or b<a or b>self.haplotype_length:
                raise ValueError('RepeatMasker interval outside reconstructed haplotype')
            j=bisect.bisect_right(self.ends,a)
            while j<len(self.segments) and self.segments[j][0]<b:
                ha,hb,ra=self.segments[j]
                lo,hi=max(a,ha),min(b,hb)
                if hi>lo:
                    result.append([ra+lo-ha,ra+hi-ha])
                j+=1
        return union(result)


def read_panel(path):
    """The metadata file's /37/ location is not a genotype-build assertion."""
    rows={}
    with Path(path).open() as stream:
        header=None
        for line in stream:
            if not line.strip():
                continue
            tokens=line.strip().lstrip('#').split()
            if header is None:
                header=tokens
                required={'sample','pop','super_pop'}
                if not required<=set(header):
                    raise ValueError(f'Sample panel needs sample/pop/super_pop columns: {path}')
                continue
            if len(tokens)!=len(header):
                raise ValueError(f'Malformed sample-panel row: {line[:100]}')
            row=dict(zip(header,tokens))
            sid=row['sample']
            if sid in rows:
                raise ValueError(f'Duplicate sample in panel: {sid}')
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',sid):
                raise ValueError(f'Unsafe/invalid sample identifier: {sid}')
            pop=row['pop'].upper()
            if pop in POPULATIONS and row['super_pop'].upper()!=POPULATIONS[pop]:
                raise ValueError(f'Inconsistent population/superpopulation for {sid}: {pop}')
            rows[sid]=dict(sample=sid,pop=pop,super_pop=row['super_pop'].upper(),
                           sex=row.get('sex','NA'),fatherID=row.get('fatherID','NA'),
                           motherID=row.get('motherID','NA'))
    if not rows:
        raise ValueError(f'Empty sample panel: {path}')
    return rows


def read_psam(path):
    ids=[]
    with Path(path).open() as stream:
        header=stream.readline().strip().lstrip('#').split()
        if 'IID' not in header:
            raise ValueError(f'PSAM requires an IID column: {path}')
        j=header.index('IID')
        for line in stream:
            if not line.strip() or line.startswith('#'):
                continue
            row=line.split()
            if len(row)!=len(header):
                raise ValueError(f'Malformed PSAM row: {path}')
            ids.append(row[j])
    if len(ids)!=len(set(ids)):
        raise ValueError(f'Duplicate sample IDs in {path}')
    return set(ids)


def pgen_directory(path):
    """Accept the GRCh38 directory itself or its pfile child; no build mixing."""
    path=Path(path).expanduser().resolve()
    candidates=[path,path/'pfile']
    if path.name=='pfile':
        candidates.append(path.parent)
    def complete(root):
        return all((root/(c+'.pgen')).is_file() and (root/(c+'.psam')).is_file()
                   and any((root/(c+x)).is_file() for x in ('.pvar','.pvar.zst')) for c in CHROMS)
    for candidate in dict.fromkeys(candidates):
        if complete(candidate):
            return candidate
    raise ValueError(f'No complete chr1–22 PGEN/PVAR/PSAM set at {path} or its pfile directory; use GRCh38 inputs')


def available_samples(directory):
    all_ids=None
    counts={}
    for chrom in CHROMS:
        ids=read_psam(Path(directory)/(chrom+'.psam'))
        counts[chrom]=len(ids)
        all_ids=ids if all_ids is None else all_ids&ids
    if not all_ids:
        raise ValueError('No samples are present in every chromosome PSAM')
    return all_ids,counts


def select_samples(panel, available, populations, per_pop, seed):
    selected=[]
    shortage=[]
    for pop in sorted(populations):
        candidates=sorted(s for s,r in panel.items() if r['pop']==pop and s in available)
        if len(candidates)<per_pop:
            shortage.append(f'{pop}: {len(candidates)} available, need {per_pop}')
            continue
        # A population-specific generator makes the draw independent of input
        # row order and of the order in which populations are requested.
        rng=random.Random(f'{seed}:{pop}')
        for sid in sorted(rng.sample(candidates,per_pop)):
            selected.append(dict(panel[sid],candidate_count=len(candidates)))
    if shortage:
        raise ValueError('Cannot make the requested balanced panel: '+'; '.join(shortage))
    if len(selected)!=len(populations)*per_pop:
        raise ValueError('Balanced-panel sample-count mismatch')
    for i,row in enumerate(selected,1):
        row['sample_index']=i
    return selected


def requested_populations(value):
    pops=sorted(POPULATIONS) if value.lower()=='all' else [x.strip().upper() for x in value.split(',') if x.strip()]
    if not pops or len(pops)!=len(set(pops)) or set(pops)-set(POPULATIONS):
        raise ValueError('--populations must be all or distinct 1000 Genomes population codes')
    return sorted(pops)


def plan(args):
    out=args.out
    ensure_layout(out)
    with file_lock(out/'tmp/plan.lock'):
        panel=read_panel(args.sample_panel)
        pfile=pgen_directory(args.pfile_dir)
        available,counts=available_samples(pfile)
        populations=requested_populations(args.populations)
        selected=select_samples(panel,available,populations,args.per_pop,args.seed)
        identity=dict(seed=args.seed,per_pop=args.per_pop,populations=populations,build='GRCh38',
                      pfile_dir=str(pfile),sample_ids=[r['sample'] for r in selected],
                      sample_populations=[r['pop'] for r in selected])
        path=out/'qc/plan.json'
        if path.is_file():
            old=read_json(path)
            if old.get('identity')!=identity:
                raise ValueError('The frozen HERV panel differs from current seed/populations/available genotypes. '
                                 'Use a new --out for a different panel; existing samples are preserved.')
            selected=old['samples']
            result=old
        else:
            result=dict(schema=SCHEMA,status='planned',created=now(),identity=identity,
                        samples=selected,panel_metadata=file_state(args.sample_panel),
                        psam_counts=counts,n_samples=len(selected),n_populations=len(populations),
                        panel_note='The /37/ samples.txt is ID/population metadata only; all sequence inputs are GRCh38.',
                        interpretation=LABEL)
            save_json(path,result)
            write_tsv(out/'results/samples.tsv',selected,
                      ['sample_index','sample','pop','super_pop','sex','fatherID','motherID','candidate_count'])
        if args.reference_bed and not (out/'results/reference_herv.bed').is_file():
            catalogue=read_catalogue(args.reference_bed)
            write_catalogue(out/'results/reference_herv.bed',catalogue)
            save_json(out/'tmp/catalogue.json',dict(schema=SCHEMA,source='provided_reference_bed',
                       input=file_state(args.reference_bed),loci=len(catalogue),coordinate='GRCh38 BED0 half-open'))
    args.pfile_dir=pfile
    print(f'[herv-scan] frozen panel: {len(selected)} samples / {len(populations)} populations / '
          f'{args.per_pop} per population; {out}',flush=True)
    return result


def complete_sample(path, sample):
    """Lightweight completeness checks; no full genome/content rehash."""
    path=Path(path)
    try:
        run=read_json(path/'run.json')
        prepared=read_json(path/'prepared.json')
        plan_info=read_json(path/'plan.json')
        if run.get('status')!='complete' or prepared.get('sample')!=sample or prepared.get('chromosomes')!=list(CHROMS):
            return False
        for name in ('analysis.done.json','genome_fraction.tsv','genome_fraction.json','variant_summary.tsv'):
            if not (path/name).is_file():
                return False
        for chrom in CHROMS:
            for name in ('hap1.json.gz','hap2.json.gz','done.json'):
                if not (path/'variants'/chrom/name).is_file():
                    return False
        for chrom,label,a,b,_ in plan_info['descriptors']:
            folder=path/'chunks'/chrom/label/f'{a:09d}-{b:09d}'
            if not (folder/'done.json').is_file() or not (folder/'unions.bed.gz').is_file():
                return False
        return True
    except (OSError,ValueError,KeyError,TypeError):
        return False


def reference_settings_match(path,args):
    saved=read_json(Path(path)/'prepared.json')
    run=read_json(Path(path)/'run.json')['settings']
    # Check scientific settings and original paths without hashing the huge
    # source files again.  Relocation is handled by the census engine itself.
    if Path(saved['fasta']).resolve()!=args.fasta.resolve():
        raise ValueError(f'Completed sample used another FASTA: {path}; use its original path or a new --out')
    if int(run['chunk_bp'])!=args.chunk_bp or int(run['halo'])!=args.halo:
        raise ValueError(f'Completed sample has different chunk/halo settings: {path}')
    if Path(run['library']).resolve()!=args.library.resolve():
        raise ValueError(f'Completed sample used another repeat library: {path}; use a new --out')


def sample_directory(out, sample):
    out = Path(out)
    return out/'tmp/census'/sample


def projection_directory(out, sample_dir):
    return Path(out)/'tmp/projection'/Path(sample_dir).name


def comparable_sequence_settings(out,samples):
    """Compare recorded scientific identities without rehashing genome files."""
    keys=('method','fasta_sha256','library_sha256','repeatmasker_sha256',
          'rmblast_version','parameters','uncertain_flank','chunk_bp','halo')
    baseline=None
    first=None
    for row in samples:
        settings=read_json(sample_directory(out,row['sample'])/'run.json')['settings']
        current={key:settings.get(key) for key in keys}
        if baseline is None:
            baseline,first=current,row['sample']
        else:
            changed=[key for key in keys if current[key]!=baseline[key]]
            if changed:
                raise ValueError('Completed samples have incompatible sequence/annotation settings: '
                                 f'{first} vs {row["sample"]} ({", ".join(changed)}). '
                                 'Use one reference and repeat library for a cohort comparison.')
    return baseline


def census_command(args,row,cpus):
    sample=row['sample']
    folder=sample_directory(args.out,sample)
    command=[sys.executable,str(HERE/'census.py'),'--sample',sample,'--gen','1kg',
             '--out',str(folder),'--fasta',str(args.fasta),'--library',str(args.library),
             '--pfile-dir',str(args.pfile_dir),'--input-format','pgen',
             '--jobs',str(args.jobs),'--cores',str(len(cpus)),'--cpu-list',','.join(map(str,cpus)),
             '--chunk-bp',str(args.chunk_bp),'--halo',str(args.halo),
             '--scratch',str(args.out/'tmp/work'/sample),
             '--reference-cache',str(args.out/'tmp/reference_cache'),
             '--replace',args.replace]
    if args.repeatmasker:
        command.extend(['--repeatmasker',str(args.repeatmasker)])
    if args.pgen_python:
        command.extend(['--pgen-python',str(args.pgen_python)])
    return command


def run_sample(args,row,cpus):
    sample=row['sample']
    folder=sample_directory(args.out,sample)
    log_dir=args.out/'logs'
    log_dir.mkdir(parents=True,exist_ok=True)
    # This lock is outside the sample directory, so --replace TRUE cannot
    # delete the lock that serializes two requests for the same individual.
    with file_lock(args.out/'tmp/sample_locks'/f'{sample}.lock'):
        if args.replace=='FALSE' and complete_sample(folder,sample):
            reference_settings_match(folder,args)
            print(f'[herv-scan] REUSE {sample}: complete sequence/annotation outputs',flush=True)
            return dict(sample=sample,pop=row['pop'],status='reused',out=str(folder))
        cmd=census_command(args,row,cpus)
        env=dict(os.environ,ALLELE_ACTIVE_HERV_BUILD='38',HERV_BUILD='38',
                 ALLELE_REPLACE=args.replace)
        print(f'[herv-scan] START {sample} ({row["pop"]}); {len(cpus)} CPUs; '
              f'{args.jobs} annotation chunks; log={log_dir/(sample+".log")}',flush=True)
        with (log_dir/(sample+'.log')).open('a') as log:
            log.write('\n'+now()+' '+json.dumps(cmd)+'\n');log.flush()
            process=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env=env)
        if process.returncode or not complete_sample(folder,sample):
            raise RuntimeError(f'{sample} did not complete (exit={process.returncode}); '
                               f'see {log_dir/(sample+".log")}. Its completed chunks are preserved.')
        print(f'[herv-scan] DONE {sample}',flush=True)
        return dict(sample=sample,pop=row['pop'],status='complete',out=str(folder))


def cpu_slots(args,number):
    available=sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else list(range(os.cpu_count() or 1))
    needed=number*args.sample_cores
    if args.cores<needed or len(available)<needed:
        raise ValueError(f'{number} simultaneous samples × {args.sample_cores} CPUs need {needed} CPUs; '
                         f'--cores={args.cores}, allocated={len(available)}. Reduce --sample-jobs or request more cores.')
    return [available[i*args.sample_cores:(i+1)*args.sample_cores] for i in range(number)]


def execute_samples(args,result,single=False):
    rows=result['samples']
    if single:
        if args.sample_index is not None:
            if not 1<=args.sample_index<=len(rows):
                raise ValueError(f'--sample-index must be 1..{len(rows)}')
            rows=[rows[args.sample_index-1]]
        elif args.sample:
            panel=read_panel(args.sample_panel)
            available,_=available_samples(args.pfile_dir)
            if args.sample not in panel or args.sample not in available:
                raise ValueError(f'--sample {args.sample} must exist in the panel and all 22 GRCh38 PSAM files')
            rows=[panel[args.sample]]
        else:
            raise ValueError('sample requires --sample ID or --sample-index N')
    count=min(args.sample_jobs,len(rows))
    slots=queue.Queue()
    for cpus in cpu_slots(args,count):
        slots.put(cpus)
    def worker(row):
        cpus=slots.get()
        try:
            return run_sample(args,row,cpus)
        finally:
            slots.put(cpus)
    outcomes=[]
    with ThreadPoolExecutor(max_workers=count) as pool:
        futures={pool.submit(worker,row):row for row in rows}
        for future in as_completed(futures):
            row=futures[future]
            try:
                outcomes.append(future.result())
            except Exception as exc:
                outcomes.append(dict(sample=row['sample'],pop=row['pop'],status='failed',error=str(exc)))
                print('[herv-scan] '+str(exc),file=sys.stderr,flush=True)
    # Per-job receipts avoid concurrent HPC workers overwriting one another.
    for outcome in outcomes:
        save_json(args.out/'tmp/status'/(outcome['sample']+'.json'),dict(outcome,updated=now()))
    failed=[x for x in outcomes if x['status']=='failed']
    if failed:
        raise RuntimeError(f'{len(failed)} sample(s) failed. Resume with --replace FALSE; '
                           'failed/missing samples are never counted as zero HERV coverage.')
    # Make the reference catalogue available as soon as one sample has finished.
    ensure_catalogue(args,result)
    return outcomes


def read_catalogue(path):
    by_family=defaultdict(list)
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as stream:
        for line in stream:
            if not line.strip() or line.startswith(('#','track','browser')):
                continue
            row=line.rstrip().split('\t')
            if len(row)<4:
                raise ValueError('--reference-bed needs BED4 with family|locus_id labels')
            chrom='chr'+row[0].removeprefix('chr')
            if chrom not in CHROMS:
                continue
            family=row[3].split('|')[0].removeprefix('family:')
            if family not in FAMILIES and len(row)>4 and row[4] in FAMILIES:
                family=row[4]
            if family not in FAMILIES:
                raise ValueError(f'Unknown HERV family in {path}: {row[3:5]}; expected {", ".join(FAMILIES)}')
            a,b=int(row[1]),int(row[2])
            if a<0 or b<=a:
                raise ValueError(f'Invalid BED0 interval in {path}: {row[:3]}')
            by_family[chrom,family].append([a,b])
    rows=[]
    for (chrom,family),intervals in by_family.items():
        for a,b in union(intervals):
            rows.append(dict(chrom=chrom,start0=a,end0=b,family=family))
    rows.sort(key=lambda r:(chrom_key(r['chrom']),r['start0'],r['end0'],r['family']))
    for i,row in enumerate(rows,1):
        row['locus_id']=f'HERV{i:07d}'
    if not rows:
        raise ValueError(f'No autosomal ERV1/ERVK/ERVL/MaLR regions in {path}')
    return rows


def write_catalogue(path,rows):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    with tmp.open('w') as stream:
        stream.write('#chrom\tstart0\tend0\tfamily|locus_id\n')
        for row in rows:
            stream.write(f'{row["chrom"]}\t{row["start0"]}\t{row["end0"]}\t{row["family"]}|{row["locus_id"]}\n')
    # Rewriting identical stability BEDs would invalidate the enrichment
    # engine's lightweight input identity on every resumed cohort summary.
    if path.is_file() and filecmp.cmp(path,tmp,shallow=False):
        tmp.unlink()
        return
    tmp.replace(path)


def annotation_intervals(sample_dir,chrom,label):
    found=defaultdict(list)
    chunks=sorted((Path(sample_dir)/'chunks'/chrom/label).glob('*/unions.bed.gz'))
    if not chunks:
        raise ValueError(f'Missing {chrom}/{label} RepeatMasker chunk annotations: {sample_dir}')
    for path in chunks:
        with gzip.open(path,'rt') as stream:
            for line in stream:
                c,a,b,group=line.rstrip().split('\t')
                family=group.removeprefix('family:')
                if group.startswith('family:') and family in FAMILIES:
                    if c!=chrom:
                        raise ValueError(f'Unexpected chromosome in {path}')
                    found[family].append([int(a),int(b)])
    return {family:union(found[family]) for family in FAMILIES}


def ensure_catalogue(args,result):
    path=args.out/'results/reference_herv.bed'
    with file_lock(args.out/'tmp/catalogue.lock'):
        if path.is_file():
            return read_catalogue(path)
        completed=[r for r in result['samples'] if complete_sample(sample_directory(args.out,r['sample']),r['sample'])]
        if not completed:
            return None
        first=completed[0]['sample']
        sample_dir=sample_directory(args.out,first)
        rows=[]
        for chrom in CHROMS:
            for family,intervals in annotation_intervals(sample_dir,chrom,'reference').items():
                for a,b in intervals:
                    rows.append(dict(chrom=chrom,start0=a,end0=b,family=family))
        rows.sort(key=lambda r:(chrom_key(r['chrom']),r['start0'],r['end0'],r['family']))
        for i,row in enumerate(rows,1):
            row['locus_id']=f'HERV{i:07d}'
        if not rows:
            raise ValueError('Completed reference annotation has no HERV-family regions')
        write_catalogue(path,rows)
        save_json(args.out/'tmp/catalogue.json',dict(schema=SCHEMA,source='RepeatMasker reference control',
                   sample=first,source_run=str(sample_dir),loci=len(rows),
                   coordinate='GRCh38 BED0 half-open; same-family contiguous annotation unions',
                   note='Regions are annotation unions, not claims of intact provirus boundaries.'))
        print(f'[herv-scan] reference catalogue: {len(rows):,} regions from {first}',flush=True)
        return rows


def noncallable_reference(args,chrom):
    """Cache reference non-ACGT spans once, shared by every personal model."""
    folder=args.out/'tmp/reference_gaps'
    folder.mkdir(parents=True,exist_ok=True)
    bed=folder/(chrom+'.bed.gz')
    receipt=folder/(chrom+'.json')
    with file_lock(folder/(chrom+'.lock')):
        state=file_state(args.fasta)
        if receipt.is_file() and bed.is_file() and read_json(receipt).get('fasta')==state:
            meta=read_json(receipt)
        else:
            command=['samtools','faidx',str(args.fasta),chrom]
            process=subprocess.run(command,capture_output=True,text=True,check=True)
            sequence=''.join(process.stdout.splitlines()[1:]).upper()
            if not sequence:
                raise ValueError(f'Empty reference contig {chrom}')
            gaps=[[match.start(),match.end()] for match in re.finditer('[^ACGT]+',sequence)]
            with gzip.open(bed,'wt') as stream:
                for a,b in gaps:
                    stream.write(f'{chrom}\t{a}\t{b}\n')
            meta=dict(fasta=state,length=len(sequence),non_acgt_bp=sum(b-a for a,b in gaps))
            save_json(receipt,meta)
        with gzip.open(bed,'rt') as stream:
            gaps=[[int(r[1]),int(r[2])] for r in (line.split() for line in stream)]
    return meta['length'],gaps


def projected_annotation(args,sample_dir,chrom,hap,length,reference_gaps):
    """Persist a small, reference-coordinate projection for summary resumption."""
    folder=projection_directory(args.out,sample_dir)
    folder.mkdir(parents=True,exist_ok=True)
    base=folder/(chrom+'.'+hap)
    bed=Path(str(base)+'.bed.gz')
    masks_path=Path(str(base)+'.uncallable.bed.gz')
    receipt=Path(str(base)+'.json')
    model_path=Path(sample_dir)/'variants'/chrom/(hap+'.json.gz')
    identity=dict(schema=SCHEMA,model=file_state(model_path),
                  annotation=file_state(Path(sample_dir)/'genome_fraction.tsv'),
                  reference_gaps=file_state(args.out/'tmp/reference_gaps'/(chrom+'.bed.gz')))
    if receipt.is_file() and bed.is_file() and masks_path.is_file() and read_json(receipt).get('identity')==identity:
        meta=read_json(receipt)
    else:
        with gzip.open(model_path,'rt') as stream:
            model=json.load(stream)
        mapper=HomologyMap(length,model['edits'])
        masks=union(reference_gaps+model['masks']+mapper.ambiguous)
        annotations=annotation_intervals(sample_dir,chrom,hap)
        projected={}
        counts=[]
        for family,hits in annotations.items():
            mapped=mapper.project(hits)
            usable=subtract(mapped,masks)
            projected[family]=usable
            counts.append(dict(family=family,raw_haplotype_annotation_bp=sum(b-a for a,b in hits),
                               mapped_reference_annotation_bp=sum(b-a for a,b in mapped),
                               callable_annotation_bp=sum(b-a for a,b in usable),
                               unprojected_annotation_bp=sum(b-a for a,b in hits)-sum(b-a for a,b in mapped)))
        with gzip.open(bed,'wt') as stream:
            for family,intervals in projected.items():
                for a,b in intervals:
                    stream.write(f'{chrom}\t{a}\t{b}\t{family}\n')
        with gzip.open(masks_path,'wt') as stream:
            for a,b in masks:
                stream.write(f'{chrom}\t{a}\t{b}\n')
        meta=dict(identity=identity,counts=counts,reference_length=length,
                  uncallable_bp=sum(b-a for a,b in masks),resolved_deletion_bp=sum(b-a for a,b in mapper.deletions),
                  ambiguous_delins_bp=sum(b-a for a,b in mapper.ambiguous),
                  interpretation='Resolved deletions count as absence; unknown/masked spans are not zero calls; '
                                 'inserted bases without homologous reference coordinates are reported separately.')
        save_json(receipt,meta)
    hits=defaultdict(list)
    with gzip.open(bed,'rt') as stream:
        for line in stream:
            c,a,b,family=line.split()
            hits[family].append([int(a),int(b)])
    with gzip.open(masks_path,'rt') as stream:
        masks=[[int(r[1]),int(r[2])] for r in (line.split() for line in stream)]
    return hits,masks,meta


def measure_sample(args,row,catalogue,chrom_groups):
    import numpy as np
    folder=sample_directory(args.out,row['sample'])
    n=len(catalogue)
    evaluable=np.zeros(n,dtype=np.uint8)
    retained=np.zeros(n,dtype=np.uint8)
    callable_bp=np.zeros(n,dtype=np.uint32)
    annotation_bp=np.zeros(n,dtype=np.uint32)
    projections=[]
    for chrom in CHROMS:
        if chrom not in chrom_groups:
            continue
        length,gaps=noncallable_reference(args,chrom)
        for hap in ('hap1','hap2'):
            hits,masks,meta=projected_annotation(args,folder,chrom,hap,length,gaps)
            for counts in meta['counts']:
                projections.append(dict(sample=row['sample'],pop=row['pop'],chrom=chrom,haplotype=hap,**counts,
                                        uncallable_reference_bp=meta['uncallable_bp'],
                                        resolved_deletion_bp=meta['resolved_deletion_bp']))
            for family,indices in chrom_groups[chrom].items():
                coords=[(catalogue[i]['start0'],catalogue[i]['end0']) for i in indices]
                lengths=np.array([b-a for a,b in coords],dtype=np.uint32)
                if any(b>length for a,b in coords):
                    raise ValueError(f'Reference HERV catalogue exceeds {chrom} length; check GRCh38')
                called=lengths-overlap_lengths(coords,masks)
                observed=overlap_lengths(coords,hits[family])
                if np.any(observed>called):
                    raise ValueError('Projected annotation exceeds the callable reference denominator')
                ok=(called>0)&(called/lengths>=args.min_callable)
                present=ok&(np.divide(observed,called,out=np.zeros(len(called),dtype=float),where=called>0)>=args.min_region_coverage)
                evaluable[indices]+=ok.astype(np.uint8)
                retained[indices]+=present.astype(np.uint8)
                callable_bp[indices]+=called
                annotation_bp[indices]+=observed
    return evaluable,retained,callable_bp,annotation_bp,projections


def bootstrap_ratio(numerator,denominator,populations,reps,seed):
    """Resample individuals within each fixed population; haplotypes stay paired."""
    import numpy as np
    numerator=np.asarray(numerator,dtype=float)
    denominator=np.asarray(denominator,dtype=float)
    if numerator.ndim!=1 or numerator.shape!=denominator.shape or len(numerator)!=len(populations):
        raise ValueError('Bootstrap arrays and population labels differ')
    estimate=float(numerator.sum()/denominator.sum()) if denominator.sum()>0 else None
    if estimate is None:
        return estimate,None,None
    groups=[[i for i,p in enumerate(populations) if p==pop] for pop in sorted(set(populations))]
    rng=np.random.default_rng(seed)
    sums_n=np.zeros(reps)
    sums_d=np.zeros(reps)
    for ids in groups:
        chosen=rng.choice(ids,size=(reps,len(ids)),replace=True)
        sums_n+=numerator[chosen].sum(axis=1)
        sums_d+=denominator[chosen].sum(axis=1)
    good=sums_d>0
    if not good.any():
        return estimate,None,None
    ci=np.quantile(sums_n[good]/sums_d[good],[.025,.975])
    return estimate,float(ci[0]),float(ci[1])


def stability_statistics(catalogue,samples,evaluable,retained,callable_bp,annotation_bp,args,planned_n,planned_pops):
    import numpy as np
    lengths=np.array([r['end0']-r['start0'] for r in catalogue],dtype=np.uint32)
    families=np.array([r['family'] for r in catalogue])
    fully=evaluable==2
    n_evaluable_samples=fully.sum(axis=0)
    n_evaluable_haps=evaluable.sum(axis=0)
    n_retained=retained.sum(axis=0)
    occupancy=np.divide(n_retained,n_evaluable_haps,out=np.zeros(len(catalogue),dtype=float),where=n_evaluable_haps>0)
    eval_pops=np.zeros(len(catalogue),dtype=np.uint8)
    for pop in planned_pops:
        ids=[i for i,r in enumerate(samples) if r['pop']==pop]
        if ids:
            eval_pops+=(fully[ids].sum(axis=0)>0).astype(np.uint8)
    sample_min=math.ceil(args.stable_min_sample_fraction*planned_n)
    pop_min=math.ceil(args.stable_min_population_fraction*len(planned_pops))
    stable=(n_evaluable_samples>=sample_min)&(eval_pops>=pop_min)&(occupancy>=args.stable_retention)
    locus_rows=[]
    for j,r in enumerate(catalogue):
        locus_rows.append(dict(r,length_bp=int(lengths[j]),planned_samples=planned_n,completed_samples=len(samples),
             evaluable_samples=int(n_evaluable_samples[j]),evaluable_haplotypes=int(n_evaluable_haps[j]),
             retained_haplotypes=int(n_retained[j]),retention_fraction=float(occupancy[j]) if n_evaluable_haps[j] else 'NA',
             carrier_samples=int(((retained[:,j]>0)&fully[:,j]).sum()),
             both_retained_samples=int(((retained[:,j]==2)&fully[:,j]).sum()),
             evaluable_populations=int(eval_pops[j]),planned_populations=len(planned_pops),
             callable_reference_fraction=float(callable_bp[:,j].sum()/(2*int(lengths[j])*len(samples))),
             stable=int(stable[j])))
    family_rows=[]
    sample_rows=[]
    population_rows=[]
    bootstrap_rows=[]
    pops=[r['pop'] for r in samples]
    for family in FAMILIES:
        keep=families==family
        if not keep.any():
            continue
        weights=lengths[keep].astype(float)
        evaluated=(evaluable[:,keep]*weights).sum(axis=1)
        kept=(retained[:,keep]*weights).sum(axis=1)
        called=callable_bp[:,keep].sum(axis=1,dtype=np.uint64)
        annotated=annotation_bp[:,keep].sum(axis=1,dtype=np.uint64)
        possible=np.repeat(2*weights.sum(),len(samples))
        metrics={'retained_region_bp_fraction':(kept,evaluated),
                 'callable_reference_bp_fraction':(called,possible),
                 'annotated_callable_bp_fraction':(annotated,called)}
        for metric,(num,den) in metrics.items():
            estimate,lo,hi=bootstrap_ratio(num,den,pops,args.bootstrap,args.seed)
            bootstrap_rows.append(dict(family=family,metric=metric,estimate=estimate,ci_low=lo,ci_high=hi,
                 bootstrap_replicates=args.bootstrap,n_samples=len(samples),n_populations=len(set(pops)),
                 method='individual resampling within each fixed population; paired haplotypes'))
        b={r['metric']:r for r in bootstrap_rows if r['family']==family}
        family_rows.append(dict(family=family,n_reference_regions=int(keep.sum()),reference_bp=int(weights.sum()),
             n_stable_regions=int(stable[keep].sum()),stable_reference_bp=int(lengths[keep&stable].sum(dtype=np.uint64)),
             stable_region_fraction=float(stable[keep].mean()),retained_fraction=b['retained_region_bp_fraction']['estimate'],
             retention_ci_low=b['retained_region_bp_fraction']['ci_low'],retention_ci_high=b['retained_region_bp_fraction']['ci_high'],
             callable_fraction=b['callable_reference_bp_fraction']['estimate'],
             annotation_fraction=b['annotated_callable_bp_fraction']['estimate'],completed_samples=len(samples),planned_samples=planned_n))
        for i,r in enumerate(samples):
            sample_rows.append(dict(sample=r['sample'],pop=r['pop'],super_pop=r['super_pop'],family=family,
                 reference_bp=int(weights.sum()),callable_bp=int(called[i]),annotated_bp=int(annotated[i]),
                 evaluable_region_bp=int(evaluated[i]),retained_region_bp=int(kept[i]),
                 callable_fraction=float(called[i]/possible[i]) if possible[i] else 'NA',
                 annotation_fraction=float(annotated[i]/called[i]) if called[i] else 'NA',
                 retention_fraction=float(kept[i]/evaluated[i]) if evaluated[i] else 'NA'))
        for pop in planned_pops:
            ids=[i for i,r in enumerate(samples) if r['pop']==pop]
            den=evaluated[ids].sum() if ids else 0
            population_rows.append(dict(pop=pop,super_pop=POPULATIONS[pop],family=family,n_samples=len(ids),
                 retained_fraction=float(kept[ids].sum()/den) if den else 'NA',
                 callable_fraction=float(called[ids].sum()/possible[ids].sum()) if ids else 'NA',
                 annotation_fraction=float(annotated[ids].sum()/called[ids].sum()) if ids and called[ids].sum() else 'NA'))
    return locus_rows,family_rows,sample_rows,population_rows,bootstrap_rows,stable


def pairwise_jaccard(catalogue,samples,evaluable,retained):
    """Jaccard over reference regions evaluable on both haplotypes of both people."""
    import numpy as np
    families=np.array([r['family'] for r in catalogue])
    lengths=np.array([r['end0']-r['start0'] for r in catalogue],dtype=np.uint64)
    result=[]
    for family in FAMILIES:
        keep=families==family
        if not keep.any():
            continue
        valid=evaluable[:,keep]==2
        present=retained[:,keep]>0
        weights=lengths[keep]
        for i in range(len(samples)):
            for j in range(i+1,len(samples)):
                both=valid[i]&valid[j]
                intersection=both&present[i]&present[j]
                joined=both&(present[i]|present[j])
                n_intersection=int(intersection.sum());n_union=int(joined.sum())
                bp_union=int(weights[joined].sum())
                result.append(dict(sample1=samples[i]['sample'],sample2=samples[j]['sample'],
                    pop1=samples[i]['pop'],pop2=samples[j]['pop'],family=family,
                    evaluable_regions=int(both.sum()),intersection_regions=n_intersection,union_regions=n_union,
                    jaccard=n_intersection/n_union if n_union else 'NA',
                    length_weighted_jaccard=float(weights[intersection].sum()/bp_union) if bp_union else 'NA',
                    status='complete' if n_union else 'no_evaluable_union'))
    return result


def sample_percentages(args,rows):
    output=[]
    for row in rows:
        path=sample_directory(args.out,row['sample'])/'genome_fraction.tsv'
        values={(r['label'],r['category']):r for r in read_tsv(path)}
        record=dict(sample=row['sample'],pop=row['pop'],super_pop=row['super_pop'],sex=row.get('sex','NA'))
        for category,short in [('ERV_with_MaLR','erv_with_malr'),('ERV_without_MaLR','erv_without_malr'),('all_LTR','all_ltr')]:
            personal=values.get(('diploid',category),{})
            reference=values.get(('reference_callable_diploid',category),{})
            for key in ('percent_genome','percent_non_N'):
                record[short+'_'+key]=float(personal[key]) if personal.get(key) else 'NA'
                record['reference_'+short+'_'+key]=float(reference[key]) if reference.get(key) else 'NA'
            if personal.get('percent_non_N') and reference.get('percent_non_N'):
                record[short+'_delta_percent_non_N']=float(personal['percent_non_N'])-float(reference['percent_non_N'])
            else:
                record[short+'_delta_percent_non_N']='NA'
        record['status']='complete'
        output.append(record)
    return output


def make_figures(out,samples,families,pairwise,percentages):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=out/'figures'
    folder.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),constrained_layout=True)
    names=[r['family'] for r in families]
    y=np.array([r['retained_fraction'] or 0 for r in families])
    low=np.array([r['retention_ci_low'] if r['retention_ci_low'] is not None else v for r,v in zip(families,y)])
    high=np.array([r['retention_ci_high'] if r['retention_ci_high'] is not None else v for r,v in zip(families,y)])
    axes[0].errorbar(range(len(names)),y,yerr=np.maximum(0,np.vstack((y-low,high-y))),fmt='o',color='#177c91',capsize=4)
    axes[0].set_xticks(range(len(names)),names)
    axes[0].set_ylabel('Retained annotation, evaluable region bp')
    axes[0].set_ylim(max(0,min(low)-.03),min(1.005,max(high)+.01))
    axes[0].set_title('Family stability with population-stratified bootstrap')
    axes[1].bar(names,[r['stable_region_fraction'] for r in families],color='#d8943d')
    axes[1].set_ylim(0,1.02)
    axes[1].set_ylabel('Fraction of reference regions meeting stability rules')
    axes[1].set_title('Stable regions')
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y',alpha=.15)
    fig.savefig(folder/'family_stability.png',dpi=180)

    plt.close(fig)
    ids={r['sample']:i for i,r in enumerate(samples)}
    fig,axes=plt.subplots(2,2,figsize=(12,10),constrained_layout=True)
    for family,ax in zip(FAMILIES,axes.flat):
        matrix=np.full((len(samples),len(samples)),np.nan)
        np.fill_diagonal(matrix,1)
        for r in pairwise:
            if r['family']==family and r['jaccard']!='NA':
                i,j=ids[r['sample1']],ids[r['sample2']]
                matrix[i,j]=matrix[j,i]=r['jaccard']
        finite=matrix[np.isfinite(matrix)]
        vmin=max(0,float(finite.min())-.002) if len(finite) else 0
        im=ax.imshow(matrix,cmap='viridis',vmin=vmin,vmax=1)
        ax.set_title(family)
        labels=[r['pop']+' / '+r['sample'] for r in samples]
        ax.set_xticks(range(len(samples)),labels,rotation=90,fontsize=4)
        ax.set_yticks(range(len(samples)),labels,fontsize=4)
        fig.colorbar(im,ax=ax,shrink=.7,label='Mutually evaluable-region Jaccard')
    fig.savefig(folder/'pairwise_jaccard.png',dpi=180)
    plt.close(fig)
    valid=[r for r in percentages if r['erv_with_malr_percent_genome']!='NA']
    if valid:
        fig,ax=plt.subplots(figsize=(12,4),constrained_layout=True)
        x=np.arange(len(valid))
        ax.scatter(x,[r['erv_with_malr_percent_genome'] for r in valid],label='Personal model',color='#177c91',s=25)
        ax.scatter(x,[r['reference_erv_with_malr_percent_genome'] for r in valid],label='Matched reference',marker='_',color='#b75d39')
        ax.set_xticks(x,[r['pop']+' / '+r['sample'] for r in valid],rotation=90,fontsize=6)
        ax.set_ylabel('Autosomal ERV + MaLR (% model length)')
        ax.set_title('Sequence census and matched reference control')
        ax.legend(frameon=False)
        ax.spines[['top','right']].set_visible(False)
        fig.savefig(folder/'sample_coverage.png',dpi=180)
        plt.close(fig)


def html_table(rows,columns):
    def fmt(v):
        if v is None:
            return 'NA'
        if isinstance(v,float):
            return f'{v:.5g}'
        return str(v)
    return '<table><thead><tr>'+''.join('<th>'+html.escape(label)+'</th>' for key,label in columns)+'</tr></thead><tbody>'+''.join(
        '<tr>'+''.join('<td>'+html.escape(fmt(row.get(key,'NA')))+'</td>' for key,label in columns)+'</tr>' for row in rows)+'</tbody></table>'


def write_report(args,result,samples,families,percentages,pairwise,missing,catalogue_meta):
    family_table=html_table(families,[('family','Family'),('n_reference_regions','Reference regions'),
        ('n_stable_regions','Stable regions'),('retained_fraction','Retained fraction'),
        ('retention_ci_low','95% CI lower'),('retention_ci_high','95% CI upper'),('callable_fraction','Callable fraction')])
    sample_table=html_table(percentages,[('sample','Sample'),('pop','Population'),
        ('erv_with_malr_percent_genome','ERV + MaLR %'),('erv_without_malr_percent_genome','ERV without MaLR %'),
        ('reference_erv_with_malr_percent_genome','Matched reference %'),
        ('erv_with_malr_delta_percent_non_N','Delta non-N percentage points')])
    missing_text=('<p class="note">Incomplete panel: '+html.escape(', '.join(missing))+'. These individuals are missing, not zero-HERV samples.</p>') if missing else ''
    source=html.escape(str(catalogue_meta.get('source','reference catalogue')))
    page=f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>HERV sequence census across populations</title><style>
body{{font:16px/1.6 system-ui,sans-serif;color:#1b303a;background:#f5f8f9;margin:0}}main{{max-width:1120px;margin:auto;padding:36px 24px}}
h1,h2{{line-height:1.2}}h1{{font-size:34px}}h2{{margin-top:42px}}.kpis{{display:flex;gap:20px;flex-wrap:wrap}}.kpi{{background:white;border-radius:10px;padding:18px 28px;box-shadow:0 2px 12px #1231}}
.kpi strong{{display:block;font-size:30px;color:#177c91}}.note{{border-left:4px solid #d8943d;background:#fff7e9;padding:15px}}table{{border-collapse:collapse;width:100%;font-size:13px;background:white}}th,td{{text-align:left;padding:8px;border-bottom:1px solid #dde4e7}}th{{background:#e9f2f5}}img{{max-width:100%;background:white;border-radius:8px}}a{{color:#126d83}}code{{background:#e4edef;padding:2px 4px}}.scroll{{overflow-x:auto}}
</style><main><h1>HERV sequence census across populations</h1>
<p>{html.escape(LABEL)}. Both haplotypes were reconstructed and annotated with RepeatMasker for every completed individual.</p>
<div class="kpis"><div class="kpi"><strong>{len(samples)} / {result['n_samples']}</strong>completed individuals</div>
<div class="kpi"><strong>{len(set(r['pop'] for r in samples))} / {result['n_populations']}</strong>populations represented</div>
<div class="kpi"><strong>{sum(r['n_stable_regions'] for r in families):,}</strong>stable family regions</div></div>{missing_text}
<h2>What was measured</h2><p>RepeatMasker calls were projected through each individual's resolved sequence edits to homologous GRCh38 positions. Equal-length substitutions retain positional correspondence. Resolved deletions count as absence. Reference non-ACGT sequence, unresolved variants and ambiguous internal homology of unequal-length replacements are excluded from the callable denominator. Inserted bases without a reference counterpart are reported separately in <a href="../results/projection_summary.tsv">projection_summary.tsv</a>.</p>
<p>Evaluability describes this reconstruction model. A sites-only genotype file does not establish independent read-level coverage at every other reference base; positions absent from the variant input retain the reference sequence.</p>
<p>Reference catalogue: {source}. Regions are contiguous same-family annotation intervals, not claims of complete or independently discovered proviruses. The <a href="../results/reference_herv.bed">full reference BED</a> and <a href="../results/stable_loci.bed">stable BED</a> use 0-based half-open coordinates and <code>family|locus_id</code> labels.</p>
<h2>Family-specific stability</h2><p>A haplotype is evaluable when at least {args.min_callable:.0%} of a reference region is callable; its annotation is retained when at least {args.min_region_coverage:.0%} of those callable bases are annotated in the same family. A stable region requires at least {args.stable_min_sample_fraction:.0%} of the planned samples evaluable on both haplotypes, representation of at least {args.stable_min_population_fraction:.0%} of the planned populations, and retention in at least {args.stable_retention:.0%} of evaluable haplotypes. Bootstrap CIs resample individuals within each fixed population; the two haplotypes of an individual stay paired.</p>
<div class="scroll">{family_table}</div><img src="../figures/family_stability.png" alt="Family stability and stable-region fractions">
<h2>Agreement between individuals</h2><p>Jaccard compares reference-region carrier sets only where both haplotypes of both individuals are evaluable. Unknown regions are excluded. Because all models start from a shared reference and known variants, high agreement is expected and does not independently establish completeness of HERV insertions.</p>
<img src="../figures/pairwise_jaccard.png" alt="Pairwise Jaccard for four HERV families">
<h2>Sequence coverage and the approximately 8% quantity</h2><p>ERV + MaLR combines ERV1, ERVK, ERVL and ERVL-MaLR. ERV without MaLR excludes the last family; all LTR is a separate, broader class. Comparisons to the matched reference use each sample's uncertainty masks.</p>
<img src="../figures/sample_coverage.png" alt="Personal sequence percentages and matched reference"><div class="scroll">{sample_table}</div>
<h2>Outputs and reproducibility</h2><p><a href="../results/samples.tsv">Frozen sample panel</a> · <a href="../results/sample_summary.tsv">Sample percentages</a> · <a href="../results/family_stability.tsv">Family stability</a> · <a href="../results/population_summary.tsv">Population summary</a> · <a href="../results/stability_bootstrap.tsv">Bootstrap results</a> · <a href="../results/pairwise_jaccard.tsv">Pairwise comparisons</a> · <a href="../results/locus_stability.tsv.gz">All region-level results</a></p>
<p>The selection seed is {args.seed}; candidates were sorted before sampling. The population file under a directory named <code>37</code> supplies identifiers and labels only. Genotypes, reference sequence and all reported genomic coordinates are GRCh38. No raw reads or independent long-read assemblies were analysed by this workflow.</p></main></html>'''
    (args.out/'reports/herv_scan.html').write_text(page)


def enrichment_windows(value):
    try:
        windows=sorted(set(int(x.strip()) for x in str(value).split(',') if x.strip()))
    except ValueError as exc:
        raise argparse.ArgumentTypeError('Use comma-separated positive window sizes') from exc
    if not windows or windows[0]<1:
        raise argparse.ArgumentTypeError('Use comma-separated positive window sizes')
    return ','.join(map(str,windows))


def positive_float(value):
    number=float(value)
    if not math.isfinite(number) or number<=0:
        raise argparse.ArgumentTypeError('Must be a finite positive number')
    return number


HIGHLIGHT_TARGETS = ('gnomad', 'allofus')


def parse_highlight_targets(value):
    value = value.strip().lower()
    if value == 'auto':
        return value
    names = value.split(',')
    names = [name.strip() for name in names]
    if any(name not in HIGHLIGHT_TARGETS for name in names):
        raise argparse.ArgumentTypeError('Supported targets: '+','.join(HIGHLIGHT_TARGETS)+
                                         '; separate selected sources with commas')
    return ','.join(dict.fromkeys(names))


def highlight_targets(args):
    """Resolve explicitly named derivative stores; never scan source genomes."""
    target = parse_highlight_targets(args.target)
    if args.scan_dir:
        if ',' in target:
            raise ValueError('--scan-dir names one target store; use --summary-root for multiple targets')
        if target == 'auto':
            manifest = data_path(args.scan_dir, 'qc')/'manifest.json'
            target = read_json(manifest).get('target') if manifest.is_file() else args.scan_dir.name.lower()
            if target not in HIGHLIGHT_TARGETS:
                raise ValueError('Cannot infer a supported target from --scan-dir; specify --target')
        return [(target, args.scan_dir)]
    targets = HIGHLIGHT_TARGETS if target == 'auto' else target.split(',')
    return [(name, args.summary_root/name) for name in targets]


def enrichment_arguments(args,target,scan,bed,stable):
    command=['enrich_test','--ref','herv','--target',target,'--scan-dir',str(scan),
             '--out',str(args.out/'enrichment'/target),'--herv-scan-dir',str(args.out),
             '--herv-bed',str(bed),'--region-windows',args.enrichment_windows,
             '--highlight-fold',str(args.enrichment_fold),'--fasta',str(args.fasta),
             '--threshold-kind',args.threshold_kind,
             '--min-ac',str(args.min_ac),'--min-an',str(args.min_an),'--min-af',str(args.min_af),
             '--seed',str(args.seed),'--replace',args.replace,
             '--jobs',str(getattr(args,'enrichment_jobs',0))]
    if stable:
        command+=['--stable-herv-bed',str(stable)]
    if args.no_context:
        command.append('--no-context')
    for excluded in args.exclude_bed:
        command+=['--exclude-bed',str(excluded)]
    return command


def highlight_retry_command(args,target,scan,bed):
    command=['bash',str(ROOT/'3.allele3_scan.sh'),'highlights','--target',target,
             '--scan-dir',str(scan),'--out',str(args.out),'--herv-bed',str(bed),
             '--enrichment-fold',str(args.enrichment_fold),'--enrichment-windows',args.enrichment_windows,
             '--threshold-kind',args.threshold_kind,
             '--min-ac',str(args.min_ac),'--min-an',str(args.min_an),'--min-af',str(args.min_af),
             '--fasta',str(args.fasta),'--seed',str(args.seed),'--replace',args.replace,
             '--submit',os.environ.get('ALLELE_EXECUTION_MODE','hpc')]
    if args.herv_input_dir:
        command+=['--dir-in-herv_scan',str(args.herv_input_dir)]
    if args.no_context:
        command.append('--no-context')
    for excluded in args.exclude_bed:
        command+=['--exclude-bed',str(excluded)]
    return shlex.join(command)


def highlight_section(args,status):
    cutoff=f'{args.enrichment_fold:g}'
    body=[f'<section id="enrichment-highlights"><h2>All regions with conditional O/E ≥ {cutoff}</h2>',
        '<p>The genome-window and individual HERV-region results use completed target scan stores. '
        'They measure multi-ALT enrichment among eligible observed variant positions, conditional on the stated '
        'background model. The full reference HERV catalogue is used; sample stability is an annotation, '
        'not a requirement for including a region.</p>',
        '<p><a href="herv_highlights.html">Dedicated enrichment report</a> · '
        '<a href="../qc/enrichment_status.json">Target availability and exact commands</a></p>']
    for row in status['targets']:
        target=html.escape(row['target'])
        body.append(f'<h3>{target}: {html.escape(row["status"])}</h3>')
        if row['status']=='complete':
            base='enrichment/'+row['target']
            count=row.get('highlight_rows')
            text=f'{count:,} exported rows satisfy O/E ≥ {cutoff}.' if count is not None else f'All O/E ≥ {cutoff} rows are exported.'
            if count==0:
                text=f'The completed analysis found no estimable region with O/E ≥ {cutoff} under these settings.'
            body.append('<p>'+text+' Low support and mask sensitivity are reported rather than silently removed. '
                        'Overlapping windows and individual HERV regions are separate rows, not a count of independent discoveries.</p>')
            if row.get('supported_highlight_regions') is not None:
                body.append(f'<p>Support-qualified rows: {row["supported_highlight_regions"]:,}; '
                            f'rows annotated as HERV-related: {row.get("herv_related_highlight_regions",0):,}. '
                            'A whole-genome window that overlaps HERV sequence does not show that HERV bases cause the window effect.</p>')
            links=[('reports/highlights.html','Full regional report'),('results/highlights.tsv.gz','All matching rows (TSV)'),
                   ('results/highlights.bed.gz','All matching regions (BED)'),
                   ('results/highlights_supported.tsv.gz','Support-qualified subset'),
                   ('results/regions.tsv.gz','All predefined region tests'),
                   (('reports/region_enrichment.html' if Path(row.get('report') or 'highlights.html').with_name('region_enrichment.html').is_file() else 'reports/enrich_test.html'),'Family enrichment')]
            body.append('<p>'+' · '.join(f'<a href="../{base}/{path}">{label}</a>' for path,label in links)+'</p>')
            for path,label in [('figures/highlights_genome.png','Genome-wide positions of enriched regions'),
                               ('figures/highlights_ranked.png','Ranked enriched regions')]:
                if (args.out/base/path).is_file():
                    body.append(f'<img src="../{base}/{path}" alt="{label}">')
        else:
            body.append('<p class="note">'+html.escape(row['reason'])+
                        f' This is not evidence that no ≥ {cutoff}-fold region exists.</p>')
            if row.get('prepare_command'):
                body.append('<p>Prepare or complete this target store:</p><pre><code>'+html.escape(row['prepare_command'])+'</code></pre>')
            body.append('<p>Run the enrichment after its inputs are available:</p><pre><code>'+html.escape(row['run_command'])+'</code></pre>')
    comparison=status.get('comparison',{})
    if comparison and comparison.get('status')!='not_requested':
        body.append('<h3>Across-target comparison</h3>')
        if comparison.get('report') and Path(comparison['report']).is_file():
            link=html.escape(os.path.relpath(comparison['report'],args.out/'reports'),quote=True)
            body.append(f'<p><a href="{link}">Complete regional concordance report</a>. '
                        'The comparison selects the union of regions reaching the fold cutoff in either cohort; '
                        'the complete region tables supply the other cohort\'s estimate even below the cutoff. '
                        'Missing or unsupported results remain NA. '
                        'Count-based and AF-only qualification may differ and are explicitly assessed before interpreting agreement.</p>')
        else:
            body.append('<p class="note">'+html.escape(comparison.get('reason','Comparison is pending.'))+'</p>')
    body.append('<p>Sequence reconstruction and the approximately 8% ERV + MaLR coverage are separate measurements '
                'from conditional multi-ALT enrichment. A large regional O/E is not the family-wide effect, '
                'an insertion-discovery result, or proof of balancing selection.</p></section>')
    return '\n'.join(body)


def update_highlight_reports(args,status):
    section=highlight_section(args,status)
    style='body{font:16px/1.6 system-ui,sans-serif;max-width:1150px;margin:36px auto;padding:0 24px;color:#173244}img{max-width:100%}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#edf3f5;padding:14px}.note{background:#fff5df;padding:14px;border-left:4px solid #d8943d}a{color:#126d83}'
    full='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HERV regional enrichment highlights</title><style>'+style+'</style><main><h1>HERV regional enrichment highlights</h1><p><a href="herv_scan.html">Sequence census and report navigation</a></p>'+section+'</main></html>'
    (args.out/'reports/herv_highlights.html').write_text(full)
    page=args.out/'reports/herv_scan.html'
    block=HIGHLIGHT_START+section+HIGHLIGHT_END
    if page.is_file():
        content=page.read_text()
        if HIGHLIGHT_START in content and HIGHLIGHT_END in content:
            content=content.split(HIGHLIGHT_START,1)[0]+block+content.split(HIGHLIGHT_END,1)[1]
        elif '</main>' in content:
            content=content.replace('</main>',block+'</main>',1)
        else:
            content+=block
    else:
        content='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HERV analysis</title><style>'+style+'</style><main><h1>HERV analysis</h1><p>The population sequence-census summary is not available here yet. The enrichment analysis below can run independently from existing reference HERV intervals and completed target stores.</p>'+block+'</main></html>'
    page.write_text(content)


def highlights(args,automatic=False):
    """Use reusable target summaries; never reconstruct a genome or run RM."""
    args.out.mkdir(parents=True,exist_ok=True)
    with file_lock(args.out/'tmp/highlights.lock'):
        source=args.herv_input_dir or args.out
        bed=args.herv_bed or source/'results/reference_herv.bed'
        stable=source/'results/stable_loci.bed'
        stable=stable if stable.is_file() else None
        rows=[]
        for target,scan in highlight_targets(args):
            scan=Path(scan)
            destination=args.out/'enrichment'/target
            source=data_path(scan, 'qc')/'manifest.json'
            row=dict(target=target,status='pending',scan_dir=str(scan),out=str(destination),
                     run_command=highlight_retry_command(args,target,scan,bed),highlight_rows=None)
            missing=[]
            if automatic:
                missing.append('Enrichment is scheduled after summary scans by 3.allele3_scan.sh; use highlights for a manual retry.')
            if not bed.is_file():
                missing.append(f'Reference HERV catalogue is missing: {bed}. Complete one HERV sample or provide --herv-bed FILE.')
            if not source.is_file():
                missing.append(f'Completed target store is missing: {source}.')
            else:
                try:
                    meta=read_json(source)
                    if meta.get('status')!='complete':
                        missing.append(f'Target scan is not complete (status={meta.get("status","unknown")}).')
                    elif meta.get('target')!=target:
                        missing.append(f'Target manifest is for {meta.get("target")}, not {target}.')
                    elif meta.get('build',meta.get('genome_build'))!='GRCh38':
                        missing.append('Target scan is not GRCh38.')
                except (OSError,ValueError) as exc:
                    missing.append(f'Target manifest cannot be read: {exc}')
            if not args.no_context and not args.fasta.is_file():
                missing.append(f'Trinucleotide context FASTA is missing: {args.fasta}; supply --fasta or explicitly choose --no-context.')
            elif not args.no_context:
                indexes=[Path(str(args.fasta)+'.fai')]
                if str(args.fasta).endswith(('.gz','.bgz')):
                    indexes.append(Path(str(args.fasta)+'.gzi'))
                for index in indexes:
                    if not index.is_file():
                        missing.append(f'Pre-existing context FASTA index is missing: {index}.')
            for excluded in args.exclude_bed:
                if not excluded.is_file():
                    missing.append(f'Exclusion BED is missing: {excluded}.')
            if missing:
                row['reason']=' '.join(missing)
                if not source.is_file() or any('Target scan is not complete' in reason for reason in missing):
                    row['prepare_command']=shlex.join(['bash',str(ROOT/'3.allele3_scan.sh'),'scan','--target',target,
                        '--out',str(scan),'--submit',os.environ.get('ALLELE_EXECUTION_MODE','hpc'),'--replace','FALSE'])
                print(f'[herv-highlights] PENDING {target}: {row["reason"]}',flush=True)
            else:
                command=enrichment_arguments(args,target,scan,bed,stable)
                row['engine_arguments']=command
                try:
                    from research_tests import main as run_enrichment
                    print(f'[herv-highlights] START {target}: {destination}',flush=True)
                    from output_layout import prepare
                    prepare(destination)
                    manifest=run_enrichment(command)
                    needed=('results/highlights.tsv.gz','results/highlights.bed.gz','results/highlights_supported.tsv.gz',
                            'results/regions.tsv.gz','reports/highlights.html','qc/region_highlights.json',
                            'figures/highlights_genome.png','figures/highlights_ranked.png')
                    absent=[name for name in needed if not (destination/name).is_file()]
                    if absent:
                        raise RuntimeError('The enrichment builder did not publish: '+', '.join(absent))
                    regional=read_json(destination/'qc/region_highlights.json')
                    if regional.get('status')!='complete':
                        raise RuntimeError('Regional enrichment report is not complete')
                    counts=regional.get('counts',{})
                    row['highlight_rows']=counts.get('highlight_regions')
                    if row['highlight_rows'] is None:
                        with gzip.open(destination/'results/highlights.tsv.gz','rt') as stream:
                            row['highlight_rows']=sum(1 for _ in csv.DictReader(stream,delimiter='\t'))
                    for key in ('tested_regions','supported_highlight_regions','herv_related_highlight_regions'):
                        if key in counts:
                            row[key]=counts[key]
                    row.update(status='complete',source_manifest=file_state(source),
                               enrichment_manifest=str(destination/'qc/manifest.json'),
                               report=str(destination/'reports/highlights.html'),
                               output_files=[str(destination/name) for name in needed],
                               engine_status=manifest.get('status') if isinstance(manifest,dict) else 'complete')
                    print(f'[herv-highlights] DONE {target}: {row["highlight_rows"]:,} exported rows with O/E >= {args.enrichment_fold:g}',flush=True)
                except (OSError,ValueError,RuntimeError,KeyError,ImportError) as exc:
                    row.update(status='failed',reason=str(exc))
                    print(f'[herv-highlights] FAILED {target}: {exc}',file=sys.stderr,flush=True)
            rows.append(row)
        count=sum(row['status']=='complete' for row in rows)
        state='complete' if count==len(rows) else 'partial' if count else 'failed' if any(row['status']=='failed' for row in rows) else 'pending'
        comparison=dict(status='not_requested')
        if set(row['target'] for row in rows)=={'gnomad','allofus'}:
            comparison=dict(status='pending',reason='Both target enrichment analyses must be complete before the across-target comparison is generated; unavailable target results are not zero.')
            if count==2:
                try:
                    from highlight_compare import compare_regions
                    from output_layout import prepare
                    prepare(args.out/'enrichment/comparison')
                    comparison=compare_regions(args.out/'enrichment/gnomad',args.out/'enrichment/allofus',
                                               args.out/'enrichment/comparison')
                except (OSError,ValueError,RuntimeError,KeyError,ImportError) as exc:
                    comparison=dict(status='failed',reason=str(exc))
                    print(f'[herv-highlights] COMPARISON FAILED: {exc}',file=sys.stderr,flush=True)
        status=dict(schema='herv-enrichment-highlights-1',status=state,updated=now(),requested_target=args.target,
                    completed_targets=count,fold_cutoff=args.enrichment_fold,windows_bp=args.enrichment_windows,
                    reference_herv_bed=str(bed),stable_herv_bed=str(stable) if stable else None,targets=rows,
                    comparison=comparison,
                    report=str(args.out/'reports/herv_highlights.html'),
                    interpretation='Conditional multi-ALT enrichment among observed eligible variant positions; all threshold-crossing rows retained with support labels.')
        save_json(args.out/'qc/enrichment_status.json',status)
        update_highlight_reports(args,status)
        if not automatic and any(row['status']=='failed' for row in rows):
            raise RuntimeError('One or more HERV enrichment analyses failed; see enrichment_status.json and herv_highlights.html')
        return status


def summarize(args,result=None):
    import numpy as np
    result=result or read_json(args.out/'qc/plan.json')
    with file_lock(args.out/'tmp/summary.lock'):
        missing=[]
        samples=[]
        completion=[]
        for row in result['samples']:
            good=complete_sample(sample_directory(args.out,row['sample']),row['sample'])
            completion.append(dict(sample=row['sample'],pop=row['pop'],status='complete' if good else 'missing_or_incomplete'))
            (samples if good else missing).append(row if good else row['sample'])
        write_tsv(args.out/'qc/sample_completion.tsv',completion,['sample','pop','status'])
        if missing and not args.allow_partial:
            raise ValueError(f'{len(missing)} of {result["n_samples"]} planned samples are incomplete: '+
                             ', '.join(missing)+'. Resume run, or explicitly use summarize --allow-partial for a labelled interim report.')
        if not samples:
            raise ValueError('No completed sequence census is available to summarize')
        scientific_settings=comparable_sequence_settings(args.out,samples)
        catalogue=ensure_catalogue(args,result)
        if catalogue is None:
            raise ValueError('Reference HERV catalogue is unavailable')
        groups=defaultdict(lambda:defaultdict(list))
        for i,row in enumerate(catalogue):
            groups[row['chrom']][row['family']].append(i)
        n,m=len(samples),len(catalogue)
        evaluable=np.zeros((n,m),dtype=np.uint8)
        retained=np.zeros((n,m),dtype=np.uint8)
        called=np.zeros((n,m),dtype=np.uint32)
        annotated=np.zeros((n,m),dtype=np.uint32)
        projections=[]
        for i,row in enumerate(samples):
            print(f'[herv-scan] compare {i+1}/{n}: {row["sample"]}',flush=True)
            ev,pr,ca,an,projection=measure_sample(args,row,catalogue,groups)
            evaluable[i],retained[i],called[i],annotated[i]=ev,pr,ca,an
            projections.extend(projection)
        locus,family,individual,population,bootstrap,stable=stability_statistics(
            catalogue,samples,evaluable,retained,called,annotated,args,result['n_samples'],result['identity']['populations'])
        pairs=pairwise_jaccard(catalogue,samples,evaluable,retained)
        percentages=sample_percentages(args,samples)
        for name,rows in [('locus_stability.tsv.gz',locus),('family_stability.tsv',family),
                          ('sample_family_stability.tsv',individual),('population_summary.tsv',population),
                          ('stability_bootstrap.tsv',bootstrap),('pairwise_jaccard.tsv',pairs),
                          ('sample_summary.tsv',percentages),('projection_summary.tsv',projections)]:
            if rows:
                write_tsv(product_path(args.out,name),rows,list(rows[0]))
        write_catalogue(args.out/'results/stable_loci.bed',[r for j,r in enumerate(catalogue) if stable[j]])
        # A compact machine-readable matrix supports alternative thresholds
        # without another sequence reconstruction/RepeatMasker run.
        np.savez_compressed(args.out/'results/locus_calls.npz',
                            samples=np.array([r['sample'] for r in samples]),
                            locus_ids=np.array([r['locus_id'] for r in catalogue]),
                            evaluable_haplotypes=evaluable,retained_haplotypes=retained,
                            callable_bp=called,annotated_bp=annotated)
        make_figures(args.out,samples,family,pairs,percentages)
        cat_meta=read_json(args.out/'tmp/catalogue.json') if (args.out/'tmp/catalogue.json').is_file() else {}
        write_report(args,result,samples,family,percentages,pairs,missing,cat_meta)
        summary=dict(schema=SCHEMA,status='partial' if missing else 'complete',completed=now(),
            planned_samples=result['n_samples'],completed_samples=n,missing_samples=missing,
            planned_populations=result['n_populations'],completed_populations=len(set(r['pop'] for r in samples)),
            reference_regions=m,stable_regions=int(stable.sum()),
            settings={k:getattr(args,k) for k in ('min_callable','min_region_coverage','stable_retention',
                       'stable_min_sample_fraction','stable_min_population_fraction','bootstrap','seed')},
            sequence_annotation_settings=scientific_settings,
            stable_bed=str(args.out/'results/stable_loci.bed'),reference_bed=str(args.out/'results/reference_herv.bed'),
            report=str(args.out/'reports/herv_scan.html'),interpretation=LABEL,
            missing_rule='Missing samples and unresolved spans are excluded from HERV presence/absence; no zero imputation.',
            bootstrap='Individuals resampled within fixed populations, paired haplotypes retained; not uncertainty over all human populations.')
        save_json(args.out/'qc/summary.json',summary)
        if not missing:
            from herv_export import export_cohort
            summary['export']=export_cohort(args.out, io_workers=args.cores)
        enrichment=highlights(args,automatic=True)
        summary['enrichment']=dict(status=enrichment['status'],completed_targets=enrichment['completed_targets'],
                                   status_file=str(args.out/'qc/enrichment_status.json'),report=enrichment['report'])
        save_json(args.out/'qc/summary.json',summary)
        if not missing:
            from herv_publish import publish
            summary['publication']=publish(args.out,args.archive)
        print(json.dumps(summary,indent=2),flush=True)
        return summary


def positive(value):
    n=int(value)
    if n<1:
        raise argparse.ArgumentTypeError('Must be a positive integer')
    return n


def fraction(value):
    x=float(value)
    if not 0<x<=1:
        raise argparse.ArgumentTypeError('Must be in (0, 1]')
    return x


def parser():
    p=argparse.ArgumentParser(description=__doc__,allow_abbrev=False,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('command',nargs='?',choices=('plan','run','sample','summarize','export') if os.environ.get('HERV_COHORT_ONLY')=='1' else ('plan','run','sample','summarize','export','highlights'),default='run')
    p.add_argument('--out',type=Path,default=Path(os.environ.get('HERV_OUT',os.environ.get('ALLELE_HERV_SCAN_STORE',configured_path('ANALYSIS_ROOT','allele3/2.herv_scan')))))
    p.add_argument('--dir-out-herv_scan',dest='archive',type=Path,default=Path(os.environ.get('ALLELE_HERV_ARCHIVE',configured_path('DIR_GEN','1kg/38/herv_scan'))),help='Permanent validated products; runtime stays under --out in analysis')
    p.add_argument('--sample-panel',type=Path,default=Path(os.environ.get('HERV_SAMPLE_PANEL',configured_path('DIR_GEN','1kg/38/samples.txt'))))
    p.add_argument('--pfile-dir',type=Path,default=Path(os.environ.get('HERV_PFILE_DIR',configured_path('DIR_GEN','1kg/38/pfile'))))
    p.add_argument('--fasta',type=Path,default=Path(os.environ.get('HERV_FASTA',configured_path('DIR_GEN','1hgp/38/GRCH38.fasta.gz'))))
    p.add_argument('--library',type=Path,default=Path(os.environ.get('HERV_LIBRARY',configured_path('DIR_ANNOT','dfam_4.0_human/human.curated.fa'))))
    p.add_argument('--reference-bed',type=Path,help='Optional GRCh38 family BED4 catalogue; otherwise same-parameter RepeatMasker reference control')
    p.add_argument('--repeatmasker',type=Path)
    p.add_argument('--pgen-python',type=Path)
    p.add_argument('--dat',choices=['1kg'],default='1kg')
    p.add_argument('--build',choices=['38'],default='38',help='Genotype/reference build; the panel metadata path does not change this')
    p.add_argument('--populations',default='all',help='all = all 26 standard populations; a named subset is an explicitly labelled pilot')
    p.add_argument('--per-pop',type=positive,default=2)
    p.add_argument('--seed',type=int,default=20261115)
    which=p.add_mutually_exclusive_group()
    which.add_argument('--sample')
    which.add_argument('--sample-index',type=positive,help='1-based row of the frozen balanced panel, for independent HPC jobs')
    p.add_argument('--enrichment-jobs',type=int,default=0,help='Highlights chromosome workers; 0 uses allocated cores')
    p.add_argument('--sample-jobs',type=positive,default=2,help='Maximum simultaneous individuals for local run')
    p.add_argument('--sample-cores',type=positive,default=16)
    p.add_argument('--cores',type=positive,default=32,help='Total allocated CPU budget for this cohort process')
    p.add_argument('--jobs',type=positive,default=8,help='Concurrent RepeatMasker chunks inside each individual')
    p.add_argument('--chunk-bp',type=positive,default=5000000)
    p.add_argument('--halo',type=int,default=20000)
    p.add_argument('--replace',choices=['TRUE','FALSE'],default='FALSE',help='FALSE resumes completed samples/chunks; TRUE explicitly recomputes selected sample outputs')
    p.add_argument('--min-callable',type=fraction,default=.95)
    p.add_argument('--min-region-coverage',type=fraction,default=.80)
    p.add_argument('--stable-retention',type=fraction,default=.95)
    p.add_argument('--stable-min-sample-fraction',type=fraction,default=.90)
    p.add_argument('--stable-min-population-fraction',type=fraction,default=1.0)
    p.add_argument('--bootstrap',type=positive,default=2000)
    p.add_argument('--allow-partial',action='store_true',help='Explicitly produce an interim summary; missing samples remain missing')
    p.add_argument('--target',type=parse_highlight_targets,default='auto',metavar='SOURCE[,SOURCE...]',help='Enrichment sources: gnomad,allofus; auto selects supported sources')
    p.add_argument('--scan-dir',type=Path,help='One completed target store; --target gnomad,allofus instead uses --summary-root')
    p.add_argument('--summary-root',type=Path,default=Path(os.environ.get('ALLELE_SUMMARY_STORE',str(__import__('product_paths').defaults()['summary']))))
    p.add_argument('--herv-input-dir',type=Path,help='Read-only HERV cohort input directory')
    p.add_argument('--herv-bed',type=Path,help='Reference HERV BED for enrichment; defaults to OUT/results/reference_herv.bed, not the stable subset')
    p.add_argument('--enrichment-fold',type=positive_float,default=3.0)
    p.add_argument('--enrichment-windows',type=enrichment_windows,default='10000,100000,1000000')
    p.add_argument('--exclude-bed',type=Path,action='append',default=[],help='Optional masks passed to the conditional enrichment analysis')
    p.add_argument('--no-context',action='store_true',help='Explicit enrichment analysis without trinucleotide context adjustment')
    p.add_argument('--threshold-kind',choices=('auto','AC','AF'),default='auto',help='Allele qualification; AF permits a shared AF-threshold sensitivity analysis across target types')
    p.add_argument('--min-ac',type=positive,default=5)
    p.add_argument('--min-an',type=positive,default=100)
    p.add_argument('--min-af',type=float,default=.0001)
    return p


def main(argv=None):
    args=parser().parse_args(argv)
    from herv_archive import try_reuse
    reused=try_reuse(args)
    if reused is not None:return reused
    args.out=args.out.expanduser().resolve()
    if args.command=='highlights' and args.herv_input_dir:
        source=args.herv_input_dir.expanduser().resolve()
        if args.out==source or source in args.out.parents or args.out in source.parents:
            raise ValueError('Highlights output must be separate from the HERV input directory')
    ensure_layout(args.out)
    if args.out == args.archive.expanduser().resolve() or args.archive.expanduser().resolve() in args.out.parents or args.out in args.archive.expanduser().resolve().parents:
        raise ValueError('HERV runtime and permanent archive must be separate directories')
    for kind in ('logs','jobs','tmp','qc'):
        (args.out/kind).mkdir(parents=True,exist_ok=True)
    for name in ('sample_panel','pfile_dir','fasta','library','reference_bed','repeatmasker','pgen_python',
                 'scan_dir','summary_root','herv_bed','herv_input_dir'):
        value=getattr(args,name)
        if value is not None:
            setattr(args,name,value.expanduser().resolve())
    args.exclude_bed=[path.expanduser().resolve() for path in args.exclude_bed]
    if not math.isfinite(args.min_af) or not 0<=args.min_af<=1:
        raise ValueError('--min-af must be between zero and one')
    if args.command=='export':
        from herv_export import export_cohort
        from herv_publish import publish
        export_cohort(args.out, io_workers=args.cores)
        return publish(args.out,args.archive)
    if args.command=='highlights':
        return highlights(args)
    if args.halo<0:
        raise ValueError('--halo must be nonnegative')
    if args.jobs>args.sample_cores:
        raise ValueError('--jobs must not exceed --sample-cores; CPU affinity caps the full sample process tree')
    if args.command=='summarize':
        if not (args.out/'qc/plan.json').is_file():
            raise ValueError('Run plan or run before summarize')
        result=read_json(args.out/'qc/plan.json')
        # A summary does not redraw the panel or require the original PGENs.
        args.seed=int(result['identity']['seed'])
        return summarize(args,result)
    result=plan(args)
    if args.command=='plan':
        return result
    execute_samples(args,result,single=args.command=='sample')
    return summarize(args,result) if args.command=='run' else dict(status='sample_complete')


def cli(argv=None):
    result = main(argv)
    if result.get('schema') == 'herv-enrichment-highlights-1' and (
            result.get('status') != 'complete' or result.get('comparison', {}).get('status') == 'failed'):
        raise RuntimeError('Highlights did not complete; see enrichment_status.json and retry highlights after its inputs are ready.')
    return result


if __name__=='__main__':
    try:
        cli()
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
        print('herv-scan: '+str(exc),file=sys.stderr)
        raise SystemExit(2)


