"""Read flat GU chromosome exports and build one reusable GRCh38 track per cohort."""
import argparse
import csv
import fcntl
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zlib
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

from research_common import Intervals, atomic_json, chrom_key, identity, write_tsv

COHORTS = ('1kg', *(f'ygb{i:03d}' for i in range(1, 50)))
GU_CHROMS = tuple('chr'+str(i) for i in range(1,23)) + ('chrX',)
READ_ERRORS = (OSError, EOFError, ValueError, csv.Error, zlib.error)


class NoGuSegments(ValueError):
    """The optional GU comparison has no usable chromosomes."""


def lift_program():
    return shutil.which(os.environ.get('ALLELE_ARCHAIC_LIFTOVER') or 'liftOver') or shutil.which('liftOver') or ''


def cohort_for(raw):
    raw = Path(raw)
    if raw.name != 'raw' or raw.parent.name not in COHORTS:
        raise ValueError('Expected GU <1kg|ygb001..ygb049>/raw: ' + str(raw))
    return raw.parent.name


def chromosome_unit(raw, chrom):
    root = raw.parent
    done = root/'jobs'/(chrom+'.complete')
    meta = root/'jobs'/(chrom+'.run.meta.tsv')
    if not done.is_file() or not meta.is_file():
        raise ValueError('Missing GU completion marker or run metadata for '+chrom)
    cfg = dict(line.split('\t',1) for line in meta.read_text().splitlines() if '\t' in line)
    refs = cfg.get('refs','').split()
    build = cfg.get('genome_build','').removeprefix('GRCh').removeprefix('b')
    if build not in ('37','38') or not refs:
        raise ValueError('Missing reference list or unsupported genome build: '+str(meta))
    descriptor = root/'metadata'/(chrom+'.json')
    if descriptor.is_file():
        info = json.loads(descriptor.read_text())
        other = str(info.get('genome_build', build)).removeprefix('GRCh').removeprefix('b')
        if (other != build or info.get('native_coordinate_system','1-based-exclusive-end') != '1-based-exclusive-end'):
            raise ValueError('Inconsistent native coordinate metadata: '+str(descriptor))
    files = []
    populations = defaultdict(set)
    for path in sorted(raw.glob(chrom+'.*.txt.gz')):
        stem = path.name[len(chrom)+1:-len('.txt.gz')]
        ref = next((r for r in sorted(refs,key=len,reverse=True) if stem==r or stem.startswith(r+'.')),None)
        if ref is None:
            raise ValueError('Unexpected reference: '+path.name)
        pop = stem[len(ref):].lstrip('.') or 'YGB'
        populations[ref].add(pop)
        files.append((path,ref,pop))
    if set(populations)!=set(refs) or len({frozenset(p) for p in populations.values()}) != 1:
        raise ValueError('Incomplete reference/population raw files for '+chrom)
    groups = root/'tmp'/(chrom+'.samples.'+('' if chrom=='chrX' else 'C'+chrom[3:]+'.')+'populations.tsv')
    if groups.is_file():
        with groups.open() as stream:
            expected={row['population'] for row in csv.DictReader(stream,delimiter='\t')}
        if expected != next(iter(populations.values())):
            raise ValueError('Raw files disagree with GU population list: '+str(groups))
    return build, files


def inventory(raw, chroms, chain, lift):
    raw = Path(raw)
    units, skipped = {}, {}
    for chrom in chroms:
        try:
            build, files = chromosome_unit(raw, chrom)
            if build=='37' and (not chain or not Path(chain).is_file() or not lift or not shutil.which(lift)):
                raise ValueError('GRCh37 calls require hg19ToHg38 chain and liftOver')
            units[chrom] = build, files
        except READ_ERRORS + (KeyError,) as exc:
            skipped[chrom] = str(exc)
    return units, skipped


def load_segments(files, chrom, build, options, temp):
    # Use the same native (start-1,end-1) conversion and LOD/length filters as GU.
    # Import only on workers; metadata-only planners need no scientific packages.
    gu = str(Path(__file__).resolve().parents[2]/'gu/f')
    if gu not in sys.path: sys.path.insert(0, gu)
    from ibdmix_raw import native_chunks
    unit = SimpleNamespace(name=chrom,config={'lod_cut':options.archaic_min_slod,'len_cut':options.archaic_min_bp})
    people, events = defaultdict(list), defaultdict(lambda:defaultdict(int))
    selected = [(p,r,pop) for p,r,pop in files if options.archaic_reference in ('all',r)]
    if not selected:
        raise ValueError('Requested archaic reference is absent')
    src, dst, reject = (temp/n for n in ('source.bed','mapped.bed','rejected.bed'))
    ids = {}
    count = unmapped = 0
    with src.open('w') as stream:
        for path,ref,pop in selected:
            with gzip.open(path,'rt') as gz:
                if not {'ID','chrom','start','end','slod'} <= set(next(csv.reader(gz,delimiter='\t'),())):
                    raise ValueError('Invalid native header: '+str(path))
            for frame in native_chunks(path,unit):
                for r in frame.itertuples(index=False):
                    key = ids.setdefault((ref,str(r.ID)),len(ids))
                    stream.write(f'{chrom}\t{r.start}\t{r.end}\t{key}\n')
                    count += 1
    mapped = src
    if count and build=='37':
        result = subprocess.run([options.liftover,str(src),str(options.archaic_chain),str(dst),str(reject)],capture_output=True,text=True)
        if result.returncode:
            raise RuntimeError('liftOver failed: '+result.stderr.strip())
        mapped = dst
        with reject.open() as stream:
            unmapped=sum(1 for line in stream if line.strip() and not line.startswith('#'))
    keys = list(ids)
    with mapped.open() as stream:
        for target,start,end,index in csv.reader(stream,delimiter='\t'):
            lo,hi,i=int(start),int(end),int(index)
            if not (0<=lo<hi and 0<=i<len(keys)):
                raise ValueError('Invalid projected interval for '+chrom)
            if target==chrom:people[keys[i]].append((lo,hi))
            else:unmapped+=1
    if count and not people:
        raise ValueError('No retained segment mapped to '+chrom)
    for (ref,_),intervals in people.items():
        for _,lo,hi,_ in Intervals((chrom,x,y,ref) for x,y in intervals).rows():
            events[ref][lo]+=1
            events[ref][hi]-=1
    rows,coverage=[],[]
    for ref in sorted({r for _,r,_ in selected}):
        active=bp=0
        last=None
        for pos,change in sorted(events[ref].items()):
            if last is not None and active>=options.archaic_min_carriers:
                rows.append((chrom,last,pos,'archaic:'+ref))
                bp+=pos-last
            active+=change
            last=pos
        coverage.append((chrom,ref,bp,options.archaic_min_carriers))
    return list(Intervals(rows).rows()),coverage,count,unmapped,selected


def build_track(raw, out, chroms, options):
    raw,out=Path(raw),Path(out)
    cohort=cohort_for(raw)
    out.mkdir(parents=True,exist_ok=True)
    units,skipped=inventory(raw,chroms,options.archaic_chain,options.liftover)
    coverage,rows,surveyed,files=[],[],{},[]
    n=unmapped=0
    with tempfile.TemporaryDirectory(dir=out) as tmp:
        for chrom,(build,parts) in units.items():
            try:
                tracts,cov,count,rejected,selected=load_segments(parts,chrom,build,options,Path(tmp))
            except READ_ERRORS as exc:
                skipped[chrom]=str(exc)
                continue
            surveyed[chrom]=build
            coverage.extend(cov)
            rows.extend(tracts)
            files.extend(identity(p) for p,_,_ in selected)
            n+=count
            unmapped+=rejected
        index=Intervals(rows)
        path=out/'annotations/archaic.bed'
        write_tsv(Path(tmp)/'track.bed',['#chrom','start0','end0','label'],index.rows())
        write_tsv(Path(tmp)/'coverage.tsv',['chrom','reference','covered_bp','min_distinct_carriers'],coverage)
        path.parent.mkdir(parents=True,exist_ok=True)
        os.replace(Path(tmp)/'track.bed',path)
        os.replace(Path(tmp)/'coverage.tsv',path.with_name('archaic_coverage.tsv'))
    meta=dict(status='complete' if surveyed else 'skipped',gu_cohort=cohort,raw=str(raw),files=files,
        input_rows=n,source_builds=surveyed,source_populations=sorted({pop for c in surveyed for _,_,pop in units[c][1]}),
        references=sorted({ref for _,ref,_,_ in coverage}),surveyed_chromosomes=sorted(surveyed,key=chrom_key),
        skipped_chromosomes=skipped,unmapped_grch37_segments=unmapped,coordinate='GRCh38 BED0',
        min_carriers=options.archaic_min_carriers,min_slod=options.archaic_min_slod,min_bp=options.archaic_min_bp,
        selected_reference=options.archaic_reference,output=str(path),
        interpretation='Overlap with unphased GU tracts; the carrier threshold is applied within this cohort and reference.')
    if not surveyed:meta['reason']='No usable GU chromosomes'
    atomic_json(out/'track.json',meta)
    return index,meta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--chain',dest='archaic_chain',required=True)
    p.add_argument('--liftover',default=lift_program())
    args=p.parse_args()
    args.archaic_min_carriers,args.archaic_min_slod,args.archaic_min_bp,args.archaic_reference=2,4,50000,'all'
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/'track.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        atomic_json(args.out/'track.json',dict(status='running',gu_cohort=cohort_for(args.raw)))
        _,meta=build_track(args.raw,args.out,GU_CHROMS,args)
        print(f"[GU] {meta['gu_cohort']}: {meta['status']}, {len(meta['surveyed_chromosomes'])} chromosomes, {meta['input_rows']} retained raw segments")

if __name__=='__main__':main()
