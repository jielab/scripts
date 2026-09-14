#!/usr/bin/env python3
"""Phased WGS -> SINGER posterior windows -> chromosome-coordinate ARG samples.

SINGER is an ARG backend, not an independent introgression caller. Keep each
retained posterior draw separate; never treat windows as posterior replicates.
"""
from __future__ import annotations
import argparse, csv, fcntl, hashlib, json, math, os, random, re, runpy, shutil, subprocess, sys, time
from itertools import groupby
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading, signal
import gzip, io, tempfile
from contextlib import contextmanager
from pathlib import Path
import numpy as np
import pysam
import tskit

HERE = Path(__file__).resolve().parent

def stamp(p):
    p=Path(p).resolve(strict=True); s=p.stat()
    return [str(p),s.st_size,s.st_mtime_ns]

def save(p,data):
    p=Path(p); p.parent.mkdir(parents=True,exist_ok=True)
    q=p.with_name(p.name+'.next'); q.write_text(json.dumps(data,indent=2)+'\n');q.replace(p)

def valid(p,request=None):
    try:
        d=json.loads(Path(p).read_text())
        return (request is None or d['request']==request) and all(stamp(x[0])==x for x in d['outputs'])
    except (OSError,ValueError,KeyError): return False

def record(p,request,files):save(p,dict(request=request,outputs=[stamp(x) for x in files]))

def run(cmd,log,stop=None):
    import shlex
    print('[SINGER] '+shlex.join(map(str,cmd)),flush=True)
    with Path(log).open('a') as f:
        f.write('\n'+shlex.join(map(str,cmd))+'\n');f.flush()
        if stop is not None and stop.is_set():raise RuntimeError('SINGER cancelled')
        proc=subprocess.Popen(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT)
        try:
            while True:
                if stop is not None and stop.is_set():raise RuntimeError('SINGER cancelled')
                try: rc=proc.wait(timeout=1);break
                except subprocess.TimeoutExpired:
                    pass  # Keep cancellation responsive without periodic console output.
            if rc:raise RuntimeError(f'command failed ({rc}); see {log}')
        finally:
            if proc.poll() is None:
                proc.terminate()
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:proc.kill();proc.wait()

def parse():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--method',choices=['singer'],default='singer')
    p.add_argument('--action',choices=['all','build','prepare','check'],default='build')
    p.add_argument('--format',choices=['native','trace'],default='native')
    p.add_argument('--dir-gen',default='/mnt/i/refGen/1kg/37')
    p.add_argument('--dir-vcf','--target-vcf-dir',dest='dir_vcf')
    p.add_argument('--arg-dir');p.add_argument('--gen4arg-dir');p.add_argument('--sample-file','--sample-panel',dest='sample_file')
    p.add_argument('--keep');p.add_argument('--chr','--chrs',dest='chr',default='1-22,X')
    p.add_argument('--grch',default='37');p.add_argument('--replace',type=str.upper,choices=['TRUE','FALSE'],default='FALSE')
    p.add_argument('--threads',type=int,default=1,help='total CPU budget; each SINGER window uses one CPU')
    p.add_argument('--jobs',type=int,default=1,help='maximum simultaneous SINGER windows, also limited by CPU and memory budgets')
    p.add_argument('--singer-worker-memory-gb',type=float,default=24,help='reserved GiB per active window including conversion (default: 24 for 5 Mb / ~1000 haplotypes); raise for larger inputs')
    p.add_argument('--features',type=str.upper,choices=['FALSE'],default='FALSE')
    p.add_argument('--sample-match-min-work',type=int,default=50000000)
    p.add_argument('--singer-populations',default='GBR,CHB,ITU,YRI,LWK')
    p.add_argument('--singer-per-pop',type=int,default=100)
    p.add_argument('--singer-samples',type=int,default=50)
    p.add_argument('--singer-burnin',type=int,default=50,help='discard this many saved draws, each after --singer-thin iterations')
    p.add_argument('--singer-thin',type=int,default=20)
    p.add_argument('--singer-window-bp',type=int,default=5000000)
    p.add_argument('--singer-min-sites',type=int,default=20)
    p.add_argument('--singer-ne',type=float,default=10000)
    p.add_argument('--mutation-rate',type=float,default=1.25e-8)
    p.add_argument('--singer-recombination-rate',type=float,default=1.25e-8)
    p.add_argument('--singer-start',type=int,default=0);p.add_argument('--singer-end',type=int)
    p.add_argument('--seed',type=int,default=20260909)
    p.add_argument('--mac-min',type=int,default=1)
    p.add_argument('--missing-max',type=float,choices=[0.0],default=0)
    a=p.parse_args()
    if os.environ.get('REFGEN_GEN4ARG_ONLY')=='1':a.action='prepare'
    for k in ['threads','jobs','singer_samples','singer_thin','singer_window_bp','singer_min_sites','mac_min']:
        if getattr(a,k)<1:p.error(k+' must be positive')
    for k in ['singer_ne','mutation_rate','singer_recombination_rate','singer_worker_memory_gb']:
        if not math.isfinite(getattr(a,k)) or getattr(a,k)<=0:p.error(k+' must be finite and positive')
    if min(a.singer_burnin,a.singer_per_pop,a.singer_start,a.seed)<0:p.error('counts, coordinates and seed must be nonnegative')
    if a.singer_end is not None and a.singer_end<=a.singer_start:p.error('end must exceed start')
    a.grch={'b37':'37','grch37':'37','hg19':'37','b38':'38','grch38':'38','hg38':'38'}.get(a.grch.lower(),a.grch)
    if a.grch not in ('37','38'):p.error('invalid genome build')
    chrs=[]
    for s in re.split('[,; ]+',a.chr):
        s=re.sub('^chr','',s,flags=re.I).upper()
        if not s:continue
        chrs.extend(map(str,range(int(s.split('-')[0]),int(s.split('-')[1])+1))) if '-' in s else chrs.append('X' if s=='23' else s)
    a.chrs=list(dict.fromkeys(chrs))
    if not a.chrs or not set(a.chrs)<=set(map(str,range(1,23)))|{'X'}:p.error('invalid chromosome selection')
    root=Path(a.dir_gen).resolve()
    if root.name in ('vcf','pfile','hap'):root=root.parent
    a.root=root;a.vcf=Path(a.dir_vcf or root/'vcf').resolve()
    a.prep=Path(a.gen4arg_dir or str(a.vcf)+'.4arg.singer').resolve()
    a.output=Path(a.arg_dir or root/'arg.singer').resolve()
    a.panel=Path(a.sample_file or root/'samples.txt').resolve()
    if len({a.vcf,a.prep,a.output})!=3:p.error('source, prepared and ARG directories must differ')
    return a

def sample_panel(a,vcf,chrom):
    available=list(vcf.header.samples)
    meta={}
    if a.panel.is_file():
        with a.panel.open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                name=r.get('sample',r.get('IID',r.get('#IID','')))
                if name:meta[name]=r
    if a.keep:
        wanted=[s.strip().split()[0] for s in Path(a.keep).read_text().splitlines() if s.strip()]
        if len(set(wanted))!=len(wanted):raise ValueError('duplicate --keep sample')
        if set(wanted)-set(available):raise ValueError('--keep contains samples absent from VCF')
        keep=set(wanted)
    else:
        pops=a.singer_populations.split(',');groups={}
        for name in available:
            pop=meta.get(name,{}).get('pop','UNKNOWN')
            if pops==['ALL'] or pop in pops:groups.setdefault(pop,[]).append(name)
        keep=set()
        for pop,names in sorted(groups.items()):
            # Hash ranking is stable across chromosomes and independent of VCF order.
            names.sort(key=lambda n:hashlib.sha256(f'{a.seed}:{n}'.encode()).hexdigest())
            keep.update(names[:a.singer_per_pop] if a.singer_per_pop else names)
    selected=[s for s in available if s in keep]
    if not selected:raise ValueError('no selected samples; provide --sample-file / --keep / --singer-populations ALL')
    if chrom=='X':
        bad=[s for s in selected if str(meta.get(s,{}).get('sex','')).lower() not in ('1','m','male')]
        if bad:raise ValueError('X requires male-only VCF/panel; nonmale or unknown sex: '+','.join(bad[:5]))
    if len(selected)*(1 if chrom=='X' else 2)<4:raise ValueError('SINGER requires at least four selected haplotypes in this workflow')
    return selected,meta

def source_vcf(a,chrom):
    return a.vcf/('chrX.male.vcf.gz' if chrom=='X' else f'chr{chrom}.vcf.gz')

def unique_snps(records,start,end,dropped):
    """Exclude every PASS biallelic SNP at a repeated position, before GT filters."""
    def candidates():
        for r in records:
            # Indexed fetch can include an earlier structural variant overlapping this window.
            if not start < r.pos <= end:continue
            reason=None
            if set(r.filter) not in (set(),{'PASS'}):reason='filter'
            elif len(r.alleles)!=2 or any(x not in ('A','C','G','T') for x in r.alleles):reason='not_biallelic_snp'
            if reason:dropped[reason]=dropped.get(reason,0)+1
            else:yield r
    for _,group in groupby(candidates(),key=lambda r:r.pos):
        first=next(group);count=1+sum(1 for _ in group)
        if count>1:
            dropped['duplicate_position']=dropped.get('duplicate_position',0)+count
        else:yield first

def prep_request(a,chrom):
    src=source_vcf(a,chrom)
    return dict(source=stamp(src),index=stamp(str(src)+('.tbi' if Path(str(src)+'.tbi').is_file() else '.csi')),
        code=stamp(__file__),panel=stamp(a.panel) if a.panel.is_file() else None,
        keep=stamp(a.keep) if a.keep else None,populations=a.singer_populations,per_pop=a.singer_per_pop,
        seed=a.seed,window=a.singer_window_bp,start=a.singer_start,end=a.singer_end,
        min_sites=a.singer_min_sites,mac=a.mac_min,build=a.grch)

@contextmanager
def compressed_vcf_output(path):
    temporary = path.with_name(path.name+'.next')
    try:
        with io.TextIOWrapper(pysam.BGZFile(str(temporary), 'w')) as out:
            yield out
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def prepared_vcf(prefix):
    compressed = Path(str(prefix)+'.vcf.gz')
    return compressed if compressed.exists() else Path(str(prefix)+'.vcf')


def remove_native_text(wd):
    for path in wd.glob('arg_*.txt*'):
        path.unlink(missing_ok=True)


def prepare(a,chrom):
    dest=a.prep/f'chr{chrom}';done=dest/'complete.json';req=prep_request(a,chrom)
    if a.replace=='FALSE' and valid(done,req):return json.loads((dest/'windows.json').read_text())
    done.unlink(missing_ok=True);dest.mkdir(parents=True,exist_ok=True)
    v=pysam.VariantFile(str(source_vcf(a,chrom)))
    matches=[c for c in v.header.contigs if c.removeprefix('chr') in (chrom,'23' if chrom=='X' else chrom)]
    if len(matches)!=1:raise ValueError(f'expected one chr{chrom} contig')
    contig=matches[0];length=v.header.contigs[contig].length
    if not length:raise ValueError('VCF contig must declare chromosome length')
    samples,meta=sample_panel(a,v,chrom);v.subset_samples(samples)
    lo=a.singer_start;hi=min(a.singer_end or length,length)
    if chrom=='X':lo=max(lo,2699520 if a.grch=='37' else 2781479);hi=min(hi,154931043 if a.grch=='37' else 155701382)
    if hi<=lo:raise ValueError('empty selected interval')
    ploidy=1 if chrom=='X' else 2;windows=[];files=[];qc=[]
    for start in range(lo,hi,a.singer_window_bp):
        end=min(start+a.singer_window_bp,hi);prefix=dest/f'w{start}_{end}';path=Path(str(prefix)+'.vcf.gz');count=0;dropped={}
        with compressed_vcf_output(path) as out:
            out.write('##fileformat=VCFv4.2\n'+f'##contig=<ID={contig},length={end-start}>\n')
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n')
            out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n')
            for r in unique_snps(v.fetch(contig,start,end),start,end,dropped):
                reason=None
                gts=[];flat=[]
                for name in samples:
                    call=r.samples[name];gt=call.get('GT')
                    if gt is None or any(x is None for x in gt):reason='missing';break
                    if len(gt)!=ploidy:raise ValueError(f'wrong ploidy {name} at {r.pos}')
                    if ploidy==2 and len(set(gt))>1 and not call.phased:raise ValueError(f'unphased heterozygote {name} at {r.pos}')
                    gts.append('|'.join(map(str,gt)));flat.extend(gt)
                if not reason and min(sum(flat),len(flat)-sum(flat))<a.mac_min:reason='mac'
                if reason:dropped[reason]=dropped.get(reason,0)+1;continue
                out.write(f'{contig}\t{r.pos-start}\t.\t{r.ref}\t{r.alts[0]}\t.\tPASS\t.\tGT\t'+'\t'.join(gts)+'\n');count+=1
        qc.append(dict(start=start,end=end,sites=count,dropped=dropped))
        print(f'[SINGER] chr{chrom} window {start}-{end}: sites={count} dropped={dropped}',flush=True)
        if count<a.singer_min_sites:path.unlink();continue
        pysam.tabix_index(str(path),preset='vcf',force=True)
        files.extend([path,Path(str(path)+'.tbi')]);windows.append(dict(start=start,end=end,sites=count,prefix=str(prefix)))
    v.close()
    if not windows:raise ValueError('no windows with sufficient polymorphic sites')
    info=dict(chrom=chrom,length=length,ploidy=ploidy,samples=samples,
        populations={s:meta.get(s,{}).get('pop','UNKNOWN') for s in samples},windows=windows,qc=qc,
        regional=a.singer_start!=0 or a.singer_end is not None)
    save(dest/'windows.json',info);record(done,req,files+[dest/'windows.json'])
    for old in dest.glob('w*.vcf*'):
        if old not in files:old.unlink()
    print(f'[SINGER] prepared chr{chrom}: samples={len(samples)} windows={len(windows)}',flush=True)
    return info

def load_expected(prefix,ploidy):
    """Keep the full expected genotype matrix compact and load it once per window."""
    positions=[];rows=[]
    with pysam.VariantFile(str(prepared_vcf(prefix))) as f:
        n=len(f.header.samples)*ploidy
        for r in f:
            positions.append(r.pos-1)
            rows.append(np.fromiter((x for c in r.samples.values() for x in c['GT']),dtype=np.int8,count=n))
    if not rows:raise ValueError('empty SINGER input')
    return np.asarray(positions,dtype=np.float64),np.stack(rows)

def validate_window(ts,prefix,ploidy,expected=None):
    """Check every leaf and site; explicit alleles make integer comparison exact."""
    positions,genotypes=load_expected(prefix,ploidy) if expected is None else expected
    if ts.num_samples!=genotypes.shape[1]:raise ValueError('SINGER sample count changed')
    if ts.num_sites!=len(positions):raise ValueError('SINGER lost input sites')
    for i,v in enumerate(ts.variants(alleles=('0','1'))):
        if v.site.position!=positions[i]:raise ValueError(f'SINGER coordinate mismatch {v.site.position}')
        if not np.array_equal(v.genotypes,genotypes[i]):
            raise ValueError(f'SINGER leaf genotype mismatch at {positions[i]}')

def merge_windows(chrom,info,sources,dest,validated=False):
    tables=tskit.TableCollection(info['length']);tables.time_units='generations'
    tables.metadata_schema=tskit.MetadataSchema.permissive_json()
    tables.metadata=dict(method='singer',chromosome=chrom,coordinate_system='0-based chromosome bp',regional=info['regional'])
    tables.nodes.metadata_schema=tskit.MetadataSchema.permissive_json()
    tables.individuals.metadata_schema=tskit.MetadataSchema.permissive_json()
    n=len(info['samples'])*info['ploidy']
    for name in info['samples']:
        ind=tables.individuals.add_row(metadata=dict(sample=name,population=info['populations'][name]))
        for h in range(info['ploidy']):tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE,time=0,individual=ind,metadata=dict(sample=name,haplotype=h+1))
    for w,source in zip(info['windows'],sources):
        ts=tskit.load(source)
        if not validated:validate_window(ts,w['prefix'],info['ploidy'])
        if ts.sequence_length!=w['end']-w['start']:raise ValueError('SINGER window length mismatch')
        mapping={int(node):j for j,node in enumerate(ts.samples())}
        for node in ts.nodes():
            if not node.is_sample():mapping[node.id]=tables.nodes.add_row(flags=node.flags,time=node.time)
        for e in ts.edges():tables.edges.add_row(e.left+w['start'],e.right+w['start'],mapping[e.parent],mapping[e.child])
        for site in ts.sites():
            s=tables.sites.add_row(site.position+w['start'],site.ancestral_state)
            for mut in site.mutations:tables.mutations.add_row(s,mapping[mut.node],mut.derived_state)
    tables.sort();tables.build_index();tables.compute_mutation_parents();tables.compute_mutation_times()
    ts=tables.tree_sequence()
    if ts.num_samples!=n or not ts.num_mutations:raise ValueError('invalid merged ARG')
    tmp=dest.with_suffix('.next')
    try:
        ts.dump(tmp);tmp.replace(dest)
    finally:
        tmp.unlink(missing_ok=True)
    return dict(sites=ts.num_sites,trees=ts.num_trees,samples=ts.num_samples,sequence_length=ts.sequence_length)

def memory_budget_gb():
    """Use available host RAM and all visible cgroup ancestor limits."""
    gib=1024**3;limits=[]
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):limits.append(int(line.split()[1])*1024*.85/gib)
    if os.environ.get('ARG_MEMORY_LIMIT_GB'):limits.append(float(os.environ['ARG_MEMORY_LIMIT_GB'])*.85)
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        if line.startswith('0::'):
            root=Path('/sys/fs/cgroup');path=root/line[3:].lstrip('/')
            while path==root or root in path.parents:
                try:
                    cap=(path/'memory.max').read_text().strip()
                    used=int((path/'memory.current').read_text())
                    if cap!='max':limits.append(max(0,int(cap)-used)*.85/gib)
                except (OSError,ValueError):pass
                if path==root:break
                path=path.parent
    return min(limits) if limits else 24.

def window_workers(a,count):
    capacity=memory_budget_gb()
    cpu=len(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else os.cpu_count() or 1
    workers=min(count,a.jobs,a.threads,cpu,max(1,int(capacity/a.singer_worker_memory_gb)))
    if capacity<a.singer_worker_memory_gb:
        print(f'[SINGER] WARNING: available budget {capacity:.1f} GiB is below the per-window reservation; one window may exceed the memory cap',flush=True)
    return workers

def infer_window(a,info,home,root,req,i,w,stop,read_ts):
    if stop.is_set():raise RuntimeError('SINGER cancelled')
    n=a.singer_burnin+a.singer_samples
    wd=root/f"w{w['start']}_{w['end']}";wd.mkdir(exist_ok=True);prefix=wd/'arg';mark=wd/'complete.json'
    request=dict(req,window=w)
    if a.replace=='FALSE' and valid(mark,request):return
    mark.unlink(missing_ok=True)
    for partial in wd.glob('posterior_*.trees'):partial.unlink()
    raw=wd/'raw.complete.json'
    if a.replace=='TRUE' or not valid(raw,request):
        raw.unlink(missing_ok=True)
        remove_native_text(wd)
        try:
            # SINGER accepts a filename prefix, so expand only this active window.
            with tempfile.TemporaryDirectory(prefix='input.', dir=wd) as temporary:
                input_prefix=Path(temporary)/'window'
                source=prepared_vcf(w['prefix'])
                with (gzip.open(source,'rb') if source.suffix=='.gz' else source.open('rb')) as inp, Path(str(input_prefix)+'.vcf').open('wb') as out:
                    shutil.copyfileobj(inp,out)
                run([home/'singer','-Ne',a.singer_ne,'-m',a.mutation_rate,'-r',a.singer_recombination_rate,
                 '-input',input_prefix,'-output',prefix,'-start',0,'-end',w['end']-w['start']+1,
                 '-polar',.5,'-n',n,'-thin',a.singer_thin,'-ploidy',info['ploidy'],'-seed',a.seed+i],wd/'run.log',stop)
            raw_files=[]
            for j in range(n):
                for kind in ('nodes','branches','muts'):
                    path=Path(str(prefix)+f'_{kind}_{j}.txt')
                    if j<a.singer_burnin:
                        path.unlink(missing_ok=True);continue
                    packed=Path(str(path)+'.gz')
                    with path.open('rb') as inp,gzip.open(packed,'wb',compresslevel=1) as out:
                        shutil.copyfileobj(inp,out)
                    path.unlink();raw_files.append(packed)
            for unused in wd.glob('arg_*.txt'):unused.unlink()
            record(raw,request,raw_files)
        except BaseException:
            remove_native_text(wd)
            raise
    # Compressed native draws allow conversion to resume after a failure.
    try:
        files=[];expected=load_expected(w['prefix'],info['ploidy'])
        for j in range(a.singer_burnin,n):
            if stop.is_set():raise RuntimeError('SINGER cancelled')
            with tempfile.TemporaryDirectory(prefix='convert.',dir=wd) as temporary:
                paths={}
                for kind in ('nodes','branches','muts'):
                    paths[kind]=Path(temporary)/kind
                    source=Path(str(prefix)+f'_{kind}_{j}.txt.gz')
                    with gzip.open(source,'rb') as inp,paths[kind].open('wb') as out:
                        shutil.copyfileobj(inp,out)
                tables=read_ts(str(paths['nodes']),str(paths['branches']))
                mutations=np.atleast_2d(np.loadtxt(paths['muts']))
            previous=None
            for row in mutations:
                if row[0]!=previous:
                    site=tables.sites.add_row(position=row[0],ancestral_state='0');previous=row[0]
                tables.mutations.add_row(site=site,node=int(row[1]),derived_state=str(int(row[3])))
            tables.sort();tables.build_index();tables.compute_mutation_parents();tables.compute_mutation_times()
            ts=tables.tree_sequence().keep_intervals([[1,w['end']-w['start']+1]],simplify=False).trim()
            validate_window(ts,w['prefix'],info['ploidy'],expected)
            file=wd/f'posterior_{j}.trees';ts.dump(file);files.append(file)
        (wd/'convert.log').write_text(f'Validated {len(files)} retained posterior draws against every input genotype.\n')
        record(mark,request,files+[wd/'run.log',wd/'convert.log'])
        remove_native_text(wd);raw.unlink(missing_ok=True)
    except BaseException:
        for partial in wd.glob('posterior_*.trees'):partial.unlink(missing_ok=True)
        raise

def infer(a,chrom,info,home):
    root=a.output/'native'/f'chr{chrom}';root.mkdir(parents=True,exist_ok=True)
    done=root/'complete.json';req=dict(prepared=stamp(a.prep/f'chr{chrom}'/'complete.json'),code=stamp(__file__),
        binary=stamp(home/'singer'),converter=stamp(home/'convert_to_tskit'),samples=a.singer_samples,burnin=a.singer_burnin,
        thin=a.singer_thin,Ne=a.singer_ne,mutation=a.mutation_rate,recombination=a.singer_recombination_rate,seed=a.seed)
    if a.replace=='FALSE' and valid(done,req):print(f'[SINGER] reuse chr{chrom}',flush=True);return
    done.unlink(missing_ok=True);n=a.singer_burnin+a.singer_samples
    workers=window_workers(a,len(info['windows']))
    print(f'[SINGER] chr{chrom}: {workers} concurrent windows; CPU budget={a.threads}, requested jobs={a.jobs}, reserve={a.singer_worker_memory_gb:g} GiB/window',flush=True)
    stop=threading.Event()
    read_ts=runpy.run_path(str(home/'convert_to_tskit'))['read_ts']
    pool=ThreadPoolExecutor(max_workers=workers)
    futures=[]
    try:
        for i,w in enumerate(info['windows']):
            futures.append(pool.submit(infer_window,a,info,home,root,req,i,w,stop,read_ts))
        for completed,future in enumerate(as_completed(futures),1):
            future.result()
            print(f'[SINGER] chr{chrom}: windows complete {completed}/{len(futures)}',flush=True)
    except BaseException:
        stop.set()
        for future in futures:future.cancel()
        raise
    finally:
        pool.shutdown(wait=True,cancel_futures=True)
    # Only merge outputs whose complete genotype validation and file stamps pass.
    for w in info['windows']:
        mark=root/f"w{w['start']}_{w['end']}"/'complete.json'
        if not valid(mark,dict(req,window=w)):raise ValueError(f'changed window outputs: {mark}')
    merged=a.output/'posterior';merged.mkdir(exist_ok=True);files=[];stats=[]
    # Remove only stale posterior exports for this chromosome on explicit rebuild.
    for p in merged.glob(f'chr{chrom}.sample*.trees'):p.unlink()
    for j in range(a.singer_burnin,n):
        sources=[root/f"w{w['start']}_{w['end']}"/f'posterior_{j}.trees' for w in info['windows']]
        dest=merged/f'chr{chrom}.sample{j-a.singer_burnin:03d}.trees'
        stats.append(merge_windows(chrom,info,sources,dest,validated=True));files.append(dest)
    smap=merged/f'chr{chrom}.sample_map.tsv'
    with smap.open('w') as f:
        f.write('tree_node_id\tsample\thaplotype\n')
        for i,s in enumerate(info['samples']):
            for h in range(info['ploidy']):f.write(f'{i*info["ploidy"]+h}\t{s}\t{h+1}\n')
    mask=merged/f'chr{chrom}.inferred.bed';mask.write_text(''.join(f'chr{chrom}\t{w["start"]}\t{w["end"]}\n' for w in info['windows']))
    save(root/'qc.json',dict(stats=stats,regional=info['regional'],retained=a.singer_samples,burnin=a.singer_burnin,
        thin=a.singer_thin,convergence='not assessed: inspect MCMC logs before scientific use'))
    record(done,req,files+[smap,mask,root/'qc.json'])
    # The validated chromosome posterior is now the durable scientific output.
    for w in info['windows']:
        wd=root/f"w{w['start']}_{w['end']}"
        for path in wd.glob('posterior_*.trees'):path.unlink()
        remove_native_text(wd)
        (wd/'raw.complete.json').unlink(missing_ok=True)

def export(a,chrom):
    dest=a.output/'trace/singer';dest.mkdir(parents=True,exist_ok=True)
    for p in dest.glob(f'chr{chrom}.sample*.trees'):p.unlink()
    files=list((a.output/'posterior').glob(f'chr{chrom}.*'))
    for p in files:
        q=dest/p.name
        if q.is_symlink() or q.exists():q.unlink()
        q.symlink_to(p.resolve())
    marker=dest/'ARG_TRACE_BUILD.tsv'
    text='method\tsinger\nformat\ttrace\ncoordinate_system\t0-based chromosome bp\ntime_units\tgenerations\n'
    if not marker.exists() or marker.read_text()!=text:marker.write_text(text)
    record(dest/f'chr{chrom}.complete.json',dict(native=stamp(a.output/'native'/f'chr{chrom}'/'complete.json')), [dest/p.name for p in files]+[dest/'ARG_TRACE_BUILD.tsv'])

def check(a,chrom):
    p=a.output/'native'/f'chr{chrom}'/'complete.json'
    if not valid(p):raise ValueError(f'incomplete or changed SINGER outputs: {p}')
    d=json.loads(p.read_text());req=d['request']
    for key in ['prepared','code','binary','converter']:
        if stamp(req[key][0])!=req[key]:raise ValueError(f'SINGER {key} changed; rebuild required')
    prep=Path(req['prepared'][0]);data=json.loads(prep.read_text())
    if not valid(prep):raise ValueError('prepared input changed')
    for key in ['source','index','panel','keep']:
        row=data['request'].get(key)
        if row and stamp(row[0])!=row:raise ValueError(f'source {key} changed')
    tsfiles=[Path(x[0]) for x in d['outputs'] if x[0].endswith('.trees')]
    if len(tsfiles)!=req['samples']:raise ValueError('posterior count mismatch')
    for f in tsfiles:
        ts=tskit.load(f)
        if ts.time_units!='generations' or not ts.num_sites or not ts.num_mutations:raise ValueError('invalid tree')
    if a.format=='trace' and not valid(a.output/'trace/singer'/f'chr{chrom}.complete.json'):raise ValueError('TRACE export missing or changed')
    print(f'[SINGER] CHECK PASS chr{chrom}: {len(tsfiles)} posterior samples',flush=True)

def main():
    def interrupted(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    for sig in (signal.SIGTERM,signal.SIGHUP):signal.signal(sig,interrupted)
    a=parse();home=Path(os.environ.get('GU_SINGER_HOME',str(Path.home()/'software/singer-0.1.9')))
    if a.action=='check':
        for chrom in a.chrs:check(a,chrom)
        return
    if not (home/'singer').is_file():raise ValueError('SINGER missing; run bash f/install_singer.sh')
    locks=[]
    for path in (a.prep,a.output):
        path.mkdir(parents=True,exist_ok=True);lock=(path/'.refgen.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);locks.append(lock)
    for chrom in a.chrs:
        info=prepare(a,chrom)
        if a.action=='prepare':continue
        infer(a,chrom,info,home)
        if a.format=='trace':export(a,chrom)
        check(a,chrom)
    print('[SINGER] COMPLETE',flush=True)

if __name__=='__main__':
    try:main()
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as e:
        print('ERROR: '+str(e),file=sys.stderr);sys.exit(2)
