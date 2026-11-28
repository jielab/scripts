#!/usr/bin/env python3
"""Streaming multi-scale AVI landscape, and exact substitution extraction.

The C++/zlib worker reads each chromosome once, using bounded memory. These are
predicted substitutions, not observed population alleles. SHAP is not a raw
molecular-effect vector. No Python row object is created for every genome base.
"""
from __future__ import annotations
from product_manifest import finish
import argparse, array, csv, fcntl, gzip, html, json, os, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from output_layout import data_path

import sys as _path_sys
_path_sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "0f"))
from workflow_paths import configured_path
from research_common import (atomic_json, identity, chromosomes, chrom_name,
    chrom_key, open_text, number, write_tsv, html_report, html_table, Intervals)

HERE=Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent/"0f"))
from workflow_cpu import core_budget, worker_cores
METRICS=('mean_avi','sd_site_mean','mean_alt_range','mean_abs_gradient','high_fraction','low_fraction')


def native_binary(cache,explicit=None):
    if explicit:
        p=Path(explicit).resolve()
        if not p.is_file() or not os.access(p,os.X_OK):raise ValueError('Native binary not executable: '+str(p))
        return p
    cache=Path(cache);cache.mkdir(parents=True,exist_ok=True);binary=cache/'avi_stream'
    source=HERE/'avi_stream.cpp'
    with (cache/'build.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if binary.is_file() and binary.stat().st_mtime_ns>=source.stat().st_mtime_ns:return binary
        # Conda g++ uses an isolated sysroot that may omit system zlib headers.
        # Honor explicit toolchains; otherwise prefer the host compiler and libraries.
        requested=os.environ.get('CXX')
        compiler=shutil.which(requested) if requested else (shutil.which('/usr/bin/g++') or shutil.which('g++'))
        if not compiler:raise ValueError('AVI streaming requires g++/C++17 and zlib headers, or --native-binary')
        tmp=binary.with_suffix('.partial')
        subprocess.run([compiler,'-O3','-std=c++17',str(source),'-lz','-o',str(tmp)],check=True)
        os.replace(tmp,binary)
    return binary


def resolve_vcf(directory,chrom):
    directory=Path(directory);chrom=chrom_name(chrom)
    stems=('chrM','M','chrMT','MT') if chrom=='chrM' else (chrom,chrom[3:])
    for stem in stems:
        for suffix in ('.vcf.bgz','.vcf.gz','.vcf'):
            p=directory/(stem+suffix)
            if p.is_file():return p
    raise ValueError('Atlas VCF missing for '+chrom+' in '+str(directory))


def query_file(source,chrom,targets,out,cache,*,native=None):
    """Targets: sorted or unsorted one-based positions, one per line.

    Output is an ordinary TSV; all available Number=A scores preserve ALT order.
    Reads once per chromosome; no millions of individual tabix invocations.
    """
    binary=native_binary(cache,native);out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    temp=out.with_suffix(out.suffix+'.partial')
    try:
        subprocess.run([str(binary),'--input',str(source),'--chrom',chrom_name(chrom),
            '--targets',str(targets),'--out',str(temp)],check=True)
        os.replace(temp,out)
    finally:
        if temp.exists():temp.unlink()
    return out


def run_chromosome(a,chrom,source,binary):
    out=Path(a.out);shard=out/'shards'/f'{chrom}.windows.tsv.gz';receipt=data_path(out, 'qc')/f'{chrom}.json'
    parameters=dict(windows=a.window_sizes,field=a.score_field,high=a.high,low=a.low)
    current=dict(source=identity(source),parameters=parameters,worker=identity(HERE/'avi_stream.cpp'),native=identity(binary))
    if a.replace=='FALSE' and shard.is_file() and receipt.is_file():
        old=json.loads(receipt.read_text())
        if old.get('status')=='complete' and old.get('output_identity')==identity(shard) and all(old.get(k)==v for k,v in current.items()):
            print('[AVI reuse] '+chrom,flush=True);return old
    log=data_path(out, 'logs')/f'{chrom}.log';raw=data_path(out, 'tmp')/f'{chrom}.windows.tsv';temp=shard.with_suffix('.partial.gz')
    for p in (shard.parent,receipt.parent,raw.parent,log.parent):p.mkdir(parents=True,exist_ok=True)
    cmd=[str(binary),'--input',str(source),'--chrom',chrom,'--out',str(raw),
        '--windows',a.window_sizes,'--field',a.score_field,'--high',str(a.high),'--low',str(a.low)]
    atomic_json(receipt,dict(current,status='running',command=cmd))
    try:
        with log.open('w') as f:subprocess.run(cmd,stderr=f,stdout=f,check=True)
        rows=0
        with raw.open('rb') as src,gzip.open(temp,'wb',compresslevel=3) as dst:
            shutil.copyfileobj(src,dst,4*1024*1024)
        with gzip.open(temp,'rt') as f:rows=sum(1 for _ in f)-1
        if rows<1:raise ValueError('No scored AVI positions for '+chrom)
        os.replace(temp,shard);raw.unlink()
        record=dict(current,status='complete',chrom=chrom,rows=rows,output=str(shard),output_identity=identity(shard))
        atomic_json(receipt,record);print('[AVI complete] '+chrom,flush=True);return record
    except Exception as e:
        atomic_json(receipt,dict(current,status='failed',reason=str(e)));raise
    finally:
        if temp.exists():temp.unlink()


def eligible(row,a):
    return (float(row['coverage_fraction'])>=a.min_coverage and
        int(row['n_positions'])>=a.min_positions and
        int(row['n_complete_positions'])/int(row['n_positions'])>=a.min_complete_fraction)


def merge(a,requested=None):
    import numpy as np
    out=Path(a.out);shards=sorted((out/'shards').glob('chr*.windows.tsv.gz'),key=lambda p:chrom_key(p.name.split('.')[0]))
    if requested is not None:
        wanted=set(requested);shards=[p for p in shards if p.name.split('.')[0] in wanted]
        if len(shards)!=len(wanted):raise ValueError('Not all requested AVI shards completed')
    if not shards:raise ValueError('No AVI shards; run scan first')
    values={};counts={};fields=None;parameters=None;worker=None;native=None
    requested_parameters=dict(windows=a.window_sizes,field=a.score_field,high=a.high,low=a.low)
    (out/'results').mkdir(parents=True,exist_ok=True)
    merged=out/'results/windows.tsv.gz'
    with gzip.open(merged.with_suffix('.partial.gz'),'wt',compresslevel=3) as dest:
        writer=None
        for shard in shards:
            rec=json.loads((data_path(out, 'qc')/(shard.name.split('.')[0]+'.json')).read_text())
            if rec.get('status')!='complete':raise ValueError('Incomplete shard: '+str(shard))
            if rec['parameters']!=requested_parameters:raise ValueError('Requested AVI score/window parameters differ from shard: '+str(shard))
            if rec.get('output_identity')!=identity(shard):raise ValueError('AVI shard changed after completion; rerun scan: '+str(shard))
            if parameters is None:parameters=rec['parameters']
            elif parameters!=rec['parameters']:raise ValueError('AVI shards use different score/window parameters')
            if worker is None:worker=rec.get('worker');native=rec.get('native')
            elif worker!=rec.get('worker') or native!=rec.get('native'):raise ValueError('AVI shards use different worker versions')
            with open_text(shard) as f:
                reader=csv.DictReader(f,delimiter='\t')
                if writer is None:fields=reader.fieldnames;writer=csv.DictWriter(dest,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
                for row in reader:
                    writer.writerow(row);w=int(row['window_bp']);counts[w]=counts.get(w,0)+1
                    if eligible(row,a):
                        for m in METRICS:
                            v=number(row[m])
                            if v is not None:values.setdefault((w,m),array.array('d')).append(v)
    os.replace(merged.with_suffix('.partial.gz'),merged)
    cutoffs=[];lookup={}
    for (w,m),v in sorted(values.items()):
        lo,hi=np.quantile(np.asarray(v),[a.low_quantile,a.high_quantile])
        row=dict(window_bp=w,metric=m,n_eligible_windows=len(v),low_cut=float(lo),high_cut=float(hi),
            low_quantile=a.low_quantile,high_quantile=a.high_quantile,status='ready' if len(v)>=a.min_windows else 'insufficient_windows')
        cutoffs.append(row);lookup[(w,m)]=row
    write_tsv(out/'results/cutoffs.tsv',['window_bp','metric','n_eligible_windows','low_cut','high_cut','low_quantile','high_quantile','status'],cutoffs)
    hot=[];top={};plotrows=[];primary=min(counts,key=lambda x:abs(x-100000))
    with open_text(merged) as f:
        for row in csv.DictReader(f,delimiter='\t'):
            w=int(row['window_bp'])
            if w==primary:plotrows.append(row)
            if not eligible(row,a):continue
            for m in METRICS:
                c=lookup.get((w,m));v=number(row[m])
                if not c or c['status']!='ready' or v is None or c['low_cut']==c['high_cut']:continue
                for direction in ('low','high'):
                    if (direction=='low' and v<=c['low_cut']) or (direction=='high' and v>=c['high_cut']):
                        label=f'AVI_{m}_{direction}_w{w}'
                        hot.append((row['chrom'],int(row['start']),int(row['end']),label))
                        candidates=top.setdefault(label,[]);candidates.append(dict(row,label=label,value=v))
                        if len(candidates)>400:
                            candidates.sort(key=lambda r:float(r['value']),reverse=direction=='high');del candidates[100:]
    index=Intervals(hot)
    write_tsv(out/'results/hotspots.bed',['#chrom','start0','end0','label'],index.rows())
    top_rows=[]
    for label,rows in top.items():
        rows.sort(key=lambda r:float(r['value']),reverse='_high_' in label);top_rows.extend(rows[:100])
    write_tsv(out/'results/top_windows.tsv',['label','chrom','start','end','window_bp','value','coverage_fraction','n_complete_positions'],top_rows)
    plot(a,plotrows)
    report=dict(schema='allele3-avi-scan-1',status='complete',genome_build='GRCh38',
        chromosomes=[p.name.split('.')[0] for p in shards],score_parameters=parameters,
        eligibility=dict(min_coverage=a.min_coverage,min_positions=a.min_positions,min_complete_fraction=a.min_complete_fraction,min_windows=a.min_windows),
        windows_by_scale=counts,semantics='All possible nonreference substitutions; prediction availability is not population observation; missing scores are never zero-filled.',
        gradient='Adjacent scored genomic positions within the same window only; gaps and window boundaries are not connected.',
        cutoffs='Empirical quantiles of eligible windows, independent of target allele status; not significance tests.',
        outputs=dict(windows=str(merged),cutoffs=str(out/'results/cutoffs.tsv'),hotspots=str(out/'results/hotspots.bed')))
    atomic_json(data_path(out, 'qc') / 'manifest.json',report)
    html_report(out/'reports/avi_scan.html','AlphaGenome AVI landscape — GRCh38',[
        ('Scope','<p>Prediction landscape over the requested chromosomes. These windows do not imply observed SNPs or balancing selection.</p>'),
        ('Genome-wide landscape','<img src="../figures/avi_landscape.png">'),
        ('Window distribution','<img src="../figures/avi_density.png">'),
        ('Quantile thresholds',html_table(cutoffs,['window_bp','metric','n_eligible_windows','low_cut','high_cut','status'])),
        ('Definitions','<p>Between-position SD, within-position ALT range and adjacent-position gradient are separate metrics. High/low substitution fractions use the configured score cutoffs. Sparse or incompletely scored windows are excluded from hotspot definitions.</p><p><a href="../results/windows.tsv.gz">All windows</a> · <a href="../results/hotspots.bed">Hotspot BED0</a> · <a href="../results/top_windows.tsv">Top windows</a></p>')])
    print(json.dumps(dict(status='complete',out=str(out),windows=counts)));return report


def plot(a,rows):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    out=Path(a.out)/'figures';out.mkdir(parents=True,exist_ok=True)
    rows=sorted(rows,key=lambda r:(chrom_key(r['chrom']),int(r['start'])))
    fig,axes=plt.subplots(3,1,figsize=(14,8),sharex=True)
    offset=0;xs=[];labels=[];ticks=[];last=None;maximum=0
    for r in rows:
        c=r['chrom']
        if c!=last:
            if last is not None:offset+=maximum+1000000
            labels.append(c.removeprefix('chr'));ticks.append(offset);maximum=0;last=c
        xs.append(offset+(int(r['start'])+int(r['end']))/2);maximum=max(maximum,int(r['end']))
    for ax,m,title in zip(axes,['mean_avi','sd_site_mean','mean_abs_gradient'],['Mean predicted impact','Between-position SD','Adjacent-position absolute gradient']):
        vals=[number(r[m]) for r in rows]
        ax.scatter(xs,[v if v is not None else np.nan for v in vals],s=3,color='#246f91',rasterized=True);ax.set_ylabel(title);ax.grid(alpha=.15)
    axes[-1].set_xticks(ticks,labels);axes[-1].set_xlabel('Chromosome (GRCh38)');fig.tight_layout();fig.savefig(out/'avi_landscape.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5));vals=[float(r['mean_avi']) for r in rows if eligible(r,a)]
    if vals:ax.hist(vals,bins=min(80,max(5,int(len(vals)**.5))),density=True,color='#298d7f',alpha=.85)
    else:ax.text(.5,.5,'No windows pass the declared coverage criteria',ha='center',transform=ax.transAxes)
    ax.set(xlabel='Window mean AVI ('+a.score_field+')',ylabel='Density of eligible windows',title='Distribution of window means; not individual-SNV score density');fig.tight_layout();fig.savefig(out/'avi_density.png',dpi=180);plt.close(fig)


def parser():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['scan','merge','plan'],nargs='?',default='scan')
    p.add_argument('--vcf-dir',default=os.environ.get('ALLELE_ATLAS_DIR',configured_path('DIR_ANNOT', 'alphagenome/38')));p.add_argument('--vcf')
    p.add_argument('--out',default=os.environ.get('ALLELE_AVI_STORE',configured_path('DIR_ANNOT', 'alphagenome/38/avi_scan')))
    p.add_argument('--chr',default='autosomes');p.add_argument('--window-sizes',default='1000,10000,100000,1000000')
    p.add_argument('--score-field',choices=['AG_AVI_RAW','AG_AVI_PHRED'],default='AG_AVI_PHRED')
    p.add_argument('--high',type=float,default=20);p.add_argument('--low',type=float,default=5)
    p.add_argument('--min-coverage',type=float,default=.8);p.add_argument('--min-complete-fraction',type=float,default=.95)
    p.add_argument('--min-positions',type=int,default=10);p.add_argument('--min-windows',type=int,default=10)
    p.add_argument('--low-quantile',type=float,default=.05);p.add_argument('--high-quantile',type=float,default=.95)
    p.add_argument('--cores',type=int);p.add_argument('--jobs',type=int);p.add_argument('--replace',choices=['TRUE','FALSE'],default='FALSE');p.add_argument('--native-binary')
    return p


def main(argv=None):
    a=parser().parse_args(argv);chosen=chromosomes(a.chr)
    budget=core_budget(a.cores)
    if a.jobs is not None and not 1 <= a.jobs <= budget: raise ValueError('--jobs must be within --cores')
    a.jobs=min(a.jobs or budget, len(chosen), worker_cores(budget))
    sizes=[int(x) for x in a.window_sizes.split(',')]
    if min(sizes)<1 or len(set(sizes))!=len(sizes):raise ValueError('Unique positive --window-sizes required')
    if a.jobs<1 or not 0<=a.low_quantile<a.high_quantile<=1:raise ValueError('Invalid jobs/quantiles')
    if not 0<=a.min_coverage<=1 or not 0<=a.min_complete_fraction<=1:raise ValueError('Coverage fractions must lie in [0,1]')
    if a.high<a.low:raise ValueError('--high must be >= --low')
    if a.vcf and len(chosen)!=1:raise ValueError('--vcf requires one --chr')
    if a.command=='plan':
        print(json.dumps(dict(chromosomes=chosen,vcf=a.vcf,vcf_dir=a.vcf_dir,out=a.out,windows=sizes,score_field=a.score_field,mode='one streaming native worker per chromosome'),indent=2));return
    from output_layout import prepare
    prepare(a.out)
    manifest=data_path(Path(a.out), 'qc') / 'manifest.json'
    atomic_json(manifest,dict(status='running',genome_build='GRCh38',command=a.command))
    try:
        if a.command=='merge':return finish(a.out,merge(a,chosen))
        sources={c:Path(a.vcf) if a.vcf else resolve_vcf(a.vcf_dir,c) for c in chosen}
        binary=native_binary(data_path(Path(a.out), 'tmp') / 'bin',a.native_binary)
        with ThreadPoolExecutor(max_workers=a.jobs) as pool:
            fs=[pool.submit(run_chromosome,a,c,p,binary) for c,p in sources.items()]
            for f in as_completed(fs):f.result()
        return finish(a.out,merge(a,chosen))
    except Exception as e:
        atomic_json(manifest,dict(status='failed',genome_build='GRCh38',reason=str(e)));raise


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as e:raise SystemExit('avi-scan: '+str(e))
