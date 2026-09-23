#!/usr/bin/env python3
"""Score synchronized PRS-CSx draws; retain individual joint covariance.

One chromosome at a time. Cross-chromosome covariances are exactly zero in
the fitted chromosome-factorized model; iteration numbers from independent
chromosomes must NOT be treated as coupled joint draws.
"""
import argparse, hashlib, json, os, subprocess
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from csx_combined import genotype
from csx_config import POPS
from pipeline_io import signature

PAIRS = [(j,k) for j in range(4) for k in range(j,4)]
COMP = str.maketrans('ACGT','TGCA')

def align_frequencies(g, table):
    """EAF refers to table A1; orient it to the scored A1, never use MAF."""
    if not {'SNP','A1','A2','EAF'} <= set(table):
        raise ValueError('Frequency table needs SNP,A1,A2,EAF')
    if table.SNP.duplicated().any(): raise ValueError('Duplicate frequency SNPs')
    ids = g['SNP'].asstr()[:]
    d = table.set_index('SNP').reindex(ids)
    af = pd.to_numeric(d.EAF, errors='coerce').to_numpy()
    a1,a2 = g['A1'].asstr()[:],g['A2'].asstr()[:]
    f1,f2 = d.A1.fillna('').str.upper().to_numpy(),d.A2.fillna('').str.upper().to_numpy()
    c1=np.array([s.translate(COMP) for s in f1]);c2=np.array([s.translate(COMP) for s in f2])
    same=((a1==f1)&(a2==f2))|((a1==c1)&(a2==c2))
    flip=((a1==f2)&(a2==f1))|((a1==c2)&(a2==c1))
    valid=(same^flip)&np.isfinite(af)&(af>0)&(af<1)
    if not valid.all():
        raise ValueError(f'{np.sum(~valid)} posterior SNPs lack unambiguous GWAS EAF; first: {ids[~valid][:5].tolist()}')
    return np.where(flip,1-af,af)

def pack_scores(sscore, dest, ndraw):
    """Convert PLINK score sums to row-chunked HDF5, with strict ID checks."""
    header=pd.read_csv(sscore,sep=r'\s+',nrows=0).columns
    idcol='IID' if 'IID' in header else '#IID'
    cols=[f'DRAW_{i:04d}_SUM' for i in range(ndraw)]
    if not set(cols+[idcol])<=set(header): raise ValueError('PLINK posterior score columns are incomplete')
    n=0;seen=set()
    with h5py.File(str(dest)+'.tmp','w') as h:
        values=h.create_dataset('scores',(0,ndraw),maxshape=(None,ndraw),chunks=(256,ndraw),dtype='f8',compression='lzf')
        ids=h.create_dataset('eid',(0,),maxshape=(None,),dtype=h5py.string_dtype())
        for d in pd.read_csv(sscore,sep=r'\s+',usecols=[idcol]+cols,dtype={idcol:str},chunksize=512):
            ix=d[idcol]
            if ix.isna().any() or (ix=='').any() or ix.duplicated().any() or seen.intersection(ix): raise ValueError('Missing/duplicate posterior IDs')
            ix=ix.astype(str);seen.update(ix);v=d[cols].to_numpy(float)
            if not np.isfinite(v).all(): raise ValueError('Nonfinite posterior scores')
            values.resize(n+len(d),axis=0);ids.resize(n+len(d),axis=0)
            values[n:n+len(d)]=v;ids[n:n+len(d)]=ix.to_numpy();n+=len(d)
        if not n:raise ValueError('No posterior score participants')
        h.attrs['complete']=True
    os.replace(str(dest)+'.tmp',dest)

def moments(pop_files, output):
    handles=[h5py.File(p,'r') for p in pop_files]
    try:
        shape=handles[0]['scores'].shape
        if shape[1]<20 or any(h['scores'].shape!=shape or not h.attrs.get('complete') for h in handles):
            raise ValueError('Need >=20 aligned draws and identical score dimensions')
        n,b=shape
        with h5py.File(str(output)+'.tmp','w') as out:
            out.create_dataset('eid',shape=(n,),dtype=h5py.string_dtype())
            out.create_dataset('mean',shape=(n,4),dtype='f8',chunks=True,compression='lzf')
            out.create_dataset('cov',shape=(n,10),dtype='f8',chunks=True,compression='lzf')
            for start in range(0,n,512):
                sl=slice(start,min(n,start+512));ids=handles[0]['eid'].asstr()[sl]
                if any(not np.array_equal(ids,h['eid'].asstr()[sl]) for h in handles[1:]):
                    raise ValueError('Population score sample order differs')
                v=np.stack([h['scores'][sl] for h in handles],axis=2)
                mean=v.mean(axis=1);u=v-mean[:,None,:]
                cov=np.einsum('nbj,nbk->njk',u,u)/(b-1)
                out['eid'][sl]=ids;out['mean'][sl]=mean
                out['cov'][sl]=np.column_stack([cov[:,j,k] for j,k in PAIRS])
            out.attrs.update(complete=True,ndraw=b)
        os.replace(str(output)+'.tmp',output)
    finally:
        for h in handles:h.close()

def combine_chromosomes(files, output):
    handles=[h5py.File(p,'r') for p in files]
    try:
        n=len(handles[0]['eid']);output=Path(output);tmp=output.with_name(output.name+'.tmp.gz')
        if any(len(h['eid'])!=n or not h.attrs.get('complete') for h in handles):raise ValueError('Incomplete chromosome moments')
        first=True
        for start in range(0,n,4096):
            sl=slice(start,min(start+4096,n));ids=handles[0]['eid'].asstr()[sl]
            if any(not np.array_equal(ids,h['eid'].asstr()[sl]) for h in handles[1:]):raise ValueError('Chromosome score sample order differs')
            mu=sum(h['mean'][sl] for h in handles)
            cov=sum(h['cov'][sl] for h in handles)
            d=pd.DataFrame({'eid':ids})
            for j,p in enumerate(POPS):d[f'csx.{p}']=mu[:,j]
            for j,(p,q) in enumerate(PAIRS):d[f'cov.{POPS[p]}.{POPS[q]}']=cov[:,j]
            d.to_csv(tmp,sep='\t',index=False,mode='w' if first else 'a',header=first,compression='gzip',float_format='%.12g');first=False
        tmp.replace(output)
    finally:
        for h in handles:h.close()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('raw-dir','sumstats-dir','target-dir','output','chrs'):p.add_argument('--'+key,required=True)
    p.add_argument('--frequency-dir',help='Optional POP.tsv.gz files with discovery SNP,A1,A2,EAF; defaults to normalized GWAS')
    p.add_argument('--threads',type=int,default=4);p.add_argument('--memory',type=int,default=8192)
    p.add_argument('--remove',default='');p.add_argument('--keep',default='')
    p.add_argument('--replace',choices=['TRUE','FALSE'],default='FALSE')
    p.add_argument('--min-draws',type=int,default=100)
    a=p.parse_args();chrs=list(map(int,a.chrs.split()))
    if not chrs or len(set(chrs))!=len(chrs) or not set(chrs)<=set(range(1,23)):raise ValueError('Invalid chromosomes')
    if min(a.threads,a.memory,a.min_draws)<1:raise ValueError('Invalid resources/draw count')
    output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
    frequency={};freq_paths=[]
    for pop in POPS:
        path=Path(a.frequency_dir or a.sumstats_dir)/f'{pop}.tsv.gz';freq_paths.append(path)
        frequency[pop]=pd.read_csv(path,sep='\t',usecols=['SNP','A1','A2','EAF'],dtype={'SNP':str,'A1':str,'A2':str})
    raw=[Path(a.raw_dir)/f'chr{c}'/'joint_posterior.h5' for c in chrs]
    gen={};gen_files=[]
    for c in chrs:gen[c],fs=genotype(Path(a.target_dir)/f'chr{c}');gen_files.extend(fs)
    sig=signature([chrs,a.min_draws],raw+freq_paths+gen_files+[Path(__file__)]+[Path(x) for x in (a.keep,a.remove) if x])
    meta_path=Path(str(output)+'.json')
    if a.replace=='FALSE' and output.is_file() and meta_path.is_file() and Path(str(output)+'.metadata.tsv').is_file() and json.loads(meta_path.read_text()).get('signature')==sig:
        print('SKIP posterior scores: matching individual covariance');return
    work=Path(a.raw_dir).parent/'posterior_scores'/sig;work.mkdir(parents=True,exist_ok=True)
    chr_files=[];draw_counts=[]
    for c,source in zip(chrs,raw):
        dst=work/f'chr{c}';dst.mkdir(exist_ok=True);moment=dst/'moments.h5';chr_files.append(moment)
        with h5py.File(source,'r') as h:
            if not h.attrs.get('complete') or h.attrs.get('schema')!='grid_csx_draws_v1':raise ValueError(f'Incomplete draws: {source}')
            ndraw=len(h['iteration']);draw_counts.append(ndraw)
            if ndraw<a.min_draws:raise ValueError(f'chr{c}: {ndraw} draws < {a.min_draws}')
            if moment.is_file() and a.replace=='FALSE':continue
            print(f'START posterior scoring chr{c}: {ndraw} joint draws',flush=True)
            score_files=[]
            for pop in POPS:
                g=h[pop];eaf=align_frequencies(g,frequency[pop]);ids=g['SNP'].asstr()[:]
                w=dst/f'{pop}.weights.tsv';af=dst/f'{pop}.afreq';pref=dst/pop;cache=dst/f'{pop}.scores.h5';score_files.append(cache)
                if cache.is_file() and a.replace=='FALSE':continue
                # Fixed discovery allele means; never re-centre within the target group.
                pd.DataFrame({'#CHROM':c,'ID':ids,'REF':g['A2'].asstr()[:],'ALT':g['A1'].asstr()[:],'ALT_FREQS':eaf,'OBS_CT':2}).to_csv(af,sep='\t',index=False)
                columns=['SNP','A1']+[f'DRAW_{i:04d}' for i in range(ndraw)]
                for st in range(0,len(ids),2048):
                    sl=slice(st,min(st+2048,len(ids)));v=pd.DataFrame(g['beta'][sl],columns=columns[2:])
                    v.insert(0,'A1',g['A1'].asstr()[sl]);v.insert(0,'SNP',ids[sl])
                    v.to_csv(w,sep='\t',index=False,mode='w' if st==0 else 'a',header=st==0,float_format='%.12g')
                cmd=['plink2',*gen[c],'--extract',str(w),'--read-freq',str(af),'--error-on-freq-calc',
                     '--score',str(w),'1','2','header-read','center','no-mean-imputation','list-variants','cols=maybefid,scoresums',
                     '--score-col-nums',f'3-{ndraw+2}','--threads',str(a.threads),'--memory',str(a.memory),'--out',str(pref)]
                # Centred missing genotypes contribute zero, mathematically equal to
                # discovery-mean imputation. Read SUM, never AVG (its denominator
                # varies with missingness). Explicit no-mean-imputation also avoids
                # old PLINK builds adding an uncentred mean for missing calls.
                # --extract must contain only IDs, without a header or effect columns.
                snps=dst/f'{pop}.snps';snps.write_text('\n'.join(ids)+'\n');cmd[cmd.index('--extract')+1]=str(snps)
                if a.keep:cmd+=['--keep',a.keep]
                if a.remove and Path(a.remove).stat().st_size:cmd+=['--remove',a.remove]
                with (dst/f'{pop}.command.json').open('w') as f:json.dump(cmd,f)
                with (dst/f'{pop}.run.log').open('w') as f:
                    subprocess.run(cmd,check=True,stdout=f,stderr=subprocess.STDOUT)
                used=Path(str(pref)+'.sscore.vars').read_text().split()
                if len(used)!=len(ids) or set(used)!=set(ids):raise ValueError(f'{pop} chr{c}: PLINK skipped posterior SNPs; cannot claim complete uncertainty')
                pack_scores(str(pref)+'.sscore',cache,ndraw)
                for tmp in (w,af,snps,Path(str(pref)+'.sscore')):tmp.unlink(missing_ok=True)
            moments(score_files,moment)
            for file in score_files:file.unlink(missing_ok=True)
        print(f'DONE posterior scoring chr{c}',flush=True)
    combine_chromosomes(chr_files,output)
    metadata={'schema':'grid_csx_moments_v1','signature':sig,'populations':POPS,'chromosomes':chrs,'draw_counts':draw_counts,
              'centering':'discovery_EAF','covariance':'joint within chromosome; summed across independent chromosomes',
              'means':'posterior mean scores centred at discovery EAF','effect_scale':'standardized_phenotype_per_allele',
              'source_draws':[str(p.resolve()) for p in raw],'frequency_files':[str(p.resolve()) for p in freq_paths]}
    temp=Path(str(meta_path)+'.tmp');temp.write_text(json.dumps(metadata,indent=2)+'\n');temp.replace(meta_path)
    audit={key:metadata[key] for key in ('schema','centering','covariance','effect_scale')}
    audit.update(chromosomes=','.join(map(str,chrs)),draw_counts=','.join(map(str,draw_counts)))
    digest=hashlib.md5()
    with output.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):digest.update(block)
    audit['md5']=digest.hexdigest()
    mt=Path(str(output)+'.metadata.tsv.tmp')
    pd.DataFrame({'field':list(audit),'value':list(audit.values())}).to_csv(mt,sep='\t',index=False)
    mt.replace(str(output)+'.metadata.tsv')
    print(f'DONE individual posterior moments: {output}',flush=True)

if __name__=='__main__':main()
