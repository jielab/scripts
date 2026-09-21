#!/usr/bin/env python3
"""Official PRS-CSx posterior meta scores, with learned or fixed phi."""
import argparse, concurrent.futures, fcntl, gzip, hashlib, json, os
from pathlib import Path
import subprocess, threading
import numpy as np
import pandas as pd
from score_output import filter_samples

ROOT = Path(__file__).resolve().parent
POPS = ['AFR','EAS','EUR','SAS']

def main():
    p=argparse.ArgumentParser()
    for key in ['trait','mode','gwas-dir','target-dir','snpinfo','ref-dir','bim','work','score-home','chrs']:
        p.add_argument('--'+key,required=True)
    for key,default in [('jobs',4),('threads',1),('iterations',4000),('burnin',2000),('thin',5),('seed',20260904)]:
        p.add_argument('--'+key,type=int,default=default)
    p.add_argument('--phi',default='1e-2');p.add_argument('--n-gwas',default='')
    p.add_argument('--remove',default='');p.add_argument('--keep',default='')
    p.add_argument('--stage',default='all');p.add_argument('--replace',default='FALSE');p.add_argument('--check',action='store_true')
    a=p.parse_args();assert a.mode in ['auto','meta']
    assert a.iterations>a.burnin>=0 and a.thin>0 and a.jobs>0 and a.threads>0
    if a.mode=='meta' and (a.phi=='auto' or float(a.phi)<=0):raise ValueError('meta requires a positive fixed --phi')
    chrs=list(map(int,a.chrs.split()));home=Path(a.score_home);work=Path(a.work)
    inputs=[]
    for pop in POPS:
        f=Path(a.gwas_dir)/f'{a.trait}.{pop}'/'gwas'/f'{a.trait}.{pop}.gz'
        if not f.exists() and a.trait=='t2dm' and pop=='AFR':f=Path(a.gwas_dir)/'t2dm.AFA/gwas/t2dm.AFA.gz'
        if not f.is_file():raise ValueError(f'Missing GWAS: {f}')
        inputs.append(f);print('input GWAS:',f,flush=True)
    ref_type='1kg' if Path(a.snpinfo).name=='snpinfo_mult_1kg_hm3' else 'ukbb'
    refs=[];ld=[]
    for pop in POPS:
        d=Path(a.ref_dir)/f'ldblk_{ref_type}_{pop.lower()}'
        if not d.is_dir():d=Path(a.ref_dir)/f'ldblk_{ref_type}_{pop}'
        refs.append(d);ld += [d/f'ldblk_{ref_type}_chr{c}.hdf5' for c in chrs]
    def stamp(f):
        f=Path(f);s=f.stat();return [str(f.resolve()),s.st_size,s.st_mtime_ns]
    settings=[a.mode,a.phi if a.mode=='meta' else 'auto',a.iterations,a.burnin,a.thin,a.seed,a.n_gwas,chrs]
    files=inputs+[Path(a.snpinfo),Path(a.bim+'.bim')]+ld+[ROOT/'csx'/n for n in ['PRScsx.py','parse_genet.py','mcmc_gtb.py','gigrnd.py']]+[ROOT/'prepare_sumstats.py',Path(__file__)]
    signature=hashlib.sha256(json.dumps([settings,[stamp(f) for f in files]]).encode()).hexdigest()
    print(f'CSx {a.mode}: official --meta=True; phi={settings[1]}',flush=True)
    print(f'output score file: {home}/csx.pgs.gz; column=csx.{a.mode}',flush=True)
    if a.check:return
    run=work/a.trait/'combined'/a.mode/signature;run.mkdir(parents=True,exist_ok=True)
    lock=(work/a.trait/'run.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    logs=work/'log'/a.trait/a.mode;logs.mkdir(parents=True,exist_ok=True)
    command_lock=threading.Lock()
    env=dict(os.environ,OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS=str(a.threads),MKL_NUM_THREADS=str(a.threads))
    def execute(cmd,name):
        cmd=list(map(str,cmd))
        with command_lock:
            with (run/'commands.jsonl').open('a') as f:f.write(json.dumps(cmd)+'\n')
        with (logs/(name+'.log')).open('w') as f:
            rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
        if rc:raise RuntimeError(f'Command exit={rc}; see {logs/(name+".log")}')
    permanent=home/'.weights';permanent.mkdir(parents=True,exist_ok=True)
    weight=permanent/f'csx.{a.mode}.gz';meta=permanent/f'csx.{a.mode}.json'
    ready=weight.is_file() and meta.is_file() and json.loads(meta.read_text()).get('signature')==signature
    if a.stage=='score' and not ready:raise ValueError('Matching combined weights missing; run --stage all or weights')
    if a.stage!='score' and (not ready or a.replace=='TRUE'):
        ref=run/'reference';ref.mkdir(exist_ok=True)
        for src,name in [(Path(a.snpinfo),Path(a.snpinfo).name)]+[(d,f'ldblk_{ref_type}_{pop.lower()}') for d,pop in zip(refs,POPS)]:
            dest=ref/name
            if not dest.exists():dest.symlink_to(src)
        n=[];std=[]
        for pop,src in zip(POPS,inputs):
            dst=run/f'{pop}.tsv.gz';info=run/f'{pop}.json'
            if not dst.exists() or not info.exists():
                # Reuse compatible normalized inputs; never infer identity from a filename alone.
                candidates=sorted((work/a.trait).glob(f'*/sumstats/{pop}.json'),key=lambda x:x.stat().st_mtime,reverse=True)
                cached=None
                for c in candidates:
                    data=json.loads(c.read_text());table=c.with_suffix('.tsv.gz')
                    if data.get('input')==str(src.resolve()) and table.exists() and table.stat().st_mtime_ns>=src.stat().st_mtime_ns:
                        cached=(table,c);break
                if cached:
                    for target,source in [(dst,cached[0]),(info,cached[1])]:
                        if not target.exists():target.symlink_to(source)
                else:execute(['python3',ROOT/'prepare_sumstats.py','--input',src,'--output',dst,'--metadata',info,'--snpinfo',a.snpinfo,'--trait',a.trait,'--pop',pop],f'prepare.{pop}')
            value=json.loads(info.read_text())['n_gwas_median']
            if a.n_gwas:
                value=dict(x.split('=') for x in a.n_gwas.split(',')).get(pop,value) if '=' in a.n_gwas else a.n_gwas
            if value is None or float(value)<=0:raise ValueError(f'Missing N for {pop}')
            n.append(str(round(float(value))));std.append(dst)
            execute(['python3',ROOT/'split_sumstats.py','--input',dst,'--out-dir',run/'sumstats','--prefix',pop,'--chrs',a.chrs],f'split.{pop}')
        raw=run/'raw';raw.mkdir(exist_ok=True)
        def infer(c):
            dest=raw/f'chr{c}';dest.mkdir(exist_ok=True);done=dest/'done'
            if done.exists() and list(dest.glob('*META*pst_eff*.txt')) and a.replace!='TRUE':return
            sst=','.join(str(run/'sumstats'/f'{pop}.chr{c}.tsv') for pop in POPS)
            cmd=['python3',ROOT/'csx/PRScsx.py',f'--ref_dir={ref}',f'--bim_prefix={a.bim}',f'--sst_file={sst}',f'--n_gwas={",".join(n)}',f'--pop={",".join(POPS)}',f'--chrom={c}',f'--n_iter={a.iterations}',f'--n_burnin={a.burnin}',f'--thin={a.thin}',f'--seed={a.seed+c}',f'--out_dir={dest}',f'--out_name={a.trait}','--meta=True']
            if a.mode=='meta':cmd.append(f'--phi={a.phi}')
            print(f'RUN {a.trait} {a.mode} chr{c}',flush=True);execute(cmd,f'infer.chr{c}');done.touch()
        with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:list(ex.map(infer,chrs))
        from normalize_csx_weights import read
        tables=[]
        for c in chrs:
            found=list((raw/f'chr{c}').glob('*META*pst_eff*.txt'))
            if len(found)!=1:raise ValueError(f'Expected one META posterior for chr{c}')
            tables.append(read(found[0]))
        z=pd.concat(tables,ignore_index=True)
        if z.empty or z.SNP.duplicated().any() or not np.isfinite(z.BETA).all():raise ValueError('Invalid combined weights')
        tmp=weight.with_suffix('.tmp');z.to_csv(tmp,sep='\t',index=False,compression='gzip');tmp.replace(weight)
        meta.write_text(json.dumps({'signature':signature,'settings':settings,'inputs':list(map(str,inputs))},indent=2)+'\n')
    if a.stage=='weights':return
    targetfiles=[]
    for c in chrs:
        prefix=Path(a.target_dir)/f'chr{c}'
        targetfiles.extend([Path(str(prefix)+s) for s in ['.pgen','.psam']])
        targetfiles.append(Path(str(prefix)+('.pvar' if Path(str(prefix)+'.pvar').exists() else '.pvar.zst')))
    scoring_sig=hashlib.sha256(json.dumps([signature,[stamp(f) for f in targetfiles+[weight]+[Path(x) for x in [a.remove,a.keep] if x]]]).encode()).hexdigest()
    score=run/'scores'/scoring_sig;score.mkdir(parents=True,exist_ok=True)
    def scoring(c):
        out=score/f'chr{c}';done=Path(str(out)+'.done')
        if done.exists() and Path(str(out)+'.sscore').exists() and a.replace!='TRUE':return
        prefix=Path(a.target_dir)/f'chr{c}';cmd=['plink2','--pfile',prefix]
        if not Path(str(prefix)+'.pvar').exists():cmd.append('vzs')
        cmd+=['--score',weight,'1','2','3','header-read','no-mean-imputation','cols=+scoresums','--threads',a.threads,'--memory','4096','--out',out]
        if a.remove:cmd+=['--remove',a.remove]
        if a.keep:cmd+=['--keep',a.keep]
        execute(cmd,f'score.chr{c}');done.touch()
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:list(ex.map(scoring,chrs))
    combined=score/'combined.gz'
    execute(['python3',ROOT/'combine_scores.py','--inputs',*[score/f'chr{c}.sscore' for c in chrs],'--name',f'csx.{a.mode}','--output',combined],'combine')
    z=pd.read_csv(combined,sep='\t',dtype={'eid':str});z=filter_samples(z,'eid',a.remove)
    output=home/'csx.pgs.gz'
    if output.exists():
        old=pd.read_csv(output,sep='\t',dtype={'eid':str});old=filter_samples(old,'eid',a.remove)
        if set(old.eid)!=set(z.eid):raise ValueError('Combined score sample set differs from existing CSx scores')
        z=old.drop(columns=[f'csx.{a.mode}'],errors='ignore').merge(z,on='eid',validate='one_to_one')
    tmp=output.with_suffix('.gz.tmp');z.to_csv(tmp,sep='\t',index=False,compression='gzip');tmp.replace(output)
    print(f'DONE {a.trait}: csx.{a.mode}; N={len(z)}; {output}',flush=True)

if __name__=='__main__':main()
