#!/usr/bin/env python3
"""Official PRS-CSx posterior meta scores; reuse the matching joint MCMC run."""
import argparse,concurrent.futures,fcntl,json,os,shutil,subprocess,threading
from pathlib import Path
import numpy as np
import pandas as pd
from score_output import filter_samples,update_csx_table
from pipeline_io import inspect,coverage,signature,stamp
from csx_config import ROOT,POPS,find_gwas,inference_signature,sample_size,validate_mcmc

def genotype(prefix):
    prefix=Path(prefix)
    if Path(str(prefix)+'.pgen').is_file():
        var=Path(str(prefix)+'.pvar');extra=[]
        if not var.is_file():var=Path(str(prefix)+'.pvar.zst');extra=['vzs']
        files=[Path(str(prefix)+'.pgen'),Path(str(prefix)+'.psam'),var]
        command=['--pfile',str(prefix),*extra]
    else:
        files=[Path(str(prefix)+ext) for ext in ('.bed','.bim','.fam')]
        command=['--bfile',str(prefix)]
    for f in files:stamp(f)
    return command,files

def main():
    p=argparse.ArgumentParser()
    for key in ['trait','mode','gwas-dir','target-dir','snpinfo','ref-dir','bim','work','score-home','chrs']:p.add_argument('--'+key,required=True)
    for key,default in [('jobs',4),('threads',1),('iterations',4000),('burnin',2000),('thin',5),('seed',20260904)]:p.add_argument('--'+key,type=int,default=default)
    p.add_argument('--phi',default='1e-2');p.add_argument('--n-gwas',default='')
    p.add_argument('--remove',default='');p.add_argument('--keep',default='')
    p.add_argument('--stage',choices=['all','weights','score'],default='all')
    p.add_argument('--replace',choices=['TRUE','FALSE'],default='FALSE');p.add_argument('--check',action='store_true')
    a=p.parse_args()
    if a.mode not in ('auto','meta') or min(a.jobs,a.threads)<1:raise ValueError('Invalid mode/jobs/threads')
    phi='auto' if a.mode=='auto' else a.phi
    if a.mode=='meta' and phi=='auto':raise ValueError('meta requires fixed phi')
    validate_mcmc(phi,a.iterations,a.burnin,a.thin,a.seed)
    chrs=list(map(int,a.chrs.split()))
    if not chrs or len(set(chrs))!=len(chrs) or not set(chrs)<=set(range(1,23)):raise ValueError('Invalid chromosomes')
    inputs=[find_gwas(a.gwas_dir,a.trait,pop) for pop in POPS]
    name=Path(a.snpinfo).name
    if name not in ('snpinfo_mult_1kg_hm3','snpinfo_mult_ukbb_hm3'):raise ValueError('Unrecognized SNPINFO reference type')
    ref_type='1kg' if name=='snpinfo_mult_1kg_hm3' else 'ukbb';refs=[];ld=[]
    for pop in POPS:
        d=Path(a.ref_dir)/f'ldblk_{ref_type}_{pop.lower()}'
        if not d.is_dir():d=Path(a.ref_dir)/f'ldblk_{ref_type}_{pop}'
        refs.append(d.resolve());ld += [d/f'ldblk_{ref_type}_chr{c}.hdf5' for c in chrs]
    sig=inference_signature(phi,a.iterations,a.burnin,a.thin,a.seed,a.n_gwas,chrs,inputs,a.snpinfo,a.bim,ld)
    inspect(a.snpinfo,list(map(str,inputs)));coverage(a.chrs,list(map(str,inputs)))
    gen={};targetfiles=[]
    if a.stage!='weights':
        if not shutil.which('plink2'):raise ValueError('plink2 not found')
        for c in chrs:
            gen[c],files=genotype(Path(a.target_dir)/f'chr{c}');targetfiles+=files
        for f in (a.remove,a.keep):
            if f and not Path(f).is_file():raise FileNotFoundError(f)
    home=Path(a.score_home);work=Path(a.work)
    weight=home/'.weights'/f'csx.{a.mode}.gz';meta=weight.with_suffix('.json')
    ready=weight.is_file() and meta.is_file() and json.loads(meta.read_text()).get('signature')==sig
    if a.stage=='score' and not ready:raise ValueError('Matching combined weights missing; run --stage all or weights')
    print(f'CSx {a.mode}: phi={phi}; shared joint-run signature={sig[:12]}',flush=True)
    if a.check:return
    run=work/a.trait/sig;run.mkdir(parents=True,exist_ok=True)
    lock=(work/a.trait/'run.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    logs=run/'logs';logs.mkdir(exist_ok=True);command_lock=threading.Lock()
    env=dict(os.environ,OMP_NUM_THREADS=str(a.threads),OPENBLAS_NUM_THREADS=str(a.threads),MKL_NUM_THREADS=str(a.threads))
    def execute(cmd,name):
        cmd=list(map(str,cmd))
        with command_lock:
            with (run/'commands.jsonl').open('a') as f:f.write(json.dumps(cmd)+'\n')
        with (logs/(name+'.log')).open('w') as f:rc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
        if rc:raise RuntimeError(f'Command exit={rc}; see {logs/(name+".log")}')
    def posterior(c):
        found=list((run/'raw'/f'chr{c}').glob('*META*pst_eff*.txt'))
        return found[0] if len(found)==1 and found[0].stat().st_size else None
    if a.stage!='score' and (not ready or a.replace=='TRUE'):
        raw=run/'raw';raw.mkdir(exist_ok=True)
        done=all(posterior(c) is not None and (raw/f'chr{c}'/'done').is_file() for c in chrs)
        if done and a.replace=='FALSE':print('SKIP inference: reuse matching joint population/META posteriors',flush=True)
        else:
            ref=run/'reference';ref.mkdir(exist_ok=True)
            for src,nm in [(Path(a.snpinfo).resolve(),name)]+[(d,f'ldblk_{ref_type}_{pop.lower()}') for d,pop in zip(refs,POPS)]:
                dest=ref/nm
                if not dest.exists():dest.symlink_to(src)
            prep=run/'sumstats';prep.mkdir(exist_ok=True);sizes=[]
            for pop,src in zip(POPS,inputs):
                dst=prep/f'{pop}.tsv.gz';info=prep/f'{pop}.json'
                execute(['python3',ROOT/'sumstats_cache.py','--input',src,'--output',dst,'--metadata',info,'--snpinfo',a.snpinfo,'--trait',a.trait,'--pop',pop,'--work',work,'--replace',a.replace],f'prepare.{pop}')
                print((logs/f'prepare.{pop}.log').read_text().splitlines()[-1],flush=True)
                sizes.append(str(sample_size(json.loads(info.read_text()),a.n_gwas,pop)))
                execute(['python3',ROOT/'split_sumstats.py','--input',dst,'--out-dir',prep,'--prefix',pop,'--chrs',a.chrs],f'split.{pop}')
            def infer(c):
                dest=raw/f'chr{c}';dest.mkdir(exist_ok=True);marker=dest/'done'
                if marker.exists() and posterior(c) is not None and a.replace=='FALSE':return
                marker.unlink(missing_ok=True)
                sst=','.join(str(prep/f'{pop}.chr{c}.tsv') for pop in POPS)
                cmd=['python3',ROOT/'csx/PRScsx.py',f'--ref_dir={ref}',f'--bim_prefix={a.bim}',f'--sst_file={sst}',f'--n_gwas={",".join(sizes)}',f'--pop={",".join(POPS)}',f'--chrom={c}',f'--n_iter={a.iterations}',f'--n_burnin={a.burnin}',f'--thin={a.thin}',f'--seed={a.seed+c}',f'--out_dir={dest}',f'--out_name={a.trait}','--meta=TRUE']
                if phi!='auto':cmd.append(f'--phi={phi}')
                print(f'RUN {a.trait} {a.mode} chr{c}',flush=True);execute(cmd,f'infer.chr{c}')
                if posterior(c) is None:raise ValueError(f'Missing META posterior chr{c}')
                marker.touch()
            with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:list(ex.map(infer,chrs))
        from normalize_csx_weights import read
        z=pd.concat([read(posterior(c)) for c in chrs],ignore_index=True)
        if z.empty or z.SNP.duplicated().any() or not np.isfinite(z.BETA).all():raise ValueError('Invalid combined weights')
        weight.parent.mkdir(parents=True,exist_ok=True)
        tmp=weight.with_suffix('.tmp');z.to_csv(tmp,sep='\t',index=False,compression='gzip');tmp.replace(weight)
        meta.write_text(json.dumps({'signature':sig,'phi':phi,'inputs':list(map(str,inputs))},indent=2)+'\n')
    if a.stage=='weights':return
    scoring_sig=signature([sig,a.mode],targetfiles+[weight]+[Path(x) for x in (a.remove,a.keep) if x]+[Path(__file__),ROOT/'combine_scores.py',ROOT/'score_output.py'])
    score=run/'combined_scores'/a.mode/scoring_sig;score.mkdir(parents=True,exist_ok=True)
    def scoring(c):
        out=score/f'chr{c}';done=Path(str(out)+'.done')
        if done.exists() and Path(str(out)+'.sscore').is_file() and a.replace=='FALSE':return
        done.unlink(missing_ok=True)
        cmd=['plink2',*gen[c],'--score',weight,'1','2','3','header-read','no-mean-imputation','list-variants','cols=+scoresums','--threads',a.threads,'--memory','4096','--out',out]
        if a.remove and Path(a.remove).stat().st_size:cmd+=['--remove',a.remove]
        if a.keep:cmd+=['--keep',a.keep]
        execute(cmd,f'{a.mode}.score.chr{c}')
        if not Path(str(out)+'.sscore').is_file():raise ValueError(f'Missing scores: {out}')
        done.touch()
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:list(ex.map(scoring,chrs))
    combined=score/'combined.gz'
    execute(['python3',ROOT/'combine_scores.py','--inputs',*[score/f'chr{c}.sscore' for c in chrs],'--name',f'csx.{a.mode}','--output',combined],f'{a.mode}.combine')
    z=filter_samples(pd.read_csv(combined,sep='\t',dtype={'eid':str}),'eid',a.remove)
    output=home/'csx.pgs.gz';update_csx_table(z,output,a.remove)
    print(f'DONE {a.trait}: csx.{a.mode}; N={len(z)}; {output}',flush=True)

if __name__=='__main__':main()
