#!/usr/bin/env python3
"""Generic plan/checkpoint and file-view helpers for gwas.sh.

GWAS options and commands are defined in the Bash entry point. Keep checkpoint
serialization stable so previously completed jobs remain resumable.
"""
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import re
from concurrent.futures import ThreadPoolExecutor

HERE=Path(__file__).resolve().parent
SELF=str(Path(__file__).resolve())
def execute(plan):
    # Serialize jobs in one plan; never let two invocations overwrite a shared cache.
    import fcntl
    lock=open(str(plan)+'.lock','w')
    try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: raise RuntimeError(f'Already running: {plan}')
    spec=json.loads(Path(plan).read_text())
    for job in spec['jobs']:
        executable=job['cmd'][0]
        if not shutil.which(executable): raise FileNotFoundError('Required executable: '+executable)
        if len(job['cmd'])>1 and job['cmd'][1].endswith('.R') and not Path(job['cmd'][1]).is_file():
            raise FileNotFoundError('Required R script: '+job['cmd'][1]+' (set --saige-dir/--saige-rscript)')
    def run_job(job):
        for f in job['inputs']:
            if not Path(f).is_file(): raise FileNotFoundError(f)
        cmd=job['cmd'][:]
        if job['covars']:
            path,method=job['covars']; names=Path(path).read_text().strip()
            if names:
                cmd+= {'plink2':['--covar',path.removesuffix('.covars'),'--covar-name',names,'--covar-variance-standardize'],
                       'regenie':['--covarFile',path.removesuffix('.covars'),'--covarColList',names],
                       'saige':['--covarColList='+names]}[method]
        signature=hashlib.sha256(json.dumps([cmd,[(f,Path(f).stat().st_size,Path(f).stat().st_mtime_ns) for f in job['inputs']]],sort_keys=True).encode()).hexdigest()
        stamp=Path(plan).parent/(job['name']+'.done.json')
        if not spec['replace'] and stamp.exists() and json.loads(stamp.read_text()).get('signature')==signature and all(Path(f).is_file() and (Path(f).stat().st_size or str(f).endswith('.covars')) for f in job['outputs']):
            print('RESUME '+job['name'],flush=True); return
        # Invalidate success before execution so failures cannot be mistaken for completion.
        stamp.unlink(missing_ok=True)
        print(shlex.join(cmd),flush=True)
        with open(Path(plan).parent/(job['name']+'.run.log'),'w') as log:
            result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:
            logpath=Path(plan).parent/(job['name']+'.run.log')
            print('\n'.join(logpath.read_text(errors='replace').splitlines()[-12:]),file=sys.stderr)
            raise RuntimeError(f'{job["name"]} failed (exit {result.returncode}); see {logpath}')
        if not all(Path(f).is_file() and (Path(f).stat().st_size or str(f).endswith('.covars')) for f in job['outputs']): raise RuntimeError('Missing output: '+job['name'])
        stamp.write_text(json.dumps({'signature':signature,'command':cmd}))

    def run_chromosome(chain):
        for job in chain:
            run_job(job)

    def flush_chromosomes(groups):
        if not groups:
            return
        # Each chromosome keeps its view/export/index/association dependency order.
        # Wait for every chromosome before merge-results or the next serial stage.
        with ThreadPoolExecutor(max_workers=spec.get('config',{}).get('jobs',1)) as pool:
            futures=[pool.submit(run_chromosome,chain) for chain in groups.values()]
            try:
                for future in futures:
                    future.result()
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        groups.clear()

    groups={}
    try:
        for job in spec['jobs']:
            match=re.match(r'^(chr(?:[0-9]+|X))\.',job['name'])
            if match:
                groups.setdefault(match[1],[]).append(job)
            else:
                flush_chromosomes(groups)
                run_job(job)
        flush_chromosomes(groups)
    finally:
        lock.close()

def save_plan(directory,jobs,a):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    plan=directory/'run.plan.json'
    spec={'config':vars(a),'replace':a.replace,'jobs':jobs}
    # Do not silently reuse results across different builds, subsets or model settings.
    if plan.exists():
        old=json.loads(plan.read_text()); oc=old.get('config',{}).copy(); nc=vars(a).copy()
        ignored=['run','replace','sparse_grm','saige_dir','saige_rscript','jobs']
        if nc.get('module') != 'prep_gwas':
            # Per-trait checkpoints validate commands and actual input files;
            # selecting more chromosomes/traits can safely reuse completed jobs.
            ignored+=['chromosomes','phenotypes']
        for key in ignored: oc.pop(key,None); nc.pop(key,None)
        if oc!=nc and not a.replace and any(directory.glob('*.done.json')): raise ValueError(f'Configuration changed at {directory}; choose a new output directory or --replace TRUE')
    plan.write_text(json.dumps(spec,indent=2)+'\n')
    lines=['#!/usr/bin/env bash','set -euo pipefail',f'# Commands and outputs: {plan}']
    lines.extend('# '+shlex.join(j['cmd']) for j in jobs)
    lines.append(shlex.join(['python3',SELF,'execute-plan',str(plan)]))
    script=directory/'run.cmd'; script.write_text('\n'.join(lines)+'\n')
    print(script,flush=True)
    if a.run: execute(plan)

def plan_node(args):
    name, covar_file, method = args[:3]
    sep = args.index('--', 3)
    out = args.index('--outputs', 3, sep)
    if args[3] != '--inputs': raise ValueError('Invalid plan-node arguments')
    job = dict(name=name, cmd=args[sep+1:], inputs=args[4:out], outputs=args[out+1:sep],
               covars=[covar_file, method] if covar_file else None)
    print(json.dumps(job))


def finalize_plan(args):
    from types import SimpleNamespace
    directory, jobs_file = args[:2]
    config = dict(item.split('=', 1) for item in args[2:])
    for key in ('threads', 'memory', 'min_mac', 'jobs'): config[key] = int(config[key])
    for key in ('run', 'replace', 'sparse_grm'): config[key] = config[key] == 'TRUE'
    for key in ('event_col', 'extract', 'keep'):
        if config[key] == '': config[key] = None
    config['chromosomes'] = config['chromosomes'].split(',')
    with open(jobs_file) as source: jobs = [json.loads(line) for line in source if line.strip()]
    save_plan(directory, jobs, SimpleNamespace(**config))


if __name__=='__main__':
    try:
        if len(sys.argv)>1 and sys.argv[1]=='execute-plan': execute(sys.argv[2])
        elif len(sys.argv)>1 and sys.argv[1]=='build-metadata': Path(sys.argv[2]+'.grch').write_text(sys.argv[3]+'\n')
        elif len(sys.argv)>1 and sys.argv[1]=='pgen-view':
            src,dst=sys.argv[2:]
            for ext in ('.pgen',):
                target=Path(dst+ext)
                if target.is_symlink(): target.unlink()
                elif target.exists(): raise FileExistsError(target)
                target.symlink_to(src+ext)
            with open(src+'.psam') as inp, open(dst+'.psam.tmp','w') as out:
                header=inp.readline()
                iid_only=header.startswith('#IID')
                if not iid_only and not header.startswith('#FID'): raise ValueError('Unsupported PSAM header')
                out.write(('#FID\t'+header.lstrip('#')) if iid_only else header)
                for line in inp: out.write(('0\t' if iid_only else '')+line)
            os.replace(dst+'.psam.tmp',dst+'.psam')
            with open(dst+'.pvar.tmp','wb') as out:
                if Path(src+'.pvar.zst').exists(): subprocess.run(['zstd','-dc',src+'.pvar.zst'],stdout=out,check=True)
                else:
                    with open(src+'.pvar','rb') as inp: shutil.copyfileobj(inp,out)
            os.replace(dst+'.pvar.tmp',dst+'.pvar')
        elif len(sys.argv)>1 and sys.argv[1]=='create-grm':
            actual,prefix=sys.argv[2:4]
            subprocess.run(sys.argv[4:],check=True)
            pairs=[(Path(actual+suffix),Path(prefix+'.sparseGRM.mtx'+suffix))
                   for suffix in ('','.sampleIDs.txt')]
            for src,dst in pairs:
                if not src.is_file() or not src.stat().st_size:
                    raise RuntimeError('Missing GRM output: '+str(src))
            for src,dst in pairs:
                os.replace(src,dst)
        elif len(sys.argv)>1 and sys.argv[1]=='grm-links':
            actual,prefix=sys.argv[2:]
            for src,dst in [(actual,prefix+'.sparseGRM.mtx'),(actual+'.sampleIDs.txt',prefix+'.sparseGRM.mtx.sampleIDs.txt')]:
                p=Path(dst)
                if p.is_symlink(): p.unlink()
                elif p.exists(): raise FileExistsError('Refusing to replace existing GRM: '+dst)
                p.symlink_to(src)
        elif len(sys.argv)>1 and sys.argv[1]=='plan-node': plan_node(sys.argv[2:])
        elif len(sys.argv)>1 and sys.argv[1]=='finalize-plan': finalize_plan(sys.argv[2:])
        else:
            sys.exit(subprocess.call(['bash', str(HERE.parent/'gwas.sh')]+sys.argv[1:]))
    except (ValueError,RuntimeError,FileNotFoundError,subprocess.CalledProcessError) as e:
        print('ERROR: '+str(e),file=sys.stderr); sys.exit(1)
