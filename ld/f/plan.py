#!/usr/bin/env python3
"""Create independent ancestry/chromosome CPU commands and their exact merge set."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import fcntl
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from ld_pipeline import Pipeline, POPS, write_json
from resources import cpu_slots
import submit

ROOT = Path(__file__).resolve().parents[1]
ENV = ('SUBMIT GRCH INPUT_FORMAT PFILE_DIR VCF_DIR PANEL THREADS MEMORY_MB LD_MEMORY_GIB LD_HOURS '
       'LD_CXX LD_ZSTD_PREFIX LD_INPUT_MANIFEST D_GRID_KB MAX_D_KB STABLE_CONFIRMATIONS '
       'REFERENCE_COUNT_POLICY REFERENCE_MAX_REL_DEVIATION COINCIDENT_POLICY LD_STORAGE '
       'MIN_FREE_GB LD_TARGETS LD_TMP_ROOT LD_REPLACE').split()


def make_plan(preview=False):
    if sys.version_info < (3,11):
        raise ValueError('Python >=3.11 is required')
    pipeline = Pipeline()
    root = pipeline.out
    mode = os.environ['SUBMIT']
    cores, memory, hours = (int(os.environ[x]) for x in ('LD_SLOTS','LD_MEMORY_GIB','LD_HOURS'))
    queue, merge_cores = os.environ['LD_QUEUE'], 1
    if mode == 'hpc':
        queue, cores = cpu_slots(queue,cores,memory,hours)
        _, merge_cores = cpu_slots(queue,1,4,1)
    if not preview:
        subprocess.run(['bash',str(ROOT/'f/check.sh')],check=True)
    run = root/'logs'/('run-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')+f'-{os.getpid()}')
    run.mkdir(parents=True)
    (run/'commands').mkdir();(run/'logs').mkdir();(run/'receipts').mkdir()
    shutil.copytree(ROOT/'f',run/'code/f',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    code=run/'code/f';plan_path=run/'plan.json';run_id=run.name
    def command(label, output, population=None, chrom=None):
        is_merge = population is None
        slots, ram, duration = (merge_cores,4,1) if is_merge else (cores,memory,hours)
        lines=['#!/usr/bin/env bash']
        if mode=='hpc':
            lines += [f'#BSUB -q {queue}',f'#BSUB -n {slots}',f'#BSUB -W {duration}:00',
                      f'#BSUB -R "span[hosts=1] rusage[mem={ram}GB/host]"',f'#BSUB -J ld{pipeline.genome_build}.{label}',
                      f'#BSUB -oo "{run}/logs/{label}.LOG"',f'#BSUB -eo "{run}/logs/{label}.ERR"']
        lines += ['set -euo pipefail','cd '+shlex.quote(str(run/'code'))]
        for name in ENV:
            if name in os.environ:lines.append('export '+name+'='+shlex.quote(os.environ[name]))
        frozen = dict(LD_PYTHON=sys.executable,LD_WORKER='1',LD_PLAN_FILE=str(plan_path),LD_PLAN_ID=run_id,
                      OUT=str(output),CHRS=','.join(map(str,pipeline.chroms)) if is_merge else str(chrom),
                      LD_TMP_DIR=str(Path(os.environ['LD_TMP_ROOT'])/'parts'/population/f'chr{chrom}'/'tmp') if not is_merge else str(root/'tmp/merge'),
                      PLINK2=shutil.which(os.environ['PLINK2']) or os.environ['PLINK2'],
                      LD_CXX=shutil.which(os.environ['LD_CXX']) or os.environ['LD_CXX'])
        if population:frozen['LD_POP']=population
        for name,value in frozen.items():lines.append('export '+name+'='+shlex.quote(value))
        lines += ['export CUDA_VISIBLE_DEVICES= WORKFLOW_DEVICE=cpu', 'source '+shlex.quote(str(code/'env.sh'))]
        if mode=='hpc':
            lines += [f'if (( ${{LSB_DJOB_NUMPROC:-0}} < {slots} )); then echo "Insufficient CPU allocation" >&2; exit 2; fi']
        else:
            lines += ['exec >'+shlex.quote(str(run/'logs'/f'{label}.LOG'))+' 2>'+shlex.quote(str(run/'logs'/f'{label}.ERR'))]
        lines += ['exec "$PYTHON" '+shlex.quote(str(code/('merge.py' if is_merge else 'ld_pipeline.py')))]
        path=run/'commands'/f'{label}.cmd';path.write_text('\n'.join(lines)+'\n');path.chmod(0o700)
        return str(path)
    jobs=[]
    for chrom in pipeline.chroms:
        for pop in POPS:
            label=f'{pop}.chr{chrom}';out=root/'parts'/pop/f'chr{chrom}'
            path=command(label,out,pop,chrom)
            jobs.append(dict(key=label,population=pop,chromosome=chrom,output=str(out),command=path,
                             resources=dict(cores=cores,memory_gib=memory,hours=hours,gpus=0),threads=pipeline.threads))
    merge=command('LD.merge',root)
    plan=dict(schema='ld.plan.1',run_id=run_id,run_dir=str(run),output_root=str(root),mode=mode,
              genome_build=pipeline.genome_build,chromosomes=pipeline.chroms,populations=list(POPS),
              grid_kb=pipeline.grid,confirmations=pipeline.confirmations,replace=os.environ['LD_REPLACE'],
              jobs=jobs,merge=merge)
    write_json(plan_path,plan)
    print(f'[LD] Plan: {plan_path}',flush=True)
    print(f'[LD] {len(jobs)} independent CPU workers; threads={pipeline.threads}; D grid={pipeline.grid}; then joint merge',flush=True)
    if not preview:
        state=root/'logs/submission-state.json'
        with state.with_suffix('.lock').open('a') as handle:
            fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
            submit.run([j['command'] for j in jobs],mode,int(os.environ['LD_JOBS']),
                       os.environ['SUBMIT_WAIT']=='TRUE',merge,state=state,plan=plan_path)
    return plan_path


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan-only',choices=('TRUE','FALSE'),default='FALSE')
    a=p.parse_args()
    try:make_plan(a.plan_only=='TRUE')
    except (ValueError,RuntimeError,OSError,subprocess.SubprocessError) as exc:p.exit(2,f'ERROR: {exc}\n')
