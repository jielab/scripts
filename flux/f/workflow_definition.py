"""LSF command files with stable names and the repository's ID-based dependencies."""
from pathlib import Path
import json
import os
import re
import shlex
import sys

from common import OUTPUT_FORMAT, ROOT, outdir
from lsf_submit import submit_graph
from submission_io import atomic_json, receipt_lock

STAGES = ('prepare', 'train', 'readout', 'evaluate', 'baseline', 'dynamics', 'blackout', 'rollout')
PARENTS = {'train': ['prepare'], 'readout': ['train'], 'evaluate': ['readout'],
           'baseline': ['readout'], 'evaluate-baseline': ['readout'], 'dynamics': ['train'], 'blackout': ['readout'], 'rollout': ['readout'],
           'ablate': ['train']}
GPU = {'dynamics', 'train', 'readout', 'blackout', 'rollout', 'ablate', 'demo'}


def selected(args):
    if args.stage != 'run':
        return [args.stage]
    stages = list(STAGES)
    if not args.no_ablation:
        stages.append('ablate')
    return stages


def require_fresh_run(cfg):
    root = Path(cfg['paths']['out'])
    # A --plan-only invocation may already have created command files and log directories.
    occupied = any((root / name).exists() for name in ('main', 'ablation', 'paper/main', 'jobs/analysis.yml'))
    if occupied or any((root / 'jobs').glob('*.json')):
        raise ValueError(f'Analysis outputs already exist in {root}. Start ./flux.sh run in a fresh analysis directory.')


def plan(cfg, args, mode):
    stages = selected(args)
    # One CPU allocation for the two endpoint-scoring passes. Keep the heavy
    # GPU analyses independent; neither CPU scoring pass gates them.
    cpu_pair = args.stage == 'run' and 'evaluate' in stages and 'baseline' in stages
    if cpu_pair:
        stages = ['evaluate-baseline' if s == 'evaluate' else s for s in stages if s != 'baseline']
    tasks = []
    for stage in stages:
        gpu = stage in GPU and args.device != 'cpu'
        r = cfg['resources']
        cores = min(r['gpu_cores'], r['cores']) if gpu else r['cores']
        queue = str(r['gpu_queue'] if gpu else r['cpu_queue'])
        memory_gib = None
        if stage == 'prepare':
            cores = min(cores, int(r.get('prepare_cores', 8)))
            queue = str(r.get('prepare_queue', queue))
            memory_gib = int(r.get('prepare_memory_gib', 32))
            if cores < 1 or memory_gib < 1:
                raise ValueError('prepare cores and memory must be positive')
        if queue == 'auto':
            queue = os.environ.get('HPC_CPU_QUEUE', '33')
        if not re.fullmatch(r'[\w.-]+', queue):
            raise ValueError('Invalid queue name')
        prefix = 'paper/' if stage == 'paper' else ''
        logs = outdir(cfg, prefix + 'logs')
        jobs = outdir(cfg, prefix + 'jobs')
        name = 'flux.' + stage
        device = ('gpu' if gpu else 'cpu') if mode == 'hpc' else args.device
        command = [sys.executable, str(ROOT/'f/flux.py'), stage, '--worker',
                   '--device', device, '--machine', 'local', '--cores', str(cores),
                   '--settings-json', json.dumps(cfg, separators=(',', ':'))]
        if args.replace:
            command.append('--replace')
        if args.resume and stage == 'train':
            command.append('--resume')
        lines = ['#!/usr/bin/env bash', f'#BSUB -J {name}', f'#BSUB -q {queue}', f'#BSUB -n {cores}',
                 '#BSUB -R "span[hosts=1]"', f'#BSUB -o "{logs}/{stage}.%J.out"',
                 f'#BSUB -e "{logs}/{stage}.%J.err"']
        if memory_gib is not None:
            lines[4] = f'#BSUB -R "span[hosts=1] rusage[mem={memory_gib}GB/host]"'
        if gpu:
            n = r['gpus'] if stage in ('train','ablate','demo') else 1
            lines.append(f'#BSUB -gpu "num={n}"')
        threads = 1 if stage == 'prepare' else cores
        lines += ['set -euo pipefail', 'cd ' + shlex.quote(str(ROOT)),
                  f'export WORKFLOW_CORES={cores} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads}',
                  'export PYTHONDONTWRITEBYTECODE=1', 'exec ' + shlex.join(command)]
        path = jobs / (stage + '.cmd')
        tmp = path.with_suffix('.tmp')
        tmp.write_text('\n'.join(lines) + '\n')
        tmp.replace(path)
        parents = [s for s in PARENTS.get(stage, []) if s in stages]
        tasks.append(dict(key=stage,path=path,depends=[(s,'done') for s in parents]))
        print(f'{name}: {path}  after {", ".join(parents) or "—"}', flush=True)
    return tasks


def schedule(cfg, args, mode):
    directory = outdir(cfg, 'paper/jobs' if args.stage == 'paper' else 'jobs')
    receipt = directory / (args.stage + '.json')
    # The shell scripts get read by bsub at submission. The resolved scientific
    # config is inside each script, so later plan-only calls cannot change it.
    with receipt_lock(outdir(cfg,'jobs') / 'dispatch.json'):
        saved = json.loads(receipt.read_text()) if receipt.is_file() else {}
        if any(row.get('output_format') != OUTPUT_FORMAT for row in saved.values()):
            raise ValueError(f'{receipt} belongs to a different flux release. Use a fresh analysis directory.')
        if args.stage == 'run' and not saved:
            require_fresh_run(cfg)
        tasks = plan(cfg, args, mode)
        if args.plan_only:
            return
        # Recover a partially acknowledged chain before starting a new submission.
        import hashlib
        for task in tasks:
            task['sha256'] = hashlib.sha256(task['path'].read_bytes()).hexdigest()
        confirmed = {k: row for k, row in saved.items() if row.get('status') == 'submitted'}
        if args.stage == 'run':
            planned = {task['key']: task['sha256'] for task in tasks}
            if any(k not in planned or planned[k] != row.get('command_sha256') for k,row in saved.items()):
                raise ValueError('A different analysis was submitted in this directory. Use a fresh analysis directory.')
            if len(confirmed) == len(tasks):
                print('Analysis already submitted. See ' + str(receipt))
                return
        elif confirmed:
            import subprocess
            ids = [row['job_id'] for row in confirmed.values()]
            check = subprocess.run(['bjobs','-a','-noheader','-o','stat', *map(str,ids)],
                                   text=True, capture_output=True)
            if set(check.stdout.split()) & {'PEND','RUN','WAIT','PSUSP','USUSP','SSUSP','UNKWN','ZOMBI'}:
                print('This workflow already has active jobs. See ' + str(receipt))
                return
            if check.returncode and 'not found' not in check.stderr.lower():
                raise RuntimeError('Cannot verify previous job status: ' + check.stderr)
            saved = {}
        def persist(key, row):
            saved[key] = dict(row, output_format=OUTPUT_FORMAT)
            atomic_json(receipt, saved)
        os.environ.pop('WORKFLOW_JOB_NAME', None)  # Keep the distinct per-stage LSF names.
        submit_graph(tasks, saved, persist=persist, receipt=receipt)
