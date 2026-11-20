"""Flux stage resources; execution and dependency handling live in scripts/0f."""
from __future__ import annotations
import os
import sys
from common import ROOT, log, outdir, save_config, validate_resources
from workflow_jobs import Job, write_plan, run_plan
from workflow_layout import new_run_dir


def plan(cfg, args, machine):
    validate_resources(cfg['resources'])
    run_dir = new_run_dir(outdir(cfg, 'logs'))
    config_path = run_dir / 'config.yml'
    save_config(cfg, config_path)
    stages = ['prepare', 'train', 'readout', 'evaluate', 'baseline', 'report'] if args.stage == 'all' else [args.stage]
    if args.from_stage:
        stages = stages[stages.index(args.from_stage):]
    dependencies = {'prepare': [], 'train': ['prepare'], 'readout': ['train'],
                    'evaluate': ['readout'], 'baseline': ['readout'],
                    'report': ['evaluate', 'baseline']}
    resource = cfg['resources']
    jobs = []
    for stage in stages:
        gpu = stage in ('train', 'readout', 'rollout', 'ablate', 'demo') and args.device != 'cpu'
        device = ('gpu' if gpu else 'cpu') if machine == 'hpc' else args.device
        cores = min(int(resource['cores']), int(resource.get('gpu_cores', resource['cores']))) if gpu else int(resource['cores'])
        command = [sys.executable, str(ROOT / 'f/flux.py'), stage, '--worker',
                   '--machine', 'local', '--device', device, '--config', str(config_path), '--cores', str(cores)]
        if args.resume: command.append('--resume')
        if args.replace: command.append('--replace')
        threads = 1 if stage == 'prepare' else cores
        env = {key: str(threads) for key in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')}
        env.update(WORKFLOW_CORES=str(cores), PYTHONDONTWRITEBYTECODE='1')
        parents = [p for p in dependencies.get(stage, []) if p in stages] if args.stage == 'all' else []
        queue = str(resource['gpu_queue' if gpu else 'cpu_queue'])
        if queue == 'auto': queue = os.environ['HPC_CPU_QUEUE']
        jobs.append(Job(stage, command, depends=parents, job_name='flux.' + stage,
            resources=dict(queue=queue, cores=cores,
                memory_gib=float(resource['gpu_memory_gib' if gpu else 'cpu_memory_gib']),
                hours=float(resource['hours']),
                gpus=(int(resource['gpus']) if stage in ('train', 'ablate', 'demo') else 1) if gpu else 0),
            env=env, cwd=str(ROOT)))
    path = write_plan(jobs, run_dir, mode=machine, name='flux')
    result = run_plan(path, submit=machine, jobs=int(os.getenv('jobs_N', '1')),
                      wait=os.getenv('SUBMIT_WAIT') == 'TRUE', dry_run=args.plan_only)
    log(f'Plan: {path}')
    return result
