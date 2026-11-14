#!/usr/bin/env python3
"""Diet entry-specific job names with the shared local/LSF plan executor."""
import argparse
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(os.getenv('DIET_HELPERS', str(ROOT.parent / '0f')))
from project_job import scientific_budget
from workflow_jobs import Job, resolve_mode, run_plan, write_plan


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--entry', choices=['disco', 'diet'], required=True)
    p.add_argument('arguments', nargs=argparse.REMAINDER)
    a = p.parse_args(argv)
    arguments = a.arguments[1:] if a.arguments[:1] == ['--'] else a.arguments
    resources = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    resources.add_argument('--plan-only', '--write-cmd', action='store_true')
    resources.add_argument('--job-cores', type=int)
    resources.add_argument('--job-memory-gb', type=float)
    resources.add_argument('--job-hours', type=float, default=120)
    resources.add_argument('--job-gpus', type=int, default=0)
    resources.add_argument('--queue')
    resources.add_argument('--dir-jobs')
    opts, scientific = resources.parse_known_args(arguments)
    # Reuse the scientific parser so typos are rejected before submitting.
    from pipeline import parser, Workflow, validate_args
    scientific_parser = parser()
    worker_args = scientific_parser.parse_args(['--entry', a.entry, *scientific])
    validate_args(scientific_parser, worker_args)
    out = Workflow(worker_args).out
    cores, memory = scientific_budget('diet', scientific)
    if opts.job_cores is None:
        opts.job_cores = cores
    elif opts.job_cores < cores:
        raise ValueError(f'--job-cores must cover the scientific worker budget ({cores})')
    if opts.job_memory_gb is None:
        opts.job_memory_gb = max(float(os.getenv('MEMORY_DEFAULT_GIB', '32')), memory)
    elif opts.job_memory_gb < memory:
        raise ValueError('--job-memory-gb is smaller than the scientific memory budget')
    mode = resolve_mode()
    entry = ROOT / ('1.disco.sh' if a.entry == 'disco' else 'diet.sh')
    env = {k: v for k, v in os.environ.items() if k.startswith('DIET_')}
    env.update(PROJECT_JOB_ACTIVE='1', WORKFLOW_RESUME_POLICY='new_plan',
               WORKFLOW_DEVICE='auto' if mode == 'local' else ('cuda' if opts.job_gpus else 'cpu'))
    queue = opts.queue or os.getenv('HPC_GPU_QUEUE' if opts.job_gpus else 'HPC_CPU_QUEUE',
                                   '4a100-40' if opts.job_gpus else '33')
    job_name = 'diet.disco' if a.entry == 'disco' else 'diet'
    job = Job(a.entry, ['bash', str(entry), *scientific], job_name=job_name,
              resources=dict(queue=queue, cores=opts.job_cores, memory_gib=opts.job_memory_gb,
                             hours=opts.job_hours, gpus=opts.job_gpus), env=env, cwd=os.getcwd())
    plan = write_plan([job], opts.dir_jobs or out / 'logs', mode=mode, name='diet')
    result = run_plan(plan, submit=mode, jobs=1, wait=os.getenv('SUBMIT_WAIT') == 'TRUE',
                      dry_run=opts.plan_only)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, RuntimeError) as exc:
        sys.exit(f'ERROR: {exc}')
