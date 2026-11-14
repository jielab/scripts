"""LSF/local plan writer. No torch, cohort reads, or numerical work here."""
from __future__ import annotations
import importlib.util
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from uuid import uuid4
from common import ROOT, log, outdir, save_config, validate_resources


def plan(cfg, args, machine):
    validate_resources(cfg['resources'])
    logdir = outdir(cfg, 'logs')
    plan_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ') + '_' + uuid4().hex[:8]
    config_path = logdir/f'plan_{plan_id}.yml'
    save_config(cfg, config_path)
    stages = ['prepare','train','readout','evaluate','baseline','report'] if args.stage == 'all' else [args.stage]
    resource = cfg['resources']
    hpc = machine == 'hpc'
    previous = None
    index = []
    # Repository-owned QiMing core/RAM policy, if this is installed in scripts/.
    cpu_slots = None
    helper = ROOT.parent/'0f/qiming_resources.py'
    if helper.exists():
        spec = importlib.util.spec_from_file_location('qiming_resources', helper)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        cpu_slots = mod.cpu_slots
    if hpc and not args.plan_only and not shutil.which('bsub'):
        raise RuntimeError('bsub is not on PATH. Load LSF, or use --plan-only to inspect the commands.')
    for stage in stages:
        gpu = stage in ('train','readout','rollout','ablate','demo') and args.device != 'cpu'
        stage_device = ('gpu' if gpu else 'cpu') if hpc else args.device
        cmd = [sys.executable, str(ROOT/'f/flux.py'), stage, '--worker', '--machine','local',
               '--device',stage_device, '--config',str(config_path)]
        if args.resume: cmd.append('--resume')
        if args.replace: cmd.append('--replace')
        commandfile = logdir/f'{stage}.cmd'
        cores = int(resource['cores'])
        memory = float(resource['gpu_memory_gib' if gpu else 'cpu_memory_gib'])
        hours = float(resource['hours'])
        queue = str(resource['gpu_queue' if gpu else 'cpu_queue'])
        if hpc and not gpu and cpu_slots: queue, slots = cpu_slots(queue, cores, memory, hours)
        else: slots = cores
        if queue == 'auto': queue = os.getenv('HPC_CPU_QUEUE','33')
        directives = []
        if hpc:
            wall_hours, wall_minutes = divmod(math.ceil(hours * 60), 60)
            directives = [f'#BSUB -J flux.{stage}', f'#BSUB -q {queue}', f'#BSUB -n {slots}',
                f'#BSUB -W {wall_hours}:{wall_minutes:02d}',
                f'#BSUB -R "span[hosts=1] rusage[mem={memory:g}GB/host]"',
                f'#BSUB -o {shlex.quote(str(logdir/(stage+".%J.out")))}',
                f'#BSUB -e {shlex.quote(str(logdir/(stage+".%J.err")))}']
            if gpu:
                count = int(resource['gpus']) if stage in ('train','ablate','demo') else 1
                directives.append(f'#BSUB -gpu "num={count}:mode=exclusive_process"')
        script = '\n'.join(['#!/usr/bin/env bash', *directives, 'set -euo pipefail',
            f'export OMP_NUM_THREADS={cores} MKL_NUM_THREADS={cores}',
            'export PYTHONDONTWRITEBYTECODE=1', shlex.join(cmd), ''])
        commandfile.write_text(script)
        commandfile.chmod(0o755)
        if hpc and not args.plan_only:
            bsub = ['bsub'] + (['-w',f'done({previous})'] if previous else [])
            with commandfile.open() as handle:
                result = subprocess.run(bsub, stdin=handle, text=True, capture_output=True, check=True)
            found = re.search(r'Job <(\d+)>', result.stdout)
            if not found: raise RuntimeError(f'Cannot parse LSF submission: {result.stdout} {result.stderr}')
            previous = found.group(1)
            index.append(f'{stage}\t{previous}\t{commandfile}')
            log(f'Submitted {stage}: job {previous}')
        else:
            index.append(f'{stage}\tplan\t{commandfile}')
            log(f'Plan: {commandfile}')
    (logdir/'jobs.tsv').write_text('stage\tjob_id\tcommand\n'+'\n'.join(index)+'\n')
    if args.plan_only:
        log('No analysis or submission performed. Stage order/dependencies are in jobs.tsv; rerun without --plan-only to execute/submit.')
