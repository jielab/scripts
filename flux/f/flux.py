#!/usr/bin/env python3
"""Public dispatcher. Imports heavy numerical modules only inside workers."""
from __future__ import annotations
import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

from common import ROOT, load_config, log, outdir, save_config, validate_resources

from workflow_mode import resolve_mode
from workflow_paths import configured_path

from submission_io import receipt_lock
from workflow_definition import STAGES, require_fresh_run, selected


def parser():
    p = argparse.ArgumentParser(description='Independent HealthFlux v3 YGB analysis')
    p.add_argument('stage', choices=(*STAGES, 'run', 'demo', 'ablate', 'paper', 'evaluate-baseline'))
    p.add_argument('--no-ablation', action='store_true', help='Skip the two architecture ablations in run')
    p.add_argument('--endpoint-file', help='Prespecified single-column ICD-10/DEATH panel')
    p.add_argument('--settings-json', help=argparse.SUPPRESS)
    p.add_argument('--config', help='YAML overrides [config/flux.yml]')
    p.add_argument('--machine', choices=('auto','local','hpc'), default='auto')
    p.add_argument('--device', choices=('auto','cpu','gpu'), default='auto')
    p.add_argument('--dir-ygb', help='YGB phenotype root, containing rap/')
    p.add_argument('--dir-out', help='One analysis output directory')
    p.add_argument('--cores', type=int, help='CPU threads per worker')
    p.add_argument('--gpus', type=int, help='GPUs for distributed world-model training [1]')
    p.add_argument('--resume', action='store_true', help='Continue training from flux.last.pt')
    p.add_argument('--replace', action='store_true', help='Rebuild the selected stage outputs')
    p.add_argument('--plan-only', action='store_true', help='Write inspectable commands without analysis/submission')
    p.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    return p


def launch_train(cfg, args):
    gpus = int(cfg['resources']['gpus'])
    if gpus > 1 and 'RANK' not in os.environ:
        if args.device == 'cpu':
            raise ValueError('--gpus > 1 requires --device gpu or auto')
        from common import choose_device
        import torch
        device = choose_device(args.device)
        if device.type != 'cuda' or torch.cuda.device_count() < gpus:
            raise ValueError(f'{gpus} GPUs requested; only {torch.cuda.device_count()} visible')
        path = outdir(cfg, 'logs/main') / 'train.resolved.yml'
        save_config(cfg, path)
        cmd = [sys.executable, '-m', 'torch.distributed.run', '--standalone',
               f'--nproc_per_node={gpus}', str(ROOT/'f/flux.py'), 'train',
               '--worker', '--machine', 'local', '--device', 'gpu', '--config', str(path)]
        if args.resume: cmd.append('--resume')
        if args.replace: cmd.append('--replace')
        subprocess.run(cmd, check=True)
    else:
        from train import train_world
        train_world(cfg, args.device, resume=args.resume, replace=args.replace)


def execute(stage, cfg, args):
    log(f'{stage}: {cfg["paths"]["out"]}')
    if stage == 'prepare':
        from prepare import prepare_cohort
        prepare_cohort(cfg, replace=args.replace)
    elif stage == 'train':
        launch_train(cfg, args)
    elif stage == 'readout':
        from train import fit_readout
        fit_readout(cfg, args.device, replace=args.replace)
    elif stage == 'evaluate-baseline':
        execute('evaluate', cfg, args)
        execute('baseline', cfg, args)
    elif stage == 'evaluate':
        from evaluate import evaluate
        evaluate(cfg)
    elif stage == 'baseline':
        from evaluate import run_baselines
        run_baselines(cfg, args.device, replace=args.replace)
    elif stage == 'rollout':
        from rollout import rollout
        rollout(cfg, args.device, replace=args.replace)
    elif stage == 'dynamics':
        from event_model import heldout_dynamics
        heldout_dynamics(cfg,args.device)
    elif stage == 'blackout':
        from blackout import blackout
        blackout(cfg, args.device)
    elif stage == 'paper':
        from paper import paper
        paper(cfg)
    elif stage == 'ablate':
        ablate(cfg, args)
    else:
        raise ValueError(stage)


def ablate(cfg, args):
    # Retrain each ablation; changing flags on an already fitted model is invalid.
    from common import torch_load
    from train import resume_contract
    reference=torch_load(Path(cfg['paths']['out'])/'main/results/models/flux.pt')
    base_cfg=reference['config']
    old_contract=resume_contract(base_cfg); new_contract=resume_contract(cfg)
    old_contract.pop('cores',None); new_contract.pop('cores',None)
    if (old_contract!=new_contract or cfg['data']!=base_cfg['data'] or cfg['model']!=base_cfg['model']
        or cfg['readout']!=base_cfg['readout'] or
        any(cfg['train'][k]!=base_cfg['train'][k] for k in ('max_steps','epochs','patience')) or
        int(cfg['resources']['gpus'])!=int(reference['world_size'])):
        raise ValueError('Ablations must use the fitted main model settings and GPU count. Pass its original --config (or logs/resolved.yml); inspect model checkpoint config for the authoritative settings.')
    source = Path(cfg['paths']['out'])/'main/results/data'
    if not (source/'cohort.sqlite').exists(): raise FileNotFoundError('Run prepare first')
    for name in ('no_ode','no_update'):
        variant = copy.deepcopy(cfg)
        variant['paths']['out'] = str(Path(cfg['paths']['out'])/'ablation'/name)
        variant['model'][name] = True
        target = outdir(variant, 'main/results/data')
        for filename in ('cohort.sqlite', 'meta.json'):
            link = target/filename
            if not link.exists(): link.symlink_to((source/filename).resolve())
        save_config(variant, outdir(variant, 'logs/main')/'resolved.yml')
        for stage in ('train','readout','evaluate','dynamics'): execute(stage, variant, args)
        from train import state_cache_key
        (outdir(variant,'main/results')/'ablation.json').write_text(json.dumps(dict(parent_risk_key=state_cache_key(cfg))))


def main():
    if any(arg.split('=', 1)[0] == '--end-date' for arg in sys.argv[1:]):
        parser().error('--end-date is no longer supported; edit 0f/phenotype.R::date_follow_end')
    args = parser().parse_args()
    if args.settings_json:
        cfg = load_config(overrides=json.loads(args.settings_json))
    elif args.stage == 'demo':
        import yaml
        overrides = None
        if args.config:
            with open(args.config) as handle:
                overrides = yaml.safe_load(handle) or {}
        cfg = load_config(ROOT/'config/demo.yml', overrides=overrides)
        if not isinstance(overrides, dict) or 'out' not in overrides.get('paths', {}):
            cfg['paths']['out'] = configured_path('ANALYSIS_ROOT', 'flux/demo')
    else:
        cfg = load_config(args.config)
    for flag, key in ((args.dir_ygb,'ygb'), (args.dir_out,'out')):
        if flag: cfg['paths'][key] = str(Path(flag).expanduser().resolve())
    if args.no_ablation and args.stage != 'run':
        raise ValueError('--no-ablation applies to run only')
    if args.resume and args.stage != 'train':
        raise ValueError('--resume applies only to an interrupted train stage of this release')
    if args.stage == 'run' and args.replace:
        raise ValueError('run starts in a fresh analysis directory. Use --replace with an individual stage.')
    if args.endpoint_file:
        import re
        lines = Path(args.endpoint_file).read_text(encoding='utf-8-sig').splitlines()
        codes = [v.strip().upper() for v in lines if v.strip() and not v.lstrip().startswith('#')]
        if codes and codes[0] in {'ENDPOINT','CODE','ICD10'}:
            codes.pop(0)
        if not codes or len(set(codes)) != len(codes) or any(not re.fullmatch(r'[A-Z][0-9]{2}|DEATH', v) for v in codes):
            raise ValueError('Panel must contain distinct three-character ICD-10 codes or DEATH, one per line')
        cfg['data']['endpoints'] = codes
    if args.stage in ('paper','ablate') and not (args.settings_json or args.config or args.endpoint_file):
        source = Path(cfg['paths']['out']) / 'jobs' / 'analysis.yml'
        if source.exists():
            cfg = load_config(source)
            if args.dir_out:
                cfg['paths']['out'] = str(Path(args.dir_out).resolve())
    if args.cores is not None: cfg['resources']['cores'] = args.cores
    elif args.stage == 'paper': cfg['resources']['cores'] = min(2, int(cfg['resources']['cores']))
    if args.gpus is not None: cfg['resources']['gpus'] = args.gpus
    validate_resources(cfg['resources'])
    total = cfg['train'].get('global_batch_size')
    if total is not None:
        micro = int(cfg['resources']['gpus'])*int(cfg['train']['batch_size'])
        if micro <= 0 or total < micro or total % micro:
            raise ValueError('global_batch_size must be divisible by GPUs times per-GPU batch_size')
        cfg['train']['accumulation_steps'] = int(total // micro)
    os.environ['WORKFLOW_CORES'] = str(cfg['resources']['cores'])
    if int(cfg['resources']['gpus']) > 1 and args.device == 'cpu' and args.stage in ('train','run','demo','ablate'):
        raise ValueError('--gpus > 1 requires --device gpu or auto for training')
    machine = args.machine if args.machine != 'auto' else resolve_mode(os.getenv('FLUX_MACHINE', 'auto'))
    # Scheduling is intentionally done before torch imports or cohort reads.
    if args.plan_only or (machine == 'hpc' and not args.worker and not os.getenv('LSB_JOBID')):
        from workflow_definition import schedule
        schedule(cfg, args, machine)
        return
    os.environ['OMP_NUM_THREADS'] = str(cfg['resources']['cores'])
    os.environ['MKL_NUM_THREADS'] = str(cfg['resources']['cores'])
    if args.stage == 'demo':
        from demo import make_demo
        cfg = make_demo(cfg)
    if args.stage == 'run':
        # Reserve the local analysis root before opening scientific outputs.
        with receipt_lock(outdir(cfg, 'jobs') / 'dispatch.json'):
            require_fresh_run(cfg)
            save_config(cfg, outdir(cfg, 'jobs') / 'analysis.yml')
    elif args.stage in ('prepare','demo') and int(os.getenv('RANK','0')) == 0:
        save_config(cfg, outdir(cfg,'jobs') / 'analysis.yml')
    if int(os.getenv('RANK','0')) == 0:
        stage_logs = 'paper/logs' if args.stage == 'paper' else 'logs/main'
        save_config(cfg, outdir(cfg, stage_logs)/f'{args.stage}.resolved.yml')
    stages = list(STAGES) if args.stage == 'demo' else selected(args)
    for stage in stages:
        execute(stage, cfg, args)
    if int(os.getenv('RANK','0')) == 0 and args.stage in ('run','demo','train','prepare'):
        save_config(cfg,outdir(cfg,'logs/main')/'resolved.yml')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main() or 0)
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(2)
