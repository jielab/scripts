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

from common import ROOT, load_config, log, merge, outdir, save_config, validate_resources

STAGES = ('prepare', 'train', 'readout', 'evaluate', 'baseline', 'report')


def parser():
    p = argparse.ArgumentParser(description='Independent HealthFlux v3 reproduction on UKB.',
        epilog='Examples: bash flux.sh demo --machine local --device cpu\n'
               '          bash flux.sh all --machine hpc --device gpu --plan-only',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('stage', nargs='?', default='all',
        choices=(*STAGES, 'all', 'demo', 'rollout', 'ablate', 'test'))
    p.add_argument('--config', help='YAML overrides [config/flux.yml]')
    p.add_argument('--machine', choices=('auto','local','hpc'), default='hpc')
    p.add_argument('--device', choices=('auto','cpu','gpu'), default='gpu')
    p.add_argument('--dir-ukb', help='UKB phenotype root, containing rap/')
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
        path = outdir(cfg, 'logs') / 'train.resolved.yml'
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
    elif stage == 'evaluate':
        from evaluate import evaluate
        evaluate(cfg)
    elif stage == 'baseline':
        from evaluate import run_baselines
        run_baselines(cfg, args.device, replace=args.replace)
    elif stage == 'rollout':
        from rollout import rollout
        rollout(cfg, args.device, replace=args.replace)
    elif stage == 'report':
        from report import report
        report(cfg)
    elif stage == 'ablate':
        ablate(cfg, args)
    else:
        raise ValueError(stage)


def ablate(cfg, args):
    # Retrain each ablation; changing flags on an already fitted model is invalid.
    from common import torch_load
    from train import resume_contract
    reference=torch_load(Path(cfg['paths']['out'])/'model/flux.pt')
    base_cfg=reference['config']
    old_contract=resume_contract(base_cfg); new_contract=resume_contract(cfg)
    old_contract.pop('cores',None); new_contract.pop('cores',None)
    if (old_contract!=new_contract or cfg['data']!=base_cfg['data'] or cfg['model']!=base_cfg['model']
        or cfg['readout']!=base_cfg['readout'] or
        any(cfg['train'][k]!=base_cfg['train'][k] for k in ('max_steps','epochs','patience')) or
        int(cfg['resources']['gpus'])!=int(reference.get('world_size',1))):
        raise ValueError('Ablations must use the fitted main model settings and GPU count. Pass its original --config (or logs/resolved.yml); inspect model checkpoint config for the authoritative settings.')
    source = Path(cfg['paths']['out'])/'data'
    if not (source/'cohort.sqlite').exists(): raise FileNotFoundError('Run prepare first')
    for name in ('no_ode','no_update'):
        variant = copy.deepcopy(cfg)
        variant['paths']['out'] = str(Path(cfg['paths']['out'])/'ablation'/name)
        variant['model'][name] = True
        target = outdir(variant, 'data')
        for filename in ('cohort.sqlite', 'meta.json'):
            link = target/filename
            if not link.exists(): link.symlink_to((source/filename).resolve())
        save_config(variant, outdir(variant, 'logs')/'resolved.yml')
        for stage in ('train','readout','evaluate','report'): execute(stage, variant, args)
    execute('report', cfg, args)


def main():
    if any(arg.split('=', 1)[0] == '--end-date' for arg in sys.argv[1:]):
        parser().error('--end-date is no longer supported; edit 0f/phenotype.R::date_follow_end')
    args = parser().parse_args()
    if args.stage == 'test':
        return subprocess.call([sys.executable, '-m', 'unittest', 'discover', '-s', str(ROOT/'f/tests'), '-v'])
    if args.stage == 'demo':
        import yaml
        overrides = None
        if args.config:
            with open(args.config) as handle:
                overrides = yaml.safe_load(handle) or {}
        cfg = load_config(ROOT/'config/demo.yml', overrides=overrides)
        if not isinstance(overrides, dict) or 'out' not in overrides.get('paths', {}):
            cfg['paths']['out'] = str(ROOT/'demo_out')
    else:
        cfg = load_config(args.config)
    for flag, key in ((args.dir_ukb,'ukb'), (args.dir_out,'out')):
        if flag: cfg['paths'][key] = str(Path(flag).expanduser().resolve())
    if args.cores is not None: cfg['resources']['cores'] = args.cores
    if args.gpus is not None: cfg['resources']['gpus'] = args.gpus
    validate_resources(cfg['resources'])
    if int(cfg['resources']['gpus']) > 1 and args.device == 'cpu' and args.stage in ('train','all','demo','ablate'):
        raise ValueError('--gpus > 1 requires --device gpu or auto for training')
    machine = args.machine if args.machine != 'auto' else os.getenv('FLUX_MACHINE', 'local')
    # Scheduling is intentionally done before torch imports or cohort reads.
    if args.plan_only or (machine == 'hpc' and not args.worker and not os.getenv('LSB_JOBID')):
        from submit import plan
        plan(cfg, args, machine)
        return
    os.environ['OMP_NUM_THREADS'] = str(cfg['resources']['cores'])
    os.environ['MKL_NUM_THREADS'] = str(cfg['resources']['cores'])
    if args.stage == 'demo':
        from demo import make_demo
        cfg = make_demo(cfg)
    if int(os.getenv('RANK','0')) == 0:
        save_config(cfg, outdir(cfg, 'logs')/f'{args.stage}.resolved.yml')
    if args.stage in ('all','demo'):
        stages = list(STAGES)
        if args.stage == 'demo': stages.insert(-1, 'rollout')
        for stage in stages: execute(stage, cfg, args)
    else:
        execute(args.stage, cfg, args)
    if int(os.getenv('RANK','0')) == 0 and args.stage in ('all','demo','train','prepare'):
        save_config(cfg,outdir(cfg,'logs')/'resolved.yml')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main() or 0)
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(2)
