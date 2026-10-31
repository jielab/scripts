"""Serial, independent PRSformer traits sharing one frozen input cache."""
import json
import os
from pathlib import Path
import subprocess
import sys
from resume import run_lock, sha


def publication_valid(base, evaluation):
    marker = base / 'PUBLISHED.json'
    if not marker.is_file(): return False
    receipt = json.loads(marker.read_text())
    files = receipt.get('files', {})
    expected = [base/name for name in ('model.pt','split.rds','performance.png','performance.xlsx','training.png','training.xlsx','predictions.rds')]
    if not {str(p.resolve()) for p in expected} <= set(files): return False
    return (receipt.get('evaluation_sha256') == sha(evaluation) and bool(files)
            and all(Path(p).is_file() and sha(p) == digest for p,digest in files.items()))


def plan(root, trait):
    from prsformer import validate_evaluation
    base = Path(root) / 'prsformer' / trait; run = base / 'run'
    evaluation = run / 'evaluation.json'
    if evaluation.is_file():
        validate_evaluation(run, Path(root) / 'prsformer/cache')
        # The training/evaluation task set is fixed; do not silently reuse another task set.
        import pandas as pd
        actual = set(pd.read_csv(run / 'metrics.tsv', sep='\t', usecols=['trait']).trait)
        if actual != {trait}:
            raise ValueError('Requested traits differ from the saved PRSformer evaluation')
        return 'skip' if publication_valid(base, evaluation) else 'report'
    model, training = run / 'model.pt', run / 'training.pt'
    trained = run / 'TRAINED.json'
    if trained.is_file():
        receipt = json.loads(trained.read_text())
        if not model.is_file() or sha(model) != receipt.get('model_sha256'):
            raise ValueError('Completed-training checkpoint changed; restore model.pt')
        return 'predict'
    if training.is_file() and model.is_file(): return 'resume'
    if model.is_file():
        import torch
        checkpoint = torch.load(model, map_location='cpu', weights_only=True)
        if checkpoint.get('history') and checkpoint.get('test_evaluated'):
            return 'predict'
        raise ValueError('Incomplete PRSformer checkpoint: model.pt exists without resumable training.pt')
    if training.exists(): raise ValueError('training.pt exists without the selected model.pt; restore model.pt')
    return 'train'


def dispatch(code_root, args, dry=False):
    from abm_core.cold_start import benchmark_parser, benchmark_config, benchmark_device, ensure_inputs
    options = benchmark_parser("prsformer", "Train or resume an independent PRSformer for each trait, sequentially.").parse_args(args)
    config = benchmark_config(options)
    root = Path(options.analysis_root)
    traits = options.traits.split(',')
    if dry or options.dry_run:
        print('[GRID] PRSformer plan: independent traits in order: ' + ' -> '.join(traits), flush=True)
        for trait in traits:
            print(f'[GRID] {trait}: skip, repair reports, resume or train at {root / "prsformer" / trait}; shared cache {root / "prsformer/cache"}', flush=True)
        return
    device = benchmark_device(options.device)
    from kinship import analysis_guard
    with analysis_guard(root), run_lock(root / 'prsformer'):
        ensure_inputs(options, config)
        for index, trait in enumerate(traits, 1):
            print(f'[GRID] PRSformer trait {index}/{len(traits)}: {trait}', flush=True)
            _run_trait(code_root, options, config, device, trait)


def _run_trait(code_root, options, config, device, trait):
    root = Path(options.analysis_root)
    base = root / 'prsformer' / trait
    action = plan(root, trait)
    if action == 'skip':
        print(f'[GRID] SKIP PRSformer {trait}: completed model, evaluation and reports verified', flush=True)
        return
    if action != 'report' and not (root / 'prsformer/cache/genotypes.npy').is_file():
        from abm_core.cold_start import preparation_command
        subprocess.run(preparation_command(options, config), check=True)
    entry = ['bash', str(code_root / 'f/prsformer.entry.sh'), '--python', sys.executable,
             '--traits', trait, '--no-replace', '--cache-dir', str(root / 'cache/abm'),
             '--out-root', str(root / 'abm'), '--prsformer-root', str(root),
             '--prsformer-run-dir', str(base / 'run'),
             '--prsformer-publication-dir', str(base), '--prsformer-score-dir', str(root / 'prsformer'),
             '--prsformer-device', device, '--prsformer-threads', str(options.threads)]
    env = dict(os.environ, PRSFORMER_KEEP_GENOTYPES='1')
    if action == 'report':
        print('[GRID] PRSformer: reuse saved evaluation; repair reports without training', flush=True)
        subprocess.run([*entry, '--prsformer-mode', 'report', '--prsformer-replace', 'TRUE'], check=True, env=env)
    elif action == 'predict':
        print('[GRID] PRSformer: reuse completed model; regenerate evaluation without training', flush=True)
        subprocess.run([*entry, '--prsformer-mode', 'predict'], check=True, env=env)
        subprocess.run([*entry, '--prsformer-mode', 'report', '--prsformer-replace', 'TRUE'], check=True, env=env)
    else:
        print('[GRID] PRSformer: ' + ('resume last saved epoch' if action == 'resume' else 'start new training'), flush=True)
        if action == 'resume':
            subprocess.run([*entry, '--prsformer-mode', 'train', '--prsformer-arg', '--resume'], check=True, env=env)
            subprocess.run([*entry, '--prsformer-mode', 'report', '--prsformer-replace', 'TRUE'], check=True, env=env)
        else:
            subprocess.run(entry, check=True, env=env)
