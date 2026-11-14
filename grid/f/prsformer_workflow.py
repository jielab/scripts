#!/usr/bin/env python3
"""Bind CPU publication to the exact outputs of its preceding model worker."""
from __future__ import annotations

# prsformer_stages

import json
import os
from pathlib import Path

from resume import atomic_json, file_record


def receipt_path(method, trait):
    directory = os.environ.get("GRID_PRSFORMER_STAGE_DIR")
    if not directory:
        raise ValueError("Split PRSformer workers require their planner's stage directory")
    return Path(directory) / f"{method}.{trait}.json"


def freeze_model(root, method, trait, run, cache):
    from prsformer import validate_evaluation
    root, run, cache = Path(root).resolve(), Path(run).resolve(), Path(cache).resolve()
    validate_evaluation(run, cache, analysis_root=root)
    paths = [run / name for name in
             ("evaluation.json", "model.pt", "metrics.tsv", "test_predictions.tsv.gz")]
    paths += [cache / name for name in ("prepare.json", "data.tsv.gz", "variants.tsv.gz")]
    if method == "prsformer2":
        paths += [run / "roles.tsv.gz", run.parent / "RUN_MANIFEST.json"]
    active = root / "local/kinship/ACTIVE.json"
    record = dict(schema="grid.prsformer.model_handoff.1", method=method, trait=trait,
                  root=str(root), run=str(run), cache=str(cache),
                  active_file_record=file_record(active) if active.is_file() else None,
                  files={str(p): {"size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns} for p in paths})
    atomic_json(receipt_path(method, trait), record)
    return run


def frozen_model(root, method, trait, cache):
    from prsformer import validate_evaluation
    root, cache = Path(root).resolve(), Path(cache).resolve()
    path = receipt_path(method, trait)
    if not path.is_file() and os.environ.get("GRID_REPORT_ONLY") == "1":
        if method == "prsformer":
            run = root / "tmp/prsformer" / trait / "run"
        else:
            base = root / "tmp/prsformer2" / trait
            current = json.loads((base / "CURRENT.json").read_text())
            run = base / current["run"] / "_LOCAL_ONLY"
        validate_evaluation(run, cache, analysis_root=root)
        freeze_model(root, method, trait, run, cache)
    record = json.loads(path.read_text())
    if record.get("method") != method or record.get("trait") != trait or record.get("root") != str(root):
        raise ValueError("Wrong PRSformer model handoff")
    if record.get("cache") != str(cache):
        raise ValueError("Wrong PRSformer cache in model handoff")
    run = Path(record["run"])
    run.relative_to(root / "tmp" / method / trait)
    for filename, expected in record["files"].items():
        source = Path(filename)
        stat = source.stat()
        if {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns} != expected:
            raise ValueError("Model/report inputs changed between scheduled stages: " + filename)
    validate_evaluation(run, cache, analysis_root=root)
    return run

# prsformer_resume
"""Serial, independent PRSformer traits sharing one frozen input cache."""
import subprocess
import sys
from resume import run_lock


def publication_valid(base, evaluation):
    marker = base / 'PUBLISHED.json'
    if not marker.is_file():
        return False
    try:
        receipt = json.loads(marker.read_text())
        return bool(receipt.get('files')) and all(Path(p).is_file() and Path(p).stat().st_size > 0 for p in receipt['files'])
    except (OSError, ValueError, TypeError):
        return False


def plan(root, trait, cache=None):
    from prsformer import validate_evaluation
    base = Path(root) / 'tmp/prsformer' / trait
    run = base / 'run'
    try:
        validate_evaluation(run, Path(cache) if cache else Path(root)/'tmp/prsformer/cache', analysis_root=root)
        import pandas as pd
        actual = set(pd.read_csv(run/'metrics.tsv', sep='\t', usecols=['trait']).trait)
        if actual != {trait}:
            return 'train'
        return 'skip' if publication_valid(base, run/'evaluation.json') else 'report'
    except (OSError, ValueError, KeyError, TypeError):
        return 'train'


def dispatch(code_root, args, dry=False):
    from abm.cold_start import benchmark_parser, benchmark_config, benchmark_device, ensure_inputs
    options = benchmark_parser("prsformer", "Train or resume an independent PRSformer for each trait, sequentially.").parse_args(args)
    config = benchmark_config(options)
    root = Path(options.analysis_root)
    traits = options.traits.split(',')
    if dry or options.dry_run:
        print('[GRID] PRSformer plan: independent traits in order: ' + ' -> '.join(traits), flush=True)
        for trait in traits:
            print(f'[GRID] {trait}: skip, repair reports, resume or train at {root / "tmp/prsformer" / trait}; shared cache {root / "tmp/prsformer/cache"}', flush=True)
        return
    stage = options.workflow_stage
    device = "cpu" if stage == "report" else benchmark_device(options.device)
    from kinship import analysis_guard
    with analysis_guard(root):
        if stage != "report":
            ensure_inputs(options, config)
        for index, trait in enumerate(traits, 1):
            print(f'[GRID] PRSformer trait {index}/{len(traits)}: {trait}', flush=True)
            with run_lock(root / 'tmp/prsformer' / trait):
                _run_trait(code_root, options, config, device, trait)
                if stage == "model":
                    from grid_paths import read_config
                    cache = read_config(config["block_config"], root)["prsformer_cache"]
                    freeze_model(root, "prsformer", trait, root / "tmp/prsformer" / trait / "run", cache)


def _run_trait(code_root, options, config, device, trait):
    root = Path(options.analysis_root)
    from grid_paths import read_config
    block = read_config(config['block_config'], root)
    cache = Path(block['prsformer_cache'])
    base = root / 'tmp/prsformer' / trait
    stage = getattr(options, 'workflow_stage', 'all')
    if stage != "report" and os.environ.get("GRID_REPLACE", "TRUE") == "TRUE":
        shutil.rmtree(base / "run", ignore_errors=True)
        (base / "PUBLISHED.json").unlink(missing_ok=True)
    action = plan(root, trait, cache)
    if action == "train" and stage != "report":
        shutil.rmtree(base / "run", ignore_errors=True)
    if stage == 'report' and action not in {'report', 'skip'}:
        raise ValueError('CPU publication requires completed model evaluation; it cannot train or predict')
    if stage == 'model' and action in {'report', 'skip'}:
        print(f'[GRID] REUSE PRSformer {trait}: model evaluation already complete', flush=True)
        return
    if stage == 'report' and action == 'skip' and os.environ.get('GRID_REPLACE', 'TRUE') == 'TRUE':
        action = 'report'
    if action == 'skip':
        print(f'[GRID] SKIP PRSformer {trait}: completed model, evaluation and reports verified', flush=True)
        return
    if action != 'report' and not (cache / 'genotypes.npy').is_file():
        if stage == 'model':
            raise FileNotFoundError('Run the CPU preparation stage before the GPU model stage')
        from copy import copy
        from shared_cohort import ensure_inputs
        preparation = copy(options)
        preparation.complete_features = False
        ensure_inputs(preparation, config)
    entry = ['bash', str(code_root / 'f/prsformer.entry.sh'), '--python', sys.executable,
             '--traits', trait, '--no-replace', '--cache-dir', str(root / 'tmp/cache/abm'),
             '--out-root', str(root / 'tmp/abm'), '--prsformer-root', str(root),
             '--prsformer-cache-dir', str(cache),
             '--prsformer-run-dir', str(base / 'run'),
             '--prsformer-publication-dir', str(base), '--prsformer-score-dir', str(root / 'tmp/prsformer'),
             '--prsformer-device', device, '--prsformer-threads', str(options.threads)]
    env = dict(os.environ, PRSFORMER_KEEP_GENOTYPES='1')
    from dosage_storage import dosage_lock
    # Different traits share read-only dosage. The wrapper takes the matching
    # shared pipeline lock; producers and release still require exclusive locks.
    with dosage_lock(block, wait=True):
        if action == 'report':
            print('[GRID] PRSformer: reuse saved evaluation; repair reports without training', flush=True)
            subprocess.run([*entry, '--prsformer-mode', 'report', '--prsformer-replace', 'TRUE'], check=True, env=env)
        elif action == 'predict':
            print('[GRID] PRSformer: reuse completed model; regenerate evaluation without training', flush=True)
            subprocess.run([*entry, '--prsformer-mode', 'predict'], check=True, env=env)
            if stage != 'model':
                subprocess.run([*entry, '--prsformer-mode', 'report', '--prsformer-replace', 'TRUE'], check=True, env=env)
        else:
            print('[GRID] PRSformer: ' + ('resume last saved epoch' if action == 'resume' else 'start new training'), flush=True)
            resume = ['--prsformer-arg', '--resume'] if action == 'resume' else []
            subprocess.run([*entry, '--prsformer-mode', 'train', *resume], check=True, env=env)
            if stage != 'model':
                subprocess.run([*entry, '--prsformer-mode', 'report', '--prsformer-replace', 'TRUE'], check=True, env=env)

# prsformer_benchmark
"""Publish a verified standalone benchmark receipt after PRSformer evaluation.

This module only reads an already validated fitted model/evaluation. It does not
calibrate, tune, retrain or change any original PRSformer performance measure.
"""
import shutil
import numpy as np
import pandas as pd


def publish(root, trait, method="prsformer", cache=None):
    from shared_cohort import load_cohort, file_record
    from target_benchmark import metrics
    from prsformer import validate_evaluation
    from kinship import analysis_guard
    root = Path(root)
    if method == "prsformer":
        run = root / "tmp/prsformer" / trait / "run"
        label = "PRSformer independent-trait adaptation"
    elif method == "prsformer2":
        base = root / "tmp/prsformer2" / trait
        current = json.loads((base / "CURRENT.json").read_text())
        run = base / current["run"] / "_LOCAL_ONLY"
        label = "PRSformer2 exploratory extension"
    else:
        raise ValueError("Unknown PRSformer benchmark")
    cache = Path(cache) if cache is not None else root / "tmp/prsformer/cache"
    out = root / "qc" / method / trait
    out.mkdir(parents=True, exist_ok=True)
    with analysis_guard(root), run_lock(out):
        validate_evaluation(run, cache, analysis_root=root)
        frame, metadata = load_cohort(root, trait, cohort=cache / "data.tsv.gz")
        test = frame.loc[frame.split.eq("test")].copy()
        saved = pd.read_csv(run / "test_predictions.tsv.gz", sep="\t", dtype={"eid":str})
        saved = saved.loc[saved.trait.eq(trait)]
        if saved.eid.duplicated().any() or set(saved.eid) != set(test.id) or not saved.split.eq("test").all():
            raise ValueError("PRSformer test prediction roster differs from the shared untouched test")
        saved = saved.set_index("eid").loc[test.id]
        if not np.allclose(saved.y.to_numpy(float), test.y.to_numpy(float), rtol=0, atol=1e-10, equal_nan=True):
            raise ValueError("PRSformer test phenotype/units differ from shared cohort")
        prediction = saved.prediction.to_numpy(float)
        if not np.isfinite(prediction).all():
            raise ValueError("PRSformer does not predict every frozen test ID")
        local = out / "_LOCAL_ONLY"
        local.mkdir(exist_ok=True)
        pred = test[["id", "group_id", "ancestry"]].copy()
        pred["role"] = "test"
        pred[method + ".prediction"] = prediction
        pred_path, truth_path = local / "final_prs.tsv.gz", local / "test_truth.tsv.gz"
        pred.to_csv(pred_path, sep="\t", index=False)
        test[["id", "y"]].to_csv(truth_path, sep="\t", index=False)
        binary = trait == "t2dm"
        pd.DataFrame([dict(trait=trait, method=label, model=method, source="4.prsformer.sh", stage="test",
            status="completed", endpoint=metadata["endpoint_definition"],
            **metrics(test.y.to_numpy(float), prediction, binary))]).to_csv(out / "performance.tsv", sep="\t", index=False)
        shutil.copyfile(run / "metrics.tsv", out / "paper_style_metrics.tsv")
        source_files = [run / name for name in ("evaluation.json", "model.pt", "test_predictions.tsv.gz", "metrics.tsv")]
        files = [*source_files, cache / "prepare.json", cache / "data.tsv.gz", pred_path, truth_path,
                 out / "performance.tsv", out / "paper_style_metrics.tsv"]
        receipt = dict(schema="grid.target_benchmark.1", status="completed", method=method, trait=trait,
            models=[method], model_labels={method:label}, method_label=label,
            task="binary" if binary else "continuous", endpoint=metadata["endpoint_definition"],
            predictions_file=str(pred_path.resolve()), test_truth_file=str(truth_path.resolve()),
            test_labels_used_only_after_prediction_freeze=True, metadata=metadata,
            training_budget="individual-genotype train/validation split; not the small PRS-score validation budget",
            roster=dict(test_N=len(test), source=str(cache / "data.tsv.gz")),
            files={str(path.resolve()):file_record(path) for path in files}, implementation_file_record=file_record(__file__))
        atomic_json(out / "COMPLETE.json", receipt)
        print(f"[{method}] verified independent benchmark: {out / 'COMPLETE.json'}", flush=True)
        from output_layout import publish_prsformer
        publish_prsformer(root, trait, method, run, out)
        return receipt

# prsformer_dispatch
#!/usr/bin/env python3
"""Portable dispatcher for the retained, explicitly labelled PRSformer adapters."""

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.append(str(_paths_root / "0f"))
from workflow_paths import configured_path
import argparse
import hashlib
from grid_paths import GRID, analysis_default, read_config, relocate


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in {"prsformer", "prsformer2"}:
        raise ValueError("Use prsformer or prsformer2")
    module, argv = argv[0], argv[1:]
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--analysis-root", "--output-root", default=analysis_default())
    p.add_argument("--config", default=str(GRID / "config/benchmarks.yml"))
    p.add_argument("--workflow-stage", choices=("all", "model", "report"), default="all")
    known, forwarded = p.parse_known_args(argv)
    root = Path(known.analysis_root)
    if "--help" in forwarded or "-h" in forwarded:
        subprocess.run([sys.executable, str(GRID / "f/twin.py"), module, "--help"], check=True)
        return
    config = read_config(known.config, root)
    block = read_config(config["block_config"], root)
    if "--dry-run" in forwarded:
        traits = forwarded[forwarded.index("--traits")+1] if "--traits" in forwarded else "height,ldl,t2dm"
        print(json.dumps(dict(module=module, traits=traits, analysis_root=str(root),
            training_mode="independent_trait_adaptation" if module == "prsformer" else "exploratory_extension",
            shared_cache=block["prsformer_cache"], owns_test_evaluation=True)))
        return
    folder = root / "tmp/cache/benchmark_config"
    folder.mkdir(parents=True, exist_ok=True)
    with run_lock(folder):
        block_digest = hashlib.sha256(json.dumps(block, sort_keys=True).encode()).hexdigest()[:20]
        block_path = folder / ("blocks." + block_digest + ".json")
        if not block_path.is_file():
            atomic_json(block_path, block)
        config["block_config"] = str(block_path)
        digest = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()[:20]
        path = folder / ("prsformer." + digest + ".json")
        if not path.is_file():
            atomic_json(path, config)
    if module == "prsformer2" and "--upstream-dir" not in forwarded:
        forwarded += ["--upstream-dir", os.environ.get("PRSFORMER_UPSTREAM_DIR", relocate(configured_path('DIR_SOFTWARE', 'PRSformer')))]
    cache = Path(block["prsformer_cache"])
    if known.workflow_stage == "model" and not (cache / "genotypes.npy").is_file():
        raise FileNotFoundError("GPU model stage requires the CPU-prepared genotype cache")
    if known.workflow_stage == "report":
        traits = forwarded[forwarded.index("--traits")+1] if "--traits" in forwarded else "height,ldl,t2dm"
        for trait in traits.split(","):
            frozen_model(root, module, trait, cache)
    if known.workflow_stage == "all" and "--check" not in forwarded and not all((cache / name).is_file() for name in
            ("prepare.json", "data.tsv.gz", "variants.tsv.gz")):
        from shared_cohort import ensure_inputs
        from types import SimpleNamespace
        threads = int(forwarded[forwarded.index("--threads")+1]) if "--threads" in forwarded else 4
        ensure_inputs(SimpleNamespace(analysis_root=str(root), config=str(path), block_config=str(block_path),
            traits="height,ldl,t2dm", threads=threads, complete_features=False, dry_run=False), config)
    command = [sys.executable, str(GRID / "f/twin.py"), module, "--analysis-root", str(root),
               "--config", str(path), "--workflow-stage", known.workflow_stage, *forwarded]
    subprocess.run(command, check=True, env=dict(os.environ, GRID_WORKFLOW_ANALYSIS_ROOT=str(root)))
    if known.workflow_stage != "model" and "--check" not in forwarded:
        traits = forwarded[forwarded.index("--traits")+1] if "--traits" in forwarded else "height,ldl,t2dm"
        for trait in traits.split(","):
            publish(root, trait, module, cache=block["prsformer_cache"])



if __name__ == "__main__":
    main()
