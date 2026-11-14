#!/usr/bin/env python3
"""Independent PRSformer2 traits sharing inputs and frozen calibration roles."""
from __future__ import annotations

# Prsformer adapter

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.append(str(_paths_root / "0f"))
from workflow_paths import configured_path
from workflow_layout import named_directory
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import pandas as pd

from abm.cold_start import (ROOT, benchmark_parser, benchmark_config,
                         benchmark_device, ensure_inputs)
from abm.grid_inputs import partition_calibration, digest
from abm.model import file_record, write_json
from abm.reporting import performance


def save_roles(path, roles):
    """Preserve exact role bytes across resumes, including gzip metadata."""
    if path.is_file():
        saved = pd.read_csv(path, sep="\t", dtype=str)
        pd.testing.assert_frame_equal(saved, roles.astype(str).reset_index(drop=True))
    else:
        roles.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})


def action_for(local):
    required = ("evaluation.json", "model.pt", "history.tsv", "metrics.tsv", "test_predictions.tsv.gz")
    return "publish" if all((local/name).is_file() and (local/name).stat().st_size > 0 for name in required) else "train"


def verify_complete(run, cache, traits):
    receipt = json.loads((run / "COMPLETE.json").read_text())
    required = {"main_performance.tsv", "RUN_MANIFEST.json", "_LOCAL_ONLY/evaluation.json",
                "_LOCAL_ONLY/roles.tsv.gz"}
    if receipt.get("status") != "completed" or not required <= set(receipt.get("files", {})):
        raise ValueError("Incomplete PRSformer2 completion receipt")
    for name, expected in receipt["files"].items():
        if not (run / name).is_file() or (run / name).stat().st_size == 0:
            raise ValueError("PRSformer2 completed output changed: " + name)
    return performance(run / "_LOCAL_ONLY", cache, traits, "PRSformer2",
                       exact_traits=True)


def publish(run, cache, traits):
    table = performance(run / "_LOCAL_ONLY", cache, traits, "PRSformer2",
                        exact_traits=True)
    table.to_csv(run / "main_performance.tsv", sep="\t", index=False)
    files = ("main_performance.tsv", "RUN_MANIFEST.json", "_LOCAL_ONLY/evaluation.json",
             "_LOCAL_ONLY/roles.tsv.gz")
    write_json(run / "COMPLETE.json", dict(status="completed",
               files={name: file_record(run / name) for name in files}))
    write_json(run.parent / "CURRENT.json", dict(run=run.name))


def main(argv=None):
    p = benchmark_parser("prsformer2", "Train or resume independent PRSformer2 traits sequentially, with separate prediction calibration.")
    p.add_argument("--upstream-dir", default=configured_path('DIR_SOFTWARE', 'PRSformer'))
    p.add_argument("--check", action="store_true", help="Validate inputs and model runtime")
    args = p.parse_args(argv)
    config = benchmark_config(args)
    if args.dry_run:
        traits = args.traits.split(",")
        print("[GRID] PRSformer2 plan: independent traits in order: " + " -> ".join(traits), flush=True)
        for trait in traits:
            print(f'[GRID] {trait}: skip, resume or train at {Path(args.analysis_root) / "tmp/prsformer2" / trait}; shared cache {Path(args.analysis_root) / "tmp/prsformer/cache"}', flush=True)
        return
    if args.workflow_stage == "report":
        from prsformer_workflow import frozen_model
        from grid_paths import read_config
        from resume import run_lock
        from kinship import analysis_guard
        root = Path(args.analysis_root)
        cache = Path(read_config(config["block_config"], root)["prsformer_cache"])
        with analysis_guard(root):
            for trait in args.traits.split(","):
                with run_lock(root / "tmp/prsformer2" / trait):
                    local = frozen_model(root, "prsformer2", trait, cache)
                    publish(local.parent, cache, [trait])
        return
    _run(args, config, benchmark_device(args.device))
    if not args.check and os.environ.get("GRID_KEEP_SHARED_INPUTS") != "1":
        from dosage_storage import release_dosages
        release_dosages(config["block_config"], traits=config["traits"],
                        analysis_root=args.analysis_root, apply=True, best_effort=True)


def _run(args, config, device):
    from twin_core import Config, assign_roles
    from twin import verify_existing_prsformer_roster
    from resume import run_lock
    root = Path(args.analysis_root)
    from kinship import analysis_guard
    with analysis_guard(root):
        ensure_inputs(args, config)
        from grid_paths import read_config
        block = read_config(config["block_config"], root)
        cache = Path(block["prsformer_cache"])
        preparation = json.loads((cache / "prepare.json").read_text())
        if file_record(cache / "data.tsv.gz") != preparation["file_record"]["data.tsv.gz"]:
            raise ValueError("PRSformer row identities changed")
        d = pd.read_csv(cache / "data.tsv.gz", sep="\t", dtype={"eid": str, "group": str})
        verify_existing_prsformer_roster(
            root, d.rename(columns={"group": "family_id"}), d.loc[d.split.eq("test"), "eid"])
        from kinship import correct_benchmark
        d, correction = correct_benchmark(d, root, 20260904, separate_calibration=True)
        roles = d[["eid", "group"]].rename(columns={"eid": "id", "group": "group_id"})
        roles["role"] = assign_roles(roles.group_id, Config())
        roles.loc[d.split.eq("test"), "role"] = "test"
        if correction:
            eligible = ~d.split.eq("unused")
            kept = partition_calibration(roles.loc[eligible].copy(), config["head_options"]["seed"])
            roles["parent_role"] = roles.role
            roles.loc[~eligible, ["role", "parent_role"]] = "excluded_development_relative"
            roles.loc[eligible, ["role", "parent_role"]] = kept[["role", "parent_role"]]
        else:
            roles = partition_calibration(roles, config["head_options"]["seed"])
        traits = args.traits.split(",")
        for index, trait in enumerate(traits, 1):
            print(f"[GRID] PRSformer2 trait {index}/{len(traits)}: {trait}", flush=True)
            with run_lock(root / "tmp/prsformer2" / trait):
                run = _run_trait(args, config, device, cache, preparation, roles, trait)
                if args.workflow_stage == "model" and not args.check:
                    from prsformer_workflow import freeze_model
                    freeze_model(root, "prsformer2", trait, run / "_LOCAL_ONLY", cache)


def _run_trait(args, config, device, cache, preparation, roles, trait):
    root = Path(args.analysis_root)
    settings = config["prsformer"]
    run = root / "tmp/prsformer2" / trait / "model"
    run.mkdir(parents=True, exist_ok=True)
    signature = "prsformer2"
    if os.environ.get("GRID_REPLACE", "TRUE") == "TRUE":
        import shutil
        shutil.rmtree(run / "_LOCAL_ONLY", ignore_errors=True)
        (run / "COMPLETE.json").unlink(missing_ok=True)
        (run / "RUN_MANIFEST.json").unlink(missing_ok=True)
    local = run / "_LOCAL_ONLY"
    local.mkdir(parents=True, exist_ok=True)
    save_roles(local / "roles.tsv.gz", roles)
    manifest = dict(signature=signature, method="PRSformer2", traits=[trait], training_mode="independent_trait",
        roles_file_record=file_record(local / "roles.tsv.gz"), settings=settings,
        genotype_QC_source="frozen_outer_development_cache_X_only",
        validation_status="exploratory_previously_inspected_test",
        information_budget="target individual genotypes/labels; no CSx/GWAS prior")
    manifest_path = run / "RUN_MANIFEST.json"
    if manifest_path.is_file():
        pass  # Complete output reuse is explicit; no code-hash gate.
    else:
        write_json(manifest_path, manifest)
    traits = [trait]
    from dosage_storage import dosage_lock, prsformer_cache_lock
    from grid_paths import read_config
    block = read_config(config["block_config"], root)
    stage = getattr(args, "workflow_stage", "all")
    if (run / "COMPLETE.json").is_file():
        try:
            with dosage_lock(block, wait=True), prsformer_cache_lock(cache):
                verify_complete(run, cache, traits)
            if stage != "model":
                write_json(run.parent / "CURRENT.json", dict(run=run.name))
            print(f"[GRID] REUSE PRSformer2 {trait}: complete outputs", flush=True)
            return run
        except (OSError, ValueError, KeyError, TypeError):
            (run / "COMPLETE.json").unlink(missing_ok=True)
    action = action_for(local)
    if action == "train":
        import shutil
        shutil.rmtree(local, ignore_errors=True)
        local.mkdir(parents=True, exist_ok=True)
        save_roles(local / "roles.tsv.gz", roles)
    if action != "publish" and not (cache / "genotypes.npy").is_file():
        if stage == "model":
            raise FileNotFoundError("Run CPU preparation before GPU model training")
        from copy import copy
        from shared_cohort import ensure_inputs as restore_inputs
        options = copy(args)
        options.complete_features = False
        restore_inputs(options, config)
    with dosage_lock(block, wait=True), prsformer_cache_lock(cache):
        if action != "publish":
            cmd = [sys.executable, str(ROOT / "f/prsformer.py"),
                   "predict" if action == "predict" else "train",
                   "--genotypes", str(cache / "genotypes.npy"), "--data", str(cache / "data.tsv.gz"),
                   "--variants", str(cache / "variants.tsv.gz"), "--upstream-dir", args.upstream_dir,
               "--analysis-root", str(root),
                   "--out-dir", str(local), "--traits", trait, "--covariates", "age,sex,PC1,PC2",
                   "--device", device, "--threads", str(args.threads),
                   "--calibration-roles", str(local / "roles.tsv.gz")]
            for key, value in settings.items():
                cmd += ["--" + key.replace("_", "-"), str(value)]
            if action == "resume":
                cmd.append("--resume")
            elif action == "predict":
                cmd += ["--checkpoint", str(local / "model.pt")]
            if args.check:
                cmd.append("--check")
            print("[GRID] PRSformer2 " + action + ": " + shlex.join(cmd), flush=True)
            subprocess.run(cmd, check=True)
            if args.check:
                return
        if stage != "model":
            publish(run, cache, traits)
    return run



if __name__ == "__main__":
    main()
