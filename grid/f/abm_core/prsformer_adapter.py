"""PRSformer2: shared backbone with separate development and calibration roles."""
import json
from pathlib import Path
import shlex
import subprocess
import sys

import pandas as pd

from .cold_start import (ROOT, benchmark_parser, benchmark_config,
                         benchmark_device, ensure_inputs)
from .grid_inputs import partition_calibration, digest
from .model import file_sha256, write_json
from .benchmark_results import performance


def save_roles(path, roles):
    """Preserve exact role bytes across resumes, including gzip metadata."""
    if path.is_file():
        saved = pd.read_csv(path, sep="\t", dtype=str)
        pd.testing.assert_frame_equal(saved, roles.astype(str).reset_index(drop=True))
    else:
        roles.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})


def action_for(local):
    if (local / "evaluation.json").is_file():
        return "publish"
    model, training, trained = (local / name for name in ("model.pt", "training.pt", "TRAINED.json"))
    if trained.is_file():
        receipt = json.loads(trained.read_text())
        if not model.is_file() or file_sha256(model) != receipt.get("model_sha256"):
            raise ValueError("PRSformer2 completed checkpoint changed")
        import torch
        checkpoint = torch.load(model, map_location="cpu", weights_only=True)
        if checkpoint.get("preprocessing", {}).get("prediction_calibration"):
            return "predict"
    if model.is_file() and training.is_file():
        return "resume"
    if model.exists() or training.exists() or trained.exists():
        raise ValueError("Incomplete PRSformer2 checkpoint; restore both model.pt and training.pt")
    return "train"


def verify_complete(run, cache, traits):
    receipt = json.loads((run / "COMPLETE.json").read_text())
    required = {"main_performance.tsv", "RUN_MANIFEST.json", "_LOCAL_ONLY/evaluation.json",
                "_LOCAL_ONLY/roles.tsv.gz"}
    if receipt.get("status") != "completed" or not required <= set(receipt.get("files", {})):
        raise ValueError("Incomplete PRSformer2 completion receipt")
    for name, expected in receipt["files"].items():
        if file_sha256(run / name) != expected:
            raise ValueError("PRSformer2 completed output changed: " + name)
    return performance(run / "_LOCAL_ONLY", cache, traits, "PRSformer2")


def publish(run, cache, traits):
    table = performance(run / "_LOCAL_ONLY", cache, traits, "PRSformer2")
    table.to_csv(run / "main_performance.tsv", sep="\t", index=False)
    files = ("main_performance.tsv", "RUN_MANIFEST.json", "_LOCAL_ONLY/evaluation.json",
             "_LOCAL_ONLY/roles.tsv.gz")
    write_json(run / "COMPLETE.json", dict(status="completed",
               files={name: file_sha256(run / name) for name in files}))
    write_json(run.parent / "CURRENT.json", dict(run=run.name))


def main(argv=None):
    p = benchmark_parser("prsformer2", "Train or resume PRSformer2 with separate prediction calibration.")
    p.add_argument("--upstream-dir", default="/mnt/f/software/PRSformer")
    p.add_argument("--check", action="store_true", help="Validate inputs and model runtime")
    args = p.parse_args(argv)
    config = benchmark_config(args)
    if args.dry_run:
        print("[GRID] PRSformer2 plan: prepare inputs; validate completed outputs or resume; publish", flush=True)
        return
    _run(args, config, benchmark_device(args.device))
    if not args.check:
        from dosage_storage import release_dosages
        release_dosages(config["block_config"], traits=config["traits"],
                        analysis_root=args.analysis_root, apply=True, best_effort=True)


def _run(args, config, device):
    from twin_core import Config, assign_roles
    from twin import verify_existing_prsformer_roster
    from resume import run_lock
    root = Path(args.analysis_root)
    with run_lock(root / "prsformer2"):
        ensure_inputs(args, config)
        cache = root / "prsformer/cache"
        preparation = json.loads((cache / "prepare.json").read_text())
        if file_sha256(cache / "data.tsv.gz") != preparation["sha256"]["data.tsv.gz"]:
            raise ValueError("PRSformer row identities changed")
        d = pd.read_csv(cache / "data.tsv.gz", sep="\t", dtype={"eid": str, "group": str})
        verify_existing_prsformer_roster(
            root, d.rename(columns={"group": "family_id"}), d.loc[d.split.eq("test"), "eid"])
        roles = d[["eid", "group"]].rename(columns={"eid": "id", "group": "group_id"})
        roles["role"] = assign_roles(roles.group_id, Config())
        roles.loc[d.split.eq("test"), "role"] = "test"
        roles = partition_calibration(roles, config["head_options"]["seed"])
        settings = config["prsformer"]
        signature = digest(dict(
            settings=settings, roles=roles.to_dict("list"), traits=args.traits, device=device,
            sources={name: file_sha256(ROOT / "f" / name) for name in
                     ("prsformer.py", "abm_core/prsformer_adapter.py", "abm_core/benchmark_results.py")},
            upstream_dir=str(Path(args.upstream_dir).resolve()), preparation=preparation))
        run = root / "prsformer2" / signature[:16]
        local = run / "_LOCAL_ONLY"
        local.mkdir(parents=True, exist_ok=True)
        save_roles(local / "roles.tsv.gz", roles)
        manifest = dict(signature=signature, method="PRSformer2",
            roles_sha256=file_sha256(local / "roles.tsv.gz"), settings=settings,
            genotype_QC_source="frozen_outer_development_cache_X_only",
            validation_status="exploratory_previously_inspected_test",
            information_budget="target individual genotypes/labels; no CSx/GWAS prior")
        manifest_path = run / "RUN_MANIFEST.json"
        if manifest_path.is_file():
            if json.loads(manifest_path.read_text()) != manifest:
                raise ValueError("PRSformer2 run manifest changed")
        else:
            write_json(manifest_path, manifest)
        traits = args.traits.split(",")
        if (run / "COMPLETE.json").is_file():
            verify_complete(run, cache, traits)
            write_json(run.parent / "CURRENT.json", dict(run=run.name))
            print("[GRID] SKIP PRSformer2: completed outputs verified", flush=True)
            return
        action = action_for(local)
        if action != "publish":
            if not (cache / "genotypes.npy").is_file():
                from .cold_start import preparation_command
                subprocess.run(preparation_command(args, config), check=True)
            cmd = [sys.executable, str(ROOT / "f/prsformer.py"),
                   "predict" if action == "predict" else "train",
                   "--genotypes", str(cache / "genotypes.npy"), "--data", str(cache / "data.tsv.gz"),
                   "--variants", str(cache / "variants.tsv.gz"), "--upstream-dir", args.upstream_dir,
                   "--out-dir", str(local), "--traits", args.traits, "--covariates", "age,sex,PC1,PC2",
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
        publish(run, cache, traits)


if __name__ == "__main__":
    main()
