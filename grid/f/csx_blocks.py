#!/usr/bin/env python3
"""Independent CSx block feature/fit/evaluation modules owned by 2.csx.sh.

blockprepare is outcome-free. blockridge and blocknet are explicitly labelled
GRID extensions trained within a small group-isolated development budget. They
are not the transparent twin-count method and are not official PRS-CSx models.
"""
from __future__ import annotations

#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path
import os

import joblib
import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.linear_model import LogisticRegression, Ridge
from grid_paths import GRID, analysis_default, read_config
from shared_cohort import frozen_roster, load_cohort, file_record
from resume import atomic_json, run_lock
from target_benchmark import fit_glm, metrics, predict_glm


def effective_config(args):
    config = read_config(args.block_config, args.analysis_root)
    if args.score_dir:
        config["score_root"] = args.score_dir
    for name in ("epochs", "batch_size"):
        if getattr(args, name, None) is not None:
            config[name] = getattr(args, name)
    config["seed"] = args.seed
    if args.ld_block_dir:
        config["ld_block_dir"] = args.ld_block_dir
    if args.blocks:
        requested = "ALL.37.bed" if args.blocks == 1 else "AFR.37.bed,EAS.37.bed,EUR.37.bed,SAS.37.bed"
        if args.ld_block_files and args.ld_block_files != requested:
            raise ValueError("--blocks and --ld-block-files disagree")
        config["ld_block_files"] = requested
    elif args.ld_block_files:
        config["ld_block_files"] = args.ld_block_files
    folder = Path(args.analysis_root) / "tmp/cache/csx_blocks/config"
    folder.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    path = folder / (digest[:20] + ".json")
    with run_lock(folder):
        if not path.is_file():
            atomic_json(path, config)
        elif json.loads(path.read_text()) != config:
            raise ValueError("Resolved block configuration identity changed")
    return path, config


def scale_blocks(array, rows, batch=256):
    """Moments from block-build genotypes only, using bounded row reads."""
    total = np.zeros(array.shape[1:], np.float64)
    squared = np.zeros_like(total)
    if not len(rows):
        raise ValueError("Empty block-build group")
    for start in range(0, len(rows), batch):
        x = np.asarray(array[rows[start:start+batch]], float)
        if not np.isfinite(x).all():
            raise ValueError("Nonfinite block feature")
        total += x.sum(0)
        squared += np.square(x).sum(0)
    mean = total / len(rows)
    sd = np.sqrt(np.maximum(squared / len(rows) - mean**2, 1e-12))
    return mean.astype(np.float32), sd.astype(np.float32)


def design(array, rows, side, mean, sd):
    x = (np.asarray(array[rows], np.float32) - mean) / sd
    return np.column_stack([x.reshape(len(rows), -1), side]).astype(np.float64, copy=False)


def ridge_fit(xbuild, ybuild, xtune, ytune, binary, penalties=(1., 10., 100., 1000.)):
    """Hyperparameters see tune labels only; this API has no test argument."""
    if not np.isfinite(ybuild).all() or not np.isfinite(ytune).all():
        raise ValueError("Build/tune labels must be observed")
    if binary and (len(np.unique(ybuild)) != 2 or len(np.unique(ytune)) != 2):
        raise ValueError("Block build and tuning require both classes")
    history, fitted = [], []
    for penalty in penalties:
        model = (LogisticRegression(C=1/penalty, solver="lbfgs", max_iter=5000, tol=1e-7)
                 if binary else Ridge(alpha=penalty, solver="lsqr", tol=1e-7))
        model.fit(xbuild, ybuild)
        pred = model.predict_proba(xtune)[:, 1] if binary else model.predict(xtune)
        loss = float(np.mean((pred-ytune)**2))
        history.append(dict(penalty=penalty, tune_loss=loss))
        fitted.append(model)
    best = int(np.argmin([r["tune_loss"] for r in history]))
    return fitted[best], pd.DataFrame(history), dict(selected_penalty=history[best]["penalty"],
        tuning_objective="brier" if binary else "mse", refit_on_tuning=False)


def fit_trait(args, trait, config_path, config):
    from abm.grid_inputs import locate_blocks
    from kinship import analysis_guard
    outer = dict(block_config=str(config_path), execution={"device": args.device})
    array, array_cohort, tokens, block_receipt, array_path = locate_blocks(outer, trait, allow_build=True)
    try:
        cohort = args.cohort
        if not cohort and not args.prepared:
            cohort = Path(config["prsformer_cache"]) / "data.tsv.gz"
        frame, cohort_meta = load_cohort(args.analysis_root, trait, cohort, args.prepared, args.group_col)
        roster, roster_meta = frozen_roster(frame, args.analysis_root, args.validation_n, args.seed)
        frame = frame.merge(roster[["id", "role", "block_role"]], on="id", validate="one_to_one", sort=False)
        lookup = pd.Series(np.arange(len(array_cohort)), index=array_cohort.eid.astype(str))
        if not set(frame.id) <= set(lookup.index) or not lookup.index.is_unique:
            raise ValueError("Block array does not cover the fixed benchmark roster uniquely")
        rows = lookup.loc[frame.id].to_numpy(int)
        covariates = [] if args.covariates.lower() in {"none", ""} else args.covariates.split(",")
        if set(covariates) - set(frame):
            raise ValueError("Declared block covariate missing")
        side = frame[covariates].to_numpy(float)
        if not np.isfinite(side).all():
            raise ValueError("Block covariates must cover the full test roster")
        role = frame.block_role.to_numpy(str)
        build, tune, test = role == "build", role == "tune_model", role == "test"
        observed = np.isfinite(frame.y.to_numpy(float))
        build &= observed
        tune &= observed
        if build.sum() < 10 or tune.sum() < 5:
            raise ValueError("Small development budget has too few observed build/tuning outcomes")
        # Missing-label development subjects stay in the roster, never becoming
        # test subjects; they simply cannot contribute to a supervised fit.
        fitting_roles = role.copy()
        fitting_roles[~observed & ~test] = "development_unused"
        root = Path(args.analysis_root) / "csx/blocks" / block_receipt["layout"]["layout_id"] / trait / args.module
        root.mkdir(parents=True, exist_ok=True)
        local = root / "_LOCAL_ONLY"
        local.mkdir(exist_ok=True)
        contract = dict(module=args.module, trait=trait, config=config, roster=roster_meta,
            cohort_source_file_record=cohort_meta["source_file_record"], block_signature=block_receipt["signature"],
            implementation_file_record=file_record(__file__), blocknet_source_file_record=file_record(Path(__file__).with_name("blocknet.py")),
            covariates=covariates, validation_n=args.validation_n)
        with analysis_guard(args.analysis_root), run_lock(root):
            completion = root / "COMPLETE.json"
            if completion.is_file() and not args.replace:
                try:
                    old = json.loads(completion.read_text())
                    if old.get("status") == "completed" and old.get("files") and all(Path(p).is_file() and Path(p).stat().st_size > 0 for p in old["files"]):
                        print(f"[{args.module}] REUSE {trait}: {root}", flush=True)
                        return
                except (OSError, ValueError, TypeError):
                    pass
            completion.unlink(missing_ok=True)
            mean, sd = scale_blocks(array, rows[build])
            side_mean, side_sd = side[build].mean(0), side[build].std(0)
            side_sd = np.where(side_sd > 1e-12, side_sd, 1.)
            side = (side-side_mean) / side_sd
            binary = trait == "t2dm"
            y_build = frame.loc[build, "y"].to_numpy(float)
            y_tune = frame.loc[tune, "y"].to_numpy(float)
            if args.module == "blockridge":
                xbuild = design(array, rows[build], side[build], mean, sd)
                xtune = design(array, rows[tune], side[tune], mean, sd)
                model, history, fit_meta = ridge_fit(xbuild, y_build, xtune, y_tune, binary)
                del xbuild, xtune
                predictions = np.empty(len(frame), float)
                for start in range(0, len(frame), 256):
                    end = min(start+256, len(frame))
                    x = design(array, rows[start:end], side[start:end], mean, sd)
                    predictions[start:end] = model.predict_proba(x)[:, 1] if binary else model.predict(x)
                joblib.dump(dict(model=model, block_mean=mean, block_sd=sd, side_mean=side_mean,
                    side_sd=side_sd, covariates=covariates), local / "frozen_model.joblib", compress=3)
                model_path = local / "frozen_model.joblib"
            else:
                import torch
                from abm.support import _backend
                from blocknet import train, predict_model
                selected, _ = _backend(args.device)
                device = torch.device(selected)
                torch.set_num_threads(args.threads)
                torch.manual_seed(args.seed)
                cfg = dict(config, _binary=binary, binary_objective="brier", initialize_prevalence=True)
                # The existing network receives no test outcomes, including in
                # early stopping or binary initialization.
                yfit = np.full(len(frame), np.nan)
                yfit[build], yfit[tune] = y_build, y_tune
                model, history, fit_meta = train("transformer", array, rows, side, yfit,
                    fitting_roles, tokens, cfg, device, mean, sd)
                raw, _ = predict_model(model, "transformer", array, rows, side,
                    np.arange(len(frame)), mean, sd, cfg, device, fit_meta["ymean"], fit_meta["ysd"])
                predictions = expit(raw) if binary else raw
                model_path = local / "frozen_model.pt"
                torch.save(dict(state_dict={k:v.cpu() for k,v in model.state_dict().items()}, config=cfg,
                    block_mean=mean, block_sd=sd, side_mean=side_mean, side_sd=side_sd,
                    covariates=covariates, outcome_scale=fit_meta), model_path)
            nuisance = fit_glm(side[build], y_build, binary)
            baseline = predict_glm(nuisance, side)
            prediction_table = frame[["id", "group_id", "split", "role", "block_role", "ancestry"]].copy()
            prediction_table["prediction"] = predictions
            prediction_table["covariate_prediction"] = baseline
            prediction_table.to_csv(local / "predictions.tsv.gz", sep="\t", index=False)
            history.to_csv(root / "training.tsv", sep="\t", index=False)
            tokens.to_csv(root / "block_manifest.tsv.gz", sep="\t", index=False)
            atomic_json(local / "frozen_protocol.json", dict(contract=contract, fit=fit_meta, nuisance=nuisance,
                learned_ancestry_mixture_in_features=False, shared_weight_fusion=block_receipt.get("shared_weight_fusion")))
            # Test labels are used only after the fitted model and all predictions
            # have been serialized. No threshold or tuning operation follows.
            outcome = frame.loc[test, "y"].to_numpy(float)
            primary, cov = metrics(outcome, predictions[test], binary), metrics(outcome, baseline[test], binary)
            label = "GRID CSx block ridge" if args.module == "blockridge" else "GRID CSx block Transformer"
            performance = dict(trait=trait, method=label, source="2.csx.sh " + args.module,
                official_prscsx=False, stage="test", status="completed", build_N=int(build.sum()),
                tune_N=int(tune.sum()), validation_budget=args.validation_n,
                grouping_status=cohort_meta["grouping_status"], endpoint=cohort_meta["endpoint_definition"],
                test_used_for_fitting=False, **primary)
            performance.update({"covariate_only_"+k:v for k,v in cov.items() if k not in {"N", "N_total", "N_missing_outcome"}})
            pd.DataFrame([performance]).to_csv(root / "performance.tsv", sep="\t", index=False)
            files = [model_path, local / "predictions.tsv.gz", local / "frozen_protocol.json",
                     root / "training.tsv", root / "block_manifest.tsv.gz", root / "performance.tsv"]
            atomic_json(completion, dict(status="completed", contract=contract, block_input=str(array_path),
                block_receipt=block_receipt, cohort_metadata=cohort_meta,
                files={str(p):file_record(p) for p in files}))
            print(f"[{args.module}] {trait}: build={int(build.sum())}, tune={int(tune.sum())}, untouched test={int(test.sum())}; {root}", flush=True)
    finally:
        array._mmap.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("module", choices=["blockprepare", "blockridge", "blocknet"])
    p.add_argument("--analysis-root", "--output-root", default=analysis_default())
    p.add_argument("--traits", "--trait", default="height,ldl,t2dm")
    p.add_argument("--block-config", "--config", default=str(GRID / "config/abm.blocks.yml"))
    p.add_argument("--blocks", type=int, choices=[1,4])
    p.add_argument("--ld-block-files")
    p.add_argument("--ld-block-dir")
    p.add_argument("--score-dir")
    p.add_argument("--cohort")
    p.add_argument("--prepared")
    p.add_argument("--family-col", "--group-id-col", dest="group_col")
    p.add_argument("--validation-n", type=int, default=5000)
    p.add_argument("--covariates", default="age,sex,PC1,PC2")
    p.add_argument("--seed", type=int, default=20261102)
    p.add_argument("--device", choices=["cpu", "cuda", "auto"], default="auto")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--epochs", type=int)
    p.add_argument("--batch-size", type=int)
    p.add_argument("--replace", action="store_true", default=__import__("os").environ.get("GRID_REPLACE", "TRUE") == "TRUE")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)
    if a.threads < 1 or a.validation_n < 4:
        p.error("positive threads and validation-n >= 4 required")
    traits = a.traits.split(",")
    if any(t not in {"height", "ldl", "t2dm"} for t in traits):
        p.error("traits must be height,ldl,t2dm")
    if a.dry_run:
        print(json.dumps(dict(module=a.module, traits=traits, analysis_root=a.analysis_root,
            blocks=a.blocks, ld_block_files=a.ld_block_files, validation_n=a.validation_n,
            features="fixed external weights; no outcome fit" if a.module == "blockprepare" else "block PRSs + declared covariates")))
        return
    config_path, config = effective_config(a)
    import block_inputs
    if a.module == "blockprepare":
        for trait in traits:
            block_inputs.prepare_trait(config_path, trait, replace=a.replace, device=a.device)
    else:
        for trait in traits:
            fit_trait(a, trait, config_path, config)



if __name__ == "__main__":
    main()
