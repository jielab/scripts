#!/usr/bin/env python3
"""Small target-validation fitting and untouched-test performance for PRS methods.

PRS-CSx population scores are combined by four explicit regression coefficients.
The same evaluation infrastructure can calibrate an already computed Disco or
official CSx meta/auto score; it never supplies features to a twin predictor.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import (average_precision_score, log_loss, mean_absolute_error,
                             mean_squared_error, r2_score, roc_auc_score)

from grid_paths import GRID, analysis_default, read_config, relocate
from shared_cohort import (atomic_json, frozen_roster, load_cohort, sha256, table)

POPS = ("AFR", "EAS", "EUR", "SAS")


def fit_glm(x, y, binary):
    """Fit only the explicitly supplied validation labels, without regularization."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(y) < max(10, x.shape[1] + 2) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Insufficient finite small-validation observations")
    if binary and (set(np.unique(y)) != {0., 1.} or min(np.bincount(y.astype(int))) < 2):
        raise ValueError("Binary validation needs cases and controls (at least two of each)")
    if not x.shape[1]:
        p = float(np.mean(y))
        intercept = float(np.log(p / (1 - p))) if binary else p
        return dict(intercept=intercept, beta=[], binary=bool(binary), n=len(y))
    if np.linalg.matrix_rank(np.column_stack([np.ones(len(x)), x])) < x.shape[1] + 1:
        raise ValueError("Target-validation design is rank deficient; coefficients are not identifiable")
    if binary:
        model = LogisticRegression(penalty=None, solver="lbfgs", max_iter=10000, tol=1e-9)
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            model.fit(x, y)
        beta, intercept = model.coef_[0], model.intercept_[0]
    else:
        model = LinearRegression().fit(x, y)
        beta, intercept = model.coef_, model.intercept_
    if not np.isfinite(beta).all() or not np.isfinite(intercept):
        raise ValueError("Nonfinite target-validation fit")
    return dict(intercept=float(intercept), beta=beta.tolist(), binary=bool(binary), n=len(y))


def predict_glm(model, x):
    raw = model["intercept"] + np.asarray(x, float) @ np.asarray(model["beta"], float)
    return expit(raw) if model["binary"] else raw


def standardize_scores(x, roles, mode="per_split"):
    """PRS-CSx's optional per-dataset standardization reads no outcome labels.

The validation-frozen option supports prospective deployment without using a
batch of query genotypes. Nuisance covariates always use validation statistics.
"""
    x, roles = np.asarray(x, float), np.asarray(roles, str)
    fit = roles == "target_validation"
    if not fit.any() or not np.isfinite(x).all():
        raise ValueError("Score matrix must be finite on the full frozen cohort")
    mean, sd = x[fit].mean(0), x[fit].std(0)
    if np.any(sd <= 1e-12):
        raise ValueError("Constant validation PRS")
    out = (x - mean) / sd
    record = {"target_validation": {"mean": mean.tolist(), "sd": sd.tolist()}}
    if mode == "per_split":
        # Only the evaluated test batch has an independent PRS distribution.
        # Excluded/unused development rows do not define a prediction batch and
        # cannot make an otherwise valid fit fail when there is just one such
        # row. Their predictions retain the validation-frozen transformation.
        for role in ("test",):
            mask = roles == role
            if not mask.any():
                continue
            mu, sigma = x[mask].mean(0), x[mask].std(0)
            if np.any(sigma <= 1e-12):
                raise ValueError(f"Constant PRS in {role}; use validation-frozen scaling if appropriate")
            out[mask] = (x[mask] - mu) / sigma
            record[role] = {"mean": mu.tolist(), "sd": sigma.tolist()}
    elif mode != "validation":
        raise ValueError("Unknown score standardization")
    return out, dict(mode=mode, labels_used=False, statistics=record)


def metrics(y, prediction, binary):
    y, prediction = np.asarray(y, float), np.asarray(prediction, float)
    observed = np.isfinite(y)
    if not np.isfinite(prediction).all():
        raise ValueError("Prediction must cover every test participant")
    y, p = y[observed], prediction[observed]
    result = dict(N_total=len(prediction), N=len(y), N_missing_outcome=int((~observed).sum()))
    if len(y) < 2:
        raise ValueError("Test has fewer than two observed outcomes")
    if binary:
        if np.any((p < 0) | (p > 1)):
            raise ValueError("Binary predictions must be probabilities")
        both = len(np.unique(y)) == 2
        result.update(cases=int(y.sum()), prevalence=float(y.mean()),
            auc=float(roc_auc_score(y, p)) if both else None,
            auprc=float(average_precision_score(y, p)) if both else None,
            brier=float(np.mean((p-y)**2)), logloss=float(log_loss(y, np.clip(p, 1e-15, 1-1e-15), labels=[0, 1])))
    else:
        mse = float(mean_squared_error(y, p))
        correlation = float(np.corrcoef(y, p)[0, 1]) if np.std(y) > 0 and np.std(p) > 0 else None
        result.update(mse=mse, rmse=float(np.sqrt(mse)), mae=float(mean_absolute_error(y, p)),
            prediction_r2=float(r2_score(y, p)), pearson_r=correlation,
            pearson_r2=correlation**2 if correlation is not None else None)
    return result


def fit_score_model(scores, side, y_validation, roles, binary, standardization="per_split"):
    """No test-y argument exists; downstream evaluation cannot change this fit."""
    fit = np.asarray(roles) == "target_validation"
    z, scale = standardize_scores(scores, roles, standardization)
    side = np.asarray(side, float)
    if not np.isfinite(side).all():
        raise ValueError("Declared covariates must be finite on the entire frozen cohort")
    mu, sd = side[fit].mean(0), side[fit].std(0)
    sd = np.where(sd > 1e-12, sd, 1.)
    zside = (side - mu) / sd
    observed = np.isfinite(y_validation)
    if len(y_validation) != int(fit.sum()):
        raise ValueError("Fit API receives validation labels only")
    full = fit_glm(np.column_stack([z[fit], zside[fit]])[observed], np.asarray(y_validation)[observed], binary)
    nuisance = fit_glm(zside[fit][observed], np.asarray(y_validation)[observed], binary)
    final_prs = z @ np.asarray(full["beta"][:z.shape[1]])
    prediction = predict_glm(full, np.column_stack([z, zside]))
    baseline = predict_glm(nuisance, zside)
    frozen = dict(model=full, nuisance=nuisance, score_scaling=scale,
                  covariate_scaling=dict(mean=mu.tolist(), sd=sd.tolist()),
                  formula="intercept + standardized_PRS @ beta_PRS + standardized_covariates @ beta_covariates")
    return final_prs, prediction, baseline, frozen


def score_columns(frame, method, models="populations,auto,meta"):
    if method == "disco":
        column = next((c for c in ("disco", "PRS", "DiscoDivas") if c in frame), None)
        if column is None:
            raise ValueError("Disco score column missing")
        return [("disco", [column], "DiscoDivas pinned interpolation + target calibration")]
    result = []
    for requested in models.split(","):
        if requested == "populations":
            columns = [next((c for c in (f"csx.{pop}", f"CSX_{pop}", f"csx_{pop}") if c in frame), None) for pop in POPS]
            if any(c is None for c in columns):
                raise ValueError("Four population-specific PRS-CSx scores are required")
            result.append(("populations", columns, "PRS-CSx validation combination"))
        elif requested in {"auto", "meta"}:
            column = next((c for c in (f"csx.{requested}", f"CSX_{requested}") if c in frame), None)
            if column is not None:
                result.append((requested, [column], f"PRS-CSx-{requested} + target calibration"))
        else:
            raise ValueError("--models accepts populations,auto,meta")
    if not result:
        raise ValueError("No requested PRS-CSx score found")
    return result


def cohort_source(args):
    """Resolve the producer cohort while preserving explicit input precedence."""
    if args.cohort or args.prepared:
        return args.cohort
    config = read_config(getattr(args, "block_config", None) or GRID / "config/abm.blocks.json",
                         args.analysis_root)
    return Path(config["prsformer_cache"]) / "data.tsv.gz"


def evaluate(args, trait, method="csx"):
    from kinship import analysis_guard
    from resume import run_lock
    root = Path(args.analysis_root)
    score_name = "1csx.scores.rds" if method == "csx" else "2disco.scores.rds"
    score_path = Path(args.scores.format(trait=trait)) if args.scores else Path(args.score_dir) / trait / score_name
    out = (root / "csx" / trait / "targetfit") if method == "csx" else (root / "disco" / trait)
    out.mkdir(parents=True, exist_ok=True)
    with analysis_guard(root), run_lock(out):
        frame, metadata = load_cohort(root, trait, cohort_source(args), args.prepared, args.group_col)
        roster, roster_meta = frozen_roster(frame, root, args.validation_n, args.seed)
        frame = frame.merge(roster[["id", "role", "block_role"]], on="id", validate="one_to_one", sort=False)
        score = table(score_path).set_index("id")
        if not set(frame.id) <= set(score.index):
            raise ValueError("Score table does not cover the complete frozen cohort")
        score = score.loc[frame.id].reset_index()
        columns = score_columns(score, method, args.models)
        covariates = [] if args.covariates.lower() in {"none", ""} else args.covariates.split(",")
        if set(covariates) - set(frame):
            raise ValueError("Missing covariates: " + ",".join(sorted(set(covariates)-set(frame))))
        roles = frame.role.to_numpy(str)
        fit, test = roles == "target_validation", roles == "test"
        binary = trait == "t2dm"
        side = frame[covariates].to_numpy(float)
        local = out / "_LOCAL_ONLY"
        local.mkdir(exist_ok=True)
        roster.to_csv(local / "roles.tsv.gz", sep="\t", index=False)
        fitted, predictions, coefficients = {}, frame[["id", "group_id", "role", "ancestry"]].copy(), []
        # Only validation outcomes are passed to the fitting function. Test
        # labels are first used below, after frozen models/predictions are saved.
        for name, features, label in columns:
            genetic, pred, baseline, model = fit_score_model(score[features].to_numpy(float), side,
                frame.loc[fit, "y"].to_numpy(float), roles, binary, args.prs_standardization)
            fitted[name] = dict(model, method=label, score_columns=features, covariates=covariates)
            predictions[f"{name}.final_prs"] = genetic
            predictions[f"{name}.prediction"] = pred
            predictions[f"{name}.covariate_prediction"] = baseline
            terms = features + covariates
            coefficients.append(dict(model=name, term="intercept", kind="intercept", beta=model["model"]["intercept"]))
            for index, (term, beta) in enumerate(zip(terms, model["model"]["beta"])):
                coefficients.append(dict(model=name, term=term, kind="population_prs" if name == "populations" and index < 4 else
                    "score" if index < len(features) else "covariate", beta=beta))
        atomic_json(local / "frozen_models.json", fitted)
        predictions.to_csv(local / "final_prs.tsv.gz", sep="\t", index=False)
        coefficient_table = pd.DataFrame(coefficients)
        coefficient_table.to_csv(out / "coefficients.tsv", sep="\t", index=False)
        coefficient_table.loc[coefficient_table.kind.eq("population_prs")].to_csv(out / "population_betas.tsv", sep="\t", index=False)
        rows, strata = [], []
        # Evaluation only: no coefficients, standardization, model selection or
        # eligibility threshold can be changed after seeing these labels.
        ytest = frame.loc[test, "y"].to_numpy(float)
        truth_path = local / "test_truth.tsv.gz"
        frame.loc[test, ["id", "y"]].to_csv(truth_path, sep="\t", index=False)
        for name, _, label in columns:
            p = predictions.loc[test, f"{name}.prediction"].to_numpy(float)
            baseline = predictions.loc[test, f"{name}.covariate_prediction"].to_numpy(float)
            primary, nuisance = metrics(ytest, p, binary), metrics(ytest, baseline, binary)
            row = dict(trait=trait, method=label, model=name, source="2.csx.sh" if method == "csx" else "3.disco.sh",
                stage="test", status="completed", training_role="target_validation", test_used_for_fitting=False,
                validation_N=int(fit.sum()), validation_observed_N=int(np.isfinite(frame.loc[fit, "y"]).sum()),
                score_standardization=args.prs_standardization, covariates=",".join(covariates),
                endpoint=metadata["endpoint_definition"], grouping_status=metadata["grouping_status"], **primary)
            row.update({"covariate_only_" + k: v for k, v in nuisance.items() if k not in {"N", "N_total", "N_missing_outcome"}})
            if not binary:
                row["incremental_prediction_r2"] = primary["prediction_r2"] - nuisance["prediction_r2"]
                row["partial_prediction_r2"] = 1 - primary["mse"] / nuisance["mse"] if nuisance["mse"] > 0 else None
            else:
                row["delta_auc"] = primary["auc"] - nuisance["auc"] if primary["auc"] is not None else None
            rows.append(row)
            test_ancestry = frame.loc[test, "ancestry"].to_numpy(str)
            for ancestry in sorted(set(test_ancestry)):
                mask = test_ancestry == ancestry
                if np.isfinite(ytest[mask]).sum() >= 2:
                    strata.append(dict(trait=trait, model=name, ancestry=ancestry, **metrics(ytest[mask], p[mask], binary)))
        pd.DataFrame(rows).to_csv(out / "performance.tsv", sep="\t", index=False)
        pd.DataFrame(strata).to_csv(out / "performance_by_ancestry.tsv", sep="\t", index=False)
        bound = [score_path, Path(metadata["source"]), Path(roster_meta["path"]), local / "frozen_models.json",
                 local / "final_prs.tsv.gz", local / "roles.tsv.gz", truth_path, out / "performance.tsv", out / "coefficients.tsv"]
        receipt = dict(schema="grid.target_benchmark.1", status="completed", method=method, trait=trait,
            metadata=metadata, roster=roster_meta, scores=str(score_path.resolve()),
            implementation_sha256=sha256(__file__), models=list(fitted),
            model_labels={name: value["method"] for name, value in fitted.items()},
            predictions_file=str((local / "final_prs.tsv.gz").resolve()), test_truth_file=str(truth_path.resolve()),
            endpoint=metadata["endpoint_definition"], task="binary" if binary else "continuous",
            test_labels_used_only_after_prediction_freeze=True,
            files={str(p.resolve()): sha256(p) for p in bound})
        atomic_json(out / "performance.json", receipt)
        atomic_json(out / "COMPLETE.json", receipt)
        print(f"[{method}] {trait}: fitted N={int(fit.sum())}; untouched test N={int(test.sum())}; {out / 'performance.tsv'}", flush=True)
        return pd.DataFrame(rows), receipt


def parser(method="csx"):
    import os
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis-root", "--output-root", default=analysis_default())
    p.add_argument("--traits", "--trait", default="height,ldl,t2dm")
    p.add_argument("--score-dir", default=os.environ.get("GRID_SCORE_DIR", relocate("/mnt/f/gen/ukb/pgs")))
    p.add_argument("--scores", "--score-file", help="TSV/RDS; {trait} is expanded per trait")
    p.add_argument("--cohort", help="Optional standalone TSV: id, group_id, split, trait/y and covariates")
    p.add_argument("--prepared", help="Existing trusted GRID prepared.joblib alternative")
    p.add_argument("--block-config", default=str(GRID / "config/abm.blocks.json"),
                   help="Outcome-free producer config; resolves the shared cohort cache when no explicit cohort/prepared input is supplied")
    p.add_argument("--family-col", "--group-id-col", dest="group_col", help="Family/component ID column (not ancestry)")
    p.add_argument("--validation-n", type=int, default=5000, help="Maximum whole-group target validation size; default 5000")
    p.add_argument("--seed", type=int, default=20261102)
    p.add_argument("--covariates", default="age,sex,PC1,PC2")
    p.add_argument("--prs-standardization", choices=["per_split", "validation"], default="per_split")
    p.add_argument("--models", default="populations,auto,meta")
    p.add_argument("--dry-run", action="store_true")
    return p


def main(argv=None, method="csx"):
    p = parser(method)
    a = p.parse_args(argv)
    for trait in a.traits.split(","):
        if trait not in {"height", "ldl", "t2dm"}:
            p.error("traits must be height,ldl,t2dm")
    if a.dry_run:
        print(json.dumps(dict(module="targetfit" if method == "csx" else "evaluate", traits=a.traits,
            analysis_root=a.analysis_root, validation_n=a.validation_n, test_used_for_fitting=False)))
        return
    if not a.cohort and not a.prepared and not Path(cohort_source(a)).is_file():
        from shared_cohort import ensure_inputs
        from types import SimpleNamespace
        ensure_inputs(SimpleNamespace(analysis_root=a.analysis_root, traits=a.traits, threads=4, block_config=a.block_config,
            complete_features=False, dry_run=False))
    for trait in a.traits.split(","):
        evaluate(a, trait, method)


if __name__ == "__main__":
    main()
