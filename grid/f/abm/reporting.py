#!/usr/bin/env python3
"""Validate complete, frozen-test PRSformer outputs before aggregation."""
from __future__ import annotations

# Benchmark results
import numpy as np
import pandas as pd

from abm.evaluation import prediction_metrics


def performance(run, cache, traits, method, *, exact_traits=False):
    from prsformer import validate_evaluation
    validate_evaluation(run, cache)
    result = pd.read_csv(run / "test_predictions.tsv.gz", sep="\t", dtype={"eid": str})
    cohort = pd.read_csv(cache / "data.tsv.gz", sep="\t", dtype={"eid": str})
    cohort = cohort.loc[cohort.split.eq("test")].set_index("eid")
    if cohort.empty or not cohort.index.is_unique:
        raise ValueError("PRSformer prepared test roster must be nonempty and unique")
    if not result.split.eq("test").all() or result.duplicated(["eid", "trait"]).any():
        raise ValueError("PRSformer predictions must be unique held-out test rows")
    if not set(traits) <= set(result.trait):
        raise ValueError("PRSformer evaluation is missing requested traits")
    if exact_traits and set(result.trait) != set(traits):
        raise ValueError("Independent benchmark evaluation must contain exactly the requested trait")
    rows = []
    for trait in traits:
        group = result.loc[result.trait.eq(trait)].set_index("eid")
        if len(group) != len(cohort) or set(group.index) != set(cohort.index):
            raise ValueError("PRSformer predictions have an incomplete test roster: " + trait)
        expected = pd.to_numeric(cohort.loc[group.index, trait]).to_numpy(float)
        if not np.allclose(group.y, expected, rtol=0, atol=1e-10, equal_nan=True):
            raise ValueError("PRSformer test outcomes differ from prepared inputs: " + trait)
        rows.append(dict(method=method, trait=trait, stage="test",
                         **prediction_metrics(group.y.to_numpy(), group.prediction.to_numpy(),
                                              "binary" if trait == "t2dm" else "continuous")))
    return pd.DataFrame(rows)

# Reporting
"""Normalized aggregate reports; participant records stay in _LOCAL_ONLY."""

import json
from sklearn.linear_model import Ridge, LinearRegression
from sklearn.preprocessing import StandardScaler

from abm.evaluation import evaluate_predictions, make_matched_controls, _effect, _bootstrap_effects
from abm.model import row_loss, assessment_loss
from abm.model import coverage_tag, write_json, COVERAGES


def save_table(frame, path):
    frame.to_csv(path, sep="\t", index=False, float_format="%.17g")


def flatten(obj, prefix=""):
    out = {}
    for key, value in obj.items():
        name = prefix + key
        if isinstance(value, dict):
            out.update(flatten(value, name + "."))
        elif not isinstance(value, (list, tuple)):
            out[name] = value
    return out


def freeze_matches(assessed, local, seed):
    for objective in ("predictability", "benefit"):
        for coverage in COVERAGES:
            key = objective + "_" + coverage_tag(coverage)
            pairs, info = make_matched_controls(assessed, seed=seed, selected_col="selected_" + key,
                                                control_eligible_col="control_eligible")
            pairs.to_parquet(local / ("matched_controls_" + key + ".parquet"), index=False)
            write_json(local / ("matching_" + key + ".json"), info)


def evaluate(assessed, truth, model, output, config, *, label="test"):
    """Verify outcome-blind assignments, then join labels by ID."""
    local = output / "_LOCAL_ONLY" / label
    frame = assessed.merge(truth[["id", "y"]], on="id", how="left", validate="one_to_one")
    tables = {}
    for objective in ("predictability", "benefit"):
        rows = []
        for coverage in COVERAGES:
            key = objective + "_" + coverage_tag(coverage)
            result = evaluate_predictions(
                frame, task=model.heads.task, loss_name=model.heads.loss_name,
                selector="count_only_" + objective, candidate_id=model.manifest["candidate_id"],
                coverage=coverage, selected_col="selected_" + key, control_eligible_col="control_eligible",
                n_boot=config["evaluation"]["bootstrap_replicates"], seed=model.heads.seed,
                min_matched_pairs=config["evaluation"]["min_matched_pairs"],
                min_matched_groups=config["evaluation"]["min_matched_groups_each_side"])
            pairs = result.pop("match_pairs")
            saved = pd.read_parquet(local / ("matched_controls_" + key + ".parquet"))
            pd.testing.assert_frame_equal(pd.DataFrame(pairs, columns=saved.columns), saved, check_dtype=False)
            result["scope"] = "deployment_policy" if objective == "benefit" else "predictability_diagnostic_counterfactual"
            result["validation_status"] = "exploratory_previously_inspected_test"
            result["grouping_status"] = model.manifest.get("grouping_status", "unverified")
            write_json(output / f"{label}_{key}.json", result)
            row = flatten(result)
            row.update(trait=model.manifest["trait"], candidate=model.manifest["candidate_id"],
                       baseline_id=model.manifest["baseline_id"],
                       stage=label, nominal_coverage=coverage, assessment_loss=model.heads.loss_name)
            rows.append(row)
            if objective == "benefit" and coverage == .5:
                policy = []
                for name, role in (("base", model.manifest["baseline_id"]),
                                   ("candidate", model.manifest["candidate_id"]),
                                   ("policy", f"ABM(base={model.manifest['baseline_id']},candidate={model.manifest['candidate_id']})")):
                    policy.append(dict(method=role, trait=model.manifest["trait"], stage=label,
                                       assessment_loss=model.heads.loss_name,
                                       n_groups=frame.loc[frame.y.notna() & (
                                           frame.p_candidate.notna() if name == "candidate" else frame.p_base.notna()),
                                           "group_id"].nunique(),
                                       n_original_groups=result["bootstrap"]["n_original_groups"],
                                       **result["whole_cohort"][name], **result["coverage"],
                                       **flatten({k: result["whole_cohort"][k] for k in
                                                   ("policy_gain_vs_base", "policy_gain_vs_candidate")})))
                tables["policy_performance"] = pd.DataFrame(policy)
        tables[objective + "_curve"] = pd.DataFrame(rows)
    for name, d in tables.items():
        save_table(d, output / (("" if label == "test" else label + "_") + name + ".tsv"))
    save_table(evaluation_calibration(frame, model, config, label),
               output / (label + "_count_calibration.tsv"))
    return frame, tables


def evaluation_calibration(frame, model, config, label):
    """Independent observed bin means; frozen predictions are never refitted."""
    groups = [("all", "all", frame)]
    groups += [("ancestry", str(name), d) for name, d in frame.groupby("ancestry", dropna=False)]
    records = []
    for kind, value, part in groups:
        for count_bin, all_rows in part.groupby("count_bin", dropna=False):
            d = all_rows.loc[all_rows.assessment_available & all_rows.y.notna()]
            base = row_loss(d.y, d.p_base, model.heads.task, model.heads.loss_name)
            gain = base - row_loss(d.y, d.p_candidate, model.heads.task, model.heads.loss_name)
            codes = pd.factorize(d.group_id)[0]
            for objective, actual, predicted_col in (
                    ("expected_loss", base, "expected_loss_base"),
                    ("expected_gain", gain, "expected_gain")):
                expected = d[predicted_col].to_numpy(float)
                finite = np.isfinite(actual) & np.isfinite(expected)
                n_groups = d.loc[finite, "group_id"].nunique()
                observed = _effect(actual, finite)
                bias = _effect(expected - actual, finite)
                enough = n_groups >= config["head_options"]["min_calibration_groups"]
                _bootstrap_effects(codes, [observed, bias],
                    config["evaluation"]["bootstrap_replicates"] if enough else 0,
                    config["head_options"]["seed"])
                records.append(dict(stage=label, stratum=kind, group=value, count_bin=count_bin,
                    objective=objective, assessment_loss=model.heads.loss_name,
                    n_query=len(all_rows), n_observed=int(finite.sum()),
                    n_groups=n_groups, events=float(d.loc[finite, "y"].sum()) if model.heads.task == "binary" else None,
                    n_assessable=int(all_rows.assessment_available.sum()),
                    expected_mean=float(expected[finite].mean()) if finite.any() else None,
                    expected_ci_low=d.loc[finite, predicted_col + "_ci_low"].mean(),
                    expected_ci_high=d.loc[finite, predicted_col + "_ci_high"].mean(),
                    **{"observed_" + k: v for k, v in observed[0].items()},
                    **{"bias_" + k: v for k, v in bias[0].items()},
                    status="available" if enough else "insufficient_evaluation_groups",
                    interval_scope="group_bootstrap_bin_mean_not_individual_error"))
    return pd.DataFrame(records)


def plot_calibration(table, path):
    """Publish only aggregate count-bin means and group uncertainty."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, objective in zip(axes, ("expected_loss", "expected_gain")):
        d = table.loc[table.stratum.eq("all") & table.objective.eq(objective)].sort_values("count_bin")
        x = d.count_bin.to_numpy(float)
        ax.plot(x, d.expected_mean, "o-", label="Frozen expected mean")
        ax.plot(x, d.observed_estimate, "s-", label="Observed mean")
        ax.fill_between(x, d.observed_ci_low.to_numpy(float), d.observed_ci_high.to_numpy(float),
                        alpha=.2, label="95% group bootstrap CI")
        ax.set(xlabel="Frozen count bin", ylabel="Native loss units", title=objective.replace("_", " "))
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def support_diagnostics(frame, binary, pc_edges):
    rows = []
    strata = [("ancestry", frame.ancestry.fillna("UNASSIGNED").astype(str))]
    for name, edges in pc_edges.items():
        strata.append((name + "_development_quintile",
                       pd.Series(np.searchsorted(edges, frame[name], side="right"), index=frame.index).astype(str)))
    for kind, values in strata:
        labels = sorted(set(values))
        if kind == "ancestry":
            labels = sorted(set(labels) | {"EUR", "AFR", "EAS", "SAS"})
        for value in labels:
            d = frame.loc[values.eq(value)]
            row = dict(diagnostic=kind, group=value, n_query=len(d), n_observed=d.y.notna().sum(),
                       n_supported=d.supported.sum(), n_assessable=d.assessment_available.sum(),
                       n_selected=d.selected.sum(), n_switched=d.switched.sum(),
                       n_groups=d.group_id.nunique(), events=d.y.sum() if binary and len(d) else None,
                       assessment_reasons=json.dumps(d.assessment_reason.value_counts().to_dict(), sort_keys=True),
                       status="available" if len(d) else "no_queries")
            for c in ("n_twins", "twins_per_10000", "ess_group", "d_nearest", "d_k_used"):
                x = d[c].replace([np.inf, -np.inf], np.nan)
                row.update({c + "_mean": x.mean(), c + "_q10": x.quantile(.1), c + "_q90": x.quantile(.9)})
            rows.append(row)
    return pd.DataFrame(rows)


def fit_diagnostics(gate, calibration, clinical, task, policy, *, loss_name=None):
    """Same Ridge budget for A/B/C/D; uses independent head calibration."""
    loss_name = assessment_loss(task, loss_name)
    identity = dict(schema_version="abm.diagnostics.1", task=task, loss_name=loss_name)
    def targets(d):
        base = row_loss(d.y, d.p_base, task, loss_name)
        gain = base - row_loss(d.y, d.p_candidate, task, loss_name)
        return np.c_[base, gain]
    gate = gate.loc[gate.supported & np.isfinite(gate.n_twins)].copy()
    calibration = calibration.loc[calibration.supported & np.isfinite(calibration.n_twins)].copy()
    if len(gate) < 30 or len(calibration) < 20:
        return dict(identity, status="unavailable_insufficient_supported_development", models={})
    policy = policy.loc[policy.supported & np.isfinite(policy.n_twins)].copy()
    bundle = dict(identity, status="fitted", clinical=clinical, models={}, thresholds={})
    for key in ("A", "B", "C", "D"):
        xg, xc = diagnostic_design(gate, key, clinical, task), diagnostic_design(calibration, key, clinical, task)
        scaler = StandardScaler().fit(xg)
        m = Ridge(alpha=100.).fit(scaler.transform(xg), targets(gate))
        raw = m.predict(scaler.transform(xc))
        calibrators = [LinearRegression().fit(raw[:, [j]], targets(calibration)[:, j]) for j in range(2)]
        bundle["models"][key] = (scaler, m, calibrators)
        if len(policy):
            raw = m.predict(scaler.transform(diagnostic_design(policy, key, clinical, task)))
            bundle["thresholds"][key] = [float(np.median(calibrators[j].predict(raw[:, [j]]))) for j in range(2)]
        else:
            bundle["thresholds"][key] = [None, None]
    return bundle


def diagnostic_design(d, key, clinical, task):
    p = d.p_base.to_numpy(float)
    parts = [p[:, None], (p ** 2)[:, None], (p ** 3)[:, None]]
    if task == "binary":
        parts.append((p * (1 - p))[:, None])
    if key in ("C", "D"):
        parts.append(d[clinical].to_numpy(float))
    if key in ("B", "D"):
        parts.append(np.log1p(d.n_twins.to_numpy(float))[:, None])
    return np.concatenate(parts, axis=1)


def diagnostic_predictions(bundle, frame):
    """Freeze diagnostic selection using policy-calibration X, before query Y."""
    out = frame.copy()
    valid = out.supported & np.isfinite(out.n_twins)
    for key, (scaler, model, calibrators) in bundle["models"].items():
        raw = model.predict(scaler.transform(diagnostic_design(out.loc[valid], key, bundle["clinical"], bundle["task"]))) if valid.any() else np.empty((0, 2))
        for j, objective in enumerate(("expected_loss", "expected_gain")):
            col = "diagnostic_" + key + "_" + objective
            out[col] = np.nan
            out.loc[valid, col] = calibrators[j].predict(raw[:, [j]]) if len(raw) else np.array([])
            threshold = bundle["thresholds"][key][j]
            out[col + "_selected"] = False if threshold is None else (
                out[col].le(threshold) if j == 0 else out[col].ge(threshold) & out[col].gt(0))
    return out


def evaluate_diagnostics(bundle, frame, config):
    if not bundle["models"]:
        return pd.DataFrame([{"status": bundle["status"]}])
    d = frame.loc[frame.assessment_available & frame.y.notna() & frame.p_candidate.notna()].copy()
    if not len(d):
        return pd.DataFrame([{"status": "unavailable_no_assessable_evaluation_rows"}])
    loss_name = diagnostic_loss(bundle)
    base = row_loss(d.y, d.p_base, bundle["task"], loss_name)
    target = np.c_[base, base - row_loss(d.y, d.p_candidate, bundle["task"], loss_name)]
    _, codes = np.unique(d.group_id, return_inverse=True)
    records, losses, selected_masks = [], {}, {}
    for key, (scaler, m, calibrators) in bundle["models"].items():
        prediction = d[["diagnostic_" + key + "_" + obj for obj in ("expected_loss", "expected_gain")]].to_numpy(float)
        losses[key] = (prediction - target) ** 2
        for j, objective in enumerate(("expected_loss", "expected_gain")):
            chosen = d["diagnostic_" + key + "_" + objective + "_selected"].to_numpy(bool)
            selected_masks[key, j] = chosen
            reduction = losses["A"][:, j] - losses[key][:, j]
            effect = _effect(reduction, np.ones(len(d), bool))
            selection = _effect(target[:, j], chosen, target[:, j], selected_masks["A", j])
            bias = _effect(prediction[:, j] - target[:, j], np.ones(len(d), bool))
            _bootstrap_effects(codes, [effect, selection, bias], config["evaluation"]["bootstrap_replicates"],
                               config["head_options"]["seed"])
            records.append(dict(diagnostic=key, objective=objective, assessment_loss=loss_name,
                                status="available", n=len(d),
                                features="risk_polynomial" + ("+clinical_PCs" if key in ("C", "D") else "") +
                                         ("+log1p_count" if key in ("B", "D") else ""),
                                target_mse=losses[key][:, j].mean(),
                                selection_threshold_role="policy_calibration_X_only",
                                selected_n=int(chosen.sum()), actual_selection_fraction=float(chosen.mean()),
                                selected_target_mean=float(target[chosen, j].mean()) if chosen.any() else None,
                                calibration_bias=float(np.mean(prediction[:, j] - target[:, j])),
                                **{"calibration_bias_" + k: v for k, v in bias[0].items()},
                                **{"selected_target_difference_vs_A_" + k: v for k, v in selection[0].items()},
                                **{"mse_reduction_vs_A_" + k: v for k, v in effect[0].items()}))
    return pd.DataFrame(records)


def diagnostic_loss(bundle):
    if bundle.get("schema_version") != "abm.diagnostics.1" or "loss_name" not in bundle:
        raise ValueError("Diagnostic bundle needs verified loss provenance; unknown schema cannot infer an objective")
    return assessment_loss(bundle["task"], bundle["loss_name"])


def count_incremental_contrasts(bundle, frame, config):
    """Direct paired assessment-target MSE and frozen-selection differences."""
    records = []
    for j, objective in enumerate(("expected_loss", "expected_gain")):
        for comparator, augmented in (("A", "B"), ("C", "D")):
            row = dict(objective=objective, contrast=f"{augmented}_vs_{comparator}",
                comparator=comparator, augmented=augmented, n=0, n_groups=0,
                comparator_target_mse=None, augmented_target_mse=None, mse_reduction=None,
                ci_low=None, ci_high=None, bootstrap_valid_replicates=0,
                assessment_loss=diagnostic_loss(bundle), n_query=len(frame),
                unit="assessment_loss_squared", interpretation="assessment target MSE reduction; positive favours count",
                selection_direction="comparator_minus_augmented" if j == 0 else "augmented_minus_comparator",
                selection_difference=None, selection_ci_low=None, selection_ci_high=None,
                selection_bootstrap_valid_replicates=0, selection_status="unavailable",
                comparator_actual_coverage=None, augmented_actual_coverage=None)
            if not bundle["models"]:
                records.append(dict(row, status=bundle["status"]))
                continue
            loss_name = diagnostic_loss(bundle)
            row["assessment_loss"] = loss_name
            ca, cb = [f"diagnostic_{key}_{objective}" for key in (comparator, augmented)]
            if not {ca, cb, ca+"_selected", cb+"_selected"} <= set(frame):
                raise ValueError("Missing frozen diagnostic predictions/selections; cannot refit on assessment data")
            row.update(comparator_query_coverage=float(frame[ca+"_selected"].mean()) if len(frame) else None,
                       augmented_query_coverage=float(frame[cb+"_selected"].mean()) if len(frame) else None,
                       nominal_selection_coverage=.5,
                       selection_threshold_role="policy_calibration_X_only")
            valid = frame.assessment_available & frame.y.notna()
            base = row_loss(frame.y, frame.p_base, bundle["task"], loss_name)
            target = base if j == 0 else base - row_loss(frame.y, frame.p_candidate, bundle["task"], loss_name)
            valid &= np.isfinite(target) & np.isfinite(frame[ca]) & np.isfinite(frame[cb])
            d = frame.loc[valid]
            target = target[valid]
            la = (target - d[ca].to_numpy(float)) ** 2
            lb = (target - d[cb].to_numpy(float)) ** 2
            ng = d.group_id.nunique()
            row.update(n=len(d), n_groups=ng,
                       events=float(d.y.sum()) if bundle["task"] == "binary" and len(d) else None)
            if len(d) < 2 or ng < 2:
                records.append(dict(row, status="insufficient_paired_rows_or_groups"))
                continue
            _, codes = np.unique(d.group_id.astype(str), return_inverse=True)
            effect = _effect(la - lb, np.ones(len(d), bool))
            a, b = d[ca+"_selected"].to_numpy(bool), d[cb+"_selected"].to_numpy(bool)
            selection = _effect(target, a if j == 0 else b, target, b if j == 0 else a)
            _bootstrap_effects(codes, [effect, selection], config["evaluation"]["bootstrap_replicates"],
                               config["head_options"]["seed"])
            selection_ok = min(d.loc[a, "group_id"].nunique(), d.loc[b, "group_id"].nunique()) >= 2
            row.update(comparator_target_mse=float(la.mean()), augmented_target_mse=float(lb.mean()),
                mse_reduction=effect[0]["estimate"], ci_low=effect[0]["ci_low"], ci_high=effect[0]["ci_high"],
                bootstrap_valid_replicates=effect[0]["bootstrap_valid_replicates"],
                status=effect[0].get("ci_status", effect[0]["status"]),
                comparator_actual_coverage=float(a.mean()), augmented_actual_coverage=float(b.mean()),
                selection_status=selection[0].get("ci_status", "unavailable") if selection_ok else "insufficient_selected_groups")
            if selection_ok:
                row.update(selection_difference=selection[0]["estimate"],
                           selection_ci_low=selection[0]["ci_low"], selection_ci_high=selection[0]["ci_high"],
                           selection_bootstrap_valid_replicates=selection[0]["bootstrap_valid_replicates"])
            records.append(row)
    return pd.DataFrame(records)



# Derived report
"""Verified, read-only reuse of frozen scientific artifacts for derived reports."""
from pathlib import Path
import shutil
import joblib

from abm.model import ABMAssessor, file_record
from abm.grid_inputs import digest
import abm.reporting as reporting

# Audited unversioned writer: its only targets were MSE/Brier. Never infer unknown writers.
UNVERSIONED_REPORTING_SHA = "3b8502df59f14b431ed284215414e48c1ee1e80580a1fb04c0f554aea0a6ecea"


def verify_run(run):
    run = Path(run).resolve()
    complete = json.loads((run / "COMPLETE.json").read_text())
    if complete.get("status") != "completed":
        raise ValueError("Parent run is incomplete")
    if file_record(run / "RUN_MANIFEST.json") != complete["manifest_file_record"]:
        raise ValueError("Parent manifest changed")
    phases = json.loads((run / "PHASES.json").read_text())
    for phase, files in phases["phases"].items():
        for name, expected in files.items():
            path = (run / name).resolve()
            from runtime_state import artifact_matches
            if isinstance(expected, dict):
                expected = dict(expected, path=str(path))  # explicit local report-copy support
            if not path.is_relative_to(run) or not artifact_matches(path, expected):
                raise ValueError("Frozen parent artifact changed: " + name)
    manifest = json.loads((run / "RUN_MANIFEST.json").read_text())
    if phases["contract"] != manifest["contract"]:
        raise ValueError("Parent phase/run contract mismatch")
    return manifest, complete


def report_from(parent, output=None):
    from abm.pipeline import source_manifest
    parent = Path(parent).resolve()
    manifest, complete = verify_run(parent)
    sources = source_manifest()
    config = manifest["contract"]["config"]
    local = parent / "_LOCAL_ONLY"
    model = ABMAssessor.load(local / "assessor.json")
    trait = manifest["contract"]["trait"]
    if config["primary_loss"][trait] != model.heads.loss_name:
        raise ValueError("Configured assessment loss differs from frozen head; refit affected heads/rules first")
    identity = dict(parent=str(parent), parent_manifest_file_record=file_record(parent / "RUN_MANIFEST.json"),
        parent_phases_file_record=file_record(parent / "PHASES.json"),
        parent_completion_file_record=file_record(parent / "COMPLETE.json"), reporting_sources=sources)
    root = Path(output).resolve() if output else parent.parents[3] / "abm_reports" / (parent.name+"_"+digest(identity)[:12])
    if root == parent or root.is_relative_to(parent):
        raise ValueError("Derived reports must be outside the frozen parent run")
    if root.exists():
        raise FileExistsError("Derived output exists; choose a fresh --output-root: " + str(root))
    new_diagnostics = not (local / "diagnostics.joblib").is_file()
    if new_diagnostics:
        # Only the original development partitions supply diagnostic labels.
        truth_dev = pd.read_parquet(local / "truth.parquet")
        features = pd.read_parquet(local / "features.parquet")
        predictions = pd.read_parquet(local / "predictor_predictions.parquet")
        support = pd.read_parquet(local / "support.parquet")
        dev = features.merge(predictions.drop(columns=["group_id", "role"]), on="id",
            validate="one_to_one").merge(support, on="id", validate="one_to_one")
        dev = dev.merge(truth_dev[["id", "y"]], on="id", validate="one_to_one")
        clinical = pd.read_csv(parent / "feature_manifest.tsv", sep="\t").feature.tolist()
        clinical = [name for name in clinical if not name.startswith("csx.")]
        diagnostics = reporting.fit_diagnostics(dev.loc[dev.role.eq("tune_gate")],
            dev.loc[dev.role.eq("head_calibration")], clinical, model.heads.task,
            dev.loc[dev.role.eq("policy_calibration")].drop(columns="y"), loss_name=model.heads.loss_name)
        identity["diagnostics_fit"] = "original_tune_gate_head_calibration_policy_X_only"
    else:
        diagnostics = joblib.load(local / "diagnostics.joblib")
    if diagnostics.get("schema_version") is None:
        prior = manifest["contract"]["sources"].get('grid/f/abm/reporting.py')
        if prior != UNVERSIONED_REPORTING_SHA or model.heads.loss_name != (
                "brier" if model.heads.task == "binary" else "mse"):
            raise ValueError("Unversioned diagnostic objective has no recognized source evidence")
        diagnostics = dict(diagnostics, schema_version="abm.diagnostics.1",
                           loss_name=model.heads.loss_name, task=model.heads.task)
        identity["diagnostic_migration"] = dict(source_file_record=prior, convention="audited_MSE_Brier_writer")
    if reporting.diagnostic_loss(diagnostics) != model.heads.loss_name:
        raise ValueError("Frozen diagnostic/head target mismatch")
    truth = pd.read_parquet(local / "truth.parquet")
    if truth.id.astype(str).duplicated().any():
        raise ValueError("Duplicate truth IDs")
    root.mkdir(parents=True)
    for label in ("test", "audit", "policy_calibration"):
        folder = local / label
        frozen = json.loads((folder / "FROZEN.json").read_text())
        if (frozen["model_file_record"] != file_record(local / "assessor.json") or
            frozen["predictions_file_record"] != file_record(folder / "predictions.parquet") or
            frozen.get("query_outcomes_used") is not False):
            raise ValueError("Frozen prediction identity failed")
        if frozen.get("truth_file_record", file_record(local / "truth.parquet")) != file_record(local / "truth.parquet"):
            raise ValueError("Frozen truth binding failed")
        assessed = pd.read_parquet(folder / "predictions.parquet")
        if new_diagnostics:
            assessed = reporting.diagnostic_predictions(diagnostics, assessed)
        target = truth.loc[truth.role.eq(label)]
        if set(assessed.id.astype(str)) != set(target.id.astype(str)):
            raise ValueError("Prediction/truth roster mismatch")
        destination = root / "_LOCAL_ONLY" / label
        destination.mkdir(parents=True)
        if new_diagnostics:
            assessed.to_parquet(destination / "predictions.parquet", index=False)
            joblib.dump(diagnostics, root / "_LOCAL_ONLY/diagnostics.joblib")
        for path in folder.glob("matched_controls_*.parquet"):
            shutil.copy2(path, destination / path.name)
        evaluated, _ = reporting.evaluate(assessed, target, model, root, config, label=label)
        if label == "test":
            reporting.save_table(reporting.evaluate_diagnostics(diagnostics, evaluated, config),
                                 root / "count_incremental_info.tsv")
            contrast = reporting.count_incremental_contrasts(diagnostics, evaluated, config)
            for key, value in dict(trait=trait, run_id=root.name, baseline_id=model.manifest["baseline_id"],
                                  candidate_id=model.manifest["candidate_id"],
                                  representation=manifest["contract"]["representation"]).items():
                contrast[key] = value
            reporting.save_table(contrast, root / "count_incremental_contrasts.tsv")
    write_json(root / "REPORT_MANIFEST.json", dict(identity, trait=trait, assessment_loss=model.heads.loss_name,
        smoke=complete.get("smoke", False), predictors_refitted=False, query_rules_changed=False,
        generated_private_file_record={str(p.relative_to(root)): file_record(p)
                                  for p in (root / "_LOCAL_ONLY").rglob("*") if p.is_file()},
        original_SOURCE_FROZEN=str(parent / "SOURCE_FROZEN.json")))
    tables = [dict(name=p.stem, path=p.name, file_record=file_record(p), trait=trait,
                   module="abm", run_id=root.name) for p in sorted(root.glob("*.tsv"))]
    write_json(root / "table_manifest.json", tables)
    write_json(root / "COMPLETE.json", dict(status="completed", smoke=complete.get("smoke", False),
        manifest_file_record=file_record(root / "REPORT_MANIFEST.json"),
        table_manifest_file_record=file_record(root / "table_manifest.json")))
    print(json.dumps(dict(output=str(root), parent=str(parent), training=False)))
    return root
