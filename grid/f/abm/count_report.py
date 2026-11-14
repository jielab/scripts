#!/usr/bin/env python3
"""Evaluation only after label-free memberships/predictions have been frozen."""
from __future__ import annotations

# Report

import html
from pathlib import Path

import numpy as np
import pandas as pd


def metrics(y, prediction, task):
    y, p = np.asarray(y, float), np.asarray(prediction, float)
    valid = np.isfinite(y) & np.isfinite(p)
    y, p = y[valid], p[valid]
    out = dict(n_evaluable=len(y), primary_loss=np.nan, mae=np.nan)
    if not len(y):
        return out
    out.update(primary_loss=float(np.mean((y-p)**2)), mae=float(np.mean(np.abs(y-p))))
    if task == "binary":
        from sklearn.metrics import roc_auc_score, average_precision_score
        both = len(np.unique(y)) == 2
        out.update(n_cases=int(y.sum()), event_fraction=float(y.mean()),
                   auroc=float(roc_auc_score(y, p)) if both else np.nan,
                   auprc=float(average_precision_score(y, p)) if both else np.nan)
    else:
        variance = float(np.mean((y-y.mean())**2))
        out.update(rmse=float(np.sqrt(out["primary_loss"])),
                   r2=1-out["primary_loss"]/variance if variance else np.nan,
                   correlation=float(np.corrcoef(y, p)[0, 1]) if np.std(y) and np.std(p) else np.nan)
    return out


def bootstrap_mean(values, groups, repetitions, seed):
    v = np.asarray(values, float)
    valid = np.isfinite(v)
    v, groups = v[valid], np.asarray(groups, str)[valid]
    if not len(v):
        return dict(estimate=np.nan, low=np.nan, high=np.nan, groups=0)
    _, inverse = np.unique(groups, return_inverse=True)
    n = np.bincount(inverse)
    s = np.bincount(inverse, weights=v)
    out = dict(estimate=float(v.mean()), low=np.nan, high=np.nan, groups=len(n))
    if len(n) < 2 or repetitions < 2:
        return out
    rng = np.random.default_rng(seed)
    boot = np.empty(repetitions)
    for i in range(repetitions):
        draws = rng.integers(0, len(n), len(n))
        boot[i] = s[draws].sum()/n[draws].sum()
    out.update(low=float(np.quantile(boot, .025)), high=float(np.quantile(boot, .975)))
    return out


def same_risk_comparison(frame, selected):
    """Overlap-standardized descriptive comparison in frozen prediction bins.

    No query outcome chooses matches or weights. Rejected controls sharing a
    group with selected queries are excluded. Coarse bins need not control
    every determinant of error; the result is a descriptive comparison.
    """
    good = frame.y.notna().to_numpy() & np.isfinite(frame.p_population)
    selected = np.asarray(selected, bool)
    selected_groups = set(frame.loc[selected, "group_id"])
    independent_controls = ~frame.group_id.isin(selected_groups).to_numpy()
    cells = []
    for risk_bin in sorted(frame.risk_bin.unique()):
        members = good & frame.risk_bin.eq(risk_bin).to_numpy()
        a, b = frame.loc[members & selected], frame.loc[members & ~selected & independent_controls]
        if not len(a) or not len(b):
            continue
        cells.append((min(len(a), len(b)), len(a), len(b),
                      float(np.mean((a.y-a.p_population)**2)),
                      float(np.mean((b.y-b.p_population)**2)),
                      a.group_id.nunique(), b.group_id.nunique()))
    n_selected = int(np.sum(good & selected)); n_other = int(np.sum(good & ~selected))
    if not cells:
        return dict(standardized_selected_loss=np.nan, standardized_other_loss=np.nan,
                    standardized_difference=np.nan, overlap_selected_fraction=0.,
                    overlap_other_fraction=0., overlap_bins=0, comparison_status="no_common_risk_bins")
    values = np.asarray(cells)
    weights = values[:, 0]/values[:, 0].sum()
    a, b = float(weights@values[:, 3]), float(weights@values[:, 4])
    enough = values[:, 0].sum() >= 20 and values[:, 5].sum() >= 10 and values[:, 6].sum() >= 10
    return dict(standardized_selected_loss=a, standardized_other_loss=b,
                standardized_difference=a-b,
                overlap_selected_fraction=float(values[:, 1].sum()/n_selected) if n_selected else 0.,
                overlap_other_fraction=float(values[:, 2].sum()/n_other) if n_other else 0.,
                overlap_bins=len(cells),
                comparison_status="descriptive_coarse_risk_standardization" if enough else "insufficient_overlap_counts")


def evaluate(frame, *, task, trait, role, repetitions=200, seed=20261102):
    if frame.id.duplicated().any():
        raise ValueError("Evaluation requires unique frozen query IDs")
    observed = frame.y.notna()
    if task == "binary" and not np.isin(frame.loc[observed, "y"], [0., 1.]).all():
        raise ValueError("Binary evaluation outcomes must be 0/1")
    rows = []
    subsets = [("all_with_explicit_fallback", np.ones(len(frame), bool), "p_population"),
               ("global_donor_prior_all", np.ones(len(frame), bool), "p_global_prior"),
               ("supported_count_only", frame.supported.to_numpy(bool), "p_abm"),
               ("unsupported_explicit_fallback", ~frame.supported.to_numpy(bool), "p_population")]
    subsets += [(state, frame.evidence_state.eq(state).to_numpy(), "p_population") for state in
                ("supported_lower_dispersion", "supported_heterogeneous")]
    for label, selected, prediction in subsets:
        d = frame.loc[selected]
        m = metrics(d.y, d[prediction], task)
        delta = (d.y-d[prediction])**2-(d.y-d.p_global_prior)**2
        ci = bootstrap_mean(delta, d.group_id, repetitions, seed)
        rows.append(dict(trait=trait, role=role, population=label, task=task,
                         loss_name="brier" if task == "binary" else "mse",
                         n_roster=len(frame), n_selected_roster=len(d),
                         actual_coverage=len(d)/len(frame), n_missing_y=int(d.y.isna().sum()),
                         delta_loss_vs_global_prior=ci["estimate"], delta_low=ci["low"], delta_high=ci["high"],
                         bootstrap_unit="declared_group", n_bootstrap_groups=ci["groups"], **m))
    risk = []
    total_cases = float(frame.loc[observed, "y"].sum()) if task == "binary" else np.nan
    total_controls = int(observed.sum())-total_cases if task == "binary" else np.nan
    for column in (c for c in frame if c.startswith("selected_")):
        selected = frame[column].to_numpy(bool)
        d = frame.loc[selected]
        nominal = float(column.removeprefix("selected_").replace("p", "."))
        m = metrics(d.y, d.p_population, task)
        ci = bootstrap_mean((d.y-d.p_population)**2, d.group_id, repetitions, seed)
        selected_cases = float(d.y.sum()) if task == "binary" else np.nan
        selected_controls = int(d.y.notna().sum())-selected_cases if task == "binary" else np.nan
        risk.append(dict(trait=trait, role=role, nominal_calibration_coverage=nominal,
                         actual_test_coverage=float(selected.mean()), n_selected_roster=int(selected.sum()),
                         primary_loss_low=ci["low"], primary_loss_high=ci["high"],
                         case_retention=selected_cases/total_cases if total_cases else np.nan,
                         control_retention=selected_controls/total_controls if total_controls else np.nan,
                         **m, **same_risk_comparison(frame, selected)))
    return pd.DataFrame(rows), pd.DataFrame(risk)


def subgroup_table(frame, *, task, trait, role):
    records = []
    for ancestry, d in frame.groupby("ancestry", dropna=False):
        records.append(dict(trait=trait, role=role, ancestry=str(ancestry), n_roster=len(d),
                            support_coverage=float(d.supported.mean()),
                            mean_twins=float(d.n_twins.mean()),
                            **metrics(d.y, d.p_population, task)))
    return pd.DataFrame(records)


def plot_risk(risk, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for role, d in risk.groupby("role"):
        axes[0].plot(d.actual_test_coverage, d.primary_loss, "o-", label=role)
        if d.case_retention.notna().any():
            axes[1].plot(d.actual_test_coverage, d.case_retention, "o-", label=role+" cases")
            axes[1].plot(d.actual_test_coverage, d.control_retention, "s--", label=role+" controls")
        else:
            axes[1].plot(d.actual_test_coverage, d.overlap_selected_fraction, "o-", label=role+" selected overlap")
            axes[1].plot(d.actual_test_coverage, d.overlap_other_fraction, "s--", label=role+" rejected overlap")
    axes[0].set(xlabel="Actual retained fraction of frozen roster", ylabel="Observed squared loss", title="Support-ranked retention")
    axes[1].set(xlabel="Actual retained fraction of frozen roster", ylabel="Fraction", ylim=(-.02, 1.02), title="Composition / common-risk overlap")
    for ax in axes:
        ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(output, dpi=160); plt.close(fig)


def write_trait_report(output, model, performance, risk, subgroups):
    caveats = (
        "The prediction is the arithmetic mean of every eligible training donor inside one frozen genetic-profile radius. "
        "Binary predictions are observed case counts divided by all donor counts. No fitted regression, CountHeads, "
        "risk-bin error head, or test-outcome donor pruning is used. "
        "Unsupported people receive no ABM prediction; an explicitly labelled global donor mean is retained for population evaluation. "
        "More support is not evidence of lower personal error. The displayed intervals and donor dispersion are descriptive; "
        "they do not guarantee individual or conditional coverage. Coarse same-predicted-risk standardization is descriptive, "
        "and its overlap fractions must be read with the difference. Previously inspected tests remain exploratory."
    )
    body = [f"<h1>GRID ABM: {html.escape(model['trait'])}</h1>",
            "<p>Source: grid.abm.1; current completion receipt only. Historical blocknet and CountHeads reports are separate.</p>",
            f"<p>{html.escape(caveats)}</p>",
            f"<p>Group evidence: {html.escape(str(model['metadata'].get('grouping_status', 'unverified')))}; "
            f"GWAS overlap: {html.escape(str(model['metadata'].get('gwas_overlap_status', 'unverified')))}.</p>",
            "<h2>Population and supported prediction</h2>", performance.to_html(index=False, border=0),
            "<h2>Frozen selection thresholds, actual coverage</h2>", risk.to_html(index=False, border=0),
            '<img src="risk_coverage.png" alt="Risk coverage and composition" style="max-width:100%">',
            "<h2>Ancestry audit</h2>", subgroups.to_html(index=False, border=0),
            "<p>Complete donor memberships, reference outcomes, query identities, and examples are in the local-only directory and are not embedded in this report.</p>"]
    Path(output).write_text("<!doctype html><html><head><meta charset='utf-8'><title>GRID ABM</title>"
                           "<style>body{font:15px system-ui;margin:32px;max-width:1500px}table{font-size:11px;border-collapse:collapse}"
                           "td,th{padding:5px;border-bottom:1px solid #ddd}h2{margin-top:32px}</style></head><body>"+
                           "\n".join(body)+"</body></html>")

# Benchmarks
"""Optional, exact-person comparison against current independent benchmarks.

Only the new completion contract is accepted. Missing/incompatible methods are
reported, never substituted by historical aggregate tables or a smaller common
intersection. This module cannot fit, calibrate or select a benchmark model.
"""

import json


from abm.count_inputs import file_record, read_table, ids_hash


def benchmark_receipts(root, trait):
    root = Path(root)
    return {"csx": root/"csx"/trait/"targetfit/COMPLETE.json",
            "disco": root/"disco"/trait/"COMPLETE.json",
            "prsformer": root/"qc/prsformer"/trait/"COMPLETE.json"}


def endpoint_identity(value):
    if value in {"t2dm", "baseline_t2dm", "baseline_t2dm_from_Yr2e_Yt2e"}:
        return "baseline_t2dm"
    return str(value)


def bound_file(receipt, key, receipt_path):
    name = receipt.get(key)
    if not isinstance(name, str) or not name:
        raise ValueError("Missing explicit frozen snapshot path: "+key)
    path = Path(name)
    if not path.is_absolute():
        path = receipt_path.parent/path
    path = path.resolve()
    expected = receipt.get("files", {}).get(str(path))
    if not expected:
        raise ValueError("Snapshot is not producer hash-bound: "+key)
    if not path.is_file() or file_record(path) != expected:
        raise ValueError("Bound snapshot is missing or changed: "+key)
    return path


def load_benchmark(path, *, method, trait, task, endpoint, count_frame, count_csx_score_file_record=None):
    path = Path(path)
    before = file_record(path)
    receipt = json.loads(path.read_text())
    if receipt.get("schema") != "grid.target_benchmark.1" or receipt.get("status") != "completed":
        raise ValueError("Unsupported/incomplete benchmark receipt; historical results are not imported")
    if receipt.get("method") != method or receipt.get("trait") != trait or receipt.get("task") != task:
        raise ValueError("Benchmark method, trait or task identity differs")
    if receipt.get("test_labels_used_only_after_prediction_freeze") is not True:
        raise ValueError("Benchmark has no verified declaration of prediction freeze before evaluation")
    if endpoint_identity(receipt.get("endpoint")) != endpoint_identity(endpoint):
        raise ValueError("Endpoint definitions differ")
    # The current producer binds coefficient/model, roster, source and results.
    # Verify every declared artifact rather than trusting an aggregate metric.
    if not receipt.get("files"):
        raise ValueError("Benchmark receipt has no artifact bindings")
    for filename, expected in receipt["files"].items():
        source = Path(filename)
        if not source.is_absolute():
            raise ValueError("Benchmark receipt artifact paths must be absolute")
        if not source.is_file() or file_record(source) != expected:
            raise ValueError("Benchmark producer-bound artifact missing or changed: "+source.name)
    if method == "csx" and count_csx_score_file_record is not None:
        score_path = receipt.get("scores")
        if not score_path or receipt["files"].get(str(Path(score_path).resolve())) != count_csx_score_file_record:
            raise ValueError("CSx score provenance differs from the count geometry's producer")
    prediction_path = bound_file(receipt, "predictions_file", path)
    truth_path = bound_file(receipt, "test_truth_file", path)
    predictions, truth = read_table(prediction_path), read_table(truth_path)
    if "y" in predictions:
        raise ValueError("Benchmark prediction snapshot includes outcomes; require the label-free frozen snapshot")
    if "role" in predictions:
        predictions = predictions.loc[predictions.role.eq("test")].copy()
    if "id" not in predictions or "id" not in truth or "y" not in truth:
        raise ValueError("Benchmark snapshots require id and separate truth y")
    if predictions.id.isna().any() or truth.id.isna().any() or predictions.id.duplicated().any() or truth.id.duplicated().any():
        raise ValueError("Benchmark snapshots contain missing or duplicate IDs")
    expected_ids = set(count_frame.id)
    if set(predictions.id) != expected_ids or set(truth.id) != expected_ids:
        raise ValueError("Test ID roster differs; no smaller intersection is compared")
    predictions = predictions.set_index("id").loc[count_frame.id]
    truth = truth.set_index("id").loc[count_frame.id]
    y = pd.to_numeric(truth.y, errors="raise").to_numpy(float)
    if not np.array_equal(y, count_frame.y.to_numpy(float), equal_nan=True):
        raise ValueError("Test labels or missing-label pattern differ")
    if "group_id" in predictions and not np.array_equal(predictions.group_id.astype(str), count_frame.group_id.astype(str)):
        raise ValueError("Declared test group mapping differs")
    if not receipt.get("models") or not isinstance(receipt["models"], list):
        raise ValueError("Benchmark must enumerate its frozen models")
    values = {}
    for model in receipt["models"]:
        column = model+".prediction"
        if column not in predictions:
            raise ValueError("Frozen benchmark prediction column missing: "+column)
        p = pd.to_numeric(predictions[column], errors="raise").to_numpy(float)
        if not np.isfinite(p).all() or (task == "binary" and np.any((p < 0) | (p > 1))):
            raise ValueError("Benchmark must predict every test person on the declared outcome scale")
        values[model] = p
    # Detect a concurrent writer between verification and reading.
    for key in ("predictions_file", "test_truth_file"):
        bound_file(receipt, key, path)
    if file_record(path) != before:
        raise ValueError("Benchmark completion receipt changed while reading")
    return values, receipt, before


def compare_current(root, traits, *, repetitions=200, seed=20261102):
    rows, statuses = [], []
    for trait in traits:
        folder = Path(root)/"tmp/abm"/trait
        complete = folder/"EVALUATE_COMPLETE.json"
        if not complete.is_file():
            for method in benchmark_receipts(root, trait):
                statuses.append(dict(trait=trait, benchmark_family=method, status="count_result_missing", reason="No completed current count evaluation"))
            continue
        # The caller verifies the complete count receipt chain before this read.
        model = json.loads((folder/"model.json").read_text())
        count = read_table(folder/"_LOCAL_ONLY/test.predictions.tsv.gz")
        truth = read_table(folder/"_LOCAL_ONLY/test.truth.tsv.gz")
        count = count.merge(truth, on="id", how="left", validate="one_to_one", sort=False)
        endpoint = model["metadata"].get("endpoint_definition", trait)
        for method, path in benchmark_receipts(root, trait).items():
            status = dict(trait=trait, benchmark_family=method, completion_receipt=str(path),
                          n_count_test_roster=len(count), count_test_ids_sha256=ids_hash(count.id))
            if not path.is_file():
                statuses.append(dict(status, status="missing_current_compatible_receipt",
                    reason="Independent benchmark has no new completion contract; legacy results were not read"))
                continue
            try:
                values, receipt, receipt_hash = load_benchmark(path, method=method, trait=trait,
                    task=model["task"], endpoint=endpoint, count_frame=count,
                    count_csx_score_file_record=model["metadata"].get("whole_csx_score_file_record"))
            except (ValueError, FileNotFoundError, KeyError, TypeError, OSError) as exc:
                statuses.append(dict(status, status="not_comparable", reason=str(exc)))
                continue
            statuses.append(dict(status, status="comparable_exact_test_roster_and_labels",
                                 reason="No intersection filtering; missing-label positions and optional group mapping also match",
                                 benchmark_receipt_file_record=receipt_hash, models=",".join(values)))
            for name, prediction in values.items():
                label = receipt.get("model_labels", {}).get(name) or receipt.get("method_label") or method+":"+name
                for population, selected in (("all_with_explicit_fallback", np.ones(len(count), bool)),
                                              ("count_supported_same_people", count.supported.to_numpy(bool))):
                    selected = np.asarray(selected, bool)
                    sample = count.loc[selected]
                    bench_prediction = prediction[selected]
                    count_metrics = metrics(sample.y, sample.p_population, model["task"])
                    benchmark_metrics = metrics(sample.y, bench_prediction, model["task"])
                    delta = (sample.y.to_numpy(float)-sample.p_population.to_numpy(float))**2-(sample.y.to_numpy(float)-bench_prediction)**2
                    uncertainty = bootstrap_mean(delta, sample.group_id, repetitions, seed)
                    record = dict(trait=trait, benchmark_family=method, benchmark_model=name, benchmark_label=label,
                        population=population, status="comparable" if len(sample) else "no_supported_queries",
                        n_outer_roster=len(count), n_same_people=len(sample),
                        actual_count_coverage=len(sample)/len(count),
                        n_evaluable=count_metrics["n_evaluable"],
                        count_donor_n=model["lineage"]["donor"]["n"],
                        csx_source_identity=("matched_count_producer_score" if method == "csx" and model["metadata"].get("whole_csx_score_file_record")
                                             else "not_applicable_or_external_export"),
                        benchmark_validation_n=receipt.get("roster", {}).get("actual_validation_N"),
                        training_budget_note="Independent declared training protocols; equal test people does not imply equal training budgets",
                        count_minus_benchmark_loss=uncertainty["estimate"],
                        paired_delta_low=uncertainty["low"], paired_delta_high=uncertainty["high"],
                        bootstrap_unit="count_declared_group", n_bootstrap_groups=uncertainty["groups"],
                        benchmark_receipt_file_record=receipt_hash, count_evaluation_receipt_file_record=file_record(complete),
                        validation_status=model["validation_status"])
                    for key, value in count_metrics.items():
                        if key != "n_evaluable":
                            record["count_"+key] = value
                    for key, value in benchmark_metrics.items():
                        if key != "n_evaluable":
                            record["benchmark_"+key] = value
                    rows.append(record)
    return (pd.DataFrame(rows) if rows else pd.DataFrame(columns=["trait", "benchmark_family", "benchmark_model", "population", "status", "n_same_people"]),
            pd.DataFrame(statuses))
