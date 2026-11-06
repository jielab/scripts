"""Optional, exact-person comparison against current independent benchmarks.

Only the new completion contract is accepted. Missing/incompatible methods are
reported, never substituted by historical aggregate tables or a smaller common
intersection. This module cannot fit, calibrate or select a benchmark model.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .inputs import file_hash, read_table, ids_hash
from .report import metrics, bootstrap_mean


def benchmark_receipts(root, trait):
    root = Path(root)
    return {"csx": root/"csx"/trait/"targetfit/COMPLETE.json",
            "disco": root/"disco"/trait/"COMPLETE.json",
            "prsformer": root/"prsformer"/trait/"benchmark/COMPLETE.json"}


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
    if not path.is_file() or file_hash(path) != expected:
        raise ValueError("Bound snapshot is missing or changed: "+key)
    return path


def load_benchmark(path, *, method, trait, task, endpoint, count_frame, count_csx_score_sha256=None):
    path = Path(path)
    before = file_hash(path)
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
        if not source.is_file() or file_hash(source) != expected:
            raise ValueError("Benchmark producer-bound artifact missing or changed: "+source.name)
    if method == "csx" and count_csx_score_sha256 is not None:
        score_path = receipt.get("scores")
        if not score_path or receipt["files"].get(str(Path(score_path).resolve())) != count_csx_score_sha256:
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
    if file_hash(path) != before:
        raise ValueError("Benchmark completion receipt changed while reading")
    return values, receipt, before


def compare_current(root, traits, *, repetitions=200, seed=20261102):
    rows, statuses = [], []
    for trait in traits:
        folder = Path(root)/"count_abm"/trait
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
                    count_csx_score_sha256=model["metadata"].get("whole_csx_score_sha256"))
            except (ValueError, FileNotFoundError, KeyError, TypeError, OSError) as exc:
                statuses.append(dict(status, status="not_comparable", reason=str(exc)))
                continue
            statuses.append(dict(status, status="comparable_exact_test_roster_and_labels",
                                 reason="No intersection filtering; missing-label positions and optional group mapping also match",
                                 benchmark_receipt_sha256=receipt_hash, models=",".join(values)))
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
                        csx_source_identity=("matched_count_producer_score" if method == "csx" and model["metadata"].get("whole_csx_score_sha256")
                                             else "not_applicable_or_external_export"),
                        benchmark_validation_n=receipt.get("roster", {}).get("actual_validation_N"),
                        training_budget_note="Independent declared training protocols; equal test people does not imply equal training budgets",
                        count_minus_benchmark_loss=uncertainty["estimate"],
                        paired_delta_low=uncertainty["low"], paired_delta_high=uncertainty["high"],
                        bootstrap_unit="count_declared_group", n_bootstrap_groups=uncertainty["groups"],
                        benchmark_receipt_sha256=receipt_hash, count_evaluation_receipt_sha256=file_hash(complete),
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
