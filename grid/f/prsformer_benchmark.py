"""Publish a verified standalone benchmark receipt after PRSformer evaluation.

This module only reads an already validated fitted model/evaluation. It does not
calibrate, tune, retrain or change any original PRSformer performance measure.
"""
from __future__ import annotations
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from shared_cohort import load_cohort, sha256
from resume import atomic_json, run_lock
from target_benchmark import metrics


def publish(root, trait, method="prsformer", cache=None):
    from prsformer import validate_evaluation
    from kinship import analysis_guard
    root = Path(root)
    if method == "prsformer":
        run = root / "prsformer" / trait / "run"
        label = "PRSformer independent-trait adaptation"
    elif method == "prsformer2":
        base = root / "prsformer2" / trait
        current = json.loads((base / "CURRENT.json").read_text())
        run = base / current["run"] / "_LOCAL_ONLY"
        label = "PRSformer2 exploratory extension"
    else:
        raise ValueError("Unknown PRSformer benchmark")
    cache = Path(cache) if cache is not None else root / "prsformer/cache"
    out = root / method / trait / "benchmark"
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
            files={str(path.resolve()):sha256(path) for path in files}, implementation_sha256=sha256(__file__))
        atomic_json(out / "COMPLETE.json", receipt)
        print(f"[{method}] verified independent benchmark: {out / 'COMPLETE.json'}", flush=True)
        return receipt
