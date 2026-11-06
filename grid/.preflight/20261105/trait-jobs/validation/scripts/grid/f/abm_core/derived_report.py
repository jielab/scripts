"""Verified, read-only reuse of frozen scientific artifacts for derived reports."""
from pathlib import Path
import json
import shutil
import joblib
import pandas as pd

from .model import ABMAssessor, file_sha256, write_json
from .grid_inputs import digest
from . import reporting

# Audited unversioned writer: its only targets were MSE/Brier. Never infer unknown writers.
UNVERSIONED_REPORTING_SHA = "3b8502df59f14b431ed284215414e48c1ee1e80580a1fb04c0f554aea0a6ecea"


def verify_run(run):
    run = Path(run).resolve()
    complete = json.loads((run / "COMPLETE.json").read_text())
    if complete.get("status") != "completed":
        raise ValueError("Parent run is incomplete")
    if file_sha256(run / "RUN_MANIFEST.json") != complete["manifest_sha256"]:
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
    from .pipeline import source_manifest
    parent = Path(parent).resolve()
    manifest, complete = verify_run(parent)
    sources = source_manifest()
    config = manifest["contract"]["config"]
    local = parent / "_LOCAL_ONLY"
    model = ABMAssessor.load(local / "assessor.json")
    trait = manifest["contract"]["trait"]
    if config["primary_loss"][trait] != model.heads.loss_name:
        raise ValueError("Configured assessment loss differs from frozen head; refit affected heads/rules first")
    identity = dict(parent=str(parent), parent_manifest_sha256=file_sha256(parent / "RUN_MANIFEST.json"),
        parent_phases_sha256=file_sha256(parent / "PHASES.json"),
        parent_completion_sha256=file_sha256(parent / "COMPLETE.json"), reporting_sources=sources)
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
        prior = manifest["contract"]["sources"].get("grid/f/abm_core/reporting.py")
        if prior != UNVERSIONED_REPORTING_SHA or model.heads.loss_name != (
                "brier" if model.heads.task == "binary" else "mse"):
            raise ValueError("Unversioned diagnostic objective has no recognized source evidence")
        diagnostics = dict(diagnostics, schema_version="abm.diagnostics.1",
                           loss_name=model.heads.loss_name, task=model.heads.task)
        identity["diagnostic_migration"] = dict(source_sha256=prior, convention="audited_MSE_Brier_writer")
    if reporting.diagnostic_loss(diagnostics) != model.heads.loss_name:
        raise ValueError("Frozen diagnostic/head target mismatch")
    truth = pd.read_parquet(local / "truth.parquet")
    if truth.id.astype(str).duplicated().any():
        raise ValueError("Duplicate truth IDs")
    root.mkdir(parents=True)
    for label in ("test", "audit", "policy_calibration"):
        folder = local / label
        frozen = json.loads((folder / "FROZEN.json").read_text())
        if (frozen["model_sha256"] != file_sha256(local / "assessor.json") or
            frozen["predictions_sha256"] != file_sha256(folder / "predictions.parquet") or
            frozen.get("query_outcomes_used") is not False):
            raise ValueError("Frozen prediction identity failed")
        if frozen.get("truth_sha256", file_sha256(local / "truth.parquet")) != file_sha256(local / "truth.parquet"):
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
        generated_private_sha256={str(p.relative_to(root)): file_sha256(p)
                                  for p in (root / "_LOCAL_ONLY").rglob("*") if p.is_file()},
        original_SOURCE_FROZEN=str(parent / "SOURCE_FROZEN.json")))
    tables = [dict(name=p.stem, path=p.name, sha256=file_sha256(p), trait=trait,
                   module="abm", run_id=root.name) for p in sorted(root.glob("*.tsv"))]
    write_json(root / "table_manifest.json", tables)
    write_json(root / "COMPLETE.json", dict(status="completed", smoke=complete.get("smoke", False),
        manifest_sha256=file_sha256(root / "REPORT_MANIFEST.json"),
        table_manifest_sha256=file_sha256(root / "table_manifest.json")))
    print(json.dumps(dict(output=str(root), parent=str(parent), training=False)))
    return root
