"""Read trusted local comparator exports; never fit or modify their source run."""
from pathlib import Path
import argparse
import json
import joblib
import numpy as np
import pandas as pd

from abm_core.model import file_sha256, write_json
from abm_core.evaluation import prediction_metrics
from twin_core import id_hash, validate_ids, state_hash


def ids(values, label):
    a = np.asarray(values)
    if a.ndim != 1 or pd.isna(a).any():
        raise ValueError(label + ": invalid IDs")
    a = a.astype(str)
    validate_ids(a)
    return a


def import_comparators(analysis_root, run, output):
    base = Path(analysis_root).resolve()
    source = Path(run)
    source = source.resolve() if source.is_absolute() else (base / source).resolve()
    if not source.is_relative_to(base):
        raise ValueError("--run must be inside the trusted local --analysis-root")
    out = Path(output).resolve()
    if out == source or out.is_relative_to(source):
        raise ValueError("Historical output must not modify the source run")
    if out.exists():
        raise FileExistsError(out)
    required = ["COMPLETE.json", "RESUME.json", "protocol.json", "SOURCE_FROZEN.json",
                "_LOCAL_ONLY/predictions.joblib", "_LOCAL_ONLY/comparator_predictions.joblib",
                "_LOCAL_ONLY/test_truth.tsv.gz", "_LOCAL_ONLY/test_X.npz",
                "_LOCAL_ONLY/development_roles.tsv.gz", "_LOCAL_ONLY/development.npz"]
    missing = [str(source / name) for name in required if not (source / name).is_file()]
    rows, methods, evidence = [], [], {}
    if missing:
        methods = [dict(old_name="Disco_reference", run=str(source),
            actual_source=str(source / "_LOCAL_ONLY/comparator_predictions.joblib"),
            method="csx-disco-legacy/raw-interpolation", status="unavailable_missing_artifacts",
            missing=json.dumps(missing))]
    else:
        complete = json.loads((source / "COMPLETE.json").read_text())
        protocol = json.loads((source / "protocol.json").read_text())
        resume = json.loads((source / "RESUME.json").read_text())
        hashes = {name: sha for files in resume["stages"].values() for name, sha in files.items()}
        if complete.get("status") != "completed":
            raise ValueError("Historical run is incomplete")
        if protocol.get("abm_version") != "GRID-ABM-FINAL-2026-10-17":
            raise ValueError("Unknown historical writer; cannot infer outcome/loss convention")
        for name in required:
            if name in ("RESUME.json", "COMPLETE.json"):
                continue
            if name not in hashes or file_sha256(source / name) != hashes[name]:
                raise ValueError("Missing/changed historical receipt binding: " + name)
        evidence = {name: file_sha256(source / name) for name in required}
        p = joblib.load(source / "_LOCAL_ONLY/predictions.joblib")
        ext = joblib.load(source / "_LOCAL_ONLY/comparator_predictions.joblib")
        if set(ext) != {"eid", "predictions"} or not isinstance(ext["predictions"], dict):
            raise ValueError("Comparator export must use the abm_resume writer schema")
        roster, external_ids = ids(p["eid"], "prediction"), ids(ext["eid"], "comparator")
        if p["frozen_hash"] != complete["model_hash"] or id_hash(roster) != complete["test_roster_hash"]:
            raise ValueError("Historical model/test roster binding mismatch")
        if set(roster) != set(external_ids):
            raise ValueError("Comparator roster mismatch")
        truth = pd.read_csv(source / "_LOCAL_ONLY/test_truth.tsv.gz", sep="\t", dtype={"eid": str})
        truth["eid"] = ids(truth.eid, "truth")
        if set(roster) != set(truth.eid):
            raise ValueError("Historical truth roster mismatch")
        truth = truth.set_index("eid").loc[roster]
        with np.load(source / "_LOCAL_ONLY/test_X.npz", allow_pickle=False) as query:
            qids = ids(query["ids"], "test_X")
            if set(qids) != set(roster) or len(query["X"]) != len(qids):
                raise ValueError("Historical feature roster mismatch")
            features = query["features"].astype(str).tolist()
            test_groups = set(query["groups"].astype(str))
        roles = pd.read_csv(source / "_LOCAL_ONLY/development_roles.tsv.gz", sep="\t",
                            dtype={"eid": str, "group": str})
        ids(roles.eid, "development roles")
        if set(roles.eid) & set(roster) or set(roles.group) & test_groups or roles.groupby("group").role.nunique().max() != 1:
            raise ValueError("Historical development/test group leakage")
        with np.load(source / "_LOCAL_ONLY/development.npz", allow_pickle=False) as dev:
            dids = ids(dev["ids"], "development")
            aligned = roles.set_index("eid").loc[dids]
            if set(dids) != set(roles.eid) or not np.array_equal(dev["roles"].astype(str), aligned.role):
                raise ValueError("Historical calibration roles mismatch")
            if not np.array_equal(dev["groups"].astype(str), aligned.group) or dev["features"].astype(str).tolist() != features:
                raise ValueError("Historical feature/group contract mismatch")
        task = protocol["config"]["trait_type"]
        trait = protocol["config"]["trait"]
        endpoint = protocol["endpoint_definition"]
        if trait == "t2dm" and endpoint != "baseline_t2dm_from_Yr2e_Yt2e":
            raise ValueError("Unsupported historical T2DM endpoint")
        if "CSx_reference" in p.get("predictions", {}):
            values = np.asarray(p["predictions"]["CSx_reference"], float)
            if values.shape != (len(roster),):
                raise ValueError("Historical baseline prediction length mismatch")
            model_path = source / "fit/_LOCAL_ONLY/model.joblib"
            record = dict(old_name="CSx_reference", method="CSx_reference", run=str(source),
                actual_source=str(source / "_LOCAL_ONLY/predictions.joblib"),
                predictor="saved_linear_model", representation="original_CSx_and_covariates",
                calibration="unverified", status="unavailable_frozen_model",
                missing=str(model_path), endpoint=endpoint,
                assessment_loss="brier" if task == "binary" else "mse", units="original_run_native_units_unverified")
            if model_path.is_file() and hashes.get("fit/_LOCAL_ONLY/model.joblib") == file_sha256(model_path):
                model = joblib.load(model_path)
                if state_hash(model) != complete["model_hash"]:
                    raise ValueError("Historical baseline fitted-state mismatch")
                if list(model.features) != features:
                    raise ValueError("Historical baseline feature names/order mismatch")
                original = np.asarray(model.prep.keep)[np.asarray(model.columns["CSx_reference"])]
                columns = [features[int(i)] for i in original]
                calibration = model.calibration["CSx_reference"]
                record.update(feature_set=json.dumps(columns), n_features=len(columns),
                    calibration="saved_calibration_fit; " + type(calibration).__name__,
                    status="available", missing=None)
                rows.append(dict(record, trait=trait,
                                 **prediction_metrics(truth.y.to_numpy(float), values, task)))
                evidence[str(model_path.relative_to(source))] = file_sha256(model_path)
            methods.append(record)
        for name, values in ext["predictions"].items():
            values = np.asarray(values, float)
            if values.ndim != 1 or len(values) != len(external_ids):
                raise ValueError("Comparator array length mismatch: " + name)
            calibrated = False
            calibration = "unverified_missing_comparator_calibration"
            columns = features
            method = name
            if name == "Disco_reference":
                method = "csx-disco-legacy/raw-interpolation"
                path = source / "_LOCAL_ONLY/disco_calibration.joblib"
                if path.is_file() and hashes.get("_LOCAL_ONLY/disco_calibration.joblib") == file_sha256(path):
                    bundle = joblib.load(path)
                    if not {"columns", "imputer", "scale", "model", "calibration"} <= set(bundle):
                        raise ValueError("Unexpected saved Disco calibration schema")
                    columns = list(bundle["columns"])
                    calibration = "saved_build_tune_model_calibration_fit; " + type(bundle["calibration"]).__name__
                    calibrated = True
                    evidence[str(path.relative_to(source))] = file_sha256(path)
            # External methods need their own recorded calibration evidence; never guess from display names.
            record = dict(old_name=name, method=method, run=str(source),
                actual_source=str(source / "_LOCAL_ONLY/comparator_predictions.joblib"),
                predictor="locally_calibrated_raw_interpolation" if name == "Disco_reference" else name,
                representation="raw_ancestry_scores_reference_centers" if name == "Disco_reference" else "source_declared_unverified",
                feature_set=json.dumps(columns), n_features=len(columns), calibration=calibration,
                endpoint=endpoint, assessment_loss="brier" if task == "binary" else "mse",
                units="original_run_native_units_unverified",
                status="available" if calibrated else "unavailable_calibration_provenance")
            methods.append(record)
            if calibrated:
                pred = pd.Series(values, index=external_ids).loc[roster].to_numpy()
                rows.append(dict(record, trait=trait, **prediction_metrics(truth.y.to_numpy(float), pred, task)))
    out.mkdir(parents=True)
    pd.DataFrame(methods).to_csv(out / "legacy_method_map.tsv", sep="\t", index=False)
    pd.DataFrame(rows, columns=None if rows else ["method", "status"]).to_csv(out / "historical_performance.tsv", sep="\t", index=False)
    write_json(out / "IMPORT_MANIFEST.json", dict(source=str(source), source_hashes=evidence,
        missing=missing, source_modified=False, training=False, comparison_scope="historical_only"))
    write_json(out / "table_manifest.json", [dict(name=p.stem, path=p.name, sha256=file_sha256(p),
        module="historical", run_id=source.name, trait="all") for p in out.glob("*.tsv")])
    write_json(out / "COMPLETE.json", dict(status="completed", available_rows=len(rows),
        manifest_sha256=file_sha256(out / "IMPORT_MANIFEST.json")))
    print(json.dumps(dict(output=str(out), available_rows=len(rows), missing=missing)))
    return pd.DataFrame(rows), pd.DataFrame(methods)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis-root", required=True)
    p.add_argument("--run", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args(argv)
    import_comparators(a.analysis_root, a.run, a.output)


if __name__ == "__main__":
    main()
