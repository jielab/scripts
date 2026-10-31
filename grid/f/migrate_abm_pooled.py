"""Verify and migrate the audited pre-pooled ABM caches; never refit a model.

Old completed runs remain immutable history. Only redundant predictor caches,
unused interpolation inputs, and obsolete CURRENT pointers are removed on --apply.
Large immutable support arrays use hard links so historical paths keep working.
"""
from __future__ import annotations
import argparse
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

import joblib
import pandas as pd

from abm_core import pipeline
from abm_core.derived_report import verify_run
from abm_core.grid_inputs import digest
from abm_core.model import file_sha256, write_json
from abm_core.storage import support_lock, validate_support, busy_reason

COMPARATOR = "grid/f/abm_core/comparators.py"
OLD_COMPARATOR = "4cab75c6d27654592b1a93ee0884089771b162a9dc431898bc3d7ff8f89d7d2e"
OLD_PIPELINE = "ace06eb21411743903e77e85a07ea13eefcc510cd7f4c4a05efc8ed54be99370"
OLD_REPRESENTATION = "5563603a892a859335f1eee14d10616fb8b8f23dbb6505c3f7fdccc9bf53d873"
NEW_REPRESENTATION = "10959859860844877977085f3c32095dd3bd3e4b15b65bd9c130323611a3715e"
SUPPORT_FUNCTIONS = {
    "prepare_support": "c8f8f5320363a8057375e8ebdd23543771e71261b996f7ffe208b510b5289d90",
    "_prepare_support": "2a6cba56b71bc4823c2f6add59f7af81d4b7e2873df799060b981ef33ee24184",
}
RETAINED = {"Clinical", "csx-moe"}


def project_identity(identity, sources):
    """Reject changed mathematical dependencies before dropping the retired fit."""
    if identity.get("schema_version") or identity["sources"].get(COMPARATOR) != OLD_COMPARATOR:
        raise ValueError("Not the audited legacy predictor cache")
    out = copy.deepcopy(identity)
    out.pop("ancestry")
    out["sources"].pop("grid/f/" + "2." + "disco.R", None)
    out["sources"].pop("grid/f/3.disco.R", None)
    for key, expected in out["sources"].items():
        if key != COMPARATOR and sources.get(key) != expected:
            raise ValueError("Predictor dependency changed: " + key)
    out["sources"][COMPARATOR] = sources[COMPARATOR]
    out["schema_version"] = "abm.predictors.pooled.1"
    return out


def verify_retained_code(frozen, sources):
    path = frozen / COMPARATOR
    if file_sha256(path) != OLD_COMPARATOR:
        raise ValueError("Legacy comparator source changed")
    def body(path):
        fn = next(n for n in ast.parse(path.read_text()).body
                  if isinstance(n, ast.FunctionDef) and n.name == "fit_comparators")
        kept = []
        for node in fn.body:
            if isinstance(node, (ast.Try, ast.Return)):
                break
            kept.append(ast.dump(node))
        return kept
    if body(path) != body(pipeline.ROOT / "f/abm_core/comparators.py"):
        raise ValueError("Retained Clinical/MoE fitting code changed")
    # Audited differences are mmap closing, scratch release and cache locking;
    # distances, scaling, radius estimation and donor eligibility are unchanged.
    for name, expected in (("pipeline.py", OLD_PIPELINE), ("representation.py", OLD_REPRESENTATION)):
        if file_sha256(frozen / "grid/f/abm_core" / name) != expected:
            raise ValueError("Unknown legacy support implementation: " + name)
    if sources["grid/f/abm_core/representation.py"] != NEW_REPRESENTATION:
        raise ValueError("Current geometry implementation needs a new audit")
    source = (pipeline.ROOT / "f/abm_core/pipeline.py").read_text()
    tree = ast.parse(source)
    for name, expected in SUPPORT_FUNCTIONS.items():
        fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
        if hashlib.sha256(ast.get_source_segment(source, fn).encode()).hexdigest() != expected:
            raise ValueError("Current support implementation needs a new audit: " + name)


def verify_files(folder, files):
    for name, expected in files.items():
        path = (folder / name).resolve()
        if not path.is_relative_to(folder.resolve()) or file_sha256(path) != expected:
            raise ValueError("Cache artifact changed: " + str(path))


def write_predictor_cache(source, destination, receipt, identity):
    """Keep primary predictions/model bytes whenever no retired model was fitted."""
    if destination.exists():
        saved = json.loads((destination / "COMPLETE.json").read_text())
        if saved["identity"] != identity:
            raise ValueError("Existing migrated identity differs")
        verify_files(destination, saved["files"])
        return
    destination.parent.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=".pooled-", dir=destination.parent.parent))
    stage = temp / "_LOCAL_ONLY"
    stage.mkdir()
    try:
        for name in receipt["files"]:
            out = stage / name
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / name, out)
        bundle_path = stage / "_LOCAL_ONLY/predictor_bundle.joblib"
        bundle = joblib.load(bundle_path)
        if set(bundle.get("comparators", {})) - RETAINED:
            bundle["comparators"] = {k:v for k,v in bundle["comparators"].items() if k in RETAINED}
            joblib.dump(bundle, bundle_path, compress=3)
        prediction_path = stage / "_LOCAL_ONLY/predictor_predictions.parquet"
        predictions = pd.read_parquet(prediction_path)
        kept = {"id", "group_id", "role", "p_base", "p_candidate"} | RETAINED
        if set(predictions) - kept:
            predictions[[c for c in predictions if c in kept]].to_parquet(prediction_path, index=False)
        status = pd.read_csv(stage / "comparator_status.tsv", sep="\t")
        status.loc[status.method.isin(RETAINED)].to_csv(stage / "comparator_status.tsv", sep="\t", index=False)
        write_json(stage / "predictor_identity.json", identity)
        saved = dict(identity=identity, executed_device=receipt["executed_device"],
                     files={n:file_sha256(stage/n) for n in receipt["files"]},
                     migration=dict(source=str(source), source_receipt_sha256=file_sha256(source/"COMPLETE.json"),
                                    retained_fits="identical pooled primary and Clinical/MoE fits"))
        write_json(stage / "COMPLETE.json", saved)
        verify_files(stage, saved["files"])
        os.replace(temp, destination.parent)
    finally:
        if temp.exists():
            shutil.rmtree(temp)


def migrate(root, trait, sources, apply=False):
    root = Path(root).resolve()
    parent = root / "abm/FINAL" / trait
    pointer = parent / "blocknet.CURRENT.json"
    if not pointer.is_file():
        return dict(trait=trait, status="no_current_legacy_run")
    original_pointer = pointer.read_bytes()
    run = (parent / json.loads(original_pointer)["run"]).resolve()
    if run.parent != parent or run.is_symlink():
        raise ValueError("Unsafe current run")
    manifest, _ = verify_run(run)
    contract = manifest["contract"]
    if contract["config"]["evaluation"]["matching"] == "frozen_risk_bin_disjoint_group_exclusion":
        return dict(trait=trait, status="already_pooled")
    verify_retained_code(run / "_LOCAL_ONLY/SOURCE_FROZEN", sources)
    for name in ("support.py", "grid_inputs.py"):
        key = "grid/f/abm_core/" + name
        if contract["sources"][key] != sources[key]:
            raise ValueError("Support dependency changed: " + name)
    old_identity = json.loads((run / "predictor_identity.json").read_text())
    new_identity = project_identity(old_identity, sources)
    fit_root = root / "cache/abm/predictors"
    old_fit = fit_root / digest(old_identity) / "_LOCAL_ONLY"
    new_fit = fit_root / digest(new_identity) / "_LOCAL_ONLY"
    receipt = json.loads((old_fit / "COMPLETE.json").read_text())
    if receipt["identity"] != old_identity:
        raise ValueError("Legacy predictor identity mismatch")
    verify_files(old_fit, receipt["files"])
    support_pointer = json.loads((run / "_LOCAL_ONLY/support_cache.json").read_text())
    old_support = Path(support_pointer["path"]).resolve()
    support_root = root / "cache/abm/support" / trait
    if old_support.parent.parent != support_root or file_sha256(old_support/"COMPLETE.json") != support_pointer["receipt_sha256"]:
        raise ValueError("Legacy support path/receipt mismatch")
    frame = pd.read_parquet(run / "_LOCAL_ONLY/features.parquet", columns=["id","group_id","role"])
    inputs = json.loads((run / "input_manifest.json").read_text())
    args = (trait, frame, contract["config"], contract["representation"], inputs["scores"], inputs["blocks"])
    if digest(pipeline.support_identity(*args, contract["sources"])) != old_support.parent.name:
        raise ValueError("Legacy support identity mismatch")
    new_support = support_root / digest(pipeline.support_identity(*args, sources)) / "_LOCAL_ONLY"
    result = dict(trait=trait, run=str(run), old_fit=str(old_fit), new_fit=str(new_fit),
                  old_support=str(old_support), new_support=str(new_support),
                  status="planned", old_results="preserved_as_history", removed=[])
    if not apply:
        return result
    with support_lock(old_support):
        reason = busy_reason([p for p in old_support.iterdir() if p.name != ".support.lock"]
                             + list(old_fit.rglob("*")))
        if reason:
            raise BlockingIOError(reason)
        print("[pooled migration] Verify persistent support: " + trait, flush=True)
        validate_support(old_support)
        write_predictor_cache(old_fit, new_fit, receipt, new_identity)
        if new_support.exists():
            if file_sha256(new_support/"COMPLETE.json") != support_pointer["receipt_sha256"]:
                raise ValueError("Existing support destination differs")
            validate_support(new_support)
        else:
            temp = Path(tempfile.mkdtemp(prefix=".pooled-", dir=support_root))
            try:
                stage = temp / "_LOCAL_ONLY"
                stage.mkdir()
                for file in old_support.iterdir():
                    if file.is_symlink() or not file.is_file():
                        raise ValueError("Unexpected support entry: " + str(file))
                    os.link(file, stage / file.name)
                os.replace(temp, new_support.parent)
            finally:
                if temp.exists():
                    shutil.rmtree(temp)
        if pointer.read_bytes() != original_pointer:
            raise ValueError("Current run changed during migration")
        # Historical run copies retain their original predictor and status files.
        shutil.rmtree(old_fit.parent)
        result["removed"].append(str(old_fit.parent))
        # These unsigned interpolation inputs are disposable; retain the R error log.
        for file in (run / "_LOCAL_ONLY/disco").glob("*.tsv"):
            if file.name in {"AFR.tsv","EAS.tsv","EUR.tsv","SAS.tsv","pca.tsv","centers.tsv"}:
                file.unlink()
                result["removed"].append(str(file))
        pointer.unlink()
        result["removed"].append(str(pointer))
    result["status"] = "migrated"
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    p.add_argument("--traits", default="height,ldl,t2dm")
    p.add_argument("--apply", action="store_true")
    p.add_argument("--report", required=True)
    args = p.parse_args()
    traits = args.traits.split(",")
    if len(set(traits)) != len(traits) or not set(traits) <= {"height","ldl","t2dm"}:
        p.error("Use unique traits from height,ldl,t2dm")
    sources = pipeline.source_manifest()
    report = Path(args.report)
    report.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for trait in traits:
        row = migrate(args.analysis_root, trait, sources, args.apply)
        records.append(row)
        write_json(report, records)
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
