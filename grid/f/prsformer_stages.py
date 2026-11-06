"""Bind CPU publication to the exact outputs of its preceding model worker."""
from __future__ import annotations

import json
import os
from pathlib import Path

from resume import atomic_json, sha


def receipt_path(method, trait):
    directory = os.environ.get("GRID_PRSFORMER_STAGE_DIR")
    if not directory:
        raise ValueError("Split PRSformer workers require their planner's stage directory")
    return Path(directory) / f"{method}.{trait}.json"


def freeze_model(root, method, trait, run, cache):
    from prsformer import validate_evaluation
    root, run, cache = Path(root).resolve(), Path(run).resolve(), Path(cache).resolve()
    validate_evaluation(run, cache, analysis_root=root)
    paths = [run / name for name in
             ("evaluation.json", "model.pt", "metrics.tsv", "test_predictions.tsv.gz")]
    paths += [cache / name for name in ("prepare.json", "data.tsv.gz", "variants.tsv.gz")]
    if method == "prsformer2":
        paths += [run / "roles.tsv.gz", run.parent / "RUN_MANIFEST.json"]
    active = root / "local/kinship/ACTIVE.json"
    record = dict(schema="grid.prsformer.model_handoff.1", method=method, trait=trait,
                  root=str(root), run=str(run), cache=str(cache),
                  active_sha256=sha(active) if active.is_file() else None,
                  files={str(p): sha(p) for p in paths})
    atomic_json(receipt_path(method, trait), record)
    return run


def frozen_model(root, method, trait, cache):
    record = json.loads(receipt_path(method, trait).read_text())
    root, cache = Path(root).resolve(), Path(cache).resolve()
    if (record.get("schema") != "grid.prsformer.model_handoff.1"
            or record.get("method") != method or record.get("trait") != trait
            or record.get("root") != str(root) or record.get("cache") != str(cache)):
        raise ValueError("PRSformer stage identity changed")
    run = Path(record["run"])
    run.relative_to(root / method / trait)
    required = {str(run / n) for n in
                ("evaluation.json", "model.pt", "metrics.tsv", "test_predictions.tsv.gz")}
    required.update(str(cache / n) for n in ("prepare.json", "data.tsv.gz", "variants.tsv.gz"))
    if method == "prsformer2":
        required.update((str(run / "roles.tsv.gz"), str(run.parent / "RUN_MANIFEST.json")))
    if not required <= set(record.get("files", {})):
        raise ValueError("Incomplete PRSformer model handoff")
    for name, expected in record["files"].items():
        if not Path(name).is_file() or sha(name) != expected:
            raise ValueError("Frozen PRSformer output/input changed: " + name)
    active = root / "local/kinship/ACTIVE.json"
    if (sha(active) if active.is_file() else None) != record.get("active_sha256"):
        raise ValueError("Kinship selection changed between PRSformer stages")
    return run
