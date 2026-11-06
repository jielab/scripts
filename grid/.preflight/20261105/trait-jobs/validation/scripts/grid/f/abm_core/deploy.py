"""Label-free deployment against a frozen bank; inputs must carry producer bindings."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import joblib
import numpy as np
import pandas as pd

from .contracts import verify_array_binding, write_feature_manifest
from .grid_inputs import digest
from .model import ABMAssessor, file_sha256, write_json
from .predictors import predict
from .support import count_twins, _backend


def assess_dataset(run, features, geometry_array, geometry_binding, *, device="auto",
                   raw_blocks=None, raw_binding=None):
    """Arrays and side features come from a declared target-cohort producer, not Y."""
    run = Path(run)
    if {"y", "outcome", "truth", "loss_base", "loss_candidate", "realized_gain"} & set(features):
        raise ValueError("Deployment cannot accept query labels or realized errors")
    local = run / "_LOCAL_ONLY"
    model = ABMAssessor.load(local / "assessor.json")
    if file_sha256(local / "predictor_bundle.joblib") != model.manifest["predictor_bundle_sha256"]:
        raise ValueError("Frozen predictor bundle changed")
    geometry = json.loads((run / "geometry_manifest.json").read_text())
    if file_sha256(run / "geometry_manifest.json") != model.manifest["geometry_sha256"]:
        raise ValueError("Frozen geometry changed")
    bound = verify_array_binding(geometry_array, features, geometry_binding)
    if bound["geometry_sha256"] != geometry["geometry_sha256"]:
        raise ValueError("Query was transformed using a different build geometry")
    pointer = json.loads((local / "support_cache.json").read_text())
    cache = Path(pointer["path"])
    if file_sha256(cache / "COMPLETE.json") != pointer["receipt_sha256"]:
        raise ValueError("Frozen donor cache changed")
    receipt = json.loads((cache / "COMPLETE.json").read_text())
    for name in ("bank.parquet", "donor_geometry.npy"):
        if file_sha256(cache / name) != receipt["files"][name]:
            raise ValueError("Frozen donor bank changed")
    bank = pd.read_parquet(cache / "bank.parquet")
    donors = np.load(cache / "donor_geometry.npy", mmap_mode="r")
    qz = np.load(geometry_array, mmap_mode="r")
    manifest = json.loads((run / "RUN_MANIFEST.json").read_text())
    cfg = manifest["contract"]["config"]
    opt = cfg["support"]
    used_device, _ = _backend(device)
    # Nonfinite geometry is an explicit QC failure. Restore those rows with NA.
    valid = np.isfinite(qz).all(1)
    result = count_twins(donors, qz[valid], donor_ids=bank.id, donor_groups=bank.group_id,
        query_ids=features.loc[valid, "id"], query_groups=features.loc[valid, "group_id"],
        radius=geometry["radius"], device=used_device, k_used=opt["k_used"],
        min_count=opt["min_count"], min_group_ess=opt["min_group_ess"],
        query_batch=cfg["execution"]["query_batch"], donor_batch=cfg["execution"]["donor_batch"],
        # The target producer's group names are bound, but that alone does not
        # establish kinship verification in a new cohort.
        groups_verified=False)
    query = features.merge(result.table, on="id", how="left", validate="one_to_one")
    bad = query.supported.isna()
    query.loc[bad, "supported"] = False
    query.loc[bad, "groups_verified"] = False
    query.loc[bad, "support_reason"] = "genotype_qc_nonfinite_blocks"
    bundle = joblib.load(local / "predictor_bundle.joblib")
    blocks = None
    if bundle["candidate_kind"] == "blocknet":
        if raw_blocks is None or raw_binding is None:
            raise ValueError("Blocknet requires bound raw CSx block inputs as well as transformed geometry")
        inp = json.loads((run / "input_manifest.json").read_text())["blocks"]
        bound = verify_array_binding(raw_blocks, features, raw_binding,
                                     tokens=bundle["tokens"], channels=inp["channels"])
        expected = digest(dict(layout=inp["layout"]["signature"], weights=inp["signature_payload"]["weights"]))
        if bound["geometry_sha256"] != expected:
            raise ValueError("Raw blocks have a different LD layout or CSx posterior")
        blocks = np.load(raw_blocks, mmap_mode="r")
    if cfg["support"]["representation"] == "csx_scores":
        from .representation import BlockGeometry
        transform = BlockGeometry(**geometry["state"])
        score_columns = [bundle["features"][i] for i in bundle["csx_columns"]]
        expected = transform.transform(features[score_columns].to_numpy(float))
        if not np.allclose(qz, expected, rtol=1e-12, atol=1e-12, equal_nan=True):
            raise ValueError("Score-only geometry does not match the aligned four score columns")
    p0, p1 = predict(bundle, query, blocks=blocks, block_rows=np.arange(len(query)), device=used_device)
    query["p_base"], query["p_candidate"] = p0, p1
    return model.assess(query, manifest=model.manifest), model


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", required=True)
    p.add_argument("--features", required=True, help="Parquet with id, group_id and all declared side columns; ancestry is optional reporting metadata; no Y")
    p.add_argument("--geometry-array", required=True)
    p.add_argument("--geometry-binding", required=True)
    p.add_argument("--raw-blocks")
    p.add_argument("--raw-binding")
    p.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    p.add_argument("--output", required=True, help="Local output directory")
    args = p.parse_args(argv)
    out, model = assess_dataset(args.run, pd.read_parquet(args.features), args.geometry_array,
                                args.geometry_binding, device=args.device,
                                raw_blocks=args.raw_blocks, raw_binding=args.raw_binding)
    root = Path(args.output)
    local = root / "_LOCAL_ONLY"
    local.mkdir(parents=True, exist_ok=True)
    out.to_parquet(local / "predictions.parquet", index=False)
    write_feature_manifest(local / "predictions.parquet", manifest=model.manifest)
    write_json(root / "dataset_assessment.json", model.dataset_assessment(out))
    write_json(local / "FROZEN.json", dict(
        predictions_sha256=file_sha256(local / "predictions.parquet"),
        model_sha256=file_sha256(Path(args.run) / "_LOCAL_ONLY/assessor.json"), query_outcomes_used=False))
    print(json.dumps({"output": str(root), "n_query": len(out), "n_assessable": int(out.assessment_available.sum())}))


if __name__ == "__main__":
    main()
