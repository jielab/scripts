"""Reference CLI for integrating frozen GRID predictions with ABM FINAL.

Individual-level TSVs and match pairs belong under _LOCAL_ONLY. The public
result JSON omits IDs. Upstream genotype/PRS/net fitting is performed by GRID.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .evaluation import evaluate_predictions
from .model import ABMAssessor, file_sha256, write_json, coverage_tag
from .support import count_twins
from .contracts import verify_array_binding, verify_feature_binding


def read_table(path):
    # Individual identifiers must not lose leading zeros during CSV inference.
    return pd.read_csv(path, sep="\t", dtype={"id": str, "group_id": str})


def _load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _evaluate(frame, model, objective, coverage, n_boot, seed):
    col = f"selected_{objective}_{coverage_tag(coverage)}"
    result = evaluate_predictions(
        frame, task=model.heads.task, loss_name=model.heads.loss_name,
        selector=f"twin_count_{objective}", candidate_id=model.manifest["candidate_id"],
        coverage=coverage, selected_col=col, control_eligible_col="control_eligible",
        n_boot=n_boot, seed=seed,
    )
    result.pop("match_pairs", None)
    result["policy_role"] = ("deployment_benefit_policy" if objective == "benefit"
                              else "diagnostic_counterfactual_switch_not_the_ABM_deployment_rule")
    result["grouping_status"] = model.manifest.get("grouping_status", "unverified")
    result["gwas_overlap_status"] = model.manifest.get("gwas_overlap_status", "unknown")
    result["validation_status"] = model.manifest.get("validation_status", "exploratory")
    return result


def fit(args):
    config = _load_json(args.config)
    gate = read_table(args.gate)
    calibration = read_table(args.head_calibration)
    policy_calibration = read_table(args.policy_calibration)
    manifest = dict(config["manifest"])
    manifest["input_file_sha256"] = {
        "gate": file_sha256(args.gate), "head_calibration": file_sha256(args.head_calibration),
        "policy_calibration": file_sha256(args.policy_calibration),
    }
    model = ABMAssessor.fit(
        gate, calibration, policy_calibration, manifest=manifest,
        task=config["task"], loss_name=config.get("loss_name"),
        head_options=config.get("head_options", {}), min_gain=config.get("min_gain", 0.0),
    )
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    model.save(target)
    model.heads.calibration_table().to_csv(target.with_suffix(".count_calibration.tsv"), sep="\t", index=False)
    if "y" in policy_calibration:
        # Labels only enter AFTER count heads and thresholds have been frozen.
        assessed = model.assess(policy_calibration.drop(columns="y"), manifest=manifest)
        assessed = assessed.merge(policy_calibration[["id", "y"]], on="id", how="left", validate="one_to_one")
        result = _evaluate(assessed, model, "benefit", model.primary_coverage,
                           args.bootstrap, model.heads.seed)
        result["stage"] = "held_out_policy_calibration_diagnostic"
        write_json(target.with_suffix(".policy_calibration.json"), result)
    print(json.dumps({"model": str(target), "sha256": file_sha256(target),
                      "count_bins": len(model.heads.bin_estimates)}, ensure_ascii=False))


def predict(args):
    model = ABMAssessor.load(args.model)
    manifest = _load_json(args.manifest)
    verify_feature_binding(args.input, args.input_provenance or str(args.input) + ".provenance.json",
                           manifest=manifest)
    source = read_table(args.input)
    # Passing a y column is an error even if every entry is missing.
    assessed = model.assess(source, manifest=manifest)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    assessed.to_csv(target, sep="\t", index=False)
    summary = model.dataset_assessment(assessed)
    summary["predictions_sha256"] = file_sha256(target)
    summary["model_sha256"] = file_sha256(args.model)
    write_json(target.with_suffix(".assessment.json"), summary)
    print(json.dumps({"predictions": str(target), "n_query": len(assessed),
                      "n_assessment_available": summary["n_assessment_available"]}, ensure_ascii=False))


def evaluate(args):
    model = ABMAssessor.load(args.model)
    sidecar = _load_json(args.predictions_manifest or str(Path(args.predictions).with_suffix(".assessment.json")))
    if sidecar["predictions_sha256"] != file_sha256(args.predictions):
        raise ValueError("Frozen predictions changed before evaluation")
    if sidecar["model_sha256"] != file_sha256(args.model):
        raise ValueError("Predictions were produced by a different frozen ABM model")
    predictions = read_table(args.predictions)
    if "y" in predictions:
        raise ValueError("Predictions must have been frozen before attaching outcomes")
    truth = read_table(args.truth)
    if not {"id", "y"}.issubset(truth.columns) or truth.id.isna().any() or truth.id.duplicated().any():
        raise ValueError("Truth must have unique id and y columns")
    unknown = set(truth.id) - set(predictions.id)
    if unknown:
        raise ValueError("Truth contains IDs outside the frozen query roster")
    frame = predictions.merge(truth[["id", "y"]], on="id", how="left", validate="one_to_one", sort=False)
    result = _evaluate(frame, model, args.objective, args.coverage, args.bootstrap, args.seed)
    result["predictions_sha256"] = file_sha256(args.predictions)
    result["truth_sha256"] = file_sha256(args.truth)
    result["model_sha256"] = file_sha256(args.model)
    write_json(args.output, result)
    print(json.dumps({"evaluation": args.output, "n_query": len(frame)}, ensure_ascii=False))


def support(args):
    donor = read_table(args.donor_metadata)
    query = read_table(args.query_metadata)
    if "y" in query:
        raise ValueError("Query support is computed without outcomes")
    donor_binding = verify_array_binding(args.donor_z, donor, args.donor_array_manifest or args.donor_z + ".rows.json")
    query_binding = verify_array_binding(args.query_z, query, args.query_array_manifest or args.query_z + ".rows.json")
    if donor_binding["geometry_sha256"] != query_binding["geometry_sha256"]:
        raise ValueError("Donor and query arrays use different frozen geometries")
    result = count_twins(
        np.load(args.donor_z, mmap_mode="r"), np.load(args.query_z, mmap_mode="r"),
        donor_ids=donor.id, query_ids=query.id, donor_groups=donor.group_id, query_groups=query.group_id,
        radius=args.radius, k_used=args.k_used, min_count=args.min_count,
        min_group_ess=args.min_group_ess, device=args.device,
        query_batch=args.query_batch, donor_batch=args.donor_batch,
        groups_verified=args.groups_verified,
    )
    result.table["group_id"] = query.group_id.to_numpy()
    result.table.to_csv(args.output, sep="\t", index=False)
    write_json(str(Path(args.output).with_suffix(".support.json")), {
        "support_sha256": file_sha256(args.output), "donor_binding": donor_binding,
        "query_binding": query_binding, "radius": args.radius, "k_used": args.k_used,
        "min_count": args.min_count, "min_group_ess": args.min_group_ess,
        "device_requested": args.device,
        "devices_used": sorted(result.table.device_used.unique().tolist()),
        "groups_verified": args.groups_verified,
    })
    # Neighbor arrays/IDs are individual-level artifacts. Do not publish them.
    if args.neighbors:
        np.savez_compressed(args.neighbors, indices=result.neighbor_indices,
                            distances=result.neighbor_distances)
    print(json.dumps({"support": args.output, "n_query": len(query)}, ensure_ascii=False))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    f = subs.add_parser("fit", help="Fit count heads and freeze selection using development roles only")
    for name in ("config", "gate", "head-calibration", "policy-calibration", "output"):
        f.add_argument("--" + name, required=True)
    f.add_argument("--bootstrap", type=int, default=1000)
    f.set_defaults(func=fit)
    p = subs.add_parser("predict", help="Label-free assessment and primary benefit policy")
    for name in ("model", "manifest", "input", "output"):
        p.add_argument("--" + name, required=True)
    p.set_defaults(func=predict)
    p.add_argument("--input-provenance", help="Default: INPUT.provenance.json, emitted by the verified feature adapter")
    e = subs.add_parser("evaluate", help="Attach outcomes only after predictions/selection are frozen")
    for name in ("model", "predictions", "truth", "output"):
        e.add_argument("--" + name, required=True)
    e.add_argument("--objective", choices=("predictability", "benefit"), default="benefit")
    e.add_argument("--coverage", type=float, choices=(0.2, 0.4, 0.5, 0.6, 0.8, 1.0), default=0.5)
    e.add_argument("--bootstrap", type=int, default=1000)
    e.add_argument("--seed", type=int, default=20261020)
    e.add_argument("--predictions-manifest", help="Default: prediction .assessment.json sidecar")
    e.set_defaults(func=evaluate)
    s = subs.add_parser("support", help="Count all fixed-radius donors with one CPU/CUDA definition")
    for name in ("donor-z", "query-z", "donor-metadata", "query-metadata", "output"):
        s.add_argument("--" + name, required=True)
    s.add_argument("--radius", type=float, required=True)
    s.add_argument("--k-used", type=int, default=128)
    s.add_argument("--min-count", type=int, default=20)
    s.add_argument("--min-group-ess", type=float, default=15)
    s.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    s.add_argument("--query-batch", type=int, default=128)
    s.add_argument("--donor-batch", type=int, default=4096)
    s.add_argument("--groups-verified", action="store_true")
    s.add_argument("--neighbors")
    s.add_argument("--donor-array-manifest")
    s.add_argument("--query-array-manifest")
    s.set_defaults(func=support)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
