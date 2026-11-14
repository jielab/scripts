#!/usr/bin/env python3
"""Standalone prepare -> count twins -> evaluate workflow with stage receipts."""
from __future__ import annotations

# Workflow
from settings import load_settings

from file_lock import flock as shared_flock
import argparse
import ast
from contextlib import contextmanager
from dataclasses import asdict
import fcntl
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from abm import SCHEMA, VERSION
from abm.count_core import count_predict, freeze_policy, apply_policy, validate_outcomes
from abm.count_inputs import (file_record, ids_hash, read_table, json_write, assign_roles,
                     load_export, load_production, Rows)

ROOT = Path(__file__).resolve().parents[2]
PHASES = ("prepare", "twins", "evaluate", "architecture", "final")


def config_for(args):
    path = Path(args.config)
    if path.is_dir():
        path = path/"abm.yml"
    config = load_settings(path)
    if config.get("schema_version") != "grid.abm.config.1":
        raise ValueError("Use config/abm.yml for ABM; benchmark methods use config/benchmarks.yml")
    block = Path(args.block_config or config["block_config"]).expanduser()
    if not args.block_config and not block.is_absolute():
        block = path.parent/block
    config["block_config"] = str(block.resolve())
    for name in ("radius_k_ref", "min_count", "min_group_ess", "query_batch", "donor_batch"):
        if getattr(args, name, None) is not None:
            config[name] = getattr(args, name)
    if args.radius is not None:
        config["radius"] = args.radius
    if args.bootstrap is not None:
        config["bootstrap_replicates"] = args.bootstrap
    for name in ("radius_k_ref", "min_count", "query_batch", "donor_batch"):
        if isinstance(config[name], bool) or not isinstance(config[name], int) or config[name] < 1:
            raise ValueError(name+" must be a positive integer")
    if not np.isfinite(config["min_group_ess"]) or config["min_group_ess"] <= 0:
        raise ValueError("min_group_ess must be positive")
    for name in ("sd_floor_quantile", "radius_quantile"):
        if not 0 <= config[name] <= 1:
            raise ValueError(name+" must be between zero and one")
    if not 0 < config["donor_fraction"] < 1:
        raise ValueError("donor_fraction must be between zero and one")
    if isinstance(config["bootstrap_replicates"], bool) or not isinstance(config["bootstrap_replicates"], int) or config["bootstrap_replicates"] < 0:
        raise ValueError("bootstrap_replicates must be a nonnegative integer")
    if not config["coverages"] or len(set(config["coverages"])) != len(config["coverages"]) or any(not 0 < c <= 1 for c in config["coverages"]):
        raise ValueError("Use distinct coverages in (0,1]")
    if "radius" in config and (not np.isfinite(config["radius"]) or config["radius"] <= 0):
        raise ValueError("radius must be positive and finite")
    if not 0 <= config["heterogeneity_quantile"] <= 1:
        raise ValueError("Invalid heterogeneity quantile")
    config["ld_block_files_override"] = args.ld_block_files
    config["ld_block_dir_override"] = args.ld_block_dir
    return config


def trait_directory(args, trait):
    return Path(args.analysis_root)/"tmp/abm"/trait


@contextmanager
def trait_lock(folder):
    ensure_output_directory(folder)
    with (folder/".workflow.lock").open("a+") as stream:
        try:
            shared_flock(stream, fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Another ABM job is using "+str(folder)) from None
        yield


def ensure_output_directory(folder, attempts=5):
    """Retry a disappearing directory entry, never overwrite a real file."""
    for attempt in range(attempts):
        try:
            folder.mkdir(parents=True, exist_ok=True)
            return
        except FileExistsError:
            if folder.is_dir():
                return
            if folder.exists() or folder.is_symlink() or attempt == attempts - 1:
                raise
            time.sleep(1)


def commit_stage(folder, stage, artifacts, *, details=None):
    record = dict(schema_version=SCHEMA, version=VERSION, stage=stage, status="complete",
                  artifacts={str(Path(p).relative_to(folder)): Path(p).stat().st_size for p in artifacts},
                  details=details or {})
    json_write(folder/(stage.upper()+"_COMPLETE.json"), record)


def verify_stage(folder, stage):
    path = folder/(stage.upper()+"_COMPLETE.json")
    if not path.is_file():
        raise FileNotFoundError(f"Run grid.sh {stage} first: {path}")
    record = json.loads(path.read_text())
    if record.get("schema_version") != SCHEMA or record.get("status") != "complete" or record.get("stage") != stage:
        raise ValueError("Incomplete ABM stage: " + str(path))
    if not record.get("artifacts"):
        raise ValueError("Empty ABM completion record: " + str(path))
    for name in record["artifacts"]:
        target = folder/name
        if not target.is_file() or target.stat().st_size == 0:
            raise ValueError("ABM artifact missing or empty: " + str(target))
    return record


def reusable_stage(args, folder, stage):
    if args.replace:
        return False
    try:
        verify_stage(folder, stage)
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def invalidate(folder, stage):
    order = ["geometry", "prepare", "twins", "evaluate"]
    for name in order[order.index(stage):]:
        (folder/(name.upper()+"_COMPLETE.json")).unlink(missing_ok=True)


def geometry(args, config, trait):
    """CPU-only cohort, geometry and normalized donor/query banks."""
    from abm.representation import BlockGeometry
    folder = trait_directory(args, trait)
    if reusable_stage(args, folder, "geometry"):
        try:
            local = folder/"_LOCAL_ONLY"
            for role in ("donor", "calibration", "test"):
                frame = read_table(local/(role+".tsv.gz"))
                values = np.load(local/(role+".z.npy"), mmap_mode="r", allow_pickle=False)
                if frame.id.duplicated().any() or len(frame) != len(values):
                    raise ValueError("Saved ABM row identities or array shape are invalid")
        except (OSError, ValueError, KeyError, AttributeError):
            pass  # An incomplete cached stage is rebuilt from its scientific inputs.
        else:
            print(f"[GRID ABM] {trait}: reuse complete geometry", flush=True)
            return
    invalidate(folder, "geometry")
    loaded = load_export(args.input_manifest, trait) if args.input_manifest else load_production(args, config, trait)
    frame, array, source_rows, features, task, metadata = loaded
    d = assign_roles(frame, donor_fraction=config["donor_fraction"], seed=config["seed"])
    donor = d.role.eq("donor") & d.y.notna()
    validate_outcomes(d.loc[donor, "y"], task)
    if donor.sum() < 2:
        raise ValueError("At least two observed training donors are needed")
    local = folder/"_LOCAL_ONLY"
    local.mkdir(parents=True, exist_ok=True, mode=0o700)
    local.chmod(0o700)
    geometry = BlockGeometry.fit(Rows(array, source_rows[donor.to_numpy()]), features.map_id,
                                 floor_quantile=config["sd_floor_quantile"])
    artifacts = []
    for role, mask in (("donor", donor), ("calibration", d.role.eq("calibration")), ("test", d.role.eq("test"))):
        rows = d.loc[mask, ["id", "group_id", "ancestry"]].copy()
        if role == "donor":
            rows["y"] = d.loc[mask, "y"].to_numpy()
        meta_path = local/(role+".tsv.gz")
        rows.to_csv(meta_path, sep="\t", index=False)
        zpath = local/(role+".z.npy")
        geometry.transform_to_npy(Rows(array, source_rows[mask.to_numpy()]), zpath, dtype="float64")
        artifacts += [meta_path, zpath]
        if role != "donor":
            truth = d.loc[mask, ["id", "y"]].copy()
            truth_path = local/(role+".truth.tsv.gz")
            truth.to_csv(truth_path, sep="\t", index=False)
            artifacts.append(truth_path)
    features = features.copy()
    features["retained"] = np.isin(np.arange(len(features)), geometry.keep)
    features.to_csv(folder/"features.tsv", sep="\t", index=False)
    json_write(folder/"geometry.json", asdict(geometry))
    lineage = {role: dict(n=int(mask.sum()), groups=d.loc[mask, "group_id"].nunique(),
                          ordered_ids_sha256=ids_hash(d.loc[mask, "id"])) for role, mask in
               (("donor", donor), ("calibration", d.role.eq("calibration")), ("test", d.role.eq("test")))}
    model = dict(schema_version=SCHEMA, version=VERSION, trait=trait, task=task,
                 geometry_fit_role="observed_training_donors_only", geometry_sha256=geometry.fingerprint(),
                 donor_filter="observed_outcome_only_no_residual_or_case_filtering",
                 prediction_rule="all_in_radius_equal_weight_mean", prior_strength=0,
                 fallback_rule="global_observed_donor_mean_explicitly_labelled",
                 global_prior=float(d.loc[donor, "y"].mean()),
                 query_outcomes_used_for_prediction=False, lineage=lineage, config=config,
                 development_missing_y_not_donors=int((d.role.eq("donor") & d.y.isna()).sum()),
                 metadata=metadata,
                 validation_status=config["validation_status"],
                 uncertainty="descriptive_local_outcomes_and_cluster_mean_interval_no_conditional_guarantee")
    json_write(folder/"geometry_model.json", model)
    artifacts += [folder/"features.tsv", folder/"geometry.json", folder/"geometry_model.json"]
    commit_stage(folder, "geometry", artifacts, details=dict(trait=trait, task=task, lineage=lineage))
    print(f"[GRID ABM] {trait}: CPU geometry frozen for {int(donor.sum()):,} donors", flush=True)


def radius(args, config, trait):
    """Estimate radius only from the frozen CPU-produced development banks."""
    from abm.support import estimate_radius
    folder = trait_directory(args, trait)
    receipt = verify_stage(folder, "geometry")
    parent = folder/"GEOMETRY_COMPLETE.json"
    if reusable_stage(args, folder, "prepare"):
        verify_stage(folder, "prepare")
        print(f"[GRID ABM] {trait}: verified frozen radius", flush=True)
        return
    invalidate(folder, "prepare")
    local = folder/"_LOCAL_ONLY"
    donors, calibration = read_table(local/"donor.tsv.gz"), read_table(local/"calibration.tsv.gz")
    dz = np.load(local/"donor.z.npy", mmap_mode="r")
    cz = np.load(local/"calibration.z.npy", mmap_mode="r")
    radius = config.get("radius")
    if radius is None:
        radius = estimate_radius(dz, cz, donor_ids=donors.id, query_ids=calibration.id,
                                 donor_groups=donors.group_id, query_groups=calibration.group_id,
                                 k_ref=config["radius_k_ref"], quantile=config["radius_quantile"],
                                 device="cpu", query_batch=config["query_batch"], donor_batch=config["donor_batch"])
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Frozen radius must be positive and finite")
    model = json.loads((folder/"geometry_model.json").read_text())
    model.update(radius=float(radius), radius_fit_role="development_calibration_X_only" if "radius" not in config else "explicit_frozen_user_radius")
    json_write(folder/"model.json", model)
    artifacts = [folder/name for name in receipt["artifacts"]]
    artifacts += [parent, folder/"model.json"]
    commit_stage(folder, "prepare", artifacts,
                 details=dict(trait=trait, task=model["task"], lineage=model["lineage"]))
    print(f"[GRID ABM] {trait}: radius={radius:.6g}; frozen geometry reused", flush=True)


def prepare(args, config, trait):
    geometry(args, config, trait)
    radius(args, config, trait)


def radius_twins(args, config, trait):
    radius(args, config, trait)
    twins(args, config, trait)


def twins(args, config, trait):
    folder = trait_directory(args, trait)
    verify_stage(folder, "prepare")
    parent = folder/"PREPARE_COMPLETE.json"
    if reusable_stage(args, folder, "twins"):
        verify_stage(folder, "twins")
        print(f"[GRID ABM] {trait}: verified frozen count predictions", flush=True)
        return
    invalidate(folder, "twins")
    model = json.loads((folder/"model.json").read_text())
    local = folder/"_LOCAL_ONLY"
    donor = read_table(local/"donor.tsv.gz")
    dz = np.load(local/"donor.z.npy", mmap_mode="r")
    artifacts = []
    for role in ("calibration", "test"):
        query = read_table(local/(role+".tsv.gz"))
        if "y" in query:
            raise ValueError("Query metadata must not include phenotypes")
        qz = np.load(local/(role+".z.npy"), mmap_mode="r")
        trace_dir = local/(role+".memberships")
        trace_dir.mkdir(exist_ok=True)
        for old_part in trace_dir.glob("part_*.npz"):
            old_part.unlink()
        def sink(start, ids, offsets, indices, distances):
            path = trace_dir/f"part_{start:09d}.npz"
            np.savez_compressed(path, query_ids=ids, indptr=offsets,
                                donor_indices=indices, distances=distances)
            artifacts.append(path)
        result = count_predict(dz, qz, donor_ids=donor.id, query_ids=query.id,
                    donor_groups=donor.group_id, query_groups=query.group_id, donor_y=donor.y,
                    task=model["task"], radius=model["radius"], min_count=config["min_count"],
                    min_group_ess=config["min_group_ess"], device="cpu",
                    query_batch=config["query_batch"], donor_batch=config["donor_batch"], trace_sink=sink,
                    groups_verified=model["metadata"].get("groups_verified", False))
        if role == "calibration":
            policy = freeze_policy(result.table, config["coverages"], config["heterogeneity_quantile"])
            json_write(folder/"policy.json", policy)
            artifacts.append(folder/"policy.json")
        predictions = apply_policy(result.table, policy)
        predictions["ancestry"] = query.ancestry.to_numpy()
        path = local/(role+".predictions.tsv.gz")
        predictions.to_csv(path, sep="\t", index=False)
        artifacts.append(path)
    commit_stage(folder, "twins", artifacts,
                 details=dict(query_outcomes_used=False, count_cap=None,
                              neighbor_weights="equal_one_per_eligible_donor", membership_scope="all_in_radius",
                              individual_artifacts="_LOCAL_ONLY", no_regression_fitted=True))
    print(f"[GRID ABM] {trait}: complete donor memberships and predictions frozen before evaluation", flush=True)


def evaluate(args, config, trait):
    from abm.count_report import evaluate as score, subgroup_table, plot_risk, write_trait_report
    folder = trait_directory(args, trait)
    verify_stage(folder, "prepare")
    verify_stage(folder, "twins")
    if reusable_stage(args, folder, "evaluate"):
        verify_stage(folder, "evaluate")
        print(f"[GRID ABM] {trait}: verified evaluation", flush=True)
        return
    invalidate(folder, "evaluate")
    model = json.loads((folder/"model.json").read_text())
    local = folder/"_LOCAL_ONLY"
    perf, risk, groups, artifacts = [], [], [], []
    for role in ("calibration", "test"):
        predictions = read_table(local/(role+".predictions.tsv.gz"))
        truth = read_table(local/(role+".truth.tsv.gz"))
        if "y" in predictions or set(truth.id) != set(predictions.id) or truth.id.duplicated().any():
            raise ValueError("Evaluation truth must match the frozen query roster exactly")
        frame = predictions.merge(truth, on="id", how="left", validate="one_to_one", sort=False)
        p, r = score(frame, task=model["task"], trait=trait, role=role,
                     repetitions=config["bootstrap_replicates"], seed=config["seed"])
        perf.append(p); risk.append(r)
        groups.append(subgroup_table(frame, task=model["task"], trait=trait, role=role))
        # Counterexamples remain local, and are explicitly selected after seeing Y.
        examples = frame.copy()
        examples["squared_error"] = (examples.y-examples.p_population)**2
        examples["selection_reason"] = "largest_observed_error_after_unblinding_not_validation_examples"
        examples = examples.sort_values("squared_error", ascending=False).groupby("evidence_state", sort=False).head(20)
        path = local/(role+".counterexamples.tsv.gz")
        examples.to_csv(path, sep="\t", index=False); artifacts.append(path)
        artifacts += reference_examples(folder, role, examples, model["task"])
    performance, coverage, ancestry = pd.concat(perf, ignore_index=True), pd.concat(risk, ignore_index=True), pd.concat(groups, ignore_index=True)
    for name, table in (("performance.tsv", performance), ("risk_coverage.tsv", coverage), ("ancestry_audit.tsv", ancestry)):
        path = folder/name; table.to_csv(path, sep="\t", index=False); artifacts.append(path)
    plot_risk(coverage, folder/"risk_coverage.png")
    write_trait_report(folder/"report.html", model, performance, coverage, ancestry)
    artifacts += [folder/"risk_coverage.png", folder/"report.html"]
    commit_stage(folder, "evaluate", artifacts,
                 details=dict(truth_attached_after_frozen_predictions=True, validation_status=model["validation_status"]))
    print(f"[GRID ABM] {trait}: evaluation/report completed, all query rows retained", flush=True)


def reference_examples(folder, role, examples, task):
    """Inspectable supporting/contrasting references; never modifies prediction.

    Example queries are explicitly post-outcome diagnostics. Each selected query
    gets nearest references on both outcome sides, when present. All references,
    not only these displayed examples, remain in the complete membership files.
    Feature distances describe matching; they are not causal outcome effects.
    """
    local = folder/"_LOCAL_ONLY"
    donor = read_table(local/"donor.tsv.gz")
    query = read_table(local/(role+".tsv.gz"))
    dz, qz = np.load(local/"donor.z.npy", mmap_mode="r"), np.load(local/(role+".z.npy"), mmap_mode="r")
    features = read_table(folder/"features.tsv")
    features = features.loc[features.retained].reset_index(drop=True)
    selected = set(examples.id)
    lookup = pd.Series(np.arange(len(query)), index=query.id)
    prediction = examples.set_index("id")
    records, blocks = [], []
    for path in sorted((local/(role+".memberships")).glob("part_*.npz")):
        with np.load(path, allow_pickle=False) as part:
            for j, query_id in enumerate(part["query_ids"]):
                if query_id not in selected:
                    continue
                start, end = part["indptr"][j:j+2]
                indices, distances = part["donor_indices"][start:end], part["distances"][start:end]
                if not len(indices):
                    continue
                outcomes = donor.iloc[indices].y.to_numpy(float)
                sides = ([("case_reference", outcomes == 1), ("control_reference", outcomes == 0)] if task == "binary" else
                         [("at_or_above_neighbor_mean", outcomes >= outcomes.mean()),
                          ("below_neighbor_mean", outcomes < outcomes.mean())])
                qindex = int(lookup.loc[query_id])
                for label, mask in sides:
                    available = np.flatnonzero(mask)
                    order = available[np.lexsort((donor.iloc[indices[available]].id.to_numpy(str), distances[available]))[:2]]
                    for index in order:
                        donor_index = int(indices[index])
                        record = dict(query_id=str(query_id), donor_id=donor.id.iloc[donor_index],
                                      donor_group_id=donor.group_id.iloc[donor_index],
                                      donor_y=float(outcomes[index]), distance=float(distances[index]),
                                      reference_type=label, query_n_twins=int(len(indices)),
                                      query_count_prediction=float(outcomes.mean()),
                                      query_evidence_state=prediction.loc[query_id, "evidence_state"],
                                      query_selection="post_outcome_largest_error_diagnostic")
                        records.append(record)
                        delta = np.asarray(qz[qindex])-np.asarray(dz[donor_index])
                        square = delta**2
                        total = square.sum()
                        for k in np.argsort(-square, kind="stable")[:6]:
                            blocks.append(dict(query_id=str(query_id), donor_id=record["donor_id"],
                                feature_id=features.feature_id.iloc[k], map_id=features.map_id.iloc[k],
                                standardized_query_coordinate=float(qz[qindex, k]),
                                standardized_donor_coordinate=float(dz[donor_index, k]),
                                squared_distance_contribution=float(square[k]),
                                fraction_of_pair_distance=float(square[k]/total) if total else 0.,
                                interpretation="matching_distance_not_causal_outcome_effect"))
    ref_path, block_path = local/(role+".reference_examples.tsv.gz"), local/(role+".block_examples.tsv.gz")
    pd.DataFrame(records, columns=["query_id", "donor_id", "donor_group_id", "donor_y", "distance", "reference_type",
        "query_n_twins", "query_count_prediction", "query_evidence_state", "query_selection"]).to_csv(ref_path, sep="\t", index=False)
    pd.DataFrame(blocks, columns=["query_id", "donor_id", "feature_id", "map_id", "standardized_query_coordinate",
        "standardized_donor_coordinate", "squared_distance_contribution", "fraction_of_pair_distance", "interpretation"]).to_csv(block_path, sep="\t", index=False)
    return [ref_path, block_path]


def report_directory(args, stage):
    root = Path(args.analysis_root)/"tmp/abm"
    trait = getattr(args, "report_trait", None)
    return root/trait/stage if trait else root/stage


def architecture(args, config, traits):
    # Separate authorized plans can finish together. Protect the aggregate
    # directory as well as the existing per-trait scientific stage locks.
    with trait_lock(report_directory(args, "architecture")):
        return _architecture(args, config, traits)


def _architecture(args, config, traits):
    """Describe current frozen representation and optional posterior architecture.

    This function never searches historical final/paper summary tables.
    """
    out = report_directory(args, "architecture")
    out.mkdir(parents=True, exist_ok=True)
    rows, sources, posterior, top = [], {}, [], []
    for trait in traits:
        folder = trait_directory(args, trait)
        verify_stage(folder, "prepare")
        sources[trait] = file_record(folder/"PREPARE_COMPLETE.json")
        geometry = json.loads((folder/"geometry.json").read_text())
        model = json.loads((folder/"model.json").read_text())
        for label in sorted(set(geometry["input_maps"])):
            rows.append(dict(trait=trait, map_id=label,
                             input_features=geometry["input_maps"].count(label),
                             retained_features=geometry["retained_maps"].count(label),
                             sd_floor=geometry["sd_floors"].get(label),
                             radius=model["radius"], geometry_fit_role=model["geometry_fit_role"]))
        if model["metadata"]["source_kind"] == "producer_bound_raw_csx_ldblocks":
            import block_inputs as b
            from regional_dimensionality import population_metrics
            bc = b.read_config(config["block_config"])
            layout = b.resolve_block_files(args.ld_block_dir or bc["ld_block_dir"], args.ld_block_files or bc["ld_block_files"])
            tokens, _, _ = b.build_manifest(layout, bc)
            files, _ = b.locate_weights(bc, trait)
            for pop in b.POPS:
                weights = b.read_weights(files[pop])
                bed = tokens[tokens.source_population.eq(pop if layout["mode"] == "population_multiview" else "ALL")].copy()
                for definition in ("beta_squared", "abs_beta"):
                    metrics, blocks = population_metrics(weights, bed, definition)
                    posterior.append(dict(trait=trait, population=pop, mass_definition=definition, **metrics))
                    top.append(blocks.assign(trait=trait, population=pop, mass_definition=definition))
    pd.DataFrame(rows).to_csv(out/"geometry_by_trait_map.tsv", sep="\t", index=False)
    files = [out/"geometry_by_trait_map.tsv"]
    if posterior:
        pd.DataFrame(posterior).to_csv(out/"posterior_regional_dimensionality.tsv", sep="\t", index=False)
        pd.concat(top, ignore_index=True).to_csv(out/"posterior_top_blocks.tsv", sep="\t", index=False)
        files += [out/"posterior_regional_dimensionality.tsv", out/"posterior_top_blocks.tsv"]
    commit_stage(out, "architecture", files,
                 details=dict(traits=traits, prepare_receipt_file_record=sources,
                              outcome_fit=False, historical_performance_tables_read=False))
    print("[GRID ABM] current representation architecture completed", flush=True)


def final(args, config, traits):
    with trait_lock(report_directory(args, "final")):
        return _final(args, config, traits)


def _final(args, config, traits):
    """Read only this pipeline's current, verified completion receipts."""
    out = report_directory(args, "final")
    out.mkdir(parents=True, exist_ok=True)
    perf, risk, statuses, sources = [], [], [], {}
    for trait in traits:
        folder = trait_directory(args, trait)
        if not (folder/"EVALUATE_COMPLETE.json").is_file():
            statuses.append(dict(trait=trait, status="missing_current_count_evaluation")); continue
        verify_stage(folder, "prepare")
        verify_stage(folder, "twins")
        verify_stage(folder, "evaluate")
        model = json.loads((folder/"model.json").read_text())
        statuses.append(dict(trait=trait, status="complete", validation_status=model["validation_status"],
                             grouping_status=model["metadata"].get("grouping_status", "unverified"),
                             gwas_overlap_status=model["metadata"].get("gwas_overlap_status", "unverified")))
        perf.append(read_table(folder/"performance.tsv")); risk.append(read_table(folder/"risk_coverage.tsv"))
        sources[trait] = file_record(folder/"EVALUATE_COMPLETE.json")
    status = pd.DataFrame(statuses)
    performance = pd.concat(perf, ignore_index=True) if perf else pd.DataFrame(columns=["trait", "role", "population", "primary_loss"])
    coverage = pd.concat(risk, ignore_index=True) if risk else pd.DataFrame(columns=["trait", "role", "actual_test_coverage", "primary_loss"])
    paths = []
    for name, table in (("status.tsv", status), ("performance.tsv", performance), ("risk_coverage.tsv", coverage)):
        path = out/name; table.to_csv(path, sep="\t", index=False); paths.append(path)
    from abm.count_report import compare_current
    comparisons, benchmark_status = compare_current(args.analysis_root, traits,
        repetitions=config["bootstrap_replicates"], seed=config["seed"])
    for name, table in (("count_vs_benchmarks.tsv", comparisons), ("benchmark_comparison_status.tsv", benchmark_status)):
        path = out/name; table.to_csv(path, sep="\t", index=False); paths.append(path)
    explanation = (
        "# GRID ABM\n\n"
        "Source family: `grid.abm.1`. Only completion-verified results under `tmp/abm/<trait>` enter these tables. "
        "The historical blocknet / CountHeads pipeline, benchmark fits, and old `final` tables are not merged.\n\n"
        "Predictions average every observed training donor inside a training-frozen genetic radius. "
        "Binary predictions are the number of donor cases divided by the number of all donors. "
        "The global donor-mean fallback is explicitly labelled and evaluated on the full query population. "
        "Support does not guarantee correctness. Read actual coverage, case/control retention, ancestry audit, "
        "and common-risk overlap together with error. All individual evidence remains local-only.\n\n"
        "The separate count_vs_benchmarks table reads only the new independent benchmark COMPLETE contracts. "
        "Comparisons require the exact same test IDs and labels (including missingness); no smaller intersection is taken. "
        "It evaluates each method on the same full roster and the same count-supported people, with paired group bootstrap. "
        "Missing/incompatible benchmarks are explicitly listed, never replaced by old aggregates. "
        "Training protocols and budgets remain separately labelled.\n\n"
    )
    (out/"README.md").write_text(explanation+"```tsv\n"+status.to_csv(sep="\t", index=False)+"```\n")
    paths.append(out/"README.md")
    # No private identities or donor edges are embedded in the public index.
    (out/"index.html").write_text("<!doctype html><html><head><meta charset='utf-8'><title>GRID ABM</title></head><body>"
        "<h1>GRID ABM</h1><p>Current grid.abm.1 outputs; no historical fitted-model rows.</p>"+
        status.to_html(index=False)+"<h2>Performance</h2>"+performance.to_html(index=False)+
        "<h2>Coverage and composition</h2>"+coverage.to_html(index=False)+
        "<h2>Independent benchmark availability / exact-person comparability</h2>"+benchmark_status.to_html(index=False)+
        "<h2>Count vs benchmark on the same people</h2>"+comparisons.to_html(index=False)+"</body></html>")
    paths.append(out/"index.html")
    commit_stage(out, "final", paths,
                 details=dict(traits=traits, evaluate_receipt_file_record=sources,
                              fitted_any_model=False, historical_results_read=False,
                              benchmark_contract="grid.target_benchmark.1_COMPLETE_exact_test_ID_and_Y"))
    print("[GRID ABM] final current-source summary: "+str(out/"index.html"), flush=True)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=(*PHASES, "geometry", "radius", "radius_twins", "all"))
    p.add_argument("--traits", default="height,ldl,t2dm")
    from grid_paths import analysis_default
    p.add_argument("--analysis-root", default=analysis_default())
    p.add_argument("--config", default=str(ROOT/"config/abm.yml"))
    p.add_argument("--block-config")
    p.add_argument("--prepared")
    p.add_argument("--input-manifest", help="Bound export JSON; optional {trait} path template; no fitted predictions")
    p.add_argument("--ld-block-dir")
    p.add_argument("--ld-block-files")
    p.add_argument("--device", choices=("cpu",), default="cpu", help=argparse.SUPPRESS)
    p.add_argument("--replace", choices=("TRUE", "FALSE"), default=os.getenv("GRID_REPLACE", "TRUE"))
    p.add_argument("--report-trait", help="Write reports inside this trait's directory")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--radius", type=float, help="Positive predeclared radius; otherwise training-calibration X only")
    p.add_argument("--radius-k-ref", type=int)
    p.add_argument("--min-count", type=int)
    p.add_argument("--min-group-ess", type=float)
    p.add_argument("--query-batch", type=int)
    p.add_argument("--donor-batch", type=int)
    p.add_argument("--bootstrap", type=int)
    p.add_argument("--dry-run", action="store_true")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    args.replace = args.replace == "TRUE"
    os.environ["GRID_OUTPUT_ROOT"] = args.analysis_root
    traits = [t.strip() for t in args.traits.split(",") if t.strip()]
    if not traits or len(traits) != len(set(traits)) or any(not t.replace("_", "").isalnum() for t in traits):
        raise ValueError("Use unique safe trait names")
    if args.report_trait and traits != [args.report_trait]:
        raise ValueError("--report-trait must match the single requested trait")
    config = config_for(args)
    phases = PHASES if args.command == "all" else (args.command,)
    if args.dry_run:
        print(json.dumps(dict(schema_version=SCHEMA, phases=phases, traits=traits,
                              output_root=str(Path(args.analysis_root)/"tmp/abm"),
                              fits_regression=False, trains_prsformer=False, trains_blocknet=False,
                              prediction_rule=config["prediction_rule"],
                              donor_policy="observed_training_outcomes_no_error_pruning",
                              source=args.input_manifest or "prepared_cohort_and_bound_CSx_blocks",
                              block_config=config["block_config"]), indent=2))
        return
    if args.threads < 1:
        raise ValueError("threads must be positive")
    os.umask(0o077)
    from threadpoolctl import threadpool_limits
    from kinship import analysis_guard
    with threadpool_limits(limits=args.threads), analysis_guard(args.analysis_root):
        for phase in phases:
            if phase in ("architecture", "final"):
                globals()[phase](args, config, traits)
                from output_layout import publish_count_report
                publish_count_report(args.analysis_root, phase, report_directory(args, phase), args.report_trait)
            else:
                for trait in traits:
                    with trait_lock(trait_directory(args, trait)):
                        globals()[phase](args, config, trait)
                        from output_layout import publish_count_stage
                        publish_count_stage(args.analysis_root, trait, phase, trait_directory(args, trait))



if __name__ == "__main__":
    main()
