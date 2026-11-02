"""Standalone prepare -> count twins -> evaluate workflow with stage receipts."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict
import fcntl
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import SCHEMA, VERSION
from .core import count_predict, freeze_policy, apply_policy, validate_outcomes
from .inputs import (file_hash, ids_hash, read_table, json_write, assign_roles,
                     load_export, load_production, Rows)

ROOT = Path(__file__).resolve().parents[2]
PHASES = ("prepare", "twins", "evaluate", "architecture", "final")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def code_identity():
    files = sorted(Path(__file__).parent.glob("*.py")) + [ROOT/"f/abm_core/representation.py", ROOT/"f/abm_core/support.py"]
    return digest({str(p.relative_to(ROOT)): file_hash(p) for p in files})


def config_for(args):
    path = Path(args.config)
    if path.is_dir():
        path = path/"count_abm.json"
    config = json.loads(path.read_text())
    if config.get("schema_version") != "grid.count_abm.config.1":
        raise ValueError("Use count_abm.json for the new counting workflow; historical abm.json belongs to the separate fitted-model workflow")
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
    return Path(args.analysis_root)/"count_abm"/trait


@contextmanager
def trait_lock(folder):
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/".workflow.lock").open("a") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Another count ABM job is using "+str(folder)) from None
        yield


def commit_stage(folder, stage, artifacts, *, config_hash=None, parent=None, details=None):
    record = dict(schema_version=SCHEMA, version=VERSION, stage=stage, status="complete",
                  code_sha256=code_identity(), config_sha256=config_hash,
                  parent_receipt_sha256=file_hash(parent) if parent else None,
                  artifacts={str(Path(p).relative_to(folder)): file_hash(p) for p in artifacts},
                  details=details or {})
    json_write(folder/(stage.upper()+"_COMPLETE.json"), record)


def verify_stage(folder, stage, *, config_hash=None, parent=None, verify_code=True):
    path = folder/(stage.upper()+"_COMPLETE.json")
    if not path.is_file():
        raise FileNotFoundError(f"Run grid.sh {stage} first: {path}")
    record = json.loads(path.read_text())
    if record.get("schema_version") != SCHEMA or record.get("status") != "complete" or record.get("stage") != stage:
        raise ValueError("Wrong or incomplete counting-workflow receipt: "+str(path))
    if verify_code and record["code_sha256"] != code_identity():
        raise ValueError("Counting code changed since this stage. Preserve the existing result and use a new analysis root: "+str(path))
    if config_hash is not None and record["config_sha256"] != config_hash:
        raise ValueError("Counting settings changed; preserve this result and use a new analysis root: "+str(path))
    if parent and record["parent_receipt_sha256"] != file_hash(parent):
        raise ValueError("Parent stage changed: "+str(path))
    for name, expected in record["artifacts"].items():
        target = folder/name
        if not target.is_file() or file_hash(target) != expected:
            raise ValueError("Frozen counting artifact missing or changed: "+str(target))
    return record


def prepare(args, config, trait):
    from abm_core.representation import BlockGeometry
    from abm_core.support import estimate_radius
    folder = trait_directory(args, trait)
    if (folder/"PREPARE_COMPLETE.json").is_file():
        verify_stage(folder, "prepare", config_hash=digest(config))
        frozen = json.loads((folder/"model.json").read_text())["metadata"]
        if args.input_manifest:
            requested = Path(str(args.input_manifest).replace("{trait}", trait)).resolve()
            if str(requested) != frozen.get("input_manifest") or file_hash(requested) != frozen.get("input_manifest_sha256"):
                raise ValueError("Requested source export differs from the frozen counting bank; use a new analysis root")
            load_export(requested, trait)  # Also recheck the explicitly requested producer-bound source files.
        elif args.prepared:
            requested = Path(args.prepared).resolve()
            if str(requested) != str(Path(frozen.get("input_source", "")).resolve()) or file_hash(requested) != frozen.get("input_sha256"):
                raise ValueError("Requested prepared cohort differs from the frozen counting bank; use a new analysis root")
        print(f"[GRID count] {trait}: verified prepared counting bank", flush=True)
        return
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
    donors, calibration = read_table(local/"donor.tsv.gz"), read_table(local/"calibration.tsv.gz")
    dz = np.load(local/"donor.z.npy", mmap_mode="r")
    cz = np.load(local/"calibration.z.npy", mmap_mode="r")
    radius = config.get("radius")
    if radius is None:
        radius = estimate_radius(dz, cz, donor_ids=donors.id, query_ids=calibration.id,
                                 donor_groups=donors.group_id, query_groups=calibration.group_id,
                                 k_ref=config["radius_k_ref"], quantile=config["radius_quantile"],
                                 device=args.device, query_batch=config["query_batch"], donor_batch=config["donor_batch"])
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Frozen radius must be positive and finite")
    features = features.copy()
    features["retained"] = np.isin(np.arange(len(features)), geometry.keep)
    features.to_csv(folder/"features.tsv", sep="\t", index=False)
    json_write(folder/"geometry.json", asdict(geometry))
    lineage = {role: dict(n=int(mask.sum()), groups=d.loc[mask, "group_id"].nunique(),
                          ordered_ids_sha256=ids_hash(d.loc[mask, "id"])) for role, mask in
               (("donor", donor), ("calibration", d.role.eq("calibration")), ("test", d.role.eq("test")))}
    model = dict(schema_version=SCHEMA, version=VERSION, trait=trait, task=task,
                 radius=float(radius), radius_fit_role="development_calibration_X_only" if "radius" not in config else "explicit_frozen_user_radius",
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
    json_write(folder/"model.json", model)
    artifacts += [folder/"features.tsv", folder/"geometry.json", folder/"model.json"]
    commit_stage(folder, "prepare", artifacts, config_hash=digest(config), details=dict(trait=trait, task=task, lineage=lineage))
    print(f"[GRID count] {trait}: {int(donor.sum()):,} donors, radius={radius:.6g}; no predictor fitted", flush=True)


def twins(args, config, trait):
    folder = trait_directory(args, trait)
    verify_stage(folder, "prepare", config_hash=digest(config))
    parent = folder/"PREPARE_COMPLETE.json"
    if (folder/"TWINS_COMPLETE.json").is_file():
        verify_stage(folder, "twins", config_hash=digest(config), parent=parent)
        print(f"[GRID count] {trait}: verified frozen count predictions", flush=True)
        return
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
        def sink(start, ids, offsets, indices, distances):
            path = trace_dir/f"part_{start:09d}.npz"
            np.savez_compressed(path, query_ids=ids, indptr=offsets,
                                donor_indices=indices, distances=distances)
            artifacts.append(path)
        result = count_predict(dz, qz, donor_ids=donor.id, query_ids=query.id,
                    donor_groups=donor.group_id, query_groups=query.group_id, donor_y=donor.y,
                    task=model["task"], radius=model["radius"], min_count=config["min_count"],
                    min_group_ess=config["min_group_ess"], device=args.device,
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
    commit_stage(folder, "twins", artifacts, config_hash=digest(config), parent=parent,
                 details=dict(query_outcomes_used=False, count_cap=None,
                              neighbor_weights="equal_one_per_eligible_donor", membership_scope="all_in_radius",
                              individual_artifacts="_LOCAL_ONLY", no_regression_fitted=True))
    print(f"[GRID count] {trait}: complete donor memberships and predictions frozen before evaluation", flush=True)


def evaluate(args, config, trait):
    from .report import evaluate as score, subgroup_table, plot_risk, write_trait_report
    folder = trait_directory(args, trait)
    verify_stage(folder, "prepare", config_hash=digest(config))
    verify_stage(folder, "twins", config_hash=digest(config), parent=folder/"PREPARE_COMPLETE.json")
    if (folder/"EVALUATE_COMPLETE.json").is_file():
        verify_stage(folder, "evaluate", parent=folder/"TWINS_COMPLETE.json")
        print(f"[GRID count] {trait}: verified evaluation", flush=True)
        return
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
    commit_stage(folder, "evaluate", artifacts, config_hash=digest(config), parent=folder/"TWINS_COMPLETE.json",
                 details=dict(truth_attached_after_frozen_predictions=True, validation_status=model["validation_status"]))
    print(f"[GRID count] {trait}: evaluation/report completed, all query rows retained", flush=True)


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


def architecture(args, config, traits):
    # Separate authorized plans can finish together. Protect the aggregate
    # directory as well as the existing per-trait scientific stage locks.
    with trait_lock(Path(args.analysis_root)/"count_abm/architecture"):
        return _architecture(args, config, traits)


def _architecture(args, config, traits):
    """Describe current frozen representation and optional posterior architecture.

    This function never searches historical final/paper summary tables.
    """
    out = Path(args.analysis_root)/"count_abm/architecture"
    out.mkdir(parents=True, exist_ok=True)
    rows, sources, posterior, top = [], {}, [], []
    for trait in traits:
        folder = trait_directory(args, trait)
        verify_stage(folder, "prepare", config_hash=digest(config))
        sources[trait] = file_hash(folder/"PREPARE_COMPLETE.json")
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
    commit_stage(out, "architecture", files, config_hash=digest(config),
                 details=dict(traits=traits, prepare_receipt_sha256=sources,
                              outcome_fit=False, historical_performance_tables_read=False))
    print("[GRID count] current representation architecture completed", flush=True)


def final(args, config, traits):
    with trait_lock(Path(args.analysis_root)/"count_abm/final"):
        return _final(args, config, traits)


def _final(args, config, traits):
    """Read only this pipeline's current, verified completion receipts."""
    out = Path(args.analysis_root)/"count_abm/final"
    out.mkdir(parents=True, exist_ok=True)
    perf, risk, statuses, sources = [], [], [], {}
    for trait in traits:
        folder = trait_directory(args, trait)
        if not (folder/"EVALUATE_COMPLETE.json").is_file():
            statuses.append(dict(trait=trait, status="missing_current_count_evaluation")); continue
        verify_stage(folder, "prepare", config_hash=digest(config))
        verify_stage(folder, "twins", parent=folder/"PREPARE_COMPLETE.json")
        verify_stage(folder, "evaluate", parent=folder/"TWINS_COMPLETE.json")
        model = json.loads((folder/"model.json").read_text())
        statuses.append(dict(trait=trait, status="complete", validation_status=model["validation_status"],
                             grouping_status=model["metadata"].get("grouping_status", "unverified"),
                             gwas_overlap_status=model["metadata"].get("gwas_overlap_status", "unverified")))
        perf.append(read_table(folder/"performance.tsv")); risk.append(read_table(folder/"risk_coverage.tsv"))
        sources[trait] = file_hash(folder/"EVALUATE_COMPLETE.json")
    status = pd.DataFrame(statuses)
    performance = pd.concat(perf, ignore_index=True) if perf else pd.DataFrame(columns=["trait", "role", "population", "primary_loss"])
    coverage = pd.concat(risk, ignore_index=True) if risk else pd.DataFrame(columns=["trait", "role", "actual_test_coverage", "primary_loss"])
    paths = []
    for name, table in (("status.tsv", status), ("performance.tsv", performance), ("risk_coverage.tsv", coverage)):
        path = out/name; table.to_csv(path, sep="\t", index=False); paths.append(path)
    from .benchmarks import compare_current
    comparisons, benchmark_status = compare_current(args.analysis_root, traits,
        repetitions=config["bootstrap_replicates"], seed=config["seed"])
    for name, table in (("count_vs_benchmarks.tsv", comparisons), ("benchmark_comparison_status.tsv", benchmark_status)):
        path = out/name; table.to_csv(path, sep="\t", index=False); paths.append(path)
    explanation = (
        "# GRID transparent count ABM\n\n"
        "Source family: `grid.count_abm.1`. Only completion-verified results under `count_abm/<trait>` enter these tables. "
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
    (out/"index.html").write_text("<!doctype html><html><head><meta charset='utf-8'><title>GRID count ABM</title></head><body>"
        "<h1>GRID transparent count ABM</h1><p>Current grid.count_abm.1 outputs; no historical fitted-model rows.</p>"+
        status.to_html(index=False)+"<h2>Performance</h2>"+performance.to_html(index=False)+
        "<h2>Coverage and composition</h2>"+coverage.to_html(index=False)+
        "<h2>Independent benchmark availability / exact-person comparability</h2>"+benchmark_status.to_html(index=False)+
        "<h2>Count vs benchmark on the same people</h2>"+comparisons.to_html(index=False)+"</body></html>")
    paths.append(out/"index.html")
    commit_stage(out, "final", paths, config_hash=digest(config),
                 details=dict(traits=traits, evaluate_receipt_sha256=sources,
                              fitted_any_model=False, historical_results_read=False,
                              benchmark_contract="grid.target_benchmark.1_COMPLETE_exact_test_ID_and_Y"))
    print("[GRID count] final current-source summary: "+str(out/"index.html"), flush=True)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=(*PHASES, "all"))
    p.add_argument("--traits", default="height,ldl,t2dm")
    from grid_paths import analysis_default
    p.add_argument("--analysis-root", default=analysis_default())
    p.add_argument("--config", default=str(ROOT/"config/count_abm.json"))
    p.add_argument("--block-config")
    p.add_argument("--prepared")
    p.add_argument("--input-manifest", help="Bound export JSON; optional {trait} path template; no fitted predictions")
    p.add_argument("--ld-block-dir")
    p.add_argument("--ld-block-files")
    p.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
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
    os.environ["GRID_OUTPUT_ROOT"] = args.analysis_root
    traits = [t.strip() for t in args.traits.split(",") if t.strip()]
    if not traits or len(traits) != len(set(traits)) or any(not t.replace("_", "").isalnum() for t in traits):
        raise ValueError("Use unique safe trait names")
    config = config_for(args)
    phases = PHASES if args.command == "all" else (args.command,)
    if args.dry_run:
        print(json.dumps(dict(schema_version=SCHEMA, phases=phases, traits=traits,
                              output_root=str(Path(args.analysis_root)/"count_abm"),
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
            else:
                for trait in traits:
                    with trait_lock(trait_directory(args, trait)):
                        globals()[phase](args, config, trait)


if __name__ == "__main__":
    main()
