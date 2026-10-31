"""Production GRID FINAL ABM: seven content-addressed, independently resumable phases."""
from __future__ import annotations

import argparse
from contextlib import ExitStack
from dataclasses import asdict
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import shutil
import sys
import time

import joblib
import numpy as np
import pandas as pd

from .grid_inputs import (load_grid, locate_blocks, block_maps, IndexedBlocks, digest,
                          smoke_subset, role_manifest, feature_manifest, validate_lineage)
from .contracts import ordered_ids_sha256, write_array_manifest, write_feature_manifest
from .losses import row_loss
from .model import ABMAssessor, file_sha256, write_json, COVERAGES
from .representation import BlockGeometry
from .support import count_twins, estimate_radius, _backend
from .predictors import fit_predictors
from . import reporting
from runtime_state import file_identity, artifact_matches, group_contract

ROOT = Path(__file__).resolve().parents[2]
PHASES = ("prepare", "fit_predictors", "build_support", "fit_assessor",
          "freeze_predictions", "evaluate", "report")


def log(message):
    print("[GRID ABM] " + message, flush=True)


def source_manifest():
    paths = [p for p in (ROOT / "f").rglob("*") if p.suffix in (".py", ".R", ".sh")]
    paths += list((ROOT / "config").glob("*.json")) + [ROOT / "grid.sh"]
    paths += [p for p in (ROOT.parent / "0f").rglob("*") if p.suffix in (".py", ".R", ".sh")]
    return {str(p.relative_to(ROOT.parent)): file_sha256(p) for p in sorted(paths)}


class Stages:
    def __init__(self, run, contract):
        self.run = run
        self.path = run / "PHASES.json"
        self.contract = contract
        self.state = json.loads(self.path.read_text()) if self.path.exists() else {"contract": contract, "phases": {}}
        if self.state["contract"] != contract:
            raise ValueError("Run contract differs from frozen dependencies")
        self.runtime = []

    def perform(self, name, fn, device, n):
        old = self.state["phases"].get(name)
        if old:
            for path, sha in old.items():
                p = self.run / path
                if not artifact_matches(p, sha):
                    raise ValueError(f"Changed phase artifact {name}/{path}; restore or use --new-run")
            log("REUSE " + name)
            self.runtime_row(name, device, n, 0., reused=True)
            return
        log("START " + name + " device=" + device)
        start = time.monotonic()
        import torch
        if device == "cuda" and torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
        files = fn()
        if not files:
            raise ValueError("Empty phase receipt")
        self.state["phases"][name] = {str(p.relative_to(self.run)): file_identity(p) for p in files}
        write_json(self.path, self.state)
        elapsed = time.monotonic() - start
        self.runtime_row(name, device, n, elapsed)
        log(f"DONE {name}: {elapsed:.1f}s")

    def runtime_row(self, name, device, n, elapsed, reused=False):
        import torch
        details = getattr(self, "execution_details", {}).get(name, {})
        rp = self.run / "runtime.tsv"
        prior = pd.read_csv(rp, sep="\t") if rp.exists() else pd.DataFrame()
        if reused and not details and len(prior) and "executed_device" in prior:
            actual = prior.loc[prior.phase.eq(name) & prior.executed_device.isin(["cpu", "cuda"]), "executed_device"]
            cached = (prior.loc[prior.phase.eq(name), "original_cache_device"].dropna()
                      if "original_cache_device" in prior else pd.Series(dtype=str))
            details = dict(original_cache_device=actual.iloc[0] if len(actual) else
                           cached.iloc[0] if len(cached) else "unrecorded")
        reused = reused or details.get("cache_reused", False)
        row = dict(phase=name, requested_device=self.contract.get("config", {}).get("execution", {}).get("device", device),
                   executed_device="none_cache_reused" if reused else device,
                   cache_reused=reused, original_cache_device=details.get("original_cache_device"),
                   seconds=elapsed, n=n, rows_per_second=None if reused else n / max(elapsed, 1e-9),
                   max_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2,
                   peak_vram_gib=torch.cuda.max_memory_allocated() / 1024 ** 3 if device == "cuda" and not reused and torch.cuda.is_available() else 0,
                   gpu=torch.cuda.get_device_name() if device == "cuda" and not reused and torch.cuda.is_available() else None,
                   dtype="float64_distance_float32_network")
        reporting.save_table(pd.DataFrame([*prior.to_dict("records"), row]), rp)


def prepare_support(frame, array, rows, maps, receipt, cache, config, device, groups_verified):
    """Support shared by blocknet and MoE; no predictor output or query Y is a key."""
    from .storage import support_lock, release_verified_scratch, SCRATCH, RELEASE_FILE, busy_reason
    cache = Path(cache)
    with support_lock(cache):
        # Legacy processes may still be using this directory without the lock.
        reason = busy_reason([cache / name for name in (*SCRATCH, "donor_geometry.npy")])
        if reason:
            raise BlockingIOError(f"Support cache in use: {cache}: {reason}")
        if not (cache / "COMPLETE.json").exists():
            (cache / RELEASE_FILE).unlink(missing_ok=True)
        try:
            with ExitStack() as mappings:
                result = _prepare_support(frame, array, rows, maps, receipt, cache, config,
                                          device, groups_verified, mappings)
        finally:
            # Close all mappings before unlink (also required on Windows/WSL).
            # Failed support has no reusable geometry; retain completed donor banks.
            if not (cache / "COMPLETE.json").exists():
                for name in (*SCRATCH, "donor_geometry.npy"):
                    (cache / name).unlink(missing_ok=True)
        release_verified_scratch(cache, result[1])
        return result


def _prepare_support(frame, array, rows, maps, receipt, cache, config, device, groups_verified, mappings):
    from .storage import validate_support
    from resume import atomic_json
    cache.mkdir(parents=True, exist_ok=True)
    pointer = cache / "COMPLETE.json"
    if pointer.is_file():
        saved = validate_support(cache)
        log("REUSE shared fixed-geometry support")
        support = pd.read_parquet(cache / "support.parquet")
        support["groups_verified"] = groups_verified
        support.loc[support.supported, "support_reason"] = (
            "supported_verified_groups" if groups_verified else "supported_provisional_groups")
        saved = dict(saved, bank=dict(saved["bank"], groups_verified=groups_verified),
                     cache_reused=True)
        return support, saved
    opts, execution = config["support"], config["execution"]
    valid = np.zeros(len(frame), bool)
    for start in range(0, len(frame), 512):
        block = np.asarray(array[rows[start:start + 512]])
        valid[start:start + len(block)] = np.isfinite(block).reshape(len(block), -1).all(1)
    donors = frame.role.eq("build").to_numpy() & valid
    tune = frame.role.eq("tune_model").to_numpy() & valid
    if not donors.any() or not tune.any():
        raise ValueError("No QC-valid build bank or tune_model geometry")
    geometry = BlockGeometry.fit(IndexedBlocks(array, rows[donors]), maps,
                                  floor_quantile=opts["sd_floor_quantile_within_map"])
    geometry.representation_id = opts["representation"]
    geometry.transform_to_npy(IndexedBlocks(array, rows[valid]), cache / "geometry.npy", dtype="float64")
    qframe = frame.loc[valid, ["id", "group_id", "role"]].reset_index(drop=True)
    z = np.load(cache / "geometry.npy", mmap_mode="r")
    mappings.callback(z._mmap.close)
    bank_ix = np.flatnonzero(qframe.role.eq("build"))
    donor_z = np.lib.format.open_memmap(cache / "donor_geometry.npy", mode="w+", dtype="float64",
                                       shape=(len(bank_ix), z.shape[1]))
    mappings.callback(donor_z._mmap.close)
    for start in range(0, len(bank_ix), 512):
        donor_z[start:start + 512] = z[bank_ix[start:start + 512]]
    donor_z.flush()
    bank = qframe.iloc[bank_ix].reset_index(drop=True)
    radius_ix = np.flatnonzero(qframe.role.eq("tune_model"))
    # The radius role is modest; still export by bounded chunks to avoid a dense copy.
    rz = np.lib.format.open_memmap(cache / "radius_geometry.npy", mode="w+", dtype="float64",
                                  shape=(len(radius_ix), z.shape[1]))
    mappings.callback(rz._mmap.close)
    for start in range(0, len(radius_ix), 512):
        rz[start:start + 512] = z[radius_ix[start:start + 512]]
    rz.flush()
    search = dict(device=device, query_batch=execution["query_batch"], donor_batch=execution["donor_batch"])
    rframe = qframe.iloc[radius_ix]
    log(f"Radius: {len(bank):,} donors × {len(rframe):,} tune_model queries; X only")
    radius = estimate_radius(donor_z, rz, donor_ids=bank.id, donor_groups=bank.group_id,
                             query_ids=rframe.id, query_groups=rframe.group_id,
                             k_ref=opts["radius_k_ref"], quantile=opts["radius_quantile"], **search)
    rz._mmap.close()
    (cache / "radius_geometry.npy").unlink()
    log(f"Frozen radius={radius:.12g}; counting all {len(qframe):,} valid queries")
    result = count_twins(donor_z, z, donor_ids=bank.id, donor_groups=bank.group_id,
                         query_ids=qframe.id, query_groups=qframe.group_id, radius=radius,
                         k_used=opts["k_used"], min_count=opts["min_count"],
                         min_group_ess=opts["min_group_ess"], groups_verified=groups_verified, **search)
    support = frame[["id"]].merge(result.table, on="id", how="left", validate="one_to_one")
    bad = support.supported.isna()
    support.loc[bad, "supported"] = False
    support.loc[bad, "groups_verified"] = groups_verified
    support.loc[bad, "support_reason"] = "genotype_qc_nonfinite_blocks"
    support["supported"] = support.supported.astype(bool)
    support["groups_verified"] = support.groups_verified.astype(bool)
    support.to_parquet(cache / "support.parquet", index=False)
    bank.to_parquet(cache / "bank.parquet", index=False)
    qframe.to_parquet(cache / "geometry_rows.parquet", index=False)
    np.save(cache / "neighbor_indices.npy", result.neighbor_indices)
    np.save(cache / "neighbor_distances.npy", result.neighbor_distances)
    qbinding = write_array_manifest(cache / "geometry.npy", qframe, geometry_sha256=geometry.fingerprint())
    dbinding = write_array_manifest(cache / "donor_geometry.npy", bank, geometry_sha256=geometry.fingerprint())
    bank_info = dict(n=len(bank), qc="finite_block_profile_and_observed_trait_build_label",
                     pruning=False, groups_verified=groups_verified,
                     ids_sha256=ordered_ids_sha256(bank.id), groups_sha256=ordered_ids_sha256(bank.group_id))
    geometry_info = dict(state=asdict(geometry), geometry_sha256=geometry.fingerprint(),
                         representation_input=receipt, radius=radius, radius_role="tune_model_X_only",
                         absent_maps=sorted(set(maps) - set(geometry.retained_maps)),
                         dropped_blocks=sorted(set(range(len(maps))) - set(geometry.keep)),
                         storage_dtype="float64", map_weighting="equal_map_mean_squared_distance")
    write_json(cache / "geometry_manifest.json", geometry_info)
    write_json(cache / "donor_bank_manifest.json", bank_info)
    names = ["support.parquet", "bank.parquet", "geometry_rows.parquet",
             "geometry.npy", "donor_geometry.npy", "neighbor_indices.npy", "neighbor_distances.npy",
             "geometry_manifest.json", "donor_bank_manifest.json", qbinding.name, dbinding.name]
    saved = dict(files={p: file_sha256(cache / p) for p in names}, bank=bank_info, geometry=geometry_info,
                 executed_device=device, cache_reused=False)
    atomic_json(pointer, saved)
    return support, saved


def predictor_identity(frame, features, csx_columns, candidate, cfg, block_receipt, sources):
    """Only mathematical fit/prediction dependencies; assessment Y/loss is excluded."""
    used = frame.role.isin(["build", "tune_model", "prediction_calibration"])
    return dict(schema_version="abm.predictors.pooled.1", candidate=candidate,
        features=list(features), csx_columns=list(csx_columns), config=cfg,
        ids=ordered_ids_sha256(frame.id), groups=ordered_ids_sha256(frame.group_id),
        roles=ordered_ids_sha256(frame.role), side_values=hashlib_array(frame[features].to_numpy(float)),
        supervised_labels=hashlib_array(frame.loc[used, "y"].to_numpy(float)),
        blocks=block_receipt["array_binding"] if candidate == "blocknet" else None,
        sources={k: v for k, v in sources.items() if Path(k).name in
                 {"predictors.py", "comparators.py", "blocknet.py", "attention.py", "abm_routes.py",
                  "twin_core.py"}},
        versions={name: importlib.metadata.version(name) for name in ("numpy", "scikit-learn", "torch")})


def support_identity(trait, frame, config, representation, score_receipt, block_receipt, sources):
    return dict(trait=trait,
                input=score_receipt if representation == "csx_scores" else block_receipt["array_binding"],
                ids=ordered_ids_sha256(frame.id), groups=ordered_ids_sha256(frame.group_id),
                roles=ordered_ids_sha256(frame.role), support=config["support"],
                representation=representation,
                source={k: v for k, v in sources.items() if k.endswith(
                    ("support.py", "representation.py", "grid_inputs.py", "pipeline.py"))})


def run_trait(args, config, trait, device):
    from resume import run_lock
    from .comparators import fit_comparators
    input_start = time.monotonic()
    frame, features, csx_columns, metadata = load_grid(args, config, trait)
    log(f"Verified {trait}: test rows={int(frame.role.eq('test').sum()):,}; "
        f"twin_core.id_hash={metadata['test_roster_hash']}; shared side columns={len(features)}")
    if args.smoke:
        frame = smoke_subset(frame, args.smoke_rows, config["head_options"]["seed"])
    needs_blocks = args.predictor == "blocknet" or args.representation == "csx_ldblocks_fixed"
    array = tokens = block_receipt = block_path = None
    rows = np.arange(len(frame))
    if needs_blocks:
        array, cohort, tokens, block_receipt, block_path = locate_blocks(
            config, trait, args.ld_block_files, args.ld_block_dir, allow_build=not args.smoke)
        rows = pd.Index(cohort.eid.astype(str)).get_indexer(frame.id)
        if (rows < 0).any():
            raise ValueError("Block cache lacks frozen roster IDs")
        if args.smoke:
            array = np.asarray(array[rows]).copy()
            rows = np.arange(len(frame))
        log(f"Verified block producer and ID alignment: {trait}; shape={array.shape}")
    sources = source_manifest()
    dependencies = {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "scipy", "scikit-learn", "torch", "pyarrow")}
    block_hash = (block_receipt["array_binding"]["array_sha256"] if needs_blocks else None)
    scores = frame[[features[i] for i in csx_columns]].to_numpy(dtype="float64")
    score_receipt = dict(kind="csx_scores", columns=[features[i] for i in csx_columns],
        values_sha256=hashlib_array(scores), shape=list(scores.shape), dtype=str(scores.dtype),
        axes=["individual", "score"], ordered_ids_sha256=ordered_ids_sha256(frame.id))
    input_seconds = time.monotonic() - input_start
    contract = dict(config=config, trait=trait, candidate=args.predictor, representation=args.representation,
                    prepared_sha256=metadata["input_sha256"], sources=sources, dependencies=dependencies,
                    ordered_ids_sha256=ordered_ids_sha256(frame.id), block_array_sha256=block_hash,
                    score_manifest=score_receipt,
                    group_contract=group_contract(frame, metadata, Path(args.analysis_root) / "local/run_groups" / trait),
                    device=device, smoke=args.smoke)
    version = config["analysis_version"] + ("-SMOKE" if args.smoke else "")
    parent = Path(args.output_root or str(Path(args.analysis_root) / "abm")) / version / trait
    run_id = args.predictor + "_" + digest(contract)[:16]
    if args.new_run:
        from datetime import datetime
        run_id += "_" + datetime.now().strftime("%Y%m%dT%H%M%S")
    run = parent / run_id
    local = run / "_LOCAL_ONLY"
    local.mkdir(parents=True, exist_ok=True)
    with run_lock(run):
        stages = Stages(run, contract)
        if not (run / "runtime.tsv").exists():
            reporting.save_table(pd.DataFrame([dict(phase="input_verification", device="cpu",
                seconds=input_seconds, n=len(frame), rows_per_second=len(frame) / max(input_seconds, 1e-9),
                max_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2)]),
                run / "runtime.tsv")
        stop = PHASES.index(args.through)
        cfg = json.loads(Path(config["block_config"]).read_text()) if needs_blocks else {}
        cfg.update(config.get("predictors", {}), _binary=trait == "t2dm", seed=config["head_options"]["seed"])
        if args.smoke:
            cfg.update(epochs=1, patience=1, d_model=16, heads=2, latent_tokens=4, layers=1,
                       batch_size=16, moe_clusters=[4])
        train_frame = frame.copy()
        train_frame.loc[train_frame.role.isin(["test", "audit", "tune_gate", "head_calibration", "policy_calibration"]), "y"] = np.nan
        roles = role_manifest(frame, trait == "t2dm")
        loss_name = config["primary_loss"][trait]

        def prepare():
            frame.drop(columns="y").to_parquet(local / "features.parquet", index=False)
            frame[["id", "group_id", "role", "y"]].to_parquet(local / "truth.parquet", index=False)
            write_json(run / "role_manifest.json", roles)
            write_json(run / "RUN_MANIFEST.json", dict(contract=contract, metadata=metadata, stages=PHASES,
                       method=f"ABM(base=csx-linear,candidate=csx-{args.predictor})",
                       source_checkout="local_without_git_metadata", source_reference_commits=json.loads(
                           (ROOT / "config/abm.source_bases.json").read_text())))
            write_json(run / "SOURCE_FROZEN.json", sources)
            write_json(run / "input_manifest.json", dict(scores=score_receipt, blocks=block_receipt,
                groups=dict(ordered_groups_sha256=ordered_ids_sha256(frame.group_id),
                    verified=metadata["groups_verified"], evidence=metadata.get("group_verification"),
                    source_sha256=metadata["input_sha256"])))
            for path in sources:
                src = ROOT.parent / path
                target = local / "SOURCE_FROZEN" / path
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, target)
            reporting.save_table(feature_manifest(features, metadata), run / "feature_manifest.tsv")
            write_json(run / "phenotype_provenance.json", dict(
                endpoint=metadata["endpoint_definition"], units="source_native_unverified",
                LDL_treatment_adjustment="source_definition_not_recoded",
                group_status=metadata.get("family_status", metadata["grouping_status"]),
                unverified_fields=["phenotype_units", "LDL_treatment_adjustment",
                                   "GWAS_DOI", "GWAS_N", "UKB_overlap"],
                gwas_overlap="unverified", sources=metadata.get("evolution", {}).get(trait, {}).get("weights_sources", [])))
            return [local / "features.parquet", local / "truth.parquet", run / "role_manifest.json",
                    run / "RUN_MANIFEST.json", run / "SOURCE_FROZEN.json", run / "input_manifest.json"]
        stages.perform("prepare", prepare, "cpu", len(frame))
        if stop == 0:
            return run

        def fit():
            identity = predictor_identity(train_frame, features, csx_columns, args.predictor, cfg, block_receipt, sources)
            fit_cache = Path(args.analysis_root) / "cache/abm/predictors" / digest(identity) / "_LOCAL_ONLY"
            fit_marker = fit_cache / "COMPLETE.json"
            outputs = ["_LOCAL_ONLY/predictor_bundle.joblib", "_LOCAL_ONLY/predictor_predictions.parquet",
                       "comparator_status.tsv", "feature_manifest.tsv", "predictor_lineage.json",
                       "training_history.json", "predictor_identity.json", "block_staging.json"]
            if fit_marker.is_file():
                saved = json.loads(fit_marker.read_text())
                if saved["identity"] != identity:
                    raise ValueError("Predictor dependency identity mismatch")
                for name, expected in saved["files"].items():
                    if file_sha256(fit_cache / name) != expected:
                        raise ValueError("Predictor cache changed: " + name)
                    dest = run / name
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(fit_cache / name, dest)
                stages.execution_details = {"fit_predictors": dict(cache_reused=True,
                    original_cache_device=saved["executed_device"])}
                return [run / name for name in outputs]
            bundle, predictions = fit_predictors(train_frame, features, csx_columns, args.predictor, cfg,
                                                  blocks=array, block_rows=rows, tokens=tokens, device=device)
            external, external_models, statuses = fit_comparators(train_frame, features, csx_columns, cfg,
                                                                  include_moe=args.predictor != "moe")
            for name, values in external.items():
                predictions[name] = values
            bundle["comparators"] = external_models
            joblib.dump(bundle, local / "predictor_bundle.joblib", compress=3)
            predictions.to_parquet(local / "predictor_predictions.parquet", index=False)
            reporting.save_table(statuses, run / "comparator_status.tsv")
            reporting.save_table(feature_manifest(features, metadata, bundle), run / "feature_manifest.tsv")
            write_json(run / "predictor_lineage.json", bundle["lineage"])
            write_json(run / "training_history.json", bundle["histories"])
            write_json(run / "predictor_identity.json", identity)
            write_json(run / "block_staging.json", bundle.get("block_staging", {"strategy":"not_required"}))
            for name in outputs:
                dest = fit_cache / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(run / name, dest)
            write_json(fit_marker, dict(identity=identity, executed_device=device if args.predictor == "blocknet" else "cpu",
                files={name:file_sha256(fit_cache / name) for name in outputs}))
            return [run / name for name in outputs]
        stages.perform("fit_predictors", fit, device if args.predictor == "blocknet" else "cpu", len(frame))
        if stop == 1:
            return run
        support_key = digest(support_identity(trait, frame, config, args.representation,
                                             score_receipt, block_receipt, sources))
        cache = Path(args.analysis_root) / "cache/abm" / ("SMOKE" if args.smoke else "support") / trait / support_key / "_LOCAL_ONLY"

        def support_phase():
            if args.representation == "csx_scores":
                support_array = scores[:, :, None]
                support_rows = np.arange(len(frame))
                maps = np.array(["aggregate_scores"] * 4)
            else:
                support_array, support_rows = array, rows
                maps = block_maps(tokens, array.shape[2])
            support, saved = prepare_support(frame, support_array, support_rows, maps,
                                              score_receipt if args.representation == "csx_scores" else block_receipt,
                                              cache, config, device, metadata["groups_verified"])
            support.to_parquet(local / "support.parquet", index=False)
            write_json(run / "geometry_manifest.json", saved["geometry"])
            write_json(run / "donor_bank_manifest.json", dict(saved["bank"],
                       group_evidence=metadata.get("group_verification")))
            stages.execution_details = {"build_support": dict(cache_reused=saved.get("cache_reused", False),
                original_cache_device=saved.get("executed_device", "unrecorded"))}
            write_json(local / "support_cache.json", {"path": str(cache), "receipt_sha256": file_sha256(cache / "COMPLETE.json")})
            return [local / "support.parquet", local / "support_cache.json",
                    run / "geometry_manifest.json", run / "donor_bank_manifest.json"]
        stages.perform("build_support", support_phase, device, len(frame))
        if stop == 2:
            return run
        predictions = pd.read_parquet(local / "predictor_predictions.parquet")
        support = pd.read_parquet(local / "support.parquet")
        joined = frame.drop(columns="y").merge(predictions.drop(columns=["group_id", "role"]),
                 on="id", validate="one_to_one").merge(support, on="id", validate="one_to_one")
        dev = joined.merge(frame[["id", "y"]], on="id", validate="one_to_one")
        manifest = dict(baseline_id="csx-linear", candidate_id="csx-" + args.predictor,
                        bank_sha256=file_sha256(run / "donor_bank_manifest.json"),
                        geometry_sha256=file_sha256(run / "geometry_manifest.json"),
                        predictor_bundle_sha256=file_sha256(local / "predictor_bundle.joblib"),
                        role_manifest_sha256=file_sha256(run / "role_manifest.json"),
                        trait=trait, loss_name=loss_name, grouping_status=metadata["grouping_status"],
                        validation_status=metadata["validation_status"], gwas_overlap_status="unverified")

        def fit_assessor():
            gate = dev.loc[dev.role.eq("tune_gate")].copy()
            head = dev.loc[dev.role.eq("head_calibration")].copy()
            policy = joined.loc[joined.role.eq("policy_calibration")].copy()
            validate_lineage(frame)
            model = ABMAssessor.fit(gate, head, policy, manifest=manifest,
                                    task="binary" if trait == "t2dm" else "continuous",
                                    loss_name=loss_name,
                                    head_options=config["head_options"], min_gain=config["min_gain"])
            model.save(local / "assessor.json")
            reporting.save_table(model.heads.calibration_table(), run / "count_calibration.tsv")
            clinical = [c for c in features if not c.startswith("csx.")]
            diagnostics = reporting.fit_diagnostics(gate, head, clinical, model.heads.task, policy, loss_name=loss_name)
            joblib.dump(diagnostics, local / "diagnostics.joblib")
            write_json(local / "pc_edges.json", {c: np.unique(np.quantile(gate[c], [.2, .4, .6, .8])).tolist()
                       for c in features if c.startswith("match_PC")})
            return [local / "assessor.json", run / "count_calibration.tsv", local / "diagnostics.joblib", local / "pc_edges.json"]
        stages.perform("fit_assessor", fit_assessor, "cpu", len(frame))
        if stop == 3:
            return run
        model = ABMAssessor.load(local / "assessor.json")

        def freeze():
            files = []
            for label in ("test", "audit", "policy_calibration"):
                folder = local / label
                folder.mkdir(exist_ok=True)
                query = joined.loc[joined.role.eq(label)].copy().reset_index(drop=True)
                assessed = model.assess(query, manifest=manifest)
                assessed = reporting.diagnostic_predictions(joblib.load(local / "diagnostics.joblib"), assessed)
                for k, v in (("trait", trait), ("baseline_id", manifest["baseline_id"]), ("candidate_id", manifest["candidate_id"]),
                              ("bank_id", manifest["bank_sha256"]), ("geometry_id", manifest["geometry_sha256"])):
                    assessed[k] = v
                assessed.to_parquet(folder / "predictions.parquet", index=False)
                write_feature_manifest(folder / "predictions.parquet", manifest=manifest)
                reporting.freeze_matches(assessed, folder, model.heads.seed)
                write_json(run / (label + "_dataset_assessment.json"), model.dataset_assessment(assessed))
                write_json(folder / "FROZEN.json", dict(predictions_sha256=file_sha256(folder / "predictions.parquet"),
                           model_sha256=file_sha256(local / "assessor.json"), loss_name=loss_name,
                           truth_sha256=file_sha256(local / "truth.parquet"), query_outcomes_used=False))
                files += list(folder.iterdir())
            return files
        stages.perform("freeze_predictions", freeze, "cpu", len(frame))
        if stop == 4:
            return run

        def evaluate_phase():
            files = []
            for label in ("test", "audit", "policy_calibration"):
                folder = local / label
                frozen = json.loads((folder / "FROZEN.json").read_text())
                if frozen["model_sha256"] != file_sha256(local / "assessor.json") or frozen["predictions_sha256"] != file_sha256(folder / "predictions.parquet"):
                    raise ValueError("Prediction/bundle binding failed before joining truth")
                assessed = pd.read_parquet(folder / "predictions.parquet")
                evaluated, tables = reporting.evaluate(assessed, frame.loc[frame.role.eq(label), ["id", "y"]],
                                                       model, run, config, label=label)
                if label == "test":
                    reporting.save_table(reporting.support_diagnostics(evaluated, trait == "t2dm",
                        json.loads((local / "pc_edges.json").read_text())), run / "support_by_group.tsv")
                    reporting.save_table(reporting.evaluate_diagnostics(joblib.load(local / "diagnostics.joblib"), evaluated, config),
                                         run / "count_incremental_info.tsv")
                    contrast = reporting.count_incremental_contrasts(
                        joblib.load(local / "diagnostics.joblib"), evaluated, config)
                    for key, value in dict(trait=trait, run_id=run_id, baseline_id=manifest["baseline_id"],
                                           candidate_id=manifest["candidate_id"], representation=args.representation).items():
                        contrast[key] = value
                    reporting.save_table(contrast, run / "count_incremental_contrasts.tsv")
                    perf = tables["policy_performance"]
                    from .evaluation import prediction_metrics
                    for name in ("Clinical", "csx-moe"):
                        if name in evaluated:
                            perf = pd.concat([perf, pd.DataFrame([dict(method=name, trait=trait, stage="test",
                                **prediction_metrics(evaluated.y.to_numpy(), evaluated[name].to_numpy(), model.heads.task))])], ignore_index=True)
                    for key, value in dict(score_source="csx", candidate=args.predictor,
                        representation=args.representation, run_id=run_id, assessment_loss=loss_name,
                        predictor_bundle_sha256=manifest["predictor_bundle_sha256"],
                        feature_manifest_sha256=file_sha256(run / "feature_manifest.tsv")).items():
                        perf[key] = value
                    reporting.save_table(perf, run / "main_performance.tsv")
                    flow = []
                    for role, part in dev.groupby("role"):
                        flow.append(dict(role=role, n=len(part), observed=part.y.notna().sum(), supported=part.supported.sum(),
                                         groups=part.group_id.nunique(), cases=part.y.sum() if trait == "t2dm" else None))
                    flow.append(dict(role="test_assessed", n=len(assessed), observed=evaluated.y.notna().sum(),
                                     supported=assessed.supported.sum(), assessable=assessed.assessment_available.sum(),
                                     selected=assessed.selected.sum(), switched=assessed.switched.sum()))
                    reporting.save_table(pd.DataFrame(flow), run / "cohort_flow.tsv")
            files = [p for p in run.glob("*.tsv") if p.name != "runtime.tsv"]
            files += list(run.glob("test_*.json")) + list(run.glob("audit_*.json")) + list(run.glob("policy_calibration_*.json"))
            return files
        stages.perform("evaluate", evaluate_phase, "cpu", len(frame))
        if stop == 5:
            return run

        def report_phase():
            from final_compare import workbook
            perf = pd.read_csv(run / "main_performance.tsv", sep="\t")
            workbook(perf, run / "main_performance.xlsx")
            calibration_plot = run / "test_count_calibration.png"
            reporting.plot_calibration(pd.read_csv(run / "test_count_calibration.tsv", sep="\t"), calibration_plot)
            table_manifest = run / "table_manifest.json"
            write_json(table_manifest, [dict(name=p.stem, path=p.name, sha256=file_sha256(p),
                trait=trait, module="abm", run_id=run_id)
                for p in sorted(run.glob("*.tsv")) if p.name != "runtime.tsv"])
            write_json(run / "COMPLETE.json", dict(status="completed", scientific_claim="exploratory",
                       smoke=args.smoke, manifest_sha256=file_sha256(run / "RUN_MANIFEST.json"),
                       model_sha256=file_sha256(local / "assessor.json")))
            return [run / "main_performance.xlsx", run / "COMPLETE.json", calibration_plot, table_manifest]
        stages.perform("report", report_phase, "cpu", len(frame))
        from kinship import active, compatible
        if not args.smoke and compatible(metadata, active(args.analysis_root)):
            write_json(parent / (args.predictor + ".CURRENT.json"), dict(run=run_id))
        log("COMPLETE " + str(run))
        return run


def hashlib_array(array):
    import hashlib
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def parser():
    p = argparse.ArgumentParser(description="ABM: fixed genetic support, expected error and absolute benefit")
    p.add_argument("--traits", default="height,ldl,t2dm")
    p.add_argument("--source", choices=["csx"], default="csx")
    p.add_argument("--predictor", choices=["blocknet", "moe"], default="blocknet")
    p.add_argument("--representation", choices=["csx_ldblocks_fixed", "csx_scores"], default="csx_ldblocks_fixed")
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--config", default=str(ROOT / "config/abm.json"))
    p.add_argument("--analysis-root", default="/mnt/d/analysis/grid")
    p.add_argument("--output-root")
    p.add_argument("--prepared")
    p.add_argument("--ld-block-files")
    p.add_argument("--ld-block-dir")
    p.add_argument("--through", choices=PHASES, default="report", help="Run/resume through this phase")
    p.add_argument("--new-run", action="store_true")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--smoke", action="store_true", help="Small real-data test in isolated SMOKE outputs")
    p.add_argument("--smoke-rows", type=int, default=768, help="Maximum per role for smoke only")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--module-only", action="store_true", help="Accepted; ABM always runs independently")
    p.add_argument("--report-from", help="Verify a completed frozen run and derive reports without training")
    return p


def configuration(args):
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    path = Path(args.config).resolve()
    if path.is_dir():
        path /= "abm.json"
    args.config = str(path)
    config = json.loads(path.read_text())
    if config["coverages"] != list(COVERAGES) or config["primary_coverage"] != .5:
        raise ValueError("ABM freezes six coverages and primary coverage 0.50")
    if config["calibration_fit_partition_by_group"] != dict(prediction_calibration=.5, head_calibration=.25, policy_calibration=.25):
        raise ValueError("A new calibration partition requires a new reviewed protocol")
    if config["residual_correction"] or not config["count_only"]:
        raise ValueError("Count-only ABM does not enable residual routes")
    if config["evaluation"].get("matching") != "frozen_risk_bin_disjoint_group_exclusion":
        raise ValueError("ABM matches controls by frozen baseline risk, without ancestry restrictions; update evaluation.matching")
    if any(config["support"][k] for k in ("include_pcs", "include_clinical", "outcome_residual_pruning")):
        raise ValueError("The main geometry excludes PCs, clinical inputs and outcome pruning")
    config["candidate_predictor"] = "csx-" + args.predictor
    config["support"]["representation"] = args.representation
    config["execution"]["device"] = args.device
    if args.predictor == "blocknet" and args.representation == "csx_scores":
        log("Sensitivity: support in four-score space; blocknet still uses LD blocks")
    traits = [trait.strip() for trait in args.traits.split(",")]
    if not traits or len(traits) != len(set(traits)) or not set(traits) <= {"height", "ldl", "t2dm"}:
        raise ValueError("Traits must be a nonempty, unique list of height,ldl,t2dm")
    args.traits = ",".join(traits)
    for trait in traits:
        if not isinstance(config["primary_loss"].get(trait), str):
            raise ValueError("An explicit assessment loss is required for " + trait)
        row_loss(np.array([], float), np.array([], float),
                 "binary" if trait == "t2dm" else "continuous", config["primary_loss"][trait])
    if args.smoke:
        config["head_options"]["n_boot"] = 30
        config["evaluation"]["bootstrap_replicates"] = 30
    return config


def main(argv=None):
    args = parser().parse_args(argv)
    config = configuration(args)
    if args.report_from:
        from .derived_report import report_from
        if args.dry_run:
            print(json.dumps(dict(modules=["abm"], phases=["derived_report"],
                                 parent=args.report_from, output=args.output_root, training=False)))
            return
        return report_from(args.report_from, args.output_root)
    traits = args.traits.split(",")
    if args.dry_run:
        print(json.dumps(dict(modules=["abm"], phases=PHASES[:PHASES.index(args.through) + 1], traits=traits,
                              predictor=args.predictor, representation=args.representation,
                              device=args.device, config_path=str(Path(args.config).resolve()),
                              after_traits=[],
                              config=config), indent=2))
        return
    log(f"Traits in order: {' -> '.join(traits)}; predictor={args.predictor}; "
        f"representation={args.representation}; device={args.device}; config={args.config}")
    device, _ = _backend(args.device)
    from threadpoolctl import threadpool_limits
    import torch
    torch.set_num_threads(args.threads)
    from kinship import analysis_guard
    with analysis_guard(args.analysis_root), threadpool_limits(limits=args.threads):
        os.environ.setdefault("GRID_CACHE_ROOT", str(Path(args.analysis_root) / "cache"))
        from .cold_start import ensure_inputs
        from .storage import clean_completed_support
        # This also reaches historical keys and runs whose build_support phase
        # is already complete, without rehashing caches with no scratch left.
        clean_completed_support(args.analysis_root, traits=traits, apply=True, best_effort=True)
        ensure_inputs(args, config)
        for i, trait in enumerate(traits, start=1):
            log(f"TRAIT {i}/{len(traits)} {trait}; stages through {args.through}")
            try:
                run_trait(args, config, trait, device)
            finally:
                # Retry anything that was busy at startup, even if a later
                # fitting/reporting phase failed or all phases were reused.
                clean_completed_support(args.analysis_root, traits=[trait], apply=True, best_effort=True)
                if args.predictor == "blocknet" or args.representation == "csx_ldblocks_fixed":
                    from dosage_storage import release_dosages
                    release_dosages(config["block_config"], traits=config["traits"],
                                    analysis_root=args.analysis_root, apply=True, best_effort=True)


if __name__ == "__main__":
    main()
