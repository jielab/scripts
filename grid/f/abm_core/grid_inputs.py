"""Real GRID input and role contracts; all joins are by string identity."""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .contracts import ordered_ids_sha256, verify_array_binding, binding_metadata
from .model import file_sha256, write_json
from .heads import hash_fraction

ROLES = ("build", "tune_model", "tune_gate", "prediction_calibration",
         "head_calibration", "policy_calibration", "audit", "test")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str,
                                     separators=(",", ":")).encode()).hexdigest()


def table(path):
    return pd.read_csv(path, sep="\t", dtype={"id": str, "eid": str, "group_id": str,
                                             "group": str, "family_id": str})


def partition_calibration(frame, seed):
    d = frame.copy()
    original = d.role.to_numpy(str)
    u = hash_fraction(d.group_id, seed)
    sub = np.where(u < .5, "prediction_calibration",
                   np.where(u < .75, "head_calibration", "policy_calibration"))
    d["parent_role"] = original
    d["role"] = np.where(original == "calibration_fit", sub, original)
    if d.groupby("group_id").role.nunique().max() != 1:
        raise ValueError("Family/group crosses development or test roles")
    return d


def role_manifest(frame, binary):
    records = {}
    for role, d in frame.groupby("role", sort=True):
        observed = d.y.notna()
        records[role] = dict(n=len(d), n_observed=int(observed.sum()),
                             cases=int(d.y.sum()) if binary else None,
                             groups=d.group_id.nunique(),
                             ids_sha256=ordered_ids_sha256(d.id),
                             groups_sha256=ordered_ids_sha256(d.group_id))
    return records


def validate_lineage(frame):
    if frame.id.isna().any() or frame.id.duplicated().any() or frame.group_id.isna().any():
        raise ValueError("Unique string IDs and explicit groups are required")
    if frame.groupby("group_id").role.nunique().max() != 1:
        raise ValueError("Upstream predictor/head roles share a group")
    if set(frame.role) - set(ROLES):
        raise ValueError("Unknown ABM roles")
    for role in ROLES:
        if not frame.role.eq(role).any():
            raise ValueError(f"Unavailable role {role}; no samples may be borrowed from test")


def load_grid(args, config, trait):
    import twin
    from twin_core import assign_roles, id_hash
    prepared = args.prepared or str(Path(args.analysis_root) / "cache/abm/prepared.joblib")
    if not Path(prepared).is_file():
        if args.smoke:
            raise FileNotFoundError("SMOKE requires existing prepared input; no full preparation: " + prepared)
        from types import SimpleNamespace
        from abm_prepare import prepare_abm
        prepare_abm(SimpleNamespace(analysis_root=args.analysis_root,
                    cache_dir=str(Path(args.analysis_root) / "cache/abm"),
                    prepared=prepared, traits=args.traits, out_root=str(Path(args.analysis_root) / "abm"),
                    prepare_config=None))
    d, groups, metadata = twin.load_source(prepared, trait)
    outer_test = d.split.eq("test")
    roster_hash = id_hash(d.loc[outer_test, "eid"])
    roster_files = [Path(args.analysis_root) / name for name in (
        "prsformer/cache/data.tsv.gz", "prsformer/split.rds",
        "prsformer/cache/roster/prsformer.split.tsv.gz", "prsformer/cache/roster/families.tsv.gz")]
    if any(p.exists() for p in roster_files):
        twin.verify_existing_prsformer_roster(args.analysis_root, d, d.loc[outer_test, "eid"])
        metadata["roster_verification"] = "prepared_cohort_and_shared_roster"
    else:
        # load_source already checks every ID, split and family against the
        # producer's original cohort. No fitted benchmark is needed.
        metadata["roster_verification"] = "prepared_original_cohort"
    features = list(dict.fromkeys(groups["csx"] + groups["ancestry"] + groups["covariates"]))
    if len(features) != 18 or len(groups["csx"]) != 4 or len(groups["ancestry"]) != 10:
        raise ValueError("Expected shared 4 CSx + 10 match_PC + 4 covariates contract")
    d = d.rename(columns={"eid": "id", "family_id": "group_id"}).copy()
    if "group_id" not in d:
        d["group_id"] = d.id
    d["group_id"] = d.group_id.astype(str)
    d["id"] = d.id.astype(str)
    cfg = twin.config_file(str(Path(__file__).resolve().parents[2] / "config/twin.json"), trait)
    d["role"] = assign_roles(d.group_id, cfg)
    d.loc[outer_test, "role"] = "test"
    # Development rows without Y are not supervised donors. All test rows survive.
    d = d.loc[d.y.notna() | d.role.eq("test")].reset_index(drop=True)
    d = partition_calibration(d, int(config["head_options"]["seed"]))
    validate_lineage(d)
    group_status = metadata.get("grouping_status", cfg.grouping_status)
    # Config cannot turn assumed singleton identities into verified kinship.
    verified = False
    group_evidence = metadata.get("group_verification")
    evidence_status = "not_provided"
    if isinstance(group_evidence, dict) and group_evidence.get("status") == "verified":
        source = Path(group_evidence.get("source", ""))
        if source.is_file() and file_sha256(source) == group_evidence.get("sha256"):
            supplied = twin.identified(twin.read_table(source))
            group_column = next((c for c in ("group_id", "family_id", "group") if c in supplied), None)
            if group_column and set(d.id) <= set(supplied.eid):
                aligned = supplied.set_index("eid").loc[d.id, group_column]
                if aligned.isna().any() or not np.array_equal(aligned.astype(str), d.group_id):
                    raise ValueError("Verified group source differs from the prepared ID/group mapping")
                verified, evidence_status = True, "source_hash_and_ID_group_mapping_verified"
            else:
                evidence_status = "source_lacks_complete_ID_group_mapping"
        else:
            evidence_status = "source_missing_or_checksum_mismatch"
    metadata.update(test_roster_hash=roster_hash, test_roster_n=int(outer_test.sum()),
                    test_roster_hash_algorithm="twin_core.id_hash",
                    groups_verified=verified,
                    grouping_status="verified_from_source" if verified else "unverified_group_identity",
                    grouping_status_declared=group_status,
                    group_evidence_status=evidence_status,
                    endpoint_definition="baseline_t2dm_from_Yr2e_Yt2e" if trait == "t2dm" else trait,
                    gwas_overlap_status="unverified", validation_status="exploratory_previously_inspected_test")
    d = d[["id", "group_id", "role", "parent_role", "y", "ancestry", *features]].copy()
    return d, features, [features.index(c) for c in groups["csx"]], metadata


def smoke_subset(frame, per_role, seed):
    """Testing-only subcohort, never exported as the production frozen roster."""
    parts = []
    for role, d in frame.groupby("role", sort=False):
        rank = hash_fraction(d.id, seed)
        parts.append(d.iloc[np.argsort(rank, kind="stable")[:per_role]])
    return pd.concat(parts, ignore_index=True)


def locate_blocks(config, trait, block_files=None, block_dir=None, *, allow_build=True):
    import block_inputs as b
    base_cfg = b.read_config(config["block_config"])
    layout = b.resolve_block_files(block_dir or base_cfg["ld_block_dir"],
                                   block_files or base_cfg["ld_block_files"])
    paths = b.cache_paths(base_cfg, layout, trait)
    if not paths[-1].is_file():
        # Existing genuinely bound outputs can be verified in place. Unbound
        # outputs are preserved; production rebuilds only this trait/layout.
        prior = Path(base_cfg["cache_root"]) / layout["layout_id"] / trait
        marker = prior / "BLOCKS_COMPLETE.json"
        if marker.is_file():
            previous = json.loads(marker.read_text())
            if previous.get("array_binding") and previous.get("producer_version") == b.PRODUCER_VERSION:
                paths = (prior, prior/"block_scores.npy", prior/"tokens.tsv.gz", prior/"cohort.tsv.gz", marker)
            elif allow_build:
                print("[GRID ABM] Unbound blocks preserved at " + str(prior) +
                      "; rebuilding only " + trait + "/" + layout["layout_id"], flush=True)
    if not paths[-1].is_file():
        if not allow_build:
            raise FileNotFoundError("SMOKE requires producer-bound blocks; no full block build: " + str(paths[-1]))
        # Prepares genotypes/CSx blocks only, never requires any fitted CPU model.
        b.prepare_trait(config["block_config"], trait, block_dir, block_files, device=config["execution"]["device"])
    receipt = json.loads(paths[-1].read_text())
    if receipt.get("status") != "complete" or receipt["trait"] != trait:
        raise ValueError("Incomplete or wrong-trait block cache")
    payload = receipt["signature_payload"]
    if receipt.get("producer_version") != b.PRODUCER_VERSION or "array_binding" not in receipt:
        raise ValueError("Missing trustworthy block producer binding; rebuild this trait/layout from genotype and weights")
    # The original producer wrote both arrays and this cohort, in this order.
    # Verify that order against the exact input whose checksum it committed.
    original = Path(base_cfg["prsformer_cache"]) / "data.tsv.gz"
    if file_sha256(original) != payload["cohort_sha256"]:
        raise ValueError("Block producer input roster changed")
    cohort = table(paths[3])
    source = table(original)
    if not cohort.eid.astype(str).equals(source.eid.astype(str)):
        raise ValueError("Block row identities differ from the original producer")
    tokens = table(paths[2])
    expected, channels, _ = b.build_manifest(layout, base_cfg)
    verify_array_binding(paths[1], binding_metadata(cohort),
                         receipt["array_binding"], tokens=tokens, channels=channels)
    pd.testing.assert_frame_equal(tokens.reset_index(drop=True), expected.reset_index(drop=True),
                                  check_dtype=False, check_exact=False, rtol=1e-14)
    if receipt["layout"]["signature"] != layout["signature"]:
        raise ValueError("LD layout mismatch")
    for pop, path in b.locate_weights(base_cfg, trait)[0].items():
        if file_sha256(path) != payload["weights"][pop]:
            raise ValueError(f"CSx posterior changed: {trait}/{pop}")
    if file_sha256(Path(base_cfg["score_root"]) / trait / "1csx.scores.rds") != payload["whole_score_sha256"]:
        raise ValueError("Full CSx scores changed")
    arr = np.load(paths[1], mmap_mode="r")
    if list(arr.shape) != receipt["shape"] or len(arr) != len(cohort):
        raise ValueError("Block array dimensions differ from producer receipt")
    for x in receipt["reconstruction"]:
        if (x["assigned_reconstruction_correlation"] < base_cfg["reconstruction_corr_min"] or
            not base_cfg["reconstruction_slope_min"] <= x["assigned_reconstruction_slope"] <= base_cfg["reconstruction_slope_max"] or
            x["assigned_reconstruction_nrmse"] > base_cfg["reconstruction_nrmse_max"]):
            raise ValueError("Block reconstruction failed")
    receipt = dict(receipt, adapter_probe=reconstruction_probe(
        arr, tokens, cohort, base_cfg, receipt, trait, b))
    return arr, cohort, tokens, receipt, paths[1]


def reconstruction_probe(array, tokens, cohort, config, receipt, trait, b):
    from dosage_storage import dosage_lock
    with dosage_lock(config, exclusive=True):
        return _reconstruction_probe(array, tokens, cohort, config, receipt, trait, b)


def _reconstruction_probe(array, tokens, cohort, config, receipt, trait, b):
    """Independent dosage scoring; retained real rows preserve audits after cleanup."""
    from block_dosage import SupplementedDosage
    from dosage_storage import probe_rows, save_probe, load_probe
    cache = Path(config["prsformer_cache"])
    variants = table(cache / "variants.tsv.gz")
    if file_sha256(cache / "variants.tsv.gz") != receipt["signature_payload"]["variants"]:
        raise ValueError("Original genotype variant manifest changed")
    if file_sha256(cache / "prepare.json") != receipt["signature_payload"]["prsformer_prepare"]:
        raise ValueError("Original genotype preparation identity changed")
    rr = probe_rows(len(cohort))
    # The actual local block path is supplied by locate_blocks, including legacy
    # bound layouts. It is never taken from an unverified receipt path.
    base = Path(array.filename).parent
    supplement = receipt.get("dosage_supplement")
    supplement_path = None
    if supplement:
        folder = Path(config["cache_root"]) / "dosage_supplements" / supplement["signature"]
        if json.loads((folder / "COMPLETE.json").read_text()) != supplement:
            raise ValueError("Supplement producer receipt changed")
        supplement_path = folder / "genotypes.npy"
        variants = pd.concat([variants, pd.DataFrame(supplement["identity"]["variants"])], ignore_index=True)
    mappings = []
    try:
        if ((cache / "genotypes.npy").is_file() and
                b.sha256_cached(cache / "genotypes.npy") != receipt["signature_payload"]["genotype_sha256"]):
            raise ValueError("Original genotype content changed; rebuild affected blocks")
        if (supplement_path is not None and supplement_path.is_file() and
                b.sha256_cached(supplement_path) != supplement["genotype_sha256"]):
            raise ValueError("Supplement genotype content changed")
        full = (cache / "genotypes.npy").is_file() and (supplement_path is None or supplement_path.is_file())
        if full:
            genotype = np.load(cache / "genotypes.npy", mmap_mode="r", allow_pickle=False)
            mappings.append(genotype)
            if len(genotype) != len(cohort):
                raise ValueError("Original dosage row count differs")
            if supplement_path is not None:
                added = np.load(supplement_path, mmap_mode="r", allow_pickle=False)
                mappings.append(added)
                genotype = SupplementedDosage(genotype, added)
            direct_rows = rr
        else:
            genotype = load_probe(base, receipt, rr)
            direct_rows = np.arange(len(rr))
        if genotype.shape[1] != len(variants):
            raise ValueError("Audit dosage variant count differs")
        checks = []
        for channel, pop in enumerate(b.POPS):
            weights = b.read_weights(b.locate_weights(config, trait)[0][pop])
            if receipt["layout"]["mode"] == "population_multiview":
                bed = tokens[tokens.source_population.eq(pop)]
                block_score = np.asarray(array[rr][:, tokens.source_population.eq(pop), 0], float).sum(1)
            else:
                bed = tokens
                block_score = np.asarray(array[rr][:, :, channel], float).sum(1)
            direct = b.direct_assigned_score(genotype, direct_rows, weights, bed, variants)
            checks.append(dict(population=pop, **b.assigned_reconstruction(block_score, direct, config)))
        if full:
            save_probe(base, receipt, rr, genotype)
        return dict(n=len(rr), row_ids_sha256=ordered_ids_sha256(cohort.iloc[rr].eid),
                    method="independent_assigned_dosage_sum_vs_cached_blocks", checks=checks)
    finally:
        for value in mappings:
            value._mmap.close()


def block_maps(tokens, channels):
    if channels == 1:
        return tokens.source_population.to_numpy(str)
    if channels == 4 and set(tokens.source_population) == {"ALL"}:
        return np.tile(["AFR", "EAS", "EUR", "SAS"], len(tokens))
    raise ValueError("Unknown block/channel layout; cannot flatten without map identities")


class IndexedBlocks:
    """Bounded row access; no whole-genome advanced-index copy in RAM."""
    def __init__(self, array, rows, flatten=True):
        self.array, self.rows, self.flatten = array, np.asarray(rows, int), flatten
        self.shape = (len(self.rows), int(np.prod(array.shape[1:]))) if flatten else (len(self.rows), *array.shape[1:])

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, item):
        v = np.asarray(self.array[self.rows[item]])
        return v.reshape(len(v), -1) if self.flatten else v


def feature_manifest(features, metadata, bundle=None):
    settings = metadata.get("prepared_settings", {})
    records = []
    for i, name in enumerate(features):
        source = (str(settings.get("score_dir")) if name.startswith("csx.") else
                  str(settings.get("pca_file")) if name.startswith("match_PC") else str(settings.get("pheno_file")))
        records.append(dict(feature=name, order=i, source=source, unit="source_native",
                            use_linear=True, use_moe=True, use_blocknet=True, use_twin_distance=False,
                            scaler_fit_role="build",
                            mean=bundle["scale"].mean_[i] if bundle else None,
                            scale=bundle["scale"].scale_[i] if bundle else None))
    return pd.DataFrame(records)
