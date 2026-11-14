#!/usr/bin/env python3
"""Bound row identities, immutable outer test roster, training-only roles."""
from __future__ import annotations

# Inputs

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def file_record(path):
    from file_state import file_record as record
    return record(path)


def ids_hash(ids):
    return hashlib.sha256(json.dumps(list(map(str, ids)), ensure_ascii=False,
                                    separators=(",", ":")).encode()).hexdigest()


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype={name: str for name in
        ("id", "eid", "group_id", "family_id", "group", "feature_id", "map_id")})


def json_write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False,
                                    default=lambda x: str(x) if isinstance(x, Path) else x.item())+"\n")
    temporary.replace(path)


def normalized_frame(frame, *, check_groups=True):
    d = frame.copy()
    if "id" not in d and "eid" in d:
        d = d.rename(columns={"eid": "id"})
    if "group_id" not in d:
        column = next((c for c in ("family_id", "group") if c in d), None)
        if column:
            d = d.rename(columns={column: "group_id"})
        else:
            d["group_id"] = d["id"]
    for name in ("id", "group_id"):
        if name not in d or d[name].isna().any() or d[name].astype(str).str.strip().eq("").any():
            raise ValueError("Explicit nonmissing identities required: "+name)
        d[name] = d[name].astype(str)
    if d.id.duplicated().any():
        raise ValueError("Duplicate IDs in source cohort")
    if not {"y", "split"} <= set(d):
        raise ValueError("Cohort requires frozen split and y columns")
    d["split"] = d.split.astype(str).replace({"training": "train", "testing": "test"})
    if not set(d.split) <= {"train", "validation", "test"} or not {"train", "test"} <= set(d.split):
        raise ValueError("Use the existing train/validation/test outer roster")
    if check_groups and d.groupby("group_id").split.nunique().max() != 1:
        raise ValueError("A declared group crosses outer source splits")
    d["y"] = pd.to_numeric(d.y, errors="raise")
    if np.isinf(d.y.to_numpy(float)).any():
        raise ValueError("Infinite outcomes are invalid")
    if "ancestry" not in d:
        d["ancestry"] = "UNASSIGNED"
    d["ancestry"] = d.ancestry.fillna("UNASSIGNED").astype(str)
    return d


def assign_roles(frame, *, donor_fraction=.8, seed=20261102):
    if not 0 < donor_fraction < 1:
        raise ValueError("donor_fraction must be in (0,1)")
    d = normalized_frame(frame)
    if d.split.eq("validation").any():
        d["role"] = d.split.map({"train": "donor", "validation": "calibration", "test": "test"})
    else:
        fractions = {g: int(hashlib.sha256(f"count-abm:{seed}:{g}".encode()).hexdigest()[:16], 16)/2**64
                     for g in d.group_id.unique()}
        d["role"] = np.where(d.split.eq("test"), "test",
                              np.where(d.group_id.map(fractions) < donor_fraction, "donor", "calibration"))
    if d.groupby("group_id").role.nunique().max() != 1:
        raise ValueError("A group crosses donor/calibration/test roles")
    for role in ("donor", "calibration", "test"):
        if not d.role.eq(role).any():
            raise ValueError("No "+role+" rows; no test rows may be borrowed")
    return d


def write_input_manifest(table, array, features, output, *, trait, task, provenance=None):
    """Bind an explicitly exported table/array at production, for tests or adapters.

    The binding checks ordered person/group IDs and array/feature dimensions.
    File records are not content-integrity hashes. It does not verify GWAS
    independence or kinship from a caller's declaration.
    """
    d = normalized_frame(read_table(table))
    f = read_table(features)
    x = np.load(array, mmap_mode="r", allow_pickle=False)
    _validate_features(x, d, f)
    binding = dict(table_file_record=file_record(table), array_file_record=file_record(array),
                   features_file_record=file_record(features), ordered_ids_sha256=ids_hash(d.id),
                   ordered_groups_sha256=ids_hash(d.group_id))
    payload = dict(schema_version="grid.abm.input.1", trait=trait, task=task,
                   table=str(Path(table).resolve()), array=str(Path(array).resolve()),
                   features=str(Path(features).resolve()), binding=binding,
                   provenance=provenance or {})
    json_write(output, payload)
    return Path(output)


def _validate_features(x, d, features):
    if x.ndim != 2 or len(x) != len(d) or x.shape[1] != len(features) or not len(features):
        raise ValueError("Input array must be N x features in exactly the bound row/feature order")
    if not {"feature_id", "map_id"} <= set(features):
        raise ValueError("Feature manifest requires feature_id and map_id")
    if features[["feature_id", "map_id"]].isna().any().any() or features.feature_id.duplicated().any():
        raise ValueError("Feature/map identities must be nonmissing; feature IDs must be unique")


def load_export(path, trait):
    path = Path(str(path).replace("{trait}", trait))
    manifest = json.loads(path.read_text())
    if manifest.get("schema_version") != "grid.abm.input.1" or manifest.get("trait") != trait:
        raise ValueError("Wrong input manifest schema or trait")
    resolved = {key: (path.parent / manifest[key]).resolve() for key in ("table", "array", "features")}
    for key, value in resolved.items():
        if file_record(value) != manifest["binding"][key+"_file_record"]:
            raise ValueError("Input binding changed: "+key)
    d = normalized_frame(read_table(resolved["table"]))
    if ids_hash(d.id) != manifest["binding"]["ordered_ids_sha256"] or ids_hash(d.group_id) != manifest["binding"]["ordered_groups_sha256"]:
        raise ValueError("Input identities do not match the producer binding")
    features = read_table(resolved["features"])
    x = np.load(resolved["array"], mmap_mode="r", allow_pickle=False)
    _validate_features(x, d, features)
    provenance = manifest.get("provenance", {})
    safe = {key: provenance[key] for key in ("gwas_overlap_status", "endpoint_definition", "source_description") if key in provenance}
    metadata = dict(safe, input_manifest=str(path.resolve()),
                    input_manifest_file_record=file_record(path), source_kind="explicit_bound_export",
                    groups_verified=False, grouping_status="provided_identity_unverified")
    return d, x, np.arange(len(d)), features, manifest["task"], metadata


def load_production(args, config, trait):
    """Use raw prepared features and verified CSx blocks, never model predictions."""
    import twin
    from abm.grid_inputs import locate_blocks, block_maps
    from kinship import correct_grid, publication_metadata
    from types import SimpleNamespace
    def finish_features():
        from shared_cohort import ensure_inputs, complete_score_features
        from grid_paths import read_config as source_config
        options = SimpleNamespace(**dict(vars(args), traits=trait, complete_features=True))
        raw = Path(source_config(config["block_config"], args.analysis_root)["prsformer_cache"])
        if (raw/"data.tsv.gz").is_file() and (raw/"prepare.json").is_file():
            complete_score_features(options, config)
        else:
            ensure_inputs(options, config)
    def current_prepared():
        if args.prepared:
            return Path(args.prepared)
        from shared_cohort import prepared_for_trait
        return Path(prepared_for_trait(args.analysis_root, trait))
    prepared = current_prepared()
    if not prepared.is_file():
        finish_features()
        prepared = current_prepared()
    try:
        frame, _, metadata = twin.load_source(prepared, trait)
    except ValueError as exc:
        # A shared cache can be prepared trait by trait. Complete the requested
        # feature table on the same roster; do not silently replace a user-supplied file.
        completable = ("Trait absent" in str(exc) or "not completed genetic feature" in str(exc)
                       or isinstance(exc, twin.PreparedSourceChanged))
        if args.prepared or not completable:
            raise
        finish_features()
        prepared = current_prepared()
        frame, _, metadata = twin.load_source(prepared, trait)
    frame = normalized_frame(frame, check_groups=False)
    original_test = frame.loc[frame.split.eq("test"), "id"].tolist()
    frame, metadata = correct_grid(frame, metadata, args.analysis_root)
    frame = normalized_frame(frame)
    if frame.loc[frame.split.eq("test"), "id"].tolist() != original_test:
        raise ValueError("Kinship correction changed the outer test roster")
    metadata.update(publication_metadata(args.analysis_root))
    metadata.setdefault("groups_verified", False)
    metadata.setdefault("grouping_status", "provided_identity_unverified")
    metadata.setdefault("gwas_overlap_status", "unverified")
    block_cfg = {"block_config": config["block_config"], "execution": {"device": args.device}}
    # Existing producer-bound blocks can be verified without rehydrating dosage.
    # Cold input preparation is shared infrastructure only, never model fitting.
    import block_inputs as b
    bc = b.read_config(config["block_config"])
    layout = b.resolve_block_files(args.ld_block_dir or bc["ld_block_dir"], args.ld_block_files or bc["ld_block_files"])
    block_ready = b.cache_paths(bc, layout, trait)[-1].is_file()
    if not block_ready and not (Path(bc["prsformer_cache"])/"genotypes.npy").is_file():
        from shared_cohort import ensure_inputs
        ensure_inputs(SimpleNamespace(**dict(vars(args), traits=trait, complete_features=False)), config)
    array, cohort, tokens, receipt, block_path = locate_blocks(
        block_cfg, trait, args.ld_block_files, args.ld_block_dir, allow_build=True)
    source_ids = cohort.eid.astype(str)
    if source_ids.duplicated().any() or not set(frame.id) <= set(source_ids):
        raise ValueError("Block producer does not uniquely cover the analysis roster")
    rows = pd.Series(np.arange(len(cohort)), index=source_ids).loc[frame.id].to_numpy(int)
    maps = block_maps(tokens, array.shape[2])
    records = []
    for index, (_, token) in enumerate(tokens.iterrows()):
        for channel in range(array.shape[2]):
            records.append(dict(feature_id=f"{token['block_id']}:{maps[index*array.shape[2]+channel]}",
                                map_id=str(maps[index*array.shape[2]+channel]),
                                block_id=str(token["block_id"]),
                                chromosome=str(token.get("CHR", "")),
                                start=int(token.get("start", 0)), end=int(token.get("end", 0))))
    features = pd.DataFrame(records)
    if features.feature_id.duplicated().any():
        raise ValueError("Nonunique block/channel feature identity")
    metadata.update(source_kind="producer_bound_raw_csx_ldblocks", block_file=str(block_path),
                    block_producer_version=receipt.get("producer_version"),
                    block_signature=receipt.get("signature"),
                    block_array_file_record=receipt["array_binding"].get("array_file_record"),
                    whole_csx_score_file_record=receipt.get("signature_payload", {}).get("whole_score_file_record"),
                    block_weight_fusion=receipt.get("shared_weight_fusion"),
                    block_map_semantics=receipt.get("shared_map_semantics"),
                    input_source=str(prepared), input_file_record=file_record(prepared),
                    endpoint_definition="baseline_t2dm_from_Yr2e_Yt2e" if trait == "t2dm" else trait)
    # Upstream metadata also contains large evolution DataFrames. Never embed
    # those or arbitrary individual tables in the aggregate model manifest.
    safe_keys = ("source_kind", "block_file", "block_producer_version", "block_signature",
                 "block_array_file_record", "whole_csx_score_file_record", "block_weight_fusion", "block_map_semantics",
                 "input_source", "input_file_record", "groups_verified", "grouping_status",
                 "gwas_overlap_status", "endpoint_definition", "t2dm_endpoint", "phenotype_definitions",
                 "kinship_version", "kinship_token", "kinship_correction_sha256", "excluded_development_relatives")
    safe_metadata = {key: metadata[key] for key in safe_keys if key in metadata}
    return frame, array, rows, features, "binary" if trait == "t2dm" else "continuous", safe_metadata


class Rows:
    """Bounded access to an array in analysis-roster order, flattened by token."""
    def __init__(self, array, rows):
        self.array, self.rows = array, np.asarray(rows, dtype=int)
        self.shape = (len(self.rows), int(np.prod(array.shape[1:])))

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, index):
        value = np.asarray(self.array[self.rows[index]])
        return value.reshape(len(value), -1)
