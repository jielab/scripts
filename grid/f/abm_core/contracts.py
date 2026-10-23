"""Bind label-free arrays/tables to their row identities and frozen producers.

GRID adapters must write these sidecars at the time they export the arrays or
features, from the SAME ordered frame. A checksum cannot prove a producer used
the correct scientific algorithm; it prevents accidental mixing/reordering of
already exported files. Source-code and role provenance are separate records.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

import numpy as np

from .model import IDENTITY_FIELDS, file_sha256, write_json


def ordered_ids_sha256(ids):
    values = [str(x) for x in ids]
    payload = json.dumps(values, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def token_sha256(tokens):
    """Canonical tabular content, including column and row order."""
    return hashlib.sha256(tokens.to_csv(index=False, float_format="%.17g",
                                       lineterminator="\n").encode()).hexdigest()

def binding_metadata(frame):
    """Use supplied group identity when available, without claiming kinship QC."""
    out = frame.copy()
    if "id" not in out and "eid" in out:
        out = out.rename(columns={"eid": "id"})
    if "group_id" not in out:
        name = next((c for c in ("family_id", "group") if c in out), None)
        if name:
            out = out.rename(columns={name: "group_id"})
    return out


def array_record(array_path, metadata, *, geometry_sha256, tokens=None,
                 channels=None, producer_version=None):
    x = np.load(array_path, mmap_mode="r")
    if x.ndim not in (2, 3) or len(x) != len(metadata):
        raise ValueError("Array and metadata dimensions disagree")
    if metadata.id.isna().any() or metadata.id.astype(str).duplicated().any():
        raise ValueError("Array producer requires unique nonmissing IDs")
    if "group_id" in metadata and (metadata.group_id.isna().any() or metadata.group_id.astype(str).str.strip().eq("").any()):
        raise ValueError("Provided producer groups cannot be missing")
    record = {"schema_version": "abm.array.2",
              "array_sha256": file_sha256(array_path), "shape": list(x.shape),
              "dtype": str(x.dtype),
              "axes": ["individual", "token", "channel"] if x.ndim == 3 else ["individual", "feature"],
              "ordered_ids_sha256": ordered_ids_sha256(metadata.id),
              "ordered_groups_sha256": (ordered_ids_sha256(metadata.group_id)
                                        if "group_id" in metadata else None),
              "group_status": "provided_identity_only" if "group_id" in metadata else "not_provided",
              "geometry_sha256": geometry_sha256}
    if x.ndim == 3:
        if tokens is None or len(tokens) != x.shape[1] or channels is None or len(channels) != x.shape[2]:
            raise ValueError("Raw blocks require ordered tokens and channels")
        if not producer_version:
            raise ValueError("Raw blocks require a declared producer scientific version")
        record.update(tokens_sha256=token_sha256(tokens), channels=list(channels),
                      producer_version=producer_version)
    return record


def write_array_manifest(array_path, metadata, *, geometry_sha256, output=None,
                         tokens=None, channels=None, producer_version=None):
    record = array_record(array_path, metadata, geometry_sha256=geometry_sha256,
                          tokens=tokens, channels=channels, producer_version=producer_version)
    path = Path(output) if output else Path(str(array_path) + ".rows.json")
    write_json(path, record)
    return path


def verify_array_binding(array_path, metadata, sidecar_path, *, tokens=None, channels=None):
    record = (sidecar_path if isinstance(sidecar_path, dict) else
              json.loads(Path(sidecar_path).read_text(encoding="utf-8")))
    if not {"array_sha256", "shape", "ordered_ids_sha256"} <= set(record):
        raise ValueError("Missing trustworthy producer binding; rebuild from verified inputs")
    if record["array_sha256"] != file_sha256(array_path):
        raise ValueError("Genetic array bytes changed after identity binding")
    if record["ordered_ids_sha256"] != ordered_ids_sha256(metadata.id):
        raise ValueError("Genetic array row IDs are reordered or mismatched")
    if record.get("ordered_groups_sha256") is not None and (
            "group_id" not in metadata or
            record["ordered_groups_sha256"] != ordered_ids_sha256(metadata.group_id)):
        raise ValueError("Genetic array group identities changed")
    x = np.load(array_path, mmap_mode="r")
    if list(x.shape) != record["shape"]:
        raise ValueError("Genetic array shape changed")
    if x.ndim == 3 and record.get("schema_version") != "abm.array.2":
        raise ValueError("Raw blocks lack a trustworthy producer binding; rebuild from genotype/weights")
    if record.get("schema_version") == "abm.array.2":
        axes = ["individual", "token", "channel"] if x.ndim == 3 else ["individual", "feature"]
        if record["dtype"] != str(x.dtype) or record["axes"] != axes:
            raise ValueError("Genetic array dtype/axes changed")
        if x.ndim == 3:
            if tokens is None or token_sha256(tokens) != record["tokens_sha256"]:
                raise ValueError("Raw block token content/order changed")
            if channels is None or list(channels) != record["channels"]:
                raise ValueError("Raw block channel order changed")
    return record


def write_feature_manifest(path, *, manifest, output=None):
    record = {k: manifest[k] for k in IDENTITY_FIELDS}
    record["features_sha256"] = file_sha256(path)
    record["query_outcomes_used"] = False
    target = Path(output) if output else Path(str(path) + ".provenance.json")
    write_json(target, record)
    return target


def verify_feature_binding(path, sidecar_path, *, manifest):
    record = json.loads(Path(sidecar_path).read_text(encoding="utf-8"))
    if record["features_sha256"] != file_sha256(path):
        raise ValueError("Query features changed after freezing")
    if record.get("query_outcomes_used") is not False:
        raise ValueError("Query feature producer did not declare label-free construction")
    for key in IDENTITY_FIELDS:
        if record.get(key) != manifest[key]:
            raise ValueError(f"Query feature producer differs from frozen model: {key}")
    return record
