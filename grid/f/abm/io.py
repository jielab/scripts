#!/usr/bin/env python3
"""Bind label-free arrays/tables to their row identities and frozen producers.

GRID adapters must write these sidecars at the time they export the arrays or
features, from the SAME ordered frame. A checksum cannot prove a producer used
the correct scientific algorithm; it prevents accidental mixing/reordering of
already exported files. Source-code and role provenance are separate records.
"""
from __future__ import annotations

# Contracts
import hashlib
import json
from pathlib import Path

import numpy as np

from abm.model import IDENTITY_FIELDS, file_record, write_json
from runtime_state import file_identity, artifact_matches


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
              "array_file_record": file_record(array_path), "file_identity": file_identity(array_path), "shape": list(x.shape),
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
    if not {"array_file_record", "shape", "ordered_ids_sha256"} <= set(record):
        raise ValueError("Missing trustworthy producer binding; rebuild from verified inputs")
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
    record["features_file_record"] = file_record(path)
    record["query_outcomes_used"] = False
    target = Path(output) if output else Path(str(path) + ".provenance.json")
    write_json(target, record)
    return target


def verify_feature_binding(path, sidecar_path, *, manifest):
    record = json.loads(Path(sidecar_path).read_text(encoding="utf-8"))
    if record["features_file_record"] != file_record(path):
        raise ValueError("Query features changed after freezing")
    if record.get("query_outcomes_used") is not False:
        raise ValueError("Query feature producer did not declare label-free construction")
    for key in IDENTITY_FIELDS:
        if record.get(key) != manifest[key]:
            raise ValueError(f"Query feature producer differs from frozen model: {key}")
    return record

# Storage
"""Release completed support scratch arrays without changing frozen receipts.

Run with PYTHONPATH=f python -m abm.io; add --apply to delete.
The donor bank, block inputs and shared genotypes are persistent dependencies.
"""

# Shared root defaults come from 0f/0os.sh.
from file_lock import flock as shared_flock
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.append(str(_paths_root / "0f"))
from workflow_paths import configured_path

import argparse
from contextlib import contextmanager
import fcntl
import os
import time

from resume import atomic_json

RELEASE_FILE = "SCRATCH_RELEASED.json"
SCRATCH = ("geometry.npy", "radius_geometry.npy")
# Give exiting processes up to ~2 seconds to disappear or become inspectable.
# WSL can keep /proc maps/fd inaccessible for longer than a few milliseconds.
_PROC_RETRY_DELAYS = (.01, .05, .1, .25, .5, 1.)


@contextmanager
def support_lock(cache):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    with (cache / ".support.lock").open("a+") as stream:
        shared_flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def busy_reason(paths):
    """Also protect legacy writers that do not hold .support.lock (including self)."""
    paths = [Path(p) for p in paths]
    # A new/fully cleaned cache has no arrays to overwrite or unlink. An
    # unrelated, uninspectable process must not block its first build; callers
    # still hold support_lock to serialize current writers. Keep dangling
    # symlinks in the conservative path as well as existing files.
    if not any(p.exists() or p.is_symlink() for p in paths):
        return None
    targets = {str(p.resolve()) for p in paths}
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit():
            continue
        for attempt in range(len(_PROC_RETRY_DELAYS) + 1):
            try:
                if proc.stat().st_uid != os.getuid():
                    break
                if (proc / "stat").read_text().rsplit(")", 1)[1].split()[0] in {"Z", "X"}:
                    break
                if (proc / "comm").read_text().strip() == "(sd-pam)":
                    break
                for line in (proc / "maps").read_text().splitlines():
                    fields = line.split(maxsplit=5)
                    if len(fields) == 6 and fields[5].removesuffix(" (deleted)") in targets:
                        return "mapped by PID " + proc.name
                for fd in (proc / "fd").iterdir():
                    try:
                        if os.readlink(fd).removesuffix(" (deleted)") in targets:
                            return "open in PID " + proc.name
                    except FileNotFoundError:
                        pass
            except (FileNotFoundError, ProcessLookupError):
                break
            except PermissionError:
                # /proc can deny maps/fd access during exit, before stat becomes
                # Z or the PID disappears. Recheck instead of treating that race
                # as a persistent, uninspectable cache user.
                try:
                    command = (proc / "cmdline").read_bytes().split(b"\0")
                except (FileNotFoundError, ProcessLookupError):
                    break
                except OSError:
                    command = []
                # The per-user service manager is nondumpable on WSL and is
                # not a GRID worker. Other persistent denials stay conservative.
                if command[:2] in ([b"/usr/lib/systemd/systemd", b"--user"],
                                   [b"/lib/systemd/systemd", b"--user"]):
                    break
                if attempt < len(_PROC_RETRY_DELAYS):
                    time.sleep(_PROC_RETRY_DELAYS[attempt])
                    continue
                return "cannot inspect PID " + proc.name
            break
    return None


def validate_support(cache):
    """Only a bound release receipt can excuse a missing scratch array."""
    cache = Path(cache)
    pointer = cache / "COMPLETE.json"
    saved = json.loads(pointer.read_text())
    released = {}
    if (cache / RELEASE_FILE).is_file():
        release = json.loads((cache / RELEASE_FILE).read_text())
        if release.get("complete_file_record") != file_record(pointer):
            raise ValueError("Support scratch release belongs to a different completion receipt")
        released = release["files"]
        if set(released) - set(SCRATCH):
            raise ValueError("Release receipt cannot remove persistent support files")
    required = {"support.parquet", "bank.parquet", "geometry_rows.parquet",
                "geometry.npy", "donor_geometry.npy", "neighbor_indices.npy",
                "neighbor_distances.npy", "geometry_manifest.json", "donor_bank_manifest.json",
                "geometry.npy.rows.json", "donor_geometry.npy.rows.json"}
    if not required <= saved["files"].keys():
        raise ValueError("Incomplete support receipt")
    for name, expected in saved["files"].items():
        path = cache / name
        if Path(name).name != name or path.is_symlink():
            raise ValueError("Unsafe support cache file: " + name)
        if not path.exists() and name in SCRATCH and released.get(name, {}).get("file_record") == expected:
            continue
        if not path.is_file() or file_record(path) != expected:
            raise ValueError("Shared support cache changed: " + name)
    return saved


def release_verified_scratch(cache, saved):
    """Caller holds support_lock and has produced or validated every saved file."""
    cache = Path(cache)
    targets = [cache / name for name in SCRATCH if (cache / name).exists()]
    if not targets:
        return 0
    if any(p.is_symlink() for p in targets):
        raise ValueError("Refusing to release symlinked support scratch")
    reason = busy_reason([*targets, cache / "donor_geometry.npy"])
    if reason:
        print(f"[GRID cleanup] RETAIN {cache}: {reason}", flush=True)
        return 0
    pointer_sha = file_record(cache / "COMPLETE.json")
    release_path = cache / RELEASE_FILE
    release = (json.loads(release_path.read_text()) if release_path.exists()
               else dict(complete_file_record=pointer_sha, files={}))
    if release.get("complete_file_record") != pointer_sha:
        raise ValueError("Support completion receipt changed during release")
    total = 0
    for path in targets:
        size = path.stat().st_size
        release["files"][path.name] = dict(bytes=size, file_record=saved["files"].get(path.name))
        total += size
    # Commit first: interruption before unlink is safe and can be retried.
    # COMPLETE.json stays byte-for-byte unchanged for existing model pointers.
    atomic_json(release_path, release)
    for path in targets:
        path.unlink()
    print(f"[GRID cleanup] Released {total / 2**30:.2f} GiB: {cache}", flush=True)
    return total


def clean_cache(cache, *, apply=False):
    cache = Path(cache)
    if cache.is_symlink() or (cache / "COMPLETE.json").is_symlink():
        raise ValueError("Refusing symlinked support cache")
    if not (cache / "COMPLETE.json").is_file():
        print(f"[GRID cleanup] SKIP incomplete support: {cache}", flush=True)
        return 0
    with support_lock(cache):
        targets = [cache / name for name in SCRATCH if (cache / name).exists()]
        if not targets:
            return 0
        reason = busy_reason([*targets, cache / "donor_geometry.npy"])
        if reason:
            print(f"[GRID cleanup] SKIP {cache}: {reason}", flush=True)
            return 0
        size = sum(p.stat().st_size for p in targets)
        if not apply:
            print(f"[GRID cleanup] CANDIDATE {size / 2**30:.2f} GiB (validation on --apply): {cache}", flush=True)
            return size
        print(f"[GRID cleanup] Verifying completed support: {cache}", flush=True)
        saved = validate_support(cache)
        return release_verified_scratch(cache, saved)


def clean_completed_support(analysis_root, *, traits=None, apply=False, best_effort=False):
    """Sweep old completion keys too, including when the runner reuses all phases."""
    total = 0
    for namespace in ("support", "SMOKE"):
        root = Path(analysis_root) / "tmp/cache/abm" / namespace
        for cache in sorted(root.glob("*/*/_LOCAL_ONLY")):
            if traits is not None and cache.parent.parent.name not in traits:
                continue
            try:
                total += clean_cache(cache, apply=apply)
            except BlockingIOError:
                print(f"[GRID cleanup] SKIP locked support: {cache}", flush=True)
            except (OSError, ValueError, KeyError, TypeError) as exc:
                if not best_effort:
                    raise
                # An unrelated historical cache must not break a new analysis
                # or mask its original exception. Failed validation retains it.
                print(f"[GRID cleanup] RETAIN {cache}: {exc}", flush=True)
    return total


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", default=configured_path('ANALYSIS_ROOT', 'grid'))
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args(argv)
    total = clean_completed_support(args.analysis_root, apply=args.apply)
    print(f"[GRID cleanup] {'Released' if args.apply else 'Candidates'} total: {total / 2**30:.2f} GiB", flush=True)



if __name__ == "__main__":
    main()
