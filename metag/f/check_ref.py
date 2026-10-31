#!/usr/bin/env python3
"""Reusable local reference preparation; cache only completed, stable checks."""
from __future__ import annotations

import fcntl
import hashlib
import json
import stat
import time
from datetime import datetime, timezone
from pathlib import Path

from common import ROOT, DIR_REF, DIR_LIB, CFG, atomic_json, assert_reference_policy, reference_files
from local_metadata import TAXONOMY_FILES, prepare_taxonomy, summary_candidates, taxonomy_path
import ref_prepare


def file_state(path):
    try:
        info = path.stat()
    except FileNotFoundError:
        return None
    return [info.st_mode, info.st_size, info.st_mtime_ns, info.st_ctime_ns, info.st_dev, info.st_ino]


def metadata_state():
    # Include missing candidates too: adding a higher-priority or supplemental
    # summary must invalidate the previous selection.
    paths = {path for group in (*ref_prepare.SUMMARY_URLS, "human")
             for path in summary_candidates(group)}
    source = taxonomy_path()
    paths.update(source/name for name in (*TAXONOMY_FILES, "delnodes.dmp", "taxdump.tar.gz"))
    paths.update(path.with_name(path.name + ".aria2") for path in list(paths))
    return {str(path): file_state(path) for path in sorted(paths)}


def reference_state():
    sizes, states, incomplete = {}, {}, []
    files = reference_files(sizes=sizes, states=states, incomplete=incomplete)
    if incomplete:
        raise RuntimeError(f"Reference download still incomplete: {incomplete[0]}; finish downloads first")
    return states, (files, sizes, incomplete)


def implementation_state():
    return {name: hashlib.sha256((ROOT/"f"/name).read_bytes()).hexdigest()
            for name in ("check_ref.py", "ref_prepare.py", "local_metadata.py", "make_refwork.py", "common.py")}


def output_state():
    out = DIR_LIB/"reference"
    paths = [out/name for name in ("inventory.tsv", "assembly_summary.txt", "reference.receipt.json",
                                  "taxid-validation.json", "taxid-excluded.tsv")]
    result = {str(path): file_state(path) for path in paths}
    for path, state in result.items():
        if state is None or not stat.S_ISREG(state[0]):
            return None
    return result


def load_cache(path):
    try:
        value = json.loads(path.read_text())
        return value if isinstance(value, dict) and value.get("status") == "complete" else {}
    except (OSError, ValueError):
        return {}


def main(force=False):
    assert_reference_policy()
    policy = CFG.get("DELETED_TAXID_POLICY", "error").lower()
    if policy not in ("error", "exclude"):
        raise ValueError("DELETED_TAXID_POLICY must be error or exclude")
    folder = DIR_LIB/"reference"
    folder.mkdir(parents=True, exist_ok=True)
    # A failed/interrupted check never publishes a reusable success receipt.
    with (folder/"check-ref.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print("[check-ref] Waiting for another reference check to finish...", flush=True)
            fcntl.flock(lock, fcntl.LOCK_EX)
        return prepare(force, policy, folder/"check-ref.receipt.json")


def prepare(force, policy, cache_path):
    started = time.monotonic()
    saved = load_cache(cache_path)
    context = {"schema": 1, "reference_root": str(DIR_REF), "library_root": str(DIR_LIB),
               "deleted_taxid_policy": policy, "implementation": implementation_state()}
    archive = taxonomy_path()/"taxdump.tar.gz"
    before_metadata = metadata_state()
    if any(value is not None and path.endswith(".aria2") for path, value in before_metadata.items()):
        raise RuntimeError("Metadata download is incomplete (.aria2 present); finish downloads before check-ref")
    old_metadata = saved.get("metadata", {})
    archive_changed = bool(old_metadata) and old_metadata.get(str(archive)) != before_metadata.get(str(archive))
    prepare_taxonomy(include_deleted=policy == "exclude", force=archive_changed or force)
    before_metadata = metadata_state()
    print("[check-ref] Compare reference file names, sizes and timestamps; no FASTA contents are read", flush=True)
    before_references, scan = reference_state()
    outputs = output_state()
    reusable = (saved.get("context") == context and saved.get("metadata") == before_metadata
                and saved.get("references") == before_references
                and outputs is not None and saved.get("outputs") == outputs)
    if reusable and not force:
        print(f"REUSE check-ref: {saved['active_references']:,} active references; "
              f"{saved['excluded_references']} excluded; metadata/taxid checks skipped "
              f"({time.monotonic() - started:.1f}s)", flush=True)
        return 0

    reason = "forced" if force else "no valid receipt or inputs/outputs changed"
    print(f"[check-ref] Full local check: {reason}", flush=True)
    cache_path.unlink(missing_ok=True)
    ref_prepare.main(include_taxonomy=True, deleted_policy=policy, reference_scan=scan)
    after_references, _ = reference_state()
    if before_references != after_references or before_metadata != metadata_state():
        raise RuntimeError("Reference files or metadata changed during check-ref; rerun after changes finish")
    outputs = output_state()
    if outputs is None:
        raise RuntimeError("Reference preparation outputs are incomplete; rerun check-ref")
    validation = json.loads((DIR_LIB/"reference/taxid-validation.json").read_text())
    atomic_json(cache_path, {
        "status": "complete", "checked_at": datetime.now(timezone.utc).isoformat(),
        "context": context, "metadata": before_metadata, "references": before_references,
        "outputs": outputs, "active_references": validation["active_references"],
        "excluded_references": validation["excluded_references"],
    })
    print(f"PASS: check-ref cached ({time.monotonic() - started:.1f}s); future builds reuse this result", flush=True)
    return 0


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    main(parser.parse_args().force)
