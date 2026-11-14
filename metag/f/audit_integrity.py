"""Optional, separate integrity seals; legacy workflow receipts are never rewritten.

Report mode reads receipts and checks recorded paths/sizes only. Seal mode reads
every selected file once to establish a SHA256 baseline. When a legacy receipt
has no checksum, that baseline cannot establish historical correctness. Quick
verification checks size/mtime/ctime; --full additionally rehashes file contents.
Run seal/full verification on a compute node for large reference collections.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import stat
import tempfile

from common import fingerprint, load_config, require_compute, sha256, utcnow

SCHEMA = "metag-integrity-seal-v1"
OWNER = "metag.audit_integrity"
BASELINE_ONLY = "baseline_only_no_historical_checksum"
HISTORICAL = "matches_recorded_sha256"


def _state(path):
    path = Path(path)
    value = path.lstat()
    if not stat.S_ISREG(value.st_mode):
        raise RuntimeError(f"Expected a regular, non-symlink file: {path}")
    return {"size": value.st_size, "mtime_ns": value.st_mtime_ns,
            "ctime_ns": value.st_ctime_ns}


def _read(path):
    before = _state(path)
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Invalid JSON receipt: {path}: {exc}") from exc
    if not isinstance(data, dict) or _state(path) != before:
        raise RuntimeError(f"Receipt is not an object or changed while reading: {path}")
    return data, before


def _receipt(path, schema=None):
    data, state = _read(path)
    context = data.get("context")
    if (not isinstance(context, dict) or data.get("status") != "complete" or
            data.get("identity") != fingerprint(context) or
            (schema is not None and data.get("schema") != schema)):
        raise RuntimeError(f"Receipt has no valid complete context: {path}")
    return data, state


def _child(root, relative):
    root = Path(root).resolve()
    relative = Path(relative)
    candidate = root / relative
    if (relative.is_absolute() or ".." in relative.parts or not relative.parts or
            root not in candidate.resolve().parents):
        raise RuntimeError(f"Unsafe recorded relative path: {root} / {relative}")
    return candidate


def _checksum(value, path):
    if value is not None and not re.fullmatch(r"[0-9a-f]{64}", str(value)):
        raise RuntimeError(f"Invalid recorded SHA256: {path}")
    return value


def _add(entries, path, role, record=None, captured_state=None):
    path = Path(path).absolute()
    state = _state(path)
    if captured_state is not None and state != captured_state:
        raise RuntimeError(f"File changed after its receipt was read: {path}")
    record = record or {}
    recorded_size = record.get("state", {}).get("size", state["size"])
    if state["size"] != recorded_size:
        raise RuntimeError(f"Recorded size mismatch: {path}")
    old_hash = _checksum(record.get("sha256"), path)
    key = str(path)
    if key in entries:
        entry = entries[key]
        if entry["state"] != state or (entry.get("recorded_sha256") and old_hash and
                                      entry["recorded_sha256"] != old_hash):
            raise RuntimeError(f"Conflicting records or changing file: {path}")
        entry["roles"] = sorted(set(entry["roles"]) | {role})
        if old_hash:
            entry["recorded_sha256"] = old_hash
    else:
        entries[key] = {"path": key, "roles": [role], "state": state,
                        "recorded_sha256": old_hash}


def _paths(cfg):
    result = {key: str(Path(cfg[key]).absolute()) for key in
              ("reference_raw", "reads_raw", "metadata_dir", "library_dir", "work_dir")}
    result.update(db_dir=str(Path(cfg.get("db_dir") or Path(result["library_dir"]) / "db").absolute()),
                  kunpeng_db_dir=str(Path(cfg.get("kunpeng_db_dir") or Path(result["library_dir"]) / "kunpeng").absolute()),
                  prepared_dir=str(Path(cfg.get("prepared_dir") or Path(result["library_dir"]) / "prepared").absolute()))
    return result


def _index(entries, inventories, directory, name):
    directory = Path(directory)
    marker = directory / "database.json"
    if not marker.exists():
        return {"status": "not_complete", "directory": str(directory)}
    receipt, state = _receipt(marker, "metag-hpc-database-v1" if name == "kraken" else 1)
    context = receipt["context"]
    if receipt.get("db_dir") != str(directory) or receipt.get("taxonomy_dir") != str(directory / "taxonomy"):
        raise RuntimeError(f"Index receipt belongs to a different directory: {marker}")
    files, taxonomy = receipt.get("files"), receipt.get("taxonomy_files")
    if (not isinstance(files, dict) or not isinstance(taxonomy, dict) or
            files != context.get("files") or taxonomy != context.get("taxonomy_files") or
            receipt.get("parameters") != context.get("parameters") or
            not {"nodes.dmp", "names.dmp", "merged.dmp"} <= set(taxonomy)):
        raise RuntimeError(f"Index receipt file lists are incomplete/inconsistent: {marker}")
    if name == "kraken":
        expected = {"hash.k2d", "taxo.k2d", "opts.k2d"}
    else:
        partitions = receipt.get("parameters", {}).get("partitions")
        if isinstance(partitions, bool) or not isinstance(partitions, int) or partitions < 1:
            raise RuntimeError(f"Invalid Kun-Peng partition count: {marker}")
        expected = {"hash_config.k2d", "taxo.k2d", "opts.k2d"} | {
            f"hash_{i}.k2d" for i in range(1, partitions + 1)}
    if set(files) != expected:
        raise RuntimeError(f"Unexpected or missing native index records: {marker}")
    observed = {p.name for p in directory.glob("hash*.k2d")} | {
        p.name for p in directory.glob("*.k2d") if p.name in ("taxo.k2d", "opts.k2d")}
    if observed != expected:
        raise RuntimeError(f"Unexpected or missing native index files: {directory}")
    inventories.append({"kind": "index", "directory": str(directory), "paths": sorted(expected)})
    _add(entries, marker, name + "_receipt", captured_state=state)
    for filename, record in files.items():
        _add(entries, _child(directory, filename), name + "_index", record)
    for filename, record in taxonomy.items():
        _add(entries, _child(directory / "taxonomy", filename), name + "_taxonomy", record)
    return {"status": "complete", "identity": receipt["identity"], "directory": str(directory)}


def _references(entries, inventories, cfg, paths):
    from reference import _inventory
    lib, meta = Path(paths["library_dir"]), Path(paths["metadata_dir"])
    metadata_path = meta / "metadata.json"
    metadata, state = _receipt(metadata_path, "metag-hpc-reference-v1")
    _add(entries, metadata_path, "metadata_receipt", captured_state=state)
    records = metadata["context"].get("files")
    required = {"taxonomy/nodes.dmp", "taxonomy/names.dmp", "taxonomy/merged.dmp", "reference_taxids.tsv"}
    if not isinstance(records, dict) or not required <= set(records):
        raise RuntimeError("Metadata receipt lacks required taxonomy/mapping files")
    for relative, record in records.items():
        _add(entries, _child(meta, relative), "reference_metadata", record)
    custom = metadata["context"].get("custom_map")
    configured_map = str(Path(cfg["reference_taxid_map"]).absolute()) if cfg.get("reference_taxid_map") else None
    if (custom.get("path") if custom else None) != configured_map:
        raise RuntimeError("Custom taxid map differs from frozen metadata")
    if custom:
        _add(entries, custom["path"], "custom_taxid_map", custom)
    plan_path = lib / "reference_plan.json"
    plan, state = _receipt(plan_path, "metag-hpc-reference-v1")
    _add(entries, plan_path, "reference_plan", captured_state=state)
    tasks = plan.get("tasks")
    if (not isinstance(tasks, list) or not tasks or tasks != plan["context"].get("tasks") or
            [task.get("task_id") for task in tasks] != list(range(1, len(tasks) + 1)) or
            plan["context"].get("metadata_identity") != metadata["identity"] or
            plan["context"].get("reference_raw") != str(Path(paths["reference_raw"]).resolve())):
        raise RuntimeError("Reference plan/tasks do not match the frozen metadata/raw directory")
    actual = _inventory(paths)
    if (fingerprint(actual) != plan["context"].get("raw_inventory") or
            fingerprint(actual) != metadata["context"].get("raw_inventory")):
        raise RuntimeError("Raw reference path/size inventory differs from its frozen plan")
    manifests_rows = []
    for task in tasks:
        record = task["manifest"]
        manifest_path = _child(lib, record["path"])
        manifest, state = _read(manifest_path)
        _add(entries, manifest_path, "reference_task_manifest", record, state)
        rows = manifest.get("files")
        if (not isinstance(rows, list) or manifest.get("task_id") != task["task_id"] or
                [row.get("index") for row in rows] != task["files"]):
            raise RuntimeError(f"Task manifest differs from plan: {manifest_path}")
        manifests_rows.extend(rows)
        for row in rows:
            _add(entries, _child(paths["reference_raw"], row["path"]), "raw_reference", row)
        prepared = Path(paths["prepared_dir"]) / f"task{task['task_id']:06d}"
        receipt_path = prepared / "receipt.json"
        receipt, state = _receipt(receipt_path, "metag-hpc-reference-v1")
        if (receipt.get("plan_identity") != plan["identity"] or
                receipt.get("task_id") != task["task_id"] or receipt["context"].get("task") != task):
            raise RuntimeError(f"Prepared task belongs to a different plan: {receipt_path}")
        _add(entries, receipt_path, "prepared_receipt", captured_state=state)
        prepared_files = receipt["context"].get("files")
        if not isinstance(prepared_files, dict) or set(prepared_files) != {"sequences.fna", "seqid2taxid.map"}:
            raise RuntimeError(f"Incomplete prepared task file records: {receipt_path}")
        for filename, record in prepared_files.items():
            _add(entries, _child(prepared, filename), "prepared_reference", record)
    rows = sorted(manifests_rows, key=lambda row: row["index"])
    if ([row["index"] for row in rows] != list(range(len(actual))) or
            [{key: row[key] for key in ("path", "state", "accession")} for row in rows] != actual):
        raise RuntimeError("Task manifests do not cover every and only frozen raw reference")
    inventories.append({"kind": "raw_reference", "directory": paths["reference_raw"],
                        "paths": [row["path"] for row in actual]})


def collect_files(cfg, scope="indexes"):
    """Return a validated inventory; this does not hash data or claim content integrity."""
    if scope not in ("indexes", "references"):
        raise ValueError("scope must be indexes or references")
    paths = _paths(cfg)
    entries, inventories, indexes = {}, [], {}
    for name, key in (("kraken", "db_dir"), ("kunpeng", "kunpeng_db_dir")):
        indexes[name] = _index(entries, inventories, paths[key], name)
    if scope == "references":
        _references(entries, inventories, cfg, paths)
    if not entries:
        raise RuntimeError("No complete recorded files available; finish a database before sealing indexes")
    return {"scope": scope, "paths": paths, "indexes": indexes,
            "inventories": inventories, "files": [entries[key] for key in sorted(entries)]}


def report(cfg, scope="indexes"):
    result = collect_files(cfg, scope)
    result.update(schema=SCHEMA, mode="report", content_hashes_checked=False,
                  assurance="receipt_structure_and_recorded_paths_sizes_only",
                  historical_integrity="not_established", files_count=len(result["files"]),
                  bytes_total=sum(row["state"]["size"] for row in result["files"]))
    return result


def _seal_destination(path, inventory):
    path = Path(path).expanduser().absolute()
    resolved = path.resolve()
    for key in ("reference_raw", "reads_raw", "metadata_dir", "db_dir", "kunpeng_db_dir", "prepared_dir"):
        protected = Path(inventory["paths"][key]).resolve()
        if resolved == protected or protected in resolved.parents or resolved in protected.parents:
            raise RuntimeError(f"Integrity seal must be outside protected input/index/metadata/prepared paths: {path}")
    if path.exists() or path.is_symlink():
        raise RuntimeError(f"Refusing to overwrite an existing seal; choose a new filename: {path}")
    return path


def _publish(path, document):
    """Atomic, no-overwrite publication on the same filesystem."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".integrity-seal-", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as out:
            json.dump(document, out, sort_keys=True, indent=2)
            out.write("\n")
            out.flush()
            os.fsync(out.fileno())
        os.link(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        Path(temporary).unlink(missing_ok=True)


def seal(cfg, destination, scope="indexes", *, allow_local=False):
    require_compute(allow_local)
    inventory = collect_files(cfg, scope)
    destination = _seal_destination(destination, inventory)
    for row in inventory["files"]:
        path = Path(row["path"])
        if _state(path) != row["state"]:
            raise RuntimeError(f"File changed before hashing: {path}")
        checksum = sha256(path)
        if _state(path) != row["state"]:
            raise RuntimeError(f"File changed while hashing: {path}")
        if row["recorded_sha256"] and checksum != row["recorded_sha256"]:
            raise RuntimeError(f"Historical recorded SHA256 mismatch: {path}; no seal was published")
        row["sha256"] = checksum
        row["historical_evidence"] = HISTORICAL if row["recorded_sha256"] else BASELINE_ONLY
    # Catch files replaced early in a long audit, and changed directory membership.
    _verify_inventory(inventory, full=False)
    context = {**inventory, "created_at": utcnow(),
               "meaning": "Current content baseline; files without historical checksums have no proof of earlier correctness."}
    document = {"schema": SCHEMA, "owner": OWNER, "context": context,
                "identity": fingerprint(context)}
    _publish(destination, document)
    return {"seal_path": str(destination), "identity": document["identity"],
            "scope": scope, "files_count": len(context["files"]),
            "baseline_only_files": sum(row["historical_evidence"] == BASELINE_ONLY for row in context["files"]),
            "content_hashes_checked": True,
            "meaning": context["meaning"]}


def _verify_inventory(context, full):
    seen = set()
    for row in context["files"]:
        path = Path(row["path"])
        if not path.is_absolute() or str(path) in seen:
            raise RuntimeError("Integrity seal has duplicate or relative file paths")
        seen.add(str(path))
        if _state(path) != row["state"]:
            raise RuntimeError(f"Sealed file size/mtime/ctime changed: {path}; perform an explicit fresh content audit before accepting new content")
        if full:
            if sha256(path) != row["sha256"]:
                raise RuntimeError(f"Sealed SHA256 mismatch: {path}")
            if _state(path) != row["state"]:
                raise RuntimeError(f"Sealed file changed during full verification: {path}")
    for inventory in context["inventories"]:
        directory = Path(inventory["directory"])
        if inventory["kind"] == "index":
            actual = {p.name for p in directory.glob("hash*.k2d")} | {
                p.name for p in directory.glob("*.k2d") if p.name in ("taxo.k2d", "opts.k2d")}
        elif inventory["kind"] == "raw_reference":
            from reference import _inventory
            actual = {row["path"] for row in _inventory(context["paths"])}
        else:
            raise RuntimeError("Unknown sealed inventory kind")
        if actual != set(inventory["paths"]):
            raise RuntimeError(f"Sealed directory membership changed: {directory}; create a new seal after explicit audit")
    if full:
        for row in context["files"]:
            if _state(row["path"]) != row["state"]:
                raise RuntimeError(f"Sealed file changed during the full audit: {row['path']}")


def verify_seal(path, full=False):
    """Read-only API. Raises RuntimeError on mutation; quick mode performs no data hashing."""
    document, original_state = _read(path)
    context = document.get("context")
    if (document.get("schema") != SCHEMA or document.get("owner") != OWNER or
            not isinstance(context, dict) or fingerprint(context) != document.get("identity") or
            context.get("scope") not in ("indexes", "references") or
            not isinstance(context.get("files"), list) or not context["files"] or
            not isinstance(context.get("inventories"), list)):
        raise RuntimeError("Integrity seal schema/owner/fingerprint is invalid")
    try:
        for row in context["files"]:
            if _checksum(row["sha256"], row["path"]) is None:
                raise RuntimeError("Missing sealed SHA256")
            if (set(row["state"]) != {"size", "mtime_ns", "ctime_ns"} or
                    row["historical_evidence"] not in (HISTORICAL, BASELINE_ONLY)):
                raise RuntimeError("Malformed sealed file state/evidence")
        _verify_inventory(context, full)
    except (KeyError, TypeError, ValueError, OSError) as exc:
        raise RuntimeError(f"Invalid or missing sealed file/record: {exc}") from exc
    if _state(path) != original_state:
        raise RuntimeError("Integrity seal changed during verification")
    return {"status": "pass", "identity": document["identity"], "scope": context["scope"],
            "paths": context["paths"], "indexes": context["indexes"],
            "files_count": len(context["files"]), "content_hashes_checked": bool(full),
            "assurance": "sha256_and_file_states_match_seal" if full else "file_states_match_seal_no_content_rehash",
            "baseline_only_files": sum(row["historical_evidence"] == BASELINE_ONLY for row in context["files"])}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config")
    parser.add_argument("--root")
    parser.add_argument("--scope", choices=("indexes", "references"), default="indexes")
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--report", action="store_true", help="Default: inspect receipts and sizes; no data content hashing")
    modes.add_argument("--seal", metavar="PATH", help="Hash selected files once and publish a new separate baseline")
    modes.add_argument("--verify-seal", metavar="PATH", help="Check an existing seal without modifying data or receipts")
    parser.add_argument("--full", action="store_true", help="With --verify-seal, also rehash every sealed file")
    parser.add_argument("--allow-local", action="store_true", help="Explicitly permit expensive local seal/full verification")
    args = parser.parse_args(argv)
    if args.full and not args.verify_seal:
        parser.error("--full requires --verify-seal")
    try:
        if args.verify_seal:
            if args.full:
                require_compute(args.allow_local)
            result = verify_seal(args.verify_seal, args.full)
        else:
            cfg = load_config(args.config, args.root)
            result = seal(cfg, args.seal, args.scope, allow_local=args.allow_local) if args.seal else report(cfg, args.scope)
        print(json.dumps(result, indent=2, sort_keys=True))
    except (RuntimeError, ValueError, KeyError, OSError) as exc:
        parser.exit(1, f"Integrity audit failed: {exc}\n")


if __name__ == "__main__":
    main()
