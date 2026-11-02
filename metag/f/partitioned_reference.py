"""Portable, content-sealed reference inputs for parallel category indexes.

Planning only reads inventories, mappings and small JSON receipts. The compute
preparation jobs stream each FASTA once to establish SHA256 baselines, validate
record identities and, only for a task spanning categories, separate complete
records. Completed old tasks are never masked again. An old size-only receipt
cannot prove historical absence of bit rot; adoption records this limitation.

All category builders MUST receive the same catalog ``seqid_map`` and
``taxonomy_dir``. Mapping a subset changes Kraken's internal taxonomy IDs and
would make evidence from independently built hashes incompatible.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil

from hpc_common import atomic_json, check_space, file_state, fingerprint, locked, sha256, utcnow
import hpc_reference as legacy

SCHEMA = "metag-partitioned-reference-v1"
CATALOG_SCHEMA = "metag-partitioned-inputs-v1"
CATEGORIES = ("microbes", "viruses", "plants", "vertebrates", "invertebrates")
ALIASES = {"microbe": "microbes", "microbes": "microbes", "microbial": "microbes",
           "virus": "viruses", "viruses": "viruses", "viral": "viruses",
           "plant": "plants", "plants": "plants",
           "vertebrate": "vertebrates", "vertebrates": "vertebrates",
           "invertebrate": "invertebrates", "invertebrates": "invertebrates"}
POLICY = "whole-record-unique-id-lowercase-to-x-v1"
BASELINE_ONLY = "legacy_size_receipt_new_sha256_baseline_not_historical_proof"
HEADER = re.compile(rb">metag_f(\d{10})_r(\d{10})\|kraken:taxid\|(\d+)\n\Z")
HEX64 = re.compile(r"[0-9a-f]{64}\Z")


def _paths(cfg):
    def absolute(value):
        value = Path(value).expanduser()
        return (value if value.is_absolute() else Path(cfg.get("root", ".")) / value).resolve()
    library = absolute(cfg["library_dir"])
    old = cfg.get("legacy_library_dir")
    old = absolute(old) if old else (library if (library / "reference_plan.json").is_file() else None)
    metadata = absolute(cfg.get("legacy_metadata_dir") or cfg["metadata_dir"])
    return {"base": absolute(cfg.get("parallel_reference_dir") or library / "parallel"),
            "raw": absolute(cfg["reference_raw"]), "metadata": metadata, "legacy": old}


def _safe_child(root, relative):
    value = PurePosixPath(str(relative))
    if value.is_absolute() or not value.parts or any(x in ("..", ".") for x in value.parts):
        raise RuntimeError(f"Unsafe relative reference path: {relative!r}")
    target = Path(root).joinpath(*value.parts)
    # Root itself is an explicitly resolved configuration path. Descendants may
    # not silently redirect an owned artifact or a cache read elsewhere.
    if any(p.is_symlink() for p in (target, *list(target.parents)[:len(value.parts) - 1])):
        raise RuntimeError(f"Reference/cache descendants must not be symlinks: {target}")
    return target


def _sealed_json(path, schema=SCHEMA):
    document = legacy._read_json(path)
    if document.get("schema") != schema or document.get("status") != "complete" or \
            fingerprint(document.get("context")) != document.get("identity"):
        raise RuntimeError(f"Invalid complete checkpoint/fingerprint: {path}")
    return document


def _stat(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Expected a regular non-symlink cache file: {path}")
    st = path.stat()
    return {"size": st.st_size, "mtime_ns": st.st_mtime_ns, "ctime_ns": st.st_ctime_ns,
            "device": st.st_dev, "inode": st.st_ino}


def _directory(path, *, create=False):
    path = Path(path)
    if path.is_symlink() or (path.exists() and not path.is_dir()):
        raise RuntimeError(f"Expected an owned, non-symlink workflow directory: {path}")
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


def _small_record(path):
    before = _stat(path)
    digest = sha256(path)
    if before != _stat(path):
        raise RuntimeError(f"File changed while hashing: {path}")
    return {"size": before["size"], "sha256": digest}


def _check_old_record(path, record, full=False):
    if _stat(path)["size"] != record.get("state", {}).get("size"):
        raise RuntimeError(f"Legacy frozen file size changed: {path}")
    if full and record.get("sha256") and sha256(path) != record["sha256"]:
        raise RuntimeError(f"Legacy recorded SHA256 mismatch: {path}")


def _metadata(paths, inventory_id=None):
    root = paths["metadata"]
    document = _sealed_json(root / "metadata.json", legacy.SCHEMA)
    context = document["context"]
    required = {"reference_taxids.tsv", *("taxonomy/" + n for n in legacy.TAX_FILES)}
    if not required.issubset(context.get("files", {})):
        raise RuntimeError("Reference metadata lacks the mapping or full taxonomy files")
    for relative, record in context["files"].items():
        _check_old_record(_safe_child(root, relative), record)
    if inventory_id is not None and context.get("raw_inventory") != inventory_id:
        raise RuntimeError("Current raw relative inventory/size differs from frozen metadata")
    return document


def _mapping(paths, inventory):
    mapping = {}
    path = paths["metadata"] / "reference_taxids.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"path", "taxid"}.issubset(reader.fieldnames or ()):
            raise RuntimeError(f"Missing path/taxid columns: {path}")
        for row in reader:
            relative = row["path"]
            _safe_child(paths["raw"], relative)
            try:
                taxid = int(row["taxid"])
            except (TypeError, ValueError) as exc:
                raise RuntimeError(f"Invalid taxid in {path}: {row}") from exc
            if relative in mapping or taxid <= 0:
                raise RuntimeError(f"Duplicate reference or invalid taxid in {path}: {relative}")
            mapping[relative] = taxid
    if set(mapping) != {row["path"] for row in inventory}:
        raise RuntimeError("Reference mapping must cover every and only current raw FASTA file")
    return mapping


def _category(relative, cfg):
    overrides = cfg.get("reference_category_map") or {}
    if not isinstance(overrides, dict):
        raise RuntimeError("reference_category_map must be a path/prefix-to-category mapping")
    matches = [(str(prefix).rstrip("/"), name) for prefix, name in overrides.items()
               if relative == str(prefix).rstrip("/") or relative.startswith(str(prefix).rstrip("/") + "/")]
    value = max(matches, key=lambda item: len(item[0]))[1] if matches else relative.split("/", 1)[0]
    name = ALIASES.get(str(value).lower())
    if name is None:
        raise RuntimeError(f"Unknown reference category for {relative!r}: {value!r}; use the five category folders or an explicit reference_category_map entry")
    return name


def _mask_policy(cfg, software):
    if software.get("source_commit") != legacy.SOURCE_COMMIT or software.get("source_tag") != "2.17.2":
        raise RuntimeError("Reference preparation does not use the pinned native source commit")
    if bool(cfg.get("masking", True)) and not HEX64.fullmatch(str(software.get("files", {}).get("k2mask", ""))):
        raise RuntimeError("Masked legacy references lack a recorded k2mask binary SHA256")
    return {"masking": bool(cfg.get("masking", True)), "record_policy": POLICY,
            "source_commit": legacy.SOURCE_COMMIT, "source_tag": "2.17.2"}


def _legacy_plan(paths, cfg, inventory, mapping, meta):
    root = paths["legacy"]
    document = _sealed_json(root / "reference_plan.json", legacy.SCHEMA)
    context = document["context"]
    if document.get("tasks") != context.get("tasks") or \
            context.get("raw_inventory") != fingerprint(inventory) or \
            context.get("metadata_identity") != meta["identity"]:
        raise RuntimeError("Legacy reference plan differs from its metadata, task list or current raw inventory")
    if context.get("record_policy") != POLICY or context.get("masking") != bool(cfg.get("masking", True)):
        raise RuntimeError("Legacy record/masking policy differs; refusing to relabel prepared references")
    _mask_policy(cfg, context.get("software", {}))
    rows_by_index, tasks = {}, []
    for old_task in context["tasks"]:
        if old_task.get("task_id") != len(tasks) + 1:
            raise RuntimeError("Legacy task IDs must be contiguous starting at one")
        record = old_task["manifest"]
        manifest = _safe_child(root, record["path"])
        _check_old_record(manifest, record, full=True)
        payload = legacy._read_json(manifest)
        rows = payload.get("files", [])
        if payload.get("task_id") != old_task["task_id"] or [r.get("index") for r in rows] != old_task["files"]:
            raise RuntimeError("Legacy task manifest does not match its frozen plan")
        for row in rows:
            index = row.get("index")
            if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(inventory) or index in rows_by_index:
                raise RuntimeError("Legacy raw indices are duplicate/out of range")
            expected = {**inventory[index], "index": index, "taxid": mapping[inventory[index]["path"]]}
            if row != expected:
                raise RuntimeError(f"Legacy task raw inventory/taxid mismatch at index {index}")
            rows_by_index[index] = row
        if not rows:
            raise RuntimeError("Empty legacy preparation task")
        receipt_path = root / "prepared" / f"task{old_task['task_id']:06d}" / "receipt.json"
        imported = None
        if receipt_path.exists():
            receipt = _sealed_json(receipt_path, legacy.SCHEMA)
            rc = receipt["context"]
            if receipt.get("plan_identity") != document["identity"] or receipt.get("task_id") != old_task["task_id"] or \
                    rc.get("plan_identity") != document["identity"] or rc.get("task") != old_task or \
                    receipt.get("execution", {}).get("masking") != context["masking"]:
                raise RuntimeError(f"Incompatible legacy prepared receipt: {receipt_path}")
            if set(rc.get("files", {})) != {"sequences.fna", "seqid2taxid.map"}:
                raise RuntimeError(f"Unexpected legacy prepared files: {receipt_path}")
            for name, record in rc["files"].items():
                _check_old_record(receipt_path.parent / name, record)
            _validate_expected_stats(rc.get("stats", {}), rows)
            imported = {"plan_identity": document["identity"], "task": old_task,
                        "receipt_identity": receipt["identity"], "receipt_record": _small_record(receipt_path),
                        "stats": rc["stats"], "files": rc["files"],
                        "mask_binary_sha256": context["software"]["files"].get("k2mask"),
                        "historical_integrity": BASELINE_ONLY}
        tasks.append((rows, imported))
    if set(rows_by_index) != set(range(len(inventory))):
        raise RuntimeError("Legacy task manifests omit current raw references")
    return document, tasks


def _validate_expected_stats(stats, rows):
    sources = stats.get("sources", [])
    if len(sources) != len(rows):
        raise RuntimeError("Prepared source-counter coverage differs from its raw task")
    for source, row in zip(sources, rows):
        if source.get("path") != row["path"] or source.get("state") != row["state"] or \
                not isinstance(source.get("records"), int) or source["records"] <= 0 or \
                not isinstance(source.get("bases"), int) or source["bases"] <= 0:
            raise RuntimeError("Prepared source counters/path/state differ from the raw manifest")
    if stats.get("records") != sum(s["records"] for s in sources) or \
            stats.get("bases") != sum(s["bases"] for s in sources) or \
            not isinstance(stats.get("longest_record"), int) or stats["longest_record"] <= 0 or \
            not HEX64.fullmatch(str(stats.get("header_sha256", ""))):
        raise RuntimeError("Invalid prepared task counters/header checksum")


def plan(cfg):
    """Create/reuse a cheap plan. Metadata bootstrap must already be complete."""
    paths = _paths(cfg)
    base = paths["base"]
    base.mkdir(parents=True, exist_ok=True)
    with locked(base / "reference.lock"):
        inventory = legacy._inventory({**cfg, "reference_raw": str(paths["raw"])})
        metadata = _metadata(paths, fingerprint(inventory))
        if (base / "plan.json").exists():
            document = _load_plan(cfg)
            if document["context"]["raw_inventory"] != fingerprint(inventory):
                raise RuntimeError("Raw inventory changed after the parallel plan; use a new library")
            return document
        mapping = _mapping(paths, inventory)
        existing = None
        if paths["legacy"] is not None:
            existing, source_tasks = _legacy_plan(paths, cfg, inventory, mapping, metadata)
            policy = _mask_policy(cfg, existing["context"]["software"])
        else:
            policy = _mask_policy(cfg, legacy._software(cfg))
            target = float(cfg.get("prepare_target_gib", 2)) * legacy.GIB
            maximum = int(cfg.get("prepare_max_files", 1000))
            if not math.isfinite(target) or target <= 0 or maximum <= 0:
                raise RuntimeError("prepare_target_gib and prepare_max_files must be positive")
            source_tasks, pending, size, current_category = [], [], 0, None
            for index, item in enumerate(inventory):
                row = {**item, "index": index, "taxid": mapping[item["path"]]}
                category = _category(row["path"], cfg)
                if pending and (category != current_category or size + row["state"]["size"] > target or len(pending) >= maximum):
                    source_tasks.append((pending, None))
                    pending, size = [], 0
                pending.append(row)
                size += row["state"]["size"]
                current_category = category
            if pending:
                source_tasks.append((pending, None))
        manifests = _directory(base / "manifests", create=True)
        tasks, categories = [], set()
        for number, (rows, imported) in enumerate(source_tasks, 1):
            rows = [{**row, "category": _category(row["path"], cfg)} for row in rows]
            groups = [name for name in CATEGORIES if name in {row["category"] for row in rows}]
            categories.update(groups)
            content = {"task_id": number, "files": rows, "legacy": imported}
            path = manifests / f"task{number:06d}.json"
            if path.is_symlink() or (path.exists() and legacy._read_json(path) != content):
                raise RuntimeError(f"Refusing an unowned/incompatible existing task manifest: {path}")
            if not path.exists():
                atomic_json(path, content)
            tasks.append({"task_id": number, "mode": "adopt" if imported else "prepare",
                          "categories": groups, "compressed_bytes": sum(row["state"]["size"] for row in rows),
                          "manifest": {"path": path.relative_to(base).as_posix(), **_small_record(path)}})
        context = {"raw_inventory": fingerprint(inventory), "metadata_identity": metadata["identity"],
                   "mapping": _small_record(paths["metadata"] / "reference_taxids.tsv"),
                   "policy": policy, "tasks": tasks,
                   "categories": [name for name in CATEGORIES if name in categories],
                   "category_overrides": cfg.get("reference_category_map") or {},
                   "legacy_plan_identity": existing["identity"] if existing else None}
        document = {"schema": SCHEMA, "status": "complete", "context": context,
                    "identity": fingerprint(context), "tasks": tasks, "categories": context["categories"],
                    "references": len(inventory), "compressed_bytes": sum(row["state"]["size"] for row in inventory),
                    "created_at": utcnow()}
        atomic_json(base / "plan.json", document)
        return document


def _load_plan(cfg):
    paths = _paths(cfg)
    document = _sealed_json(paths["base"] / "plan.json")
    context = document["context"]
    if document.get("tasks") != context.get("tasks") or document.get("categories") != context.get("categories"):
        raise RuntimeError("Parallel reference plan aliases differ from the fingerprint")
    if context["policy"].get("masking") != bool(cfg.get("masking", True)) or \
            context.get("category_overrides", {}) != (cfg.get("reference_category_map") or {}):
        raise RuntimeError("Reference masking/category policy changed; use a new parallel reference directory")
    metadata = _metadata(paths)
    if metadata["identity"] != context["metadata_identity"]:
        raise RuntimeError("Parallel plan metadata differs from the current metadata receipt")
    mapping = paths["metadata"] / "reference_taxids.tsv"
    if _small_record(mapping) != context["mapping"]:
        raise RuntimeError("Resolved reference mapping changed after planning")
    if context["legacy_plan_identity"]:
        if paths["legacy"] is None:
            raise RuntimeError("This plan needs legacy_library_dir for its adopted cache")
        old = _sealed_json(paths["legacy"] / "reference_plan.json", legacy.SCHEMA)
        if old["identity"] != context["legacy_plan_identity"]:
            raise RuntimeError("Configured legacy cache has a different frozen reference plan")
    return document


def load_plan(cfg):
    """Validate the existing small plan without rescanning the raw inventory."""
    paths = _paths(cfg)
    with locked(paths["base"] / "reference.lock", shared=True):
        return _load_plan(cfg)


def _task_manifest(paths, task):
    path = _safe_child(paths["base"], task["manifest"]["path"])
    if _small_record(path) != {key: task["manifest"][key] for key in ("size", "sha256")}:
        raise RuntimeError(f"Parallel task manifest checksum changed: {path}")
    result = legacy._read_json(path)
    if result.get("task_id") != task["task_id"]:
        raise RuntimeError("Parallel task manifest ID mismatch")
    for row in result["files"]:
        if _stat(_safe_child(paths["raw"], row["path"]))["size"] != row["state"]["size"]:
            raise RuntimeError(f"Raw reference size changed: {row['path']}")
    return result


def _source_path(paths, task_id, source, folder=None):
    if source.get("kind") == "legacy":
        if paths["legacy"] is None:
            raise RuntimeError("Missing legacy_library_dir for an adopted prepared file")
        root = paths["legacy"] / "prepared" / f"task{int(source['task_id']):06d}"
    elif source.get("kind") == "local":
        root = folder or paths["base"] / "prepared" / f"task{int(task_id):06d}"
    else:
        raise RuntimeError("Unknown prepared source kind")
    return _safe_child(root, source["name"])


def _scan_prepared(sequence_path, mapping_path, rows, stage, *, split, masking):
    """Validate/map whole records and hash source/output bytes in one read pass."""
    before = {"sequence": _stat(sequence_path), "mapping": _stat(mapping_path)}
    groups = {name: {"records": 0, "bases": 0, "longest_record": 0, "size": 0,
                     "digest": hashlib.sha256()} for name in CATEGORIES if name in {r["category"] for r in rows}}
    outputs = {name: (stage / f"{name}.fna").open("wb") for name in groups} if split else {}
    source_hash, map_hash, headers = hashlib.sha256(), hashlib.sha256(), hashlib.sha256()
    sources = [{"path": r["path"], "state": r["state"], "records": 0, "bases": 0} for r in rows]
    positions = {r["index"]: i for i, r in enumerate(rows)}
    current, position, record_bases, longest = None, -1, 0, 0
    try:
        with Path(sequence_path).open("rb", buffering=8 << 20) as fasta, Path(mapping_path).open("rb", buffering=1 << 20) as mapping:
            for line in fasta:
                source_hash.update(line)
                if line.startswith(b">"):
                    if current is not None:
                        if not record_bases:
                            raise RuntimeError("Empty prepared FASTA record")
                        groups[current]["longest_record"] = max(groups[current]["longest_record"], record_bases)
                        longest = max(longest, record_bases)
                    match = HEADER.fullmatch(line)
                    if match is None:
                        raise RuntimeError("Prepared header is not a canonical metag source/record/taxid ID")
                    index, record, taxid = map(int, match.groups())
                    new_position = positions.get(index, -1)
                    if new_position < 0 or new_position not in (position, position + 1):
                        raise RuntimeError("Prepared source IDs are missing, repeated or out of order")
                    if new_position != position and position >= 0 and not sources[position]["records"]:
                        raise RuntimeError("Prepared reference contains an empty source")
                    position = new_position
                    row = rows[position]
                    if taxid != row["taxid"] or record != sources[position]["records"] + 1:
                        raise RuntimeError("Prepared taxid/record counter differs from the frozen raw manifest")
                    mapped = mapping.readline()
                    map_hash.update(mapped)
                    if mapped != line[1:-1] + b"\t" + str(taxid).encode() + b"\n":
                        raise RuntimeError("Prepared seqid-to-taxid map differs from its FASTA headers")
                    current = row["category"]
                    groups[current]["records"] += 1
                    sources[position]["records"] += 1
                    headers.update(line)
                    record_bases = 0
                else:
                    bases = line.rstrip(b"\r\n")
                    allowed = (legacy.DNA - frozenset(b"acgturyswkmbdhvn")) if masking else legacy.DNA
                    if current is None or not bases or not set(bases).issubset(allowed):
                        raise RuntimeError("Invalid/empty prepared FASTA sequence line")
                    groups[current]["bases"] += len(bases)
                    sources[position]["bases"] += len(bases)
                    record_bases += len(bases)
                    if masking and record_bases > legacy.MAX_MASK_RECORD:
                        raise RuntimeError("Legacy masked FASTA record exceeds the pinned k2mask signed-position limit; its masking correctness cannot be assumed")
                groups[current]["digest"].update(line)
                groups[current]["size"] += len(line)
                if split:
                    outputs[current].write(line)
            if current is None or not record_bases or position != len(rows) - 1 or any(not x["records"] for x in sources):
                raise RuntimeError("Prepared FASTA is empty/truncated or omits a raw source")
            groups[current]["longest_record"] = max(groups[current]["longest_record"], record_bases)
            longest = max(longest, record_bases)
            if mapping.read(1):
                raise RuntimeError("Prepared sequence map contains extra rows")
        for handle in outputs.values():
            legacy._sync_file(handle)
    finally:
        for handle in outputs.values():
            handle.close()
    if before != {"sequence": _stat(sequence_path), "mapping": _stat(mapping_path)}:
        raise RuntimeError("Prepared cache changed during checksum/record validation")
    for value in groups.values():
        value["sha256"] = value.pop("digest").hexdigest()
    stats = {"records": sum(s["records"] for s in sources), "bases": sum(s["bases"] for s in sources),
             "longest_record": longest, "header_sha256": headers.hexdigest(), "sources": sources}
    observed = {"source_sequence": before["sequence"], "mapping": before["mapping"],
                "groups": {name: _stat(stage / f"{name}.fna") if split else before["sequence"] for name in groups}}
    return groups, stats, {"size": before["sequence"]["size"], "sha256": source_hash.hexdigest()}, \
        {"size": before["mapping"]["size"], "sha256": map_hash.hexdigest()}, observed


def _verified_files(paths, task_id, receipt, folder, previous=None, *, full=False):
    """Quick local stat verification; relocation/change requires a full rehash.

    Verification states are separate from portable content identities. A caller
    can request full=True to check unchanged-stat storage corruption as well.
    The native builder always checks the SHA256 supplied in inputs.tsv itself.
    """
    previous = previous or {}
    result = {}
    records = {"mapping": receipt["context"]["mapping"], **receipt["context"]["inputs"]}
    for label, record in records.items():
        path = _source_path(paths, task_id, record["source"], folder)
        state = _stat(path)
        if state["size"] != record["size"]:
            raise RuntimeError(f"Prepared content size changed: {path}")
        observed = {"path": str(path), "state": state, "sha256": record["sha256"]}
        if full or previous.get(label) != observed:
            if sha256(path) != record["sha256"] or _stat(path) != state:
                raise RuntimeError(f"Prepared content SHA256 mismatch/change: {path}")
        result[label] = observed
    return result


def _receipt(paths, document, task, *, full=False):
    prepared = _directory(paths["base"] / "prepared")
    folder = _directory(prepared / f"task{task['task_id']:06d}")
    receipt = _sealed_json(folder / "receipt.json")
    context = receipt["context"]
    if context.get("plan_identity") != document["identity"] or context.get("task") != task:
        raise RuntimeError(f"Prepared receipt has a different plan/task: {folder}")
    previous = legacy._read_json(folder / "verification.json") if (folder / "verification.json").exists() else {}
    observed = _verified_files(paths, task["task_id"], receipt, folder, previous, full=full)
    if observed != previous:
        atomic_json(folder / "verification.json", observed)
    return receipt


def prepare_task(cfg, task_id):
    """Compute job: adopt/mask one task, seal it, and commit atomically."""
    paths = _paths(cfg)
    base = paths["base"]
    with locked(base / "reference.lock", shared=True):
        document = _load_plan(cfg)
        if isinstance(task_id, bool) or not 1 <= int(task_id) <= len(document["tasks"]):
            raise RuntimeError(f"Prepare task ID must be in 1..{len(document['tasks'])}")
        task_id = int(task_id)
        task = document["tasks"][task_id - 1]
        prepared = _directory(base / "prepared", create=True)
        final = _directory(prepared / f"task{task_id:06d}")
        with locked(prepared / f"task{task_id:06d}.lock"):
            manifest = _task_manifest(paths, task)
            if final.exists():
                receipt = _receipt(paths, document, task)
                print(f"REUSE sealed reference task {task_id}", flush=True)
                return receipt
            check_space(prepared, float(cfg.get("min_free_gib", 100)))
            owner = {"schema": SCHEMA, "kind": "parallel-reference-task", "plan_identity": document["identity"], "task_id": task_id}
            stage = prepared / f".task{task_id:06d}.partial"
            legacy._owned_stage(stage, owner)
            rows, imported = manifest["files"], manifest.get("legacy")
            if imported:
                old_id = imported["task"]["task_id"]
                old_folder = paths["legacy"] / "prepared" / f"task{old_id:06d}"
                if _small_record(old_folder / "receipt.json") != imported["receipt_record"]:
                    raise RuntimeError("Legacy prepared receipt changed after parallel planning")
                for name, record in imported["files"].items():
                    _check_old_record(old_folder / name, record)
                sequence, mapping = old_folder / "sequences.fna", old_folder / "seqid2taxid.map"
                sequence_source = {"kind": "legacy", "task_id": old_id, "name": "sequences.fna"}
                mapping_source = {"kind": "legacy", "task_id": old_id, "name": "seqid2taxid.map"}
                expected_stats = imported["stats"]
                execution = {"mode": "adopt", "mask_binary_sha256": imported["mask_binary_sha256"],
                             "historical_integrity": BASELINE_ONLY, "masking_repeated": False}
            else:
                software = legacy._software(cfg)
                if _mask_policy(cfg, software) != document["context"]["policy"]:
                    raise RuntimeError("Native preparation source/policy changed after planning")
                original_rows = [{key: value for key, value in row.items() if key != "category"} for row in rows]
                expected_stats = legacy._prepare_native({**cfg, "reference_raw": str(paths["raw"])}, original_rows, stage)
                sequence, mapping = stage / "sequences.fna", stage / "seqid2taxid.map"
                sequence_source = {"kind": "local", "name": "sequences.fna"}
                mapping_source = {"kind": "local", "name": "seqid2taxid.map"}
                execution = {"mode": "prepare", "mask_binary_sha256": software["files"].get("k2mask"),
                             "threads": int(cfg.get("prepare_threads", 8)), "historical_integrity": "sha256_at_creation"}
            split = len(task["categories"]) > 1
            groups, stats, source_record, map_record, scan_states = _scan_prepared(sequence, mapping, rows, stage, split=split,
                                                                                  masking=document["context"]["policy"]["masking"])
            if stats != expected_stats:
                raise RuntimeError("Prepared stream counters/header checksum differ from the preparation receipt")
            if imported:
                for name, observed in (("sequences.fna", source_record), ("seqid2taxid.map", map_record)):
                    previous_sha = imported["files"][name].get("sha256")
                    if previous_sha is not None and previous_sha != observed["sha256"]:
                        raise RuntimeError(f"Legacy historical SHA256 mismatch: {name}")
            for name, value in groups.items():
                value["source"] = {"kind": "local", "name": f"{name}.fna"} if split else sequence_source
            if split and not imported:
                sequence.unlink()  # Owned intermediate; complete category records have been synced.
            context = {"plan_identity": document["identity"], "task": task, "stats": stats,
                       "inputs": groups, "mapping": {**map_record, "source": mapping_source},
                       "source_sequence": source_record, "legacy_receipt_identity": imported["receipt_identity"] if imported else None}
            receipt = {"schema": SCHEMA, "status": "complete", "identity": fingerprint(context), "context": context,
                       "task_id": task_id, "execution": execution, "created_at": utcnow()}
            atomic_json(stage / "receipt.json", receipt)
            # Each byte was hashed above. Record post-close states without an
            # unnecessary second multi-terabyte read of newly adopted inputs.
            verified = {}
            for label, record in {"mapping": context["mapping"], **groups}.items():
                path = _source_path(paths, task_id, record["source"], stage)
                state = _stat(path)
                expected_state = scan_states["mapping"] if label == "mapping" else scan_states["groups"][label]
                if state != expected_state or state["size"] != record["size"]:
                    raise RuntimeError("Prepared file changed before publication")
                final_path = _source_path(paths, task_id, record["source"], final)
                verified[label] = {"path": str(final_path), "state": state, "sha256": record["sha256"]}
            atomic_json(stage / "verification.json", verified)
            legacy._sync_dir(stage)
            stage.rename(final)
            legacy._sync_dir(prepared)
            return receipt


def _copy_hash(source, target):
    before = _stat(source)
    digest = hashlib.sha256()
    with Path(source).open("rb", buffering=8 << 20) as src, Path(target).open("wb", buffering=8 << 20) as dst:
        for block in iter(lambda: src.read(8 << 20), b""):
            digest.update(block)
            dst.write(block)
        legacy._sync_file(dst)
    if before != _stat(source):
        raise RuntimeError(f"Source changed while freezing metadata: {source}")
    return {"size": before["size"], "sha256": digest.hexdigest()}


def finalize_preparation(cfg):
    """Compute job after all tasks: publish global map/taxonomy and input TSVs.

    Content identity excludes deployment paths. Regenerating input TSVs after
    copying a completed cache changes only paths, never reference identities.
    """
    paths = _paths(cfg)
    base = paths["base"]
    with locked(base / "reference.lock", shared=True), locked(base / "finalize.lock"):
        document = _load_plan(cfg)
        receipts, groups, identifiers = [], {}, []
        for task in document["tasks"]:
            with locked(base / "prepared" / f"task{task['task_id']:06d}.lock"):
                _task_manifest(paths, task)
                receipt = _receipt(paths, document, task)
            receipts.append((task, receipt))
            identifiers.append(receipt["identity"])
            for name, item in receipt["context"]["inputs"].items():
                group = groups.setdefault(name, {"inputs": [], "records": 0, "bases": 0, "longest_record": 0})
                source = _source_path(paths, task["task_id"], item["source"])
                if any(c in str(source) for c in ("\t", "\r", "\n")):
                    raise RuntimeError("Native input TSV paths must not contain tab/newline characters")
                group["inputs"].append({"path": str(source), **{k: item[k] for k in ("size", "sha256", "records", "bases")}})
                group["records"] += item["records"]
                group["bases"] += item["bases"]
                group["longest_record"] = max(group["longest_record"], item["longest_record"])
        source_id = fingerprint({"plan_identity": document["identity"], "prepared_identities": identifiers})
        destination = _directory(base / "global")
        if destination.exists():
            frozen = _sealed_json(destination / "receipt.json")
            if frozen["context"].get("source_identity") != source_id:
                raise RuntimeError("Global reference snapshot belongs to different prepared inputs")
            # These metadata files are small compared with references. Rehash
            # them at every finalize, including copied deployments.
            for relative, record in frozen["context"]["files"].items():
                if _small_record(_safe_child(destination, relative)) != record:
                    raise RuntimeError(f"Global frozen mapping/taxonomy SHA256 changed: {relative}")
        else:
            check_space(base, float(cfg.get("min_free_gib", 100)))
            owner = {"schema": SCHEMA, "kind": "parallel-global-reference", "source_identity": source_id}
            stage = base / ".global.partial"
            legacy._owned_stage(stage, owner)
            (stage / "taxonomy").mkdir()
            files = {}
            for name in legacy.TAX_FILES + ("delnodes.dmp",):
                source = paths["metadata"] / "taxonomy" / name
                if source.is_file():
                    files["taxonomy/" + name] = _copy_hash(source, stage / "taxonomy" / name)
            # Taxids must still resolve to exactly the IDs stored in the raw
            # manifests. This catches changed/inconsistent taxonomy metadata.
            resolve = legacy._taxonomy_ids(stage / "taxonomy")
            taxids = set()
            for task in document["tasks"]:
                taxids.update(row["taxid"] for row in _task_manifest(paths, task)["files"])
            for taxid in taxids:
                if resolve(taxid) != taxid:
                    raise RuntimeError(f"Prepared taxid {taxid} changed in the frozen taxonomy")
            map_hash, map_size = hashlib.sha256(), 0
            with (stage / "seqid2taxid.map").open("wb", buffering=8 << 20) as out:
                for task, receipt in receipts:
                    record = receipt["context"]["mapping"]
                    source = _source_path(paths, task["task_id"], record["source"])
                    before, digest, length = _stat(source), hashlib.sha256(), 0
                    with source.open("rb", buffering=8 << 20) as handle:
                        for block in iter(lambda: handle.read(8 << 20), b""):
                            digest.update(block)
                            map_hash.update(block)
                            length += len(block)
                            out.write(block)
                    if digest.hexdigest() != record["sha256"] or length != record["size"] or before != _stat(source):
                        raise RuntimeError("Prepared map changed during global concatenation")
                    map_size += length
                legacy._sync_file(out)
            files["seqid2taxid.map"] = {"size": map_size, "sha256": map_hash.hexdigest()}
            context = {"source_identity": source_id, "files": files}
            frozen = {"schema": SCHEMA, "status": "complete", "context": context,
                      "identity": fingerprint(context), "created_at": utcnow()}
            atomic_json(stage / "receipt.json", frozen)
            legacy._sync_dir(stage)
            stage.rename(destination)
            legacy._sync_dir(base)
        category_root = _directory(base / "categories", create=True)
        for name in document["categories"]:
            group = groups[name]
            target = category_root / f"{name}.inputs.tsv"
            temporary = category_root / f".{name}.inputs.tsv.partial"
            if target.is_symlink() or temporary.is_symlink():
                raise RuntimeError("Refusing a symlink native input list")
            with temporary.open("w", encoding="utf-8", newline="\n") as handle:
                for item in group["inputs"]:
                    handle.write("\t".join(str(item[key]) for key in ("path", "size", "sha256", "records", "bases")) + "\n")
                legacy._sync_file(handle)
            temporary.replace(target)
            group["inputs_tsv"] = str(target)
        legacy._sync_dir(category_root)
        portable_groups = {name: {**{key: value for key, value in group.items() if key not in ("inputs", "inputs_tsv")},
                                 "inputs": [{k: v for k, v in item.items() if k != "path"} for item in group["inputs"]]}
                           for name, group in groups.items()}
        context = {"plan_identity": document["identity"], "source_identity": source_id,
                   "global_identity": frozen["identity"], "global_files": frozen["context"]["files"],
                   "groups": portable_groups}
        catalog = {"schema": CATALOG_SCHEMA, "status": "complete", "context": context,
                   "identity": fingerprint(context), "groups": groups,
                   "taxonomy_dir": str(destination / "taxonomy"), "seqid_map": str(destination / "seqid2taxid.map"),
                   "created_at": utcnow()}
        atomic_json(base / "catalog.json", catalog)
        return catalog


def catalog_snapshot(cfg, *, full=False):
    """Validate a completed catalog; ``full`` additionally rehashes all inputs.

    This is a compute-side integrity operation. It never creates a missing
    prepared task, changes a reference, or trusts a stale absolute deployment
    path. After relocating a cache run finalize_preparation once to regenerate
    native TSV paths; catalog content identity remains unchanged.
    """
    paths = _paths(cfg)
    base = paths["base"]
    with locked(base / "reference.lock", shared=True), locked(base / "finalize.lock"):
        document = _load_plan(cfg)
        catalog = _sealed_json(base / "catalog.json", CATALOG_SCHEMA)
        if catalog["context"].get("plan_identity") != document["identity"]:
            raise RuntimeError("Catalog belongs to a different reference plan")
        expected_taxonomy = str(base / "global" / "taxonomy")
        expected_mapping = str(base / "global" / "seqid2taxid.map")
        if catalog.get("taxonomy_dir") != expected_taxonomy or catalog.get("seqid_map") != expected_mapping:
            raise RuntimeError("Reference catalog was relocated; run finalize_preparation to refresh paths")
        receipts, actual_groups = [], {}
        for task in document["tasks"]:
            with locked(base / "prepared" / f"task{task['task_id']:06d}.lock"):
                _task_manifest(paths, task)
                receipt = _receipt(paths, document, task, full=full)
            receipts.append(receipt["identity"])
            for name, item in receipt["context"]["inputs"].items():
                group = actual_groups.setdefault(name, {"inputs": [], "records": 0, "bases": 0, "longest_record": 0})
                group["inputs"].append({"path": str(_source_path(paths, task["task_id"], item["source"])),
                                        **{k: item[k] for k in ("size", "sha256", "records", "bases")}})
                group["records"] += item["records"]
                group["bases"] += item["bases"]
                group["longest_record"] = max(group["longest_record"], item["longest_record"])
        source_identity = fingerprint({"plan_identity": document["identity"], "prepared_identities": receipts})
        if source_identity != catalog["context"].get("source_identity"):
            raise RuntimeError("Catalog prepared input identities changed")
        for name, group in actual_groups.items():
            path = base / "categories" / f"{name}.inputs.tsv"
            group["inputs_tsv"] = str(path)
            expected = "".join("\t".join(str(row[k]) for k in ("path", "size", "sha256", "records", "bases")) + "\n"
                               for row in group["inputs"])
            if path.is_symlink() or path.read_text(encoding="utf-8") != expected:
                raise RuntimeError(f"Native category input list changed: {path}")
        if actual_groups != catalog.get("groups"):
            raise RuntimeError("Catalog group paths or counters differ from its prepared receipts; finalize after relocation")
        global_root = base / "global"
        frozen = _sealed_json(global_root / "receipt.json")
        if frozen["identity"] != catalog["context"].get("global_identity") or \
                frozen["context"].get("files") != catalog["context"].get("global_files"):
            raise RuntimeError("Catalog global mapping/taxonomy identity differs")
        previous_path = global_root / "verification.json"
        previous = legacy._read_json(previous_path) if previous_path.exists() else {}
        observations = {}
        for relative, record in frozen["context"]["files"].items():
            path = _safe_child(global_root, relative)
            state = _stat(path)
            observed = {"path": str(path), "state": state, "sha256": record["sha256"]}
            if state["size"] != record["size"]:
                raise RuntimeError(f"Global reference file size changed: {path}")
            if full or previous.get(relative) != observed:
                if sha256(path) != record["sha256"] or _stat(path) != state:
                    raise RuntimeError(f"Global reference SHA256 mismatch/change: {path}")
            observations[relative] = observed
        if observations != previous:
            atomic_json(previous_path, observations)
        return catalog
