#!/usr/bin/env python3
"""Build a capped, competing all-group Kraken index from existing prepared FASTA.

The five full indexes are reused by stage 2 and are never rewritten here. This
new single index retains a deterministic hash-selected fraction of minimizers
from ALL prepared references. A smaller hash is a sensitivity tradeoff, not an
equivalent representation of the full library, nor a hard bound on total RAM.
No prepared FASTA, mapping, or large index is copied or rehashed.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import shutil
import struct
import subprocess
import sys

from common import atomic_json, check_space, fingerprint, load_config, locked, sha256, utcnow
from ordinal_query import NativeTelemetry, stop_process
from reference import GIB, PARAMETERS, SOURCE_COMMIT, _interruptible, _native_index_structure, _software
from result_output import replace_bool
from storage_layout import category_jobs, library_jobs, library_root

SCHEMA = "metag-compact-group-library-v1"
DEFAULT_MAX_GIB = 128
NATIVE_FILES = ("hash.k2d", "taxo.k2d", "opts.k2d")
CAUTION = ("Hash-based minimizer downsampling reduces sensitivity, especially for short, "
           "divergent, masked or rare-target fragments. The retention fraction is a hash "
           "selection fraction, not a measured read recall. All groups compete together; "
           "a group call is provisional and subsequent species classification is conditional "
           "on routing. Validate retained targets and quantify missed-reference controls.")


def _read(path):
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Cannot read metadata {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise RuntimeError(f"Expected a metadata object: {path}")
    return data


def _absolute(cfg, path):
    path = Path(path).expanduser()
    return (path if path.is_absolute() else Path(cfg.get("root", ".")) / path).absolute()


def group_directory(cfg):
    return _absolute(cfg, cfg.get("group_library_dir") or library_root(cfg) / "group")


def native_executable(cfg):
    path = Path(cfg["tools_dir"]) / "kraken" / "bin" / "classify"
    if not path.is_file() or not os.access(path, os.X_OK):
        raise RuntimeError(f"Missing native Kraken classifier: {path}; run setup")
    return path.resolve()


def _observation(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Expected a regular non-symlink file: {path}")
    value = path.stat()
    return {"size": value.st_size, "mtime_ns": value.st_mtime_ns,
            "ctime_ns": value.st_ctime_ns, "inode": value.st_ino}


def _check_observations(items):
    for item in items:
        if _observation(item["path"]) != item["observation"]:
            raise RuntimeError(f"Input changed during the group build: {item['path']}")


def source_catalog(cfg):
    """Read the sealed small catalog and stat inputs; never rehash their payloads."""
    base = _absolute(cfg, library_root(cfg))
    from library_paths import prepared_catalog
    catalog = prepared_catalog(base)
    if cfg.get('source_preparation_identity',catalog['identity']) != catalog['identity']:
        raise RuntimeError('Frozen bundle preparation changed')
    context = catalog.get("context", {})
    from partitioned_reference import CATALOG_SCHEMA
    if catalog.get("schema") != CATALOG_SCHEMA or \
            catalog.get("status") != "complete" or fingerprint(context) != catalog.get("identity"):
        raise RuntimeError("Complete partitioned prepared catalog is required; finalize the existing library first")
    groups = catalog.get("groups", {})
    portable = {name: {**{k: v for k, v in group.items() if k not in ("inputs", "inputs_tsv")},
                      "inputs": [{k: v for k, v in row.items() if k != "path"}
                                 for row in group.get("inputs", [])]}
                for name, group in groups.items()}
    if not groups or portable != context.get("groups"):
        raise RuntimeError("Prepared groups do not match the sealed catalog")
    sources, seen = [], set()
    for category, group in sorted(groups.items()):
        for row in group["inputs"]:
            path = Path(row["path"])
            if not path.is_absolute():
                raise RuntimeError(f"Prepared catalog needs finalized absolute paths: {path}")
            canonical = path.resolve(strict=True)
            if canonical in seen:
                raise RuntimeError(f"Prepared FASTA is duplicated across catalog inputs: {path}")
            seen.add(canonical)
            observed = _observation(path)
            if observed["size"] != int(row["size"]) or observed["size"] <= 0:
                raise RuntimeError(f"Prepared FASTA size differs from catalog: {path}")
            sources.append({"path": str(path), "category": category, "observation": observed,
                            "records": int(row["records"]), "bases": int(row["bases"])})
    taxonomy = Path(catalog["taxonomy_dir"])
    mapping = Path(catalog["seqid_map"])
    if taxonomy != base / "global" / "taxonomy" or mapping != base / "global" / "seqid2taxid.map":
        raise RuntimeError("Prepared catalog paths were relocated; finalize paths before building the group index")
    metadata = []
    global_files = context.get("global_files", {})
    if not {"taxonomy/nodes.dmp", "taxonomy/names.dmp", "seqid2taxid.map"} <= set(global_files):
        raise RuntimeError("Prepared catalog lacks shared taxonomy or sequence mapping records")
    for relative, record in sorted(global_files.items()):
        parts = PurePosixPath(relative)
        if parts.is_absolute() or ".." in parts.parts:
            raise RuntimeError("Unsafe relative path in the global metadata catalog")
        path = base / "global" / relative
        observed = _observation(path)
        if observed["size"] != int(record["size"]):
            raise RuntimeError(f"Frozen global metadata size differs: {path}")
        metadata.append({"path": str(path), "observation": observed})
    if not sources:
        raise RuntimeError("No prepared FASTA inputs in catalog")
    return catalog, sources, metadata


def _reference_table(cfg, catalog, destination):
    """Export raw-reference metadata using the very same frozen taxonomy."""
    from taxonomy import Taxonomy
    base, raw = _absolute(cfg, library_root(cfg)), _absolute(cfg, cfg["reference_raw"])
    plan = _read(base / "plan.json")
    if plan.get("status") != "complete" or plan.get("identity") != catalog["context"].get("plan_identity") or \
            plan.get("identity") != fingerprint(plan.get("context")):
        raise RuntimeError("Reference plan does not match the prepared catalog")
    tasks = plan.get("tasks") or plan.get("context", {}).get("tasks")
    if not isinstance(tasks, list) or not tasks or tasks != plan["context"].get("tasks"):
        raise RuntimeError("Reference plan has no tasks for references.tsv")
    taxonomy = Taxonomy(catalog["taxonomy_dir"])
    fields = ("path", "relative_path", "source_taxid", "taxid", "species_taxid", "name",
              "group", "category", "accession", "records", "bases", "compressed_bytes", "task_id")
    seen, count, receipt_ids = set(), 0, []
    with Path(destination).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        for task in tasks:
            task_id = int(task["task_id"])
            relative_manifest = PurePosixPath(task["manifest"]["path"])
            if relative_manifest.is_absolute() or ".." in relative_manifest.parts:
                raise RuntimeError("Unsafe reference manifest path")
            manifest_path = base / relative_manifest
            if manifest_path.stat().st_size != task["manifest"]["size"] or \
                    sha256(manifest_path) != task["manifest"]["sha256"]:
                raise RuntimeError(f"Frozen reference manifest changed: {manifest_path}")
            manifest = _read(manifest_path)
            receipt = _read(base / "prepared" / f"task{task_id:06d}" / "receipt.json")
            if manifest.get("task_id") != task_id or receipt.get("status") != "complete" or \
                    receipt.get("identity") != fingerprint(receipt.get("context")):
                raise RuntimeError(f"Invalid prepared metadata receipt for task {task_id}")
            receipt_ids.append(receipt["identity"])
            rows = manifest.get("files", [])
            counters = receipt.get("context", {}).get("stats", {}).get("sources", [])
            if not rows or len(rows) != len(counters):
                raise RuntimeError(f"Prepared source counters differ from manifest for task {task_id}")
            for row, counter in zip(rows, counters):
                relative = PurePosixPath(row["path"])
                if relative.is_absolute() or ".." in relative.parts or str(relative) in seen or \
                        counter.get("path") != str(relative):
                    raise RuntimeError(f"Invalid/duplicate source in reference plan: {relative}")
                seen.add(str(relative))
                taxid = taxonomy.resolve(int(row["taxid"]))
                info = taxonomy.describe(taxid)
                writer.writerow({"path": str(raw / relative), "relative_path": str(relative),
                                 "source_taxid": row["taxid"], "taxid": taxid,
                                 "species_taxid": info["species_taxid"], "name": info["name"],
                                 "group": info["group"], "category": row.get("category", ""),
                                 "accession": row.get("accession", ""), "records": counter["records"],
                                 "bases": counter["bases"], "compressed_bytes": row.get("state", {}).get("size", 0),
                                 "task_id": task_id})
                count += 1
    if fingerprint({"plan_identity": plan["identity"], "prepared_identities": receipt_ids}) != \
            catalog["context"].get("source_identity"):
        raise RuntimeError("Prepared receipt identities differ from the catalog")
    return count


def _run_stream(command, sources, log_path, stdout_path):
    """Pipe existing bytes once; do not create a concatenated FASTA or checksum."""
    process, monitor = None, None
    with Path(log_path).open("w", encoding="utf-8") as log, Path(stdout_path).open("wb") as out, _interruptible():
        log.write("COMMAND " + shlex.join(str(x) for x in command) + "\n")
        log.flush()
        try:
            process = subprocess.Popen([str(x) for x in command], stdin=subprocess.PIPE,
                                       stdout=out, stderr=log, start_new_session=True)
            monitor = NativeTelemetry(process.pid)
            for item in sources:
                _check_observations([item])
                with Path(item["path"]).open("rb", buffering=8 << 20) as handle:
                    for block in iter(lambda: handle.read(8 << 20), b""):
                        process.stdin.write(block)
                # A FASTA may lack its final newline; separating files must
                # never append the next header to the previous sequence line.
                process.stdin.write(b"\n")
                _check_observations([item])
            process.stdin.close()
            status = process.wait()
            if status:
                raise RuntimeError(f"Native group command exited {status}; see {log_path}")
        except BaseException as exc:
            stop_process(process)
            if isinstance(exc, BrokenPipeError):
                raise RuntimeError(f"Native group command stopped reading; see {log_path}") from exc
            raise
        finally:
            if process is not None and process.stdin is not None and not process.stdin.closed:
                try:
                    process.stdin.close()
                except BrokenPipeError:
                    pass
            if monitor is not None:
                log.write("NATIVE_METRICS " + json.dumps(monitor.finish(), sort_keys=True) + "\n")


def _native_options(path):
    data = Path(path).read_bytes()
    if len(data) != 64 or sys.byteorder != "little" or struct.calcsize("P") != 8:
        raise RuntimeError("Unsupported/truncated pinned native options ABI (requires 64-bit little-endian Linux)")
    k, length, seed, toggle = struct.unpack_from("<QQQQ", data)
    minimum = struct.unpack_from("<Q", data, 40)[0]
    if k != PARAMETERS["k"] or length != PARAMETERS["l"] or not data[32]:
        raise RuntimeError("Group index does not have the requested DNA k/l options")
    return {"k": k, "l": length, "spaced_seed_mask": seed, "toggle_mask": toggle,
            "minimum_acceptable_hash_value": minimum,
            "selected_hash_fraction": ((1 << 64) - minimum) / (1 << 64)}


def _saved_estimate(cfg, catalog):
    estimates, receipts = [], []
    for category in sorted(catalog["groups"]):
        path = category_jobs(cfg, category) / "capacity.json"
        row = _read(path)
        expected = {"build_identity": row.get("build_identity"), "estimate": row.get("estimate")}
        if row.get("identity") != fingerprint(expected) or int(row.get("estimate", 0)) <= 0:
            raise RuntimeError(f"Invalid saved category capacity estimate: {path}")
        category_marker = _read(_absolute(cfg, library_root(cfg)) / "index" / category / "database.json")
        if category_marker.get("status") != "complete" or \
                category_marker.get("identity") != fingerprint(category_marker.get("context")) or \
                category_marker.get("context", {}).get("preparation_identity") != catalog["identity"] or \
                category_marker.get("context", {}).get("build_identity") != row.get("build_identity"):
            raise RuntimeError(f"Saved category estimate is not tied to the current prepared catalog: {category}")
        estimates.append(int(row["estimate"]))
        receipts.append({"category": category, "path": str(path), "estimate": int(row["estimate"]),
                         "identity": row["identity"]})
    return sum(estimates), {"method": "sum_saved_category_estimates", "receipts": receipts,
                            "caution": "Shared minimizers are counted in multiple categories; this can overestimate "
                                       "union uniqueness and retain fewer minimizers than native union estimation."}


def _owned(path):
    path = Path(path)
    return not path.is_symlink() and path.is_dir() and \
        _read(path / "owner.json").get("schema") == SCHEMA


def compact_snapshot(cfg, directory=None):
    folder = _absolute(cfg, directory) if directory is not None else group_directory(cfg)
    marker = _read(folder / "database.json")
    context = marker.get("context", {})
    if marker.get("schema") != SCHEMA or marker.get("status") != "complete" or \
            marker.get("layout") != "single" or \
            marker.get("identity") != fingerprint(context):
        raise RuntimeError(f"Incomplete or invalid compact group index: {folder}")
    for key in ("files", "taxonomy_files", "taxonomy_dir", "downsampling", "parameters", "preparation_identity"):
        if marker.get(key) != context.get(key):
            raise RuntimeError(f"Compact snapshot {key} differs from its completion identity")
    # Require the current library paths recorded at publication.
    if Path(marker['taxonomy_dir']).parts[-2:] != ('global','taxonomy'):
        raise RuntimeError('Compact catalog taxonomy is outside the prepared global layout')
    # Named indexes can be nested below group/; taxonomy belongs to the common
    # library root, not necessarily to the index directory's immediate parent.
    if Path(marker['taxonomy_dir']) != _absolute(cfg, library_root(cfg)) / 'global/taxonomy' or Path(marker['db_dir']) != folder:
        raise RuntimeError('Index paths differ from the current publication; rebuild in the intended library')
    if set(marker["files"]) != set(NATIVE_FILES) or \
            not {"nodes.dmp", "names.dmp"} <= set(marker["taxonomy_files"]):
        raise RuntimeError("Compact snapshot lacks native files or frozen taxonomy")
    for base, table in ((folder, marker["files"]), (Path(marker["taxonomy_dir"]), marker["taxonomy_files"])):
        for name, record in table.items():
            if Path(name).name != name or _observation(base / name) != record["observation"]:
                raise RuntimeError(f"Compact index file or taxonomy changed: {base / name}")
    actual = _native_options(folder / "opts.k2d")
    if actual != marker["downsampling"]["native_options"]:
        raise RuntimeError("Compact native downsampling options changed")
    _native_index_structure(folder, expected_capacity=context["capacity"])
    refs = folder / "references.tsv"
    if _observation(refs) != context["references"]["observation"]:
        raise RuntimeError("Compact reference metadata table changed")
    return marker


def reusable_union_estimate(cfg, catalog, software, parameters):
    base = _absolute(cfg, library_root(cfg))
    candidates = [Path(cfg['group_estimate_source'])] if cfg.get('group_estimate_source') else []
    candidates += [base/'group/optimized']
    reasons=[]
    for directory in dict.fromkeys(candidates):
        if not (directory/'database.json').is_file():continue
        try:
            prior=compact_snapshot(cfg,directory)
            if (prior['preparation_identity'] != catalog['identity'] or prior['parameters'] != parameters
                    or prior['context']['software'] != software):
                raise RuntimeError('inputs, parameters or native program differ')
            info=prior['downsampling']
            if info['estimate_basis'].get('method')!='native_union_estimate_capacity':
                raise RuntimeError('estimate was not computed over the native union')
            return int(info['distinct_minimizers_estimate']),dict(info['estimate_basis'],reused_from=str(directory/'database.json'))
        except (OSError,ValueError,KeyError,RuntimeError) as exc:
            reasons.append(str(directory)+': '+str(exc))
    raise RuntimeError('No matching saved native union estimate'+('; '+'; '.join(reasons) if reasons else ''))


def build_group(cfg, out=None, threads=None, max_gib=None, estimate=None, replace=False):
    folder = _absolute(cfg, out) if out is not None else group_directory(cfg)
    threads = int(threads if threads is not None else cfg.get("group_build_threads", cfg.get("build_threads", 8)))
    max_gib = float(max_gib if max_gib is not None else cfg.get("group_max_gib", DEFAULT_MAX_GIB))
    method = estimate or cfg.get("group_estimate", "native")
    if threads < 1 or not math.isfinite(max_gib) or max_gib <= 0 or method not in ("auto", "native", "saved", "reuse"):
        raise ValueError("Invalid group build threads, hash cap, or estimate mode")
    max_bytes = int(max_gib * GIB)
    max_capacity = (max_bytes - 32) // 4
    if max_capacity < 1024:
        raise ValueError("Group hash cap must allow at least 1024 cells plus its 32-byte header")
    base = _absolute(cfg, library_root(cfg))
    protected = [base / name for name in ("prepared", "index", "global", "categories", "manifests")]
    protected += [_absolute(cfg, cfg[name]) for name in ("reference_raw", "metadata_dir", "tools_dir")]
    if folder == base or any(folder == p or p in folder.parents for p in protected):
        raise RuntimeError("Group output must be separate from existing references, indexes, metadata and tools")
    folder.parent.mkdir(parents=True, exist_ok=True)
    stage = folder.with_name("." + folder.name + ".building")
    backup = folder.with_name("." + folder.name + ".previous")
    lock = folder.with_name("." + folder.name + ".lock")
    with locked(library_jobs(cfg) / "reference.lock", shared=True), locked(lock):
        if folder.is_symlink() or stage.is_symlink() or backup.is_symlink():
            raise RuntimeError("Compact output/staging directories must not be symlinks")
        if backup.exists():
            if not _owned(backup):
                raise RuntimeError(f"Unowned interrupted publication backup: {backup}")
            if not folder.exists():
                backup.rename(folder)
            elif _owned(folder) and (folder / "database.json").is_file():
                compact_snapshot(cfg, folder)
                shutil.rmtree(backup)
            else:
                raise RuntimeError("Interrupted compact publication requires an intact completed destination")
        catalog, sources, metadata = source_catalog(cfg)
        software = _software(cfg)
        build_context = {"schema": SCHEMA, "source_selection": "all_prepared_sequences",
                         "preparation_identity": catalog["identity"], "sources": sources,
                         "global_metadata": metadata, "max_hash_bytes": max_bytes,
                         "threads": threads, "estimate_method": method, "software": software,
                         "parameters": {k: PARAMETERS[k] for k in ("k", "l", "seed", "load_factor", "block_size")}}
        build_id = fingerprint(build_context)
        if folder.exists():
            if not _owned(folder):
                raise RuntimeError(f"Refusing to replace an unowned group directory: {folder}")
            try:
                complete = compact_snapshot(cfg, folder)
            except (RuntimeError, OSError):
                if not replace:
                    raise
                complete = None
            if complete is not None and complete["context"]["build_identity"] == build_id:
                print(f"Reusing completed compact group index: {folder}", flush=True)
                return complete
            if not replace:
                raise RuntimeError("Group build settings/inputs changed; use --replace or a different --out")
        if stage.exists():
            if not _owned(stage):
                raise RuntimeError(f"Refusing an unowned group staging directory: {stage}")
            old_plan = _read(stage / "plan.json")
            if old_plan.get("identity") != build_id:
                if not replace:
                    raise RuntimeError("Interrupted group build has different inputs/settings; use --replace")
                shutil.rmtree(stage)
        stage.mkdir(exist_ok=True)
        atomic_json(stage / "owner.json", {"schema": SCHEMA, "destination": str(folder)})
        atomic_json(stage / "plan.json", {"status": "planned", "identity": build_id,
                                         "context": build_context, "created_at": utcnow(), "caution": CAUTION})
        binary = Path(cfg["tools_dir"]) / "kraken" / "bin"
        flags = ["-k", "35", "-l", "31", "-S", PARAMETERS["seed"], "-p", str(threads)]
        capacity_file = stage / "capacity.json"
        if capacity_file.is_file():
            saved = _read(capacity_file)
            if saved.get("build_identity") != build_id or saved.get("identity") != fingerprint(saved.get("context")):
                raise RuntimeError("Interrupted union capacity estimate does not match this build")
            estimate_value, basis = int(saved["context"]["estimate"]), saved["context"]["basis"]
        else:
            reused = None
            if method in ('auto','reuse'):
                try:reused=reusable_union_estimate(cfg,catalog,software,build_context['parameters'])
                except RuntimeError:
                    if method=='reuse':raise
                    print('No reusable union estimate; estimating from prepared FASTA',flush=True)
            if reused is not None:
                estimate_value,basis=reused
                print('Reusing the matching native union estimate; no estimation FASTA pass',flush=True)
            elif method in ('auto','native'):
                print(f"Estimating compact union from all {len(sources)} prepared FASTA files; one sequential pass", flush=True)
                _run_stream([binary / "estimate_capacity", *flags], sources,
                            stage / "estimate.log", stage / "estimate.stdout")
                text = (stage / "estimate.stdout").read_text().strip()
                if not text.isdigit() or int(text) <= 0:
                    raise RuntimeError(f"Invalid native union estimate: {text!r}")
                estimate_value = int(text)
                basis = {"method": "native_union_estimate_capacity", "sampled_hash_sections": 4,
                         "total_hash_sections": 1024, "note": "Native estimate from all sequences; not a sum of overlapping category estimates."}
            else:
                estimate_value, basis = _saved_estimate(cfg, catalog)
            _check_observations(sources + metadata)
            context = {"estimate": estimate_value, "basis": basis}
            atomic_json(capacity_file, {"build_identity": build_id, "context": context,
                                       "identity": fingerprint(context), "created_at": utcnow()})
        capacity = (estimate_value + 8192) * 10 // 7
        actual_capacity = min(capacity, max_capacity)
        hash_bytes = 32 + 4 * actual_capacity
        longest = max(int(group.get("longest_record", 0)) for group in catalog["groups"].values())
        resource = {"capacity_before_cap": capacity, "capacity": actual_capacity,
                    "hash_bytes": hash_bytes, "max_hash_bytes": max_bytes,
                    "estimated_build_ram_bytes": math.ceil(hash_bytes * 1.15) + 8 * GIB + 2 * longest,
                    "note": "Hash capacity is bounded; total process RAM is not. Native estimation uses sampled sets; "
                            "building also needs global sequence mapping, taxonomy and complete-record buffers."}
        atomic_json(stage / "resource_estimate.json", resource)
        budget = float(cfg.get("group_build_memory_gib", 192)) * GIB
        if not math.isfinite(budget) or resource["estimated_build_ram_bytes"] > budget:
            raise RuntimeError("Estimated group build RAM exceeds group_build_memory_gib; "
                               "increase the job memory or lower group_max_gib. The estimate is saved.")
        check_space(stage, float(cfg.get("min_free_gib", 0)) + hash_bytes * 1.01 / GIB)
        print(f"Building group hash {hash_bytes / GIB:.3f} GiB; selected capacity fraction "
              f"{actual_capacity / capacity:.6f}. {CAUTION}", flush=True)
        references = _reference_table(cfg, catalog, stage / "references.tsv")
        _check_observations(sources + metadata)
        for name in NATIVE_FILES:
            path = stage / name
            if path.is_symlink():
                raise RuntimeError(f"Unsafe native staging file: {path}")
            path.unlink(missing_ok=True)
        command = [binary / "build_db", *flags, "-H", stage / "hash.k2d", "-t", stage / "taxo.k2d",
                   "-o", stage / "opts.k2d", "-n", catalog["taxonomy_dir"], "-m", catalog["seqid_map"],
                   "-c", str(capacity), "-B", "16384", "-b", str(math.ceil(16384 / threads)), "-r", "0"]
        if actual_capacity < capacity:
            command += ["-M", str(actual_capacity)]
        _run_stream(command, sources, stage / "build.log", stage / "build.stdout")
        _check_observations(sources + metadata)
        counters = re.findall(r"Completed processing of (\d+) sequences, (\d+) bp",
                              (stage / "build.log").read_text(errors="replace"))
        expected_counts = (sum(row["records"] for row in sources), sum(row["bases"] for row in sources))
        if len(counters) != 1 or tuple(map(int, counters[0])) != expected_counts:
            raise RuntimeError("Native processed sequence/base counts differ from all prepared catalog inputs; "
                               "no group index was published. Inspect build.log and the frozen mapping.")
        options = _native_options(stage / "opts.k2d")
        expected_fraction = actual_capacity / capacity
        if not math.isclose(options["selected_hash_fraction"], expected_fraction, rel_tol=0, abs_tol=1e-14):
            raise RuntimeError("Native hash downsampling threshold differs from the requested cap")
        structure = _native_index_structure(stage, expected_capacity=actual_capacity)
        def record(path):
            obs = _observation(path)
            return {"state": {"size": obs["size"]}, "observation": obs}
        files = {name: record(stage / name) for name in NATIVE_FILES}
        taxonomy_files = {Path(row["path"]).name: record(row["path"]) for row in metadata
                          if Path(row["path"]).parent == Path(catalog["taxonomy_dir"])}
        downsampling = {"enabled": actual_capacity < capacity, "source_selection": "all_prepared_sequences",
                        "capacity_fraction": expected_fraction, "native_options": options,
                        "distinct_minimizers_estimate": estimate_value, "estimate_basis": basis,
                        "caution": CAUTION, "measured_read_recall": None}
        context = {"build_identity": build_id, "preparation_identity": catalog["identity"],
                   "capacity": actual_capacity, "resource_estimate": resource,
                   "files": files, "taxonomy_dir": catalog["taxonomy_dir"], "taxonomy_files": taxonomy_files,
                   "downsampling": downsampling, "parameters": build_context["parameters"],
                   "native_structure": structure, "software": software,
                   "processed_records": expected_counts[0], "processed_bases": expected_counts[1],
                   "references": {"rows": references, "observation": _observation(stage / "references.tsv"),
                                  "sha256": sha256(stage / "references.tsv")},
                   "input_verification": "Existing catalog identity and portable size/mtime/ctime/inode checks; no new FASTA/index SHA256 pass."}
        marker = {"schema": SCHEMA, "status": "complete", "layout": "single", "db_dir": str(folder),
                  "context": context, "identity": fingerprint(context), "created_at": utcnow(),
                  **{key: context[key] for key in ("files", "taxonomy_dir", "taxonomy_files", "downsampling",
                                                  "parameters", "preparation_identity")}}
        for name in NATIVE_FILES:
            with (stage / name).open("rb") as handle:
                os.fsync(handle.fileno())
        atomic_json(stage / "database.json", marker)
        moved_old = False
        try:
            if folder.exists():
                folder.rename(backup)
                moved_old = True
            stage.rename(folder)
        except BaseException:
            if moved_old and not folder.exists():
                backup.rename(folder)
            raise
        complete = compact_snapshot(cfg, folder)
        if backup.exists():
            shutil.rmtree(backup)
        print(f"Completed compact group index: {folder}", flush=True)
        return complete


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config")
    parser.add_argument("--root")
    parser.add_argument("--out")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--max-gib", type=float, help="Maximum hash.k2d GiB (default config or 128); not total RAM")
    parser.add_argument("--estimate", choices=("auto", "native", "saved", "reuse"), help="reuse validates the existing union estimate; native scans prepared references; saved uses category estimates")
    parser.add_argument("--replace", type=replace_bool, default=False, metavar="TRUE|FALSE")
    parser.add_argument("--allow-local", action="store_true")
    args = parser.parse_args(argv)
    from common import require_compute
    require_compute(args.allow_local)
    cfg = load_config(args.config, args.root)
    build_group(cfg, out=args.out, threads=args.threads, max_gib=args.max_gib,
                estimate=args.estimate, replace=args.replace)


if __name__ == "__main__":
    main()
