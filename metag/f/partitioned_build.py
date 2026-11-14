"""Five independent resumable hashes, one per-minimizer union classifier.

The original installed Kraken binaries are never patched in place. The private
builder persists its mutable native table. Completed category files are only
published after structural validation and the native payload SHA receipt.
"""
from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from common import atomic_json, check_space, file_state, fingerprint, locked, sha256, utcnow
from reference import PARAMETERS, _native_index_structure, _native_opts, _run_feed

from storage_layout import library_root, library_jobs, category_jobs, category_logs, separated

SCHEMA = "metag-partitioned-index-v1"
CATEGORIES = ("microbes", "viruses", "plants", "invertebrates", "vertebrates")
GIB = 2 ** 30


def root(cfg):
    value = library_root(cfg).expanduser()
    return (value if value.is_absolute() else Path(cfg["root"]) / value).resolve()


def index_root(cfg):
    return root(cfg) / "index"


def runtime_state(path):
    value = Path(path).stat()
    return {"size": value.st_size, "mtime_ns": value.st_mtime_ns, "ctime_ns": value.st_ctime_ns}


def record(path, digest=None):
    return {"state": file_state(path), "runtime_state": runtime_state(path), "sha256": digest or sha256(path)}


def check_record(path, item, *, full=False):
    path = Path(path)
    if not path.is_file() or path.is_symlink() or file_state(path) != item["state"]:
        raise RuntimeError(f"Missing or changed index file: {path}")
    if full or runtime_state(path) != item.get("runtime_state"):
        if sha256(path) != item["sha256"]:
            raise RuntimeError(f"Index SHA256 mismatch: {path}")


def resources(cfg, category):
    if category not in CATEGORIES:
        raise ValueError(f"Unknown reference category: {category}")
    result = dict(cfg)
    settings = cfg.get("category_resources", {})
    if not isinstance(settings, dict) or set(settings) - set(CATEGORIES):
        raise ValueError("category_resources must map the five recognized categories to resource mappings")
    values = settings.get(category, {})
    if not isinstance(values, dict):
        raise ValueError(f"{category}: resources must be a mapping")
    for field in ("queue", "cores", "threads", "memory_gib", "hours"):
        if field in values:
            result["build_" + field] = values[field]
    for key in ("build_cores", "build_threads"):
        if isinstance(result[key], bool) or not isinstance(result[key], int) or result[key] < 1:
            raise ValueError(f"{category}: {key} must be a positive integer")
    if result["build_threads"] > result["build_cores"]:
        raise ValueError(f"{category}: build_threads exceeds build_cores")
    for key in ("build_memory_gib", "build_hours"):
        if not isinstance(result[key], (int, float)) or isinstance(result[key], bool) or not math.isfinite(result[key]) or result[key] <= 0:
            raise ValueError(f"{category}: {key} must be positive and finite")
    return result


def _read(path):
    return json.loads(Path(path).read_text())


def committed_receipt(checkpoints):
    """Validate both small digest links before trusting the native payload hash."""
    checkpoints = Path(checkpoints)
    pointer = (checkpoints / "CURRENT").read_text().strip().split("\t")
    if len(pointer) != 2 or not re.fullmatch(r"cp-[A-Za-z0-9_-]+", pointer[0]) or not re.fullmatch(r"[0-9a-f]{64}", pointer[1]):
        raise RuntimeError("Invalid native checkpoint CURRENT pointer")
    folder = checkpoints / pointer[0]
    if folder.is_symlink() or sha256(folder / "receipt.tsv") != pointer[1]:
        raise RuntimeError("Native checkpoint receipt SHA256 mismatch")
    fields = {}
    for line in (folder / "receipt.tsv").read_text().splitlines():
        row = line.split("\t", 1)
        if len(row) != 2 or row[0] in fields:
            raise RuntimeError("Invalid canonical native receipt")
        fields[row[0]] = row[1]
    if sha256(folder / "receipt.json") != fields.get("receipt_json_sha256"):
        raise RuntimeError("Native checkpoint JSON SHA256 mismatch")
    return folder, _read(folder / "receipt.json")


def _preparation(cfg):
    from partitioned_reference import catalog_snapshot
    return catalog_snapshot(cfg)


def _category_snapshot(cfg, category, *, full=False):
    folder = index_root(cfg) / category
    value = _read(folder / "database.json")
    if value.get("schema") != SCHEMA or value.get("status") != "complete" or fingerprint(value.get("context")) != value.get("identity"):
        raise RuntimeError(f"Invalid completed category receipt: {category}")
    if value["context"].get("category") != category or value["context"].get("db_dir") != str(folder):
        raise RuntimeError(f"Category database has moved or belongs to another category: {folder}")
    for name, item in value["context"]["files"].items():
        if name not in ("hash.k2d", "taxo.k2d", "opts.k2d"):
            raise RuntimeError("Unexpected category index filename")
        check_record(folder / name, item, full=full)
    _native_index_structure(folder, expected_capacity=value["context"]["resource"]["capacity"])
    _native_opts(folder / "opts.k2d")
    return value


def _copy_checked(source, destination):
    source, destination = Path(source), Path(destination)
    digest = sha256(source)
    if destination.exists():
        if destination.is_symlink() or sha256(destination) != digest:
            raise RuntimeError(f"Frozen file differs: {destination}")
        return
    temporary = destination.with_name(destination.name + ".partial")
    shutil.copyfile(source, temporary)
    with temporary.open("rb") as handle:
        os.fsync(handle.fileno())
    temporary.replace(destination)


def _native_run(command, log, state, max_seconds):
    """Forward scheduler warnings to the checkpoint-aware child, never kill it on a timer."""
    child = None
    old = {}
    def request_pause(signum, frame):
        if child is not None and child.poll() is None:
            try:
                child.send_signal(signal.SIGUSR1)
            except ProcessLookupError:
                pass
        print(f"Signal {signum}: requested durable checkpoint at next complete input file", flush=True)
    with Path(log).open("ab", buffering=0) as output:
        for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1, signal.SIGUSR2):
            old[sig] = signal.signal(sig, request_pause)
        try:
            child = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT)
            atomic_json(state, {"status": "building", "pid": child.pid, "started_at": utcnow(),
                                "native_soft_seconds": max_seconds, "command": command})
            return child.wait()
        finally:
            for sig, handler in old.items():
                signal.signal(sig, handler)


def build_category(cfg, category):
    """Return 0 complete, 75 deliberately paused. Queue/RAM changes retain checkpoints."""
    started = time.monotonic()
    cfg = resources(cfg, category)
    prepared = _preparation(cfg)
    if category not in prepared["groups"]:
        raise RuntimeError(f"Category is absent from the frozen input plan: {category}")
    group = prepared["groups"][category]
    folder = index_root(cfg) / category
    folder.mkdir(parents=True, exist_ok=True)
    state = category_jobs(cfg, category)
    logs = category_logs(cfg, category)
    state.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    with locked(state / "build.lock"):
        # Published indexes are independent of build-time paths and thread count.
        if (folder / "database.json").exists():
            complete = _category_snapshot(cfg, category)
            if complete["context"]["preparation_identity"] != prepared["identity"]:
                raise RuntimeError("Completed category belongs to different prepared inputs")
            print(f"Reusing completed {category} index", flush=True)
            return 0
        from native_checkpoint import install
        installed = install(cfg)
        executable = Path(installed["executable"])
        threads = int(cfg["build_threads"])
        params = {**PARAMETERS, "threads": threads, "subblock_size": math.ceil(16384 / threads)}
        build_context = {"schema": SCHEMA, "category": category,
                         "preparation_identity": prepared["identity"],
                         "inputs_sha256": sha256(group["inputs_tsv"]),
                         "global_mapping_sha256": sha256(prepared["seqid_map"]),
                         "taxonomy": {p.name: sha256(p) for p in sorted(Path(prepared["taxonomy_dir"]).glob("*.dmp"))},
                         "parameters": params, "native_sha256": sha256(executable)}
        identity = fingerprint(build_context)
        state_path = state / "context.json"
        if state_path.exists():
            if _read(state_path).get("identity") != identity:
                raise RuntimeError(f"{category}: content/build thread identity changed; use a new parallel_reference_dir. Existing checkpoints were preserved.")
        else:
            allowed = {"build.lock"}
            if any(p.name not in allowed for p in folder.iterdir()):
                raise RuntimeError(f"Unowned index folder, refusing to overwrite: {folder}")
            atomic_json(state_path, {"identity": identity, "context": build_context})
        base = ["-k", "35", "-l", "31", "-S", PARAMETERS["seed"], "-p", str(threads)]
        capacity_file = state / "capacity.json"
        if capacity_file.exists():
            estimate_value = _read(capacity_file)
            expected = {"build_identity": identity, "estimate": estimate_value.get("estimate")}
            if estimate_value.get("identity") != fingerprint(expected) or estimate_value.get("build_identity") != identity:
                raise RuntimeError("Invalid category capacity checkpoint")
            estimate = int(estimate_value["estimate"])
        else:
            sources = [(Path(item["path"]), {"state": {"size": item["size"]}, "sha256": item["sha256"]}) for item in group["inputs"]]
            binary = Path(cfg["tools_dir"]) / "kraken/bin/estimate_capacity"
            print(f"{category}: estimating uncapped capacity from {len(sources)} prepared files", flush=True)
            _run_feed([str(binary), *base], sources, logs / "estimate.log", logs / "estimate.stdout")
            output = (logs / "estimate.stdout").read_text().strip()
            if not output.isdigit() or int(output) <= 0:
                raise RuntimeError(f"Invalid capacity estimate: {output!r}")
            estimate = int(output)
            estimate_value = {"build_identity": identity, "estimate": estimate}
            atomic_json(capacity_file, {**estimate_value, "identity": fingerprint(estimate_value), "created_at": utcnow()})
        capacity = (estimate + 8192) * 10 // 7
        hash_bytes = 32 + capacity * 4
        ram_bytes = math.ceil(hash_bytes * 1.15) + 8 * GIB + group["longest_record"] * (threads + 2)
        resource = {"capacity": capacity, "estimated_hash_bytes": hash_bytes,
                    "estimated_ram_bytes": ram_bytes, "budget_bytes": int(cfg["build_memory_gib"] * GIB),
                    "checkpoint_peak_hash_copies": 2, "no_capacity_cap": True}
        atomic_json(state / "resource_estimate.json", resource)
        print(f"{category}: hash~{hash_bytes/GIB:.1f} GiB, RAM estimate~{ram_bytes/GIB:.1f} GiB, budget={cfg['build_memory_gib']} GiB", flush=True)
        if ram_bytes > resource["budget_bytes"]:
            raise RuntimeError(f"{category}: uncapped RAM estimate exceeds this job. Increase this category's memory/queue and regenerate; capacity estimate and prepared references are reusable.")
        # Allow the next whole table generation in addition to the existing one.
        # When starting from scratch both generations may coexist at the next checkpoint.
        if not (state / "checkpoints").exists():
            check_space(folder, cfg["min_free_gib"] + 2 * hash_bytes / GIB)
        # On recovery, the native transaction validates CURRENT then reclaims
        # its abandoned stages before checking space for the next generation.
        # A Python preflight here could block that cleanup after a disk-full
        # interrupted write. Native errors always preserve the valid CURRENT.
        soft = min(float(cfg.get("build_soft_hours", 110)), float(cfg["build_hours"]) - float(cfg.get("checkpoint_headroom_hours", 10)))
        if soft <= 0:
            raise ValueError("build_hours must exceed checkpoint_headroom_hours")
        seconds = max(1, int(soft * 3600 - (time.monotonic() - started)))
        interval = int(cfg.get("checkpoint_hours", 4) * 3600)
        if interval < 1:
            raise ValueError("checkpoint_hours must be positive")
        command = [str(executable), *base, "-H", str(folder / "hash.k2d"),
                   "-t", str(folder / "taxo.k2d"), "-o", str(folder / "opts.k2d"),
                   "-n", str(prepared["taxonomy_dir"]), "-m", str(prepared["seqid_map"]),
                   "-c", str(capacity), "-B", "16384", "-b", str(params["subblock_size"]), "-r", "0",
                   "--metag-inputs", str(group["inputs_tsv"]), "--metag-context", identity,
                   "--metag-checkpoint-dir", str(state / "checkpoints"),
                   "--metag-checkpoint-seconds", str(interval), "--metag-max-seconds", str(seconds)]
        if cfg.get("native_stop_after_files"):
            command += ["--metag-stop-after-files", str(cfg["native_stop_after_files"])]
        code = _native_run(command, logs / "build.log", state / "process.json", seconds)
        if code == 75:
            atomic_json(state / "process.json", {"status": "paused", "finished_at": utcnow(), "exit_code": 75})
            print(f"{category}: durable pause; resubmit the same content/thread configuration to resume", flush=True)
            return 75
        if code:
            atomic_json(state / "process.json", {"status": "failed", "finished_at": utcnow(), "exit_code": code})
            raise RuntimeError(f"{category}: native builder exited {code}; inspect {logs / 'build.log'}. Existing committed checkpoints remain.")
        _native_index_structure(folder, expected_capacity=capacity)
        options = _native_opts(folder / "opts.k2d")
        generation, receipt = committed_receipt(state / "checkpoints")
        if not receipt.get("complete") or receipt.get("context") != identity or receipt.get("completed_files") != len(group["inputs"]):
            raise RuntimeError("Native builder returned success without a complete checkpoint")
        if receipt.get("sequences") != group["records"] or receipt.get("bases") != group["bases"]:
            raise RuntimeError("Completed native counters differ from the frozen category inputs")
        # Native computes hash digest after fsync; the completed -H is its hardlink.
        checkpoint_hash = generation / "hash.k2d"
        if not os.path.samefile(folder / "hash.k2d", checkpoint_hash):
            raise RuntimeError("Published hash is not the committed native payload")
        files = {name: record(folder / name, receipt["hash_sha256"] if name == "hash.k2d" else None)
                 for name in ("hash.k2d", "taxo.k2d", "opts.k2d")}
        context = {"category": category, "db_dir": str(folder), "build_identity": identity,
                   "preparation_identity": prepared["identity"], "parameters": params,
                   "files": files, "resource": resource, "native_options": options}
        atomic_json(folder / "database.json", {"schema": SCHEMA, "status": "complete",
                    "identity": fingerprint(context), "context": context, "completed_at": utcnow()})
        atomic_json(state / "process.json", {"status": "complete", "finished_at": utcnow(), "exit_code": 0,
                                             "records": group["records"], "bases": group["bases"]})
        print(f"Completed {category}: {folder}", flush=True)
        return 0


def publish(cfg):
    prepared = _preparation(cfg)
    folder = index_root(cfg)
    folder.mkdir(parents=True, exist_ok=True)
    with locked((library_jobs(cfg) if separated(cfg) else folder) / "catalog.lock"):
        if (folder / "database.json").exists():
            saved = database_snapshot(cfg)
            if saved["context"]["preparation_identity"] != prepared["identity"]:
                raise RuntimeError("Existing catalog belongs to another prepared input set")
            return saved
        categories = [name for name in CATEGORIES if name in prepared["groups"]]
        values = {name: _category_snapshot(cfg, name) for name in categories}
        if not values:
            raise RuntimeError("No completed category indexes")
        for name, value in values.items():
            if value["context"]["preparation_identity"] != prepared["identity"]:
                raise RuntimeError(f"Mismatched category preparation: {name}")
        first = folder / categories[0]
        for name in categories[1:]:
            for item in ("taxo.k2d", "opts.k2d"):
                if values[name]["context"]["files"][item]["sha256"] != values[categories[0]]["context"]["files"][item]["sha256"]:
                    raise RuntimeError(f"Category {name} does not share exactly the same taxonomy/options")
        for name in ("taxo.k2d", "opts.k2d"):
            _copy_checked(first / name, folder / name)
        hashes = [f"{name}/hash.k2d" for name in categories]
        hash_list = folder / "hashes.txt"
        text = "".join(str(folder / name) + "\n" for name in hashes)
        if hash_list.exists() and hash_list.read_text() != text:
            raise RuntimeError("Existing union hash list differs")
        if not hash_list.exists():
            with (folder / "hashes.txt.partial").open("w") as handle:
                handle.write(text)
                handle.flush()
                os.fsync(handle.fileno())
            (folder / "hashes.txt.partial").replace(hash_list)
        files = {f"{name}/{item}": record(folder / name / item, details["sha256"])
                 for name in categories for item, details in values[name]["context"]["files"].items()}
        files.update({name: record(folder / name) for name in ("taxo.k2d", "opts.k2d", "hashes.txt")})
        taxonomy = Path(prepared["taxonomy_dir"])
        taxonomy_files = {p.name: record(p) for p in sorted(taxonomy.glob("*.dmp"))}
        context = {"layout": "partitioned", "preparation_identity": prepared["identity"],
                   "categories": {name: values[name]["identity"] for name in categories},
                   "absent_categories": [name for name in CATEGORIES if name not in categories],
                   "db_dir": str(folder), "taxonomy_dir": str(taxonomy),
                   "hashes": hashes, "files": files, "taxonomy_files": taxonomy_files,
                   "parameters": dict(PARAMETERS), "classification": "per-minimizer-LCA-before-thresholds"}
        if separated(cfg):
            context["verification_path"] = str(library_jobs(cfg) / "verified-runtime.json")
        value = {"schema": SCHEMA, "status": "complete", "identity": fingerprint(context), "context": context,
                 **{key: context[key] for key in ("layout", "db_dir", "taxonomy_dir", "hashes", "files", "taxonomy_files", "parameters")},
                 "completed_at": utcnow()}
        if separated(cfg):
            value["verification_path"] = context["verification_path"]
        atomic_json(folder / "database.json", value)
        return value


def validate_snapshot(value, *, full=False):
    """Check published structure cheaply; read all payload bytes only with full=True.

    The build catalog remains the size/SHA baseline. File timestamps and the
    optional audit receipt are not prerequisites for using a completed library.
    """
    if value.get("schema") != SCHEMA or value.get("status") != "complete" or fingerprint(value.get("context")) != value.get("identity"):
        raise RuntimeError("No valid complete partitioned database catalog")
    context = value["context"]
    for key in ("layout", "db_dir", "taxonomy_dir", "hashes", "files", "taxonomy_files", "parameters"):
        if value.get(key) != context.get(key):
            raise RuntimeError(f"Catalog field is inconsistent: {key}")
    if value.get("verification_path") != context.get("verification_path"):
        raise RuntimeError("Catalog verification path is inconsistent")
    db = Path(value["db_dir"])
    categories = list(context.get("categories", {}))
    expected_hashes = [f"{name}/hash.k2d" for name in CATEGORIES if name in categories]
    if not categories or any(name not in CATEGORIES for name in categories) or value["hashes"] != expected_hashes:
        raise RuntimeError("Catalog must contain exactly one hash for every recorded reference category")
    if sorted(context.get("absent_categories", [])) != sorted(set(CATEGORIES) - set(categories)):
        raise RuntimeError("Catalog must explicitly record every absent category")
    required = {"taxo.k2d", "opts.k2d", "hashes.txt"} | {f"{name}/{part}" for name in categories for part in ("hash.k2d", "taxo.k2d", "opts.k2d")}
    if set(value["files"]) != required or not {"nodes.dmp", "names.dmp", "merged.dmp"} <= set(value["taxonomy_files"]):
        raise RuntimeError("Catalog lacks exact category index or shared taxonomy files")
    saved = _read(db / "database.json")
    if saved.get("identity") != value["identity"] or saved != value:
        raise RuntimeError("Partitioned catalog changed")
    observed = {}
    for directory, records in ((db, value["files"]), (Path(value["taxonomy_dir"]), value["taxonomy_files"])):
        for name, item in records.items():
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise RuntimeError("Unsafe path in partitioned catalog")
            path = directory / name
            if not path.is_file() or path.is_symlink() or file_state(path) != item["state"]:
                raise RuntimeError(f"Missing or changed index file: {path}")
            state = runtime_state(path)
            if full:
                check_record(path, item, full=True)
                after = runtime_state(path)
                if any(after[key] != state[key] for key in ("size", "mtime_ns")):
                    raise RuntimeError(f"Index changed during full verification: {path}")
                observed[str(path)] = {"state": after, "sha256": item["sha256"]}
    if (db / "hashes.txt").read_text() != "".join(str(db / name) + "\n" for name in value["hashes"]):
        raise RuntimeError("Union hash list differs from committed catalog")
    for name in value["hashes"]:
        _native_index_structure((db / name).parent)
    if full:
        verification_path = Path(value.get("verification_path") or db / "verified-runtime.json")
        try:
            atomic_json(verification_path, {"identity": value["identity"], "files": observed, "verified_at": utcnow()})
        except OSError as exc:
            print(f"Warning: SHA256 verification passed; could not write optional audit record {verification_path}: {exc}", file=sys.stderr)
    return value


def database_snapshot(cfg, *, full=False):
    return validate_snapshot(_read(index_root(cfg) / "database.json"), full=full)
