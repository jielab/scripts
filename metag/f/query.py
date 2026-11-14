"""Sample discovery, isolated queries and atomic publication for one HPC job.

Workers read the same database paths. Work is deliberately staged on the
configured shared filesystem; no node-local capacity or TMPDIR is assumed.
"""

import sys as _layout_sys
from pathlib import Path as _LayoutPath
_code_root = _LayoutPath(__file__).resolve().parents[1]
_layout_root = _code_root.parent
_layout_sys.path.append(str(_layout_root / "0f"))
from workflow_layout import named_directory

from collections import deque
from contextlib import contextmanager, ExitStack
import json
import multiprocessing
from multiprocessing.connection import wait
import os
from pathlib import Path
import re
import shutil
import signal
import tarfile
import threading
import time
import traceback
import tempfile

from common import atomic_json, check_space, file_state, fingerprint, locked, sha256, utcnow
from ordinal_query import input_states
import ordinal_query
import kunpeng_query
from species_split import split_context, split_originals, verify_host_context

SCHEMA = "metag-hpc-sample-v3"
OWNER = ".metag-owner.json"
FASTQ_RE = re.compile(r"(.+?)(?:_R?([12])(?:_001)?|\.([12]))?\.fastq\.gz$", re.I)


def _engine(cfg):
    value = cfg.get("classifier_engine", "kraken")
    if value not in ("kraken", "kunpeng_full", "kunpeng_memory"):
        raise ValueError(f"Unsupported classifier_engine: {value}")
    return value


def resolve_native(cfg):
    return (ordinal_query if _engine(cfg) == "kraken" else kunpeng_query).resolve_native(cfg)


def validate_snapshot(snapshot):
    module = kunpeng_query if "hash_config.k2d" in snapshot.get("files", {}) else ordinal_query
    return module.validate_snapshot(snapshot)


def database_snapshot(cfg):
    from common import engine_snapshot
    return engine_snapshot(cfg)


def query_once(cfg, sample, snapshot, work, limit=0, taxonomy=None):
    is_kunpeng = "hash_config.k2d" in snapshot.get("files", {})
    if is_kunpeng != (_engine(cfg) != "kraken"):
        raise RuntimeError("Selected classifier and database format differ")
    module = kunpeng_query if is_kunpeng else ordinal_query
    if hasattr(taxonomy, "validate_current"):
        taxonomy.validate_current()
    result = module.query_once(cfg, sample, snapshot, work, limit, taxonomy=taxonomy)
    if hasattr(taxonomy, "validate_current"):
        taxonomy.validate_current()
    return result


def _name(value):
    if not isinstance(value, str) or not value or value in (".", "..") or re.search(r"[/\\\x00-\x1f]", value):
        raise ValueError(f"Expected one safe sample/dataset component: {value!r}")
    return value


def _selectors(values):
    if values is None:
        return []
    if isinstance(values, str):
        return [x.strip() for x in values.split(",") if x.strip()]
    return list(values)


def discover_samples(cfg, datasets=None, only=None, ready_only=False):
    """Snapshot ready paired .fastq.gz inputs; warn and skip incomplete arrivals.

    Incomplete pairs are always skipped. Workers verify the frozen input
    states before processing.
    """
    from project_cache import inventory, probe_samples
    raw = Path(cfg["reads_raw"]).resolve(strict=True)
    wanted = _selectors(datasets)
    selected = _selectors(only)
    for selector in selected:
        parts = selector.split("/")
        if len(parts) not in (1, 2):
            raise ValueError(f"Use sample or dataset/sample selector: {selector!r}")
        for part in parts:
            _name(part)
    wanted = wanted or sorted(p.name for p in raw.iterdir() if p.is_dir() and not p.is_symlink())
    if not wanted or len(wanted) != len(set(wanted)):
        raise RuntimeError("No datasets found, or duplicate datasets selected")
    result, candidates = [], set()
    for dataset in wanted:
        _name(dataset)
        folder = raw / dataset
        if not folder.is_dir():
            raise RuntimeError(f"Dataset not found: {folder}")
        groups = {}
        for path in inventory(cfg, dataset):
            match = FASTQ_RE.fullmatch(path.name)
            if not match or not path.is_file():
                continue
            stem, mate = match[1], match[2] or match[3] or "single"
            _name(stem)
            if selected and stem not in selected and f"{dataset}/{stem}" not in selected:
                continue
            candidates.add(f"{dataset}/{stem}")
            files = groups.setdefault(stem, {})
            if mate in files:
                raise RuntimeError(f"Ambiguous FASTQ copies for {dataset}/{stem}: {files[mate]} and {path}")
            if files and any(previous.parent != path.parent for previous in files.values()):
                raise RuntimeError(f"FASTQ mates are in different folders: {dataset}/{stem}")
            files[mate] = path
        if not groups and not selected:
            print(f"WARNING: no .fastq.gz candidates in {folder}; skipping dataset", flush=True)
        for stem, files in sorted(groups.items()):
            if set(files) != {"1", "2"}:
                missing = ", ".join("R" + mate for mate in ("1", "2") if mate not in files)
                reason = f"missing mate(s): {missing}" if missing else "mixed paired/unpaired filenames"
                print(f"WARNING: skipping {dataset}/{stem}: {reason}; requires paired .fastq.gz files", flush=True)
                continue
            paths = []
            for mate in ("1", "2"):
                path = files[mate]
                try:
                    real = path.resolve(strict=True)
                    if not real.stat().st_size or Path(str(path) + ".aria2").exists() or Path(str(real) + ".aria2").exists():
                        print(f"WARNING: skipping {dataset}/{stem}: empty or downloading FASTQ: {path}", flush=True)
                        break
                except FileNotFoundError:
                    print(f"WARNING: skipping {dataset}/{stem}: FASTQ disappeared during discovery: {path}", flush=True)
                    break
                paths.append(str(real))
            if len(paths) != 2:
                continue
            if Path(paths[0]).samefile(paths[1]):
                raise RuntimeError(f"R1 and R2 resolve to the same file: {dataset}/{stem}")
            result.append({"sample_id": f"{dataset}/{stem}", "dataset": dataset,
                           "sample": stem, "r1": paths[0], "r2": paths[1]})
    for selector in selected:
        matches = [name for name in candidates if selector in (name, name.split("/", 1)[1])]
        if len(matches) != 1:
            raise RuntimeError(f"Selector {selector!r} matched {len(matches)} samples; use dataset/sample")
    if not result:
        print("WARNING: no complete .fastq.gz pairs are ready; nothing to process", flush=True)
    else:
        probe_samples(cfg, result)
    return result


@contextmanager
def _signals():
    previous = {}
    if threading.current_thread() is threading.main_thread():
        def interrupted(signum, frame):
            raise InterruptedError(f"HPC sample interrupted by signal {signum}")
        for signum in (signal.SIGTERM, signal.SIGINT):
            previous[signum] = signal.getsignal(signum)
            signal.signal(signum, interrupted)
    try:
        yield
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def _directory(path):
    path = Path(path)
    if path.is_symlink():
        raise RuntimeError(f"Refusing symlink at writable path: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _owned(path, owner, reset=False):
    path = Path(path)
    if path.is_symlink():
        raise RuntimeError(f"Refusing symlink at private stage: {path}")
    if path.exists():
        try:
            saved = json.loads((path / OWNER).read_text())
        except (OSError, ValueError) as exc:
            raise RuntimeError(f"Refusing unowned staging directory: {path}") from exc
        if saved != owner and not (reset and saved.get("target") == owner.get("target")):
            raise RuntimeError(f"Stage belongs to a different output: {path}")
        if reset:
            shutil.rmtree(path)
    if not path.exists():
        path.mkdir(parents=True)
        atomic_json(path / OWNER, owner)
    return path


def _remove_owned(path, owner):
    if Path(path).exists():
        _owned(path, owner)
        shutil.rmtree(path)


def _sync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _runtime_state(path):
    """Cheap local change detection, separate from portable size-only receipts.

    This is conservative about metadata changes and is not a content checksum.
    Keep reference/preparation receipt identities unchanged; query publication
    additionally records these guards for its own inputs and output files.
    """
    value = Path(path).stat()
    return {"size": value.st_size, "mtime_ns": value.st_mtime_ns,
            "ctime_ns": value.st_ctime_ns}


def _input_runtime_states(sample):
    result = []
    for mate in ("r1", "r2"):
        if sample.get(mate) is not None:
            path = Path(sample[mate]).resolve(strict=True)
            result.append({"path": str(path), "state": _runtime_state(path)})
    return result


def _database_runtime_states(snapshot):
    """Detect writes during a sample without treating chmod as a library edit."""
    result = {}
    for directory, records in ((snapshot["db_dir"], snapshot["files"]),
                               (snapshot["taxonomy_dir"], snapshot["taxonomy_files"])):
        for name in records:
            path = (Path(directory) / name).resolve(strict=True)
            state = _runtime_state(path)
            if snapshot.get("schema") == "metag-partitioned-index-v1":
                state = {key: state[key] for key in ("size", "mtime_ns")}
            result[str(path)] = state
    return result


def _manifest(folder, exclude=()):
    """Record file sizes without another read of every FASTQ."""
    result = {}
    for path in sorted(Path(folder).rglob("*")):
        name = path.relative_to(folder).as_posix()
        if path.is_symlink():
            raise RuntimeError(f"Symlink appeared in private output: {path}")
        if path.is_file() and name not in exclude:
            result[name] = {"state": file_state(path)}
    return result


def _complete(folder, identity=None, verify_hashes=False):
    """Cheap explicit reuse: complete receipt, file sizes and conserved counts.

    --replace FALSE deliberately reuses an existing complete result. It is not
    a parameter/content identity cache. Missing or incomplete results rebuild.
    """
    folder = Path(folder)
    if folder.is_symlink():
        raise RuntimeError(f"Refusing symlink output: {folder}")
    marker = folder / "qc/complete.json"
    if not marker.is_file():
        return None
    try:
        saved = json.loads(marker.read_text())
        if saved.get("schema") != SCHEMA or saved.get("status") != "complete":
            return None
        expected = saved.get("outputs", {})
        if not expected:
            return None
        for name, item in expected.items():
            rel = Path(name)
            path = folder / rel
            if rel.is_absolute() or ".." in rel.parts or path.is_symlink() or path.resolve() != path or not path.is_file():
                return None
            if path.stat().st_size != item["state"]["size"]:
                return None
        summary = saved["summary"]
        if not summary.get("count_conservation") or int(summary["fragments"]) <= 0:
            return None
        if int(summary["reads"]) != int(summary["fragments"]) * (2 if summary["paired"] else 1):
            return None
        if int(summary.get("requested_limit", 0)) == 0:
            from split_summary import load_sample
            load_sample(folder, summary["sample_id"])
        return saved
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _context(cfg, sample, snapshot, limit):
    native = resolve_native(cfg)
    source = Path(__file__).resolve().parent
    return {"schema": SCHEMA, "sample": sample, "inputs": input_states(sample),
            "source_release": str(source.parent),
            "input_runtime_states": _input_runtime_states(sample),
            "database_runtime_states": _database_runtime_states(snapshot),
            "classifier_engine": _engine(cfg),
            "engine_parameters": kunpeng_query.parameters(cfg) if _engine(cfg) != "kraken" else {},
            "database_identity": snapshot["identity"], "database_files": snapshot["files"],
            "source_database_identity": snapshot.get("source_database_identity", snapshot["identity"]),
            "taxonomy_files": snapshot["taxonomy_files"],
            "db_dir": str(Path(snapshot["db_dir"]).resolve()),
            "taxonomy_dir": str(Path(snapshot["taxonomy_dir"]).resolve()),
            "native": {"path": str(native), "state": _runtime_state(native)},
            "confidence": float(cfg.get("confidence", 0.1)),
            "minimum_hit_groups": int(cfg.get("minimum_hit_groups", 2)),
            "minimum_base_quality": int(cfg.get("minimum_base_quality", 0)),
            "split_policy": split_context(cfg),
            "output_format": cfg.get("output_format", "directory"), "limit": limit,
            "code_release": str(source.parent)}


def _verify_sources(sample, snapshot, context):
    verify_host_context(context["split_policy"]["human_host"])
    if input_states(sample) != context["inputs"]:
        raise RuntimeError("Original FASTQ files changed since this run began")
    if _input_runtime_states(sample) != context["input_runtime_states"]:
        raise RuntimeError("Original FASTQ runtime state changed since this run began")
    validate_snapshot(snapshot)
    if _database_runtime_states(snapshot) != context["database_runtime_states"]:
        raise RuntimeError("Database or taxonomy runtime state changed since this run began")
    if _runtime_state(context["native"]["path"]) != context["native"]["state"]:
        raise RuntimeError("Native executable changed during this sample")


def _roots(cfg):
    out, work = Path(cfg["output_dir"]).resolve(), Path(cfg["work_dir"]).resolve()
    if cfg.get("storage_layout") == "separated":
        work = work / "tmp"
    # work_dir also contains scheduler plans, metadata and library directories
    # in the default HPC layout.  Only its queries/ subtree is sample scratch
    # that this module can create or recursively remove.
    query_root = work / "queries"
    if query_root.resolve() != query_root:
        raise RuntimeError(f"Symlink/alias in query work root: {query_root}")
    if out == query_root or out in query_root.parents or query_root in out.parents:
        raise RuntimeError(f"Query output and scratch roots overlap: {out} / {query_root}")
    for key in ("reads_raw", "reference_raw", "metadata_dir", "library_dir", "tools_dir"):
        protected = Path(cfg[key]).resolve()
        for path in (out, query_root):
            if path == protected or protected in path.parents or path in protected.parents:
                raise RuntimeError(f"Writable query root overlaps protected {key}: {path}")
    return _directory(out), _directory(work)


def _publish_sample(cfg, sample, snapshot, limit, prepared_query=None, prepared_taxonomy=None):
    """Build a fresh result, then replace the sample directory under its lock."""
    started = time.monotonic()
    _name(sample["dataset"])
    _name(sample["sample"])
    label = _name(cfg.get("run_label", "split"))
    if sample["sample_id"] != f"{sample['dataset']}/{sample['sample']}":
        raise ValueError("sample_id must equal dataset/sample")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")
    if cfg.get("output_format", "directory") not in ("tar", "directory"):
        raise ValueError("output_format must be tar or directory")
    validate_snapshot(snapshot)
    output_root, work_root = _roots(cfg)
    context = _context(cfg, sample, snapshot, limit)
    if hasattr(prepared_taxonomy, "validate_current"):
        prepared_taxonomy.validate_current()
    identity = fingerprint(context)
    target = (named_directory(output_root / "pilots" / sample["dataset"] / sample["sample"], identity, "pilot")
              if limit else output_root / sample["dataset"] / sample["sample"])
    if target.resolve() != target:
        raise RuntimeError(f"Symlink/alias in sample output path: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    lock = Path(cfg["work_dir"]) / "jobs/locks" / label / sample["dataset"] / (sample["sample"] + ".lock")
    with locked(lock):
        completed = _complete(target) if not cfg.get("replace", True) else None
        if completed:
            return {"status": "skipped", "sample_id": sample["sample_id"],
                    "identity": completed["identity"], "output": str(target), "summary": completed["summary"]}
        owner = {"schema": SCHEMA, "identity": identity, "target": str(target)}
        if target.exists() and (target.is_symlink() or not target.is_dir()):
            raise RuntimeError(f"Expected a sample directory, not a link/file: {target}")
        reserve = float(cfg.get("min_free_gib", 100))
        check_space(work_root, reserve)
        check_space(target.parent, reserve)
        work = work_root / "queries" / label
        if limit:
            work /= "pilots"
        work = work / sample["dataset"] / sample["sample"]
        if limit:
            work /= target.name
        if work.resolve() != work:
            raise RuntimeError(f"Symlink/alias in sample work path: {work}")
        _owned(work, owner, reset=True)
        query_dir = work / "query"
        if prepared_query is None:
            query_dir.mkdir()
            query = query_once(cfg, sample, snapshot, query_dir, limit, taxonomy=prepared_taxonomy)
        else:
            query = prepared_query
            if query["sample_id"] != sample["sample_id"] or query["requested_limit"] != limit:
                raise RuntimeError("Prepared classification belongs to another sample/scope")
        if prepared_taxonomy is not None and hasattr(prepared_taxonomy, "preparation"):
            query["project_preparation"] = prepared_taxonomy.preparation
        _verify_sources(sample, snapshot, context)
        stage = work / "publish"
        stage.mkdir()
        payload = stage / "data"
        summary = split_originals(cfg, sample, query, snapshot, payload, host_work=stage / "host",
                                  frozen_host_context=context["split_policy"]["human_host"], taxonomy=prepared_taxonomy)
        metadata = payload / "qc"
        atomic_json(metadata / "provenance.json", {
            "schema": SCHEMA, "identity": identity, "created_utc": utcnow(), "context": context,
            "query": query, "classifier": _engine(cfg), "classify_threads": cfg.get("classify_threads", 8),
            "publication": "fresh stage; complete receipt; directory replacement under sample lock",
            "prefix_note": "prefix is a throughput pilot, not a recall validation" if limit else ""})
        log_root = Path(cfg["work_dir"]) / "logs" / label / "samples" / sample["dataset"] / sample["sample"]
        if limit:
            log_root /= target.name
        log_root.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(query.get("log_path", query_dir / "query.log"), log_root / "classifier.log")
        archive_started = time.monotonic()
        if cfg.get("output_format", "directory") == "tar":
            bins = [path for path in payload.iterdir() if path.name != "qc"]
            files = [file for item in bins for file in ([item] if item.is_file() else item.rglob("*")) if file.is_file()]
            check_space(stage, reserve + sum(path.stat().st_size for path in files) / 2**30 + 0.1)
            with tarfile.open(payload / "reads.tar", "w", format=tarfile.PAX_FORMAT) as tar:
                for path in sorted(files):
                    if path.is_symlink():
                        raise RuntimeError("Symlink appeared during archive publication")
                    tar.add(path, arcname=path.relative_to(payload).as_posix(), recursive=False)
            for item in bins:
                shutil.rmtree(item) if item.is_dir() else item.unlink()
        summary.update(archive_seconds=time.monotonic() - archive_started,
                       query_seconds_this_run=query["seconds"],
                       total_seconds_before_commit=time.monotonic() - started,
                       artifact_bytes=sum(path.stat().st_size for path in payload.rglob("*")
                                          if path.is_file() and "qc" not in path.relative_to(payload).parts),
                       timing_note=query.get("classification_timing", "sample wall time"),
                       classifier_log=str(log_root / "classifier.log"))
        atomic_json(metadata / "summary.json", summary)
        _verify_sources(sample, snapshot, context)
        # Rename avoids a second full write when scratch/output share a filesystem.
        # Across filesystems copy once into a private sibling before publishing.
        publication_started = time.monotonic()
        destination = Path(tempfile.mkdtemp(prefix=f".{sample['sample']}.publish-", dir=target.parent))
        backup = None
        try:
            same_fs = payload.stat().st_dev == target.parent.stat().st_dev
            if same_fs:
                destination.rmdir()
                os.replace(payload, destination)
            else:
                check_space(target.parent, reserve + sum(p.stat().st_size for p in payload.rglob("*") if p.is_file()) / 2**30)
                shutil.copytree(payload, destination, dirs_exist_ok=True)
            summary.update(publication_seconds=time.monotonic() - publication_started,
                           publication_method="same-filesystem rename" if same_fs else "cross-filesystem copy",
                           total_seconds_before_commit=time.monotonic() - started)
            atomic_json(destination / "qc/summary.json", summary)
            atomic_json(destination / "qc/owner.json", owner)
            outputs = _manifest(destination, exclude=("qc/owner.json", "qc/complete.json"))
            atomic_json(destination / "qc/complete.json", {
                "schema": SCHEMA, "status": "complete", "identity": identity,
                "completed_utc": utcnow(), "summary": summary, "outputs": outputs})
            _verify_sources(sample, snapshot, context)
            if target.exists():
                backup = Path(tempfile.mkdtemp(prefix=f".{sample['sample']}.old-", dir=target.parent))
                backup.rmdir()
                os.replace(target, backup)
            try:
                os.replace(destination, target)
            except BaseException:
                if backup is not None and backup.exists() and not target.exists():
                    os.replace(backup, target)
                raise
            _sync_directory(target.parent)
            if backup is not None:
                shutil.rmtree(backup)
        finally:
            if destination.exists():
                shutil.rmtree(destination)
        _remove_owned(work, owner)
        print(f"[{sample['sample_id']}] PASS {summary['fragments']:,} fragments -> {target}", flush=True)
        return {"status": "complete", "sample_id": sample["sample_id"], "identity": identity,
                "output": str(target), "summary": summary}


def run_sample(cfg, sample, snapshot=None, limit=0, prepared_query=None, prepared_taxonomy=None):
    if snapshot is None:
        snapshot = database_snapshot(cfg)
    with _signals(), ExitStack() as stack:
        if prepared_taxonomy is None:
            from project_cache import load_taxonomy
            prepared_taxonomy = stack.enter_context(load_taxonomy(cfg, snapshot, [sample["dataset"]]))
        return _publish_sample(cfg, sample, snapshot, limit, prepared_query=prepared_query, prepared_taxonomy=prepared_taxonomy)


class BatchFailure(RuntimeError):
    def __init__(self, summary):
        self.summary = summary
        super().__init__(f"{len(summary['failed'])}/{summary['samples']} samples failed; completed samples remain reusable")


def _worker(connection, cfg, sample, snapshot):
    try:
        connection.send({"ok": True, "result": run_sample(cfg, sample, snapshot)})
    except BaseException as exc:
        connection.send({"ok": False, "sample_id": sample["sample_id"],
                         "error": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()[-8000:]})
    finally:
        connection.close()


def run_batch(cfg, samples, snapshot=None):
    """Bound worker processes; preserve successes and raise after all failures are reported."""
    samples = list(samples)
    if not samples or len({s["sample_id"] for s in samples}) != len(samples):
        raise ValueError("Batch requires a nonempty list of distinct sample_id values")
    workers = int(cfg.get("classify_workers", 4))
    if workers < 1:
        raise ValueError("classify_workers must be positive")
    if snapshot is None:
        snapshot = database_snapshot(cfg)
    validate_snapshot(snapshot)
    if workers == 1 and _engine(cfg) == "kraken" and cfg.get("batch_query", True):
        return _resident_batch(cfg, samples, snapshot)
    ctx = multiprocessing.get_context("spawn")
    pending, active, successes, failures = deque(samples), {}, [], []
    started = time.monotonic()
    with _signals():
        try:
            while pending or active:
                while pending and len(active) < workers:
                    sample = pending.popleft()
                    receive, send = ctx.Pipe(duplex=False)
                    process = ctx.Process(target=_worker, args=(send, cfg, sample, snapshot),
                                          name="metag-" + sample["sample_id"].replace("/", "_"))
                    process.start()
                    send.close()
                    active[receive] = (process, sample)
                ready = wait(list(active), timeout=0.2)
                for receive in list(active):
                    process, sample = active[receive]
                    if receive not in ready and process.is_alive():
                        continue
                    try:
                        message = receive.recv()
                    except EOFError:
                        message = {"ok": False, "sample_id": sample["sample_id"],
                                   "error": f"Worker exited without result, exitcode={process.exitcode}"}
                    process.join(timeout=10)
                    if process.is_alive():
                        process.terminate()
                        process.join(timeout=10)
                        message = {"ok": False, "sample_id": sample["sample_id"], "error": "Worker did not exit after reporting"}
                    receive.close()
                    del active[receive]
                    if message.get("ok") and process.exitcode == 0:
                        successes.append(message["result"])
                    else:
                        failures.append(message)
                        print(f"FAIL {sample['sample_id']}: {message.get('error', 'worker exit failure')}", flush=True)
        finally:
            # Workers handle SIGTERM and terminate their own native process groups.
            for process, sample in active.values():
                if process.is_alive():
                    process.terminate()
            for receive, (process, sample) in active.items():
                process.join(timeout=15)
                if process.is_alive():
                    process.kill()
                    process.join(timeout=5)
                receive.close()
    summary = {"samples": len(samples), "workers": workers, "succeeded": successes,
               "failed": failures, "seconds": time.monotonic() - started}
    if failures:
        raise BatchFailure(summary)
    return summary


def _resident_batch(cfg, samples, snapshot):
    """One native invocation/index load for all pending samples in this job."""
    started = time.monotonic()
    successes, pending, failures = [], [], []
    output_root, work_root = _roots(cfg)
    for sample in samples:
        target = output_root / sample["dataset"] / sample["sample"]
        completed = _complete(target) if not cfg.get("replace", True) else None
        if completed:
            successes.append({"status": "skipped", "sample_id": sample["sample_id"],
                              "identity": completed["identity"], "output": str(target),
                              "summary": completed["summary"]})
        else:
            pending.append(sample)
    native_invocations = 0
    if pending:
        batch_root = work_root / "queries" / "resident-batches"
        batch_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="batch-", dir=batch_root) as temporary, _signals(), ExitStack() as stack:
            try:
                native_invocations = 1
                from project_cache import load_taxonomy
                taxonomy_started = time.monotonic()
                taxonomy = stack.enter_context(load_taxonomy(cfg, snapshot, [sample["dataset"] for sample in pending]))
                taxonomy_seconds = time.monotonic() - taxonomy_started
                queries = ordinal_query.query_many(cfg, pending, snapshot, Path(temporary), taxonomy=taxonomy)
                taxonomy.validate_current()
            except BaseException as exc:
                if isinstance(exc, (KeyboardInterrupt, SystemExit, InterruptedError)):
                    raise
                queries = []
                failures.extend({"ok": False, "sample_id": sample["sample_id"],
                                 "error": f"Batch classification: {type(exc).__name__}: {exc}"} for sample in pending)
            for sample, prepared in zip(pending, queries):
                prepared["batch_taxonomy_load_seconds"] = taxonomy_seconds
                try:
                    successes.append(run_sample(cfg, sample, snapshot, prepared_query=prepared, prepared_taxonomy=taxonomy))
                except BaseException as exc:
                    if isinstance(exc, (KeyboardInterrupt, SystemExit, InterruptedError)):
                        raise
                    failures.append({"ok": False, "sample_id": sample["sample_id"],
                                     "error": f"{type(exc).__name__}: {exc}",
                                     "traceback": traceback.format_exc()[-8000:]})
                    print(f"FAIL {sample['sample_id']}: {exc}", flush=True)
    summary = {"samples": len(samples), "workers": 1, "succeeded": successes,
               "failed": failures, "seconds": time.monotonic() - started,
               "native_invocations": native_invocations, "resident_index_batch": True}
    if failures:
        raise BatchFailure(summary)
    return summary
