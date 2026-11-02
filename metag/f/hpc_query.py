"""Sample discovery, isolated queries and atomic publication for one HPC job.

Workers mmap the same database paths. Work is deliberately staged on the
configured shared filesystem; no node-local capacity or TMPDIR is assumed.
"""
from collections import deque
from contextlib import contextmanager
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

from hpc_common import atomic_json, check_space, file_state, fingerprint, locked, sha256, utcnow
from ordinal_query import input_states
import ordinal_query
import kunpeng_query
from species_split import split_originals

SCHEMA = "metag-hpc-sample-v1"
OWNER = ".metag-owner.json"
FASTQ_RE = re.compile(r"(.+?)(?:_R?([12])(?:_001)?|\.([12]))?\.(?:fastq|fq)(?:\.gz)?$", re.I)


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
    from hpc_common import engine_snapshot
    return engine_snapshot(cfg)


def query_once(cfg, sample, snapshot, work, limit=0):
    is_kunpeng = "hash_config.k2d" in snapshot.get("files", {})
    if is_kunpeng != (_engine(cfg) != "kraken"):
        raise RuntimeError("Selected classifier and database format differ")
    module = kunpeng_query if is_kunpeng else ordinal_query
    return module.query_once(cfg, sample, snapshot, work, limit)


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
    """Find actual files, rejecting ambiguous copies, missing mates and mixed sets."""
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
    result = []
    for dataset in wanted:
        _name(dataset)
        folder = raw / dataset
        if not folder.is_dir():
            raise RuntimeError(f"Dataset not found: {folder}")
        groups = {}
        for path in sorted(folder.rglob("*")):
            if not path.is_file():
                continue
            match = FASTQ_RE.fullmatch(path.name)
            if not match:
                continue
            stem, mate = match[1], match[2] or match[3] or "single"
            _name(stem)
            if selected and stem not in selected and f"{dataset}/{stem}" not in selected:
                continue
            files = groups.setdefault(stem, {})
            if mate in files:
                raise RuntimeError(f"Ambiguous FASTQ copies for {dataset}/{stem}: {files[mate]} and {path}")
            if files and any(previous.parent != path.parent for previous in files.values()):
                raise RuntimeError(f"FASTQ mates are in different folders: {dataset}/{stem}")
            files[mate] = path
        if not groups and not selected and ready_only:
            print(f"SKIP no completed FASTQ candidates: {folder}", flush=True)
            continue
        if not groups and not selected:
            raise RuntimeError(f"No FASTQs in requested dataset: {folder}")
        for stem, files in sorted(groups.items()):
            if ready_only and not selected:
                paired = set(files) in ({"1", "2"}, {"single"})
                incomplete = any(not p.stat().st_size or Path(str(p) + ".aria2").exists()
                                 or Path(str(p.resolve()) + ".aria2").exists() for p in files.values())
                if not paired or incomplete:
                    print(f"SKIP incomplete sample: {dataset}/{stem}", flush=True)
                    continue
            if set(files) == {"1", "2"}:
                r1, r2 = files["1"], files["2"]
            elif set(files) == {"single"}:
                r1, r2 = files["single"], None
            else:
                raise RuntimeError(f"Missing mate or mixed single/paired FASTQs: {dataset}/{stem}")
            paths = []
            for path in (r1, r2):
                if path is None:
                    paths.append(None)
                    continue
                real = path.resolve(strict=True)
                if not real.stat().st_size or Path(str(path) + ".aria2").exists() or Path(str(real) + ".aria2").exists():
                    raise RuntimeError(f"FASTQ empty or download incomplete: {path}")
                paths.append(str(real))
            if paths[1] and Path(paths[0]).samefile(paths[1]):
                raise RuntimeError(f"R1 and R2 resolve to the same file: {dataset}/{stem}")
            result.append({"sample_id": f"{dataset}/{stem}", "dataset": dataset,
                           "sample": stem, "r1": paths[0], "r2": paths[1]})
    for selector in selected:
        matches = [s for s in result if selector in (s["sample_id"], s["sample"])]
        if len(matches) != 1:
            raise RuntimeError(f"Selector {selector!r} matched {len(matches)} samples; use dataset/sample")
    if not result:
        raise RuntimeError("No samples selected")
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
        if saved != owner:
            raise RuntimeError(f"Stage belongs to different inputs: {path}")
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
    result = {}
    for directory, records in ((snapshot["db_dir"], snapshot["files"]),
                               (snapshot["taxonomy_dir"], snapshot["taxonomy_files"])):
        for name in records:
            path = (Path(directory) / name).resolve(strict=True)
            result[str(path)] = _runtime_state(path)
    return result


def _manifest(folder, exclude=()):
    result = {}
    for path in sorted(Path(folder).rglob("*")):
        name = path.relative_to(folder).as_posix()
        if path.is_symlink():
            raise RuntimeError(f"Symlink appeared in private output: {path}")
        if path.is_file() and name not in exclude:
            with path.open("rb") as handle:
                os.fsync(handle.fileno())
            result[name] = {"state": file_state(path), "runtime_state": _runtime_state(path),
                            "sha256": sha256(path)}
    return result


def _complete(folder, identity, verify_hashes=False):
    folder = Path(folder)
    if folder.is_symlink():
        raise RuntimeError(f"Refusing symlink output: {folder}")
    marker = folder / "complete.json"
    if not marker.is_file():
        return None
    saved = json.loads(marker.read_text())
    if saved.get("schema") != SCHEMA or saved.get("identity") != identity or saved.get("status") != "complete":
        raise RuntimeError(f"Output completion identity conflicts: {folder}")
    expected = saved.get("outputs", {})
    if not expected:
        raise RuntimeError(f"Empty output completion receipt: {folder}")
    actual = set()
    for path in folder.rglob("*"):
        if path.is_symlink():
            raise RuntimeError(f"Symlink appeared in completed output: {path}")
        if path.is_file():
            actual.add(path.relative_to(folder).as_posix())
    if actual != set(expected) | {"complete.json", OWNER}:
        raise RuntimeError(f"Completed output file inventory differs: {folder}")
    for name, item in expected.items():
        path = folder / name
        if Path(name).is_absolute() or ".." in Path(name).parts or file_state(path) != item["state"]:
            raise RuntimeError(f"Completed output changed: {path}")
        runtime_state = item.get("runtime_state")
        if runtime_state is not None and _runtime_state(path) != runtime_state:
            raise RuntimeError(f"Completed output changed: {path}")
        # Legacy output receipts have no local change guard; size alone cannot
        # justify reuse.  Validate their recorded content hashes instead.
        if (verify_hashes or runtime_state is None) and sha256(path) != item["sha256"]:
            raise RuntimeError(f"Completed output hash mismatch: {path}")
    return saved


def _context(cfg, sample, snapshot, limit):
    native = resolve_native(cfg)
    source = Path(__file__).resolve().parent
    return {"schema": SCHEMA, "sample": sample, "inputs": input_states(sample),
            "input_runtime_states": _input_runtime_states(sample),
            "database_runtime_states": _database_runtime_states(snapshot),
            "classifier_engine": _engine(cfg),
            "engine_parameters": kunpeng_query.parameters(cfg) if _engine(cfg) != "kraken" else {},
            "database_identity": snapshot["identity"], "database_files": snapshot["files"],
            "source_database_identity": snapshot.get("source_database_identity", snapshot["identity"]),
            "taxonomy_files": snapshot["taxonomy_files"],
            "db_dir": str(Path(snapshot["db_dir"]).resolve()),
            "taxonomy_dir": str(Path(snapshot["taxonomy_dir"]).resolve()),
            "native": {"path": str(native), "sha256": sha256(native)},
            "confidence": float(cfg.get("confidence", 0.1)),
            "minimum_hit_groups": int(cfg.get("minimum_hit_groups", 2)),
            "minimum_base_quality": int(cfg.get("minimum_base_quality", 0)),
            "output_format": cfg.get("output_format", "tar"), "limit": limit,
            "code": {name: sha256(source / name) for name in
                     ("hpc_query.py", "ordinal_query.py", "kunpeng_query.py", "species_split.py", "fastq_stream.py",
                      "taxonomy.py", "hpc_common.py", *(("partitioned_build.py",) if snapshot.get("layout") == "partitioned" else ()))}}


def _verify_sources(sample, snapshot, context):
    if input_states(sample) != context["inputs"]:
        raise RuntimeError("Original FASTQ files changed since this run began")
    if _input_runtime_states(sample) != context["input_runtime_states"]:
        raise RuntimeError("Original FASTQ runtime state changed since this run began")
    validate_snapshot(snapshot)
    if _database_runtime_states(snapshot) != context["database_runtime_states"]:
        raise RuntimeError("Database or taxonomy runtime state changed since this run began")
    if sha256(context["native"]["path"]) != context["native"]["sha256"]:
        raise RuntimeError("Native executable changed during this sample")
    source = Path(__file__).resolve().parent
    if any(sha256(source / name) != value for name, value in context["code"].items()):
        raise RuntimeError("Query code changed during this sample; do not edit a running release")


def _roots(cfg):
    out, work = Path(cfg["output_dir"]).resolve(), Path(cfg["work_dir"]).resolve()
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


def _publish_sample(cfg, sample, snapshot, limit):
    started = time.monotonic()
    _name(sample["dataset"])
    _name(sample["sample"])
    if sample["sample_id"] != f"{sample['dataset']}/{sample['sample']}":
        raise ValueError("sample_id must equal dataset/sample")
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")
    if cfg.get("output_format", "tar") not in ("tar", "directory"):
        raise ValueError("output_format must be tar or directory")
    validate_snapshot(snapshot)
    output_root, work_root = _roots(cfg)
    context = _context(cfg, sample, snapshot, limit)
    identity = fingerprint(context)
    if limit:
        target = output_root / "pilots" / sample["dataset"] / sample["sample"] / identity[:20]
    else:
        target = output_root / sample["dataset"] / sample["sample"]
    # Resolve only the configured root, rejecting aliases below it.
    if target.resolve() != target:
        raise RuntimeError(f"Symlink/alias in sample output path: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    lock = output_root / ".locks" / sample["dataset"] / (sample["sample"] + ".lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    with locked(lock):
        completed = _complete(target, identity, bool(cfg.get("verify_output_hashes", False)))
        if completed:
            _verify_sources(sample, snapshot, context)
            return {"status": "skipped", "sample_id": sample["sample_id"], "identity": identity,
                    "output": str(target), "summary": completed["summary"]}
        if target.exists() or target.is_symlink():
            raise RuntimeError(f"Output exists without a matching complete receipt: {target}")
        reserve = float(cfg.get("min_free_gib", 100))
        check_space(work_root, reserve)
        check_space(target.parent, reserve)
        owner = {"schema": SCHEMA, "identity": identity, "target": str(target)}
        work = work_root / "queries" / sample["dataset"] / sample["sample"] / identity
        if work.resolve() != work:
            raise RuntimeError(f"Symlink/alias in sample work path: {work}")
        _owned(work, owner)
        stage = target.with_name(target.name + ".partial." + identity[:16])
        if stage.exists() and (stage / "complete.json").is_file():
            _owned(stage, owner)
            recovered = _complete(stage, identity, bool(cfg.get("verify_output_hashes", False)))
            _verify_sources(sample, snapshot, context)
            stage.rename(target)
            _sync_directory(target.parent)
            _remove_owned(work, owner)
            return {"status": "recovered", "sample_id": sample["sample_id"], "identity": identity,
                    "output": str(target), "summary": recovered["summary"]}
        query_dir = work / "query"
        query_owner = dict(owner, kind="ordinal-query")
        query = None
        if (query_dir / "query.json").is_file():
            _owned(query_dir, query_owner)
            saved = json.loads((query_dir / "query.json").read_text())
            candidate = saved.get("query", {})
            if (saved.get("identity") == identity and saved.get("status") == "complete" and
                    Path(candidate.get("taxids", "")).resolve() == (query_dir / "taxids.bin").resolve() and
                    file_state(query_dir / "taxids.bin") == candidate.get("taxids_state") and
                    sha256(query_dir / "taxids.bin") == candidate.get("taxids_sha256")):
                query = candidate
        query_reused = query is not None
        if query is None:
            _owned(query_dir, query_owner, reset=True)
            query = query_once(cfg, sample, snapshot, query_dir, limit)
            _verify_sources(sample, snapshot, context)
            atomic_json(query_dir / "query.json", {"status": "complete", "identity": identity, "query": query})
        else:
            print(f"[{sample['sample_id']}] REUSE validated ordinal classifications", flush=True)
        _owned(stage, owner, reset=True)
        payload = stage / "species"
        summary = split_originals(cfg, sample, query, snapshot, payload)
        provenance = {"schema": SCHEMA, "identity": identity, "created_utc": utcnow(),
                      "context": context, "query": query,
                      "classifier": {"kraken": ("CPU mmap category hashes, per-minimizer union LCA before paired/confidence/hit-group decision"
                                                if snapshot.get("layout") == "partitioned" else "Kraken native CPU, mmap of one shared database"),
                                     "kunpeng_full": "Kun-peng native direct CPU, private fully loaded index",
                                     "kunpeng_memory": "Kun-peng native splitr/annotate/resolve CPU, verified zero-hit completion"}[_engine(cfg)],
                      "classify_threads": cfg.get("classify_threads", 8),
                      "scratch_policy": "private stages on configured shared work/output filesystems",
                      "prefix_note": "prefix is a throughput pilot, not a recall validation" if limit else ""}
        atomic_json(payload / "provenance.json", provenance)
        shutil.copyfile(query_dir / "query.log", payload / "classifier.log")
        atomic_json(payload / "manifest.json", {"files": _manifest(payload)})
        for name in ("summary.json", "counts.tsv", "provenance.json", "classifier.log"):
            shutil.copyfile(payload / name, stage / name)
        archive_started = time.monotonic()
        if cfg.get("output_format", "tar") == "tar":
            payload_bytes = sum(path.stat().st_size for path in payload.rglob("*") if path.is_file())
            payload_files = sum(1 for path in payload.rglob("*") if path.is_file())
            # Include tar/PAX headers and per-member padding, not only gzip sizes.
            check_space(stage, reserve + (payload_bytes + payload_files * 4096 + (1 << 20)) / 2**30)
            archive = stage / "species.tar"
            with tarfile.open(archive, "w", format=tarfile.PAX_FORMAT) as tar:
                for path in sorted(payload.rglob("*")):
                    if path.is_symlink():
                        raise RuntimeError("Symlink appeared during archive publication")
                    if path.is_file():
                        tar.add(path, arcname=path.relative_to(payload).as_posix(), recursive=False)
            with archive.open("rb") as handle:
                os.fsync(handle.fileno())
            # stage ownership was verified, and this exact payload was created above.
            shutil.rmtree(payload)
        else:
            for folder in sorted((p for p in payload.rglob("*") if p.is_dir()),
                                 key=lambda p: len(p.parts), reverse=True):
                _sync_directory(folder)
            _sync_directory(payload)
        archive_seconds = time.monotonic() - archive_started
        outputs = _manifest(stage, exclude=(OWNER, "complete.json"))
        _verify_sources(sample, snapshot, context)
        check_space(stage, reserve)
        summary.update(archive_seconds=archive_seconds, query_reused=query_reused,
                       query_seconds_this_run=0.0 if query_reused else query["seconds"],
                       total_seconds_before_commit=time.monotonic() - started,
                       artifact_bytes=sum(item["state"]["size"] for name, item in outputs.items()
                                          if name == "species.tar" or name.startswith("species/")),
                       timing_note="Outer summary includes archiving and validation; archive-internal summary precedes archiving")
        atomic_json(stage / "summary.json", summary)
        outputs["summary.json"] = {"state": file_state(stage / "summary.json"),
                                   "runtime_state": _runtime_state(stage / "summary.json"),
                                   "sha256": sha256(stage / "summary.json")}
        atomic_json(stage / "complete.json", {"schema": SCHEMA, "status": "complete",
                    "identity": identity, "completed_utc": utcnow(), "summary": summary, "outputs": outputs})
        if target.exists() or target.is_symlink():
            raise RuntimeError(f"Output appeared before publication: {target}")
        stage.rename(target)
        _sync_directory(target.parent)
        # A complete, published output is the prerequisite for deleting private work.
        _remove_owned(work, owner)
        print(f"[{sample['sample_id']}] PASS {summary['fragments']:,} fragments -> {target}", flush=True)
        return {"status": "complete", "sample_id": sample["sample_id"], "identity": identity,
                "output": str(target), "summary": summary}


def run_sample(cfg, sample, snapshot=None, limit=0):
    if snapshot is None:
        snapshot = database_snapshot(cfg)
    with _signals():
        return _publish_sample(cfg, sample, snapshot, limit)


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
