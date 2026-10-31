"""Durable receipts, verified libraries, helper build, and process supervision."""
from __future__ import annotations
from contextlib import contextmanager
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from common import ROOT, DIR_LIB, DIR_REF, DIR_RAW, THREADS
from library import completed_library, file_state
from metacache_build import complete_part

VERSION = "metag-compact-v2.0"


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def sync_directory(path):
    try:
        fd = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
    except OSError as exc:
        if exc.errno not in (errno.EINVAL, errno.ENOTSUP):
            raise


def durable_unlink(path):
    path = Path(path)
    path.unlink()
    sync_directory(path.parent)


def durable_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".writing")
    with temp.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)
    sync_directory(path.parent)


@contextmanager
def locked(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another process owns this run/output: {path}") from exc
        yield


def space_guard(path, reserve, extra=0):
    path = Path(path)
    while not path.exists():
        if path.parent == path:
            raise RuntimeError(f"Cannot locate filesystem for {path}")
        path = path.parent
    free = shutil.disk_usage(path).free
    if free < reserve + extra:
        raise RuntimeError(f"Insufficient free space on {path}: {free / 2**30:.1f} GiB; "
                           f"need {(reserve+extra) / 2**30:.1f} GiB including reserve")


def safe_destination(path):
    path = Path(path).resolve()
    for protected in (DIR_LIB, DIR_REF, DIR_RAW):
        root = protected.resolve()
        if path == root or root in path.parents or path in root.parents:
            raise RuntimeError(f"Work/output must be separate from source/library paths: {path}, protected {root}")
    return path


def library_snapshot():
    root = DIR_LIB / "metacache/db"
    saved = completed_library(root)
    if not saved.get("identity") or not saved.get("results"):
        raise RuntimeError(f"No completed MetaCache libraries in {root}")
    parts = []
    seen = set()
    for row in saved["results"]:
        db = Path(row["db"]).resolve()
        identity = row.get("identity", saved["identity"])
        if str(db) in seen or not complete_part(db, identity):
            raise RuntimeError(f"Duplicate, missing or changed completed database: {db}")
        seen.add(str(db))
        receipt = json.loads(Path(str(db) + ".complete.json").read_text())
        parts.append({"db": str(db), "identity": identity, "outputs": receipt["outputs"]})
    tax = DIR_LIB / "metacache/taxonomy"
    taxonomy = {name: file_state(tax / name) for name in ("nodes.dmp", "names.dmp", "merged.dmp")}
    return {"identity": saved["identity"], "parts": parts, "taxonomy": taxonomy, "taxonomy_path": str(tax)}


def choose_backend(request="auto"):
    directory = Path(os.environ.get("METAG_TOOLS", Path.home() / ".cache/metag/tools")) / "metacache"
    gpu, cpu = directory / "metacache_gpu", directory / "metacache"
    available = False
    if os.access(gpu, os.X_OK) and shutil.which("nvidia-smi"):
        try:
            available = subprocess.run(["nvidia-smi", "-L"], stdout=subprocess.DEVNULL,
                                       stderr=subprocess.DEVNULL, timeout=15).returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            pass
    if request == "gpu" and not available:
        raise RuntimeError("GPU backend requested, but metacache_gpu / NVIDIA GPU is unavailable")
    if request != "cpu" and available:
        return gpu.resolve(), "gpu"
    if not os.access(cpu, os.X_OK):
        raise RuntimeError(f"MetaCache executable not found: {cpu}; run ./install.sh")
    return cpu.resolve(), "cpu"


def build_helper():
    source = ROOT / "f/compact_merge.cpp"
    base = Path(os.environ.get("METAG_TOOLS", Path.home() / ".cache/metag/tools")) / "compact"
    base.mkdir(parents=True, exist_ok=True)
    binary, receipt = base / "compact_merge", base / "build.json"
    identity = sha256(source)
    with locked(base / "build.lock"):
        if receipt.exists() and os.access(binary, os.X_OK):
            data = json.loads(receipt.read_text())
            if data.get("source_sha256") == identity and data.get("binary") == file_state(binary):
                return binary
        compiler = shutil.which("g++")
        if not compiler:
            raise RuntimeError("g++ is required for the compact helper (the existing MetaCache installer installs it)")
        temp = binary.with_name(binary.name + ".partial")
        temp.unlink(missing_ok=True)
        print("Compile compact merge helper (existing MetaCache databases are reused)", flush=True)
        subprocess.run([compiler, "-O3", "-std=c++17", "-DNDEBUG", str(source), "-o", str(temp)], check=True)
        temp.replace(binary)
        durable_json(receipt, {"source_sha256": identity, "binary": file_state(binary)})
    return binary


def helper_call(helper, *args):
    result = subprocess.run([str(helper), *map(str, args)], text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"Compact helper failed ({args[0]}): {result.stderr.strip()}")
    try:
        return json.loads(result.stdout)
    except ValueError as exc:
        raise RuntimeError(f"Invalid compact helper completion: {result.stdout[:300]}") from exc


def _stop(process):
    if process is None or process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=5)
    except ProcessLookupError:
        pass
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def query_to_delta(exe, backend, helper, db, prep, state, delta, log, batch_size, reserve):
    """Separate mapping FD from logs; no FIFO filesystem support is needed on WSL."""
    delta, log = Path(delta), Path(log)
    delta.unlink(missing_ok=True)  # caller recovered any committed pending delta first
    consumer_log = log.with_name(log.stem + ".merge.log")
    summary_path = log.with_name(log.stem + ".merge.json")
    read_fd, write_fd = os.pipe()
    producer = consumer = None
    started = time.monotonic()
    command = [str(exe), "query", str(db), prep["query_r1"]]
    if prep.get("query_r2"):
        command += [prep["query_r2"], "-pairfiles"]
    command += ["-threads", str(THREADS), "-batch-size", str(batch_size),
                "-tophits", "-queryids", "-lowest", "species", "-highest", "domain",
                "-maxcand", "2", "-hitmin", "5", "-hitdiff", "100",
                "-taxids-only", "-separate-cols", "-out", f"/proc/self/fd/{write_fd}"]
    try:
        with log.open("w") as native_log, consumer_log.open("w") as merge_log, summary_path.open("w") as summary:
            native_log.write("COMMAND " + json.dumps(command) + "\n")
            native_log.flush()
            consumer = subprocess.Popen([str(helper), "update", "--state", str(state),
                "--input", "-", "--delta", str(delta), "--paired", str(int(prep["paired"])),
                "--headers", prep["candidate_headers"]], stdin=read_fd,
                stdout=summary, stderr=merge_log, start_new_session=True)
            os.close(read_fd); read_fd = -1
            producer = subprocess.Popen(command, stdout=native_log, stderr=subprocess.STDOUT,
                                        pass_fds=(write_fd,), start_new_session=True)
            os.close(write_fd); write_fd = -1
            checked = 0.0
            while producer.poll() is None or consumer.poll() is None:
                now = time.monotonic()
                if now - checked >= 5:
                    space_guard(delta.parent, reserve)
                    checked = now
                if consumer.poll() not in (None, 0):
                    _stop(producer)
                    break
                if producer.poll() not in (None, 0):
                    _stop(consumer)
                    break
                time.sleep(0.1)
            pcode, ccode = producer.wait(), consumer.wait()
        if pcode or ccode:
            details = consumer_log.read_text(errors="replace")[-1500:]
            raise RuntimeError(f"Query/merge failed for {db}: MetaCache={pcode}, merge={ccode}. "
                               f"No checkpoint committed. See {log} and {consumer_log}. {details}")
        try:
            details = json.loads(summary_path.read_text())
        except ValueError as exc:
            raise RuntimeError(f"Missing compact query completion: {summary_path}") from exc
        return {"db": str(db), "backend": backend, "batch_size": batch_size,
                "seconds": time.monotonic()-started, "delta_bytes": delta.stat().st_size,
                "native_log": str(log), "compact": details}
    finally:
        if read_fd >= 0: os.close(read_fd)
        if write_fd >= 0: os.close(write_fd)
        _stop(producer); _stop(consumer)
