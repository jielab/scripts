"""Small shared helpers. All learned objects are local, trusted research artifacts."""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import json
import os
import time
import numpy as np

VERSION = "4.0.0"

def words(value):
    return [s.strip() for s in (value or "").split(",") if s.strip()]

def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value

def dump(path, value):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(json_safe(value), ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)

def digest(value):
    return hashlib.sha256(json.dumps(json_safe(value), sort_keys=True).encode()).hexdigest()

def log(action, step, detail=""):
    print(f"[PANOME] {action} {step}" + (f" | {detail}" if detail else ""), flush=True)

@contextmanager
def stage_log(step):
    start = time.monotonic()
    log("START", step)
    yield
    log("DONE", step, f"{time.monotonic() - start:.1f}s")

def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    ix = ok[np.argsort(p[ok])]
    q[ix] = np.minimum(1, np.minimum.accumulate(
        (p[ix] * len(ix) / np.arange(1, len(ix) + 1))[::-1])[::-1])
    return q

def fingerprints(paths, full=False):
    rows = []
    for path in sorted(set(str(Path(p).resolve()) for p in paths if p)):
        p = Path(path)
        if not p.is_file():
            raise FileNotFoundError(p)
        stat = p.stat()
        row = dict(path=path, bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
        if full or stat.st_size < 2_000_000:
            h = hashlib.sha256()
            with p.open("rb") as handle:
                for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                    h.update(block)
            row["sha256"] = h.hexdigest()
        rows.append(row)
    return rows

@contextmanager
def run_lock(root):
    path = root / ".lock"
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f"Run is locked: {path}. Inspect its PID before removing a stale lock.") from exc
    with os.fdopen(fd, "w") as f:
        f.write(str(os.getpid()))
    try:
        yield
    finally:
        path.unlink(missing_ok=True)
