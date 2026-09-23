"""Small shared helpers. All learned objects are local, trusted research artifacts."""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import json
import os
import shutil
import time
import numpy as np

VERSION = "5.1.0"

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


def run_complete(root, train_only=False):
    """A frozen model alone is not a completed run: audits may still fail."""
    root = Path(root)
    required = ["manifest.json", "MODEL_FROZEN.json", "model_bundle.joblib"]
    marker = "TRAIN_DONE.json" if train_only else "DONE.json"
    if train_only and (root/"DONE.json").is_file():
        return run_complete(root)
    if not train_only:
        required += ["REPORT.md", "test_metrics.csv"]
    if any(not (root/name).is_file() or (root/name).stat().st_size == 0
           for name in required + [marker]):
        return False
    try:
        record = json.loads((root/marker).read_text())
        return isinstance(record, dict) and bool(record.get("version"))
    except (ValueError, OSError):
        return False


def prepare_run_directory(root, resume=False, replace=False, train_only=False):
    """Called under run_lock; keep that lock held while clearing old outputs."""
    root = Path(root)
    entries = [p for p in root.iterdir() if p.name != ".lock"]
    if not entries:
        return True
    manifest = root/"manifest.json"
    if not manifest.is_file():
        raise ValueError("Refusing to overwrite a nonempty directory without a Panome manifest")
    try:
        record = json.loads(manifest.read_text())
    except (ValueError, OSError) as exc:
        raise ValueError("Refusing to overwrite a directory with an unreadable Panome manifest") from exc
    if not isinstance(record, dict) or not {"version", "config", "signature"} <= record.keys():
        raise ValueError("Refusing to overwrite a directory without a valid Panome manifest")
    if resume:
        return True
    if not replace and run_complete(root, train_only):
        log("SKIP", "completed", str(root))
        return False
    log("RESTART", "replace" if replace else "incomplete_or_failed", str(root))
    for path in entries:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        else:
            path.unlink()
    return True
