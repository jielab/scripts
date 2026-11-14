"""Simple complete-output reuse for optional analysis commands, without hashes."""
from contextlib import ExitStack
from functools import wraps
import json
from pathlib import Path
import shutil

from common import atomic_json, locked

_ACTIVE = []
_OWNER = ".metag-result.json"


def replace_bool(value):
    if value not in ("TRUE", "FALSE"):
        raise ValueError("--replace must be TRUE or FALSE")
    return value == "TRUE"


def managed_outputs(function):
    """Hold result-directory locks until the CLI finishes, including exceptions."""
    @wraps(function)
    def wrapper(*args, **kwargs):
        with ExitStack() as stack:
            _ACTIVE.append(stack)
            try:
                return function(*args, **kwargs)
            finally:
                _ACTIVE.pop()
    return wrapper


def prepare(path, *, replace=True, cfg=None, protected=()):
    """Return True to compute, False to reuse a complete current-format result."""
    path = Path(path).expanduser().absolute()
    if path.is_symlink() or path.resolve() != path:
        raise ValueError(f"Refusing a symlink/alias result directory: {path}")
    if not _ACTIVE:
        raise RuntimeError("Analysis output must be managed by managed_outputs")
    cfg = cfg or {}
    for key in ("reads_raw", "reference_raw", "metadata_dir", "library_dir", "output_dir", "tools_dir"):
        if cfg.get(key):
            source = Path(cfg[key]).resolve()
            if path == source or path in source.parents or source in path.parents:
                raise ValueError(f"Analysis output overlaps protected {key}: {path}")
    for source in protected:
        if source:
            source = Path(source).expanduser().resolve()
            if path == source or path in source.parents:
                raise ValueError(f"Analysis output contains its own input: {source}")
    if path == Path(path.anchor):
        raise ValueError("Result directory cannot be a filesystem root")
    path.parent.mkdir(parents=True, exist_ok=True)
    _ACTIVE[-1].enter_context(locked(path.parent / ("." + path.name + ".lock")))
    marker = path / "complete.json"
    if not replace and marker.is_file():
        try:
            saved = json.loads(marker.read_text())
            records = saved.get("files", {})
            complete = saved.get("status") == "complete" and bool(records)
            for name, size in records.items():
                rel = Path(name)
                file = path / rel
                if rel.is_absolute() or ".." in rel.parts or file.is_symlink() or not file.is_file() or file.stat().st_size != size:
                    complete = False
                    break
            if complete:
                print(f"REUSE complete analysis: {path}", flush=True)
                return False
        except (OSError, ValueError, TypeError):
            pass
    if path.exists():
        if not path.is_dir():
            raise ValueError(f"Expected an analysis directory: {path}")
        if any(path.iterdir()):
            owner = path / _OWNER
            try:
                owned = json.loads(owner.read_text()).get("format") == "metag-analysis-v1"
            except (OSError, ValueError):
                owned = False
            if not owned:
                raise ValueError(f"Nonempty output is not a current METAG analysis directory; choose a dedicated empty --out: {path}")
            shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    atomic_json(path / _OWNER, {"format": "metag-analysis-v1"})
    return True


def finish(path):
    """Commit only after all required work succeeded; store sizes, not checksums."""
    path = Path(path)
    marker = path / "complete.json"
    saved = json.loads(marker.read_text()) if marker.is_file() else {}
    records = {}
    for file in path.rglob("*"):
        if file.is_symlink():
            raise RuntimeError(f"Symlink in analysis output: {file}")
        if file.is_file() and file.name not in (_OWNER, "complete.json"):
            records[file.relative_to(path).as_posix()] = file.stat().st_size
    if not records:
        raise RuntimeError("Analysis did not produce any files")
    saved.update(status="complete", files=records)
    atomic_json(marker, saved)
