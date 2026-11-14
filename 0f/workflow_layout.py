"""Small, collision-safe directories for command plans and support files."""
from pathlib import Path
import errno
import os
import shutil
import tempfile


def new_run_dir(parent):
    """Create run, run2, ... without overwriting an active or completed plan."""
    parent = Path(parent).expanduser().resolve()
    parent.mkdir(parents=True, exist_ok=True)
    revision = 1
    while True:
        path = parent / ("run" if revision == 1 else f"run{revision}")
        try:
            path.mkdir()
            return path
        except FileExistsError:
            revision += 1


def publish_result(source, destination, *, immutable=False):
    """Expose an explicitly selected, completed result at a short stable path.

    Copy mutable work files so a later failed run cannot change public results.
    For explicitly immutable, receipt-bound sources, prefer a hard link (no
    duplicate data) and fall back to an atomic copy across filesystems.
    Call only after the owning scientific stage validates its completion record.
    Model caches, provenance and individual-level tables are not auto-discovered.
    """
    source, destination = Path(source), Path(destination)
    if not source.is_file():
        raise FileNotFoundError(source)
    if source.resolve() == destination.resolve():
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".publish-", dir=destination.parent)
    os.close(fd)
    try:
        os.unlink(temporary)
        if not immutable:
            shutil.copy2(source, temporary)
        else:
            try:
                os.link(source, temporary)
            except OSError as exc:
                if exc.errno not in {errno.EXDEV, errno.EPERM, errno.EACCES, errno.ENOTSUP}:
                    raise
                shutil.copy2(source, temporary)
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return destination


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent")
    print(new_run_dir(parser.parse_args().parent))


import fcntl
import json
import os
from pathlib import Path
import re
import tempfile
from datetime import datetime, timezone

from workflow_lock import flock


def run_name():
    return datetime.now(timezone.utc).strftime("run-%Y%m%dT%H%M%S-%fZ") + f"-p{os.getpid()}"


def named_directory(parent, identity, label):
    """Map a complete fingerprint to a stable readable, collision-safe name.

    Reserve names under a shared lock without creating the result directory:
    callers can still publish a completed staging directory atomically.
    """
    parent = Path(parent)
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*", label):
        raise ValueError("Invalid output directory label")
    qc = parent / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    registry = qc / "directory_names.json"
    with (qc / "directory_names.lock").open("a") as lock:
        flock(lock, fcntl.LOCK_EX)
        entries = json.loads(registry.read_text()) if registry.exists() else {}
        key = label + ":" + identity
        if key not in entries:
            used = set(entries.values()) | {p.name for p in parent.iterdir()}
            number = 1
            while f"{label}-{number:03d}" in used:
                number += 1
            entries[key] = f"{label}-{number:03d}"
            with tempfile.NamedTemporaryFile(mode="w", dir=qc, delete=False) as stream:
                json.dump(entries, stream, indent=2, sort_keys=True)
                temporary = stream.name
            os.replace(temporary, registry)
        name = entries[key]
        if not re.fullmatch(re.escape(label) + r"-\d+", name):
            raise ValueError("Invalid output directory registry entry")
        return parent / name
