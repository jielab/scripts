"""A result tree belongs to exactly one current LE8 source release."""
from pathlib import Path
import hashlib
import json
import os
import tempfile
from file_lock import file_lock

SOURCE = Path(__file__).resolve().parents[1]
SCHEMA = "le8-current-only-1"


def source_identity(source=SOURCE):
    digest = hashlib.sha256()
    for p in sorted(source.rglob("*")):
        rel = p.relative_to(source)
        if p.is_file() and p.suffix in {".py", ".R", ".sh", ".json", ".r"} and not set(rel.parts) & {"tests", "__pycache__", ".git"}:
            digest.update(str(rel).replace("\\", "/").encode())
            digest.update(b"\0")
            digest.update(p.read_bytes())
    return dict(schema=SCHEMA, source_sha256=digest.hexdigest())


def require_current(root, *, create=False, identity=None):
    root = Path(root).resolve()
    expected = source_identity() if identity is None else identity
    marker = root / "RELEASE.json"
    if not root.exists() and not create:
        raise ValueError("No current LE8 analysis: " + str(root))
    root.mkdir(parents=True, exist_ok=True)
    with file_lock(root / ".release.lock"):
        if marker.is_file():
            if json.loads(marker.read_text()) != expected:
                raise ValueError("LE8 source/result release mismatch; start in a clean output directory. No migration or repair is supported.")
            return
        if not create or any(p.name != ".release.lock" for p in root.iterdir()):
            raise ValueError("Missing current RELEASE.json; existing outputs are not accepted. Start in a clean output directory.")
        fd, name = tempfile.mkstemp(dir=root, prefix=".release-")
        try:
            with os.fdopen(fd, "w") as f:
                json.dump(expected, f, sort_keys=True)
            os.replace(name, marker)
        finally:
            if os.path.exists(name):
                os.unlink(name)
