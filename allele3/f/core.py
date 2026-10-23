"""Small, dependency-free utilities for the local allele3 workflow."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

VERSION = "3.2.0-final-20261022"


def open_text(path, mode="rt"):
    path = Path(path)
    if path.name.endswith((".gz", ".bgz", ".bgzip")):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8", newline="")


def dump_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(value, f, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")
    os.replace(tmp, path)


def fingerprint(path, sample_bytes=65536):
    """Detect ordinary input changes without hashing hundreds of GB twice.

    This is deliberately labelled a metadata/edge fingerprint, NOT a full-file
    cryptographic integrity check. Download checksums remain the user's source
    integrity check. Small files are hashed completely.
    """
    p = Path(path).resolve()
    s = p.stat()
    h = hashlib.sha256()
    with p.open("rb") as f:
        h.update(f.read(sample_bytes))
        if s.st_size > sample_bytes:
            f.seek(max(sample_bytes, s.st_size - sample_bytes))
            h.update(f.read(sample_bytes))
    return {"path": str(p), "size": s.st_size, "mtime_ns": s.st_mtime_ns,
            "edge_sha256": h.hexdigest()}


def stable_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def code_hash():
    h = hashlib.sha256()
    for path in sorted(Path(__file__).parent.glob("*.py")):
        h.update(path.name.encode())
        h.update(path.read_bytes())
    return h.hexdigest()


def tsv_writer(path, fields):
    """Return a handle/writer; caller closes handle before committing a shard."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.name.endswith((".gz", ".bgz")):
        f = gzip.open(p, "wt", encoding="utf-8", newline="", compresslevel=2)
    else:
        f = p.open("w", encoding="utf-8", newline="")
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n",
                            extrasaction="raise", restval=".")
    writer.writeheader()
    return f, writer


def fmt(value):
    if value is None:
        return "."
    if isinstance(value, float):
        return format(value, ".10g")
    return value
