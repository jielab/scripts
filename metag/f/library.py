"""Permanent library inputs and completion receipts."""
import csv
import hashlib
import json
from pathlib import Path
from common import DIR_REF, DIR_LIB


def active_references():
    inventory = DIR_LIB/"reference/inventory.tsv"
    if not inventory.is_file():
        raise RuntimeError("Validated inventory missing; run ./1.library.sh --check-only")
    with inventory.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    seen = set()
    for row in rows:
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts or relative.name != row["filename"]:
            raise RuntimeError(f"Unsafe inventory path: {relative}")
        if relative.name in seen:
            raise RuntimeError(f"Duplicate reference filename: {relative.name}")
        seen.add(relative.name)
        row["path"] = str(DIR_REF/relative)
    if not rows:
        raise RuntimeError("Validated reference inventory is empty")
    return rows


def file_state(path):
    info = Path(path).stat()
    return [info.st_size, info.st_mtime_ns, info.st_ctime_ns]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def reference_identity():
    root = DIR_LIB/"reference"
    saved = json.loads((root/"check-ref.receipt.json").read_text())
    if saved.get("status") != "complete":
        raise RuntimeError("Reference check incomplete; run ./1.library.sh --check-only")
    return digest({"references": saved["references"], "metadata": saved["metadata"],
                   "mapping": hashlib.sha256((root/"assembly_summary.txt").read_bytes()).hexdigest()})


def completed_library(folder):
    path = Path(folder)/"build.json"
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Database completion receipt missing/invalid: {path}; run ./1.library.sh") from exc
    if value.get("status") != "complete":
        raise RuntimeError(f"Database build is incomplete: {path}; run ./1.library.sh")
    return value
