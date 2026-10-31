"""Bind sample completion to input reads, library versions and output files."""
import hashlib
import json
from pathlib import Path

from common import ROOT, DIR_LIB, SOURMASH_K, SOURMASH_SCALED
from library import completed_library, file_state


def context(r1, r2):
    reads = []
    for path in (r1, r2):
        if path is not None:
            path = Path(path).resolve()
            if not path.is_file():
                raise RuntimeError(f"Input FASTQ missing: {path}")
            reads.append({"path": str(path), "state": file_state(path)})
    sm = completed_library(DIR_LIB/"sourmash")
    mc = completed_library(DIR_LIB/"metacache/db")
    if not sm.get("identity") or not mc.get("identity"):
        raise RuntimeError("Library identity missing; run ./1.library.sh")
    return {"schema": 1, "reads": reads,
            "libraries": {"sourmash": sm["identity"], "metacache": mc["identity"]},
            "taxonomy": {name: file_state(DIR_LIB/"metacache/taxonomy"/name)
                         for name in ("nodes.dmp", "names.dmp", "merged.dmp")},
            "params": {"k": SOURMASH_K, "scaled": SOURMASH_SCALED},
            "code": {name: hashlib.sha256((ROOT/"f"/name).read_bytes()).hexdigest()
                     for name in ("pipeline.py", "sourmash_screen.py", "metacache_query.py",
                                  "split_metacache.py", "sample_state.py", "sample_report.py", "reference_evidence.py")}}


def output_state(path):
    info = Path(path).stat()
    return [info.st_size, info.st_mtime_ns]


def is_complete(target, identity):
    try:
        saved = json.loads((target/"split.done").read_text())
        if not isinstance(saved, dict) or saved.get("input_identity") != identity:
            return False
        if not isinstance(saved.get("outputs"), dict) or "counts.tsv" not in saved["outputs"]:
            return False
        for relative, expected in saved["outputs"].items():
            path = Path(relative)
            if path.is_absolute() or ".." in path.parts or output_state(target/path) != expected:
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False
