"""Resolve only independent benchmark outputs; joint layouts are never reused."""
import json
from pathlib import Path

TRAITS = ("height", "ldl", "t2dm")


def prsformer_directory(root, trait):
    return Path(root) / "prsformer" / trait


def prsformer_publication(root, trait):
    base = prsformer_directory(root, trait)
    return base / "split.rds", base / "predictions.rds"


def prsformer2_current(root, trait):
    base = Path(root) / "prsformer2" / trait
    pointer = base / "CURRENT.json"
    if not pointer.is_file():
        return None
    name = json.loads(pointer.read_text())["run"]
    if not isinstance(name, str) or name in ("", ".", "..") or Path(name).name != name:
        raise ValueError("PRSformer2 current run must stay inside its trait directory")
    run = base / name
    if run.is_symlink():
        raise ValueError("PRSformer2 current run cannot link to another model directory")
    return run


def benchmark_runs(root):
    """Include legacy and all independent runs when protecting unfinished work."""
    root = Path(root)
    runs = [root / "prsformer/run"]
    runs.extend(root / "prsformer" / trait / "run" for trait in TRAITS)
    runs.extend((root / "prsformer2").glob("*/_LOCAL_ONLY"))
    runs.extend((root / "prsformer2").glob("*/*/_LOCAL_ONLY"))
    return runs
