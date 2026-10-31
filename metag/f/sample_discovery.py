"""Reject ambiguous FASTQ layouts before the public split workflow starts.

This validation leaves sample names and classification cache identities unchanged.
Mixed singleton/pair data can be processed with explicit ``classify`` inputs and
distinct names; the batch workflow must never silently omit one of those files.
"""
from pathlib import Path

from common import FASTQ_RE


def validate_dataset(folder):
    folder = Path(folder)
    if not folder.is_dir():
        raise RuntimeError(f"Dataset not found: {folder}")
    groups = {}
    for path in sorted(folder.rglob("*")):
        if not path.is_file():
            continue
        match = FASTQ_RE.match(path.name)
        if not match:
            continue
        stem = match.group(1)
        mate = match.group(2) or match.group(3) or "single"
        key = (path.parent.relative_to(folder), stem)
        files = groups.setdefault(key, {})
        if mate in files:
            raise RuntimeError(
                f"Ambiguous FASTQ files for {stem!r}, mate {mate}: "
                f"{files[mate]} and {path}. The batch splitter would select only "
                "one. Choose the intended file with --r1/--r2 and a distinct "
                "--name, or keep only one intended input copy in this dataset."
            )
        files[mate] = path
    if not groups:
        raise RuntimeError(f"No FASTQ files found in dataset: {folder}")
    names = {}
    runs = []
    for (relative, stem), files in sorted(groups.items()):
        if stem in names:
            raise RuntimeError(
                f"Repeated sample stem {stem!r} in {names[stem]} and "
                f"{folder / relative}. Both would use the same split output. "
                "Process each input set with explicit --r1/--r2 and a distinct "
                "--name, or use separate dataset folders."
            )
        names[stem] = folder / relative
        mates = set(files)
        if mates == {"1", "2"}:
            runs.append((stem, files["1"], files["2"]))
        elif mates == {"single"}:
            runs.append((stem, files["single"], None))
        elif "single" in mates:
            raise RuntimeError(
                f"Mixed single/paired FASTQ files for {stem!r}: "
                + ", ".join(str(files[key]) for key in sorted(files))
                + ". The batch splitter would omit at least one file. "
                "For an SRA split-3 layout, process the pair with --r1/--r2 "
                "and process the singleton separately with --r1 and a distinct "
                "--name (for example <run>__single)."
            )
        else:
            raise RuntimeError(f"Unpaired FASTQ: {next(iter(files.values()))}")
    return runs


def validate_split_datasets(raw, dat_list=""):
    """Validate every requested dataset before any sample can be classified."""
    raw = Path(raw)
    wanted = [value.strip() for value in dat_list.split(",") if value.strip()]
    datasets = wanted or sorted(
        path.name for path in raw.iterdir() if path.is_dir() and not path.is_symlink()
    )
    if not datasets:
        raise RuntimeError(f"No datasets under {raw}")
    return [(name, validate_dataset(raw / name)) for name in datasets]
