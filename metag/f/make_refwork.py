#!/usr/bin/env python3
from pathlib import Path
import csv, os, shutil
from common import DIR_REF, DIR_LIB, assert_reference_policy

def main(sources_checked=False):
    assert_reference_policy()
    inv = DIR_LIB/"reference"/"inventory.tsv"
    mapping = DIR_LIB/"reference"/"assembly_summary.txt"
    if not inv.is_file() or not mapping.is_file():
        raise SystemExit("Reference mapping missing; run ./1.library.sh")
    root = DIR_LIB/"metacache"/"genomes"
    if root.exists() and not root.is_dir():
        raise SystemExit(f"Invalid {root}")
    root.mkdir(parents=True, exist_ok=True)

    # This inventory is the validated active set, including explicit exclusions.
    # Rescanning DIR_REF here would silently reintroduce excluded references.
    with inv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError("Validated reference inventory is empty")
    refs = []
    for row in rows:
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise RuntimeError(f"Invalid reference inventory path: {relative}")
        p = DIR_REF/relative
        if p.name != row["filename"] or (not sources_checked and not p.is_file()):
            raise RuntimeError(f"Invalid/missing inventory reference: {p}")
        refs.append(p)
    print(f"[refwork] Prepare {len(refs):,} reference links at {root}", flush=True)
    with os.scandir(root) as entries:
        existing = {entry.name: entry.is_symlink() for entry in entries}
    expected=set()
    for number, p in enumerate(refs, 1):
        if p.name in expected:
            raise RuntimeError(f"Duplicate reference filename: {p.name}")
        dst=root/p.name
        expected.add(dst.name)
        if p.name in existing:
            target = os.readlink(dst) if existing[p.name] else None
            target_path = Path(os.path.abspath(dst.parent/target)) if target is not None else None
            if target_path != p.absolute() and (target_path is None or dst.resolve() != p.resolve()):
                raise SystemExit(f"Refwork collision: {dst}")
        else:
            dst.symlink_to(p.absolute())
        if number % 5000 == 0:
            print(f"[refwork] {number:,}/{len(refs):,} links ready", flush=True)

    # Remove only stale symlinks created in this derived refwork, NEVER source FASTA.
    for name, is_symlink in existing.items():
        if name=="assembly_summary.txt": continue
        if is_symlink and name not in expected:
            (root/name).unlink()
    shutil.copyfile(mapping, root/"assembly_summary.txt")
    print(f"PASS: {len(expected):,} symlinked references at {root}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
