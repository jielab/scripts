#!/usr/bin/env python3
"""Receipts for native stages validated in this invocation only.

Receipts live in a dispatcher-owned temporary directory. They never substitute
for input validation on the next invocation. Result hashes also prevent a
changed or missing upstream result from being skipped between transactions.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path


RESULTS = {
    "c1_correlate": ["c1_correlate/c1.res.rds"],
    "c2_cause": ["c2_cause/c2.res.rds"],
    "c3_coloc": ["c3_coloc/c3.res.rds"],
    "c4_connect": ["c4_connect/" + name for name in
                   ("c4.res.rds", "c4.interactions.res.rds", "c4.nonlin.res.rds", "c4.penalty.res.rds")],
    "c4_panel_validation": ["c4_connect/c4.validation.rds"],
}


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def snapshot(root, trait, layer, module):
    paths = [root / trait / layer / name for name in RESULTS[module]]
    if any(not p.is_file() or not p.stat().st_size for p in paths):
        raise ValueError(f"Missing result for {trait}/{layer}/{module}")
    return {str(p.relative_to(root)): digest(p) for p in paths}


def receipts(action, root, traits, layers, modules):
    directory = os.getenv("LE8_RUN_STATE_DIR")
    if not directory:
        return []
    directory = Path(directory)
    identity = str(Path(os.getenv("LE8_PUBLISHED_ROOT", root)).resolve())
    found = []
    for trait in traits:
        for layer in layers:
            for module in modules:
                if module not in RESULTS:
                    continue
                key = f"{trait}|{module}|{layer}"
                file = directory / (hashlib.sha256((identity + key).encode()).hexdigest() + ".json")
                if action == "read" and not file.is_file():
                    continue
                try:
                    current = dict(root=identity, key=key, outputs=snapshot(root, trait, layer, module))
                    if action == "record":
                        temporary = file.with_suffix(".tmp")
                        temporary.write_text(json.dumps(current, sort_keys=True))
                        temporary.replace(file)
                    elif json.loads(file.read_text()) == current:
                        found.append(key)
                except (ValueError, OSError):
                    if action == "record":
                        raise
                    # Invalid receipts require ordinary dependency/input validation.
    return found


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=["read", "record"])
    p.add_argument("root", type=Path)
    p.add_argument("traits")
    p.add_argument("layers")
    p.add_argument("modules")
    a = p.parse_args()
    for key in receipts(a.action, a.root, a.traits.split(","), a.layers.split(","), a.modules.split(",")):
        print(key)
