#!/usr/bin/env python3
"""Portable dispatcher for the retained, explicitly labelled PRSformer adapters."""
from __future__ import annotations

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from grid_paths import GRID, analysis_default, read_config, relocate
from resume import atomic_json, run_lock


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] not in {"prsformer", "prsformer2"}:
        raise ValueError("Use prsformer or prsformer2")
    module, argv = argv[0], argv[1:]
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--analysis-root", "--output-root", default=analysis_default())
    p.add_argument("--config", default=str(GRID / "config/abm.json"))
    known, forwarded = p.parse_known_args(argv)
    root = Path(known.analysis_root)
    if "--help" in forwarded or "-h" in forwarded:
        subprocess.run([sys.executable, str(GRID / "f/twin.py"), module, "--help"], check=True)
        return
    config = read_config(known.config, root)
    block = read_config(config["block_config"], root)
    if "--dry-run" in forwarded:
        traits = forwarded[forwarded.index("--traits")+1] if "--traits" in forwarded else "height,ldl,t2dm"
        print(json.dumps(dict(module=module, traits=traits, analysis_root=str(root),
            training_mode="independent_trait_adaptation" if module == "prsformer" else "exploratory_extension",
            shared_cache=block["prsformer_cache"], owns_test_evaluation=True)))
        return
    folder = root / "cache/benchmark_config"
    folder.mkdir(parents=True, exist_ok=True)
    with run_lock(folder):
        block_digest = hashlib.sha256(json.dumps(block, sort_keys=True).encode()).hexdigest()[:20]
        block_path = folder / ("blocks." + block_digest + ".json")
        if not block_path.is_file():
            atomic_json(block_path, block)
        config["block_config"] = str(block_path)
        digest = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()[:20]
        path = folder / ("prsformer." + digest + ".json")
        if not path.is_file():
            atomic_json(path, config)
    if module == "prsformer2" and "--upstream-dir" not in forwarded:
        forwarded += ["--upstream-dir", os.environ.get("PRSFORMER_UPSTREAM_DIR", relocate(configured_path('DIR_SOFTWARE', 'PRSformer')))]
    cache = Path(block["prsformer_cache"])
    if "--check" not in forwarded and not all((cache / name).is_file() for name in
            ("prepare.json", "data.tsv.gz", "variants.tsv.gz")):
        from shared_cohort import ensure_inputs
        from types import SimpleNamespace
        threads = int(forwarded[forwarded.index("--threads")+1]) if "--threads" in forwarded else 4
        ensure_inputs(SimpleNamespace(analysis_root=str(root), config=str(path), block_config=str(block_path),
            traits="height,ldl,t2dm", threads=threads, complete_features=False, dry_run=False), config)
    command = [sys.executable, str(GRID / "f/twin.py"), module, "--analysis-root", str(root),
               "--config", str(path), *forwarded]
    subprocess.run(command, check=True, env=dict(os.environ, GRID_WORKFLOW_ANALYSIS_ROOT=str(root)))
    if "--check" not in forwarded:
        from prsformer_benchmark import publish
        traits = forwarded[forwarded.index("--traits")+1] if "--traits" in forwarded else "height,ldl,t2dm"
        for trait in traits.split(","):
            publish(root, trait, module, cache=block["prsformer_cache"])


if __name__ == "__main__":
    main()
