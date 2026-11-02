"""Relocate repository defaults without rewriting explicit user paths.

The WSL paths in historical JSON files are defaults, not a requirement to mount
those drives on a cluster. The common workflow helper owns platform detection;
these overrides also make direct Python scientific workers portable.
"""
from __future__ import annotations
import importlib.util
import json
import os
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
GRID = REPOSITORY / "grid"


def analysis_default():
    return os.environ.get("GRID_OUTPUT_ROOT", str(Path(os.environ["ANALYSIS_ROOT"]) / "grid")
                          if os.environ.get("ANALYSIS_ROOT") else "/mnt/d/analysis/grid")


def relocate(value, analysis_root=None):
    helper = REPOSITORY / "0f/workflow_paths.py"
    platform = lambda item: item
    if helper.is_file():
        spec = importlib.util.spec_from_file_location("grid_workflow_paths", helper)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        platform = module.relocate_defaults
    mappings = [("/mnt/d/scripts", str(REPOSITORY)),
                ("/mnt/d/analysis/grid", str(analysis_root or analysis_default()))]
    if os.environ.get("DIR_DATA"):
        mappings.append(("/mnt/f", os.environ["DIR_DATA"]))
    if os.environ.get("DIR_WORK"):
        mappings.append(("/mnt/d/data", str(Path(os.environ["DIR_WORK"]) / "data")))
    if os.environ.get("OS_FILES_ROOT"):
        mappings.append(("/mnt/d/files", os.environ["OS_FILES_ROOT"]))

    def walk(item):
        if isinstance(item, dict):
            return {key: walk(val) for key, val in item.items()}
        if isinstance(item, list):
            return [walk(val) for val in item]
        if isinstance(item, str):
            # Explicit/current analysis roots take priority over generic HPC
            # relocation; do not remap them a second time.
            for old, new in mappings[:2]:
                if item == old or item.startswith(old + "/"):
                    return new.rstrip("/") + item[len(old):]
            relocated = platform(item)
            if relocated != item:
                return relocated
            for old, new in mappings[2:]:
                if item == old or item.startswith(old + "/"):
                    return new.rstrip("/") + item[len(old):]
        return item
    return walk(value)


def read_config(path, analysis_root=None):
    source = Path(relocate(str(path), analysis_root))
    raw = json.loads(source.read_text())
    value = relocate(raw, analysis_root)
    # Only replace known repository defaults. Custom config values remain exact.
    aliases = {"score_root": ("GRID_SCORE_DIR", "/mnt/f/gen/ukb/pgs"),
               "ld_block_dir": ("GRID_LD_BLOCK_DIR", "/mnt/f/LD/block"),
               "dosage_supplement_dir": ("GRID_TARGET_DIR", "/mnt/f/gen/ukb/37/typ")}
    for key, (env, default) in aliases.items():
        if raw.get(key) == default and os.environ.get(env):
            value[key] = os.environ[env]
    return value
