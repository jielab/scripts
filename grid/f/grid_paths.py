#!/usr/bin/env python3
"""Relocate repository defaults without rewriting explicit user paths.

The WSL paths in historical JSON files are defaults, not a requirement to mount
those drives on a cluster. The common workflow helper owns platform detection;
these overrides also make direct Python scientific workers portable.
"""
from __future__ import annotations

from settings import load_settings
import importlib.util
import json
import os
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
GRID = REPOSITORY / "grid"


# Use the shared adapter; this module only adds GRID's analysis-root override.
import sys as _paths_sys
_paths_sys.path.append(str(REPOSITORY / "0f"))
from workflow_paths import configured_path, relocate_defaults


def analysis_default():
    return os.environ.get("GRID_OUTPUT_ROOT", configured_path("ANALYSIS_ROOT", "grid"))


def relocate(value, analysis_root=None):
    if isinstance(value, dict):
        return {key: relocate(item, analysis_root) for key, item in value.items()}
    if isinstance(value, list):
        return [relocate(item, analysis_root) for item in value]
    if not isinstance(value, str):
        return value
    resolved = relocate_defaults(value)
    base = configured_path("ANALYSIS_ROOT", "grid")
    if resolved == base or resolved.startswith(base + "/"):
        return str(analysis_root or analysis_default()).rstrip("/") + resolved[len(base):]
    return resolved


def read_config(path, analysis_root=None):
    source = Path(relocate(str(path), analysis_root))
    raw = load_settings(source)
    value = relocate(raw, analysis_root)
    # Only replace known repository defaults. Custom config values remain exact.
    aliases = {"score_root": ("GRID_SCORE_DIR", configured_path("DIR_GEN", "ukb/pgs")),
               "ld_block_dir": ("GRID_LD_BLOCK_DIR", configured_path("DIR_LD", "block")),
               "dosage_supplement_dir": ("GRID_TARGET_DIR", configured_path("DIR_GEN", "ukb/37/typ"))}
    for key, (env, default) in aliases.items():
        if relocate_defaults(raw.get(key)) == default and os.environ.get(env):
            value[key] = os.environ[env]
    return value
