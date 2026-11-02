"""Relocate repository defaults to the 0os.sh path table.

Use this only while loading built-in defaults, before applying explicit user
CLI values. An explicit path is never guessed or searched for elsewhere.
"""
from __future__ import annotations

import os
from pathlib import Path


def relocate_defaults(value):
    if isinstance(value, dict):
        return {k: relocate_defaults(v) for k, v in value.items()}
    if isinstance(value, list):
        return [relocate_defaults(v) for v in value]
    if isinstance(value, tuple):
        return tuple(relocate_defaults(v) for v in value)
    if not isinstance(value, str):
        return value
    root = str(Path(__file__).resolve().parents[1])
    mappings = [("/mnt/d/scripts", os.environ.get("SCRIPTS_ROOT", root))]
    if os.getenv("submit", os.getenv("SUBMIT", "local")) == "hpc":
        work = os.getenv("DIR_WORK", "/work/sph-huangj")
        data = os.getenv("DIR_DATA", "/data/sph-huangj")
        scratch = os.getenv("DIR_SCRATCH", f"/scratch/{os.getenv('DATE_SCRATCH', '2026-10-31')}/sph-huangj")
        mappings += [("/mnt/d/analysis", os.getenv("ANALYSIS_ROOT", scratch + "/analysis")),
                     ("/mnt/f/software", os.getenv("DIR_SOFTWARE", work + "/software")),
                     ("/mnt/d/files", os.getenv("OS_FILES_ROOT", data + "/files")),
                     ("/mnt/f", data), ("/mnt/d", work)]
    for old, new in mappings:
        if value == old or value.startswith(old + "/"):
            return new.rstrip("/") + value[len(old):]
    return value
