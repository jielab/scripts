"""Use the shell path table's environment detection from Python planners."""
from __future__ import annotations
import os
from pathlib import Path
import subprocess


def resolve_mode(mode=None):
    """Resolve auto without importing scientific runtimes or inspecting data."""
    requested = mode or os.environ.get("submit") or os.environ.get("SUBMIT") or "auto"
    if requested in {"local", "hpc"}:
        return requested
    if requested != "auto":
        raise ValueError("submit must be auto, local or hpc")
    env = dict(os.environ, submit="auto", SUBMIT="auto")
    value = subprocess.check_output(
        ["bash", "-c", 'source "$1" || exit; printf "%s" "$submit"',
         "workflow-mode", str(Path(__file__).with_name("0os.sh"))],
        env=env, text=True,
    ).strip()
    if value not in {"local", "hpc"}:
        raise ValueError("0os.sh returned an invalid execution mode")
    return value
