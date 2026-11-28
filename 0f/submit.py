#!/usr/bin/env python3
"""Shared project submission entry; project arguments are passed through unchanged."""
from __future__ import annotations

import argparse
from pathlib import Path
import importlib.util
import subprocess
from workflow_entry import execute_workflow
import sys

SHARED = Path(__file__).resolve().parent
PLANNERS = {
    "ems120": ("ems120", "workflow_definition.py"),
    "le8": ("le8", "workflow_definition.py"),
    "metag-qc": ("metag", "qc_workflow.py"),
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--project", required=True,
                        help="ems120, le8, metag-qc, or a project with --entry")
    # Stop parsing at the project name: --help and every later option belong
    # to that project's parser (including its own positional arguments).
    tokens = list(sys.argv[1:] if argv is None else argv)
    if tokens[:1] == ["--project"] and len(tokens) >= 2:
        project, arguments = tokens[1], tokens[2:]
    elif tokens and tokens[0].startswith("--project="):
        project, arguments = tokens[0].split("=", 1)[1], tokens[1:]
    else:
        parser.parse_args(tokens)
        parser.error("Put --project NAME first")
    if project in PLANNERS:
        folder, planner = PLANNERS[project]
        sys.path.insert(0, str(SHARED.parent / folder / "f"))
        target = SHARED.parent / folder / "f" / planner
    else:
        target = SHARED / "project_job.py"
        arguments = ["--project", project, *arguments]
    module_spec = importlib.util.spec_from_file_location('project_definition', target)
    definition = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(definition)
    if arguments[:1] == ['--verify-config'] and hasattr(definition, 'verify_config'):
        definition.verify_config(arguments[1:])
        return 0
    spec = definition.build_workflow(arguments)
    if spec is None: return 0  # Scientific summary worker; no scheduling.
    result = execute_workflow(spec)
    return int(result.get('exit_code', 0))


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        sys.exit(f"ERROR: {exc}")
