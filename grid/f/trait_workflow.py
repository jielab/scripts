#!/usr/bin/env python3
"""Run a trait's native stages sequentially inside its one allocated job."""
from __future__ import annotations

import json
import os
import subprocess
import sys


def run_steps(steps, environment):
    if not isinstance(steps, list) or not steps:
        raise ValueError("A trait workflow requires at least one stage")
    for step in steps:
        if (not isinstance(step, dict) or not isinstance(step.get("key"), str)
                or not isinstance(step.get("argv"), list) or not step["argv"]
                or not all(isinstance(arg, str) for arg in step["argv"])
                or not isinstance(step.get("env", {}), dict)
                or not all(isinstance(k, str) and isinstance(v, str)
                           for k, v in step.get("env", {}).items())):
            raise ValueError("Invalid trait workflow stage")
    for index, step in enumerate(steps, 1):
        env = {**environment, **step.get("env", {})}
        # Children run native modules, without recursively dispatching this
        # bundle or submitting additional scheduler jobs.
        env.pop("GRID_TRAIT_STEPS", None)
        env["WORKFLOW_WORKER"] = "1"
        print(f"[trait] START {index}/{len(steps)} {step['key']}", flush=True)
        result = subprocess.run(step["argv"], env=env)
        if result.returncode:
            print(f"[trait] FAILED {step['key']}: exit {result.returncode}; later stages skipped",
                  file=sys.stderr, flush=True)
            return result.returncode if result.returncode > 0 else 128 - result.returncode
        print(f"[trait] DONE {index}/{len(steps)} {step['key']}", flush=True)
    return 0


def main():
    if os.environ.get("WORKFLOW_WORKER") != "1":
        raise ValueError("Use the public entry to schedule a trait workflow")
    if os.environ.get("submit") == "hpc" and (
            not os.environ.get("LSB_JOBID") or os.environ.get("WORKFLOW_RESOURCE_MANAGER") != "lsf"):
        raise ValueError("HPC trait workflows require an allocated LSF job")
    return run_steps(json.loads(os.environ["GRID_TRAIT_STEPS"]), os.environ)


if __name__ == "__main__":
    sys.exit(main())
