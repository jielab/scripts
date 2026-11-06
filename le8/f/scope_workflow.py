#!/usr/bin/env python3
"""Run a resource-compatible group of trait/layer stages in one allocation.

Expected scientific unavailability is handled by the individual analyses.
Split allocations stop on failure so success-dependent jobs cannot consume
partial results. Historical single-allocation plans retain their continuation policy.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


@contextmanager
def aggregate_lock(root, enabled):
    if not enabled:
        yield
        return
    path = root / "logs/locks/scope-aggregate.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def save_status(path, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    os.replace(handle.name, path)


def claim_aggregates(root, scope, batch_id, scopes):
    """Let the last scientific scope own aggregation, without a barrier job."""
    directory = root / "logs/scopes/batches" / batch_id
    save_status(directory / (scope + ".ready.json"), dict(scope=scope))
    with aggregate_lock(root, True):
        if not all((directory / (s + ".ready.json")).is_file() for s in scopes):
            return False
        owner = directory / "aggregate_owner.json"
        if owner.exists():
            return json.loads(owner.read_text())["scope"] == scope
        save_status(owner, dict(scope=scope))
        return True


def run_steps(root, scope, steps, *, batch_id=None, scopes=None, phase=None, fail_fast=False):
    if os.getenv("submit") == "hpc" and not os.getenv("LSB_JOBID"):
        raise RuntimeError("Scope worker requires an LSF allocation")
    # Every scientific child inherits WORKFLOW_WORKER, preventing nested bsub.
    env = dict(os.environ, WORKFLOW_WORKER="1", LE8_SCOPE_WORKER="1")
    status_key = f"{scope}.{phase}" if phase else scope
    status = Path(root) / "logs/scopes" / f"{status_key}.json"
    result = dict(scope=scope, job_id=os.getenv("LSB_JOBID"),
                  phase=phase, started=datetime.now(timezone.utc).isoformat(), status="RUN", steps=[])
    save_status(status, result)
    aggregate_owner = None
    for step in steps:
        if step.get("aggregate") and batch_id:
            if aggregate_owner is None:
                aggregate_owner = claim_aggregates(Path(root), scope, batch_id, scopes or [scope])
            if not aggregate_owner:
                result["steps"].append(dict(key=step["key"], status="DEFERRED", exit_code=0,
                    reason="Aggregates run inside the last finished scope job"))
                save_status(status, result)
                continue
        print(f"[LE8] START {step['key']}: {shlex.join(step['argv'])}", flush=True)
        try:
            with aggregate_lock(Path(root), step.get("aggregate", False)):
                completed = subprocess.run(step["argv"], env={**env, **step.get("env", {})})
            code = completed.returncode
            if code < 0 or code in (130, 143):
                # Cancellation must stop the job, not start another analysis.
                result["status"] = "INTERRUPTED"
                save_status(status, result)
                return 128 - code if code < 0 else code
            record = dict(key=step["key"], exit_code=code, status="DONE" if code == 0 else "FAILED")
        except OSError as error:
            record = dict(key=step["key"], exit_code=1, status="FAILED", error=str(error))
        result["steps"].append(record)
        save_status(status, result)
        if fail_fast and record["status"] == "FAILED":
            print(f"[LE8] FAILED {step['key']}; stopping this phase and blocking dependent jobs", flush=True)
            break
        print(f"[LE8] {record['status']} {step['key']} exit={record['exit_code']}; "
              "continuing with remaining analyses", flush=True)
    failed = [s["key"] for s in result["steps"] if s["status"] == "FAILED"]
    result.update(status="FAILED" if failed else "DONE", failed_steps=failed,
                  finished=datetime.now(timezone.utc).isoformat())
    save_status(status, result)
    print(f"[LE8] {result['status']} {scope}; stage status: {status}", flush=True)
    return int(bool(failed))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis-root", required=True, type=Path)
    p.add_argument("--scope", required=True)
    p.add_argument("--steps-json", required=True)
    p.add_argument("--batch-id")
    p.add_argument("--scopes")
    p.add_argument("--phase")
    p.add_argument("--fail-fast", action="store_true")
    a = p.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9_]+\.(prot|met)", a.scope):
        p.error("Expected one trait.prot or trait.met scope")
    scopes = a.scopes.split(",") if a.scopes else [a.scope]
    if any(not re.fullmatch(r"[A-Za-z0-9_]+\.(prot|met)", s) for s in scopes):
        p.error("Invalid scope list")
    if a.batch_id and not re.fullmatch(r"[a-f0-9]{32}", a.batch_id):
        p.error("Invalid batch ID")
    if a.phase and not re.fullmatch(r"[A-Za-z0-9_]+", a.phase):
        p.error("Invalid phase label")
    return run_steps(a.analysis_root, a.scope, json.loads(a.steps_json), batch_id=a.batch_id,
                     scopes=scopes, phase=a.phase, fail_fast=a.fail_fast)


if __name__ == "__main__":
    raise SystemExit(main())
