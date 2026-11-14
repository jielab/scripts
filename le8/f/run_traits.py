#!/usr/bin/env python3
"""Run the sequential steps of one LE8 allocation; all dependencies are in its DAG."""
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
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT.parent/"0f"))
from file_lock import flock as shared_flock
from workflow_jobs import atomic_json


@contextmanager
def aggregate_lock(root, enabled):
    if not enabled:
        yield
        return
    path=Path(root)/"logs/locks/scope-aggregate.lock"
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("a") as lock:
        shared_flock(lock,fcntl.LOCK_EX)
        yield


def claim_aggregates(root,scope,batch_id,scopes):
    directory=Path(root)/"qc/batches"/batch_id
    atomic_json(directory/(scope+".ready.json"),dict(scope=scope))
    with aggregate_lock(root,True):
        if not all((directory/(s+".ready.json")).is_file() for s in scopes):
            return False
        owner=directory/"aggregate_owner.json"
        if owner.is_file():return json.loads(owner.read_text())["scope"]==scope
        atomic_json(owner,dict(scope=scope))
        return True


def run_steps(root,scope,steps,*,batch_id=None,scopes=None,phase=None):
    if os.getenv("submit")=="hpc" and not os.getenv("LSB_JOBID"):
        raise RuntimeError("Scope worker requires an LSF allocation")
    env=dict(os.environ,WORKFLOW_WORKER="1",LE8_SCOPE_WORKER="1")
    trait,layer=scope.split(".")
    path=Path(root)/trait/layer/"logs"/f"{phase or 'analysis'}.status.json"
    state=dict(scope=scope,phase=phase,status="RUNNING",job_id=os.getenv("LSB_JOBID"),
               started=datetime.now(timezone.utc).isoformat(),steps=[])
    atomic_json(path,state)
    skip_file=path.with_name(path.stem+".skip.json")
    skip_file.unlink(missing_ok=True)
    skip_reason=None
    # The preparation allocation may skip a whole scope. Carry that decision
    # through this plan's CPU/GPU allocations, never an earlier run's QC file.
    batch_skip=next((Path(s["env"]["LE8_C1_BATCH_DIR"])/"scope.skip.json"
                     for s in steps if s.get("env",{}).get("LE8_C1_BATCH_DIR")),None)
    if batch_skip is not None and batch_skip.is_file():
        skipped=json.loads(batch_skip.read_text())
        if skipped.get("scope")!=scope or skipped.get("status")!="SKIPPED":
            raise ValueError("Invalid batch scope skip record")
        skip_reason=skipped["reason"]
    aggregate_owner=None
    for index,step in enumerate(steps):
        if skip_reason and not step.get("aggregate"):
            state["steps"].append(dict(key=step["key"],status="SKIPPED",exit_code=0,reason=skip_reason))
            atomic_json(path,state)
            continue
        if step.get("aggregate") and batch_id:
            if aggregate_owner is None:
                aggregate_owner=claim_aggregates(root,scope,batch_id,scopes or [scope])
            if not aggregate_owner:
                state["steps"].append(dict(key=step["key"],status="DEFERRED",exit_code=0,
                    reason="The last completed scope runs the common summary"))
                atomic_json(path,state)
                continue
        print(f"[LE8] START {step['key']}: {shlex.join(step['argv'])}",flush=True)
        try:
            with aggregate_lock(root,step.get("aggregate",False)):
                code=subprocess.run(step["argv"],env={**env,**step.get("env",{}),"LE8_SCOPE_SKIP_FILE":str(skip_file)}).returncode
            if code<0 or code in (130,143):
                state["status"]="INTERRUPTED";atomic_json(path,state)
                return 128-code if code<0 else code
            record=dict(key=step["key"],status="DONE" if code==0 else "FAILED",exit_code=code)
            if code==0 and not step.get("aggregate") and skip_file.is_file():
                skipped=json.loads(skip_file.read_text())
                if skipped.get("scope")!=scope or skipped.get("status")!="SKIPPED":
                    raise ValueError("Invalid scope skip record")
                skip_reason=skipped["reason"]
                if batch_skip is not None:
                    atomic_json(batch_skip,skipped)
                record.update(status="SKIPPED",reason=skip_reason)
        except OSError as error:
            record=dict(key=step["key"],status="FAILED",exit_code=1,error=str(error))
        state["steps"].append(record);atomic_json(path,state)
        print(f"[LE8] {record['status']} {step['key']}",flush=True)
        if record["status"]=="FAILED":
            state["steps"].extend(dict(key=s["key"],status="BLOCKED",blocked_by=step["key"])
                                  for s in steps[index+1:])
            break
    failed=[r["key"] for r in state["steps"] if r["status"]=="FAILED"]
    state.update(status="FAILED" if failed else "SKIPPED" if skip_reason else "DONE",failed_steps=failed,
                 finished=datetime.now(timezone.utc).isoformat())
    atomic_json(path,state)
    return int(bool(failed))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__,allow_abbrev=False)
    p.add_argument("--analysis-root",type=Path,required=True)
    p.add_argument("--scope",required=True)
    p.add_argument("--steps-json",required=True)
    p.add_argument("--batch-id")
    p.add_argument("--scopes")
    p.add_argument("--phase")
    a=p.parse_args(argv)
    scopes=a.scopes.split(",") if a.scopes else [a.scope]
    if any(not re.fullmatch(r"[A-Za-z0-9_]+\.(prot|met)",s) for s in [a.scope,*scopes]):
        p.error("Expected trait.prot or trait.met scopes")
    if a.batch_id and not re.fullmatch(r"(?:run-[0-9TZp-]+|[a-f0-9]{32})",a.batch_id):
        p.error("Invalid batch ID")
    if a.phase and not re.fullmatch(r"[A-Za-z0-9_]+",a.phase):p.error("Invalid phase")
    return run_steps(a.analysis_root,a.scope,json.loads(a.steps_json),batch_id=a.batch_id,
                     scopes=scopes,phase=a.phase)


if __name__=="__main__":
    raise SystemExit(main())
