"""LE8 job plans with CPU/GPU reservations only, using the shared DAG dispatch and receipts."""
import argparse
import fcntl
import json
import math
import os
from contextlib import redirect_stdout
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import time

SHARED = Path(__file__).resolve().parents[2] / "0f"
sys.path.append(str(SHARED))
import workflow_jobs as common
from workflow_jobs import (KEY, atomic_json, check_binding, command_digest, digest, new_run_dir,
                           plan_digest, positive, resolve_mode, shared_flock, source_available,
                           source_identities, text, utcnow)
from file_lock import file_lock

VERSION = "le8-dag-cpu-gpu"
_STATE = common._state_path


def normalized(jobs, env=None):
    # Reuse shared argv, environment, DAG and CPU/GPU validation. Its default
    # memory and walltime belong to other workflows and are never serialized.
    items=common.normalized(jobs,env)
    for job in items:
        for key in ("memory_gib","memory_scope","hours"):
            job["resources"].pop(key,None)
    return items


def batch_number(root):
    root.mkdir(parents=True,exist_ok=True)
    counter=root/"batch_sequence.json"
    with file_lock(root/"batch_sequence.lock"):
        value=json.loads(counter.read_text())["last"]+1 if counter.is_file() else 1
        atomic_json(counter,{"last":value})
    return value


def share_allocation(first, second, key, mode):
    """Keep the two scientific commands and their per-scope environment overrides intact."""
    shared = {k:v for k,v in first["env"].items() if second["env"].get(k) == v}
    def command(job):
        overrides = [f"{k}={v}" for k,v in job["env"].items() if k not in shared]
        return ["env", *overrides, *job["argv"]]
    argv=[os.getenv("WORKFLOW_PYTHON","python3"), str(Path(__file__).resolve()), mode,
          *command(first),"--then",*command(second)]
    first.update(key=key,env=shared,argv=argv)


def join_endpoint_fits(jobs):
    """Pack terminal default CPU ICD and Any1 fits sharing the same trait, assay and prerequisites."""
    by_key = {job["key"]:job for job in jobs}
    parents = {p for job in jobs for p in job["depends"]}
    consumed=set()
    for sensitivity in jobs:
        key=sensitivity["key"]
        if not key.startswith("any1.") or not key.endswith(".analysis") or key in parents:
            continue
        tail=key.removeprefix("any1.")
        primary=by_key.get("icd."+tail)
        if not primary or primary["key"] in parents or primary["depends"] != sensitivity["depends"]:
            continue
        if primary["resources"]!=sensitivity["resources"] or primary["resources"]["gpus"] or primary["cwd"]!=sensitivity["cwd"]:
            continue
        share_allocation(primary,sensitivity,"icd-any1."+tail.removesuffix(".analysis"),"endpoint-fits")
        consumed.add(key)
    return [job for job in jobs if job["key"] not in consumed]


def join_aod_focus(jobs):
    """Share one LSF allocation for each AOD and its terminal timing/twins study."""
    by_key = {job["key"]: job for job in jobs}
    parents = {p for job in jobs for p in job["depends"]}
    renamed = {}
    consumed = set()
    for focus in jobs:
        key = focus["key"]
        if not key.startswith("focus.") or key in parents:
            continue
        trait = key.removeprefix("focus.")
        source = "aod." + trait
        aod = by_key.get(source)
        if not aod or focus["depends"] != [source] or aod["resources"] != focus["resources"] or aod["cwd"] != focus["cwd"]:
            continue
        joined = "aod-focus." + trait
        share_allocation(aod,focus,joined,"aod-focus")
        renamed[source] = joined
        consumed.add(key)
    result = [job for job in jobs if job["key"] not in consumed]
    for job in result:
        job["depends"] = [renamed.get(p, p) for p in job["depends"]]
    return result


def job_label(job,number,modules):
    key=job["key"]; parts=key.split(".")
    paper=all(m in ("paper_summary","paper","share","shiny") for m in modules)
    if key=="paper": parts=["paper","summary"]
    elif key=="dx.compare": parts=["dates","compare"]
    elif parts[0]=="aod": parts=[parts[1],"aod"]
    elif parts[0]=="aod-focus": parts=[parts[1],"aod","timing-twins"]
    elif parts[0]=="icd-any1": parts=[*parts[1:],"icd-any1"]
    elif parts[0]=="focus": parts=[parts[1],"c1","timing-twins"]
    elif parts[0]=="timing": parts=[*parts[1:],"any1","timing"]
    elif len(parts)==4:
        date,trait,layer,phase=parts
        if phase=="analysis":
            steps=json.loads(job["argv"][job["argv"].index("--steps-json")+1])
            labels=list(dict.fromkeys(s["key"].split(".")[2].split("_")[0] for s in steps))
            if paper: phase="-".join(modules)
            else: phase="-".join(labels)
        parts=(["paper"] if paper else [])+[trait,layer,date,phase]
    return ".".join(["le8",*parts,str(number)])


def dependency_condition(job, parent):
    """Reports consume available results after science ends, including failures.

    AOD remains a required input; scientific jobs keep success-only dependencies.
    """
    # For a merged job, the downstream worker checks the AOD step's own receipt.
    if parent.startswith("aod-focus."):
        return "ended"
    optional = job["key"] == "dx.compare" or job["key"].startswith("timing.")
    return "ended" if optional and not parent.startswith("aod.") else "done"


def check_parent(plan, data, job, parent, record):
    if record:
        check_binding(record, data, "Dependency status")
    if parent.startswith("aod-focus."):
        path = _STATE(plan, parent).with_suffix(".aod-ok")
        ready = json.loads(path.read_text()) if path.is_file() else {}
        if not ready or not record or ready.get("job_id") != record.get("job_id"):
            raise RuntimeError(f"AOD did not complete for this allocation: {parent}")
        check_binding(ready, data, "AOD completion")
        # LSF releases the downstream job after this allocation ends. Its only
        # scientific prerequisite was AOD, even if the subsequent focus fit failed.
        return
    if record.get("status") == "DONE":
        return
    if dependency_condition(job, parent) == "ended":
        terminal = record.get("status") in {"EXIT", "FAILED", "INTERRUPTED", "SKIPPED", "BLOCKED"}
        # An LSF kill may prevent the worker from writing its final receipt.
        # A missing receipt or a still-running parent is never treated as success.
        if not terminal and data["mode"] == "hpc":
            journal = plan.parent / "submission.json"
            saved = json.loads(journal.read_text()) if journal.is_file() else {}
            if saved:
                check_binding(saved, data, "Submission journal")
            job_id = saved.get("jobs", {}).get(parent, {}).get("job_id") or record.get("job_id")
            terminal = bool(job_id and common._probe([str(job_id)]).get(str(job_id)) in {"DONE", "EXIT"})
        if terminal:
            print(f"[LE8] {job['key']}: dependency {parent} ended without successful results; unavailable results will be skipped", flush=True)
            return
    raise RuntimeError(f"Dependency did not complete successfully: {parent}")


def write_plan(jobs,out_dir,*,mode,name,analysis_root,modules,announce=True):
    mode=resolve_mode(mode)
    shared={k:v for k,v in os.environ.items() if k.startswith(("DIR_","HPC_","_OS_DEFAULT_")) or k in ("ANALYSIS_ROOT","DATE_SCRATCH")}
    ordered=normalized(jobs,shared)
    sources=source_identities(ordered)
    from current_release import source_identity
    release=source_identity()
    if mode=="hpc":
        ordered=join_endpoint_fits(join_aod_focus(ordered))
    numbering=Path(analysis_root)
    if mode=="hpc" and os.getenv("ANALYSIS_ROOT"):
        numbering=Path(os.environ["ANALYSIS_ROOT"])/"le8"
    number=batch_number(numbering/"qc/scheduler")
    for job in ordered:
        r=job["resources"]
        job["env"].setdefault("OMP_NUM_THREADS",str(r["cores"]))
        job["job_name"]=job_label(job,number,modules)
        if not KEY.fullmatch(job["job_name"]): raise ValueError("Invalid LE8 job label: "+job["job_name"])
        if mode=="hpc" and not r["gpus"]:
            from qiming_resources import CPU_LIMITS, EXACT_CORE_QUEUES
            if r["queue"] in CPU_LIMITS:
                slots=CPU_LIMITS[r["queue"]][0]
                if r["cores"]>slots: raise ValueError("Requested cores exceed the selected single-node queue")
                if r["cores"] > 1 and r["queue"] not in EXACT_CORE_QUEUES:
                    r["cores"]=math.ceil(slots*math.ceil(4*r["cores"]/slots)/4)
    if len({j["job_name"] for j in ordered})!=len(ordered): raise ValueError("Duplicate LE8 job name")
    root=new_run_dir(Path(out_dir).resolve())
    for folder in ("commands","logs","status"): (root/folder).mkdir()
    path=root/"plan.json"
    helper=Path(__file__).resolve()
    sources[str(helper)]={"size":helper.stat().st_size}
    value=dict(schema=VERSION,created=utcnow(),name=name,mode=mode,configuration_sha256=digest(ordered),
               jobs=ordered,integrity_version=1,sources=sources,commands={},release=release)
    for job in ordered:
        key=job["key"]; r=job["resources"]
        lines=["#!/usr/bin/env bash",f"#BSUB -J {job['job_name']}",f"#BSUB -q {r['queue']}",
               f"#BSUB -n {r['cores']}",'#BSUB -R "span[hosts=1]"',
               "#BSUB -o "+shlex.quote(str(root/"logs"/(key+".%J.out"))),
               "#BSUB -e "+shlex.quote(str(root/"logs"/(key+".%J.err")))]
        if r["gpus"]: lines.append(f'#BSUB -gpu "num={r["gpus"]}/host:mode=exclusive_process:j_exclusive=yes"')
        lines.extend(["set -euo pipefail","umask 077"])
        if mode=="hpc":
            lines.extend(['[[ -n "${LSB_JOBID:-}" ]] || exit 2',
                          f'[[ ${{LSB_DJOB_NUMPROC:-0}} -ge {r["cores"]} ]] || exit 2',
                          "source "+shlex.quote(str(SHARED/"0os.sh"))])
        lines.extend(["# Dependencies: "+(", ".join(job["depends"]) or "(none)"),
                      "# Scientific command: "+shlex.join(job["argv"]),
                      shlex.join(["exec",os.getenv("WORKFLOW_PYTHON","python3"),str(helper),"worker","--plan",str(path),"--key",key])])
        cmd=root/"commands"/(key+".cmd")
        cmd.write_text("\n".join(lines)+"\n");cmd.chmod(0o750)
        value["commands"][key]=command_digest(cmd)
    value["plan_sha256"]=plan_digest(value)
    atomic_json(path,value)
    if announce:
        print(f"[PLAN] {path} ({len(ordered)} jobs; {mode})",flush=True)
    return path


def load_plan(path):
    path = Path(path).resolve()
    d = json.loads(path.read_text())
    if d.get("schema") != VERSION:
        raise ValueError("Unsupported workflow plan")
    from current_release import source_identity
    if d.get("integrity_version") != 1 or d.get("release") != source_identity():
        raise ValueError("Plan is not from the current LE8 release; create a new plan")
    d["jobs"] = normalized(d["jobs"])
    if d.get("configuration_sha256") != digest(d["jobs"]):
        raise ValueError("Workflow configuration changed; create a new plan instead of editing a frozen plan")
    if d.get("mode") not in {"local", "hpc"}:
        raise ValueError("Invalid frozen workflow mode")
    if d.get("integrity_version") == 1:
        if d.get("plan_sha256") != plan_digest(d):
            raise ValueError("Workflow plan identity changed")
        for source, expected in d["sources"].items():
            if not source_available(source):
                raise ValueError("Workflow source is missing: " + source)
        for key, expected in d["commands"].items():
            command = path.parent / "commands" / (key + ".cmd")
            if not command.is_file() or command_digest(command) != expected:
                raise ValueError("Workflow command changed; generate a new plan: " + str(command))
    elif d.get("integrity_version") is not None:
        raise ValueError("Unsupported workflow integrity version")
    return path, d


def worker(plan, key):
    plan, d = load_plan(plan)
    job = next((j for j in d["jobs"] if j["key"] == key), None)
    if job is None:
        raise ValueError("Unknown job key")
    if d["mode"] == "hpc" and not os.environ.get("LSB_JOBID"):
        raise RuntimeError("HPC worker requires an LSF allocation")
    state = _STATE(plan, key)
    with state.with_suffix(".lock").open("a+") as lock:
        try:
            shared_flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("This plan/job is already running") from exc
        for parent in job["depends"]:
            # Existing LE8 plans may still list the optional deep dependency.
            if (job.get("job_name", "").startswith("le8.")
                    and job["key"].endswith(".evaluation")
                    and parent == job["key"].removesuffix(".evaluation") + ".abm_deep"):
                continue
            parent_state = _STATE(plan, parent)
            parent_record = json.loads(parent_state.read_text()) if parent_state.is_file() else {}
            check_parent(plan, d, job, parent, parent_record)
        env = {**os.environ, **job["env"], "WORKFLOW_WORKER": "1", "SCRIPT_CONSOLE_ACTIVE": "1",
               "submit": d["mode"], "SUBMIT": d["mode"]}
        if key.startswith("aod-focus."):
            marker = state.with_suffix(".aod-ok")
            marker.unlink(missing_ok=True)
            env.update(LE8_AOD_STEP_RECEIPT=str(marker), LE8_AOD_PLAN_SHA=d["plan_sha256"])
        r = job["resources"]
        env.update(WORKFLOW_CORES=str(min(int(job["env"].get("WORKFLOW_CORES", r["cores"])), r["cores"])),
                   WORKFLOW_GPUS=str(r["gpus"]), WORKFLOW_RESOURCE_MANAGER="lsf" if d["mode"] == "hpc" else "local")
        env.update(OMP_NUM_THREADS=job["env"].get("OMP_NUM_THREADS", str(r["cores"])), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                   NUMEXPR_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1")
        # Processor selection is resolved before allocation, for both machines.
        # Zero requested GPUs means CPU, including local runs on a GPU machine.
        if not r["gpus"]:
            env["CUDA_VISIBLE_DEVICES"] = ""
            env["WORKFLOW_DEVICE"] = "cpu"
        temporary = (Path(env["WORKFLOW_TMP_ROOT"]) / plan.parent.name / key
                     if env.get("WORKFLOW_TMP_ROOT") else plan.parent / "tmp" / key)
        temporary.mkdir(parents=True, exist_ok=True)
        env.update(TMPDIR=str(temporary), TMP=str(temporary), TEMP=str(temporary))
        record = dict(key=key, status="RUN", started=utcnow(), job_id=os.getenv("LSB_JOBID"), pid=os.getpid(),
                      plan_sha256=d.get("plan_sha256"), completion_contract="process_exit_only")
        atomic_json(state, record)
        process = None
        previous = {}
        def stop(signum, _frame):
            if process is not None and process.poll() is None:
                os.killpg(process.pid, signum)
            raise KeyboardInterrupt
        try:
            for sig in (signal.SIGTERM, signal.SIGINT):
                previous[sig] = signal.signal(sig, stop)
            print("[EXEC] " + shlex.join(job["argv"]), flush=True)
            process = subprocess.Popen(job["argv"], cwd=job["cwd"], env=env, start_new_session=True)
            code = process.wait()
        except BaseException:
            if process is not None and process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            record.update(status="EXIT", ended=utcnow(), returncode=130)
            atomic_json(state, record)
            raise
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)
        record.update(status="DONE" if code == 0 else "EXIT", ended=utcnow(), returncode=code)
        atomic_json(state, record)
        return code


def _wait_hpc(path, plan, receipts, poll=5.0):
    remaining = dict(receipts)
    while remaining:
        local = common._read_states(path,plan)
        ids = [v["job_id"] for v in remaining.values()]
        observed = common._probe(ids)
        for key, row in list(remaining.items()):
            stat = local.get(key, {}).get("status")
            cluster = observed.get(row["job_id"])
            if stat == "EXIT" or cluster == "EXIT":
                raise RuntimeError(f"LSF job failed: {key} ({row['job_id']}); success-dependent jobs remain blocked")
            if stat == "DONE" and cluster in {"DONE", None}:
                remaining.pop(key)
            elif cluster is None:
                raise RuntimeError(f"Cannot verify LSF job {row['job_id']}; disappearance is not success. See submission.json")
        if remaining:
            time.sleep(min(max(poll, .05), 30.0))
    return {"status": "DONE", "jobs": list(receipts), "plan": str(path)}


class SubmissionOutput:
    def __init__(self, stream, jobs):
        self.stream = stream
        self.names = {job["key"]:job["job_name"] for job in jobs}

    def write(self, value):
        parts = value.split()
        if len(parts)==3 and parts[0]=="[SUBMITTED]" and parts[1] in self.names:
            value = f"[SUBMITTED] {self.names[parts[1]]} {parts[2]}"
        return self.stream.write(value)

    def flush(self):
        self.stream.flush()


def print_result(result):
    if result.get("status") != "SUBMITTED":
        print(json.dumps(result,indent=2))


def run_plan(plan_path, *, submit=None, jobs=128, wait=False, dry_run=False, resume=True):
    path, d = load_plan(plan_path)
    # A saved plan carries the resolved mode; never reinterpret it on resume.
    submit = d["mode"] if submit is None else resolve_mode(submit)
    jobs = positive(jobs, "jobs", True)
    if submit != d["mode"]:
        raise ValueError("Invalid concurrency or plan/submit mode mismatch")
    if dry_run:
        return {"status": "PLANNED", "plan": str(path), "jobs": len(d["jobs"])}
    with (path.parent / "dispatch.lock").open("a+") as lock:
        try:
            shared_flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another dispatcher is using this plan") from exc
        if submit == "local":
            return common._run_local(path, d, jobs, resume)
        inherited = os.environ.pop("WORKFLOW_JOB_NAME",None)
        try:
            with redirect_stdout(SubmissionOutput(sys.stdout,d["jobs"])):
                result=common._run_hpc(path,d,jobs,False,resume,dependency_condition=dependency_condition)
        finally:
            if inherited is not None:
                os.environ["WORKFLOW_JOB_NAME"] = inherited
        if wait: return _wait_hpc(path,d,result["jobs"])
        return result


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if argv[:1]==["endpoint-fits"]:
        split=argv.index("--then")
        primary=subprocess.run(argv[1:split]).returncode
        if primary<0 or primary in (130,143):
            return 128-primary if primary<0 else primary
        print(f"[LE8] ICD exit={primary}. Starting Any1 in the same allocation.",flush=True)
        sensitivity=subprocess.run(argv[split+1:]).returncode
        if sensitivity<0:
            return 128-sensitivity
        return primary or sensitivity
    if argv[:1]==["aod-focus"]:
        split = argv.index("--then")
        result = subprocess.run(argv[1:split])
        if result.returncode:
            return result.returncode if result.returncode > 0 else 128 - result.returncode
        atomic_json(Path(os.environ["LE8_AOD_STEP_RECEIPT"]),
                    dict(plan_sha256=os.environ["LE8_AOD_PLAN_SHA"], job_id=os.environ["LSB_JOBID"]))
        print("[LE8] AOD complete. Starting C1 timing/twins.", flush=True)
        os.execvpe(argv[split+1], argv[split+1:], os.environ)
    if argv[:1]==["aod"]:
        out,traits,queue,slots,dry,entry=argv[1:7]
        if argv[7]!="--": raise ValueError("Missing standalone AOD arguments")
        root=Path(out).resolve()
        key="aod."+(traits if "," not in traits else "multiple")
        job=common.Job(key,["bash",entry,"--execute",*argv[8:]],[],
                       dict(cores=int(slots),queue=queue,gpus=0),{},str(Path(entry).resolve().parents[1]))
        mode=resolve_mode()
        path=write_plan([job],root/"jobs",mode=mode,name="le8",analysis_root=root,modules=["aod"],announce=dry=="1" or mode!="hpc")
        print_result(run_plan(path,submit=mode,dry_run=dry=="1",wait=os.getenv("SUBMIT_WAIT","").upper()=="TRUE"))
        return 0
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    w = sub.add_parser("worker")
    w.add_argument("--plan", required=True)
    w.add_argument("--key", required=True)
    r = sub.add_parser("run")
    r.add_argument("plan")
    r.add_argument("--submit", choices=["auto", "local", "hpc"], default=None)
    r.add_argument("--jobs-N", type=int, default=128)
    r.add_argument("--wait", action="store_true")
    r.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)
    try:
        if a.command == "worker":
            return worker(a.plan, a.key)
        result = run_plan(a.plan, submit=a.submit, jobs=a.jobs_N, wait=a.wait, dry_run=a.dry_run)
        print_result(result)
        return 0
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        p.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
