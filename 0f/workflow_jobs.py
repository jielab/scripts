#!/usr/bin/env python3
"""Standard-library DAG executor for reviewable local and LSF .cmd jobs.

Plans contain argv arrays, never shell snippets. Every new plan has a fresh run
directory; resume operates on that exact plan and never guesses that a new run
has unchanged scientific inputs. Scientific stages retain their own validated
checkpoints. Submission receipts are flushed after every accepted LSF job.
"""
from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait as futures_wait
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import fcntl
from workflow_lock import flock as shared_flock
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import uuid
from workflow_mode import resolve_mode
from workflow_layout import new_run_dir

from submission_io import atomic_json, utcnow
from lsf_submit import submit_graph

VERSION = "workflow-dag-1"
KEY = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$")


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def command_digest(path):
    # Presentation comments may expand a long scientific argv for review.
    # Executable lines and every LSF directive remain immutable.
    lines = [line for line in Path(path).read_text().splitlines()
             if line.strip() and (not line.lstrip().startswith("#") or line.lstrip().startswith("#BSUB"))]
    return digest(lines)


def plan_digest(value):
    return digest({key: value[key] for key in
        ("schema", "name", "mode", "configuration_sha256", "jobs", "sources", "commands")})


def source_identities(jobs):
    shared = Path(__file__).resolve().parent
    paths = {shared / name for name in
             ("workflow_jobs.py", "workflow_lock.py", "workflow_mode.py", "workflow_layout.py",
              "0os.sh", "qiming_resources.py", "project_entry.sh", "project_job.py",
              "submit.py", "workflow_entry.py", "lsf_submit.py", "submission_io.py", "command_jobs.py", "submit_jobs.py")}
    for job in jobs:
        definition = Path(job['cwd']) / 'f' / 'workflow_definition.py'
        if definition.is_file(): paths.add(definition.resolve())
        # Bind direct entry/worker code, without opening scientific datasets.
        # Imported scientific modules and data use the owning project's receipts.
        for value in job["argv"][:2]:
            path = Path(value)
            if path.suffix not in {".sh", ".py", ".R"}:
                continue
            if not path.is_absolute():
                path = Path(job["cwd"]) / path
            if path.is_file():
                paths.add(path.resolve())
    return {str(path): {"size": path.stat().st_size} for path in sorted(paths)}


def check_binding(record, plan, label):
    if plan.get("integrity_version") == 1 and record.get("plan_sha256") != plan["plan_sha256"]:
        raise ValueError(f"{label} belongs to another or unverified plan; preserve it and create a new plan")


@dataclass
class Job:
    key: str
    argv: list[str]
    depends: list[str] = field(default_factory=list)
    resources: dict = field(default_factory=dict)
    env: dict = field(default_factory=dict)
    cwd: str | None = None
    job_name: str | None = None




def text(value, name):
    value = str(value)
    if any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError(f"Control character in {name}")
    return value


def positive(value, name, integral=False, allow_zero=False):
    if isinstance(value, bool):
        raise ValueError(f"Invalid {name}")
    v = float(value)
    if not math.isfinite(v) or v < 0 or (not allow_zero and v == 0):
        raise ValueError(f"Invalid {name}")
    if integral and v != int(v):
        raise ValueError(f"{name} must be an integer")
    return int(v) if integral else v


def resources(values):
    unknown = set(values) - {"queue", "cores", "memory_gib", "hours", "gpus", "gpu_model", "memory_scope"}
    if unknown:
        raise ValueError(f"Unknown resource keys: {sorted(unknown)}")
    r = dict(queue=os.environ.get("HPC_CPU_QUEUE", "33"), cores=8, memory_gib=32, hours=24, gpus=0, memory_scope="host")
    r.update(values)
    r["queue"] = text(r["queue"], "queue")
    if not KEY.fullmatch(r["queue"]):
        raise ValueError("Invalid LSF queue name")
    r["cores"] = positive(r["cores"], "cores", True)
    r["memory_gib"] = positive(r["memory_gib"], "memory_gib")
    r["hours"] = positive(r["hours"], "hours")
    r["gpus"] = positive(r["gpus"], "gpus", True, True)
    if r["memory_scope"] not in {"host", "site_host", "task"}:
        raise ValueError("memory_scope must be host, site_host or task")
    if r.get("gpu_model") and not KEY.fullmatch(str(r["gpu_model"])):
        raise ValueError("Invalid GPU model")
    return r


def normalized(jobs, common_env=None):
    env0 = dict(common_env or {})
    result = []
    for item in jobs:
        j = asdict(item) if isinstance(item, Job) else dict(item)
        if not KEY.fullmatch(str(j.get("key", ""))):
            raise ValueError("Job keys must be safe, unique labels")
        if j.get("job_name") is not None and not KEY.fullmatch(str(j["job_name"])):
            raise ValueError("Job names must be safe labels")
        argv = j.get("argv")
        if not isinstance(argv, (list, tuple)) or not argv:
            raise ValueError("A job needs a nonempty argv list")
        j["argv"] = [text(a, "argv") for a in argv]
        j["depends"] = list(dict.fromkeys(j.get("depends", [])))
        j["resources"] = resources(j.get("resources", {}))
        env = {**env0, **j.get("env", {})}
        for name in env:
            if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", str(name)):
                raise ValueError(f"Invalid environment name: {name}")
        j["env"] = {name: text(value, "environment value") for name, value in env.items()}
        j["cwd"] = str(Path(j.get("cwd") or os.getcwd()).expanduser().resolve())
        text(j["cwd"], "cwd")
        result.append(j)
    by_key = {j["key"]: j for j in result}
    if len(by_key) != len(result) or not result:
        raise ValueError("Plan requires at least one job and unique keys")
    for j in result:
        if set(j["depends"]) - set(by_key):
            raise ValueError(f"Missing dependencies for {j['key']}")
    ordered, pending, seen = [], list(result), set()
    while pending:
        ready = [j for j in pending if set(j["depends"]) <= seen]
        if not ready:
            raise ValueError("Cyclic workflow dependencies")
        for j in ready:
            pending.remove(j)
            ordered.append(j)
            seen.add(j["key"])
    return ordered


def write_plan(jobs, out_dir, *, mode=None, name="workflow", common_env=None):
    mode = resolve_mode(mode)
    if mode not in {"local", "hpc"} or not KEY.fullmatch(name):
        raise ValueError("Invalid mode or plan name")
    shared = {k:v for k,v in os.environ.items() if k.startswith(('DIR_', 'HPC_', '_OS_DEFAULT_')) or k in ('ANALYSIS_ROOT', 'DATE_SCRATCH')}
    ordered = normalized(jobs, {**shared, **(common_env or {})})
    for job in ordered:
        # Reserve LSF slots independently of scientific thread budgets.
        job["env"].setdefault("OMP_NUM_THREADS", str(job["resources"]["cores"]))
    if mode == "hpc":
        from qiming_resources import cpu_slots
        for job in ordered:
            r = job["resources"]
            if not r["gpus"]:
                r["queue"], r["cores"] = cpu_slots(r["queue"], r["cores"], r["memory_gib"], r["hours"])
    root = Path(out_dir).expanduser().resolve()
    text(root, "plan directory")
    fingerprint = digest(ordered)
    # Short immutable run directories prevent an active plan or its receipts
    # being overwritten. Results themselves retain their stable project paths.
    root = new_run_dir(root)
    for folder in ("commands", "logs", "status"):
        (root / folder).mkdir()
    path = root / "plan.json"
    value = dict(schema=VERSION, created=utcnow(), name=name, mode=mode,
                 configuration_sha256=fingerprint, jobs=ordered, integrity_version=1,
                 sources=source_identities(ordered), commands={})
    helper = str(Path(__file__).resolve())
    worker_python = os.environ.get("WORKFLOW_PYTHON", "python3")
    for j in ordered:
        r = j["resources"]
        if r["memory_scope"] == "task":
            # Round per-slot MiB, not GiB: 160 GiB / 39 slots must not
            # become 5 GiB/slot (195 GiB) on a 192 GiB node.
            mem = math.ceil(r["memory_gib"] * 1024 / r["cores"])
            reservation = f"{mem}MB/task"
        else:
            mem = math.ceil(r["memory_gib"])
            suffix = "/host" if r["memory_scope"] == "host" else ""
            reservation = f"{mem}GB{suffix}"
        job_name = os.getenv("WORKFLOW_JOB_NAME") or j.get("job_name") or f"{name}.{j['key']}"
        lines = ["#!/usr/bin/env bash", f"#BSUB -J {job_name}", f"#BSUB -q {r['queue']}",
                 f"#BSUB -n {r['cores']}", f'#BSUB -R "span[hosts=1] rusage[mem={reservation}]"',
                 f"#BSUB -o {shlex.quote(str(root / 'logs' / (j['key'] + '.%J.out')))}",
                 f"#BSUB -e {shlex.quote(str(root / 'logs' / (j['key'] + '.%J.err')))}"]
        if r["gpus"]:
            gpu = f"num={r['gpus']}/host:mode=exclusive_process:j_exclusive=yes"
            if r.get("gpu_model"):
                gpu += f":gmodel={r['gpu_model']}"
            lines.append(f'#BSUB -gpu "{gpu}"')
        lines.extend(["", "set -euo pipefail", "umask 077"])
        if mode == "hpc":
            lines.extend(['[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }',
                          f'[[ ${{LSB_DJOB_NUMPROC:-0}} -ge {r["cores"]} ]] || {{ echo "Insufficient allocated CPU slots" >&2; exit 2; }}'])
        if mode == "hpc":
            lines.append("source " + shlex.quote(str(Path(__file__).resolve().with_name("0os.sh"))))
        lines.append(f"# Dependencies: {', '.join(j['depends']) or '(none)'}")
        lines.append("# Scientific command: " + shlex.join(j["argv"]))
        lines.append(shlex.join(["exec", worker_python, helper, "worker", "--plan", str(path), "--key", j["key"]]))
        command = root / "commands" / (j["key"] + ".cmd")
        command.write_text("\n".join(lines) + "\n", encoding="utf-8")
        command.chmod(0o750)
        value["commands"][j["key"]] = command_digest(command)
    value["plan_sha256"] = plan_digest(value)
    atomic_json(path, value)
    print(f"[PLAN] {path} ({len(ordered)} jobs; {mode})", flush=True)
    return path




def source_available(source):
    path = Path(source)
    if path.is_file(): return True
    # Old frozen plans name the former shared project modules. Their worker
    # commands are unchanged; resolve only these explicit source relocations.
    moved = {'ems120_jobs.py': '../ems120/f/workflow_definition.py',
             'le8_jobs.py': '../le8/f/workflow_definition.py',
             'flux_jobs.py': '../flux/f/workflow_definition.py',
             'metag_qc_jobs.py': '../metag/f/qc_workflow.py',
             'grid_jobs.py': 'submit_jobs.py'}
    return path.parent.name == '0f' and path.name in moved and (path.parent / moved[path.name]).is_file()


def load_plan(path):
    path = Path(path).resolve()
    d = json.loads(path.read_text())
    if d.get("schema") != VERSION:
        raise ValueError("Unsupported workflow plan")
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


def _state_path(plan, key):
    return plan.parent / "status" / f"{key}.json"


def worker(plan, key):
    plan, d = load_plan(plan)
    job = next((j for j in d["jobs"] if j["key"] == key), None)
    if job is None:
        raise ValueError("Unknown job key")
    if d["mode"] == "hpc" and not os.environ.get("LSB_JOBID"):
        raise RuntimeError("HPC worker requires an LSF allocation")
    state = _state_path(plan, key)
    with state.with_suffix(".lock").open("a+") as lock:
        try:
            shared_flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("This plan/job is already running") from exc
        for parent in job["depends"]:
            # Existing LE8 plans may still list the optional deep dependency.
            if ((job.get("job_name") or "").startswith("le8.")
                    and job["key"].endswith(".evaluation")
                    and parent == job["key"].removesuffix(".evaluation") + ".abm_deep"):
                continue
            parent_state = _state_path(plan, parent)
            parent_record = json.loads(parent_state.read_text()) if parent_state.is_file() else {}
            if parent_record:
                check_binding(parent_record, d, "Dependency status")
            if parent_record.get("status") != "DONE":
                raise RuntimeError(f"Dependency did not complete successfully: {parent}")
        env = {**os.environ, **job["env"], "WORKFLOW_WORKER": "1", "SCRIPT_CONSOLE_ACTIVE": "1",
               "submit": d["mode"], "SUBMIT": d["mode"]}
        r = job["resources"]
        env.update(WORKFLOW_MEMORY_GIB=str(r["memory_gib"]), WORKFLOW_CORES=str(min(int(job["env"].get("WORKFLOW_CORES", r["cores"])), r["cores"])),
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


def _read_states(plan, definition=None):
    if definition is None:
        _, definition = load_plan(plan)
    values = {p.stem: json.loads(p.read_text()) for p in (plan.parent / "status").glob("*.json")}
    for record in values.values():
        check_binding(record, definition, "Worker status")
    return values


def check_completed_entries(plan, states):
    for job in plan["jobs"]:
        if (states.get(job["key"], {}).get("status") == "DONE" and
                job["env"].get("WORKFLOW_RESUME_POLICY") == "new_plan"):
            raise RuntimeError("This sequential entry completed. Invoke its public entry again to create a new plan "
                               "and recheck project outputs; the scheduler's DONE records process exit only.")


def _local_one(path, key):
    root = path.parent
    with (root / "logs" / f"{key}.local.out").open("a") as out, (root / "logs" / f"{key}.local.err").open("a") as err:
        result = subprocess.run(["bash", str(root / "commands" / f"{key}.cmd")], stdout=out, stderr=err)
    return key, result.returncode


def relay_local_logs(plan):
    """Keep interactive CLI results/errors visible without a second log copy."""
    path = Path(plan)
    data = json.loads(path.read_text())
    if data["mode"] != "local":
        return
    for job in data["jobs"]:
        for suffix, stream in (("out", sys.stdout), ("err", sys.stderr)):
            log = path.parent / "logs" / f"{job['key']}.local.{suffix}"
            if log.is_file():
                with log.open() as handle:
                    shutil.copyfileobj(handle, stream)


def _run_local(path, plan, jobs, resume):
    from workflow_cpu import local_jobs
    jobs = local_jobs(jobs)
    states = _read_states(path, plan) if resume else {}
    check_completed_entries(plan, states)
    done = {k for k, v in states.items() if v.get("status") == "DONE"}
    pending = {j["key"]: j for j in plan["jobs"] if j["key"] not in done}
    active, failed = {}, {}
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        while pending or active:
            if not failed:
                ready = [j for j in pending.values() if set(j["depends"]) <= done]
                for j in ready[:jobs - len(active)]:
                    active[pool.submit(_local_one, path, j["key"])] = j["key"]
                    pending.pop(j["key"])
            if not active:
                break
            finished, _ = futures_wait(active, return_when=FIRST_COMPLETED)
            for f in finished:
                active.pop(f)
                key, code = f.result()
                if code:
                    failed[key] = code
                else:
                    done.add(key)
                print(f"[{'EXIT' if code else 'DONE'}] {key}", flush=True)
    if failed or pending:
        raise RuntimeError(f"Local workflow failed: {failed}; not started: {list(pending)}. See {path.parent / 'logs'}")
    return {"status": "DONE", "jobs": sorted(done), "plan": str(path),
            "completion_contract": "process_exit_only"}


def _probe(ids):
    result = {}
    for begin in range(0, len(ids), 100):
        batch = ids[begin:begin + 100]
        p = subprocess.run(["bjobs", "-a", "-noheader", "-o", "jobid stat", *batch],
                           capture_output=True, text=True, timeout=30)
        for line in p.stdout.splitlines():
            fields = line.split()
            if len(fields) == 2 and fields[0] in batch:
                result[fields[0]] = fields[1]
    return result


def _wait_hpc(path, receipts, poll=5.0):
    remaining = dict(receipts)
    while remaining:
        local = _read_states(path)
        ids = [v["job_id"] for v in remaining.values()]
        observed = _probe(ids)
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


def _run_hpc(path, plan, jobs, wait, resume, *, dependency_condition=None):
    # Busy LSF controllers may take longer than one minute to acknowledge bsub.
    # A timeout still leaves an uncertain receipt; never retry it automatically.
    submit_timeout = positive(os.getenv("WORKFLOW_BSUB_TIMEOUT_SECONDS", "300"),
                              "WORKFLOW_BSUB_TIMEOUT_SECONDS")
    journal = path.parent / "submission.json"
    previous = json.loads(journal.read_text()) if resume and journal.is_file() else {"jobs": {}}
    if previous["jobs"]:
        check_binding(previous, plan, "Submission journal")
    if resume:
        check_completed_entries(plan, _read_states(path, plan))
    saved = previous["jobs"]
    if set(saved) - {j["key"] for j in plan["jobs"]}:
        raise RuntimeError("Submission journal contains jobs outside this plan")
    # Restore the DAG order when loading the sorted submission journal.
    receipts = {j["key"]: saved[j["key"]] for j in plan["jobs"] if j["key"] in saved}
    if not resume and journal.exists():
        raise RuntimeError("This plan has already been submitted; resume it or create a new plan")
    if any(r.get("status") == "submitting" for r in receipts.values()):
        raise RuntimeError("An earlier bsub result is unknown. Inspect bjobs and submission.json before retrying; no duplicate was submitted")
    tasks = [dict(key=j['key'], path=path.parent/'commands'/f"{j['key']}.cmd",
                  depends=[(key, dependency_condition(j, key) if dependency_condition else 'done')
                           for key in j['depends']]) for j in plan['jobs']]
    def persist(key, row):
        receipts[key] = dict(row, submitted=row['submitted_at'],
                             status='submitted' if row['status'] == 'submitted' else 'submitting')
        atomic_json(journal, dict(plan=str(path), plan_sha256=plan.get('plan_sha256'), jobs=receipts))
    submit_graph(tasks, receipts, persist=persist, receipt=journal, timeout=submit_timeout)
    if wait:
        return _wait_hpc(path, receipts)
    return {"status": "SUBMITTED", "jobs": receipts, "plan": str(path)}


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
            return _run_local(path, d, jobs, resume)
        return _run_hpc(path, d, jobs, wait, resume)


def main(argv=None):
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
        print(json.dumps(result, indent=2))
        return 0
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        p.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
