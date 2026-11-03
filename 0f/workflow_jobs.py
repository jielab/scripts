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
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import tempfile
import time
import uuid

VERSION = "workflow-dag-1"
KEY = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}$")


@dataclass
class Job:
    key: str
    argv: list[str]
    depends: list[str] = field(default_factory=list)
    resources: dict = field(default_factory=dict)
    env: dict = field(default_factory=dict)
    cwd: str | None = None


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, sort_keys=True, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


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


def write_plan(jobs, out_dir, *, mode="hpc", name="workflow", common_env=None):
    if mode not in {"local", "hpc"} or not KEY.fullmatch(name):
        raise ValueError("Invalid mode or plan name")
    shared = {k:v for k,v in os.environ.items() if k.startswith(('DIR_', 'HPC_', '_OS_DEFAULT_')) or k in ('ANALYSIS_ROOT', 'DATE_SCRATCH')}
    ordered = normalized(jobs, {**shared, **(common_env or {})})
    root = Path(out_dir).expanduser().resolve()
    text(root, "plan directory")
    fingerprint = hashlib.sha256(json.dumps(ordered, sort_keys=True).encode()).hexdigest()
    root = root / f"{name}-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')}-{uuid.uuid4().hex[:8]}"
    root.mkdir(parents=True, exist_ok=False)
    for folder in ("commands", "logs", "status"):
        (root / folder).mkdir()
    path = root / "plan.json"
    value = dict(schema=VERSION, created=utcnow(), name=name, mode=mode,
                 configuration_sha256=fingerprint, jobs=ordered)
    atomic_json(path, value)
    helper = str(Path(__file__).resolve())
    worker_python = os.environ.get("WORKFLOW_PYTHON", "python3")
    for j in ordered:
        r = j["resources"]
        minutes = math.ceil(r["hours"] * 60)
        mem = math.ceil(r["memory_gib"] / r["cores"] if r["memory_scope"] == "task" else r["memory_gib"])
        suffix = "/host" if r["memory_scope"] == "host" else ""
        lines = ["#!/usr/bin/env bash", f"#BSUB -J {name}.{j['key']}", f"#BSUB -q {r['queue']}",
                 f"#BSUB -n {r['cores']}", f"#BSUB -W {minutes//60}:{minutes%60:02d}",
                 f"#BSUB -M {math.ceil(r['memory_gib'])}GB",
                 f'#BSUB -R "span[hosts=1] rusage[mem={mem}GB{suffix}]"',
                 f"#BSUB -o {shlex.quote(str(root / 'logs' / (j['key'] + '.%J.out')))}",
                 f"#BSUB -e {shlex.quote(str(root / 'logs' / (j['key'] + '.%J.err')))}"]
        if r["gpus"]:
            gpu = f"num={r['gpus']}:mode=exclusive_process"
            if r.get("gpu_model"):
                gpu += f":gmodel={r['gpu_model']}"
            lines.append(f'#BSUB -gpu "{gpu}"')
        lines.extend(["", "set -euo pipefail", "umask 077"])
        if mode == "hpc":
            lines.extend(['[[ -n "${LSB_JOBID:-}" ]] || { echo "Submit this .cmd with bsub; do not run on a login node." >&2; exit 2; }',
                          f'[[ ${{LSB_DJOB_NUMPROC:-0}} -ge {r["cores"]} ]] || {{ echo "Insufficient allocated CPU slots" >&2; exit 2; }}'])
        lines.append(f"# Dependencies: {', '.join(j['depends']) or '(none)'}")
        lines.append("# Scientific command: " + shlex.join(j["argv"]))
        lines.append(shlex.join(["exec", worker_python, helper, "worker", "--plan", str(path), "--key", j["key"]]))
        command = root / "commands" / (j["key"] + ".cmd")
        command.write_text("\n".join(lines) + "\n", encoding="utf-8")
        command.chmod(0o750)
    print(f"[PLAN] {path} ({len(ordered)} jobs; {mode})", flush=True)
    return path


def load_plan(path):
    path = Path(path).resolve()
    d = json.loads(path.read_text())
    if d.get("schema") != VERSION:
        raise ValueError("Unsupported workflow plan")
    d["jobs"] = normalized(d["jobs"])
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
    with state.with_suffix(".lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("This plan/job is already running") from exc
        for parent in job["depends"]:
            parent_state = _state_path(plan, parent)
            if not parent_state.is_file() or json.loads(parent_state.read_text()).get("status") != "DONE":
                raise RuntimeError(f"Dependency did not complete successfully: {parent}")
        env = {**os.environ, **job["env"], "WORKFLOW_WORKER": "1", "submit": d["mode"], "SUBMIT": d["mode"]}
        r = job["resources"]
        env.update(WORKFLOW_MEMORY_GIB=str(r["memory_gib"]), WORKFLOW_CORES=str(r["cores"]),
                   WORKFLOW_GPUS=str(r["gpus"]), WORKFLOW_RESOURCE_MANAGER="lsf" if d["mode"] == "hpc" else "local")
        env.update(OMP_NUM_THREADS=str(r["cores"]), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                   NUMEXPR_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1")
        # A CPU allocation must never borrow GPUs from another LSF job. Locally,
        # explicit auto can discover a GPU without pretending LSF allocated one.
        if not r["gpus"] and not (d["mode"] == "local" and env.get("WORKFLOW_DEVICE") == "auto"):
            env["CUDA_VISIBLE_DEVICES"] = ""
        temporary = plan.parent / "tmp" / key
        temporary.mkdir(parents=True, exist_ok=True)
        env.update(TMPDIR=str(temporary), TMP=str(temporary), TEMP=str(temporary))
        record = dict(key=key, status="RUN", started=utcnow(), job_id=os.getenv("LSB_JOBID"), pid=os.getpid())
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


def _read_states(plan):
    return {p.stem: json.loads(p.read_text()) for p in (plan.parent / "status").glob("*.json")}


def _local_one(path, key):
    root = path.parent
    with (root / "logs" / f"{key}.local.out").open("a") as out, (root / "logs" / f"{key}.local.err").open("a") as err:
        result = subprocess.run(["bash", str(root / "commands" / f"{key}.cmd")], stdout=out, stderr=err)
    return key, result.returncode


def _run_local(path, plan, jobs, resume):
    states = _read_states(path) if resume else {}
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
    return {"status": "DONE", "jobs": sorted(done), "plan": str(path)}


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


def _run_hpc(path, plan, jobs, wait, resume):
    journal = path.parent / "submission.json"
    previous = json.loads(journal.read_text()) if resume and journal.is_file() else {"jobs": {}}
    receipts = previous["jobs"]
    if not resume and journal.exists():
        raise RuntimeError("This plan has already been submitted; resume it or create a new plan")
    if any(r.get("status") == "submitting" for r in receipts.values()):
        raise RuntimeError("An earlier bsub result is unknown. Inspect bjobs and submission.json before retrying; no duplicate was submitted")
    for j in plan["jobs"]:
        key = j["key"]
        if key in receipts:
            continue
        condition = [f"done({receipts[k]['job_id']})" for k in j["depends"]]
        # Submission lanes cap simultaneous jobs without running a waiting job
        # on compute resources. Lane dependencies are ended(), scientific DAG
        # dependencies are done(); a failed unrelated job is not a success gate.
        if len(receipts) >= jobs:
            lane = list(receipts.values())[-jobs]
            condition.append(f"ended({lane['job_id']})")
        argv = ["bsub"]
        if condition:
            argv += ["-w", " && ".join(dict.fromkeys(condition))]
        receipts[key] = dict(status="submitting", argv=argv, submitted=utcnow())
        atomic_json(journal, dict(plan=str(path), jobs=receipts))
        try:
            with (path.parent / "commands" / f"{key}.cmd").open() as handle:
                p = subprocess.run(argv, stdin=handle, capture_output=True, text=True,
                                   timeout=60, env={**os.environ, "LC_ALL": "C"})
        except BaseException:
            # Keep the unknown submission receipt; never retry blindly.
            raise
        found = re.search(r"Job <([0-9]+)>", p.stdout + p.stderr)
        receipts[key].update(returncode=p.returncode, stdout=p.stdout, stderr=p.stderr)
        if p.returncode or not found:
            # A server might have accepted a job even if the client failed.
            atomic_json(journal, dict(plan=str(path), jobs=receipts))
            raise RuntimeError(f"LSF submission result unknown/failed for {key}; inspect {journal}")
        receipts[key].update(status="submitted", job_id=found[1])
        atomic_json(journal, dict(plan=str(path), jobs=receipts))
        print(f"[SUBMITTED] {key} {found[1]}", flush=True)
    if wait:
        return _wait_hpc(path, receipts)
    return {"status": "SUBMITTED", "jobs": receipts, "plan": str(path)}


def run_plan(plan_path, *, submit="hpc", jobs=128, wait=False, dry_run=False, resume=True):
    path, d = load_plan(plan_path)
    jobs = positive(jobs, "jobs", True)
    if jobs > 10000 or submit != d["mode"]:
        raise ValueError("Invalid concurrency or plan/submit mode mismatch")
    if dry_run:
        return {"status": "PLANNED", "plan": str(path), "jobs": len(d["jobs"])}
    with (path.parent / "dispatch.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
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
    r.add_argument("--submit", choices=["local", "hpc"], required=True)
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
