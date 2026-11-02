#!/usr/bin/env python3
"""Run an explicit DAG locally or submit reviewable per-job LSF command files.

Planning/submission never runs a scientific payload in HPC mode. A durable
pre-submission intent prevents ambiguous scheduler responses from being retried.
The local parallelism setting is independent of the LSF concurrency ceiling.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, wait as futures_wait, FIRST_COMPLETED
from contextlib import contextmanager
from decimal import Decimal, InvalidOperation, ROUND_CEILING
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone

SCHEMA = "allele3-batch-plan-1"
ID_RE = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]*\Z")
ACTIVE = {"PEND", "RUN", "WAIT", "PSUSP", "USUSP", "SSUSP"}
METADATA_JSON_LIMIT = 16 * 1024 * 1024
STAGE_VALIDATION_SCHEMA = "allele3-scientific-stage-1"


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n"); handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def text_value(value, label):
    if not isinstance(value, (str, os.PathLike)):
        raise ValueError(f"{label} must be text")
    value = os.fspath(value)
    if any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError(f"Control character in {label}")
    return value


def integer(value, label, minimum=1, maximum=None):
    if isinstance(value, bool): raise ValueError(f"Invalid {label}")
    try:
        result = int(value)
    except (TypeError, ValueError): raise ValueError(f"Invalid {label}") from None
    if str(result) != str(value) or result < minimum or (maximum is not None and result > maximum):
        raise ValueError(f"Invalid {label}: {value}")
    return result


def positive(value, label):
    if isinstance(value, bool): raise ValueError(f"Invalid {label}")
    try: result = Decimal(str(value))
    except (InvalidOperation, TypeError): raise ValueError(f"Invalid {label}") from None
    if not result.is_finite() or result <= 0: raise ValueError(f"Invalid {label}")
    return result


def absolute(value, label):
    return str(Path(text_value(value, label)).expanduser().absolute())


def path_identity(path, full=False):
    path = Path(path)
    if not path.is_file(): return {"path": str(path), "status": "missing_or_not_file"}
    stat = path.stat()
    result = {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if full:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024*1024), b""): digest.update(chunk)
        result["sha256"] = digest.hexdigest()
    if (path.stat().st_size, path.stat().st_mtime_ns) != (stat.st_size, stat.st_mtime_ns):
        raise ValueError(f"Input changed while reading identity: {path}")
    return result


def normalize_resources(raw, plan):
    resource = dict(plan.get("resources", {})); resource.update(raw or {})
    queue = text_value(resource.get("queue", "ser"), "queue")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.+-]*", queue): raise ValueError("Invalid queue")
    requested = integer(resource.get("cores", 4), "cores")
    cores = max(requested, 40) if queue == "short" else requested
    memory = positive(resource.get("memory_gib", 16), "memory_gib")
    minutes = int((positive(resource.get("hours", 24), "hours")*60).to_integral_value(rounding=ROUND_CEILING))
    gpu = integer(resource.get("gpu", 0), "gpu", minimum=0)
    scope = plan.get("lsf_memory_scope", "host")
    per_task = plan.get("lsf_resource_reserve_per_task")
    if scope == "host": unit = "GB/host"; amount = memory
    elif scope == "site_host" and per_task is False: unit = "GB"; amount = memory
    elif scope == "site_task" and per_task is True: unit = "GB"; amount = memory/cores
    else:
        raise ValueError("Use lsf_memory_scope=host, or confirm site_host/site_task with lsf_resource_reserve_per_task=false/true")
    return dict(queue=queue, cores=cores, requested_cores=requested, memory_gib=str(memory),
                minutes=minutes, gpu=gpu, lsf_memory=f"{int(amount.to_integral_value(rounding=ROUND_CEILING))}{unit}")


def normalize_plan(plan, mode):
    if mode not in {"hpc", "local"}: raise ValueError("mode must be local or hpc")
    if not isinstance(plan, dict): raise ValueError("Plan must be a JSON object")
    name = text_value(plan.get("name", "allele3"), "name")
    if not ID_RE.fullmatch(name): raise ValueError("Invalid plan name")
    jobs = plan.get("jobs")
    if not isinstance(jobs, list) or not jobs: raise ValueError("Plan requires at least one job")
    normalized = []
    for raw in jobs:
        if not isinstance(raw, dict): raise ValueError("Job must be an object")
        job_id = text_value(raw.get("id", ""), "job id")
        if not ID_RE.fullmatch(job_id): raise ValueError("Invalid job id")
        argv = raw.get("argv")
        if not isinstance(argv, list) or not argv: raise ValueError("argv must be a nonempty argument array")
        argv = [text_value(v, "argv") for v in argv]
        if not argv[0]: raise ValueError("Empty executable")
        deps = raw.get("deps", [])
        if not isinstance(deps, list) or any(not isinstance(d, str) for d in deps): raise ValueError("deps must be job IDs")
        if len(set(deps)) != len(deps): raise ValueError("Duplicate dependency")
        env = raw.get("env", {})
        if not isinstance(env, dict): raise ValueError("env must be an object")
        for key in env:
            if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key): raise ValueError("Invalid environment key")
        env = {key:text_value(value, "env value") for key,value in env.items()}
        outputs = raw.get("outputs", [])
        if not isinstance(outputs, list): raise ValueError("outputs must be a list of completed files")
        normalized.append(dict(id=job_id, argv=argv, deps=deps, cwd=absolute(raw.get("cwd", os.getcwd()), "cwd"),
            env=env, resources=normalize_resources(raw.get("resources"), plan),
            outputs=[absolute(p,"output") for p in outputs]))
    lookup = {j["id"]:j for j in normalized}
    if len(lookup) != len(normalized): raise ValueError("Duplicate job IDs")
    for job in normalized:
        if any(d not in lookup or d == job["id"] for d in job["deps"]): raise ValueError("Unknown or self dependency")
    ordered, remaining = [], list(normalized)
    completed = set()
    while remaining:
        job = next((j for j in remaining if set(j["deps"]) <= completed), None)
        if job is None: raise ValueError("Cyclic dependency graph")
        ordered.append(job); completed.add(job["id"]); remaining.remove(job)
    max_running = integer(plan.get("max_running", 200), "max_running")
    for i, job in enumerate(ordered):
        # Artificial lane edges only limit running concurrency; they are not
        # data dependencies, and an EXIT releases its slot without satisfying
        # any downstream scientific dependency.
        job["throttle_after"] = ordered[i-max_running]["id"] if i >= max_running else None
    inputs = []
    for item in plan.get("inputs", []):
        item = {"path":item} if isinstance(item,str) else item
        if not isinstance(item, dict) or "path" not in item: raise ValueError("Invalid declared input")
        path = absolute(item["path"], "input")
        inputs.append(path_identity(path, bool(item.get("sha256", False))))
    result = dict(schema=SCHEMA, name=name, mode=mode, max_running=max_running,
        local_memory_gib=str(positive(plan.get("local_memory_gib", 32), "local_memory_gib")),
        driver_python=text_value(plan.get("driver_python", sys.executable), "driver_python"),
        input_identity=inputs, provenance=plan.get("provenance", {}), jobs=ordered,
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    result["identity"] = sha(result)
    return result


def validate_frozen(plan):
    content = {k:v for k,v in plan.items() if k != "identity"}
    if plan.get("schema") != SCHEMA or sha(content) != plan.get("identity"):
        raise ValueError("Frozen plan identity failed")


def verify_inputs(plan):
    if plan["runner_sha256"] != hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
        raise ValueError("Batch runner changed since the plan was frozen; use a new submission directory")
    for expected in plan.get("input_identity", []):
        if path_identity(expected["path"], "sha256" in expected) != expected:
            raise ValueError("Declared source input changed since planning: " + expected["path"] + "; use a new submission directory")


@contextmanager
def locked(out):
    out.mkdir(parents=True, exist_ok=True)
    with (out/".dispatch.lock").open("a") as handle:
        try: fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError("Another dispatcher owns this plan directory") from None
        try: yield
        finally: fcntl.flock(handle, fcntl.LOCK_UN)


def shell_text(plan, job, out):
    resources = job["resources"]; stem = job["id"]
    label = f"{plan['name'][:30]}.{plan['identity'][:10]}.{stem[:65]}"
    lines = ["#!/usr/bin/env bash", f'#BSUB -J "{label}"', f"#BSUB -q {resources['queue']}",
        f"#BSUB -n {resources['cores']}", f"#BSUB -W {resources['minutes']//60}:{resources['minutes']%60:02d}",
        f'#BSUB -R "span[hosts=1] rusage[mem={resources["lsf_memory"]}]"',
        "#BSUB -o " + shlex.quote(str(out/"logs"/(stem+".%J.out"))),
        "#BSUB -e " + shlex.quote(str(out/"logs"/(stem+".%J.err")))]
    if resources["gpu"]: lines.append(f'#BSUB -gpu "num={resources["gpu"]}:mode=exclusive_process"')
    lines += ["", "set -euo pipefail"]
    if plan["mode"] == "hpc":
        lines += ["[[ ${LSB_JOBID:-} =~ ^[1-9][0-9]*$ ]] || { echo 'Submit this .cmd with bsub; no login-node payload' >&2; exit 2; }",
                  f"[[ ${{LSB_DJOB_NUMPROC:-0}} =~ ^[0-9]+$ && ${{LSB_DJOB_NUMPROC:-0}} -ge {resources['cores']} ]] || {{ echo 'Insufficient allocated CPU slots' >&2; exit 2; }}"]
    else:
        lines += ["if [[ -z ${LSB_JOBID:-} ]] && command -v bsub >/dev/null 2>&1; then echo 'Local payload refused on an LSF submission host without an allocation' >&2; exit 2; fi"]
    for key, value in sorted(job["env"].items()): lines.append(f"export {key}={shlex.quote(value)}")
    lines += [f"export submit={plan['mode']} SUBMIT={plan['mode']} ALLELE_EXECUTION_MODE={plan['mode']}",
              "export WORKFLOW_WORKER=1 ALLELE_IN_COMPUTE=1 ALLELE_BATCH_WORKER=1",
              "export ALLELE_ALLOCATED_MEMORY_GIB=" + shlex.quote(resources["memory_gib"]),
              f"export ALLELE_ALLOCATED_CORES={resources['cores']}",
              "export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1",
              "cd -- " + shlex.quote(job["cwd"]),
              "# Payload (argument quoting shown for review): " + shlex.join(job["argv"]),
              "exec " + shlex.join([plan["driver_python"], str(Path(__file__).resolve()), "worker",
                      "--frozen-plan",str(out/"plan.json"),"--job-id",stem,"--out",str(out)]), ""]
    return "\n".join(lines)


def prepare(plan, out):
    out = Path(out).absolute()
    frozen_path = out/"plan.json"
    if frozen_path.exists():
        current = json.loads(frozen_path.read_text()); validate_frozen(current)
        if current["identity"] != plan["identity"]:
            raise ValueError("Plan/input identity changed; use a different submission directory; old job IDs are not reusable")
    else: atomic_json(frozen_path, plan)
    for folder in ("cmd", "logs", "runs"): (out/folder).mkdir(exist_ok=True)
    commands = {}
    for job in plan["jobs"]:
        path = out/"cmd"/(job["id"]+".cmd")
        text = shell_text(plan, job, out)
        if path.exists() and path.read_text() != text:
            raise ValueError(f"Generated command changed: {path}; use a new submission directory")
        if not path.exists(): path.write_text(text); path.chmod(0o750)
        commands[job["id"]] = {"path":str(path),"sha256":hashlib.sha256(text.encode()).hexdigest()}
    dispatcher = "\n".join(["#!/usr/bin/env bash","set -euo pipefail",
        "# Re-enters the receipt-aware dispatcher; never blindly resubmits each .cmd.",
        "exec " + shlex.join([plan["driver_python"], str(Path(__file__).resolve()), "dispatch",
            "--frozen-plan",str(frozen_path),"--out",str(out)]) + ' "$@"',""])
    path = out/"submit.sh"
    if path.exists() and path.read_text() != dispatcher: raise ValueError("Dispatcher changed; use a new submission directory")
    if not path.exists(): path.write_text(dispatcher); path.chmod(0o750)
    return commands


def read_state(out, plan):
    path=out/"submission.json"
    if path.exists():
        state=json.loads(path.read_text())
        if state.get("plan_identity") != plan["identity"]: raise ValueError("Submission receipt differs from frozen plan")
        return state
    return dict(schema=SCHEMA, plan_identity=plan["identity"], mode=plan["mode"], created_at=utcnow(), jobs={})


def save_state(out, state):
    state["updated_at"]=utcnow(); atomic_json(out/"submission.json",state)
    lines=["job\tjob_id\tstatus\tdependency\n"]
    for name,row in state["jobs"].items():
        lines.append("\t".join(str(row.get(k,"")) for k in ("name","job_id","status","dependency"))+"\n")
    path=out/"job_ids.tsv"; temporary=path.with_suffix(".tsv.partial")
    temporary.write_text("".join(lines)); os.replace(temporary,path)


def dependency(job, state):
    terms=[]
    for name in job["deps"]:
        record=state["jobs"].get(name,{})
        if record.get("status") != "submitted" or not isinstance(record.get("job_id"),int):
            raise RuntimeError("Dependency has no confirmed scheduler ID: "+name)
        terms.append(f"done({record['job_id']})")
    lane=job["throttle_after"]
    if lane and lane not in job["deps"]:
        record=state["jobs"].get(lane,{})
        if record.get("status") != "submitted" or not isinstance(record.get("job_id"),int):
            raise RuntimeError("Concurrency predecessor has no scheduler ID: "+lane)
        terms.append(f"ended({record['job_id']})")
    return " && ".join(terms)


def decoded(value):
    return value.decode(errors="replace") if isinstance(value,bytes) else value or ""


def submit_hpc(plan,out,commands,state,bsub="bsub",timeout=30):
    for job in plan["jobs"]:
        name=job["id"]; old=state["jobs"].get(name)
        if old and old.get("status") == "submitted":
            if old.get("command_sha256") != commands[name]["sha256"]: raise ValueError("Recorded command identity changed")
            print(f"[batch] reuse submitted job={name} id={old['job_id']}",flush=True)
            continue
        if old and old.get("status") in {"submitting","submission_unknown"}:
            raise RuntimeError(f"Submission outcome unknown for {name}; inspect scheduler and receipt before resolving; never auto-resubmitted")
        condition=dependency(job,state)
        argv=[bsub]
        if condition: argv += ["-w",condition]
        record=dict(name=name,status="submitting",job_id=None,dependency=condition,argv=argv,
            command_sha256=commands[name]["sha256"],attempt=(old or {}).get("attempt",0)+1,
            started_at=utcnow(),history=(old or {}).get("history",[]) +
            ([{k:v for k,v in old.items() if k!="history"}] if old else []))
        state["jobs"][name]=record; save_state(out,state)
        result=None
        try:
            with Path(commands[name]["path"]).open("rb") as handle:
                result=subprocess.run(argv,stdin=handle,capture_output=True,timeout=timeout,
                    env={**os.environ,"LC_ALL":"C"})
            record.update(returncode=result.returncode,stdout=decoded(result.stdout),stderr=decoded(result.stderr))
        except subprocess.TimeoutExpired as exc:
            record.update(returncode=None,stdout=decoded(exc.stdout),stderr=decoded(exc.stderr),error="bsub timed out; acceptance unknown")
        except OSError as exc:
            record.update(returncode=None,stdout="",stderr=str(exc),error="bsub could not be started",status="submit_failed")
        ids=set(re.findall(r"\bJob\s*<([1-9][0-9]*)>",record.get("stdout","")+"\n"+record.get("stderr","")))
        if len(ids)==1: record["job_id"]=int(next(iter(ids)))
        if result is not None and result.returncode==0 and len(ids)==1:
            record["status"]="submitted"
        elif result is not None and result.returncode!=0 and not ids:
            record["status"]="submit_failed"
        elif record["status"] != "submit_failed": record["status"]="submission_unknown"
        record["finished_at"]=utcnow(); save_state(out,state)
        if record["status"] != "submitted":
            raise RuntimeError(f"{name}: {record['status']}; inspect {out/'submission.json'}; accepted job IDs were retained")
        print(f"[batch] submitted job={name} id={record['job_id']} dependency={condition or 'none'}",flush=True)
    return state


def read_metadata_json(path):
    """Bound metadata reads; never inspect the contents of genomic products."""
    path=Path(path)
    if path.stat().st_size>METADATA_JSON_LIMIT:
        raise ValueError("Metadata JSON exceeds the bounded receipt check: "+str(path))
    with path.open() as handle:
        text=handle.read(METADATA_JSON_LIMIT+1)
    if len(text)>METADATA_JSON_LIMIT:
        raise ValueError("Metadata JSON grew during the receipt check: "+str(path))
    return json.loads(text)


def verify_stat_identity(expected,label="Scientific file identity changed: ",allow_missing=False):
    """Verify saved size/mtime identities without rehashing large files."""
    if not isinstance(expected,dict) or not isinstance(expected.get("path"),str):
        raise ValueError("Invalid recorded scientific file identity")
    actual=path_identity(expected["path"])
    if allow_missing and expected.get("status") in {"missing","missing_or_not_file"}:
        if "size" not in actual: return
    elif all(isinstance(expected.get(k),int) and not isinstance(expected[k],bool) for k in ("size","mtime_ns")) and all(actual.get(k)==expected[k] for k in ("size","mtime_ns")):
        return
    raise ValueError(label+expected["path"]+"; use a new submission directory for recomputation")


def validate_stage_receipt(path,seen):
    saved=read_metadata_json(path)
    snapshot=saved.get("identity",{})
    inputs=snapshot.get("inputs")
    outputs=saved.get("outputs")
    if not isinstance(inputs,list) or not isinstance(outputs,dict) or not outputs:
        raise ValueError("Scientific stage receipt lacks an input/output inventory: "+str(path))
    for expected in inputs:
        verify_stat_identity(expected,allow_missing=True)
    for output,expected in outputs.items():
        if not isinstance(expected,dict) or expected.get("path")!=output:
            raise ValueError("Inconsistent scientific output identity: "+str(path))
        verify_stat_identity(expected)
        validate_output_file(output,seen)


def validate_output_file(path,seen=None):
    """Recognize scientific manifests/ledgers; ordinary markers remain markers.

    This may run on the dispatcher: only bounded JSON and file stats are read.
    Current analysis_integrity inventories contain size/mtime. A checksum-only
    inventory is refused here rather than silently rehashing a huge product.
    """
    path=Path(path)
    if path.suffix.lower()!=".json": return
    seen=set() if seen is None else seen
    key=str(path.resolve())
    if key in seen: return
    seen.add(key)
    if path.stat().st_size>METADATA_JSON_LIMIT:
        if "manifest" in path.name or ".batch_stage_status" in path.parts:
            raise ValueError("Scientific metadata exceeds the bounded receipt check: "+str(path))
        return
    try: saved=read_metadata_json(path)
    except (json.JSONDecodeError,UnicodeError): return
    if not isinstance(saved,dict): return
    if "output_inventory" in saved:
        from analysis_integrity import validate_manifest
        inventory=saved["output_inventory"]
        if not isinstance(inventory,dict) or any(not isinstance(value,dict) for value in inventory.values()):
            raise ValueError("Invalid scientific output inventory: "+str(path))
        if any(not all(isinstance(value.get(k),int) and not isinstance(value[k],bool) for k in ("size","mtime_ns"))
               for value in inventory.values()):
            raise ValueError("Scientific inventory lacks verifiable size/mtime identities: "+str(path))
        if any("sha256" in value for value in inventory.values()):
            raise ValueError("Scientific inventory requires checksum verification in compute; no large-file hashing on dispatcher: "+str(path))
        validate_manifest(path.parent,path.name)
        inputs=saved.get("inputs",{})
        if not isinstance(inputs,dict): raise ValueError("Invalid scientific input inventory: "+str(path))
        for expected in inputs.values(): verify_stat_identity(expected,allow_missing=True)
    if saved.get("schema")=="allele3-batch-stage-1":
        proof=saved.get("scientific_validation",{})
        if saved.get("status") not in {"complete","reused"} or proof.get("schema")!=STAGE_VALIDATION_SCHEMA:
            raise ValueError("Stage ledger lacks verified scientific outputs; use a new submission directory: "+str(path))
        snapshot=proof.get("snapshot",{})
        if not isinstance(snapshot.get("inputs"),list):
            raise ValueError("Stage ledger lacks its source snapshot: "+str(path))
        for expected in snapshot["inputs"]: verify_stat_identity(expected,allow_missing=True)
        receipt=proof.get("receipt")
        products=proof.get("outputs",[])
        if receipt:
            verify_stat_identity(receipt)
            validate_stage_receipt(receipt["path"],seen)
        elif not products:
            raise ValueError("Stage ledger has no independently verifiable products: "+str(path))
        for expected in [*products,*proof.get("dependencies",[])]:
            verify_stat_identity(expected)
            validate_output_file(expected["path"],seen)


def run_receipt(out,name,identity):
    path=out/"runs"/(name+".json")
    if not path.exists(): return None
    row=json.loads(path.read_text())
    if row.get("plan_identity") != identity: raise ValueError("Worker receipt differs from plan")
    if row.get("status")=="succeeded":
        for expected in row.get("outputs",[]):
            verify_stat_identity(expected,"Completed worker output changed: ")
            validate_output_file(expected["path"])
        row["output_verification"]="declared_output_metadata_validated" if row.get("outputs") else "unavailable_no_declared_outputs"
    return row


def receipt_for_submission(receipt, record):
    """A previous attempt's receipt cannot complete or fail its replacement."""
    return receipt if receipt and str(receipt.get("lsf_job_id")) == str(record.get("job_id")) else None


def _worker_payload(frozen_plan,job_id,out):
    plan=json.loads(Path(frozen_plan).read_text()); validate_frozen(plan); verify_inputs(plan)
    jobs={j["id"]:j for j in plan["jobs"]}
    if job_id not in jobs: raise ValueError("Unknown worker job")
    job=jobs[job_id]; out=Path(out).absolute()
    if plan["mode"]=="hpc" and not re.fullmatch(r"[1-9][0-9]*",os.environ.get("LSB_JOBID","")):
        raise RuntimeError("HPC payload requires an LSF allocation")
    if plan["mode"]=="local" and not os.environ.get("LSB_JOBID") and shutil.which("bsub"):
        raise RuntimeError("Local payload refused on an LSF submission host without an allocation")
    for name in job["deps"]:
        previous=run_receipt(out,name,plan["identity"])
        if previous is None or previous.get("status")!="succeeded":
            raise RuntimeError("Scientific dependency has no successful worker receipt: "+name)
    path=out/"runs"/(job_id+".json")
    prior=run_receipt(out,job_id,plan["identity"])
    if prior and prior.get("status")=="succeeded":
        if not prior.get("outputs"):
            raise ValueError("Successful execution has no declared outputs to verify for reuse; use a new submission directory")
        prior.update(reused_verified_outputs=True,previous_lsf_job_id=prior.get("lsf_job_id"),
                     lsf_job_id=os.environ.get("LSB_JOBID"),revalidated_at=utcnow())
        atomic_json(path,prior)
        return 0
    row=dict(schema=SCHEMA,plan_identity=plan["identity"],name=job_id,status="running",started_at=utcnow(),
             lsf_job_id=os.environ.get("LSB_JOBID"),argv=job["argv"],cwd=job["cwd"])
    atomic_json(path,row)
    try:
        environment={**os.environ,**job["env"],"submit":plan["mode"],"SUBMIT":plan["mode"],
            "ALLELE_EXECUTION_MODE":plan["mode"],"WORKFLOW_WORKER":"1","ALLELE_BATCH_WORKER":"1", "ALLELE_IN_COMPUTE":"1",
            "ALLELE_ALLOCATED_MEMORY_GIB":job["resources"]["memory_gib"],"ALLELE_ALLOCATED_CORES":str(job["resources"]["cores"])}
        process=subprocess.run(job["argv"],cwd=job["cwd"],env=environment)
        row.update(returncode=process.returncode,status="succeeded" if process.returncode==0 else "failed")
        if process.returncode==0:
            identities=[path_identity(p) for p in job["outputs"]]
            if any("size" not in p for p in identities):
                raise ValueError("Payload exited zero but a declared completed output is missing")
            for expected in identities: validate_output_file(expected["path"])
            row["outputs"]=identities
            row["output_verification"]="declared_output_metadata_validated" if identities else "unavailable_no_declared_outputs"
    except (OSError,ValueError) as exc:
        row.update(returncode=1,status="failed",error=str(exc))
    finally:
        row["finished_at"]=utcnow(); atomic_json(path,row)
    return row.get("returncode",1)


def worker(frozen_plan,job_id,out):
    if not ID_RE.fullmatch(job_id): raise ValueError("Invalid worker job ID")
    out=Path(out).absolute()
    if Path(frozen_plan).absolute() != out/"plan.json": raise ValueError("Worker plan must belong to its output directory")
    (out/"runs").mkdir(parents=True,exist_ok=True)
    with (out/"runs"/(job_id+".lock")).open("a") as handle:
        try: fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError("Another worker owns this exact job; no duplicate payload started") from None
        try: return _worker_payload(frozen_plan,job_id,out)
        finally: fcntl.flock(handle,fcntl.LOCK_UN)


def local_budget(plan):
    limits=[float(plan["local_memory_gib"])]
    try:
        available=os.sysconf("SC_AVPHYS_PAGES")*os.sysconf("SC_PAGE_SIZE")/1024**3
        if available>0: limits.append(available)
    except (OSError,ValueError): pass
    for name in ("/sys/fs/cgroup/memory.max","/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            value=Path(name).read_text().strip()
            if value.isdigit() and 0<int(value)<2**60: limits.append(int(value)/1024**3)
        except OSError: pass
    return min(limits)


def local_command(command,out,name):
    with (out/"logs"/(name+".local.out")).open("ab") as stdout, (out/"logs"/(name+".local.err")).open("ab") as stderr:
        return subprocess.run(["bash",command],stdout=stdout,stderr=stderr).returncode


def execute_local(plan,out,commands,state,jobs_n=1):
    jobs_n=integer(jobs_n,"jobs_N",maximum=8)
    if not os.environ.get("LSB_JOBID") and shutil.which("bsub"):
        raise RuntimeError("Local mode refused on an LSF submission host without an allocation; select hpc")
    budget=local_budget(plan); by_id={j["id"]:j for j in plan["jobs"]}
    pending=set(by_id); finished={}; running={}; occupied=0.
    for name in list(pending):
        prior=run_receipt(out,name,plan["identity"])
        if prior and prior.get("status")=="succeeded":
            if not prior.get("outputs"):
                raise ValueError("Successful execution has no declared outputs to verify for reuse; use a new submission directory")
            finished[name]=True; pending.remove(name)
        elif prior and prior.get("status")=="running":
            raise RuntimeError("A prior worker has an unresolved running receipt: "+name)
    with ThreadPoolExecutor(max_workers=jobs_n) as pool:
        while pending or running:
            for job in plan["jobs"]:
                name=job["id"]
                if name not in pending: continue
                if any(d in finished and not finished[d] for d in job["deps"]):
                    state["jobs"][name]=dict(name=name,status="blocked_dependency",job_id=None,dependency=",".join(job["deps"]))
                    finished[name]=False; pending.remove(name); save_state(out,state); continue
                if not all(d in finished for d in job["deps"]): continue
                memory=float(job["resources"]["memory_gib"])
                if memory>budget:
                    raise RuntimeError(f"Local job {name} requests {memory:g} GiB, exceeding available/configured budget {budget:g} GiB")
                if len(running)>=jobs_n or occupied+memory>budget: continue
                state["jobs"][name]=dict(name=name,status="running",job_id=None,dependency=",".join(job["deps"]),command_sha256=commands[name]["sha256"])
                save_state(out,state); pending.remove(name); occupied+=memory
                future=pool.submit(local_command,commands[name]["path"],out,name); running[future]=(name,memory)
            if not running:
                if pending: raise RuntimeError("No runnable local jobs")
                break
            done,_=futures_wait(running,return_when=FIRST_COMPLETED)
            for future in done:
                name,memory=running.pop(future); occupied-=memory
                try: code=future.result()
                except OSError as exc: code=1; state["jobs"][name]["error"]=str(exc)
                receipt=run_receipt(out,name,plan["identity"])
                succeeded=code==0 and receipt is not None and receipt.get("status")=="succeeded"
                finished[name]=succeeded
                state["jobs"][name].update(status="succeeded" if succeeded else "failed",returncode=code)
                save_state(out,state)
    state["status"]="complete" if all(finished.values()) else "failed"; save_state(out,state)
    if state["status"]!="complete": raise RuntimeError("One or more local jobs failed; downstream jobs were blocked")
    return state


def scheduler_states(state,bjobs):
    ids=[str(row["job_id"]) for row in state["jobs"].values() if row.get("status")=="submitted"]
    if not ids: return {}
    result=subprocess.run([bjobs,"-a","-noheader","-o","jobid stat",*ids],capture_output=True,text=True,timeout=30)
    parsed={}
    for line in result.stdout.splitlines():
        fields=line.split()
        if len(fields)==2 and fields[0].isdigit() and fields[1] in ACTIVE|{"DONE","EXIT","UNKWN","ZOMBI"}:
            parsed[int(fields[0])]=fields[1]
    if result.returncode and not parsed: raise RuntimeError("Cannot establish scheduler state: "+result.stderr.strip())
    return parsed


def change_pending_dependency(state,out,context,name,condition,phase,bjobs,bmod,timeout):
    record=state["jobs"][name]; job_id=record["job_id"]
    observed=scheduler_states(state,bjobs).get(job_id)
    if observed!="PEND":
        raise RuntimeError(f"Refusing bmod for {name}: scheduler state={observed}; only a confirmed PEND job may be changed")
    action=dict(name=name,job_id=job_id,phase=phase,condition=condition,status="changing",started_at=utcnow())
    context.setdefault("actions",[]).append(action); save_state(out,state)
    try:
        process=subprocess.run([bmod,"-w",condition,str(job_id)],capture_output=True,text=True,timeout=timeout,
                               env={**os.environ,"LC_ALL":"C"})
        action.update(returncode=process.returncode,stdout=process.stdout,stderr=process.stderr,
                      status="applied" if process.returncode==0 else "failed")
    except (OSError,subprocess.TimeoutExpired) as exc:
        action.update(status="outcome_unknown",error=str(exc)); save_state(out,state)
        raise RuntimeError(f"bmod outcome unresolved for {name}; accepted retry IDs remain recorded") from exc
    action["finished_at"]=utcnow(); save_state(out,state)
    if action["status"]!="applied": raise RuntimeError(f"bmod failed for {name}; inspect the retained retry actions")
    if phase=="hold" and scheduler_states(state,bjobs).get(job_id)!="PEND":
        raise RuntimeError(f"{name} did not remain PEND after hold; no new retry job will be submitted")
    record["dependency"]=condition
    context["held" if phase=="hold" else "released"][name]=condition
    save_state(out,state)


def retry_failed_hpc(plan,out,commands,state,bsub,bjobs,bmod,timeout):
    """Explicitly retry confirmed EXIT jobs and rewire their own pending edges.

    Pending successors are held before a failed lane is restarted, preserving
    concurrency and preventing them from retaining old numeric dependency IDs.
    No job is killed. Accepted retry IDs and each bmod action are durable.
    """
    if not shutil.which(bmod): raise RuntimeError("--retry-failed needs bmod; nothing was resubmitted")
    context=state.get("retry")
    if context is None or context.get("status")=="complete":
        if any(r.get("status") in {"submitting","submission_unknown"} for r in state["jobs"].values()):
            raise RuntimeError("Resolve unknown submission outcomes before retrying; no job was resubmitted")
        observed=scheduler_states(state,bjobs)
        confirmed={n:r for n,r in state["jobs"].items() if r.get("status")=="submitted"}
        for name,row in confirmed.items():
            status=observed.get(row["job_id"])
            if status not in {"PEND","DONE","EXIT"}:
                raise RuntimeError(f"Retry requires settled known states; {name} is {status}; running/unknown jobs are untouched")
            if status=="DONE":
                receipt=receipt_for_submission(run_receipt(out,name,plan["identity"]),row)
                if receipt is None or receipt.get("status")!="succeeded":
                    raise RuntimeError("DONE lacks a valid successful receipt: "+name)
        roots={name:{"old_id":r["job_id"]} for name,r in confirmed.items() if observed.get(r["job_id"])=="EXIT"}
        if not roots: return
        affected=[j["id"] for j in plan["jobs"] if j["id"] in confirmed and observed[confirmed[j["id"]]["job_id"]]=="PEND"
                  and (set(j["deps"])|({j["throttle_after"]} if j["throttle_after"] else set())) & set(roots)]
        if context: state.setdefault("retry_history",[]).append(context)
        context=dict(status="holding",started_at=utcnow(),roots=roots,affected=affected,held={},released={},actions=[])
        state["retry"]=context; save_state(out,state)
    by_id={j["id"]:j for j in plan["jobs"]}
    hold_condition=f"done({next(iter(context['roots'].values()))['old_id']})"
    for name in context["affected"]:
        if name not in context["held"]:
            change_pending_dependency(state,out,context,name,hold_condition,"hold",bjobs,bmod,timeout)
    context["status"]="resubmitting"; save_state(out,state)
    for job in plan["jobs"]:
        name=job["id"]
        if name not in context["roots"]: continue
        root=context["roots"][name]; record=state["jobs"][name]
        if record.get("status")=="submitted" and record["job_id"]!=root["old_id"]:
            root["new_id"]=record["job_id"]; continue
        if record.get("status")=="submitted":
            # The hold stage finished before changing any of the old IDs.
            if scheduler_states(state,bjobs).get(root["old_id"])!="EXIT":
                raise RuntimeError("Previously failed job is no longer confirmed EXIT: "+name)
            record["status"]="retry_pending"; save_state(out,state)
        subset=dict(plan,jobs=[job])
        submit_hpc(subset,out,commands,state,bsub,timeout)
        root["new_id"]=state["jobs"][name]["job_id"]; save_state(out,state)
    context["status"]="releasing"; save_state(out,state)
    for name in context["affected"]:
        if name in context["released"]: continue
        condition=dependency(by_id[name],state)
        change_pending_dependency(state,out,context,name,condition,"release",bjobs,bmod,timeout)
    context.update(status="complete",finished_at=utcnow()); save_state(out,state)


def wait_for_plan(plan,out,state,bjobs="bjobs",bwait="bwait",timeout=None):
    started=time.monotonic()
    while True:
        receipts=[receipt_for_submission(run_receipt(out,j["id"],plan["identity"]),state["jobs"][j["id"]]) for j in plan["jobs"]]
        if all(r is not None and r.get("status")=="succeeded" for r in receipts): return state
        statuses=scheduler_states(state,bjobs); outstanding=[]
        for job in plan["jobs"]:
            name=job["id"]; row=state["jobs"][name]
            receipt=receipt_for_submission(run_receipt(out,name,plan["identity"]),row)
            scheduler=statuses.get(row["job_id"])
            row["scheduler_status"]=scheduler or "not_observed"
            if receipt and receipt.get("status")=="succeeded": continue
            if scheduler=="EXIT" or receipt and receipt.get("status")=="failed":
                save_state(out,state); raise RuntimeError(f"Job {name} failed; done dependencies remain unsatisfied")
            if scheduler=="DONE":
                raise RuntimeError(f"Job {name} reports DONE without a successful worker receipt")
            if scheduler not in ACTIVE:
                raise RuntimeError(f"Job {name} state unavailable; no completion inferred")
            outstanding.append(row["job_id"])
        save_state(out,state)
        if not outstanding: return state
        if timeout is not None and time.monotonic()-started>=timeout: raise RuntimeError("Wait timeout; jobs remain submitted")
        # One-minute scheduler waits are interruptible; reevaluate failures after
        # each wait rather than waiting forever for a failed dependency's child.
        result=subprocess.run([bwait,"-w",f"ended({outstanding[0]})","-t","1"],capture_output=True,text=True,timeout=65)
        if result.returncode not in (0,1): raise RuntimeError("bwait failed: "+result.stderr.strip())


def dispatch_frozen(plan,out,write_only=False,jobs_n=1,wait=False,bsub="bsub",bjobs="bjobs",bwait="bwait",bmod="bmod",retry_failed=False,submit_timeout=30,wait_timeout=None):
    validate_frozen(plan); verify_inputs(plan); out=Path(out).absolute()
    with locked(out):
        commands=prepare(plan,out)
        if write_only:
            return dict(status="planned",mode=plan["mode"],plan_identity=plan["identity"],jobs=len(plan["jobs"]),out=str(out))
        state=read_state(out,plan)
        if plan["mode"]=="local": return execute_local(plan,out,commands,state,jobs_n)
        if state.get("retry",{}).get("status") not in {None,"complete"} and not retry_failed:
            raise RuntimeError("An explicit retry is incomplete; resume it with --retry-failed; no IDs were discarded")
        if retry_failed: retry_failed_hpc(plan,out,commands,state,bsub,bjobs,bmod,submit_timeout)
        state=submit_hpc(plan,out,commands,state,bsub,submit_timeout)
    if wait: return wait_for_plan(plan,out,state,bjobs,bwait,wait_timeout)
    return state


def run_plan(plan,out,mode="hpc",write_only=False,jobs_n=1,wait=False,**kwargs):
    return dispatch_frozen(normalize_plan(plan,mode),out,write_only,jobs_n,wait,**kwargs)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("command",nargs="?",choices=("run","dispatch","worker"),default="run")
    p.add_argument("--plan"); p.add_argument("--frozen-plan"); p.add_argument("--job-id")
    p.add_argument("--mode",choices=("local","hpc"),default="hpc")
    p.add_argument("--out",required=True); p.add_argument("--write-cmd",action="store_true")
    p.add_argument("--jobs-N",type=int,default=1); p.add_argument("--wait",action="store_true")
    p.add_argument("--bsub",default="bsub"); p.add_argument("--bjobs",default="bjobs"); p.add_argument("--bwait",default="bwait")
    p.add_argument("--bmod",default="bmod"); p.add_argument("--retry-failed",action="store_true")
    p.add_argument("--submit-timeout",type=float,default=30); p.add_argument("--wait-timeout",type=float)
    return p


def main(argv=None):
    a=parser().parse_args(argv)
    try:
        if a.command=="worker":
            if not a.frozen_plan or not a.job_id: raise ValueError("worker requires --frozen-plan and --job-id")
            return worker(a.frozen_plan,a.job_id,a.out)
        kwargs=dict(write_only=a.write_cmd,jobs_n=a.jobs_N,wait=a.wait,bsub=a.bsub,bjobs=a.bjobs,bwait=a.bwait,bmod=a.bmod,retry_failed=a.retry_failed,
                    submit_timeout=float(positive(a.submit_timeout,"submit_timeout")),wait_timeout=a.wait_timeout)
        if a.command=="dispatch":
            if not a.frozen_plan: raise ValueError("dispatch requires --frozen-plan")
            plan=json.loads(Path(a.frozen_plan).read_text()); result=dispatch_frozen(plan,a.out,**kwargs)
        else:
            if not a.plan: raise ValueError("run requires --plan")
            result=run_plan(json.loads(Path(a.plan).read_text()),a.out,mode=a.mode,**kwargs)
        print(json.dumps(dict(status=result.get("status","submitted"),out=str(Path(a.out).absolute()),
                              plan_identity=result.get("plan_identity")),sort_keys=True))
        return 0
    except (OSError,ValueError,RuntimeError,subprocess.TimeoutExpired) as exc:
        print("batch-runner: "+str(exc),file=sys.stderr); return 2


if __name__=="__main__": raise SystemExit(main())
