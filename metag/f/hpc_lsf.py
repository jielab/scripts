"""Generate single-node LSF jobs and submit only when explicitly called.

doctor() returns its report; the CLI owns saving doctor.json.  Resource output
is retained as evidence, never treated as proof that a future job will fit.
"""
from __future__ import annotations

from decimal import Decimal, InvalidOperation, ROUND_CEILING
import hashlib
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import uuid

from hpc_common import PATH_KEYS, atomic_json, utcnow


_DEFAULTS = {
    "setup": ("ser", 8, 16, 24),
    "prepare": ("ser", 8, 32, 72),
    "build": ("smp", 32, 4096, 168),
    "pilot": ("smp", 32, 4096, 168),
    "classify": ("smp", 32, 4096, 168),
    "convert": ("short", 40, None, 72),
    "compare_summary": ("ser", 2, 16, 24),
}
_FALLBACK = {"setup": "prepare", "pilot": "classify"}
_FIELDS = (
    "RESOURCE_RESERVE_PER_TASK", "LSF_UNIT_FOR_LIMITS", "MAX_JOB_ARRAY_SIZE",
    "MAX_PEND_JOBS", "MAX_PEND_SLOTS", "LSF_LINUX_CGROUP_ACCT",
    "LSB_CGROUP_MEM_INCLUDE_CACHE", "LSB_JOB_MEMLIMIT", "LSB_RESOURCE_ENFORCE",
    "LSF_REPLACE_PIM_WITH_LINUX_CGROUP", "EGO_PIM_SWAP_REPORT", "MEMLIMIT",
    "RUNLIMIT", "SWAPLIMIT", "RES_REQ", "RESRSV_LIMIT", "UJOB_LIMIT",
    "QJOB_LIMIT", "HJOB_LIMIT", "ABS_RUNLIMIT",
)
_INDEX = "__LSB_JOBINDEX__"


def _text(value, name):
    value = str(value)
    if any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError(f"Control character in {name}")
    return value


def _number(value, name):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive number")
    try:
        number = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"{name} must be a positive number") from exc
    if not number.is_finite() or number <= 0:
        raise ValueError(f"{name} must be a positive finite number")
    return number


def _settings(cfg, kind):
    if kind not in _DEFAULTS:
        raise ValueError(f"Unknown LSF job kind: {kind}")
    values = {}
    for field, default in zip(("queue", "cores", "memory_gib", "hours"), _DEFAULTS[kind]):
        key = f"{kind}_{field}"
        fallback = f"{_FALLBACK.get(kind, kind)}_{field}"
        values[field] = cfg.get(key, cfg.get(fallback, default))
    queue = _text(values["queue"], "queue")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.+-]*", queue):
        raise ValueError(f"Invalid queue name: {queue!r}")
    cores = _number(values["cores"], "cores")
    if cores != cores.to_integral_value():
        raise ValueError("cores must be a positive integer")
    memory = _number(values["memory_gib"] if values["memory_gib"] is not None else os.environ["MEMORY_DEFAULT_GIB"], "memory_gib")
    minutes = int((_number(values["hours"], "hours") * 60).to_integral_value(rounding=ROUND_CEILING))
    return queue, int(cores), memory, f"{minutes // 60}:{minutes % 60:02d}"


def _memory(cfg, memory, cores):
    scope = cfg.get("lsf_memory_scope", "host")
    confirmation = cfg.get("lsf_resource_reserve_per_task")
    if scope == "host":
        amount = int(memory.to_integral_value(rounding=ROUND_CEILING))
        return f"{amount}GB/host"
    if scope == "site_host" and confirmation is False:
        amount = int(memory.to_integral_value(rounding=ROUND_CEILING))
        return f"{amount}GB"
    if scope == "site_task" and confirmation is True:
        amount = int((memory / cores).to_integral_value(rounding=ROUND_CEILING))
        return f"{amount}GB"
    raise ValueError("lsf_memory_scope must be host, or site_host/site_task with "
                     "lsf_resource_reserve_per_task explicitly false/true respectively")


def _path(value, name):
    return Path(_text(os.fspath(value), name)).expanduser().absolute()


def _logs(cfg):
    logs = _path(cfg["work_dir"], "work_dir") / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    return logs


def _job_label(name, array):
    name = _text(name, "job_name")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name):
        raise ValueError("job_name must contain only letters, digits, _, . and -")
    if array is None:
        return name
    value = _text(array, "array")
    if value.count("%") > 1:
        raise ValueError("Invalid array concurrency")
    indices, sep, limit = value.partition("%")
    if sep and (not limit.isdigit() or int(limit) < 1):
        raise ValueError("Array concurrency must be positive")
    for item in indices.split(","):
        match = re.fullmatch(r"([1-9]\d*)(?:-([1-9]\d*)(?::([1-9]\d*))?)?", item)
        if not match or (match[2] is not None and int(match[2]) < int(match[1])):
            raise ValueError(f"Invalid array index range: {item!r}")
    return f"{name}[{indices}]" + (f"%{limit}" if sep else "")


def _shell_arg(value):
    value = _text(os.fspath(value) if isinstance(value, os.PathLike) else value, "command argument")
    if _INDEX not in value:
        return shlex.quote(value)
    parts = value.split(_INDEX)
    return '"${LSB_JOBINDEX:?}"'.join(shlex.quote(part) if part else "" for part in parts)


def write_cmd(cfg, path, kind, command, job_name, array=None):
    """Write a reviewable script.  This function never invokes bsub."""
    if isinstance(command, (str, bytes)) or not command:
        raise ValueError("command must be a nonempty list of arguments")
    queue, cores, memory, walltime = _settings(cfg, kind)
    resource = _memory(cfg, memory, cores)
    label = _job_label(job_name, array)
    rendered = " ".join(_shell_arg(arg) for arg in command)
    destination = _path(path, "command path")
    code_root = _path(cfg.get("code_root", Path(__file__).resolve().parent.parent), "code_root")
    if not code_root.is_dir():
        raise ValueError(f"Code directory does not exist: {code_root}")
    logs = _logs(cfg)
    stem = destination.stem
    text = "\n".join([
        "#!/usr/bin/env bash",
        f'#BSUB -J "{label}"',
        f"#BSUB -q {queue}",
        f"#BSUB -n {cores}",
        f"#BSUB -W {walltime}",
        f'#BSUB -R "span[hosts=1] rusage[mem={resource}]"',
        f"#BSUB -o {shlex.quote(str(logs / (stem + '.%J_%I.out')))}",
        f"#BSUB -e {shlex.quote(str(logs / (stem + '.%J_%I.err')))}",
        "",
        "set -euo pipefail",
        'export submit="${submit:-hpc}" WORKFLOW_WORKER=1',
        "export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1",
        "export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1",
        f"cd -- {shlex.quote(str(code_root))}",
        f"exec {rendered}",
        "",
    ])
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix="." + destination.name + ".", dir=destination.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, 0o750)
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return destination


def _decoded(value):
    return value.decode("utf-8", "replace") if isinstance(value, bytes) else (value or "")


def _probe(argv, timeout):
    result = {"argv": argv, "stdout": "", "stderr": "", "returncode": None}
    try:
        process = subprocess.run(argv, capture_output=True, text=True, encoding="utf-8",
                                 errors="replace", timeout=timeout,
                                 env={**os.environ, "LC_ALL": "C"})
        result.update(stdout=process.stdout, stderr=process.stderr, returncode=process.returncode,
                      status="ok" if process.returncode == 0 else "failed")
    except FileNotFoundError as exc:
        result.update(status="missing", stderr=str(exc))
    except subprocess.TimeoutExpired as exc:
        result.update(status="timeout", stdout=_decoded(exc.stdout), stderr=_decoded(exc.stderr))
    except OSError as exc:
        result.update(status="error", stderr=str(exc))
    return result


def _field_evidence(text, source):
    lines = text.splitlines()
    evidence = []
    for index, line in enumerate(lines):
        for field in _FIELDS:
            if re.search(r"\b" + field + r"\b", line):
                match = re.match(r"\s*" + field + r"\s*[:=]\s*(.*?)\s*$", line)
                evidence.append({"field": field, "source": source,
                                 "value": match[1] if match else None,
                                 "context": "\n".join(lines[index:index + 3])})
    return evidence


def _site_config_paths():
    paths = {Path("/etc/lsf.conf")}
    for key in ("LSF_ENVDIR", "LSF_CONFDIR"):
        if os.environ.get(key):
            paths.add(Path(os.environ[key]) / "lsf.conf")
    return sorted(paths)


def doctor(cfg):
    """Inspect local LSF commands and files only; return report without writing it."""
    timeout = min(60.0, float(_number(cfg.get("doctor_timeout_seconds", 15), "doctor_timeout_seconds")))
    queues = {_settings(cfg, kind)[0] for kind in ("setup", "prepare", "build", "classify")}
    if cfg.get("database_layout") == "partitioned":
        queues |= {str(value["queue"]) for value in cfg.get("category_resources", {}).values() if "queue" in value}
    else:
        queues |= {_settings(cfg, "convert")[0]} | {str(cfg.get(f"{mode}_queue", "ser")) for mode in ("kunpeng_full", "kunpeng_memory")}
    queues = sorted(queues)
    commands = [["lsid"], ["bsub", "-V"], ["bparams", "-a"]]
    commands += [["bqueues", "-l", queue] for queue in queues]
    commands += [["bhosts", "hg_" + queue] for queue in queues]
    commands += [["lshosts", "-w"], ["busers"], ["blimits"]]
    paths = []
    df_targets = set()
    for key in (*PATH_KEYS, "root", "code_root", "db_dir"):
        if not cfg.get(key):
            continue
        path = _path(cfg[key], key)
        ancestor = path
        while not ancestor.exists() and ancestor != ancestor.parent:
            ancestor = ancestor.parent
        paths.append({"key": key, "path": str(path), "exists": path.exists(),
                      "df_path": str(ancestor)})
        df_targets.add(str(ancestor))
    commands += [["df", "-h", "--", path] for path in sorted(df_targets)]
    checks, evidence = [], []
    print("LSF doctor: collecting local scheduler and filesystem evidence.", flush=True)
    for command in commands:
        check = _probe(command, timeout)
        checks.append(check)
        evidence.extend(_field_evidence(check["stdout"] + "\n" + check["stderr"], shlex.join(command)))
        print(f"  {shlex.join(command)}: {check['status']}", flush=True)
    # Only selected scheduling keys are copied, not arbitrary site configuration.
    configuration_files = []
    for path in _site_config_paths():
        item = {"path": str(path)}
        try:
            with path.open(encoding="utf-8", errors="replace") as handle:
                selected = [line.rstrip("\r\n") for line in handle
                            if any(re.match(r"\s*" + field + r"\s*=", line) for field in _FIELDS)]
            item.update(status="ok", selected_lines=selected)
            evidence.extend(_field_evidence("\n".join(selected), str(path)))
        except FileNotFoundError:
            item["status"] = "missing"
        except OSError as exc:
            item.update(status="error", error=str(exc))
        configuration_files.append(item)
    observed = {item["value"].upper() for item in evidence
                if item["field"] == "RESOURCE_RESERVE_PER_TASK" and item["value"] is not None}
    per_task = next(iter(observed)) == "Y" if len(observed) == 1 and observed <= {"Y", "N"} else None
    warnings = [
        "Queue output is evidence, not confirmation of available memory, permission, quota or dispatch time.",
        "Review MEMLIMIT, RUNLIMIT, RES_REQ, SWAPLIMIT, cgroup accounting and job/array limits before submission.",
        "df shows filesystem capacity; project/user quotas and scratch expiry require separate verification.",
    ]
    scope = cfg.get("lsf_memory_scope", "host")
    if scope == "host":
        warnings.append("Explicit /host requires LSF 10.1 Fix Pack 9 or newer; a bsub rejection is never retried with weaker resource syntax.")
    else:
        try:
            _memory(cfg, Decimal(1), 1)
        except ValueError as exc:
            warnings.append(str(exc))
        confirmation = cfg.get("lsf_resource_reserve_per_task")
        if per_task is not None and confirmation is not per_task:
            warnings.append("Observed RESOURCE_RESERVE_PER_TASK conflicts with the configured site confirmation.")
    if per_task is None:
        warnings.append("RESOURCE_RESERVE_PER_TASK was not unambiguously observed; its site default remains unverified.")
    for warning in warnings:
        print(f"  NOTE: {warning}", flush=True)
    return {"schema": 1, "created_at": utcnow(), "commands": checks, "paths": paths,
            "configuration_files": configuration_files, "field_evidence": evidence,
            "observed_resource_reserve_per_task": per_task, "lsf_memory_scope": scope,
            "resources_verified": False, "warnings": warnings}


class SubmissionError(RuntimeError):
    """The result/receipt can include a job ID even if reporting failed afterwards."""
    def __init__(self, message, result):
        super().__init__(message)
        self.result = result


def local_submit(cfg, path, expected_sha256=None):
    raw = Path(path).read_bytes()
    command_sha256 = hashlib.sha256(raw).hexdigest()
    if expected_sha256 is not None and command_sha256 != expected_sha256:
        result = {'schema': 1, 'status': 'failed', 'job_id': None,
                  'command_path': str(path), 'command_sha256': command_sha256,
                  'error': 'Generated command changed before local execution'}
        raise SubmissionError(result['error'], result)
    text = raw.decode('utf-8')
    match = re.search(r'^#BSUB -J .*\[([0-9,:-]+)\](?:%([1-9][0-9]*))?', text, re.M)
    indices = None
    array_limit = None
    if match:
        indices = []
        array_limit = int(match[2]) if match[2] else None
        for token in match[1].split(','):
            span, _, step = token.partition(':')
            bounds = list(map(int, span.split('-')))
            indices.extend(range(bounds[0], bounds[-1]+1, int(step or 1)))
    # Full databases can require terabytes. Fail before work if this host
    # cannot supply the same algorithm's requested memory; never downsample.
    required = re.search(r'rusage\[mem=([0-9.]+)([KMGT]?B)(/host)?\]', text)
    if required is None:
        raise RuntimeError('Cannot establish the local memory request from the generated .cmd')
    units = {'B': 1, 'KB': 1024, 'MB': 1024**2, 'GB': 1024**3, 'TB': 1024**4}
    job_bytes = float(required[1]) * units[required[2]]
    if not required[3] and cfg.get('lsf_memory_scope') == 'site_task':
        cores = re.search(r'^#BSUB -n ([1-9][0-9]*)$', text, re.M)
        if cores is None:
            raise RuntimeError('Per-task memory needs an explicit CPU request')
        job_bytes *= int(cores[1])
    jobs = int(os.environ.get('jobs_N', '8'))
    if not 1 <= jobs <= 8:
        raise ValueError('jobs_N must be 1..8')
    concurrent = min(jobs, len(indices), array_limit or jobs) if indices is not None else 1
    available = os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE')
    if job_bytes * concurrent > available:
        raise RuntimeError('Requested metag memory exceeds local RAM; use submit=hpc or an appropriate explicit config')
    results = _run_local([path], concurrent, indices)
    result = {'schema': 1, 'job_id': 'local-' + uuid.uuid4().hex[:12],
              'status': 'submitted' if all(r['returncode'] == 0 for r in results) else 'failed',
              'command_path': str(path), 'command_sha256': command_sha256,
              'results': results, 'submitted_at': utcnow()}
    receipt = _logs(cfg) / (Path(path).stem + '.local.' + uuid.uuid4().hex + '.json')
    atomic_json(receipt, result)
    if result['status'] != 'submitted':
        raise SubmissionError('Local metag stage failed: ' + str(receipt), result)
    return result


def _run_local(paths, jobs=1, indices=None):
    """Explicit developer-only local runner; no external repository dependency."""
    from concurrent.futures import ThreadPoolExecutor
    def one(item):
        path, index = item
        env = {**os.environ, "submit": "local", "WORKFLOW_WORKER": "1"}
        if index is not None:
            env["LSB_JOBINDEX"] = str(index)
        suffix = f".{index}" if index is not None else ""
        with Path(str(path) + suffix + ".out").open("a") as out, Path(str(path) + suffix + ".err").open("a") as err:
            process = subprocess.run(["bash", str(path)], stdout=out, stderr=err, env=env)
        return {"cmd": str(path), "index": index, "returncode": process.returncode}
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(one, [(path, index) for path in paths for index in (indices or [None])]))


def submit_cmd(cfg, path, dependency=None, expected_sha256=None):
    """Submit one script via stdin, record the real ID, and never auto-retry."""
    path = _path(path, "command path")
    if os.environ.get('submit') == 'local':
        return local_submit(cfg, path, expected_sha256=expected_sha256)
    logs = _logs(cfg)
    argv = ["bsub", "-o", str(logs / f"{path.stem}.%J.%I.out"),
            "-e", str(logs / f"{path.stem}.%J.%I.err")]
    if dependency is not None:
        if isinstance(dependency, int) and not isinstance(dependency, bool) and dependency > 0:
            dependency = f"done({dependency})"
        dependency = _text(dependency, "dependency")
        if not dependency.strip():
            raise ValueError("dependency must be nonempty")
        argv += ["-w", dependency]
    timeout = min(60.0, float(_number(cfg.get("submit_timeout_seconds", 30), "submit_timeout_seconds")))
    result = {"schema": 1, "submitted_at": utcnow(), "command_path": str(path),
              "dependency": dependency, "argv": argv, "job_id": None,
              "returncode": None, "stdout": "", "stderr": "", "status": "unknown"}
    # Keep this descriptor open for both the recorded hash and bsub's stdin.
    with path.open("rb") as handle:
        result["command_sha256"] = hashlib.sha256(handle.read()).hexdigest()
        if expected_sha256 is not None and result["command_sha256"] != expected_sha256:
            result.update(status="failed", error="Generated command changed before submission; not sent to bsub")
            raise SubmissionError(result["error"], result)
        handle.seek(0)
        try:
            process = subprocess.run(argv, stdin=handle, capture_output=True, text=True,
                                     encoding="utf-8", errors="replace", timeout=timeout,
                                     env={**os.environ, "LC_ALL": "C"})
            result.update(stdout=process.stdout, stderr=process.stderr, returncode=process.returncode)
        except subprocess.TimeoutExpired as exc:
            result.update(stdout=_decoded(exc.stdout), stderr=_decoded(exc.stderr),
                          error="bsub timed out; scheduler acceptance is unknown")
        except OSError as exc:
            result.update(error=str(exc), status="failed")
    ids = set(re.findall(r"\bJob\s*<([1-9]\d*)>", result["stdout"] + "\n" + result["stderr"]))
    if len(ids) == 1:
        result["job_id"] = int(next(iter(ids)))
    if result["returncode"] == 0 and result["job_id"] is not None:
        result["status"] = "submitted"
    elif result["returncode"] not in (None, 0):
        result["status"] = "failed"
    receipt = logs / f"{path.stem}.submission.{uuid.uuid4().hex}.json"
    result["receipt_path"] = str(receipt)
    try:
        atomic_json(receipt, result)
    except OSError as exc:
        raise SubmissionError(f"bsub returned job_id={result['job_id']}; submission receipt could not be saved: {exc}. "
                              "Inspect scheduler state before submitting again.", result) from exc
    if result["status"] != "submitted":
        message = (f"bsub submission status={result['status']}, returncode={result['returncode']}, "
                   f"job_id={result['job_id']}; receipt={receipt}\n"
                   f"{result.get('error', '')}\n{result['stdout']}{result['stderr']}")
        raise SubmissionError(message, result)
    return result
