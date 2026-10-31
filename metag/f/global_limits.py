"""Run the entire global pipeline in its own verified cgroup v2 memory scope."""
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

from global_resources import GIB, memory_status

MARKER = "METAG_GLOBAL_MEMORY_SCOPE"


def stop_scope(unit, process):
    """Stop only the transient scope created by this launcher."""
    try:
        subprocess.run(["systemctl", "--user", "stop", "--no-block", unit],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=5)
    finally:
        if process.poll() is None:
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                subprocess.run(["systemctl", "--user", "kill", "--kill-whom=all",
                                "--signal=SIGKILL", unit],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=5)
                process.kill()
                process.wait(timeout=5)


def current_group():
    for line in Path("/proc/self/cgroup").read_text().splitlines():
        if line.startswith("0::/"):
            relative = line[4:]
            if ".." in Path(relative).parts:
                break
            return Path("/sys/fs/cgroup") / relative
    raise RuntimeError("A cgroup v2 memory controller is required; refusing an unlimited run")


def verify_scope(unit, limit):
    group = current_group()
    if group.name != unit or not unit.startswith("metag-global-"):
        raise RuntimeError("Global memory scope does not match the current cgroup")
    observed = int((group / "memory.max").read_text().strip())
    swap = (group / "memory.swap.max").read_text().strip()
    oom_group = (group / "memory.oom.group").read_text().strip()
    if observed != limit or swap != "0" or oom_group != "1":
        raise RuntimeError("Global memory/swap/OOM group limits were not applied; refusing to run")
    # Prefer this job if the machine encounters unrelated global memory pressure.
    Path("/proc/self/oom_score_adj").write_text("500")
    return {"cgroup": str(group), "memory_max_bytes": observed,
            "swap_max_bytes": 0, "oom_group": True}


def ensure_memory_scope(args, argv):
    """Return a child exit status, or None once the current scope is verified.

    All native children inherit this cgroup, including new process sessions.
    No process outside this scope is terminated or reconfigured.
    """
    limit = math.ceil(args.memory_limit_gib * GIB)
    if args.command in ("all", "build") and (
            args.db_gib + args.memory_reserve_gib + 1 > args.memory_limit_gib):
        raise RuntimeError("--memory-limit-gib must exceed --db-gib + --memory-reserve-gib "
                           "by at least 1 GiB for the Python controller")
    if os.environ.get(MARKER):
        args.resource_limits = verify_scope(os.environ[MARKER], limit)
        print(f"Memory scope verified: {args.memory_limit_gib:g} GiB total for this job and "
              "all children; swap disabled; OOM confined to this job", flush=True)
        return None
    state = memory_status()
    # Leave headroom outside the job as well as the internal build/query reserve.
    if state["effective_available_bytes"] < limit + 8 * GIB:
        raise RuntimeError(f"Insufficient headroom to start a {args.memory_limit_gib:g} GiB job: "
                           f"{state['effective_available_bytes']/GIB:.1f} GiB available; "
                           "8 GiB must remain outside the job. Other tasks are left running.")
    runner = shutil.which("systemd-run")
    if not runner:
        raise RuntimeError("systemd-run --user is required for a process-tree memory limit; "
                           "refusing to run without it")
    unit = "metag-global-" + uuid.uuid4().hex + ".scope"
    command = [runner, "--user", "--scope", "--quiet", "--unit", unit,
               "-p", f"MemoryMax={limit}", "-p", f"MemoryHigh={int(limit * 0.9)}",
               "-p", "MemorySwapMax=0", "-p", "OOMPolicy=kill",
               sys.executable, str(Path(__file__).with_name("global_pipeline.py")), *argv]
    process = subprocess.Popen(command, env={**os.environ, MARKER: unit})
    try:
        while process.poll() is None:
            if memory_status()["effective_available_bytes"] < 8 * GIB:
                raise RuntimeError("Host available memory fell below 8 GiB; stopping only this "
                                   "global job to leave other programs running. Retry later.")
            time.sleep(1)
        return process.returncode if process.returncode >= 0 else 128 - process.returncode
    finally:
        if process.poll() is None:
            stop_scope(unit, process)
