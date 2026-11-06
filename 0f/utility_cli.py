"""Shared local/HPC options for standalone maintenance commands.

Paths come only from 0os.sh. HPC work uses the existing .cmd planner and LSF
receipt handling; local execution preserves foreground progress and exit codes.
"""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys


def add_options(parser):
    group = parser.add_argument_group("local / HPC")
    group.add_argument("--submit", choices=["local", "hpc", "auto"],
                       default=os.environ.get("submit", os.environ.get("SUBMIT", "hpc")))
    for name in ("dat", "big", "tmp"):
        group.add_argument("--dir-" + name, help="Override the corresponding 0os.sh root")
    group.add_argument("--write-cmd", action="store_true", help="Write .cmd/plan only; do not execute or submit")
    group.add_argument("--dir-jobs", help="Plan/log parent directory")
    group.add_argument("--queue", default=os.environ.get("HPC_CPU_QUEUE", "33"))
    group.add_argument("--slots", type=int, default=1)
    group.add_argument("--memory-gb", type=float, default=4)
    group.add_argument("--walltime", default="24:00")
    group.add_argument("--wait", choices=["TRUE", "FALSE"], default=os.environ.get("SUBMIT_WAIT", "FALSE"))
    group.add_argument("--execute", action="store_true", help=argparse.SUPPRESS)


def initialize(options):
    environment = dict(os.environ, submit=options.submit, SUBMIT=options.submit)
    for suffix in ("dat", "big", "tmp"):
        value = getattr(options, "dir_" + suffix)
        if value is not None:
            if not Path(value).is_absolute():
                raise ValueError("--dir-" + suffix + " requires an absolute path")
            environment["DIR_" + suffix.upper()] = value.rstrip("/") or "/"
            if suffix == "big":
                environment["DIR_DATA"] = environment["DIR_BIG"]
            elif suffix == "tmp":
                environment["DIR_SCRATCH"] = environment["DIR_TMP"]
    command = ('source "$1" || exit; '
               'for name in "${OS_PATH_KEYS[@]}" submit SUBMIT SHARED_PYTHON HPC_CPU_QUEUE; do '
               'printf "%s\\0%s\\0" "$name" "${!name-}"; done')
    raw = subprocess.check_output(["bash", "-c", command, "utility-paths",
        str(Path(__file__).with_name("0os.sh"))], env=environment)
    fields = raw.decode().split("\0")[:-1]
    os.environ.update(dict(zip(fields[::2], fields[1::2])))
    options.submit = os.environ["submit"]
    if options.execute and options.submit == "hpc" and not os.environ.get("LSB_JOBID"):
        raise ValueError("HPC execution requires an LSF allocation; submit the .cmd with bsub")
    return os.environ


def dispatch(options, name, argv, *, preview=False, extra_env=None, job_name=None):
    if options.execute or (preview and not options.write_cmd):
        return False
    if options.submit == "local" and not options.write_cmd:
        return False
    match = re.fullmatch(r"([0-9]+):([0-5][0-9])", options.walltime)
    if not match:
        raise ValueError("--walltime must be HHH:MM")
    hours = int(match[1]) + int(match[2]) / 60
    from workflow_jobs import Job, write_plan, run_plan
    python = os.environ.get("UTILITY_PYTHON") or os.environ["SHARED_PYTHON"]
    os.environ["WORKFLOW_PYTHON"] = python
    env = {k: v for k, v in os.environ.items() if k.startswith("GITHUB_COPY_") or k == "UTILITY_PYTHON"}
    env.update(extra_env or {})
    env.update(submit=options.submit, SUBMIT=options.submit, WORKFLOW_WORKER="1")
    root = options.dir_jobs or str(Path(os.environ["ANALYSIS_ROOT"]) / name / "jobs")
    plan = write_plan([Job(name, [*argv, "--submit", options.submit, "--execute"],
        resources={"queue": options.queue, "cores": options.slots,
                   "memory_gib": options.memory_gb, "hours": hours},
        env=env, cwd=os.getcwd(), job_name=job_name)], root, mode=options.submit, name=name)
    if not options.write_cmd:
        run_plan(plan, submit=options.submit, jobs=1, wait=options.wait == "TRUE")
    else:
        print("Command file only: " + str(plan.parent / "commands" / (name + ".cmd")), flush=True)
    return True
