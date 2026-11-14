#!/usr/bin/env python3
"""Plan one existing sequential project entry without loading its data/runtime."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import sys
from workflow_jobs import Job, resolve_mode, run_plan, write_plan


def scientific_budget(project, arguments):
    """Honor exposed worker counts without rewriting scientific arguments."""
    cores = int(os.getenv("DOSE_THREADS", "4")) if project == "dose" else 4
    memory = float(os.getenv("DOSE_MEMORY_LIMIT_GB", "16")) if project == "dose" else 0
    if project == "dose" and os.getenv("DOSE_SHAP_JOBS"):
        cores = max(cores, int(os.environ["DOSE_SHAP_JOBS"]))
    if project == "pfi":
        cores = max(cores, *(int(os.getenv(k, "4")) for k in
            ("N_JOBS", "PMC_PARSE_WORKERS", "PFI_FEATURE_WORKERS", "PFI_DUCKDB_THREADS")),
            int(os.getenv("PFI_BUILD_JOBS", "2")) *
            int(os.getenv("PFI_BUILD_DUCKDB_THREADS", os.getenv("PFI_DUCKDB_THREADS", "4"))))
    if project == "refGen" and arguments[:1] == ["download_archaic"]:
        cores = 8
    for i, token in enumerate(arguments):
        if token == "--":
            break
        key, sep, value = token.partition("=")
        if key not in {"--jobs", "--threads", "--bgzip-threads", "--memory-limit-gb"}:
            continue
        if not sep:
            if i + 1 >= len(arguments):
                raise ValueError(f"{key} requires a value")
            value = arguments[i + 1]
        if key == "--memory-limit-gb":
            memory = max(memory, float(value))
        else:
            cores = max(cores, int(value))
    return cores, memory


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--entry", required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    a = parser.parse_args(argv)
    arguments = a.arguments[1:] if a.arguments[:1] == ["--"] else a.arguments
    resource = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    resource.add_argument("--plan-only", "--write-cmd", action="store_true")
    resource.add_argument("--job-cores", type=int)
    resource.add_argument("--job-memory-gb", type=float)
    resource.add_argument("--job-hours", type=float, default=120)
    resource.add_argument("--queue")
    resource.add_argument("--dir-jobs")
    if any(token == "--job-gpus" or token.startswith("--job-gpus=") for token in arguments):
        raise ValueError("This sequential project entry is CPU-only; no --job-gpus option is provided")
    opts, scientific = resource.parse_known_args(arguments)
    required_cores, required_memory = scientific_budget(a.project, scientific)
    if opts.job_cores is None:
        opts.job_cores = required_cores
    elif opts.job_cores < required_cores:
        raise ValueError(f"--job-cores must cover the scientific worker budget ({required_cores})")
    if opts.job_memory_gb is None:
        opts.job_memory_gb = max(float(os.getenv("MEMORY_DEFAULT_GIB", "32")), required_memory)
    elif opts.job_memory_gb < required_memory:
        raise ValueError(f"--job-memory-gb must cover --memory-limit-gb ({required_memory:g})")
    mode = resolve_mode()
    root = Path(os.environ["ANALYSIS_ROOT"]) / a.project
    if a.project == "pfi" and os.getenv("PFI_ANALYSIS_ROOT"):
        root = Path(os.environ["PFI_ANALYSIS_ROOT"])
    elif a.project == "dose":
        root = Path(os.getenv("DOSE_OUTDIR") or os.getenv("DOSE_PROJECT") or root)
    # Existing custom output roots remain authoritative, not silently ignored.
    aliases = {"--dir-output", "--output-dir", "--dir-analysis-root", "--analysis-root",
               "--dir-out", "--out-dir", "--outdir", "--dir-analysis", "--analysis-dir"}
    for i, token in enumerate(scientific):
        if token == "--":
            break
        key, sep, value = token.partition("=")
        if key in aliases:
            if not sep:
                if i + 1 >= len(scientific):
                    raise ValueError(f"{key} requires a path")
                value = scientific[i + 1]
            root = Path(value).expanduser().resolve()
    queue = opts.queue or os.getenv("HPC_CPU_QUEUE", "33")
    env = {"PROJECT_JOB_ACTIVE": "1", "WORKFLOW_RESUME_POLICY": "new_plan",
           "WORKFLOW_DEVICE": "cpu"}
    # Freeze only workflow/project settings, never the full environment.
    prefix = a.project.upper().replace("-", "_") + "_"
    env.update({k: v for k, v in os.environ.items() if k.startswith(prefix)})
    if a.project == "pfi":
        env.update({k: os.environ[k] for k in ("N_JOBS", "PMC_PARSE_WORKERS") if k in os.environ})
    job = Job(a.project, ["bash", str(Path(a.entry).resolve()), *scientific],
              resources={"queue": queue, "cores": opts.job_cores, "memory_gib": opts.job_memory_gb,
                         "hours": opts.job_hours, "gpus": 0}, env=env, cwd=os.getcwd())
    plan = write_plan([job], opts.dir_jobs or root / "logs", mode=mode, name=a.project)
    result = run_plan(plan, submit=mode, jobs=1, wait=os.getenv("SUBMIT_WAIT") == "TRUE",
                      dry_run=opts.plan_only)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, RuntimeError) as exc:
        sys.exit(f"ERROR: {exc}")
