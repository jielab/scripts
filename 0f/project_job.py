#!/usr/bin/env python3
"""Plan one existing sequential project entry without loading its data/runtime."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import sys
from workflow_cpu import core_budget
from workflow_jobs import Job, resolve_mode
from workflow_entry import Workflow, execute_workflow


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


def build_workflow(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True)
    parser.add_argument("--entry", required=True)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    a = parser.parse_args(argv)
    arguments = a.arguments[1:] if a.arguments[:1] == ["--"] else a.arguments
    resource = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    resource.add_argument("--plan-only", "--write-cmd", action="store_true")
    resource.add_argument("--cores", type=int)
    resource.add_argument("--job-cores", type=int)
    resource.add_argument("--job-memory-gb", type=float)
    resource.add_argument("--job-hours", type=float, default=float(os.getenv("JOB_HOURS", "120")))
    resource.add_argument("--queue")
    resource.add_argument("--dir-jobs")
    resource.add_argument("--job-gpus", type=int, default=0)
    opts, scientific = resource.parse_known_args(arguments)
    entry_name = Path(a.entry).name
    gpu_stage = a.project == "pfi" and (
        entry_name in {"09b_train_pubmedbert_pfi_model.sh", "09c_predict_pubmedbert_pfi.sh"}
        or (entry_name == "pfi.sh" and scientific[:1] in (["train"], ["predict"])))
    explicit_gpu = any(t == "--job-gpus" or t.startswith("--job-gpus=") for t in arguments)
    if explicit_gpu and not gpu_stage and a.project != "diet":
        raise ValueError("--job-gpus is supported only by PFI train/predict and stages 09b/09c")
    if opts.job_gpus < 0:
        raise ValueError("--job-gpus must be nonnegative")
    budget = core_budget(opts.cores)
    os.environ['WORKFLOW_CORES'] = str(budget)
    if a.project == 'pfi':
        for key in ('N_JOBS', 'PMC_PARSE_WORKERS', 'PFI_FEATURE_WORKERS', 'PFI_DUCKDB_THREADS'):
            os.environ.setdefault(key, str(budget))
        os.environ.setdefault('PFI_BUILD_JOBS', str(min(2, budget)))
        os.environ.setdefault('PFI_BUILD_DUCKDB_THREADS', str(max(1, budget // int(os.environ['PFI_BUILD_JOBS']))))
    diet_workflow = None
    if a.project == "diet":
        sys.path.insert(0, str(Path(a.entry).resolve().parent / "f"))
        from pipeline import parser as diet_parser, Workflow as DietWorkflow, validate_args
        diet_entry = "disco" if entry_name == "1.disco.sh" else "diet"
        if opts.cores is not None:
            scientific += ["--cores", str(opts.cores)]
        dp = diet_parser()
        worker_args = dp.parse_args(["--entry", diet_entry, *scientific])
        validate_args(dp, worker_args)
        if not any(v.split("=", 1)[0] in {"--cores", "--jobs"} for v in scientific):
            scientific += ["--cores", str(worker_args.jobs)]
        diet_workflow = DietWorkflow(worker_args)
        if not opts.plan_only:
            diet_workflow.preflight(diet_workflow.selected())
        required_cores, required_memory = worker_args.jobs, 0
    else:
        required_cores, required_memory = scientific_budget(a.project, scientific)
    if a.project in {'ckm', 'pfi', 'nextstrain'}:
        # Only explicit worker settings constrain the shared budget.
        for i, token in enumerate(scientific):
            key, sep, value = token.partition('=')
            if key in {'--threads', '--jobs', '--bgzip-threads'} and int(value if sep else scientific[i+1]) > budget:
                raise ValueError(f'{key} exceeds --cores')
        if a.project == 'pfi' and (any(int(os.environ[k]) > budget for k in ('N_JOBS', 'PMC_PARSE_WORKERS', 'PFI_FEATURE_WORKERS', 'PFI_DUCKDB_THREADS')) or int(os.environ['PFI_BUILD_JOBS']) * int(os.environ['PFI_BUILD_DUCKDB_THREADS']) > budget):
            raise ValueError('PFI worker settings exceed --cores')
        required_cores = budget
        if a.project == 'ckm': os.environ['CKM_THREADS'] = str(budget)
        if a.project == 'nextstrain': os.environ['NEXTSTRAIN_THREADS'] = str(budget)
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
    if diet_workflow is not None:
        root = diet_workflow.out
    if entry_name.endswith(".paper.sh") and (a.project != "gu" or not any(t.split("=", 1)[0] == "--dir-out" for t in scientific)):
        root = root / ("paper" if a.project == "diet" else Path(entry_name).stem)
    queue = opts.queue or (os.getenv("HPC_GPU_QUEUE", "4a100-40") if opts.job_gpus else os.getenv("HPC_CPU_QUEUE", "33"))
    env = {"WORKFLOW_CORES": str(budget), "PROJECT_JOB_ACTIVE": "1", "WORKFLOW_RESUME_POLICY": "new_plan",
           "WORKFLOW_DEVICE": "gpu" if opts.job_gpus else "cpu"}
    # Freeze only workflow/project settings, never the full environment.
    prefix = a.project.upper().replace("-", "_") + "_"
    env.update({k: v for k, v in os.environ.items() if k.startswith(prefix)})
    if a.project == "pfi":
        env.update({k: os.environ[k] for k in ("N_JOBS", "PMC_PARSE_WORKERS") if k in os.environ})
    if a.project == "diet":
        env["WORKFLOW_DEVICE"] = "auto" if mode == "local" else ("cuda" if opts.job_gpus else "cpu")
    job_name = "diet.disco" if a.project == "diet" and entry_name == "1.disco.sh" else a.project
    job = Job(a.project, ["bash", str(Path(a.entry).resolve()), *scientific], job_name=job_name,
              resources={"queue": queue, "cores": opts.job_cores, "memory_gib": opts.job_memory_gb,
                         "hours": opts.job_hours, "gpus": opts.job_gpus}, env=env, cwd=os.getcwd())
    return Workflow([job], opts.dir_jobs or root / 'logs', a.project, mode, plan_only=opts.plan_only)


if __name__ == '__main__':
    result = execute_workflow(build_workflow())
    sys.exit(int(result.get('exit_code', 0)))
