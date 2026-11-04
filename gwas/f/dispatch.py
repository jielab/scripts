#!/usr/bin/env python3
"""Schedule the existing GWAS pipeline as a bounded local or LSF compute job.

Only standard-library planning occurs here. No genotype/reference/GWAS scans
or statistical work are performed on a login node. The scientific driver keeps
its chromosome/trait order and its existing validated completion markers.
"""
from __future__ import annotations
import math
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "0f"))
import workflow_jobs as workflow
from qiming_resources import cpu_queue


def option(argv, name, default=None):
    value = default
    for i, token in enumerate(argv):
        if token == name:
            if i + 1 >= len(argv) or argv[i + 1].startswith("--"):
                raise ValueError(f"{name} requires a value")
            value = argv[i + 1]
        elif token.startswith(name + "="):
            value = token.split("=", 1)[1]
    return value


def set_option(argv, name, value):
    out, i = [], 0
    while i < len(argv):
        if argv[i] == name:
            i += 2
        elif argv[i].startswith(name + "="):
            i += 1
        else:
            out.append(argv[i]); i += 1
    return out + [name, str(value)]


def integer(value, name):
    if not str(value).isdigit() or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def build(entry, argv):
    if entry not in {"format.sh", "gwas.sh", "compare.sh"}:
        raise ValueError("Unknown GWAS entry")
    if entry == "format.sh":
        # Its historical parser permits the module before or after valued options.
        modules, i = [], 0
        while i < len(argv):
            if argv[i].startswith("--"):
                if "=" in argv[i]:
                    i += 1
                else:
                    if i + 1 >= len(argv) or argv[i + 1].startswith("--"):
                        raise ValueError(f"{argv[i]} requires a value")
                    i += 2
            else:
                modules.append(argv[i]); i += 1
        if len(modules) > 1:
            raise ValueError("Specify one format module or one comma-separated module list")
        module = modules[0] if modules else "all"
    else:
        module = argv[0] if argv and not argv[0].startswith("-") else ""
    allowed = {"format.sh": {"format", "thin", "magma", "liftover", "cis", "lead", "mplot", "h2", "pgs", "all"},
               "gwas.sh": {"sparseGRM", "acount", "regenie", "run_plink2", "run_saige", "prep_gwas", "account", "prep_regenie", "filter_variants", "run_regenie"},
               "compare.sh": {"compare", "ldsc", "shiny"}}
    if not module or set(module.split(",")) - allowed[entry] or (entry != "format.sh" and "," in module):
        raise ValueError(f"Invalid {entry} module: {module}")
    mode = os.environ.get("submit", "local")
    if mode not in {"local", "hpc"}:
        raise ValueError("submit must be local or hpc")
    run_default = "FALSE" if entry == "format.sh" and "magma" not in module.split(",") else "TRUE"
    run = option(argv, "--run-cmd", run_default).upper()
    if entry == "compare.sh" and module == "ldsc":
        run = option(argv, "--run", "TRUE").upper()
    if run not in {"TRUE", "FALSE"}:
        raise ValueError("--run-cmd/--run must be TRUE or FALSE")
    dry = os.environ.get("GWAS_PLAN_ONLY", "FALSE").upper() == "TRUE" or run == "FALSE"
    count = integer(option(argv, "--jobs", "8" if entry == "gwas.sh" else "4" if entry == "format.sh" else "1"), "--jobs")
    per_process = 1
    if entry == "gwas.sh":
        per_process = max(integer(option(argv, "--threads", "4"), "--threads"),
                          integer(option(argv, "--step2-threads", option(argv, "--threads", "4")), "--step2-threads"))
    elif entry == "format.sh":
        per_process = max(2, integer(option(argv, "--pgs-threads", "1"), "--pgs-threads"),
                          integer(os.environ.get("GWAS_TOOL_THREADS", "2"), "--tool-threads"),
                          integer(os.environ.get("GWAS_EAF_THREADS", "2"), "GWAS_EAF_THREADS"),
                          integer(os.environ.get("GWAS_CLEAN_COMP_THREADS", "1"), "GWAS_CLEAN_COMP_THREADS"))
    required_cores = count * per_process
    queue = os.environ.get("GWAS_QUEUE", "auto" if mode == "hpc" else os.environ.get("HPC_CPU_QUEUE", "33"))
    cores = integer(os.environ.get("GWAS_CORES", str(max(40 if mode == "hpc" and queue == "short" else 1, required_cores))), "--cores")
    if cores < required_cores or (mode == "hpc" and queue == "short" and cores < 40):
        raise ValueError(f"--cores must cover {required_cores} concurrent worker threads; short requires >=40")
    minimum_memory = 1
    if entry == "gwas.sh":
        minimum_memory = math.ceil(integer(option(argv, "--memory", "24000"), "--memory") * count / 1024 + 4)
    elif entry == "format.sh":
        plink_mb = integer(os.environ.get("GWAS_PLINK_MEMORY_MB", "8192"), "--tool-memory-mb")
        if plink_mb < 640:
            raise ValueError("--tool-memory-mb must be at least 640 for PLINK")
        minimum_memory = math.ceil(max(plink_mb,
            integer(os.environ.get("PGS_PLINK_MEMORY_MB", "8192"), "PGS_PLINK_MEMORY_MB"),
            integer(os.environ.get("GWAS_EAF_MEMORY_MB", "2048"), "GWAS_EAF_MEMORY_MB")) * count / 1024 + 4)
    memory = float(os.environ.get("GWAS_MEMORY_GIB", str(max(minimum_memory, int(os.environ.get("MEMORY_DEFAULT_GIB", "160" if mode == "hpc" else "32"))))))
    if memory < minimum_memory:
        raise ValueError(f"--memory-gb must be at least {minimum_memory} for the requested PLINK workers; reduce --jobs/--memory")
    hours = float(os.environ.get("GWAS_HOURS", "120"))
    if mode == "hpc":
        queue = cpu_queue(queue, cores, memory, hours)
    elif queue == "auto":
        queue = os.environ.get("HPC_CPU_QUEUE", "33")
    resources = workflow.resources(dict(queue=queue, cores=cores, memory_gib=memory,
                                       hours=hours,
                                       memory_scope=os.environ.get("GWAS_MEMORY_SCOPE", "site_host")))
    if os.environ.get("GWAS_RUNTIME_ENV") and not Path(os.environ["GWAS_RUNTIME_ENV"]).is_file():
        raise ValueError("--runtime-env must be a readable file on the shared filesystem")
    worker_args = list(argv)
    if entry in {"format.sh", "gwas.sh"}:
        worker_args = set_option(worker_args, "--run-cmd", "TRUE")
    if entry == "format.sh":
        worker_args = set_option(worker_args, "--foreground", "TRUE")
    if entry == "compare.sh" and module == "ldsc":
        worker_args = set_option(worker_args, "--run", "TRUE")
    keep = {"DIR_WORK", "DIR_DATA", "DIR_SCRATCH", "DIR_SOFTWARE", "DATE_SCRATCH", "SHARED_PYTHON",
            "WORKFLOW_PYTHON", "PATH", "CONDA_PREFIX", "CONDA_DEFAULT_ENV", "CONDA_EXE", "LD_LIBRARY_PATH",
            "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "SAIGE_ENV_DIR", "CHECK_GRCH_SNP_LIST", "PYTHONPATH", "PGS_PLINK_MEMORY_MB"}
    env = {k: v for k, v in os.environ.items() if k in keep or k.startswith("GWAS_")}
    env.update(GWAS_WORKER="1", GWAS_PLAN_ONLY="FALSE", GWAS_PROCESS_THREADS=str(per_process), SCRIPT_CONSOLE_PASSTHROUGH="1",
               # Bound implicit OpenMP/R workers as well as explicit PLINK threads.
               OMP_NUM_THREADS=str(per_process), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               NUMEXPR_NUM_THREADS="1", R_DATATABLE_NUM_THREADS=str(per_process))
    job = workflow.Job(key=module.replace(",", "_"), argv=["bash", str(ROOT / "f/worker.sh"), str(ROOT / entry), *worker_args],
                       resources=resources, env=env, cwd=str(Path.cwd()))
    return mode, dry, job


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    in_allocation = bool(args and args[0] == "--in-allocation")
    if in_allocation:
        args.pop(0)
    if not args:
        raise ValueError("An entry script is required")
    entry = args.pop(0)
    mode, dry, job = build(entry, args)
    if in_allocation and not dry:
        if mode != "hpc" or not os.environ.get("LSB_JOBID"):
            raise ValueError("--in-allocation requires an LSF job")
        slots = int(os.environ.get("LSB_DJOB_NUMPROC", "0"))
        if slots < job.resources["cores"]:
            raise ValueError(f"Insufficient allocated CPU slots ({slots}); need {job.resources['cores']}. Reduce --jobs/--threads or request a larger job.")
        if os.environ.get("WORKFLOW_MEMORY_GIB") and float(os.environ["WORKFLOW_MEMORY_GIB"]) < job.resources["memory_gib"]:
            raise ValueError("Allocated workflow memory is smaller than the requested GWAS budget")
        os.execvpe(job.argv[0], job.argv, {**os.environ, **job.env})
    directory = Path(os.environ.get("GWAS_PLAN_DIR", str(Path(os.environ["GWAS_ANALYSIS_ROOT"]) / "jobs")))
    path = workflow.write_plan([job], directory, mode=mode, name="gwas-" + entry.removesuffix(".sh"))
    print(f"[GWAS] {job.resources['cores']} CPU slots; {job.resources['memory_gib']:g} GiB total; {job.resources['hours']:g} h", flush=True)
    if dry:
        print(f"[GWAS] Plan only. Execute the reviewed plan with: {sys.executable} {ROOT.parent / '0f/workflow_jobs.py'} run {path} --submit {mode} --wait", flush=True)
        return 0
    workflow.run_plan(path, submit=mode, jobs=1, wait=os.environ.get("SUBMIT_WAIT") == "TRUE")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
