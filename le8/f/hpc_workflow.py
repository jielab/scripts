#!/usr/bin/env python3
"""Lightweight LE8 local/LSF scope planner: standard library only, no cohort reads.

All fitting, native R preparation, dependency checks and aggregation are workers.
One trait/layer is one job; its internal stages run sequentially in that allocation.
"""
from __future__ import annotations

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
# Keep project modules (notably abm.py) ahead of shared helpers.
_paths_sys.path.append(str(_paths_root / "0f"))
from workflow_paths import configured_path

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT.parent / "0f"))
from workflow_jobs import Job, write_plan, run_plan

DEFAULT = "c1_correlate,c2_cause,c3_coloc,c4_connect,c5_cellulation"
MODULES = {"c1_correlate", "c1_reliability", "c2_cause", "c3_coloc", "c4_connect",
           "c5_cellulation", "pgs_focus", "final", "share", "shiny", "install"}


def bool_word(value):
    if value.lower() in {"true", "1", "yes", "y"}:
        return "TRUE"
    if value.lower() in {"false", "0", "no", "n"}:
        return "FALSE"
    raise argparse.ArgumentTypeError("requires TRUE or FALSE")


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.set_defaults(modules=DEFAULT)
    p.add_argument("--Y", "--trait", default=os.getenv("Y", "cvd_cad,ra"))
    p.add_argument("--biom", default=os.getenv("BIOM", "prot,met"))
    p.add_argument("--analysis-root", "--dir-analysis-root", type=Path,
                   default=Path(os.environ.get("LE8_ANALYSIS_ROOT", configured_path('ANALYSIS_ROOT', 'le8'))))
    p.add_argument("--cores", type=int, default=int(os.getenv("N_CORES", "8")))
    p.add_argument("--memory-limit-gb", type=float, default=float(os.getenv("LE8_MEMORY_LIMIT_GB", "32")))
    p.add_argument("--queue", default=os.getenv("LE8_QUEUE", "33"))
    p.add_argument("--gpu-queue", default=os.getenv("LE8_GPU_QUEUE", os.getenv("HPC_GPU_QUEUE", "4a100-40")))
    p.add_argument("--gpu-memory-limit-gb", type=float, default=float(os.getenv("HPC_GPU_MEMORY_GIB", "160")))
    p.add_argument("--gpu-cores", type=int, default=int(os.getenv("LE8_GPU_CORES", "20")))
    p.add_argument("--hours", type=float, default=float(os.getenv("LE8_HOURS", "72")))
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--eigh-backend", choices=["auto", "cpu", "cuda"], default=os.getenv("MRLINK2_EIGH_BACKEND", "auto"))
    p.add_argument("--run-gpu-coloc", type=bool_word, choices=["TRUE", "FALSE"],
                   help="Explicit C3 GPU switch; overrides --device and RUN_GPU_COLOC for C3")
    p.add_argument("--python", default=os.getenv("LE8_REPORT_PYTHON", os.getenv("PYTHON_BIN", "python3")))
    p.add_argument("--r-bin", default=os.getenv("R_BIN", "Rscript"))
    p.add_argument("--methods")
    p.add_argument("--config", type=Path, default=ROOT / "config/c1_methods.json")
    p.add_argument("--skip-abm", action="store_true")
    p.add_argument("--run-abm", action="store_true", help="Legacy alias: include C1 in the selected stages")
    p.add_argument("--plan-dir", type=Path)
    p.add_argument("--resume-plan", type=Path)
    p.add_argument("--plan-only", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--preflight", action="store_true")
    for name in ("fit-joint", "fit-genetic", "fit-reference", "details", "index-only"):
        p.add_argument("--" + name, action="store_true")
    p.add_argument("--summary-worker", action="store_true", help=argparse.SUPPRESS)
    return p


def parse_cli(argv=None):
    # An optional argparse positional can steal an unknown scientific option's
    # value (e.g. --prepared-source /path). Only the first token can be modules.
    tokens = list(sys.argv[1:] if argv is None else argv)
    modules = tokens.pop(0) if tokens and not tokens[0].startswith("-") else DEFAULT
    args, extra = parser().parse_known_args(tokens)
    args.modules = modules
    return args, extra


def selected_methods(a):
    if a.skip_abm:
        return ["glmnet", "lightgbm"]
    if a.methods:
        methods = a.methods.split(",")
    else:
        config = json.loads(a.config.read_text())
        methods = config.get("methods", ["glmnet", "lightgbm", "abm_count", "abm_deep"])
    allowed = {"glmnet", "lightgbm", "abm_count", "abm_legacy", "abm_case_preserving", "abm_profile",
               "abm_lenient", "abm_prevalent", "abm_deep", "molecular_deviation"}
    aliases = {"pradeep": "glmnet", "lgbm": "lightgbm", "abm_multiview": "abm_profile"}
    methods = list(dict.fromkeys(aliases.get(x, x) for x in methods))
    if not methods or not set(methods) <= allowed:
        raise ValueError("Invalid C1 methods")
    return methods


def build_jobs(a, extra=()):
    """Return Job records without importing scientific code or reading UKB data."""
    modules = list(dict.fromkeys(a.modules.split(",")))
    old_aliases = {"c1_abm": ["abm_legacy", "molecular_deviation"], "c1_abm_GPU": ["abm_deep"]}
    legacy = [m for m in modules if m in old_aliases]
    inherited = [method for module in legacy for method in old_aliases[module]]
    had_c1 = "c1_correlate" in modules
    if legacy:
        modules = [m for m in modules if m not in old_aliases]
        if "c1_correlate" not in modules:
            modules.insert(0, "c1_correlate")
    if a.run_abm and "c1_correlate" not in modules:
        modules.insert(0, "c1_correlate")
    if set(modules) - MODULES:
        raise ValueError("Unknown modules: " + ",".join(sorted(set(modules) - MODULES)))
    traits, layers = list(dict.fromkeys(a.Y.split(","))), list(dict.fromkeys(a.biom.split(",")))
    if not traits or any(not re.fullmatch(r"[A-Za-z0-9_]+", t) for t in traits):
        raise ValueError("Invalid trait list")
    if not layers or not set(layers) <= {"prot", "met"}:
        raise ValueError("Invalid biomarker layers")
    if a.cores < 1 or a.gpu_cores < 1 or a.memory_limit_gb <= 0 or a.gpu_memory_limit_gb <= 0 or a.hours <= 0:
        raise ValueError("Resource requests must be positive")
    if not float(a.memory_limit_gb).is_integer():
        raise ValueError("LE8 --memory-limit-gb uses whole GiB, matching the native R launcher")
    if "install" in modules:
        if modules != ["install"]:
            raise ValueError("Install environments separately from scientific analysis")
        env = {k: v for k, v in os.environ.items() if k.startswith(("DIR_", "LE8_", "OS_"))
               or k in {"PATH", "CONDA_EXE", "CONDA_SH", "ANALYSIS_ROOT", "DATE_SCRATCH", "SHARED_PYTHON"}}
        env.update(PYTHON_BIN=a.python, LE8_REPORT_PYTHON=a.python, R_BIN=a.r_bin,
                   SCRIPT_LOG_DIR=str(a.analysis_root.resolve() / "logs/install"))
        return [Job("install", ["bash", str(ROOT / "install.sh"), "--python", a.python,
                    "--r-bin", a.r_bin, *map(str, extra)], [],
                    dict(queue=a.queue, cores=a.cores, memory_gib=a.memory_limit_gb, hours=a.hours, gpus=0),
                    env, str(ROOT))]
    methods = selected_methods(a)
    if legacy and not a.methods and not a.skip_abm:
        methods = list(dict.fromkeys((methods if had_c1 else []) + inherited))
    root = str(a.analysis_root.resolve())
    common = ["--analysis-root", root, "--cores", str(a.cores),
              "--memory-limit-gb", f"{a.memory_limit_gb:g}", "--r-bin", a.r_bin,
              "--device", a.device, "--config", str(a.config.resolve()), *map(str, extra)]
    if a.preflight:
        common.append("--preflight")
    jobs = []
    scientific_env = {k: v for k, v in os.environ.items() if
        k.startswith(("LE8_", "C1_", "C2_", "C3_", "C4_", "C5_", "FINAL_", "MRLINK2_", "PGS_", "RUN_"))
        or k in {"DIR_WORK", "DIR_DATA", "DIR_SCRATCH", "DIR_SOFTWARE", "ANALYSIS_ROOT", "DATE_SCRATCH",
                 "UKB_PHE", "UKB_COMMON", "UKB_BGEN_DIR", "PHE_F", "CONDA_SH", "CONDA_EXE", "PATH", "LD_LIBRARY_PATH"}}
    scientific_env.update(PYTHON_BIN=a.python, ABM_PYTHON=a.python, LE8_REPORT_PYTHON=a.python,
        R_BIN=a.r_bin, LE8_TOTAL_CORES=str(a.cores), LE8_RESOURCE_MEMORY_GIB=f"{a.memory_limit_gb:g}",
        LE8_MEMORY_LIMIT_GB=f"{a.memory_limit_gb:g}", LE8_RUNTIME_READY="1",
        LE8_COHORT_CACHE_ROOT=str(a.analysis_root.resolve() / "logs/cohort_cache"),
        SCRIPT_LOG_DIR=str(a.analysis_root.resolve() / "logs/runs"),
        SCRIPT_FOREGROUND="1", SCRIPT_CONSOLE_ACTIVE="1")

    def add(key, module, deps, trait=None, layer=None, tail=(), stage_env=None, gpu=False):
        env = dict(scientific_env, RUN_GPU_COLOC="FALSE", MRLINK2_EIGH_BACKEND="cpu")
        allocated_cores = max(a.cores, a.gpu_cores) if gpu else a.cores
        env["LE8_TOTAL_CORES"] = str(allocated_cores)
        if trait and layer:
            env["LE8_CACHE_ROOT"] = str(a.analysis_root.resolve() / "logs/cache" / f"{trait}.{layer}")
        if gpu and module == "c3_coloc":
            env["RUN_GPU_COLOC"] = "TRUE"
        if gpu and module == "c2_cause":
            env["MRLINK2_EIGH_BACKEND"] = "cuda"
        memory = min(a.memory_limit_gb, a.gpu_memory_limit_gb) if gpu else a.memory_limit_gb
        env.update(LE8_RESOURCE_MEMORY_GIB=f"{memory:g}", LE8_MEMORY_LIMIT_GB=f"{memory:g}")
        env.update(stage_env or {})
        argv = ["bash", str(ROOT / "le8.sh"), module, *common,
                "--Y", trait or a.Y, "--biom", layer or a.biom, "--cores", str(allocated_cores), "--memory-limit-gb", f"{memory:g}", *tail]
        if module == "c1_correlate":
            argv += ["--methods", ",".join(methods)]
            if gpu:
                # An allocated deep-model GPU must be usable; a broken CUDA
                # runtime must not silently turn the requested GPU fit into CPU.
                argv += ["--device", "cuda"]
        jobs.append(Job(key, argv, list(deps), dict(queue=a.gpu_queue if gpu else a.queue,
            cores=allocated_cores, memory_gib=memory, hours=a.hours, gpus=int(gpu)), env, str(ROOT)))
        return key

    # Build the inspectable internal commands first, then bundle them by scope.
    # They are subprocesses inside one allocation, never additional LSF jobs.
    scopes = []
    batch_id = uuid.uuid4().hex
    scope_keys = [f"{trait}.{layer}" for trait in traits for layer in layers]
    for trait in traits:
        for layer in layers:
            scope = f"{trait}.{layer}"
            start = len(jobs)
            deps = []
            if "c1_correlate" in modules:
                for stage in ("prepare", "association", "models", "evaluation"):
                    deps = [add(f"{scope}.c1.{stage}", "c1_correlate", deps, trait, layer,
                        ["--c1-stage", stage, *(["--replace-models", "FALSE"] if stage == "evaluation" else [])],
                        gpu=stage == "models" and "abm_deep" in methods and
                            (a.device == "cuda" or (a.device == "auto" and os.getenv("submit") == "hpc")))]
            if "c1_reliability" in modules or ("c1_correlate" in modules and "glmnet" in methods):
                deps = [add(f"{scope}.c1.reliability", "c1_reliability", deps, trait, layer)]
            for module in ("c2_cause", "c3_coloc", "c4_connect", "c5_cellulation", "pgs_focus"):
                if module not in modules:
                    continue
                stages = ("c4_connect", "c4_panel_validation", "c4_explain") if module == "c4_connect" else (module,)
                for stage in stages:
                    # Explicit C3 switch wins. Otherwise an explicit CPU device
                    # must override a site default that enables GPU colocalization.
                    coloc_gpu = (a.run_gpu_coloc == "TRUE" if a.run_gpu_coloc is not None else
                                 (a.device != "cpu" and
                                  (a.device == "cuda" or os.getenv("RUN_GPU_COLOC") == "TRUE")))
                    eigh_gpu = (a.eigh_backend == "cuda" or
                                (a.eigh_backend == "auto" and a.device != "cpu" and
                                 (a.device == "cuda" or os.getenv("submit") == "hpc")))
                    gpu = ((module == "c3_coloc" and coloc_gpu)
                           or (module == "c2_cause" and eigh_gpu))
                    deps = [add(f"{scope}.{stage}", module, deps, trait, layer,
                        stage_env={"LE8_C4_STAGE": stage} if module == "c4_connect" else None, gpu=gpu)]
            if any(m in modules for m in ("c1_correlate", "c1_reliability")):
                argv = [a.python, str(Path(__file__).resolve()), "--summary-worker", "--analysis-root", root,
                        "--Y", a.Y, "--biom", a.biom]
                jobs.append(Job(f"{scope}.c1.summary", argv, [],
                    dict(queue=a.queue, cores=1, memory_gib=min(a.memory_limit_gb, 8), hours=a.hours, gpus=0),
                    scientific_env, str(ROOT)))
            for module in ("final", "share", "shiny"):
                if module in modules:
                    tail = ["--prepare-only"] if module == "shiny" and os.getenv("submit", "local") == "hpc" else []
                    if module == "final":
                        tail += ["--" + flag for flag in ("fit-joint", "fit-genetic", "fit-reference", "details", "index-only")
                                 if getattr(a, flag.replace("-", "_"))]
                    # The last finished scope performs the original aggregate
                    # request, including optional joint-layer fits, in its job.
                    add(f"{scope}.{module}", module, [], tail=tail)
            stages = jobs[start:]
            gpu = any(j.resources["gpus"] for j in stages)
            allocated_cores = max(j.resources["cores"] for j in stages)
            memory = max(j.resources["memory_gib"] for j in stages)
            if gpu and memory > a.gpu_memory_limit_gb:
                raise ValueError("The single scope job needs more host RAM than --gpu-memory-limit-gb; "
                                 "increase that cap or select CPU backends")
            env = dict(scientific_env, LE8_TOTAL_CORES=str(allocated_cores),
                       LE8_RESOURCE_MEMORY_GIB=f"{memory:g}", LE8_MEMORY_LIMIT_GB=f"{memory:g}",
                       RUN_GPU_COLOC="TRUE" if any(j.env.get("RUN_GPU_COLOC") == "TRUE" for j in stages) else "FALSE",
                       MRLINK2_EIGH_BACKEND="cuda" if any(j.env.get("MRLINK2_EIGH_BACKEND") == "cuda" for j in stages) else "cpu")
            steps = [dict(key=j.key, argv=j.argv, env={k:v for k,v in j.env.items() if env.get(k) != v},
                          aggregate=j.key.endswith((".c1.summary", ".final", ".share", ".shiny"))) for j in stages]
            argv = [a.python, str(ROOT / "f/scope_workflow.py"), "--analysis-root", root,
                    "--scope", scope, "--batch-id", batch_id, "--scopes", ",".join(scope_keys),
                    "--steps-json", json.dumps(steps, separators=(",", ":"))]
            scopes.append(Job(scope, argv, [], dict(queue=a.gpu_queue if gpu else a.queue,
                cores=allocated_cores, memory_gib=memory, hours=a.hours, gpus=int(gpu)), env, str(ROOT)))
    return scopes


def main(argv=None):
    a, extra = parse_cli(argv)
    if a.summary_worker:
        if os.getenv("submit") == "hpc" and not os.getenv("LSB_JOBID"):
            raise RuntimeError("Summary worker requires an LSF allocation")
        from c1_report import write_reliability_cross_trait_summary, write_pgs_cross_trait_summary
        for layer in a.biom.split(","):
            write_reliability_cross_trait_summary(a.analysis_root, layer)
            write_pgs_cross_trait_summary(a.analysis_root, layer)
        return 0
    mode = os.getenv("submit", "local")
    if a.resume_plan:
        path = a.resume_plan
    else:
        jobs = build_jobs(a, extra)
        path = write_plan(jobs, a.plan_dir or a.analysis_root / "logs/workflows", mode=mode, name="le8")
        # Keep the actual .cmd reviewable without expanding a huge quoted JSON
        # argument in its comment. Execution still uses the immutable plan.
        for job in jobs:
            if "--steps-json" not in job.argv:
                continue
            steps = json.loads(job.argv[job.argv.index("--steps-json") + 1])
            command = path.parent / "commands" / (job.key + ".cmd")
            lines = []
            for line in command.read_text().splitlines():
                if line.startswith("# Scientific command:"):
                    lines.append("# Sequential stages inside this single allocation:")
                    for step in steps:
                        lines.append("# " + step["key"] + (" (last scope aggregates)" if step["aggregate"] else ""))
                        lines.append("# " + shlex.join(["env", *[f"{k}={v}" for k,v in sorted(step["env"].items())], *step["argv"]]))
                else:
                    lines.append(line)
            command.write_text("\n".join(lines) + "\n")
    if a.plan_only or a.dry_run:
        print(f"[LE8] Reviewable plan written; no analysis or submission: {path}", flush=True)
        return 0
    result = run_plan(path, submit=mode, jobs=int(os.getenv("jobs_N", "1")),
                      wait=os.getenv("SUBMIT_WAIT", "FALSE") == "TRUE")
    print(json.dumps(result, indent=2), flush=True)
    return int(result.get("exit_code", 0))


if __name__ == "__main__":
    sys.exit(main())
