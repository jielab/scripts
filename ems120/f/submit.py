#!/usr/bin/env python3
"""Lightweight EMS120 planner; scientific libraries load only in workers."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "0f"))
from workflow_jobs import Job, run_plan, write_plan
from workflow_layout import run_name

STEPS = ("all", "fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "figs")
METHODS = ("phone", "geo", "dx")
GPU_LIMITS = {"4a100-40": (56, 384), "4a100-80": (32, 128), "2a100-40": (40, 192),
              "2a100-80": (40, 192), "2h100-80": (192, 768), "hgx": (64, 1024),
              "8v100-32-sc": (48, 384), "4v100-16-e5": (24, 64),
              "4v100-16-sc": (40, 192), "2v100-16-sc": (40, 384),
              "2v100-32-e5": (24, 256), "2v100-32-sc": (40, 192)}
CPU_QUEUES = {"33", "38", "ot38", "2t50c", "1t75c", "1t88c", "1t81c", "73x", "52", "63",
              "v4-3t", "sk-6t-104", "75m", "v3-6t", "v3-64", "spec-cpu"}


def boolean(value):
    if str(value).upper() in {"TRUE", "FALSE"}:
        return str(value).upper() == "TRUE"
    raise argparse.ArgumentTypeError("expected TRUE or FALSE")


def absolute(value):
    value = os.fspath(value)
    if not value or any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError("A path must be nonempty and contain no control characters")
    return str(Path(value).expanduser().resolve())


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument("--plan-only", action="store_true", help="write plan/.cmd without running or submitting")
    action = p.add_mutually_exclusive_group()
    action.add_argument("--report-only", action="store_true", help="CPU redraw from complete data; never prepare or fit models")
    action.add_argument("--preflight", action="store_true", help="allocated checks using input headers and synthetic training; no production analysis")
    for method in METHODS:
        p.add_argument(f"--{method}-run", choices=("cpu", "gpu", "auto"), default="cpu",
                       help=f"{method} processor; cpu by default, auto resolves to cpu")
    p.add_argument("--replace", type=boolean, default=True, help="TRUE overwrites requested outputs; FALSE reuses complete outputs")
    local = os.getenv("submit") == "local"
    p.add_argument("--cpu-queue", default=os.getenv("EMS120_CPU_QUEUE", "33"))
    p.add_argument("--cpu-cores", type=int, default=int(os.getenv("EMS120_CPU_CORES", str(min(8, os.cpu_count() or 1)) if local else "20")))
    p.add_argument("--cpu-threads", type=int, default=os.getenv("EMS120_CPU_THREADS"), help="default min(8, cpu-cores)")
    p.add_argument("--cpu-memory-gb", type=float, default=float(os.getenv("EMS120_CPU_MEMORY_GB", "32" if local else "96")))
    p.add_argument("--cpu-hours", type=float, default=float(os.getenv("EMS120_CPU_HOURS", "120")))
    p.add_argument("--gpu-queue", default=os.getenv("EMS120_GPU_QUEUE", "2v100-32-sc"))
    p.add_argument("--gpu-cores", type=int, default=int(os.getenv("EMS120_GPU_CORES", "8")))
    p.add_argument("--gpu-memory-gb", type=float, default=float(os.getenv("EMS120_GPU_MEMORY_GB", "32" if local else "96")))
    p.add_argument("--gpu-hours", type=float, default=float(os.getenv("EMS120_GPU_HOURS", "120")))
    p.add_argument("--gpu-model", default=os.getenv("EMS120_GPU_MODEL", ""), help="optional LSF model constraint; one GPU per selected method")
    for flag, key, description in (
        ("dir-data", "DATA_DIR", "EMS input directory containing clean/ and szu_xia/"),
        ("dir-analysis", "ANALYSIS_DIR", "analysis outputs"),
        ("dir-base-model", "BASE_MODEL_DIR", "offline pretrained HFL/MacBERT directory"),
        ("dir-model", "MODEL_DIR", "trained-model root, default <analysis>/ml"),
        ("dir-files", "FILES_DIR", "figure assets"),
        ("dir-jobs", "JOBS_DIR", "plans/logs/status, default <analysis>/logs"),
    ):
        p.add_argument(f"--{flag}", default=os.getenv(f"EMS120_{key}"), help=description)
    p.add_argument("--conda-env", default=os.getenv("EMS120_CONDA_ENV", "ai"))
    p.add_argument("--conda-exe", default=os.getenv("EMS120_CONDA_EXE", os.getenv("CONDA_EXE", "")))
    p.add_argument("--rscript", default=os.getenv("EMS120_RSCRIPT", ""))
    p.add_argument("--python-lib", default=os.getenv("EMS120_PYTHON_LIB", ""))
    p.add_argument("--start-step", type=str.lower, choices=STEPS, default="all")
    p.add_argument("--end-step", type=str.lower, choices=STEPS[1:], default="")
    p.add_argument("--run-cv", type=boolean, default=True)
    g = p.add_argument_group("shared options (use ems120.sh)")
    g.add_argument("--submit", choices=("auto", "hpc", "local"))
    g.add_argument("--dir-dat")
    g.add_argument("--dir-big")
    g.add_argument("--wait", choices=("TRUE", "FALSE"))
    return p


def configuration(a):
    mode = os.environ.get("submit", "hpc")
    if mode not in {"hpc", "local"}:
        raise ValueError("--submit must resolve to hpc or local")
    if a.submit or a.dir_dat or a.dir_big or a.wait:
        raise ValueError("Use ems120.sh so shared options are resolved through 0os.sh")
    for method in METHODS:
        if getattr(a, method + "_run") == "auto":
            setattr(a, method + "_run", "cpu")
    if a.report_only and any(getattr(a, m + "_run") != "cpu" for m in METHODS):
        raise ValueError("--report-only uses CPU and cannot request a model GPU")
    if a.cpu_threads is None:
        a.cpu_threads = min(8, a.cpu_cores)
    for prefix in ("cpu", "gpu"):
        cores, memory, hours, queue = (getattr(a, prefix + "_" + name) for name in ("cores", "memory_gb", "hours", "queue"))
        if not 1 <= cores <= 192:
            raise ValueError(f"--{prefix}-cores must be in 1..192")
        if not math.isfinite(memory) or memory <= 0:
            raise ValueError(f"--{prefix}-memory-gb must be positive and finite")
        if not math.isfinite(hours) or not 0 < hours <= 336:
            raise ValueError(f"--{prefix}-hours must be in (0, 336]")
        if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}", queue):
            raise ValueError(f"Invalid {prefix} queue name")
    if not 1 <= a.cpu_threads <= a.cpu_cores:
        raise ValueError("--cpu-threads must be within allocated --cpu-cores")
    if a.gpu_model and not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}", a.gpu_model):
        raise ValueError("Invalid GPU model constraint")
    if mode == "hpc":
        if a.cpu_queue in GPU_LIMITS:
            raise ValueError("Select a CPU queue for --cpu-queue")
        if a.cpu_queue == "33" and (a.cpu_cores not in {10, 20, 30, 40} or a.cpu_memory_gb > 192 * a.cpu_cores / 40):
            raise ValueError("Queue 33 needs 10/20/30/40 slots with at most 48/96/144/192 GiB respectively")
        if any(getattr(a, m + "_run") == "gpu" for m in METHODS):
            if a.gpu_queue in CPU_QUEUES:
                raise ValueError("Select a GPU queue for --gpu-queue")
            if a.gpu_queue in GPU_LIMITS:
                cores, ram = GPU_LIMITS[a.gpu_queue]
                if a.gpu_cores > cores or a.gpu_memory_gb > ram * .94:
                    raise ValueError(f"Request exceeds conservative capacity for {a.gpu_queue}")
    if a.end_step == "figs" and a.start_step in {"fig2", "fig3", "fig4", "fig5", "fig6"}:
        raise ValueError("--end-step figs requires --start-step all, fig1, or figs so supplementary dependencies exist")
    if a.end_step and STEPS.index(a.start_step) > STEPS.index(a.end_step):
        raise ValueError("--end-step cannot precede --start-step")
    if not a.conda_env.strip():
        raise ValueError("--conda-env must name an existing environment or prefix")
    project = Path(os.environ["EMS120_SCRIPT_DIR"]).resolve()
    analysis_root = os.environ["PERSISTENT_ANALYSIS_ROOT" if mode == "hpc" else "ANALYSIS_ROOT"]
    analysis = absolute(a.dir_analysis or Path(analysis_root) / "ems120")
    config = {
        "EMS120_SCRIPT_DIR": str(project), "EMS120_SHARED_DIR": absolute(os.environ["EMS120_SHARED_DIR"]),
        "EMS120_DATA_DIR": absolute(a.dir_data or Path(os.environ["DIR_DAT"]) / "phe/ems120"),
        "EMS120_ANALYSIS_DIR": analysis,
        "EMS120_BASE_MODEL_DIR": absolute(a.dir_base_model or Path(os.environ["DIR_BIG"]) / "AI_model/hfl"),
        "EMS120_MODEL_DIR": absolute(a.dir_model or Path(analysis) / "ml"),
        "EMS120_FILES_DIR": absolute(a.dir_files or os.environ["DIR_FILES"]),
        "EMS120_JOBS_DIR": absolute(a.dir_jobs or Path(analysis) / "logs"),
        "EMS120_TMP_ROOT": absolute(Path(analysis) / "tmp"), "WORKFLOW_TMP_ROOT": absolute(Path(analysis) / "tmp"),
        "EMS120_CONDA_ENV": a.conda_env, "EMS120_CONDA_EXE": a.conda_exe,
        "EMS120_RSCRIPT": a.rscript, "EMS120_PYTHON_LIB": absolute(a.python_lib) if a.python_lib else "",
        "EMS120_REPLACE": "TRUE" if a.replace else "FALSE",
        "EMS120_PREFLIGHT": "1" if a.preflight else "0", "EMS120_RUN_DX_CV": "1" if a.run_cv else "0",
        "START_STEP": a.start_step, "END_STEP": a.end_step,
        "RETICULATE_AUTOCONFIGURE": "FALSE", "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1",
        "TOKENIZERS_PARALLELISM": "false", "MPLBACKEND": "Agg", "PYTHONNOUSERSITE": "1",
    }
    for method in METHODS:
        config[f"EMS120_{method.upper()}_RUN"] = getattr(a, method + "_run")
    for key in ("CHROME_BIN", "CHROME_PATH", "EDGE_BIN", "FIG2D_PYTHON", "R_LIBS_USER", "R_LIBS_SITE", "HF_HOME",
                "EMS120_DX_BATCH_SIZE", "EMS120_GEO_BATCH_SIZE", "EMS120_PHONE_BATCH_SIZE", "EMS120_SENS_N_PERM"):
        if key in os.environ:
            config[key] = os.environ[key]
    if mode == "local" and "CUDA_VISIBLE_DEVICES" in os.environ:
        config["CUDA_VISIBLE_DEVICES"] = os.environ["CUDA_VISIBLE_DEVICES"]
    for key in ("EMS120_CONDA_ENV", "EMS120_CONDA_EXE", "EMS120_RSCRIPT", "CHROME_BIN", "CHROME_PATH", "EDGE_BIN", "FIG2D_PYTHON"):
        if config.get(key) and ("/" in config[key] or config[key].startswith("~")):
            config[key] = absolute(config[key])
    return mode, config


def build_jobs(a, mode, env):
    cpu = dict(queue=a.cpu_queue, cores=a.cpu_cores, memory_gib=a.cpu_memory_gb,
               hours=a.cpu_hours, gpus=0, memory_scope="host")
    gpu = dict(queue=a.gpu_queue, cores=a.gpu_cores, memory_gib=a.gpu_memory_gb,
               hours=a.gpu_hours, gpus=1, memory_scope="host")
    if a.gpu_model:
        gpu["gpu_model"] = a.gpu_model
    env = {**env, "EMS120_HANDOFF_DIR": str(Path(env["EMS120_JOBS_DIR"]) / "handoffs" / run_name())}
    phases = ["report"] if a.report_only else list(METHODS) if a.preflight else ["prepare", *METHODS, "report"]
    jobs = []
    for phase in phases:
        device = getattr(a, phase + "_run") if phase in METHODS else "cpu"
        settings = {**env, "EMS120_PHASE": phase, "EMS120_DEVICE": device,
                    "EMS120_REPORT_ONLY": "1" if a.report_only else "0",
                    "EMS120_CORES": str(a.gpu_cores if device == "gpu" else a.cpu_threads)}
        if device == "cpu":
            settings["CUDA_VISIBLE_DEVICES"] = ""
        key = ("preflight." if a.preflight else "") + phase
        jobs.append(Job(key=key, argv=["bash", str(Path(env["EMS120_SCRIPT_DIR"]) / "f/analyse.sh")],
                        resources=dict(gpu if device == "gpu" else cpu), env=settings,
                        cwd=env["EMS120_SCRIPT_DIR"], depends=[jobs[-1].key] if jobs else []))
    return jobs


def main(argv=None):
    p = parser()
    a = p.parse_args(argv)
    try:
        mode, env = configuration(a)
        jobs = build_jobs(a, mode, env)
        os.environ["WORKFLOW_PYTHON"] = sys.executable
        path = write_plan(jobs, env["EMS120_JOBS_DIR"], mode=mode, name="ems120")
        print("[EMS120] " + json.dumps({"mode": mode, "replace": a.replace,
              "stages": [{"key": j.key, "device": j.env["EMS120_DEVICE"], "depends": j.depends,
                          "resources": j.resources} for j in jobs],
              "data": env["EMS120_DATA_DIR"], "analysis": env["EMS120_ANALYSIS_DIR"],
              "base_model": env["EMS120_BASE_MODEL_DIR"], "model": env["EMS120_MODEL_DIR"]}, ensure_ascii=False), flush=True)
        result = run_plan(path, submit=mode, jobs=1,
                          wait=os.environ.get("SUBMIT_WAIT", "FALSE").upper() == "TRUE",
                          dry_run=a.plan_only, resume=False)
        print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
        print(f"[LOGS] {path.parent / 'logs'}", flush=True)
        return 0
    except (ValueError, KeyError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        p.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
