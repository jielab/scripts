#!/usr/bin/env python3
"""Lightweight EMS120 planner. No Torch, Conda, R, or patient data on login nodes."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

from workflow_jobs import Job, run_plan, write_plan

STEPS = ("all", "fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "figs")


def boolean(value):
    value = str(value).lower()
    if value in {"1", "true", "yes", "on"}:
        return True
    if value in {"0", "false", "no", "off"}:
        return False
    raise argparse.ArgumentTypeError("expected TRUE or FALSE")


def absolute(value):
    value = os.fspath(value)
    if not value or any(ord(c) < 32 or ord(c) == 127 for c in value):
        raise ValueError("A path must be nonempty and contain no control characters")
    return str(Path(value).expanduser().resolve())


def parser():
    p = argparse.ArgumentParser(
        description="EMS120 at scripts/ems120: default QiMing HPC, explicit local execution.",
        epilog="Common root options are resolved by ../0f/0os.sh. No Conda activation or GPU initialization occurs while planning.",
    )
    p.add_argument("--plan-only", "--dry-run", action="store_true", help="write plan/.cmd without submission or analysis")
    p.add_argument("--check-gpu", action="store_true", help="run only R/reticulate + GPU environment checks, inside an LSF GPU job by default")
    p.add_argument("--queue", default=os.getenv("EMS120_QUEUE", "4a100-40"))
    default_cores = min(8, os.cpu_count() or 1) if os.getenv("submit") == "local" else 8
    p.add_argument("--cores", "--threads", type=int, default=os.getenv("EMS120_CORES", str(default_cores)))
    p.add_argument("--memory-gb", type=float, default=os.getenv("EMS120_MEMORY_GB", "96" if os.getenv("submit") != "local" else "32"), help="total host RAM in GiB; not GPU VRAM (local is advisory)")
    p.add_argument("--hours", type=float, default=os.getenv("EMS120_HOURS"), help="wall limit: 120 hours for analysis, 1 hour for --check-gpu")
    p.add_argument("--gpus", type=int, default=1, help="exactly one GPU; this pipeline is not distributed")
    p.add_argument("--gpu", type=boolean, default=True, help="TRUE by default; FALSE is allowed only with --submit local")
    p.add_argument("--phone-gpu", choices=("required", "auto", "off"), default=os.getenv("EMS120_PHONE_GPU", "required"), help="XGBoost candidate GPU policy; cached CPU winners remain valid")
    p.add_argument("--dir-data", default=os.getenv("EMS120_DATA_DIR"), help="EMS input directory, containing clean/ and szu_xia/")
    p.add_argument("--dir-analysis", "--analysis-dir", default=os.getenv("EMS120_ANALYSIS_DIR"))
    p.add_argument("--dir-base-model", default=os.getenv("EMS120_BASE_MODEL_DIR"), help="local pretrained HFL/MacBERT model, no automatic download")
    p.add_argument("--dir-model", default=os.getenv("EMS120_MODEL_DIR"), help="trained-model root; default <analysis>/ml")
    p.add_argument("--dir-files", default=os.getenv("EMS120_FILES_DIR"))
    p.add_argument("--dir-jobs", default=os.getenv("EMS120_JOBS_DIR"), help="persistent plans/logs/status; default <analysis>/.jobs")
    p.add_argument("--conda-env", default=os.getenv("EMS120_CONDA_ENV", "ai"), help="existing Conda env name or absolute prefix")
    p.add_argument("--conda-exe", default=os.getenv("EMS120_CONDA_EXE", os.getenv("CONDA_EXE", "")))
    p.add_argument("--rscript", default=os.getenv("EMS120_RSCRIPT", ""), help="optional Rscript executable; default from activated environment")
    p.add_argument("--start-step", type=str.lower, choices=STEPS, default="all")
    p.add_argument("--end-step", type=str.lower, choices=STEPS[1:], default="")
    p.add_argument("--run-cv", type=boolean, default=True)
    # These flags have already been consumed by 0os.sh. Including their help
    # here documents the public shell interface without a second path table.
    g = p.add_argument_group("shared options (use ems120.sh)")
    g.add_argument("--submit", choices=("hpc", "local"), help="default hpc; missing bsub never triggers local execution")
    for flag in ("--dir-dat", "--dir-big", "--dir-tmp", "--date-scratch"):
        g.add_argument(flag, help="shared 0os.sh override")
    g.add_argument("--wait", nargs="?", const="TRUE", help="wait for LSF completion; default FALSE")
    return p


def configuration(a):
    mode = os.environ.get("submit", "hpc")
    if mode not in {"hpc", "local"}:
        raise ValueError("--submit must be hpc or local")
    if a.submit or a.dir_dat or a.dir_big or a.dir_tmp or a.date_scratch or a.wait:
        raise ValueError("Use ems120.sh so shared options are resolved through 0os.sh")
    if a.gpus != 1:
        raise ValueError("EMS120 supports --gpus 1. Use --submit local --gpu FALSE for an explicit local CPU run")
    if not a.gpu and mode == "hpc":
        raise ValueError("HPC EMS120 requires a GPU allocation; --gpu FALSE is local-only")
    if not a.gpu and a.check_gpu:
        raise ValueError("--check-gpu requires --gpu TRUE")
    if a.cores < 1 or a.cores > 192:
        raise ValueError("--cores must be in 1..192")
    if not math.isfinite(a.memory_gb) or a.memory_gb <= 0:
        raise ValueError("--memory-gb must be positive and finite")
    if a.hours is None:
        a.hours = 1.0 if a.check_gpu else 120.0
    if not math.isfinite(a.hours) or not 0 < a.hours <= 336:
        raise ValueError("--hours must be in (0, 336] for the supplied QiMing GPU queues")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}", a.queue):
        raise ValueError("Invalid queue name")
    # Use the smaller node variant for mixed GPU queues from the supplied table.
    limits = {"4a100-40": (56, 384), "4a100-80": (32, 128), "2a100-40": (40, 192), "2a100-80": (40, 192),
              "2h100-80": (192, 768), "hgx": (64, 1024), "8v100-32-sc": (48, 384),
              "4v100-16-e5": (24, 64), "4v100-16-sc": (40, 192),
              "2v100-16-sc": (40, 384), "2v100-32-e5": (24, 256), "2v100-32-sc": (40, 192)}
    if mode == "hpc" and a.queue in {"33", "38", "ot38", "2t50c", "1t75c", "1t88c", "1t81c",
                                     "73x", "52", "63", "v4-3t", "sk-6t-104", "75m", "v3-6t", "v3-64", "spec-cpu"}:
        raise ValueError("Select a GPU queue for EMS120")
    if mode == "hpc" and a.queue in limits:
        cpu, ram = limits[a.queue]
        if a.cores > cpu or a.memory_gb > ram * 0.94:
            raise ValueError(f"Request exceeds conservative capacity for {a.queue}: <= {cpu} cores, <= {ram * .94:g} GiB RAM. Select another GPU queue or lower the request")
    if a.end_step and STEPS.index(a.start_step) > STEPS.index(a.end_step):
        raise ValueError("--end-step cannot precede --start-step")
    if not a.conda_env.strip():
        raise ValueError("--conda-env must name an existing environment or prefix")

    project = Path(os.environ["EMS120_SCRIPT_DIR"]).resolve()
    analysis_root = os.environ["PERSISTENT_ANALYSIS_ROOT" if mode == "hpc" else "ANALYSIS_ROOT"]
    analysis = absolute(a.dir_analysis or Path(analysis_root) / "ems120")
    config = {
        "EMS120_SCRIPT_DIR": str(project),
        "EMS120_SHARED_DIR": absolute(os.environ["EMS120_SHARED_DIR"]),
        "EMS120_DATA_DIR": absolute(a.dir_data or Path(os.environ["DIR_DAT"]) / "ems120"),
        "EMS120_ANALYSIS_DIR": analysis,
        "EMS120_BASE_MODEL_DIR": absolute(a.dir_base_model or Path(os.environ["DIR_BIG"]) / "AI_model" / "hfl"),
        "EMS120_MODEL_DIR": absolute(a.dir_model or Path(analysis) / "ml"),
        "EMS120_FILES_DIR": absolute(a.dir_files or os.environ["DIR_FILES"]),
        "EMS120_JOBS_DIR": absolute(a.dir_jobs or Path(analysis) / ".jobs"),
        "EMS120_TMP_ROOT": absolute(Path(analysis) / "tmp"),
        "EMS120_CONDA_ENV": a.conda_env,
        "EMS120_CONDA_EXE": a.conda_exe,
        "EMS120_RSCRIPT": a.rscript,
        "EMS120_REQUIRE_GPU": "1" if a.gpu else "0",
        "EMS120_PHONE_GPU": a.phone_gpu if a.gpu else "off",
        "EMS120_CHECK_GPU_ONLY": "1" if a.check_gpu else "0",
        "EMS120_RUN_DX_CV": "1" if a.run_cv else "0",
        "EMS120_CORES": str(a.cores),
        "START_STEP": a.start_step,
        "END_STEP": a.end_step,
        "RETICULATE_AUTOCONFIGURE": "FALSE",
        "HF_HUB_OFFLINE": "1", "TRANSFORMERS_OFFLINE": "1",
        "TOKENIZERS_PARALLELISM": "false",
        "MPLBACKEND": "Agg", "PYTHONNOUSERSITE": "1",
    }
    # Preserve selected runtime/render settings, not the submit host's CUDA
    # visibility, active Conda prefix, or full ambient environment.
    names = ("CHROME_BIN", "CHROME_PATH", "EDGE_BIN", "FIG2D_PYTHON", "R_LIBS_USER", "R_LIBS_SITE", "HF_HOME",
             "EMS120_DX_BATCH_SIZE", "EMS120_GEO_BATCH_SIZE", "EMS120_PHONE_BATCH_SIZE", "EMS120_SENS_N_PERM")
    for key in names:
        if key in os.environ:
            config[key] = os.environ[key]
    if mode == "local" and "CUDA_VISIBLE_DEVICES" in os.environ:
        config["CUDA_VISIBLE_DEVICES"] = os.environ["CUDA_VISIBLE_DEVICES"]
    # Freeze an explicit executable/prefix path before worker cwd changes.
    for key in ("EMS120_CONDA_ENV", "EMS120_CONDA_EXE", "EMS120_RSCRIPT", "CHROME_BIN", "CHROME_PATH", "EDGE_BIN", "FIG2D_PYTHON"):
        if config.get(key) and ("/" in config[key] or config[key].startswith("~")):
            config[key] = absolute(config[key])
    return mode, config


def main(argv=None):
    p = parser()
    a = p.parse_args(argv)
    try:
        mode, env = configuration(a)
        resources = {"queue": a.queue, "cores": a.cores, "memory_gib": a.memory_gb,
                     "hours": a.hours, "gpus": 1 if a.gpu else 0, "memory_scope": "host"}
        os.environ["WORKFLOW_PYTHON"] = sys.executable
        job = Job(
            key="check_gpu" if a.check_gpu else "analysis",
            argv=["bash", str(Path(env["EMS120_SCRIPT_DIR"]) / "f" / "run_worker.sh")],
            resources=resources, env=env, cwd=env["EMS120_SCRIPT_DIR"],
        )
        path = write_plan([job], env["EMS120_JOBS_DIR"], mode=mode, name="ems120")
        print("[EMS120] " + json.dumps({"mode": mode, "resources": resources,
              "data": env["EMS120_DATA_DIR"], "analysis": env["EMS120_ANALYSIS_DIR"],
              "base_model": env["EMS120_BASE_MODEL_DIR"], "model": env["EMS120_MODEL_DIR"],
              "phone_gpu": env["EMS120_PHONE_GPU"]}, ensure_ascii=False), flush=True)
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
