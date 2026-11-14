#!/usr/bin/env python3
"""Plan LE8 jobs without loading participant data or scientific runtimes.

All methods default to CPU. A selected deep model can explicitly use one GPU;
its input preparation and common evaluation remain CPU jobs in the same DAG.
"""
from __future__ import annotations
import argparse
import json
import math
import os
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT.parent / "0f"))
from workflow_paths import configured_path
from workflow_jobs import Job, write_plan, run_plan, resolve_mode
from outputs import run_name

DEFAULT = "c1_correlate,c2_cause,c3_coloc,c4_connect,c5_cellulation"
MODULES = {"c1_correlate", "c1_reliability", "c2_cause", "c3_coloc", "c4_connect",
           "c5_cellulation", "pgs_focus", "final", "share", "shiny", "install"}
METHODS = ("glmnet", "lightgbm", "abm", "abm_deep")
REMOVED = {"--schedule", "--allocation", "--device", "--abm-run", "--abm-count-run",
           "--gpu-abm-deep", "--gpu-coloc", "--run-gpu-coloc", "--coloc-validation",
           "--history-config", "--eigh-backend", "--skip-abm", "--run-abm", "--resume-plan", "--dry-run", "--preflight",
           "--dir-analysis-root", "--dir-ukb-phe", "--trait", "-Y", "-b", "--min_cases_with_biom",
           "--abm-gpu-config", "--abm-host-args",
           "--gpu-jobs", "--group-column", "--Y-date", "--vars.adj",
           "--gpu-coloc-bin",
           "--replace-association", "--replace-models", "--replace-evaluation", "--replace-reliability"}


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument("--Y", default=os.getenv("Y", "cvd_cad,ra"))
    p.add_argument("--biom", default=os.getenv("BIOM", "prot,met"))
    p.add_argument("--analysis-root", type=Path,
                   default=Path(os.getenv("LE8_ANALYSIS_ROOT", configured_path("ANALYSIS_ROOT", "le8"))))
    p.add_argument("--cores", type=int, default=int(os.getenv("N_CORES", "8")))
    p.add_argument("--memory-limit-gb", type=float, default=float(os.getenv("LE8_MEMORY_LIMIT_GB", "32")))
    p.add_argument("--queue", default=os.getenv("LE8_QUEUE", "33"))
    p.add_argument("--gpu-queue", default=os.getenv("LE8_GPU_QUEUE", os.getenv("HPC_GPU_QUEUE", "2v100-32-e5")))
    p.add_argument("--gpu-cores", type=int, default=int(os.getenv("LE8_GPU_CORES", "8")))
    p.add_argument("--gpu-memory-limit-gb", type=float, default=float(os.getenv("HPC_GPU_MEMORY_GIB", "160")))
    p.add_argument("--hours", type=float, default=float(os.getenv("LE8_HOURS", "72")))
    p.add_argument("--abm-deep-run", choices=("cpu", "gpu", "auto"), default=None,
                   help="Execution of the selected abm_deep; auto uses its CPU default")
    p.add_argument("--methods", help="Comma-separated glmnet,lightgbm,abm,abm_deep")
    p.add_argument("--config", type=Path, default=ROOT / "config/c1_methods.json")
    p.add_argument("--replace", choices=("TRUE", "FALSE"), default="TRUE")
    p.add_argument("--python", default=os.getenv("LE8_REPORT_PYTHON", os.getenv("PYTHON_BIN", "python3")))
    p.add_argument("--r-bin", default=os.getenv("R_BIN", "Rscript"))
    p.add_argument("--plan-dir", type=Path)
    p.add_argument("--plan-only", action="store_true")
    for flag in ("fit-joint", "fit-genetic", "fit-reference", "details", "index-only"):
        p.add_argument("--" + flag, action="store_true")
    p.add_argument("--summary-worker", action="store_true", help=argparse.SUPPRESS)
    return p


def parse_cli(argv=None):
    tokens = list(sys.argv[1:] if argv is None else argv)
    modules = tokens.pop(0) if tokens and not tokens[0].startswith("-") else DEFAULT
    removed = {t.split("=", 1)[0] for t in tokens} & REMOVED
    if removed:
        raise ValueError("Removed option(s): " + ", ".join(sorted(removed)) +
                         ". Select methods with --methods; only abm_deep accepts --abm-deep-run cpu|gpu|auto.")
    a, extra = parser().parse_known_args(tokens)
    a.modules = modules
    a.deep_run_explicit = a.abm_deep_run is not None
    a.abm_deep_run = "cpu" if a.abm_deep_run in (None, "auto") else a.abm_deep_run
    return a, extra


def selected_modules(a):
    result = list(dict.fromkeys(m.strip() for m in a.modules.split(",")))
    if not result or set(result) - MODULES:
        raise ValueError("Unknown modules: " + ", ".join(sorted(set(result) - MODULES)))
    return result


def selected_methods(a):
    raw = a.methods.split(",") if a.methods else json.loads(a.config.read_text()).get("methods", list(METHODS[:3]))
    result = list(dict.fromkeys(str(m).strip() for m in raw))
    if not result or set(result) - set(METHODS):
        raise ValueError("C1 methods must be selected from " + ",".join(METHODS))
    return result


def build_jobs(a, extra=()):
    modules = selected_modules(a)
    traits, layers = list(dict.fromkeys(a.Y.split(","))), list(dict.fromkeys(a.biom.split(",")))
    if not traits or any(not re.fullmatch(r"[A-Za-z0-9_]+", t) for t in traits):
        raise ValueError("Invalid outcome list")
    if not layers or not set(layers) <= {"prot", "met"}:
        raise ValueError("Biomarker layers must be prot and/or met")
    if any(not math.isfinite(v) or v <= 0 for v in (a.cores, a.gpu_cores, a.memory_limit_gb, a.gpu_memory_limit_gb, a.hours)):
        raise ValueError("Resource requests must be positive")
    if not all(float(v).is_integer() for v in (a.memory_limit_gb,a.gpu_memory_limit_gb)):
        raise ValueError("--memory-limit-gb must be a whole number of GiB")
    if a.methods and "c1_correlate" not in modules:
        raise ValueError("--methods requires c1_correlate")
    methods = selected_methods(a) if "c1_correlate" in modules else []
    deep = "abm_deep" in methods
    if a.deep_run_explicit and not deep:
        raise ValueError("--abm-deep-run requires selecting abm_deep with --methods or --config")
    primary_methods = [m for m in methods if m != "abm_deep"]
    deep_gpu = deep and a.abm_deep_run == "gpu"
    root = a.analysis_root.resolve()
    env = {k:v for k,v in os.environ.items() if
           k.startswith(("LE8_", "C1_", "C2_", "C3_", "C4_", "C5_", "FINAL_", "MRLINK2_", "PGS_"))
           or k in {"DIR_WORK", "DIR_DATA", "DIR_SCRATCH", "DIR_SOFTWARE", "ANALYSIS_ROOT", "DATE_SCRATCH",
                    "UKB_PHE", "UKB_COMMON", "UKB_BGEN_DIR", "PHE_F", "CONDA_SH", "CONDA_EXE", "PATH", "LD_LIBRARY_PATH"}}
    env.update(MRLINK2_EIGH_BACKEND="cpu", LE8_GPU_JOBS="1", PYTHON_BIN=a.python,
               ABM_PYTHON=a.python, LE8_REPORT_PYTHON=a.python, R_BIN=a.r_bin,
               LE8_RUNTIME_READY="1", LE8_COHORT_CACHE_ROOT=str(root / "tmp/cohort_cache"),
               WORKFLOW_TMP_ROOT=str(root / "tmp/workflows"), SCRIPT_LOG_DIR=str(root / "logs/runs"),
               SCRIPT_FOREGROUND="1", SCRIPT_CONSOLE_ACTIVE="1", WORKFLOW_DEVICE="cpu")
    def resources(gpu=False):
        return dict(queue=a.gpu_queue if gpu else a.queue, cores=a.gpu_cores if gpu else a.cores,
                    memory_gib=min(a.memory_limit_gb,a.gpu_memory_limit_gb) if gpu else a.memory_limit_gb,
                    hours=a.hours,gpus=int(gpu))
    def stage_env(gpu=False):
        r=resources(gpu)
        return dict(env, LE8_TOTAL_CORES=str(r["cores"]), LE8_RESOURCE_MEMORY_GIB=str(r["memory_gib"]),
                    LE8_MEMORY_LIMIT_GB=str(r["memory_gib"]),WORKFLOW_DEVICE="cuda" if gpu else "cpu")
    if "install" in modules:
        if modules != ["install"]:
            raise ValueError("Run installation separately from analysis")
        return [Job("install", ["bash",str(ROOT/"install.sh"),"--python",a.python,"--r-bin",a.r_bin,*extra],
                    [],resources(),stage_env(),str(ROOT))]
    batch=run_name()
    scope_keys=[f"{t}.{b}" for t in traits for b in layers]
    jobs=[]
    for trait in traits:
        for layer in layers:
            scope=f"{trait}.{layer}"
            handoff=str(root/trait/layer/"logs/c1_stages"/batch)
            deep_root=root/"exploratory/abm_deep"
            parts=(["primary"] if primary_methods else [])+(["deep"] if deep else [])
            def step(module,key,tail=(),chosen=None,gpu=False,part=None,output=None,aggregate=False):
                r=resources(gpu); e=stage_env(gpu)
                destination=output or root
                e.update(LE8_C1_BATCH_DIR=handoff, LE8_C1_MODEL_PARTS=",".join(parts),
                         LE8_ANALYSIS_ROOT=str(destination), LE8_PRIMARY_ANALYSIS_ROOT=str(root),
                         LE8_CACHE_ROOT=str(destination/trait/layer/"tmp/cache"),
                         SCRIPT_LOG_DIR=str(destination/trait/layer/"logs/runs"))
                if part: e["LE8_C1_MODEL_PART"]=part
                argv=["bash",str(ROOT/"le8.sh"),module,"--analysis-root",str(destination),
                      "--Y",trait,"--biom",layer,"--cores",str(r["cores"]),
                      "--memory-limit-gb",f"{r['memory_gib']:g}","--r-bin",a.r_bin,
                      "--replace",a.replace,"--config",str(a.config.resolve()),*map(str,extra),*tail]
                if module=="c1_correlate":
                    argv += ["--methods",",".join(methods if chosen is None else chosen),
                             "--device","cuda" if gpu else "cpu"]
                return dict(key=f"{scope}.{key}",argv=argv,env=e,aggregate=aggregate)
            def pack(phase,steps,deps=(),gpu=False,output=None):
                key=f"{scope}.{phase}"
                argv=[a.python,str(ROOT/"f/run_traits.py"),"--analysis-root",str(output or root),
                      "--scope",scope,"--phase",phase,"--batch-id",batch,"--scopes",",".join(scope_keys),
                      "--steps-json",json.dumps(steps,separators=(",",":"))]
                jobs.append(Job(key,argv,list(deps),resources(gpu),stage_env(gpu),str(ROOT),"le8."+key))
                return key
            prepare=[]; work=[]; finish=[]
            if methods:
                prepare.append(step("c1_correlate","c1.prepare",["--c1-stage","prepare"]))
                work.append(step("c1_correlate","c1.association",["--c1-stage","association"]))
                if primary_methods:
                    work.append(step("c1_correlate","c1.primary.models",["--c1-stage","models"],
                                     chosen=primary_methods,part="primary"))
                if deep and not deep_gpu:
                    work.append(step("c1_correlate","c1.deep.models",["--c1-stage","models"],
                                     chosen=["abm_deep"],part="deep"))
                finish.append(step("c1_correlate","c1.evaluation",["--c1-stage","evaluation"]))
            if "c1_reliability" in modules or "glmnet" in primary_methods:
                finish.append(step("c1_reliability","c1.reliability"))
            c1_work_count=len(work)
            c1_finish_count=len(finish)
            for module in ("c2_cause","c3_coloc","c4_connect","c5_cellulation","pgs_focus"):
                if module not in modules: continue
                stages=("c4_connect","c4_panel_validation","c4_explain") if module=="c4_connect" else (module,)
                for phase in stages:
                    item=step(module,phase)
                    if module=="c4_connect": item["env"]["LE8_C4_STAGE"]=phase
                    work.append(item)
            if methods or "c1_reliability" in modules:
                finish.append(dict(key=f"{scope}.c1.summary",argv=[a.python,str(Path(__file__).resolve()),
                    "--summary-worker","--analysis-root",str(root),"--Y",a.Y,"--biom",a.biom],
                    env=stage_env(),aggregate=True))
            for module in ("final","share","shiny"):
                if module in modules:
                    tail=["--"+f for f in ("fit-joint","fit-genetic","fit-reference","details","index-only")
                          if getattr(a,f.replace("-","_"))] if module=="final" else []
                    if module=="shiny" and resolve_mode()=="hpc": tail += ["--prepare-only"]
                    finish.append(step(module,module,tail,aggregate=True))
            if deep_gpu:
                prepared=pack("prepare",prepare)
                primary=pack("primary",work,[prepared]) if work else prepared
                deep_job=pack("abm_deep",[step("c1_correlate","c1.deep.models",["--c1-stage","models","--association","skip"],
                    chosen=["abm_deep"],gpu=True,part="deep",output=deep_root)], [prepared],True,deep_root)
                pack("evaluation",finish,list(dict.fromkeys([primary,deep_job])))
            else:
                pack("analysis",prepare+work[:c1_work_count]+finish[:c1_finish_count]+work[c1_work_count:]+finish[c1_finish_count:])
    return jobs


def main(argv=None):
    a,extra=parse_cli(argv)
    if a.summary_worker:
        if resolve_mode()=="hpc" and not os.getenv("LSB_JOBID"):
            raise RuntimeError("Summary worker requires an LSF allocation")
        from c1_report import write_reliability_cross_trait_summary,write_pgs_cross_trait_summary
        for layer in a.biom.split(","):
            write_reliability_cross_trait_summary(a.analysis_root,layer)
            write_pgs_cross_trait_summary(a.analysis_root,layer)
        return 0
    jobs=build_jobs(a,extra);mode=resolve_mode()
    for job in jobs:
        job.job_name = "le8." + job.key.removesuffix(".analysis")
    plan=write_plan(jobs,a.plan_dir or a.analysis_root/"logs/workflows",mode=mode,name="le8")
    print(f"[LE8] {len(jobs)} jobs; {sum(j.resources['gpus'] for j in jobs)} explicitly selected GPU allocations; plan={plan}",flush=True)
    if a.plan_only:
        print("[LE8] Plan only; no analysis or submission",flush=True)
        return 0
    result=run_plan(plan,submit=mode,jobs=int(os.getenv("jobs_N","1")),wait=os.getenv("SUBMIT_WAIT","FALSE")=="TRUE")
    print(f"[LE8] {result.get('status', 'FINISHED')}; plan: {plan}", flush=True)
    return int(result.get("exit_code",0))


if __name__=="__main__":
    raise SystemExit(main())
