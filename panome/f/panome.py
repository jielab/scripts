#!/usr/bin/env python3
"""Panome 4: reference-based fixed-horizon disease risk research pipeline."""
from pathlib import Path
import argparse
import os
import sys
import shutil
from types import SimpleNamespace


def parser():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Usage examples:
  ./panome.sh --Y cvd_cad --biom prot --dry-run
  ./panome.sh --Y cvd_cad --biom prot --preflight
  ./panome.sh --Y cvd_cad --biom prot --cores 16
  ./panome.sh evaluate --run-dir /mnt/d/analysis/panome/cvd_cad/prot/v4_reference
  ./panome.sh --demo --tree hist --max-samples 1600 --run-name synthetic_check --analysis-root /tmp/panome_check
""")
    p.add_argument("command", nargs="?", choices=["run", "evaluate", "project"], default="run")
    p.add_argument("--Y", dest="trait", default="cvd_cad")
    p.add_argument("--biom", choices=["prot", "met"], default="prot")
    p.add_argument("--ukb-phe", default=os.environ.get("UKB_PHE", "/mnt/d/data/ukb/phe"))
    p.add_argument("--phe-file")
    p.add_argument("--omics-file")
    p.add_argument("--input-source", choices=["raw", "cleaned"], default="raw")
    p.add_argument("--met-input", choices=["raw", "named"], default=None)
    p.add_argument("--met-map")
    p.add_argument("--id-col", default="eid")
    p.add_argument("--baseline-col", default="date_attend")
    p.add_argument("--diagnosis-col")
    p.add_argument("--death-col", default="date_death")
    p.add_argument("--lost-col", default="date_lost")
    p.add_argument("--end-date", default="2023-04-01")
    p.add_argument("--disease-evidence-col", default="")
    p.add_argument("--healthy-date-cols", default="")
    p.add_argument("--covariates", default="age,sex,tdi,PC1,PC2,center")
    p.add_argument("--residualize", default=None)
    p.add_argument("--categorical", default=None)
    p.add_argument("--group-col", default="")
    p.add_argument("--split-file", default="")
    p.add_argument("--module-file", default="")
    p.add_argument("--transform", choices=["none", "log1p"], default=None)
    p.add_argument("--exclude-features", default="")
    p.add_argument("--feature-missing", type=float, default=.2)
    p.add_argument("--sample-missing", type=float, default=.2)
    p.add_argument("--horizon", type=float, default=10)
    p.add_argument("--min-censor-survival", type=float, default=.05)
    p.add_argument("--min-events", type=int, default=20)
    p.add_argument("--panel-size", type=int, default=100)
    p.add_argument("--panel-sizes", default="100,300,1000")
    p.add_argument("--match-k", default="5,10,20")
    p.add_argument("--all-k", default="30,100,300")
    p.add_argument("--donor-neighbors", type=int, default=100)
    p.add_argument("--prior-strength", type=float, default=10)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--oof-repeats", type=int, default=3)
    p.add_argument("--teacher-c", type=float, default=.01)
    p.add_argument("--quality-teacher", choices=["elasticnet", "ensemble"], default="ensemble")
    p.add_argument("--quality-trees", type=int, default=100)
    p.add_argument("--c-grid", default="0.001,0.01,0.1,1")
    p.add_argument("--max-iter", type=int, default=3000)
    p.add_argument("--fit-stability", type=float, default=2/3)
    p.add_argument("--min-gain", type=float, default=0.)
    p.add_argument("--dimensions", type=int, default=32)
    p.add_argument("--accept-quantile", type=float, default=.95)
    p.add_argument("--min-match-ess", type=float, default=3.)
    p.add_argument("--min-coverage", type=float, default=.5)
    p.add_argument("--tree", choices=["lightgbm", "hist", "none"], default="lightgbm")
    p.add_argument("--tree-estimators", type=int, default=300)
    p.add_argument("--bootstrap", type=int, default=200)
    p.add_argument("--pair-caliper", type=float, default=.01)
    p.add_argument("--max-pairs", type=int, default=100)
    p.add_argument("--explanation-samples", type=int, default=1000)
    p.add_argument("--mask-fraction", type=float, default=.1)
    p.add_argument("--mask-repeats", type=int, default=3)
    p.add_argument("--mask-models", default="panome,copy1_fullproteome,random_panel,diversity_panel")
    p.add_argument("--analysis-root", default="/mnt/d/analysis/panome")
    p.add_argument("--run-name", default="v4_reference")
    p.add_argument("--run-dir")
    p.add_argument("--output", help="CSV output for project")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--cores", type=int, default=4)
    p.add_argument("--r-bin", default="Rscript")
    p.add_argument("--max-samples", type=int, default=0)
    p.add_argument("--demo-features", type=int, default=80)
    p.add_argument("--demo", action="store_true")
    p.add_argument("--train-only", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--preflight", action="store_true")
    p.add_argument("--full-input-hash", action="store_true")
    p.add_argument("--replace", action="store_true")
    return p


def configure(a):
    for field, cast in [("panel_sizes", int), ("match_k", int), ("all_k", int), ("c_grid", float)]:
        setattr(a, field, sorted(set(cast(v) for v in getattr(a, field).split(","))))
        if not getattr(a, field) or min(getattr(a, field)) <= 0:
            raise ValueError(f"Invalid {field}")
    a.panel_sizes = sorted(set(a.panel_sizes+[a.panel_size]))
    a.outcome_type, a.target_col = "survival", a.trait
    a.diagnosis_col = a.diagnosis_col or "fod_icd10_"+a.trait
    base = Path(a.ukb_phe)
    a.phe_file = a.phe_file or str(base/"Rdata/all.rds")
    a.omics_file = a.omics_file or str(base/(f"Rdata/{a.biom}.rds" if a.input_source == "cleaned" else
                                           "rap/raw/prot.tab.gz" if a.biom == "prot" else "rap/met.tab.gz"))
    a.met_map = a.met_map or str(base/"common/met.lst")
    a.met_input = a.met_input or ("raw" if a.biom == "met" and a.input_source == "raw" else "named")
    a.residualize = f"age,sex,{a.biom}.plate" if a.residualize is None else a.residualize
    a.categorical = f"sex,center,{a.biom}.plate" if a.categorical is None else a.categorical
    a.transform = a.transform or ("log1p" if a.biom == "met" and a.input_source == "raw" and not a.demo else "none")
    if a.horizon <= 0 or a.panel_size < 3 or a.dimensions < 1 or a.folds < 2 or a.oof_repeats < 2:
        raise ValueError("Invalid horizon, panel size, dimensions, folds, or OOF repeats (minimum 2)")
    if not 0 < a.accept_quantile <= 1 or not 0 < a.fit_stability <= 1:
        raise ValueError("Invalid acceptance/stability quantile")
    if not 0 < a.mask_fraction < 1 or a.mask_repeats < 1 or a.explanation_samples < 0:
        raise ValueError("Invalid masking sensitivity settings")
    for field in ["feature_missing", "sample_missing", "min_censor_survival"]:
        if not 0 < getattr(a, field) < 1:
            raise ValueError(f"Invalid {field}")
    if any(v in a.run_name for v in ["/", "\\"]) or a.run_name in ["", ".", ".."]:
        raise ValueError("Use a simple nonempty --run-name")
    return a


def main():
    a = configure(parser().parse_args())
    if a.cores < 1:
        raise ValueError("--cores must be positive")
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        os.environ[name] = str(a.cores)
    # Imports occur after thread limits have been set.
    from common import log, run_lock, dump
    from pipeline import train, evaluate, project
    from threadpoolctl import threadpool_limits
    out = Path(a.run_dir) if a.run_dir else Path(a.analysis_root)/a.trait/a.biom/a.run_name
    if a.command == "project":
        if not a.run_dir or not a.output:
            raise ValueError("project requires --run-dir and --output")
        with threadpool_limits(a.cores):
            project(out, a.phe_file, a.omics_file, a.output, a.r_bin, a.met_input)
        return
    if a.dry_run:
        import json
        print(json.dumps(dict(output=str(out), config=vars(a)), indent=2))
        return
    if a.preflight:
        import importlib.util
        paths = [] if a.demo else [a.phe_file, a.omics_file]
        if a.module_file:
            paths.append(a.module_file)
        if a.split_file:
            paths.append(a.split_file)
        if a.biom == "met" and a.met_input == "raw" and not a.demo:
            paths.append(a.met_map)
        missing = [v for v in paths if not Path(v).is_file()]
        if missing:
            raise FileNotFoundError("Missing inputs: "+", ".join(missing))
        if any(str(v).lower().endswith(".rds") for v in paths) and not importlib.util.find_spec("pyreadr"):
            raise RuntimeError("RDS needs pyreadr for lossless binary reading; install requirements.txt")
        if a.tree == "lightgbm" and not importlib.util.find_spec("lightgbm"):
            raise RuntimeError("LightGBM missing; install requirements.txt or explicitly choose --tree hist")
        log("DONE", "preflight", "paths and required reader/tree dependencies available; data values not yet read")
        return
    if a.command == "evaluate":
        if not (out/"MODEL_FROZEN.json").is_file():
            raise ValueError("No frozen model found")
        with run_lock(out), threadpool_limits(a.cores):
            evaluate(out)
        return
    if a.tree == "lightgbm":
        import importlib.util
        if not importlib.util.find_spec("lightgbm"):
            raise RuntimeError("LightGBM missing; install requirements.txt or explicitly choose --tree hist before training")
    if out.exists() and any(out.iterdir()):
        if not a.replace:
            raise FileExistsError(f"Existing run: {out}. Use evaluate, a new run-name, or --replace")
        if not (out/"manifest.json").is_file() or (out/".lock").exists():
            raise ValueError("Refusing replacement without Panome manifest or while locked")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    with run_lock(out), threadpool_limits(a.cores):
        train(a, out)
        if not a.train_only:
            log("START", "test_evaluation")
            evaluate(out)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"[PANOME] ERROR {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        raise
