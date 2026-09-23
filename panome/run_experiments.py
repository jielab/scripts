#!/usr/bin/env python3
"""Repeated training on a fixed outer split; evaluate only after all fits freeze.

Example: python run_experiments.py --output-root /mnt/d/analysis/panome_suite \
  --seeds 2026,2027,2028 --residualization-sensitivity -- --Y cvd_cad --biom prot --device cuda
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import pandas as pd
import numpy as np


def collect(paths, output):
    scores, agreement = [], []
    for run in paths:
        frame = pd.read_csv(run/"test_metrics.csv")
        frame["run"] = run.name
        frame["regime"] = run.name.split("_",2)[-1]
        scores.append(frame)
    frame = pd.concat(scores,ignore_index=True)
    frame.to_csv(output/"all_run_metrics.csv",index=False)
    frame.groupby(["regime","model","subset"])[["AUC_IPCW","Brier_IPCW","AUPRC_IPCW","coverage"]].agg(
        ["mean","std","min","max"]).to_csv(output/"repeat_summary.csv")
    for i,left in enumerate(paths):
        for right in paths[i+1:]:
            a = pd.read_csv(left/"panel_panome_transformer.csv",dtype=str).iloc[:,0]
            b = pd.read_csv(right/"panel_panome_transformer.csv",dtype=str).iloc[:,0]
            sa,sb = set(a),set(b)
            pa = pd.read_csv(left/"test_individuals.csv",dtype={a.name:str})
            pb = pd.read_csv(right/"test_individuals.csv",dtype={a.name:str})
            joined = pa[[a.name,"panome_transformer"]].merge(pb[[a.name,"panome_transformer"]],on=a.name)
            corr = joined[["panome_transformer_x","panome_transformer_y"]].corr().iloc[0,1] if len(joined)>2 else np.nan
            agreement.append(dict(run_A=left.name,run_B=right.name,panel_Jaccard=len(sa&sb)/len(sa|sb),
                common_test_people=len(joined),risk_correlation=corr))
    pd.DataFrame(agreement).to_csv(output/"reference_and_risk_stability.csv",index=False)
    (output/"INTERPRETATION.txt").write_text(
        "Repeated runs share people and are not independent cohorts. SD is training/split sensitivity, not a standard error. "
        "All model choices were frozen before test evaluation. Choosing a best approach from these test results requires external confirmation.\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root",required=True)
    parser.add_argument("--seeds",default="2026,2027,2028")
    parser.add_argument("--residualization-sensitivity",action="store_true")
    parser.add_argument("--vary-splits",action="store_true")
    parser.add_argument("--include-null",action="store_true")
    parser.add_argument("--dry-run",action="store_true")
    parser.add_argument("args",nargs=argparse.REMAINDER)
    args = parser.parse_args(); rest = args.args[1:] if args.args[:1]==["--"] else args.args
    prohibited = ["--run-dir","--run-name","--seed","--analysis-root","--split-file","--train-only","--replace","--resume"]
    if any(v in rest for v in prohibited):
        raise ValueError("Suite controls run paths, seed, split and training phases; remove conflicting options")
    seeds = [int(v) for v in args.seeds.split(",")]
    if len(set(seeds))!=len(seeds): raise ValueError("Duplicate seeds")
    if len(seeds)>1 and ("--demo" in rest or ("--max-samples" in rest and rest[rest.index("--max-samples")+1]!="0")):
        raise ValueError("Repeated fits require one fixed input cohort. Export synthetic/subsampled inputs once before using multiple seeds.")
    output = Path(args.output_root).resolve(); output.mkdir(parents=True,exist_ok=True)
    root = Path(__file__).resolve().parent
    biom = rest[rest.index("--biom")+1] if "--biom" in rest else "prot"
    modes = ["technical_only","age_sex_residualized"] if args.residualization_sensitivity else ["base"]
    if args.include_null: modes += ["shuffled_outcomes"]
    if args.residualization_sensitivity and "--residualize" in rest:
        raise ValueError("Suite sensitivity controls --residualize")
    commands, runs = [], []
    first = output/f"seed_{seeds[0]}_{modes[0]}"
    for mode in modes:
        for seed in seeds:
            run = output/f"seed_{seed}_{mode}"; runs.append(run)
            cmd = [sys.executable,str(root/"f/panome.py"),*rest,"--seed",str(seed),"--run-dir",str(run),"--train-only"]
            if run!=first and not args.vary_splits:
                cmd += ["--split-file",str(first/"split_before_qc.csv")]
            if mode in ["technical_only","age_sex_residualized"]:
                cmd += ["--residualize",f"{biom}.plate" if mode=="technical_only" else f"age,sex,{biom}.plate"]
            if mode=="shuffled_outcomes": cmd += ["--shuffle-development-outcomes"]
            commands.append(cmd)
    (output/"commands.json").write_text(json.dumps(commands,indent=2)+"\n")
    if args.dry_run:
        print(output/"commands.json"); return
    for cmd in commands:
        subprocess.run(cmd,check=True)
    for run in runs:
        subprocess.run([sys.executable,str(root/"f/panome.py"),"evaluate","--run-dir",str(run)],check=True)
    collect(runs,output)
    print(output/"repeat_summary.csv")


if __name__=="__main__": main()
