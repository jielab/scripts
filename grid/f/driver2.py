#!/usr/bin/env python3
"""GRID FINAL: consensus LD-block PRS, published comparators, CPU twins and optional deep learning."""
from __future__ import annotations
import argparse
import contextlib
import fcntl
import glob
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

GRID=Path(__file__).resolve().parents[1]
SHARED=GRID.parent/"0f"
sys.path.insert(0,str(SHARED))

def boolean(value):
    if str(value).lower() in {"true","1","yes"}:return True
    if str(value).lower() in {"false","0","no"}:return False
    raise argparse.ArgumentTypeError("Use TRUE or FALSE")

ENTRY_MODULES = {
    "grid": ["run","paper-run","regions","block_models","abm","twins","abm_deep","evaluate","architecture","robustness","final","check"],
    "pca": ["all","projection","pca"],
    "csx": ["all","prepare","inputs","infer","score","fit","infer_prepare","infer_chr"],
    "disco": ["all","disco"],
    "prsformer": ["all","prsformer"],
}

def parser(entry="grid"):
    p=argparse.ArgumentParser(description=f"{entry}: independent analysis workflow",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("module",nargs="?",default="run" if entry=="grid" else "all",choices=ENTRY_MODULES[entry])
    p.set_defaults(entry=entry,with_prsformer=False,prsformer_stage="all")
    p.add_argument("--grch", type=int, choices=[38], default=38, help="Target genome build; example: --grch 38")
    p.add_argument("--config",default=str(GRID/"config/grid.yml"))
    p.add_argument("--analysis-root");p.add_argument("--traits",help="comma-separated height,ldl,t2dm")
    p.set_defaults(train_n=None,label_budget=None,phi_grid=None,geometry=None,
                   bootstrap=None,with_deep=False,abm_deep_run=None)
    p.add_argument("--genotype-dir","--dir-imp",dest="genotype_dir")
    p.add_argument("--threads",type=int)
    if entry=="csx":
        p.add_argument("--chr","--chrs",dest="chromosomes",help="autosomes 1-22 or a subset, e.g. 1,2")
        p.add_argument("--worker-chromosome",type=int,help=argparse.SUPPRESS)
        p.add_argument("--infer-jobs",type=int,help="concurrent model chains inside each chromosome worker [4]")
        p.add_argument("--blas-threads",type=int,help="BLAS threads per MCMC chain; default 1")
        p.add_argument("--infer-memory-gb",type=float,default=16,help="host RAM per chromosome worker [16 GiB]")
    if entry=="grid":
        p.add_argument("--consensus-dir",help="existing GRCh38 chr*.bed catalog; never rebuilds LD")
        p.add_argument("--geometry",choices=["consensus","regional","global_prs","legacy_unit_block"])
        p.add_argument("--bootstrap",type=int)
        deep=p.add_mutually_exclusive_group()
        deep.add_argument("--with-deep",action="store_true",help="include deep (enabled in default config)")
        deep.add_argument("--without-deep",action="store_true",help="omit deep from run")
        p.add_argument("--gpu-queue",help="deep GPU queue; overrides config/environment")
        p.add_argument("--gpu-cores",type=int,default=8,help="CPU slots per deep GPU job")
        p.add_argument("--gpu-memory-limit-gb",type=float,default=float(os.environ.get("HPC_GPU_MEMORY_GIB","160")),help="deep GPU host RAM cap; also bounded by --job-memory-gb")
        p.add_argument("--abm-deep-run",choices=["cpu","gpu","auto"],help="default config uses GPU; auto means GPU; cpu opts into CPU")
    else:
        for name in ("cohort-file","pca-file","ancestry-file","group-file","relationship-file","split-file","gwas-manifest"):
            p.add_argument("--"+name)
        p.add_argument("--phenotype-file","--pheno-file",dest="phenotype_file")
        p.add_argument("--train-n",type=int)
        p.add_argument("--label-budget",type=int)
        p.add_argument("--phi-grid")
    if entry=="prsformer":p.add_argument("--prsformer-stage",choices=["all","plan","prepare","train","predict"],default="all")
    p.add_argument("--submit",choices=["auto","local","hpc"],default=os.environ.get("submit","auto"))
    p.add_argument("--queue",default="auto");p.add_argument("--cores","--job-cores",dest="job_cores",type=int,
        help="CPU budget from shared 0f/0os.sh WORKFLOW_CORES; --job-cores is a compatibility alias")
    p.add_argument("--job-memory-gb",type=float,default=128 if os.environ.get("submit")=="hpc" else 32);p.add_argument("--job-hours",type=float,default=120)
    p.add_argument("--jobs",type=int,default=1,help="concurrent workflow jobs")
    p.add_argument("--plan-only",action="store_true");p.add_argument("--wait",type=boolean,default=False)
    p.add_argument("--replace",type=boolean,default=False,help="recompute requested stages; otherwise reuse matching completed outputs")
    p.add_argument("--worker",action="store_true",help=argparse.SUPPRESS)
    return p

def shared_cpu_budget(mode):
    """Use the shell's single source of defaults, including direct Python use."""
    env=dict(os.environ,submit=mode,SUBMIT=mode)
    value=subprocess.check_output(["bash","-c",
        'source "$1" || exit; printf "%s" "$WORKFLOW_CORES"',
        "grid-cpu-budget",str(SHARED/"0os.sh")],env=env,text=True).strip()
    try:cores=int(value)
    except ValueError:raise ValueError("0f/0os.sh returned an invalid WORKFLOW_CORES")
    if not 1<=cores<=10000:raise ValueError("WORKFLOW_CORES must be in 1..10000")
    return cores


def configure(args):
    from core2 import load_config
    cfg=load_config(args.config)
    if args.job_cores is None:
        if args.worker:
            args.job_cores=int(os.environ.get("WORKFLOW_CORES",cfg["runtime"]["threads"]))
        else:
            args.job_cores=shared_cpu_budget(args.submit)
    if not 1<=args.job_cores<=10000:raise ValueError("--cores must be in 1..10000")
    if args.entry=="csx":
        from genomics import parse_chromosomes
        cfg["chromosomes"]=parse_chromosomes(getattr(args,"chromosomes",None) or ",".join(map(str,cfg["chromosomes"])))
        chrom=getattr(args,"worker_chromosome",None)
        if args.module=="infer_chr":
            if not args.worker or chrom not in cfg["chromosomes"]:
                raise ValueError("infer_chr is an internal worker; use 2.csx.sh infer --chr 1-22")
            cfg["chromosomes"]=[chrom]
        if args.module=="infer_prepare" and not args.worker:
            raise ValueError("infer_prepare is an internal worker; use 2.csx.sh infer")
    cfg.setdefault("worker_python",os.environ.get("GRID_WORKER_PYTHON",sys.executable))
    for name in ("cohort_file","pca_file","ancestry_file","group_file","relationship_file","split_file","gwas_manifest","phenotype_file","genotype_dir","consensus_dir"):
        if getattr(args,name,None):cfg["inputs"][name]=None if getattr(args,name).lower()=="none" else str(Path(getattr(args,name)).expanduser().resolve())
    if args.analysis_root:cfg["analysis_root"]=str(Path(args.analysis_root).expanduser().resolve())
    cfg["genome_build"] = args.grch
    for key in ("genotype_dir", "pca_file", "ancestry_file", "consensus_dir"):
        value = str(cfg["inputs"].get(key, ""))
        if "/37/" in value or "/hg19/" in value:
            raise ValueError(f"{key} conflicts with --grch 38: {value}")
    if str(cfg["inputs"].get("consensus_build", "GRCh38")) not in ("38", "GRCh38"):
        raise ValueError("consensus_build must match --grch 38")
    if args.traits:cfg["traits"]=[t.strip() for t in args.traits.split(",") if t.strip()]
    if not cfg["traits"] or set(cfg["traits"])-{"height","ldl","t2dm"}:raise ValueError("Traitsmustbeheight,ldl,t2dm")
    if args.train_n is not None:cfg["cohort"]["train_n"]=args.train_n
    if args.label_budget is not None:cfg["cohort"]["label_budget"]=args.label_budget
    if args.phi_grid:cfg["inference"]["phi_grid"]=[float(x) for x in args.phi_grid.split(",")]
    import math
    grid=[float(x) for x in cfg["inference"]["phi_grid"]]
    if not grid or any(not math.isfinite(x) or x<=0 for x in grid):raise ValueError("inference.phi_grid must contain positive finite values")
    fixed=cfg["benchmark"].get("input_phi",0.01)
    if str(fixed).lower()!="auto" and not any(math.isclose(float(fixed),x,rel_tol=1e-10,abs_tol=0) for x in grid):
        raise ValueError("benchmark.input_phi must occur in inference.phi_grid (or be auto); it fixes the block/Disco input independently of grid order")
    if args.geometry:cfg["abm"]["geometry"]=args.geometry
    if args.bootstrap is not None:cfg["evaluation"]["bootstrap"]=args.bootstrap
    if args.threads is not None:cfg["runtime"]["threads"]=args.threads
    if args.threads is None and not args.worker:cfg["runtime"]["threads"]=args.job_cores
    for option,key in (("infer_jobs","jobs"),("blas_threads","blas_threads")):
        if getattr(args,option,None) is not None:cfg["inference"][key]=getattr(args,option)
    inference_resources(cfg)
    from workflow_mode import resolve_mode
    if resolve_mode(args.submit) == "hpc":
        cfg["runtime"]["rscript"] = str(SHARED / "0os.r.sh")
    os.environ["GRID_RSCRIPT"]=cfg["runtime"].get("rscript","Rscript")
    if args.entry=="prsformer":cfg["prsformer"]["enabled"]=True
    if args.with_deep or args.module=="abm_deep":cfg["block_deep"]["enabled"]=True
    if getattr(args,"without_deep",False):
        if args.module=="abm_deep":raise ValueError("--without-deep cannot be used with abm_deep")
        cfg["block_deep"]["enabled"]=False
    if args.abm_deep_run:cfg["block_deep"]["device"]="cpu" if args.abm_deep_run=="cpu" else "cuda"
    if args.entry=="grid":
        if cfg["block_deep"]["device"]=="auto":cfg["block_deep"]["device"]="cuda"
        if args.gpu_cores<1 or not math.isfinite(args.gpu_memory_limit_gb) or args.gpu_memory_limit_gb<=0:
            raise ValueError("GPU cores and host RAM cap must be positive")
    if cfg["inputs"].get("split_file")=="auto":
        old=cfg["inputs"].get("legacy_split_file")
        candidates=[old,*cfg["inputs"].get("legacy_split_fallbacks",[])]
        cfg["inputs"]["split_file"]=next((p for p in candidates if p and Path(p).is_file()),old)
        # A missing preserved TEST roster must not silently create a new split.
    if min(cfg["cohort"].get("train_n",0),cfg["cohort"].get("label_budget",0))<0 or cfg["runtime"]["threads"]<1:raise ValueError("Invalid training budget/threads")
    return cfg

def inference_resources(cfg):
    """Bound chain concurrency separately from each chain's BLAS thread count."""
    options=cfg["inference"]
    jobs=int(options.get("jobs",4))
    blas=int(options.get("blas_threads",1))
    budget=int(cfg["runtime"]["threads"])
    if min(jobs,blas,budget)<1 or blas>budget:
        raise ValueError("Inference jobs/BLAS threads must be positive and fit the CPU budget")
    return min(jobs,budget//blas),blas

def manifest_path(cfg,root):
    value=cfg["inputs"].get("gwas_manifest")
    return Path(root)/"tmp/gwas_manifest.tsv" if value in (None,"auto") else Path(value)

def dependencies(cfg,root,stage,trait=None):
    """Readable stat records. No whole-genome hashing and no test truth in fit inputs."""
    from core2 import stamp,VERSION
    root=Path(root);work=root/f"tmp/{trait}";inp=cfg["inputs"]
    cohort=root/"tmp/cohort/participants.tsv.gz";dev=root/f"tmp/cohort/{trait}.development.tsv.gz"
    files=[];settings={"version":VERSION,"stage":stage,"trait":trait}
    def add(*values):files.extend(Path(v) for v in values if v)
    if stage=="pca":
        settings.update(pca=cfg["pca"],distance_pcs=cfg["disco"]["distance_pcs"])
        add(inp["pca_weight"],inp["pca_centers"],inp["phenotype_file"],inp.get("rsid_map"))
    elif stage=="prepare":
        settings.update(cohort=cfg["cohort"],covariates=cfg["benchmark"]["covariates"],traits=cfg["traits"])
        for k in ("cohort_file","phenotype_file","pca_file","ancestry_file","group_file","relationship_file","split_file"):add(inp.get(k))
    elif stage=="inputs":
        settings.update(n_gwas=cfg["inference"].get("n_gwas"),effective_n_columns=cfg["inference"].get("effective_n_columns",{}),overlap=cfg["inference"].get("discovery_target_overlap"),traits=cfg["traits"])
        if inp.get("gwas_manifest") not in (None,"auto"):add(inp["gwas_manifest"])
        else:
            from inputs2 import find_gwas
            from core2 import POPS
            for t in cfg["traits"]:
                for pop in POPS:add(find_gwas(inp["gwas_dir"],t,pop))
            add(Path(inp["ld_ref"])/"snpinfo_mult_1kg_hm3")
    elif stage in {"infer","infer_prepare","infer_chr"}:
        settings.update(inference=cfg["inference"],chromosomes=cfg["chromosomes"])
        add(manifest_path(cfg,root),inp.get("rsid_map"))
        from infer import read_manifest
        for row in read_manifest(manifest_path(cfg,root),trait).values():add(*row["files"])
        for p in Path(inp["ld_ref"]).glob("**/*"):
            if p.is_file() and (p.name.startswith("snpinfo") or p.suffix==".hdf5"):add(p)
        add(GRID/"f/infer.py",GRID/"f/csx_parallel.py",GRID/"f/genomics.py")
        if stage in {"infer","infer_chr"}:add(work/"inference/prepared.json")
        if stage=="infer":
            for c in cfg["chromosomes"]:add(work/f"inference/partitions/chr{c}/partition.complete.json")
    elif stage=="score":
        settings.update(scoring=cfg["scoring"]);add(cohort,work/"inference/weights/manifest.tsv")
        add(*sorted((work/"inference/weights").glob("*.chr*.tsv")))
    elif stage=="fit":
        settings.update(benchmark=cfg["benchmark"],seed=cfg["cohort"]["seed"]);add(cohort,dev,work/"scores.tsv.gz",work/"inference/weights/manifest.tsv")
    elif stage=="disco":
        settings.update(disco=cfg["disco"],benchmark=cfg["benchmark"]);add(cohort,dev,work/"scores.tsv.gz",work/"baseline.predictions.tsv.gz",work/"benchmark.json")
    elif stage=="regions":
        settings.update(scoring=cfg["scoring"],regions=cfg["regions"],consensus_build=inp["consensus_build"])
        add(cohort,work/"benchmark.json",work/"inference/weights/manifest.tsv")
        source=Path(inp["consensus_dir"])
        add(*sorted(source.glob("chr*.bed")) if source.is_dir() else [source])
        add(*sorted((work/"inference/weights").glob("*.chr*.tsv")))
    elif stage in {"block_models","abm_deep"}:
        section="block_deep" if stage=="abm_deep" else "block_models"
        settings.update(method={k:v for k,v in cfg[section].items() if k!="enabled"},seed=cfg["cohort"]["seed"])
        add(cohort,dev,work/"baseline.predictions.tsv.gz",work/"regions/scores.npy",work/"regions/ids.tsv",work/"regions/features.tsv")
    elif stage in {"abm","robustness"}:
        settings.update(abm=cfg["abm"],seed=cfg["cohort"]["seed"]);add(cohort,dev,work/"predictions.tsv.gz")
        if stage=="robustness" or cfg["abm"]["geometry"]!="global_prs":add(work/"regions/scores.npy",work/"regions/ids.tsv",work/"regions/features.tsv")
        if stage=="robustness":add(root/f"tmp/cohort/_private/{trait}.test.tsv.gz")
    elif stage=="evaluate":
        settings.update(evaluation=cfg["evaluation"],seed=cfg["cohort"]["seed"])
        add(cohort,root/f"tmp/cohort/_private/{trait}.test.tsv.gz",work/"predictions.tsv.gz",work/"abm/abm.predictions.tsv.gz",work/"abm/model.json")
        for path in verified_optional_predictions(cfg,root,trait).values():add(path)
    elif stage=="twins":
        add(cohort,work/"abm/abm.predictions.tsv.gz",work/"abm/model.json")
    elif stage=="architecture":add(cohort,work/"regions/scores.npy",work/"regions/features.tsv",work/"regions/ids.tsv")
    elif stage=="final":
        settings.update(traits=cfg["traits"])
        for t in cfg["traits"]:
            add(root/f"final/{t}.performance.tsv",root/f"final/{t}.ancestry.tsv",root/f"tmp/{t}/benchmark.json",root/f"tmp/{t}/abm/model.json")
            for path in (root/f"final/{t}.twins.html",root/f"tmp/{t}/disco/model.json",root/f"tmp/{t}/block_models/model.json",root/f"tmp/{t}/deep/model.json"):
                if path.is_file():add(path)
    elif stage=="prsformer":
        settings.update(prsformer=cfg["prsformer"],traits=cfg["traits"]);add(cohort)
        for t in cfg["traits"]:add(root/f"tmp/cohort/{t}.development.tsv.gz")
        add(root/f"tmp/{cfg['traits'][0]}/inference/target/manifest.tsv")
    if stage in {"pca","prepare","infer","infer_prepare","infer_chr","score","regions","prsformer"}:
        for c in cfg["chromosomes"]:
            prefix=Path(inp["genotype_dir"])/f"chr{c}"
            add(str(prefix)+".psam")
            pv=Path(str(prefix)+".pvar");add(pv if pv.is_file() else str(prefix)+".pvar.zst")
            if stage in {"pca","score","regions","prsformer"}:add(str(prefix)+".pgen")
    # When a ready cohort supplies families, unused relationship/phenotype paths
    # are not mandatory. Missing relevant files are raised by the actual worker.
    return {"settings":settings,"files":[stamp(p) if p.is_file() else {"path":str(p.resolve()),"missing":True} for p in dict.fromkeys(files)]}

def stage_outputs(cfg,root,stage,trait=None):
    root=Path(root);work=root/f"tmp/{trait}"
    if stage=="pca":return [Path(cfg["inputs"]["pca_file"]),Path(cfg["inputs"]["ancestry_file"])]
    if stage=="infer_prepare":
        from csx_parallel import prepared_outputs
        return prepared_outputs(work/"inference",cfg["chromosomes"])
    if stage=="infer_chr":
        from infer import model_manifest,parse_phi_grid
        c=cfg["chromosomes"][0];part=work/f"inference/partitions/chr{c}"
        models=model_manifest(part,parse_phi_grid(",".join(map(str,cfg["inference"]["phi_grid"]))),cfg["inference"].get("prscs_fixed",False))
        return [part/"partition.complete.json",part/"provenance.json",part/"weights/manifest.tsv",part/"weights/snp_counts.tsv",
                *[part/f"weights/{m}.chr{c}.tsv" for m in models.score_name]]
    if stage=="infer":
        from infer import model_manifest,parse_phi_grid
        models=model_manifest(work/"inference",parse_phi_grid(",".join(map(str,cfg["inference"]["phi_grid"]))),cfg["inference"].get("prscs_fixed",False))
        return [work/"inference/weights/manifest.tsv",work/"inference/provenance.json",work/"inference/target/manifest.tsv",
            *[work/f"inference/target/chr{c}.tsv" for c in cfg["chromosomes"]],
            *[work/f"inference/weights/{m}.chr{c}.tsv" for m in models.score_name for c in cfg["chromosomes"]]]
    if stage=="evaluate":
        return [root/f"final/{trait}.performance.tsv",root/f"final/{trait}.ancestry.tsv",root/f"final/{trait}.report.html",root/f"final/{trait}.evaluation.json",
                *([root/f"final/{trait}.paired_discrimination.tsv"] if trait=="t2dm" else [])]
    return {
      "prepare":[root/"tmp/cohort/participants.tsv.gz",root/"tmp/cohort/cohort.json",*[root/f"tmp/cohort/{t}.development.tsv.gz" for t in cfg["traits"]],*[root/f"tmp/cohort/_private/{t}.test.tsv.gz" for t in cfg["traits"]]],
      "inputs":[manifest_path(cfg,root)],
      "score":[work/"scores.tsv.gz"],"fit":[work/"baseline.predictions.tsv.gz",work/"benchmark.json",work/"benchmark_models.joblib"],
      "disco":[work/"predictions.tsv.gz",work/"disco/model.json"],
      "regions":[work/"regions"/name for name in ("scores.npy","ids.tsv","features.tsv","catalog.tsv","sum_check.tsv","region_counts.tsv","covered_global.tsv.gz","unassigned_variants.tsv.gz","metadata.json")],
      "block_models":[work/"block_models.predictions.tsv.gz",work/"block_models/model.json",work/"block_models/model.joblib"],
      "abm_deep":[work/"deep"/name for name in ("predictions.tsv.gz","model.json","model.pt","bundle.joblib","training.tsv","tuning.tsv")],
      "abm":[work/"abm/abm.predictions.tsv.gz",work/"abm/model.json",work/"abm/abm.model.joblib",
             *([work/"abm/geometry.npy",work/"abm/twins.tsv.gz"] if cfg["abm"]["geometry"]=="consensus" else [])],
      "twins":[root/f"qc/{trait}.twins.tsv",root/f"final/{trait}.twins.html"],
      "architecture":[root/f"qc/{trait}.regional_architecture.tsv"],
      "robustness":[root/f"final/{trait}.geometry_sensitivity.tsv"],
      "final":[root/"final/index.html",root/"final/all.performance.tsv",root/"final/all.ancestry.tsv",
               *[root/f"final/{stem}.{suffix}" for stem in ("overall_performance","support_risk_coverage") for suffix in ("png","pdf")]],
      "prsformer":[root/f"tmp/{t}/prsformer.predictions.tsv.gz" for t in cfg["traits"]]
    }[stage]

def read_manifest_weights(work):
    import pandas as pd
    return pd.read_csv(Path(work)/"inference/weights/manifest.tsv",sep="\t",dtype=str)

def regional_inputs(root,trait,cohort):
    import numpy as np
    from core2 import read_table
    from abm2 import IndexedFeatures
    folder=Path(root)/f"tmp/{trait}/regions"
    data=np.load(folder/"scores.npy",mmap_mode="r");ids=read_table(folder/"ids.tsv").id.astype(str)
    lookup={v:i for i,v in enumerate(ids)}
    if len(lookup)!=len(ids) or any(x not in lookup for x in cohort.id):raise ValueError("RegionalfeaturesIDalignmentfailed")
    index=np.array([lookup[x] for x in cohort.id],dtype=int)
    features=data if np.array_equal(index,np.arange(len(data))) else IndexedFeatures(data,index)
    manifest=read_table(folder/"features.tsv")
    key=next(c for c in ("map","population","pop") if c in manifest)
    return features,manifest[key].astype(str).to_numpy()

def run_abm(cfg,root,trait,geometry=None,output=None):
    from core2 import read_table,identified
    from abm2 import fit_abm
    root=Path(root);work=root/f"tmp/{trait}";o=cfg["abm"]
    cohort=identified(read_table(root/"tmp/cohort/participants.tsv.gz"));mode=geometry or o["geometry"]
    features,maps=(None,None) if mode=="global_prs" else regional_inputs(root,trait,cohort)
    return fit_abm(cohort,read_table(root/f"tmp/cohort/{trait}.development.tsv.gz"),read_table(work/"predictions.tsv.gz"),output or work/"abm",
        trait=trait,features=features,feature_maps=maps,feature_manifest=None if mode=="global_prs" else read_table(work/"regions/features.tsv"),
        geometry=mode,dimensions=o["dimensions"],k_values=o["k_values"],
        radius_k=o["radius_k"],radius_quantile=o["radius_quantile"],seed=cfg["cohort"]["seed"],workers=cfg["runtime"]["threads"],
        batch_size=o["batch_size"],min_count=o["min_count"],min_group_ess=o["min_group_ess"],shrinkage=o["shrinkage"],
        donor_batch=o.get("donor_batch",4096),evidence_queries=o.get("evidence_queries",200),evidence_neighbors=o.get("evidence_neighbors",16))

def fit_regional(cfg,root,trait,deep=False):
    from core2 import read_table,identified
    root=Path(root);work=root/f"tmp/{trait}"
    cohort=identified(read_table(root/"tmp/cohort/participants.tsv.gz"))
    features,_=regional_inputs(root,trait,cohort)
    common=(cohort,read_table(root/f"tmp/cohort/{trait}.development.tsv.gz"),
            read_table(work/"baseline.predictions.tsv.gz"),features,read_table(work/"regions/features.tsv"))
    options=dict(trait=trait,seed=cfg["cohort"]["seed"],threads=cfg["runtime"]["threads"])
    if deep:
        from block_deep import fit_block_deep
        return fit_block_deep(*common,work/"deep",config=cfg["block_deep"],**options)
    from block_models import fit_block_models
    return fit_block_models(*common,work/"block_models",config=cfg["block_models"],rscript=cfg["runtime"]["rscript"],**options)

def verified_optional_predictions(cfg,root,trait,*,load=False):
    """Include separately completed optional stages only on their original inputs.

    Enabling controls scheduling, not evaluation membership. Light file records
    prevent a same-ID prediction from an older training split being mixed in.
    Device/settings are taken from the fitted model; GPU fit -> CPU evaluation
    is allowed. No genotype bytes or TEST outcomes are read by this guard.
    """
    from core2 import VERSION,stamp,read_table,identified
    root=Path(root);work=root/f"tmp/{trait}";result={}
    specs=(("block_models",work/"block_models.predictions.tsv.gz",f"{trait}.block_models"),
           ("abm_deep",work/"deep/predictions.tsv.gz",f"{trait}.abm_deep"),
           ("prsformer",work/"prsformer.predictions.tsv.gz","prsformer"))
    roster=set(identified(read_table(root/"tmp/cohort/participants.tsv.gz")).id) if load else None
    for stage,path,name in specs:
        if not path.is_file():continue
        marker=root/"tmp/state"/f"{name}.json"
        command="4.prsformer.sh" if stage=="prsformer" else f"grid.sh {stage}"
        problem=f"{stage}: stale or unverified optional predictions; rerun {command} --replace TRUE on this analysis root"
        if not marker.is_file():raise ValueError(problem+" (missing producer receipt)")
        previous=json.loads(marker.read_text());record=previous.get("inputs",{})
        producer=dict(cfg)
        if stage=="prsformer":producer["traits"]=record.get("settings",{}).get("traits",cfg["traits"])
        current=dependencies(producer,root,stage,None if stage=="prsformer" else trait)
        if (record.get("settings",{}).get("version")!=VERSION or record.get("files")!=current["files"]
                or previous.get("outputs",{}).get(str(path))!=stamp(path)):
            raise ValueError(problem+" (producer inputs or output changed)")
        if load:
            frame=identified(read_table(path),f"{stage} predictions")
            if set(frame.id)!=roster:raise ValueError(problem+" (prediction roster differs)")
            result[stage]=frame
        else:result[stage]=path
    return result

def evaluate(cfg,root,trait,abm_path=None,output=None):
    from core2 import read_table,identified
    from evaluate2 import evaluate_predictions
    root=Path(root);work=root/f"tmp/{trait}"
    benchmark=read_table(work/"predictions.tsv.gz")
    abm=read_table(abm_path or work/"abm/abm.predictions.tsv.gz")
    optional=verified_optional_predictions(cfg,root,trait,load=True)
    if "block_models" in optional:
        benchmark=benchmark.merge(optional["block_models"],on="id",how="left",validate="one_to_one")
    if "abm_deep" in optional:
        deep=optional["abm_deep"]
        benchmark=benchmark.merge(deep[["id","abm.deep.prediction"]],on="id",how="left",validate="one_to_one")
        columns=[c for c in deep if c.startswith("abm.deep.selected_")]
        masks=deep[["id",*columns]].rename(columns={c:c.replace("abm.deep.selected_","selected_deep_") for c in columns})
        abm=abm.merge(masks,on="id",how="left",validate="one_to_one")
    if "prsformer" in optional:
        benchmark=benchmark.merge(optional["prsformer"],on="id",how="left",validate="one_to_one")
    return evaluate_predictions(identified(read_table(root/"tmp/cohort/participants.tsv.gz")),read_table(root/f"tmp/cohort/_private/{trait}.test.tsv.gz"),benchmark,
        abm,output or root/"final",trait=trait,bootstrap=cfg["evaluation"]["bootstrap"],seed=cfg["cohort"]["seed"],workers=cfg["runtime"]["threads"])

def execute(cfg,root,stage,trait,args):
    from core2 import atomic_json,atomic_table,read_table
    root=Path(root);work=root/f"tmp/{trait}";inp=cfg["inputs"]
    if stage=="pca":
        from pca2 import project;project(cfg,root)
    elif stage=="prepare":
        from cohort2 import prepare;prepare(cfg,root,cfg["traits"])
    elif stage=="inputs":
        from inputs2 import ensure_manifest;ensure_manifest(cfg,root,cfg["traits"])
    elif stage in {"infer","infer_prepare","infer_chr"}:
        o=cfg["inference"]
        jobs,blas=inference_resources(cfg)
        print(f"[infer resources] concurrent chains={jobs}; BLAS threads/chain={blas}; CPU budget={cfg['runtime']['threads']}",flush=True)
        command=[sys.executable,str(GRID/"f/infer.py"),"--manifest",str(manifest_path(cfg,root)),"--trait",trait,"--ld-ref",inp["ld_ref"],"--genotype-dir",inp["genotype_dir"],"--output",str(work/"inference"),"--chrs",",".join(map(str,cfg["chromosomes"])),"--phi-grid",",".join(map(str,o["phi_grid"])),"--jobs",str(jobs),"--threads",str(blas),"--backend",o["backend"],"--n-iter",str(o["n_iter"]),"--n-burnin",str(o["n_burnin"]),"--thin",str(o["thin"]),"--seed",str(o["seed"])]
        if o.get("prscs_fixed"):command.append("--prscs-fixed")
        if inp.get("rsid_map"):command += ["--rsid-map",inp["rsid_map"]]
        if args.replace:command.append("--replace")
        if stage=="infer_prepare":command.append("--prepare-only")
        elif stage=="infer_chr":
            command[command.index("--output")+1]=str(work/f"inference/partitions/chr{cfg['chromosomes'][0]}")
            command += ["--prepared-inputs",str(work/"inference")]
        else:command.append("--merge-only")
        subprocess.run(command,check=True)
    elif stage=="score":
        from scoring_parallel import score_models
        m=read_manifest_weights(work)
        score_models(root/"tmp/cohort/participants.tsv.gz",dict(zip(m.score_name,m.weight_pattern)),inp["genotype_dir"],work/"scores.tsv.gz",
            common_variants=cfg["scoring"]["common_variants"],plink=cfg["runtime"]["plink"],threads=cfg["runtime"]["threads"],memory_mb=cfg["runtime"]["plink_memory_mb"],batch_size=cfg["scoring"]["batch_size"])
    elif stage=="fit":
        from benchmark2 import fit_models;fit_models(cfg,root,trait)
    elif stage=="disco":
        from benchmark2 import fit_disco;fit_disco(cfg,root,trait)
    elif stage=="regions":
        from scoring_parallel import regional_scores
        m=read_manifest_weights(work);meta=json.loads((work/"benchmark.json").read_text())
        mode=cfg["regions"]["representation"]
        if mode not in {"meta","ancestry_channels"}:raise ValueError("regions.representation must be meta or ancestry_channels")
        family="csx_meta" if mode=="meta" else "csx"
        fixed=m[m.family.eq(family)&m.phi.eq(meta["geometry_phi"])]
        if len(fixed)!=(1 if mode=="meta" else 4):raise ValueError("Incomplete fixed-phi regional weights")
        regional_scores(root/"tmp/cohort/participants.tsv.gz",dict(zip(fixed.population,fixed.weight_pattern)),inp["consensus_dir"],inp["genotype_dir"],work/"regions",plink=cfg["runtime"]["plink"],threads=cfg["runtime"]["threads"],memory_mb=cfg["runtime"]["plink_memory_mb"],batch_size=cfg["scoring"]["batch_size"],block_build=inp["consensus_build"])
    elif stage=="block_models":fit_regional(cfg,root,trait)
    elif stage=="abm_deep":fit_regional(cfg,root,trait,deep=True)
    elif stage=="abm":run_abm(cfg,root,trait)
    elif stage=="twins":
        from report2 import twin_report;twin_report(cfg,root,trait)
    elif stage=="evaluate":evaluate(cfg,root,trait)
    elif stage=="architecture":
        from report2 import architecture;architecture(cfg,root,trait)
    elif stage=="robustness":
        import pandas as pd
        rows=[]
        for mode in ("consensus","global_prs","regional","legacy_unit_block"):
            location=work/"geometry_sensitivity"/mode
            run_abm(cfg,root,trait,mode,location)
            performance,_,_=evaluate(cfg,root,trait,location/"abm.predictions.tsv.gz",location/"evaluation")
            rows.append(performance.assign(geometry=mode,interpretation="predeclared exploratory sensitivity; no test selection of primary geometry"))
        atomic_table(pd.concat(rows,ignore_index=True),root/f"final/{trait}.geometry_sensitivity.tsv")
    elif stage=="final":
        from report2 import final_report;final_report(cfg,root)
    elif stage=="prsformer":
        from copy import deepcopy
        from prsformer2 import run
        worker_cfg=deepcopy(cfg)
        if args.replace and args.prsformer_stage in {"all","train"}:
            worker_cfg["prsformer"]["force_retrain"]=True
        if args.replace and args.prsformer_stage in {"all","prepare"}:
            worker_cfg["prsformer"].setdefault("prepare",{})["replace"]=True
        run(worker_cfg,root,cfg["traits"],stage=args.prsformer_stage)
    else:raise ValueError(stage)

def run_stage(cfg,root,stage,trait,args):
    from core2 import atomic_json,stamp
    root=Path(root);name=f"{trait}.{stage}" if trait else stage
    if stage=="infer_chr":name=f"{trait}.chr{cfg['chromosomes'][0]}.infer"
    state=root/"tmp/state";state.mkdir(parents=True,exist_ok=True)
    with (state/f"{name}.lock").open("a+") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        signature=dependencies(cfg,root,stage,trait);marker=state/f"{name}.json"
        outputs=stage_outputs(cfg,root,stage,trait)
        # Inference workers validate their detailed preparation/chain receipts
        # themselves, including raw posteriors and code identities.
        if marker.is_file() and not args.replace and stage not in {"infer","infer_prepare","infer_chr"}:
            previous=json.loads(marker.read_text())
            same_outputs=all(p.is_file() and previous.get("outputs",{}).get(str(p))==stamp(p) for p in outputs)
            if previous.get("inputs")==signature and same_outputs:
                print(f"[reuse] {name}",flush=True);return
            if previous.get("inputs")!=signature:
                raise ValueError(f"{name}: inputs/settings changed; use --replace TRUE for this stage, or a new --analysis-root")
            print(f"[repair] {name}: missing or changed output; reconstruct from matching inputs",flush=True)
        print(f"[run] {name}",flush=True)
        execute(cfg,root,stage,trait,args)
        if stage=="prsformer" and args.prsformer_stage in {"plan","prepare"}:return
        if not all(p.is_file() for p in outputs):raise RuntimeError(f"{name}:expectedoutputmissing")
        atomic_json({"inputs":signature,"outputs":{str(p):stamp(p) for p in outputs}},marker)

def plan_jobs(cfg,args):
    from core2 import atomic_json
    from workflow_jobs import Job,write_plan,run_plan
    root=Path(cfg["analysis_root"]).resolve();root.mkdir(parents=True,exist_ok=True)
    # Every plan has a separate immutable config snapshot, shared by its workers.
    from tempfile import mkstemp
    folder=root/"tmp/config";folder.mkdir(parents=True,exist_ok=True)
    fd,name=mkstemp(prefix="run-",suffix=".json",dir=folder);os.close(fd)
    config=Path(name);atomic_json(cfg,config)
    jobs=[];last=[]
    worker_python=cfg.get("worker_python",sys.executable)
    def add(stage,trait=None,depends=(),chromosome=None):
        key=f"{trait}.{stage}" if trait else stage
        if chromosome is not None:key=f"{trait}.chr{chromosome}.infer"
        argv=[worker_python,str(GRID/f"f/{args.entry}_driver.py"),stage,"--grch","38","--config",str(config),"--worker","--replace",str(args.replace)]
        if trait:argv += ["--traits",trait]
        if chromosome is not None:argv += ["--worker-chromosome",str(chromosome)]
        if stage=="prsformer":argv += ["--prsformer-stage",args.prsformer_stage]
        resources={"queue":args.queue,"cores":args.job_cores,"memory_gib":args.job_memory_gb,"hours":args.job_hours}
        blas=cfg["runtime"]["threads"]
        if args.entry=="csx" and stage in {"infer_prepare","infer_chr","infer","inputs"}:
            _,blas=inference_resources(cfg)
            # Reserve the requested job CPUs independently of BLAS parallelism.
            # The HPC default is 40 slots, including chromosome/preparation jobs.
            resources.update(memory_gib=args.infer_memory_gb)
            if stage!="infer_chr":blas=1
        if stage=="prsformer" and args.prsformer_stage in {"all","train","predict"}:
            resources.update(queue=os.environ.get("HPC_GPU_QUEUE","4a100-40"),gpus=1)
        if stage=="abm_deep" and cfg["block_deep"]["device"] in {"cuda","gpu"}:
            resources.update(queue=args.gpu_queue or cfg["runtime"].get("gpu_queue",os.environ.get("HPC_GPU_QUEUE","4a100-40")),
                             cores=args.gpu_cores,memory_gib=min(args.job_memory_gb,args.gpu_memory_limit_gb),gpus=1)
            blas=min(cfg["runtime"]["threads"],args.gpu_cores)
            argv += ["--threads",str(blas)]
        jobs.append(Job(key=key,argv=argv,depends=list(depends),resources=resources,cwd=str(GRID),
                        job_name="grid.pca" if args.entry=="pca" and stage=="pca" else (f"grid.disco.{trait}" if args.entry=="disco" and stage=="disco" else None),env={"GRID_MEMORY_CAP_GB":str(resources["memory_gib"]),"PYTHONDONTWRITEBYTECODE":"1","LD_LIBRARY_PATH":os.environ.get("LD_LIBRARY_PATH",""),**{var:str(blas) for var in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS")}}))
        return key
    if args.entry=="grid":
        from grid_workflow import build_jobs
        build_jobs(cfg,args,add)
    else:
        from standalone_workflows import build_jobs
        build_jobs(cfg,args,add)
    os.environ["WORKFLOW_PYTHON"]=worker_python
    plan=write_plan(jobs,root/"log",mode=args.submit,name=args.entry,common_env={"WORKFLOW_PYTHON":worker_python})
    print(f"Plan:{plan}",flush=True)
    for j in jobs:print(f"{j.key}: after {','.join(j.depends) or '-'}",flush=True)
    print(run_plan(plan,submit=args.submit,jobs=args.jobs,wait=args.wait,dry_run=args.plan_only),flush=True)

def main(argv=None,entry="grid"):
    args=parser(entry).parse_args(argv);cfg=configure(args)
    if args.module in {"all","projection"} and entry in {"pca","disco","prsformer"}:args.module=entry
    if entry=="grid" and not args.worker and not args.plan_only:
        from grid_workflow import validate_inputs
        validate_inputs(cfg,args.module)
    if entry!="grid" and not args.worker and not args.plan_only:
        from standalone_workflows import validate_inputs
        validate_inputs(cfg,args)
    if args.module=="check":
        from report2 import check_inputs;return check_inputs(cfg)
    if args.worker:
        stage=args.module
        trait=None if stage in {"pca","prepare","inputs","final","prsformer"} else cfg["traits"][0]
        # Shared 0f intentionally starts BLAS at one thread. Apply this stage's
        # explicit limit after NumPy import so full-space ABM can use its CPUs.
        from threadpoolctl import threadpool_limits
        blas=1 if entry=="csx" and stage in {"infer","infer_prepare","infer_chr","inputs"} else cfg["runtime"]["threads"]
        allocated=int(os.environ.get("WORKFLOW_CORES",os.environ.get("LSB_DJOB_NUMPROC",cfg["runtime"]["threads"])))
        cfg["runtime"]["threads"]=min(cfg["runtime"]["threads"],allocated)
        blas=min(blas,allocated)
        # The shared scheduler starts BLAS at one thread. Set the environment
        # too: SciPy/R/other libraries imported later are outside the existing
        # threadpoolctl snapshot, and subprocesses inherit only the environment.
        for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
            os.environ[name]=str(blas)
        print(f"[resources] stage={stage}; allocated CPUs={allocated}; CPU budget={cfg['runtime']['threads']}; native threads={blas}",flush=True)
        with threadpool_limits(limits=blas):
            run_stage(cfg,cfg["analysis_root"],stage,trait,args)
    else:plan_jobs(cfg,args)

def entry_main(entry="grid"):
    try:main(entry=entry)
    except (ValueError,FileNotFoundError,RuntimeError,BlockingIOError) as exc:
        print(f"GRID ERROR: {exc}",file=sys.stderr);raise SystemExit(2)


if __name__=="__main__":
    entry_main()
