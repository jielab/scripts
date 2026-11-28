#!/usr/bin/env python3
"""GRID FINAL: consensus LD-block PRS, published comparators, CPU twins and optional deep learning."""
from __future__ import annotations
from output_layout import project_path, trait_path, module_root, upstream_root
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
    "grid": ["prepare","run","direct_blocknet","blocknet_ablation","blocknet_twins","regions","block_models","abm","twins","abm_deep","architecture","robustness","check"],
    "pca": ["all","projection","pca"],
    "csx": ["all","prepare","inputs","infer","score","fit","infer_prepare","infer_chr"],
    "disco": ["all","disco"],
    "prsformer": ["all","prsformer"],
}

def method_list(value):
    from methods import selection
    try:
        return selection(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def parser(entry="grid", *, paper=False):
    p=argparse.ArgumentParser(description=f"{entry}: independent analysis workflow",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("module",nargs="?",default="run" if entry=="grid" else "all",choices=["run"] if paper else ENTRY_MODULES[entry])
    p.set_defaults(entry=entry,paper=paper,method=None,with_prsformer=False,prsformer_stage="all")
    p.add_argument("--grch", type=int, choices=[38], default=38, help="Target genome build; example: --grch 38")
    p.add_argument("--config",default=str(GRID/"config/grid.yml"))
    p.add_argument("--analysis-root");p.add_argument("--traits",help="comma-separated height,ldl,t2dm")
    p.set_defaults(train_n=None,label_budget=None,phi_grid=None,geometry=None,
                   bootstrap=None,with_deep=False,abm_deep_run=None)
    p.add_argument("--genotype-dir","--dir-imp",dest="genotype_dir")
    p.add_argument("--threads",type=int)
    if entry=="pca":
        p.add_argument("--dir-out-pca_projection", "--dir-out-pca-projection",dest="pca_projection")
    if entry=="csx":
        p.add_argument("--chr","--chrs",dest="chromosomes",help="autosomes 1-22 or a subset, e.g. 1,2")
        p.add_argument("--worker-chromosome",type=int,help=argparse.SUPPRESS)
        p.add_argument("--infer-jobs",type=int,help="concurrent model chains inside each chromosome worker [4]")
        p.add_argument("--blas-threads",type=int,help="BLAS threads per MCMC chain; default 1")
        p.add_argument("--infer-memory-gb",type=float,default=16,help="host RAM per chromosome worker [16 GiB]")
    if entry=="prsformer":
        p.add_argument("--gpu-queue", help="PRSformer GPU queue; overrides config/environment")
        p.add_argument("--gpu-cores", type=int, default=8, help="CPU slots per PRSformer GPU job")
    if entry=="grid":
        if not paper:
            p.add_argument("--method",type=method_list,metavar="METHOD[,METHOD...]",help="Shared BED default: universal.auto and universal.fixed. Use all for all comparisons.")
        p.add_argument("--dir-ld-bed", help="directory containing the selected regional BED file(s)")
        p.add_argument("--ld-bed", help="one shared BED or four comma-separated ancestry BEDs")
        p.add_argument("--ld-layout", choices=["ancestry","uniform","consensus"])
        p.add_argument("--representation", choices=["ancestry_channels","meta"])
        p.add_argument("--ld-block-dir", help="directory with matching EUR/AFR/EAS/SAS.<build>.bed files")
        p.add_argument("--ld-bed-build",choices=["37","38"],help="BED coordinate system (target PGEN is GRCh38)")
        p.add_argument("--run-name", help="separate GRID experiment directory")
        p.add_argument("--consensus-dir",help="existing GRCh38 chr*.bed catalog; never rebuilds LD")
        p.add_argument("--geometry",choices=["consensus","regional","global_prs","unit_block"])
        p.add_argument("--bootstrap",type=int)
        p.add_argument("--ldl-threshold",type=float,help="Observed high-LDL cutoff in mmol/L")
        if paper:
            p.set_defaults(without_deep=False,gpu_cores=8,gpu_memory_limit_gb=160)
        else:
            deep=p.add_mutually_exclusive_group()
            deep.add_argument("--with-deep",action="store_true",help="include deep (enabled in default config)")
            deep.add_argument("--without-deep",action="store_true",help="omit deep from run")
            p.add_argument("--gpu-queue",help="deep GPU queue; overrides config/environment")
            p.add_argument("--gpu-cores",type=int,default=8,help="CPU slots per deep GPU job")
            p.add_argument("--gpu-memory-limit-gb",type=float,default=float(os.environ.get("HPC_GPU_MEMORY_GIB","160")),help="deep GPU host RAM cap; also bounded by --job-memory-gb")
            p.add_argument("--blocknet-run","--abm-deep-run",dest="abm_deep_run",choices=["cpu","gpu","auto"],help="CPU by default for both BlockNet fits; gpu requests CUDA. The prior auto value explicitly requests GPU.")
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
    p.add_argument("--job-memory-gb",type=float,default=(128 if entry=="prsformer" else 80) if os.environ.get("submit")=="hpc" else 32);p.add_argument("--job-hours",type=float,default=120)
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
            from workflow_mode import resolve_mode
            if resolve_mode(args.submit)=="hpc":args.job_cores=min(20,args.job_cores)
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
    for name in ("cohort_file","pca_file","ancestry_file","group_file","relationship_file","split_file","gwas_manifest","phenotype_file","genotype_dir","consensus_dir","ld_block_dir"):
        if getattr(args,name,None):cfg["inputs"][name]=None if getattr(args,name).lower()=="none" else str(Path(getattr(args,name)).expanduser().resolve())
    if getattr(args,"pca_projection",None):
        directory=Path(args.pca_projection).expanduser().resolve()
        for key in ("pca_file","ancestry_file"):
            cfg["inputs"][key]=str(directory/Path(cfg["inputs"][key]).name)
    if args.analysis_root:cfg["analysis_root"]=str(Path(args.analysis_root).expanduser().resolve())
    cfg["genome_build"] = args.grch
    for key in ("genotype_dir", "pca_file", "ancestry_file", "consensus_dir"):
        value = str(cfg["inputs"].get(key, ""))
        if "/37/" in value or "/hg19/" in value:
            raise ValueError(f"{key} conflicts with --grch 38: {value}")
    if str(cfg["inputs"].get("consensus_build", "GRCh38")) not in ("38", "GRCh38"):
        raise ValueError("consensus_build must match --grch 38")
    if args.traits:cfg["traits"]=[t.strip() for t in args.traits.split(",") if t.strip()]
    if not cfg["traits"] or len(set(cfg["traits"])) != len(cfg["traits"]) or set(cfg["traits"])-{"height","ldl","t2dm"}:
        raise ValueError("--traits must be a unique, nonempty subset of height,ldl,t2dm")
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
    if args.entry == "grid":
        configure_region_run(cfg, args)
    if getattr(args,"ldl_threshold",None) is not None:
        if not math.isfinite(args.ldl_threshold) or args.ldl_threshold<=0:raise ValueError("LDL cutoff must be finite and positive")
        cfg["evaluation"]["ldl_threshold"]=args.ldl_threshold
    if args.bootstrap is not None:cfg["evaluation"]["bootstrap"]=args.bootstrap
    if cfg["evaluation"]["bootstrap"] < 0: raise ValueError("--bootstrap must be nonnegative")
    if args.threads is not None:cfg["runtime"]["threads"]=args.threads
    if args.threads is None and not args.worker:cfg["runtime"]["threads"]=args.job_cores
    for option,key in (("infer_jobs","jobs"),("blas_threads","blas_threads")):
        if getattr(args,option,None) is not None:cfg["inference"][key]=getattr(args,option)
    inference_resources(cfg)
    from workflow_mode import resolve_mode
    if resolve_mode(args.submit) == "hpc":
        cfg["runtime"]["rscript"] = str(GRID / "f/rscript.sh")
        os.environ["WORKFLOW_RSCRIPT"]=cfg["runtime"]["rscript"]
    os.environ["GRID_RSCRIPT"]=cfg["runtime"].get("rscript","Rscript")
    if args.entry=="prsformer":cfg["prsformer"]["enabled"]=True
    if args.with_deep or args.module=="abm_deep":cfg["block_deep"]["enabled"]=True
    if getattr(args,"without_deep",False):
        if args.module=="abm_deep":raise ValueError("--without-deep cannot be used with abm_deep")
        cfg["block_deep"]["enabled"]=False
    if args.abm_deep_run:cfg["block_deep"]["device"]="cpu" if args.abm_deep_run=="cpu" else "cuda"
    if args.entry=="grid":
        if cfg["block_deep"]["device"] in {"gpu","auto"}:cfg["block_deep"]["device"]="cuda"
        if args.gpu_cores<1 or not math.isfinite(args.gpu_memory_limit_gb) or args.gpu_memory_limit_gb<=0:
            raise ValueError("GPU cores and host RAM cap must be positive")
    if cfg["inputs"].get("split_file")=="auto":
        cfg["inputs"]["split_file"]=None
    if min(cfg["cohort"].get("train_n",0),cfg["cohort"].get("label_budget",0))<0 or cfg["runtime"]["threads"]<1:raise ValueError("Invalid training budget/threads")
    return cfg

def configure_region_run(cfg, args):
    import re
    from ancestry_regions import configure_beds
    if args.worker:
        return  # BED paths and output location are frozen in the job snapshot.
    configure_beds(cfg, args)
    base = upstream_root(cfg["analysis_root"])
    name = args.run_name or ("ancestry" if cfg["regions"]["layout"]=="ancestry" else "")
    if name and (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or name in {"ldl","height","t2dm","tables","figures","main"}):
        raise ValueError("Choose an experiment name distinct from trait and report directories")
    cfg["analysis_root"] = str(base / "main" / name)


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
    if value not in (None,"auto"):return Path(value)
    scope=cfg.get("gwas_manifest_scope")
    return project_path(root, f"tmp/gwas_inputs/{scope}/manifest.tsv") if scope else project_path(root, "tmp/gwas_manifest.tsv")

def dependencies(cfg,root,stage,trait=None):
    """Readable stat records. No whole-genome hashing and no test truth in fit inputs."""
    from core2 import stamp,VERSION
    root=Path(root);inp=cfg["inputs"]
    cohort=project_path(root, "tmp/cohort/participants.tsv.gz");dev=project_path(root, f"tmp/cohort/{trait}.development.tsv.gz")
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
            from infer import select_ld_reference
            add(select_ld_reference(inp["ld_ref"],cfg["chromosomes"])[1])
    elif stage in {"infer","infer_prepare","infer_chr"}:
        settings.update(inference=cfg["inference"],chromosomes=cfg["chromosomes"])
        add(manifest_path(cfg,root),inp.get("rsid_map"))
        from infer import read_manifest
        for row in read_manifest(manifest_path(cfg,root),trait).values():add(*row["files"])
        for p in Path(inp["ld_ref"]).glob("**/*"):
            if p.is_file() and (p.name.startswith("snpinfo") or p.suffix==".hdf5"):add(p)
        add(GRID/"f/infer.py",GRID/"f/csx_parallel.py",GRID/"f/genomics.py")
        if stage in {"infer","infer_chr"}:add(trait_path(root, trait, "inference/prepared.json"))
        if stage=="infer":
            for c in cfg["chromosomes"]:add(trait_path(root, trait, f"inference/partitions/chr{c}/partition.complete.json"))
    elif stage=="score":
        settings.update(scoring=cfg["scoring"]);add(cohort,trait_path(root, trait, "inference/weights/manifest.tsv"))
        add(*sorted((trait_path(root, trait, "inference/weights")).glob("*.chr*.tsv")))
    elif stage=="fit":
        settings.update(benchmark=cfg["benchmark"],seed=cfg["cohort"]["seed"]);add(cohort,dev,trait_path(root, trait, "scores.tsv.gz"),trait_path(root, trait, "inference/weights/manifest.tsv"))
    elif stage=="disco":
        settings.update(disco=cfg["disco"],benchmark=cfg["benchmark"]);add(cohort,dev,trait_path(root, trait, "scores.tsv.gz"),trait_path(root, trait, "baseline.predictions.tsv.gz"),trait_path(root, trait, "benchmark.json"))
    elif stage=="regions":
        settings.update(scoring=cfg["scoring"],regions=cfg["regions"],consensus_build=inp["consensus_build"])
        add(cohort,trait_path(root, trait, "benchmark.json"),trait_path(root, trait, "inference/weights/manifest.tsv"))
        from ancestry_regions import bed_paths
        from consensus3 import consensus_files
        for source in bed_paths(cfg).values():
            add(*consensus_files(source))
        add(GRID/"f/ancestry_regions.py",GRID/"f/scoring_parallel.py",GRID/"f/consensus3.py")
        weights,_ = regional_weights(cfg,root,trait)
        for paths in weights.values():
            add(*paths)
    elif stage in {"direct_blocknet","blocknet_ablation"}:
        settings.update(model=cfg["direct_blocknet"],covariates=cfg["benchmark"]["covariates"],device=cfg["block_deep"]["device"])
        add(cohort,dev,trait_path(root,trait,"baseline.predictions.tsv.gz"),GRID/"f/direct_blocknet.py",GRID/"f/blocknet.py",
            *[trait_path(root,trait,"regions")/f for f in ("scores.npy","features.tsv","ids.tsv")])
        if stage=="blocknet_ablation":add(GRID/"f/blocknet_ablation.py")
    elif stage=="blocknet_twins":
        settings.update(twins=cfg["ldl_twins"],seed=cfg["cohort"]["seed"])
        add(cohort,dev,GRID/"f/blocknet_twins.py",trait_path(root,trait,"direct_blocknet/embedding.npy"),
            trait_path(root,trait,"direct_blocknet/predictions.tsv.gz"),trait_path(root,trait,"abm/abm.predictions.tsv.gz"))
    elif stage in {"block_models","abm_deep"}:
        section="block_deep" if stage=="abm_deep" else "block_models"
        settings.update(method={k:v for k,v in cfg[section].items() if k!="enabled"},seed=cfg["cohort"]["seed"])
        add(cohort,dev,trait_path(root, trait, "baseline.predictions.tsv.gz"),trait_path(root, trait, "regions/scores.npy"),trait_path(root, trait, "regions/ids.tsv"),trait_path(root, trait, "regions/features.tsv"))
    elif stage in {"abm","robustness"}:
        settings.update(abm=cfg["abm"],seed=cfg["cohort"]["seed"]);add(cohort,dev,trait_path(root, trait, "predictions.tsv.gz"))
        add(GRID/"f/abm2.py")
        if stage=="robustness" or cfg["abm"]["geometry"]!="global_prs":add(trait_path(root, trait, "regions/scores.npy"),trait_path(root, trait, "regions/ids.tsv"),trait_path(root, trait, "regions/features.tsv"))
        if stage=="robustness":add(project_path(root, f"tmp/cohort/_private/{trait}.test.tsv.gz"))
    elif stage=="evaluate":
        settings.update(evaluation=cfg["evaluation"],seed=cfg["cohort"]["seed"])
        add(cohort,project_path(root, f"tmp/cohort/_private/{trait}.test.tsv.gz"),trait_path(root, trait, "predictions.tsv.gz"),trait_path(root, trait, "abm/abm.predictions.tsv.gz"),trait_path(root, trait, "abm/model.json"))
        add(GRID/"f/driver2.py",GRID/"f/evaluate2.py", GRID/"f/manuscript_stats.py",GRID/"f/ldl_report.py",trait_path(root,trait,"benchmark.json"),
            trait_path(root,trait,"disco/model.json"),trait_path(root,trait,"disco/interpolation.tsv.gz"),
            trait_path(root,trait,"baseline.predictions.tsv.gz"))
        skipped = {}
        ready=verified_optional_predictions(cfg,root,trait,skipped=skipped)
        for path in ready.values():add(path)
        if trait=="ldl":
            from ldl_report import boundary_inputs
            peers,status=boundary_inputs(cfg,root)
            add(*peers.values())
            settings["boundary_comparison"]=status
        settings.update(optional_ready=sorted(ready),skipped_optional=skipped, include_deep=cfg["block_deep"]["enabled"])
    elif stage=="twins":
        add(cohort,trait_path(root, trait, "abm/abm.predictions.tsv.gz"),trait_path(root, trait, "abm/model.json"))
    elif stage=="architecture":add(cohort,trait_path(root, trait, "regions/scores.npy"),trait_path(root, trait, "regions/features.tsv"),trait_path(root, trait, "regions/ids.tsv"))
    elif stage=="paper":
        from manuscript import context_files
        settings.update(traits=cfg["traits"], validation_status=cfg["cohort"].get("validation_status"))
        add(GRID/"f/manuscript.py",GRID/"f/report2.py",GRID/"f/paper_compare.py",GRID/"f/comparison.html")
        if "ldl" in cfg["traits"]:
            add(GRID/"f/ldl_report.py",project_path(root,"grid.paper/ldl.blocknet_twins.tsv"),project_path(root,"grid.paper/ldl.boundaries.tsv"))
        for t in cfg["traits"]:
            add(*context_files(root,t))
            for suffix in ("contrasts.tsv","support.tsv","anchor_associations.tsv","evaluation.json","report.html"):
                add(project_path(root,f"grid.paper/{t}.{suffix}"))
            if t=="t2dm": add(project_path(root,f"grid.paper/{t}.paired_discrimination.tsv"))
            add(project_path(root, f"grid.paper/{t}.performance.tsv"),project_path(root, f"grid.paper/{t}.ancestry.tsv"),project_path(root, f"tmp/{t}/benchmark.json"),project_path(root, f"tmp/{t}/abm/model.json"))
            for path in (project_path(root, f"grid.paper/{t}.twins.html"),project_path(root, f"tmp/{t}/disco/model.json"),project_path(root, f"tmp/{t}/block_models/model.json"),project_path(root, f"tmp/{t}/deep/model.json")):
                if path.is_file():add(path)
    elif stage=="prsformer":
        from prsformer2 import PREPARED_FORMAT
        from prsformer_inputs import input_contract
        settings.update(prsformer=cfg["prsformer"],traits=cfg["traits"],prepared_format=PREPARED_FORMAT,
                        independent_inputs=input_contract(cfg,cfg["traits"]))
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
    root=Path(root)
    if stage=="pca":return [Path(cfg["inputs"]["pca_file"]),Path(cfg["inputs"]["ancestry_file"]),module_root(root,"pca")/"main/projection.json"]
    if stage=="infer_prepare":
        from csx_parallel import prepared_outputs
        return prepared_outputs(trait_path(root, trait, "inference"),cfg["chromosomes"])
    if stage=="infer_chr":
        from infer import model_manifest,parse_phi_grid
        c=cfg["chromosomes"][0];part=trait_path(root, trait, f"inference/partitions/chr{c}")
        models=model_manifest(part,parse_phi_grid(",".join(map(str,cfg["inference"]["phi_grid"]))),cfg["inference"].get("cs_fixed",False))
        return [part/"partition.complete.json",part/"provenance.json",part/"weights/manifest.tsv",part/"weights/snp_counts.tsv",
                *[part/f"weights/{m}.chr{c}.tsv" for m in models.score_name]]
    if stage=="infer":
        from infer import model_manifest,parse_phi_grid
        models=model_manifest(trait_path(root, trait, "inference"),parse_phi_grid(",".join(map(str,cfg["inference"]["phi_grid"]))),cfg["inference"].get("cs_fixed",False))
        return [trait_path(root, trait, "inference/weights/manifest.tsv"),trait_path(root, trait, "inference/provenance.json"),trait_path(root, trait, "inference/target/manifest.tsv"),
            *[trait_path(root, trait, f"inference/target/chr{c}.tsv") for c in cfg["chromosomes"]],
            *[trait_path(root, trait, f"inference/weights/{m}.chr{c}.tsv") for m in models.score_name for c in cfg["chromosomes"]]]
    if stage=="evaluate":
        return [*([project_path(root,"grid.paper/ldl.blocknet_twins.tsv")] if trait=="ldl" else []),*([project_path(root,"grid.paper/ldl.boundaries.tsv")] if trait=="ldl" else []),*[project_path(root,f"grid.paper/{trait}.{name}.tsv") for name in ("contrasts","support","anchor_associations")],project_path(root, f"grid.paper/{trait}.performance.tsv"),project_path(root, f"grid.paper/{trait}.ancestry.tsv"),project_path(root, f"grid.paper/{trait}.report.html"),project_path(root, f"grid.paper/{trait}.evaluation.json"),
                *([project_path(root, f"grid.paper/{trait}.paired_discrimination.tsv")] if trait=="t2dm" else [])]
    if stage=="paper":
        from manuscript import TABLES, FIGURES
        folder=project_path(root,"grid.paper")
        return [project_path(root,"grid.paper/"+n) for n in ("tables/table_method_comparison.tsv","comparison.html","index.html","all.performance.tsv","all.ancestry.tsv","manuscript.html")] + [
            folder/"tables"/f"table_{name}.tsv" for name in TABLES] + [
            folder/"figures"/f"{stem}.png" for stem in (*FIGURES,"fig5_method_comparison","overall_performance","support_risk_coverage")] + (
            [folder/"figures/fig4_ldl_blocknet.png"] if "ldl" in cfg["traits"] else [])
    return {
      "prepare":[project_path(root, "tmp/cohort/participants.tsv.gz"),project_path(root, "tmp/cohort/cohort.json"),*[project_path(root, f"tmp/cohort/{t}.development.tsv.gz") for t in cfg["traits"]],*[project_path(root, f"tmp/cohort/_private/{t}.test.tsv.gz") for t in cfg["traits"]]],
      "inputs":[manifest_path(cfg,root)],
      "score":[trait_path(root, trait, "scores.tsv.gz")],"fit":[trait_path(root, trait, "baseline.predictions.tsv.gz"),trait_path(root, trait, "benchmark.json"),trait_path(root, trait, "benchmark_models.joblib")],
      "disco":[trait_path(root, trait, "predictions.tsv.gz"),trait_path(root, trait, "disco/model.json")],
      "direct_blocknet":[trait_path(root,trait,"direct_blocknet")/f for f in ("model.json","model.joblib","predictions.tsv.gz","training.tsv","embedding.npy")],
      "blocknet_ablation":[trait_path(root,trait,"blocknet_ablation")/f for f in ("model.json","model.joblib","predictions.tsv.gz","training.tsv")],
      "blocknet_twins":[trait_path(root,trait,"blocknet_twins")/f for f in ("model.json","model.joblib","predictions.tsv.gz")],
      "regions":[trait_path(root, trait, "regions")/name for name in ("scores.npy","ids.tsv","features.tsv","catalog.tsv","sum_check.tsv","region_counts.tsv","covered_global.tsv.gz","unassigned_variants.tsv.gz","metadata.json")],
      "block_models":[trait_path(root, trait, "block_models.predictions.tsv.gz"),trait_path(root, trait, "block_models/model.json"),trait_path(root, trait, "block_models/model.joblib")],
      "abm_deep":[trait_path(root, trait, "deep")/name for name in ("predictions.tsv.gz","model.json","model.pt","bundle.joblib","training.tsv","tuning.tsv")],
      "abm":[trait_path(root, trait, "abm/abm.predictions.tsv.gz"),trait_path(root, trait, "abm/model.json"),trait_path(root, trait, "abm/abm.model.joblib"),
             trait_path(root, trait, "abm/twins.tsv.gz"),
             *([trait_path(root, trait, "abm/geometry.npy")] if cfg["abm"]["geometry"]=="consensus" else [])],
      "twins":[project_path(root, f"qc/{trait}.twins.tsv"),project_path(root, f"grid.paper/{trait}.twins.html")],
      "architecture":[project_path(root, f"qc/{trait}.regional_architecture.tsv")],
      "robustness":[project_path(root, f"grid.paper/{trait}.geometry_sensitivity.tsv")],
      "prsformer":[project_path(root, f"tmp/{t}/prsformer.predictions.tsv.gz") for t in cfg["traits"]]
    }[stage]

def read_manifest_weights(root, trait, chromosomes, *, family=None, phi=None):
    """Bind each model to its published, per-chromosome posterior files."""
    import pandas as pd
    folder = trait_path(root, trait, "inference/weights").resolve()
    manifest = pd.read_csv(folder/"manifest.tsv", sep="\t", dtype=str)
    required = {"score_name", "family", "phi", "population", "target_build", "weight_pattern"}
    if manifest.empty or not required <= set(manifest) or manifest.score_name.duplicated().any():
        raise ValueError(f"Missing or ambiguous posterior models in {folder}/manifest.tsv")
    if not manifest.target_build.eq("GRCh38").all():
        raise ValueError(f"Posterior weights in {folder} must be harmonized to GRCh38")
    if family is not None:
        manifest = manifest[manifest.family.eq(family) & manifest.phi.eq(str(phi))].copy()
        from core2 import POPS
        expected = {"META"} if family=="csx_meta" else set(POPS)
        if set(manifest.population)!=expected or manifest.population.duplicated().any():
            raise ValueError(f"Incomplete {family} posterior channels for {trait}, phi={phi}: {folder}")
    paths = []
    for row in manifest.itertuples():
        pattern = Path(row.weight_pattern)
        if (pattern.name != f"{row.score_name}.chr*.tsv" or
                (pattern.parent != Path(".") and pattern.parent.resolve() != folder)):
            raise ValueError(f"Posterior manifest points outside the expected model export: {pattern}; expected {folder}")
        selected = [folder/f"{row.score_name}.chr{chrom}.tsv" for chrom in chromosomes]
        for path in selected:
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing chromosome posterior for {trait}: {path}")
        paths.append(selected)
    return manifest.assign(weight_files=paths)


def regional_weights(cfg,root,trait):
    metadata = json.loads(trait_path(root,trait,"benchmark.json").read_text())
    family = "csx_meta" if cfg["regions"]["representation"]=="meta" else "csx"
    phi = metadata["geometry_phi"]
    table = read_manifest_weights(root,trait,cfg["chromosomes"],family=family,phi=phi)
    return dict(zip(table.population,table.weight_files)),phi


def regional_inputs(root,trait,cohort):
    import numpy as np
    from core2 import read_table
    from abm2 import IndexedFeatures
    folder=project_path(root, f"tmp/{trait}/regions")
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
    root=Path(root);o=cfg["abm"]
    cohort=identified(read_table(project_path(root, "tmp/cohort/participants.tsv.gz")));mode=geometry or o["geometry"]
    features,maps=(None,None) if mode=="global_prs" else regional_inputs(root,trait,cohort)
    return fit_abm(cohort,read_table(project_path(root, f"tmp/cohort/{trait}.development.tsv.gz")),read_table(trait_path(root, trait, "predictions.tsv.gz")),output or trait_path(root, trait, "abm"),
        trait=trait,features=features,feature_maps=maps,feature_manifest=None if mode=="global_prs" else read_table(trait_path(root, trait, "regions/features.tsv")),
        geometry=mode,dimensions=o["dimensions"],k_values=o["k_values"],
        radius_k=o["radius_k"],radius_quantile=o["radius_quantile"],seed=cfg["cohort"]["seed"],workers=cfg["runtime"]["threads"],
        batch_size=o["batch_size"],min_count=o["min_count"],min_group_ess=o["min_group_ess"],shrinkage=o["shrinkage"],
        donor_batch=o.get("donor_batch",4096),evidence_queries=o.get("evidence_queries",200),evidence_neighbors=o.get("evidence_neighbors",16))

def fit_regional(cfg,root,trait,deep=False):
    from core2 import read_table,identified
    root=Path(root)
    cohort=identified(read_table(project_path(root, "tmp/cohort/participants.tsv.gz")))
    features,_=regional_inputs(root,trait,cohort)
    common=(cohort,read_table(project_path(root, f"tmp/cohort/{trait}.development.tsv.gz")),
            read_table(trait_path(root, trait, "baseline.predictions.tsv.gz")),features,read_table(trait_path(root, trait, "regions/features.tsv")))
    options=dict(trait=trait,seed=cfg["cohort"]["seed"],threads=cfg["runtime"]["threads"])
    if deep:
        from block_deep import fit_block_deep
        return fit_block_deep(*common,trait_path(root, trait, "deep"),config=cfg["block_deep"],**options)
    from block_models import fit_block_models
    return fit_block_models(*common,trait_path(root, trait, "block_models"),config=cfg["block_models"],rscript=cfg["runtime"]["rscript"],**options)

def verified_optional_predictions(cfg,root,trait,*,load=False,skipped=None,stages=None):
    """Read completed optional exports without waiting for their producer."""
    from core2 import VERSION,stamp,read_table,identified
    import numpy as np
    import pandas as pd
    root=Path(root); result={}
    specs=(("block_models","block_models.predictions.tsv.gz"),("abm_deep","deep/predictions.tsv.gz"))
    if trait=="ldl" and cfg["block_deep"]["enabled"] and cfg.get("direct_blocknet"):
        specs += tuple((name,name+"/predictions.tsv.gz") for name in ("direct_blocknet","blocknet_ablation","blocknet_twins"))
    roster = None
    for stage, relative in specs:
        if stages is not None and stage not in stages:continue
        path = trait_path(root,trait,relative)
        try:
            if stage=="abm_deep" and not cfg["block_deep"]["enabled"]:
                raise ValueError("omitted by configuration")
            owner = root
            marker=project_path(owner,"tmp/state")/f"{trait}.{stage}.json"
            if not path.is_file():
                raise ValueError("final predictions not yet published")
            if not marker.is_file():
                raise ValueError("final producer receipt not yet published")
            previous=json.loads(marker.read_text())
            record=previous["inputs"]
            producer=dict(cfg)
            current=dependencies(producer,owner,stage,trait)
            expected=previous["outputs"].get(str(path))
            if (record["settings"].get("version")!=VERSION or record.get("files")!=current["files"]
                    or expected!=stamp(path)):
                raise ValueError("producer inputs or final output do not match receipt")
            if stage in {"direct_blocknet","blocknet_ablation","blocknet_twins"}:
                if ({k:v for k,v in current["settings"].items() if k!="device"}!={k:v for k,v in record["settings"].items() if k!="device"}) or any(
                        not Path(p).is_file() or stamp(Path(p))!=sig for p,sig in previous["outputs"].items()):
                    raise ValueError("LDL producer settings or artifacts differ from the completed fit")
            if not load:
                result[stage]=path
                continue
            if roster is None:
                roster=identified(read_table(project_path(root,"tmp/cohort/participants.tsv.gz"),columns=["id","role"]))
            frame=identified(read_table(path),f"{stage} predictions")
            if set(frame.id)!=set(roster.id):raise ValueError("prediction roster differs")
            if stamp(path)!=expected: raise ValueError("prediction file changed during read")
            result[stage]=frame
        except (OSError,EOFError,ValueError,KeyError,TypeError) as exc:
            if skipped is None:
                if path.is_file():raise
            else:
                skipped[stage]=str(exc)
    return result

def evaluate(cfg,root,trait,abm_path=None,output=None):
    from core2 import read_table,identified
    from evaluate2 import evaluate_predictions
    root=Path(root)
    benchmark=read_table(trait_path(root, trait, "predictions.tsv.gz"))
    abm=read_table(abm_path or trait_path(root, trait, "abm/abm.predictions.tsv.gz"))
    skipped={}
    optional=verified_optional_predictions(cfg,root,trait,load=True,skipped=skipped)
    if output is None and trait=="ldl" and cfg["block_deep"]["enabled"] and cfg.get("direct_blocknet"):
        missing={name:skipped.get(name,"missing") for name in ("direct_blocknet","blocknet_ablation","blocknet_twins") if name not in optional}
        if missing:
            raise ValueError(f"LDL report needs completed BlockNet, ablation and twin fits for this BED layout: {missing}")
    for method,reason in skipped.items():print(f"[evaluate] omitted optional {method}: {reason}",flush=True)
    if "block_models" in optional:
        benchmark=benchmark.merge(optional["block_models"],on="id",how="left",validate="one_to_one")
    for stage in ("direct_blocknet","blocknet_ablation"):
        if stage in optional:
            benchmark=benchmark.merge(optional[stage],on="id",how="left",validate="one_to_one")
    if "abm_deep" in optional:
        deep=optional["abm_deep"]
        benchmark=benchmark.merge(deep[["id","abm.deep.prediction"]],on="id",how="left",validate="one_to_one")
        columns=[c for c in deep if c.startswith("abm.deep.selected_")]
        masks=deep[["id",*columns]].rename(columns={c:c.replace("abm.deep.selected_","selected_deep_") for c in columns})
        abm=abm.merge(masks,on="id",how="left",validate="one_to_one")
    from manuscript_stats import manuscript_evaluation
    from core2 import atomic_json
    import numpy as np
    from evaluate2 import aligned
    cohort=identified(read_table(project_path(root,"tmp/cohort/participants.tsv.gz")))
    truth=read_table(project_path(root,f"tmp/cohort/_private/{trait}.test.tsv.gz"))
    folder=Path(output) if output is not None else project_path(root,f"grid.paper/traits/{trait}")
    reps=cfg["evaluation"]["bootstrap"]; seed=cfg["cohort"]["seed"]
    disco_path=trait_path(root,trait,"disco/model.json")
    disco=json.loads(disco_path.read_text()) if disco_path.is_file() else {}
    has_disco="disco.prediction" in benchmark
    if disco.get("status") in {"complete","unavailable"} and (disco["status"]=="complete")!=has_disco:
        raise ValueError(f"{trait}: Disco model status and predictions disagree; restore matching 3.disco outputs")
    source=trait_path(root,trait,"baseline.predictions.tsv.gz")
    if source.is_file():
        cols=["id","csx.population.prediction"]
        baseline=aligned(read_table(source,columns=cols),benchmark.id,"CSx baseline")
        if not np.allclose(baseline[cols[1]].to_numpy(float),benchmark[cols[1]].to_numpy(float),rtol=1e-12,atol=1e-12):
            raise ValueError(f"{trait}: CSx baseline and combined Disco snapshot differ; restore matching outputs")
    abm_model=Path(abm_path).parent/"model.json" if abm_path else trait_path(root,trait,"abm/model.json")
    if abm_model.is_file() and json.loads(abm_model.read_text()).get("alpha")==0:
        same=aligned(benchmark[["id","csx.population.prediction"]],abm.id,"ABM baseline check")
        if not np.allclose(abm["abm.prediction"].to_numpy(float),same["csx.population.prediction"].to_numpy(float),
                           rtol=1e-12,atol=1e-12):
            raise ValueError(f"{trait}: ABM alpha is zero but frozen predictions differ from CSx; check snapshot versions")
    result=evaluate_predictions(cohort,truth,benchmark,abm,folder,trait=trait,bootstrap=reps,seed=seed,
                                workers=cfg["runtime"]["threads"])
    if output is None:
        if trait=="ldl":
            from ldl_report import twin_audit,boundary_audit
            boundary_audit(cfg,root,cohort,truth,folder)
            twin_audit(cohort,truth,benchmark,abm,optional.get("blocknet_twins"),folder,reps,seed,cfg["runtime"]["threads"],trait)
        model_path=trait_path(root,trait,"benchmark.json")
        frozen=json.loads(model_path.read_text()) if model_path.is_file() else {}
        manuscript_evaluation(trait,folder,cohort,truth,benchmark,abm,trait_path(root,trait,"disco/interpolation.tsv.gz"),
                              frozen.get("covariates",cfg["benchmark"]["covariates"]),reps,seed)
        result[2]["optional_methods_included"]=list(optional)
        result[2]["optional_methods_omitted"]=skipped
        result[2]["disco_status"]=disco.get("status","present" if has_disco else "unavailable")
        result[2]["disco_reason"]=disco.get("reason","")
        atomic_json(result[2],folder/f"{trait}.evaluation.json")
    return result

def execute(cfg,root,stage,trait,args):
    from core2 import atomic_json,atomic_table,read_table
    root=Path(root);inp=cfg["inputs"]
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
        command=[sys.executable,str(GRID/"f/infer.py"),"--manifest",str(manifest_path(cfg,root)),"--trait",trait,"--ld-ref",inp["ld_ref"],"--genotype-dir",inp["genotype_dir"],"--output",str(trait_path(root, trait, "inference")),"--chrs",",".join(map(str,cfg["chromosomes"])),"--phi-grid",",".join(map(str,o["phi_grid"])),"--jobs",str(jobs),"--threads",str(blas),"--backend",o["backend"],"--n-iter",str(o["n_iter"]),"--n-burnin",str(o["n_burnin"]),"--thin",str(o["thin"]),"--seed",str(o["seed"])]
        if o.get("cs_fixed"):command.append("--cs-fixed")
        if inp.get("rsid_map"):command += ["--rsid-map",inp["rsid_map"]]
        if args.replace:command.append("--replace")
        if stage=="infer_prepare":command.append("--prepare-only")
        elif stage=="infer_chr":
            command[command.index("--output")+1]=str(trait_path(root, trait, f"inference/partitions/chr{cfg['chromosomes'][0]}"))
            command += ["--prepared-inputs",str(trait_path(root, trait, "inference"))]
        else:command.append("--merge-only")
        subprocess.run(command,check=True)
    elif stage=="score":
        from scoring_parallel import score_models
        m=read_manifest_weights(root,trait,cfg["chromosomes"])
        score_models(project_path(root, "tmp/cohort/participants.tsv.gz"),dict(zip(m.score_name,m.weight_files)),inp["genotype_dir"],trait_path(root, trait, "scores.tsv.gz"),
            common_variants=cfg["scoring"]["common_variants"],plink=cfg["runtime"]["plink"],threads=cfg["runtime"]["threads"],memory_mb=cfg["runtime"]["plink_memory_mb"],batch_size=cfg["scoring"]["batch_size"])
    elif stage=="fit":
        from benchmark2 import fit_models;fit_models(cfg,root,trait)
    elif stage=="disco":
        from benchmark2 import fit_disco;fit_disco(cfg,root,trait)
    elif stage=="regions":
        from ancestry_regions import score
        weights,phi = regional_weights(cfg,root,trait)
        score(cfg,project_path(root,"tmp/cohort/participants.tsv.gz"),weights,
              trait_path(root,trait,"regions"),phi)
    elif stage in {"direct_blocknet","blocknet_ablation"}:
        from direct_blocknet import fit
        for method in (("direct.global_mlp","direct.regional_ridge") if stage=="blocknet_ablation" else ("direct.linear","direct.blocknet")):
            fit(cfg,root,trait,method)
    elif stage=="blocknet_twins":
        from blocknet_twins import fit
        fit(cfg,root,trait)
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
        for mode in ("consensus","global_prs","regional","unit_block"):
            location=trait_path(root, trait, "geometry_sensitivity")/mode
            run_abm(cfg,root,trait,mode,location)
            performance,_,_=evaluate(cfg,root,trait,location/"abm.predictions.tsv.gz",location/"evaluation")
            rows.append(performance.assign(geometry=mode,interpretation="predeclared exploratory sensitivity; no test selection of primary geometry"))
        atomic_table(pd.concat(rows,ignore_index=True),project_path(root, f"grid.paper/{trait}.geometry_sensitivity.tsv"))
    elif stage=="paper":
        from report2 import final_report
        from manuscript import publish
        final_report(cfg,root)
        publish(cfg,root)
    elif stage=="prsformer":
        from copy import deepcopy
        from prsformer2 import run
        worker_cfg=deepcopy(cfg)
        if args.replace and args.prsformer_stage in {"all","train"}:
            worker_cfg["prsformer"]["force_retrain"]=True
        if args.replace and args.prsformer_stage in {"all","prepare"}:
            from prsformer2 import _section
            worker_cfg["prsformer"]["prepare"]={**_section(worker_cfg["prsformer"],"prepare"),"replace":True}
        run(worker_cfg,root,cfg["traits"],stage=args.prsformer_stage)
    else:raise ValueError(stage)

def run_stage(cfg,root,stage,trait,args):
    from core2 import atomic_json,stamp
    root=Path(root);name=f"{trait}.{stage}" if trait else stage
    if stage=="infer_chr":name=f"{trait}.chr{cfg['chromosomes'][0]}.infer"
    if stage=="prsformer" and args.prsformer_stage in {"prepare","plan"}:name+=f".{args.prsformer_stage}"
    if stage=="inputs" and cfg.get("gwas_manifest_scope"):name+="."+cfg["gwas_manifest_scope"]
    if args.entry in {"pca","csx","disco","prsformer"}:
        state = module_root(root, args.entry) / "qc" / "state" if stage!="prsformer" else module_root(root, args.entry)/"qc"
    elif stage in {"paper","evaluate"}:
        state = project_path(root,"grid.paper/logs")/"state"
    else:
        state = project_path(root,"tmp/state")
    state.mkdir(parents=True,exist_ok=True)
    with (state/f"{name}.lock").open("a+") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        signature=dependencies(cfg,root,stage,trait);marker=state/f"{name}.json"
        outputs=stage_outputs(cfg,root,stage,trait)
        # Inference workers validate their detailed preparation/chain receipts
        # themselves, including raw posteriors and code identities.
        if marker.is_file() and not args.replace and stage not in {"infer","infer_prepare","infer_chr"}:
            previous=json.loads(marker.read_text())
            if stage in {"prepare","inputs"} and previous.get("inputs")!=signature:
                from copy import deepcopy
                old_signature=deepcopy(previous.get("inputs",{}))
                new_signature=deepcopy(signature)
                old_traits=old_signature.get("settings",{}).pop("traits",[])
                new_traits=new_signature["settings"].pop("traits",[])
                # A subset can reuse the original shared manifest. New traits
                # get a separate manifest in plan_jobs, never overwrite it.
                if stage=="inputs" and set(new_traits)<=set(old_traits):
                    expanded=deepcopy(cfg);expanded["traits"]=old_traits
                    new_signature=dependencies(expanded,root,stage)
                    new_signature["settings"].pop("traits",None)
                intact=all(Path(p).is_file() and stamp(Path(p))==record
                           for p,record in previous.get("outputs",{}).items())
                if old_signature==new_signature and intact:
                    missing=[t for t in new_traits if t not in old_traits]
                    if stage=="prepare" and missing:
                        from cohort2 import extend_traits
                        extend_traits(cfg,root,missing)
                        expanded=deepcopy(cfg);expanded["traits"]=old_traits+missing
                        signature=dependencies(expanded,root,stage)
                        outputs=stage_outputs(expanded,root,stage)
                        if not all(p.is_file() for p in outputs):raise RuntimeError("prepare: expected output missing")
                        atomic_json({"inputs":signature,"outputs":{str(p):stamp(p) for p in outputs}},marker)
                        print(f"[extend] prepare: {','.join(missing)}; existing cohort and labels preserved",flush=True)
                        return
                    if not missing:
                        print(f"[reuse] {name}: requested trait subset",flush=True);return
            same_outputs=all(p.is_file() and previous.get("outputs",{}).get(str(p))==stamp(p) for p in outputs)
            if previous.get("inputs")==signature and same_outputs:
                print(f"[reuse] {name}",flush=True);return
            if previous.get("inputs")!=signature and stage not in {"evaluate","paper"}:
                raise ValueError(f"{name}: inputs/settings changed; use --replace TRUE for this stage, or a new --analysis-root")
            raise ValueError(f"{name}: result receipt does not match current inputs/outputs; rerun this stage with --replace TRUE or use a new --analysis-root")
        if stage=="prsformer" and args.prsformer_stage in {"train","predict","all"}:
            marker.unlink(missing_ok=True)
        print(f"[run] {name}",flush=True)
        execute(cfg,root,stage,trait,args)
        if stage=="prsformer" and args.prsformer_stage in {"plan","prepare"}:return
        if not all(p.is_file() for p in outputs):raise RuntimeError(f"{name}:expectedoutputmissing")
        atomic_json({"inputs":signature,"outputs":{str(p):stamp(p) for p in outputs}},marker)

def plan_jobs(cfg,args):
    from core2 import atomic_json
    from workflow_jobs import Job,write_plan,run_plan
    root=Path(cfg["analysis_root"]).resolve();root.mkdir(parents=True,exist_ok=True)
    if args.entry=="csx" and cfg["inputs"].get("gwas_manifest") in (None,"auto"):
        marker=module_root(root,"csx")/"qc/state/inputs.json"
        if marker.is_file() and not cfg.get("gwas_manifest_scope"):
            previous=json.loads(marker.read_text()).get("inputs",{}).get("settings",{})
            if not set(cfg["traits"])<=set(previous.get("traits",[])):
                # Freeze the selection in the plan: per-trait workers must use
                # this same path even after --traits narrows their selection.
                cfg["gwas_manifest_scope"]="-".join(sorted(set(cfg["traits"])))
    # Every plan has a separate immutable config snapshot, shared by its workers.
    from workflow_layout import new_run_dir
    import shutil
    paper = args.entry=="grid" and args.module in {"paper-run","paper","evaluate"}
    component = module_root(root,args.entry)
    if paper:
        workspace=project_path(root,"grid.paper/logs")/"plans"
    elif args.entry=="grid":
        workspace=project_path(root,"tmp/plans")
    else:
        workspace=component/"tmp"/"plans"
    snapshot=new_run_dir(workspace)
    code=GRID
    if args.entry=="prsformer":
        # Pending jobs and predict subprocesses must use the same scientific code as prepare.
        code=snapshot/"grid"
        shutil.copytree(GRID/"f",code/"f",ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        import importlib.util
        shared_copy = snapshot / "0f"
        shutil.copytree(SHARED, shared_copy, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        spec = importlib.util.spec_from_file_location("grid_snapshot_jobs", shared_copy/"workflow_jobs.py")
        scheduler = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = scheduler
        spec.loader.exec_module(scheduler)
        Job, write_plan, run_plan = scheduler.Job, scheduler.write_plan, scheduler.run_plan
    config=snapshot/"config.json"
    atomic_json(cfg,config)
    jobs=[];last=[]
    worker_python=cfg.get("worker_python",sys.executable)
    def prsformer_device(action):
        from prsformer2 import _options
        return _options(cfg, "predict" if action=="predict" else "train").get("device","cuda")

    def add(stage,trait=None,depends=(),chromosome=None,action=None):
        key=f"{trait}.{stage}" if trait else stage
        action=action or args.prsformer_stage
        if stage=="prsformer":key=f"{trait}.{action}"
        if chromosome is not None:key=f"{trait}.chr{chromosome}.infer"
        argv=[worker_python,str(code/f"f/{args.entry}_driver.py"),stage,"--grch","38","--config",str(config),"--worker","--replace",str(args.replace)]
        if trait:argv += ["--traits",trait]
        if chromosome is not None:argv += ["--worker-chromosome",str(chromosome)]
        if stage=="prsformer":argv += ["--prsformer-stage",action]
        resources={"queue":args.queue,"cores":args.job_cores,"memory_gib":args.job_memory_gb,"hours":args.job_hours}
        blas=cfg["runtime"]["threads"]
        if args.entry=="csx" and stage in {"infer_prepare","infer_chr","infer","inputs"}:
            _,blas=inference_resources(cfg)
            # Preparation/aggregation is serial; chromosome BLAS uses its own budget.
            resources["cores"] = min(args.job_cores, blas) if stage=="infer_chr" else 1
            resources.update(memory_gib=args.infer_memory_gb)
            if stage!="infer_chr":blas=1
        if stage=="prsformer" and action in {"all","train","predict"} and prsformer_device(action) != "cpu":
            resources.update(queue=args.gpu_queue or cfg["prsformer"].get("gpu_queue") or os.environ.get("HPC_GPU_QUEUE","4a100-40"),
                             cores=args.gpu_cores,gpus=1)
            blas=min(cfg["runtime"]["threads"],args.gpu_cores)
            argv += ["--threads",str(blas)]
        if stage in {"abm_deep","direct_blocknet","blocknet_ablation"} and cfg["block_deep"]["device"] in {"cuda","gpu"}:
            resources.update(queue=args.gpu_queue or cfg["runtime"].get("gpu_queue",os.environ.get("HPC_GPU_QUEUE","4a100-40")),
                             cores=args.gpu_cores,memory_gib=min(args.job_memory_gb,args.gpu_memory_limit_gb),gpus=1)
            blas=min(cfg["runtime"]["threads"],args.gpu_cores)
            argv += ["--threads",str(blas)]
        job_name = "grid.pca" if args.entry=="pca" and stage=="pca" else (f"grid.disco.{trait}" if args.entry=="disco" and stage=="disco" else (f"prsformer.{trait}.{action}" if stage=="prsformer" and trait else None))
        if paper and stage == "evaluate" and trait:
            job_name = f"grid.paper.{trait}"
        elif paper and stage == "paper":
            job_name = "grid.paper.paper"
        if args.entry=="grid" and root.parent.name=="main":
            job_name=(job_name or key)+"."+root.name
        jobs.append(Job(key=key,argv=argv,depends=list(depends),resources=resources,cwd=str(code),
                        job_name=job_name,env={"GRID_MEMORY_CAP_GB":str(resources["memory_gib"]),"PYTHONDONTWRITEBYTECODE":"1","WORKFLOW_TMP_ROOT":str((workspace.parent if not paper else workspace)/"workers"),"LD_LIBRARY_PATH":os.environ.get("LD_LIBRARY_PATH",""),**{var:str(blas) for var in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS")}}))
        return key
    if args.entry=="grid":
        from grid_workflow import build_jobs
        build_jobs(cfg,args,add)
    else:
        from standalone_workflows import build_jobs
        build_jobs(cfg,args,add)
    os.environ["WORKFLOW_PYTHON"]=worker_python
    log_root = (project_path(root,"grid.paper/logs") if args.entry == "grid" and args.module in {"paper-run","paper","evaluate"}
                else project_path(root,"jobs") if args.entry=="grid" else component/"jobs")
    from workflow_mode import resolve_mode
    if args.entry in {"pca","csx","disco"} and resolve_mode(args.submit)=="hpc":
        from cluster import numbered
        return numbered(jobs,component,args,worker_python)
    plan=write_plan(jobs,log_root,mode=args.submit,name=args.entry,common_env={"WORKFLOW_PYTHON":worker_python})
    if not paper:
        plan_dir=Path(plan).parent
        log_parent=project_path(root,"logs") if args.entry=="grid" else component/"logs"
        log_dir=log_parent/plan_dir.name
        log_dir.mkdir(parents=True,exist_ok=True)
        (plan_dir/"logs").rmdir()  # empty, freshly created by write_plan()
        (plan_dir/"logs").symlink_to(log_dir.resolve(),target_is_directory=True)
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
        trait=None if stage in {"pca","prepare","inputs","paper"} else cfg["traits"][0]
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
