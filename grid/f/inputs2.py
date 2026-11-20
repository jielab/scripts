"""Discover the existing common GWAS layout and record an explicit sample-size manifest."""
from pathlib import Path
import numpy as np
import pandas as pd
from core2 import POPS,read_table,atomic_table,atomic_json

def find_gwas(directory,trait,pop):
    tags=[f"{trait}.{pop}"]+(["t2dm.AFA"] if (trait,pop)==("t2dm","AFR") else [])
    for tag in tags:
        for file in (Path(directory)/tag/"gwas"/f"{tag}.gz",Path(directory)/f"{tag}.gz",Path(directory)/f"{tag}.tsv.gz"):
            if file.is_file():return file.resolve()
    raise FileNotFoundError(f"No GWAS for {trait}.{pop}; setinputs.gwas_manifestexplicitly")

def sample_size(path,trait,snps,effective_n_column=None):
    header=pd.read_csv(path,sep=r"\s+",nrows=0).columns.tolist()
    lower={c.lower().lstrip("#"):c for c in header}
    def col(*names):return next((lower[n] for n in names if n in lower),None)
    sid=col("snp","rsid","rs_id","markername","rsids")
    neff=col("neff","n_eff","effective_n","n_effective")
    if effective_n_column:
        if effective_n_column not in header:
            raise ValueError(f"{path}: declared effective-N column {effective_n_column!r} is missing")
        neff=effective_n_column
    cases=col("ncase","n_case","n_cases","ncases","cases")
    controls=col("ncontrol","n_control","n_controls","nctrl","n_ctrl","ncontrols","controls")
    total=col("n","n_total","n_samples","samplesize","totaln")
    if neff:cols=[neff];rule=f"median SNP effective N (column {neff}; " + ("explicit verified mapping" if effective_n_column else "header") + ")"
    elif trait=="t2dm" and cases and controls:cols=[cases,controls];rule="median SNP 4/(1/Ncase+1/Ncontrol)"
    elif trait!="t2dm" and total:cols=[total];rule="median SNP N"
    else:raise ValueError(f"{path}: supplyexplicitn_gwas; T2DMrequiresNefforNcase+Ncontrol, notunverifiedtotalN")
    if not sid:raise ValueError(f"{path}: automaticNrequiresSNP/rsID; provideexplicitn_gwasinmanifest")
    values=[]
    for chunk in pd.read_csv(path,sep=r"\s+",usecols=[sid,*cols],chunksize=250000):
        chunk=chunk[chunk[sid].astype(str).isin(snps)]
        if len(cols)==1:n=pd.to_numeric(chunk[cols[0]],errors="coerce").to_numpy(float)
        else:
            ca=pd.to_numeric(chunk[cases],errors="coerce").to_numpy(float);co=pd.to_numeric(chunk[controls],errors="coerce").to_numpy(float)
            with np.errstate(divide="ignore",invalid="ignore"):n=4/(1/ca+1/co)
            n[(ca<=0)|(co<=0)]=np.nan
        values.append(n[np.isfinite(n)&(n>0)])
    values=np.concatenate(values)
    if not len(values):raise ValueError(f"{path}: no valid N among reference HM3 SNPs")
    return float(np.median(values)),{"rule":rule,"snp_count":len(values),"min":float(values.min()),"max":float(values.max())}

def ensure_manifest(cfg,root,traits):
    declared=cfg["inputs"].get("gwas_manifest","auto")
    if declared and declared!="auto":return str(Path(declared).resolve())
    root=Path(root);destination=root/"tmp/gwas_manifest.tsv"
    snpinfo=Path(cfg["inputs"]["ld_ref"])/"snpinfo_mult_1kg_hm3"
    reference=read_table(snpinfo,columns=["SNP"])
    if "SNP" not in reference:raise ValueError("LD snpinfo must declare SNP rsIDs")
    snps=set(reference.SNP)
    records=[];provenance=[]
    for trait in traits:
        for pop in POPS:
            path=find_gwas(cfg["inputs"]["gwas_dir"],trait,pop)
            override=cfg["inference"].get("n_gwas",{}).get(trait,{}).get(pop)
            if override:
                n=float(override);detail={"rule":"explicit configuration override"}
            else:n,detail=sample_size(path,trait,snps,cfg["inference"].get("effective_n_columns",{}).get(trait))
            if not np.isfinite(n) or n<=0:raise ValueError(f"InvalidNfor{trait}.{pop}")
            records.append({"trait":trait,"population":pop,"path":str(path),"n_gwas":n,
                "sample_overlap":cfg["inference"].get("discovery_target_overlap","unknown")})
            provenance.append({"trait":trait,"population":pop,"path":str(path),"n_gwas":n,**detail})
    atomic_table(pd.DataFrame(records),destination)
    atomic_json({"sample_size_resolution":provenance,"overlap_status":cfg["inference"].get("discovery_target_overlap","unknown"),
        "note":"GWAS-target and cross-GWAS sample overlap cannot be inferred from summary files; retain source-study documentation"},root/"qc/gwas_inputs.json")
    return str(destination)
