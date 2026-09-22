#!/usr/bin/env python3
"""Project new people with frozen models; outcome columns are not required."""
from pathlib import Path
import argparse
import json
import joblib
import numpy as np
import pandas as pd
import torch
from common import words, log, STAGES
from io_data import read_table,numeric,map_metabolites
from representation import load_encoder,encode
from prediction import make_blocks
from neural import NeuralModel,random_neighbors

def project_people(root,phenotype,omics,batch=512,met_input="named",r_bin="Rscript"):
    root=Path(root)
    manifest=json.loads((root/"manifest.json").read_text())
    cfg=manifest["config"]
    from panome import completed
    for stage in STAGES[:5]:
        if not completed(root/stage,manifest["signature"]):
            raise ValueError(f"Projection requires completed {stage}")
    idcol=cfg["id_col"]
    cols=list(dict.fromkeys([idcol]+words(cfg["covariates"])+words(cfg["residualize"])+
                           ([cfg["group_col"]] if cfg["group_col"] else [])))
    p=read_table(phenotype,idcol,cols,r_bin)
    raw=read_table(omics,idcol,r_bin=r_bin)
    if not len(raw):
        raise ValueError("No new people supplied")
    if cfg["biom"]=="prot":
        raw.columns=[c if c==idcol else str(c).upper() for c in raw]
        if raw.columns.duplicated().any():
            raise ValueError("Protein names collide after uppercasing")
    elif met_input=="raw":
        raw,_=map_metabolites(raw,cfg["met_map"],idcol)
    if not raw[idcol].isin(p[idcol]).all():
        raise ValueError("Some molecular IDs lack baseline metadata")
    p=p.set_index(idcol).loc[raw[idcol]].reset_index()
    pre=joblib.load(root/"s2_preprocess/preprocessor.joblib")
    original=(root/"s1_prepare/features.txt").read_text().splitlines()
    needed=[original[j] for j in pre.keep]
    absent=set(needed)-set(raw)
    if absent:
        raise ValueError("Missing required assays: "+", ".join(sorted(absent)[:20]))
    values=numeric(raw,needed,"projection")
    if np.any(np.mean(~np.isfinite(values),axis=1)>=cfg["sample_missing"]):
        raise ValueError("New people fail the training-defined sample-missingness threshold")
    full=np.full((len(values),len(original)),np.nan,dtype="float32")
    full[:,pre.keep]=values
    x,observed=pre.transform(full,p)
    torch.set_num_threads(cfg["cores"])
    encoder=load_encoder(root/"s3_representation/ae.pt")
    z=encode(encoder,x,observed,batch=batch)
    pca=joblib.load(root/"s3_representation/pca.joblib").transform(x)
    atlas=joblib.load(root/"s4_graph/atlas.joblib")
    if cfg["group_col"] and (p[cfg["group_col"]].isna().any() or
                            p[cfg["group_col"]].astype(str).str.strip().eq("").any()):
        raise ValueError("New people require complete family/group IDs")
    groups=p[cfg["group_col"]].to_numpy(str) if cfg["group_col"] else None
    graph=atlas.project(z,p[idcol].to_numpy(str),groups)
    clinical=joblib.load(root/"s5_predict/clinical_preprocessor.joblib")
    c=clinical.transform(p)
    pw=joblib.load(root/"s5_predict/pwas_score.joblib")
    score=x[:,pw["indices"]]@pw["beta"]
    blocks=make_blocks(c,x,z,pca,graph,score,varying=True)
    output=p[[idcol]].copy()
    log("AUDIT","projection categories",str(clinical.unknown_categories(p)))
    output["state"]=graph["state"]+1
    output["novel"]=graph["novelty"]
    output["observed_fraction"]=observed.mean(1)
    output["effective_neighbors"]=graph["effective_neighbors"]
    for j in range(z.shape[1]):
        output[f"AE{j+1}"]=z[:,j]
    for j in range(graph["soft"].shape[1]):
        output[f"state_weight_{j+1}"]=graph["soft"][:,j]
    artifacts=json.loads((root/"s5_predict/models.json").read_text())
    u=bank=random=None
    if any(row["type"]=="neural" for row in artifacts):
        token_encoder=load_encoder(root/"s3_representation/transformer.pt")
        u=encode(token_encoder,x,observed,batch=batch)
        u=joblib.load(root/"s5_predict/transformer_scaler.joblib").transform(u).astype("float32")
        bank=np.load(root/"s5_predict/reference_bank.npy")
        random=random_neighbors(atlas,p[idcol].to_numpy(str),groups,cfg["seed"])
    train_time=None
    if cfg["outcome_type"]=="survival":
        cohort=pd.read_csv(root/"s2_preprocess/cohort.csv",usecols=["split","time"])
        train_time=float(cohort.loc[cohort.split=="train","time"].max())
    for row in artifacts:
        name=row["name"]
        path=root/"s5_predict"/row["artifact"]
        if row["type"]=="classical":
            fit=joblib.load(path)
            b=blocks[row["inputs"]]
            predicted=fit.predict(b)
            survival=lambda h,fit=fit,b=b:fit.survival_at(b,h)
        else:
            fit=NeuralModel.load(path)
            inputs=x if row["inputs"]=="molecular" else u
            neighbors=random if row["inputs"]=="random" else graph["neighbors"]
            predicted=fit.predict(inputs,c,neighbors,bank,batch=batch)
            survival=lambda h,fit=fit,s=predicted:fit.survival(s,h)
        output[f"{name}_prediction"]=predicted
        if cfg["outcome_type"]=="survival":
            for horizon in cfg["horizons"]:
                if horizon<train_time:
                    try:
                        output[f"{name}_net_risk_{horizon:g}y"]=1-survival(horizon)
                    except ValueError:
                        pass
    return output

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-dir",required=True,type=Path)
    ap.add_argument("--phenotype",required=True)
    ap.add_argument("--omics",required=True)
    ap.add_argument("--output",required=True,type=Path)
    ap.add_argument("--batch-size",type=int,default=512)
    ap.add_argument("--met-input",choices=["named","raw"],default="named")
    ap.add_argument("--r-bin",default="Rscript")
    a=ap.parse_args()
    if a.batch_size<1:
        raise ValueError("batch-size must be positive")
    if a.output.exists():
        raise FileExistsError(f"Output already exists: {a.output}")
    result=project_people(a.run_dir,a.phenotype,a.omics,a.batch_size,a.met_input,a.r_bin)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    result.to_csv(a.output,index=False)
    log("DONE","frozen projection",f"N={len(result)}; {a.output}")

if __name__=="__main__":
    main()
