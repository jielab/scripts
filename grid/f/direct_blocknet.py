"""Direct trait predictors sharing the same TRAIN, tuning families and calibration families."""
import hashlib

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge

from abm2 import RAW_COLUMNS, check_roster, development_outcomes, regional_feature_schema
from core2 import NumericDesign, aligned, atomic_json, atomic_table, read_table
from output_layout import project_path, trait_path


def fit(cfg, root, trait, method):
    binary = trait == "t2dm"
    folder=trait_path(root,trait,method)
    folder.mkdir(parents=True,exist_ok=True)
    region=trait_path(root,trait,"regions")
    cohort=read_table(project_path(root,"tmp/cohort/participants.tsv.gz"))
    check_roster(cohort)
    ids=read_table(region/"ids.tsv").id
    ix=pd.Index(ids).get_indexer(cohort.id)
    array=np.load(region/"scores.npy",mmap_mode="r")[:,:,None]
    if np.any(ix<0) or len(array)!=len(ids):
        raise ValueError("Regional matrix differs from its participant manifest")
    manifest=read_table(region/"features.tsv")
    suffix="37" if "chrom37" in manifest else "38"
    start=manifest["start"+suffix].to_numpy(float)
    end=manifest["end"+suffix].to_numpy(float)
    tokens=pd.DataFrame(dict(CHR=manifest["chrom"+suffix].to_numpy(int),
                             source_population=manifest["map"].replace({"META":"ALL"}),
                             midpoint=(start+end)/2,length=end-start))
    scores=aligned(read_table(trait_path(root,trait,"baseline.predictions.tsv.gz")),cohort.id)
    for col in RAW_COLUMNS:cohort[col]=pd.to_numeric(scores[col],errors="raise")
    columns=RAW_COLUMNS+cfg["benchmark"]["covariates"]
    dev=read_table(project_path(root,f"tmp/cohort/{trait}.development.tsv.gz"))
    y=development_outcomes(dev,cohort)
    build=cohort.role.eq("train").to_numpy() & np.isfinite(y)
    tune=cohort.role.eq("tune").to_numpy() & np.isfinite(y)
    groups=cohort.group_id.to_numpy(str)
    tune_groups=sorted(set(groups[tune]),key=lambda g:hashlib.sha256((f"blocknet-{trait}|"+g).encode()).digest())
    tune_model=tune & np.isin(groups,tune_groups[::2])
    calibrate=tune & ~tune_model
    if min(build.sum(),tune_model.sum(),calibrate.sum())<2:
        raise ValueError("Direct models need train families and disjoint tuning/calibration families")
    design=NumericDesign(columns).fit(cohort.loc[build])
    side=design.transform(cohort)
    roles=np.full(len(cohort),"query",dtype=object)
    roles[build]="build"
    roles[tune_model]="tune_model"
    opt=dict(cfg["direct_blocknet"],_binary=binary, binary_objective="brier", initialize_prevalence=True)
    if binary and any(len(np.unique(y[mask])) != 2 for mask in (build, tune_model, calibrate)):
        raise ValueError("Binary fits require cases and controls in build, tuning and calibration families")
    bundle=dict(design=design,config=opt,feature_schema=regional_feature_schema(manifest,len(tokens)))
    device="cpu"
    best_epoch=None
    if method=="direct.linear":
        models=[Ridge(alpha=a).fit(side[build],y[build]) for a in (1.,10.,100.,1000.)]
        ridge=min(models,key=lambda m: np.mean((m.predict(side[tune_model])-y[tune_model])**2))
        raw=ridge.predict(side);bundle["ridge"]=ridge
    elif method=="direct.regional_ridge":
        from blocknet import feature_scaler
        from blocknet_ablation import regional_ridge
        mu,sd=feature_scaler(array,ix[build])
        raw,ridge=regional_ridge(array,ix,side,y,build,tune_model,mu,sd)
        bundle.update(ridge=ridge,block_mean=mu,block_sd=sd)
    else:
        import torch
        from blocknet import feature_scaler,train,predict_model
        device=torch.device(cfg["block_deep"]["device"])
        torch.set_num_threads(cfg["runtime"]["threads"])
        opt["amp"]="bf16" if device.type=="cuda" and torch.cuda.is_bf16_supported() else "off"
        kind="transformer" if method=="direct.blocknet" else "side_mlp"
        if kind=="side_mlp":
            array=np.empty((len(array),0,1),np.float32)
            mu=sd=np.empty((0,1),np.float32)
        else:mu,sd=feature_scaler(array,ix[build])
        torch.manual_seed(opt["seed"])
        arr_side=side.astype(np.float32)
        net,history,info=train(kind,array,ix,arr_side,y,roles,tokens,opt,device,mu,sd)
        raw,emb=predict_model(net,kind,array,ix,arr_side,np.arange(len(y)),mu,sd,opt,device,
                              info["ymean"],info["ysd"],embedding=kind=="transformer")
        if emb is not None:np.save(folder/"embedding.npy",emb,allow_pickle=False)
        atomic_table(history,folder/"training.tsv")
        best_epoch=info["best_epoch"]
        bundle.update(tokens=tokens,info=info,block_mean=mu,block_sd=sd,
                      weights={k:v.detach().cpu() for k,v in net.state_dict().items()})
    cal=(LogisticRegression(C=1e6, max_iter=1000) if binary else LinearRegression()).fit(raw[calibrate,None],y[calibrate])
    calibrated = cal.predict_proba(raw[:,None])[:,1] if binary else cal.predict(raw[:,None])
    bundle["calibration"]={method:cal}
    pred=pd.DataFrame({"id":cohort.id,method+".prediction":calibrated})
    joblib.dump(bundle,folder/"model.joblib",compress=3)
    atomic_table(pred,folder/"predictions.tsv.gz")
    atomic_json(dict(trait=trait, binary=binary,
                     methods=[method],
                     side_columns=columns,n_build=int(build.sum()),n_tune_model=int(tune_model.sum()),
                     n_prediction_calibration=int(calibrate.sum()),best_epoch=best_epoch,
                     device=str(device),threads=cfg["runtime"]["threads"],test_labels_used=False,calibration_labels_used=False),folder/"model.json")
