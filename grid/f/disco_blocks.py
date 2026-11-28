"""Disco anchor interpolation inside regional prediction; all labels come from development."""
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.linear_model import Ridge, LogisticRegression

from abm2 import RAW_COLUMNS
from block_models import BlockTransform, combine, corrections, fit_glmnet, inputs
from core2 import NumericDesign, aligned, atomic_json, atomic_table, read_table
from driver2 import regional_inputs
from output_layout import project_path, trait_path

POPS=("EUR","AFR","EAS","SAS")


class InterpolatedRegions:
    """One token per common block; separate annotated tokens for distinct ancestry BEDs."""
    def __init__(self,src,manifest,weights,shared):
        self.src,self.weights,self.shared=src,weights,shared
        pop=manifest["map"].to_numpy(str)
        if set(pop)!=set(POPS):
            raise ValueError("Disco regional methods require all four CSx channels")
        self.cols=[np.flatnonzero(pop==a) for a in POPS]
        if shared:
            token=manifest.iloc[self.cols[0]].copy()
            for ix in self.cols[1:]:
                if not np.array_equal(token.block_id,manifest.iloc[ix].block_id):
                    raise ValueError("Shared BED channel coordinates differ")
            token["map"]="ALL"; channels=5
        else:
            token=manifest;channels=2
        self.pop_index=np.array([POPS.index(a) for a in pop])
        suffix="37" if "chrom37" in token else "38"
        start=token["start"+suffix].to_numpy(float);end=token["end"+suffix].to_numpy(float)
        self.tokens=pd.DataFrame(dict(CHR=token["chrom"+suffix].astype(int).to_numpy(),
                  source_population=token["map"].to_numpy(),midpoint=(start+end)/2,length=end-start))
        self.shape=(src.shape[0],len(token),channels)
    def __len__(self):return self.shape[0]
    def __getitem__(self,rows):
        raw=np.asarray(self.src[rows],np.float32)
        w=self.weights[rows]
        if self.shared:
            x=np.stack([raw[:,j] for j in self.cols],axis=-1)
            mix=np.einsum("nbp,np->nb",x,w)
            return np.concatenate((x,mix[:,:,None]),axis=2)
        return np.stack([raw,raw*w[:,self.pop_index]],axis=2)


class MixedFeatures:
    def __init__(self,tokens):
        self.tokens=tokens;self.shape=tokens.shape[:2]
    def __getitem__(self,rows):return self.tokens[rows][:,:,-1]


def smooth_offset(cohort,cols,y,train,fold,binary):
    """Cross-fit age curvature, sex-specific age effects and global PRS effects with family-held-out folds."""
    def matrix(design,frame):
        x=design.transform(frame);age=x[:,cols.index("age")];sex=x[:,cols.index("sex")]
        spline=np.stack([age,age**2,np.maximum(age,0)**2],axis=1)
        return np.column_stack([x,spline[:,1:],sex[:,None]*spline,x[:,:4]*sex[:,None],x[:,:4]*age[:,None]])
    def fit(rows):
        design=NumericDesign(cols).fit(cohort.iloc[rows])
        model=LogisticRegression(C=1,max_iter=500) if binary else Ridge(alpha=10.)
        model.fit(matrix(design,cohort.iloc[rows]),y[rows])
        return design,model
    def predict(design,model,rows):
        x=matrix(design,cohort.iloc[rows])
        return model.predict_proba(x)[:,1] if binary else model.predict(x)
    held=np.full(len(y),np.nan)
    for k in np.unique(fold[train]):
        sel=train[fold[train]==k]
        held[sel]=predict(*fit(train[fold[train]!=k]),sel)
    fitted=fit(train)
    return predict(*fitted,np.arange(len(y))),held,fitted


def blend_with_disco(raw,baseline,gates,y,cal,strength):
    """Estimate a bounded blend at each anchor on CAL, shrinking small strata toward the global blend."""
    from scipy.optimize import lsq_linear
    if not np.isfinite(strength) or strength<=0:raise ValueError("anchor prior size must be positive")
    # Signed interpolation is preserved inside the model. A convex deployment blend uses positive memberships.
    membership=np.maximum(gates.astype(float),0.)
    membership/=membership.sum(1,keepdims=True)
    delta=raw-baseline
    d=delta[cal];residual=y[cal]-baseline[cal];total=d@d
    global_blend=float(np.clip(d@residual/total,0,1)) if total>0 else 0.
    prior=np.sqrt(max(float(np.mean(d*d)),1e-12)*strength)
    x=np.vstack((d[:,None]*membership[cal],np.eye(gates.shape[1])*prior))
    z=np.r_[residual,np.full(gates.shape[1],prior*global_blend)]
    result=lsq_linear(x,z,bounds=(0,1))
    if not result.success:raise ValueError("Disco blend did not converge")
    return baseline+(membership@result.x)*delta,global_blend,result.x


def fit(cfg,root,trait,method,out):
    cohort=read_table(project_path(root,"tmp/cohort/participants.tsv.gz"))
    dev=read_table(project_path(root,f"tmp/cohort/{trait}.development.tsv.gz"))
    scores=aligned(read_table(trait_path(root,trait,"predictions.tsv.gz")),cohort.id)
    baseline=aligned(read_table(trait_path(root,trait,"baseline.predictions.tsv.gz")),cohort.id)
    cohort,y,_,_,train,tune,binary=inputs(cohort,dev,baseline,trait)
    if "disco.genetic_score" not in scores or "disco.prediction" not in scores:
        raise ValueError("Completed Disco anchor interpolation is required for "+method)
    details=aligned(read_table(root/trait/"disco.weights"/"weights.tsv.gz"),cohort.id)
    gates=details[["gate_"+p for p in POPS]].to_numpy(float)
    slopes=details[["slope_"+p for p in POPS]].to_numpy("float32")
    if not np.isfinite(gates).all() or not np.allclose(gates.sum(1),1):
        raise ValueError("Invalid Disco anchor weights")
    features,_=regional_inputs(root,trait,cohort)
    manifest=read_table(trait_path(root,trait,"regions/features.tsv"))
    tokens=InterpolatedRegions(features,manifest,slopes,cfg["regions"]["layout"]=="consensus")
    for c in RAW_COLUMNS:cohort[c]=baseline[c].to_numpy(float)
    cols=RAW_COLUMNS+cfg["benchmark"]["covariates"]
    seed=cfg["cohort"]["seed"];threads=cfg["runtime"]["threads"]
    store=dict(token_catalog=manifest,anchors=list(POPS),shared=tokens.shared)
    if method=="block.glmnet.disco":
        base,oof,offset=smooth_offset(cohort,cols,y,train,baseline.oof_fold.to_numpy(int),binary)
        mixed=MixedFeatures(tokens);transform=BlockTransform().fit(mixed,train)
        opt=cfg["block_models"]
        model=fit_glmnet(mixed,transform,y,base,oof,train,tune,binary,out/"glmnet",opt,seed,cfg["runtime"]["rscript"],threads)
        raw=base.copy()
        for start in range(0,len(y),512):
            r=slice(start,start+512)
            raw[r]=combine(base[r],corrections(model,transform.transform(mixed[r])),binary,model["blend"],binary)
        store.update(offset=offset,columns=cols,transform=transform,model=model)
    else:
        import torch
        from blocknet import feature_scaler,predict_model,train as fit_net
        # 3.disco exports family-out-of-fold genetic predictions for TRAIN, frozen predictions elsewhere.
        cohort["disco.genetic_score"]=scores["disco.genetic_score"].to_numpy(float)
        for i,a in enumerate(POPS):cohort["weight_"+a]=gates[:,i]
        cols=cols+["disco.genetic_score"]+["weight_"+a for a in POPS]
        design=NumericDesign(cols).fit(cohort.iloc[train]);side=design.transform(cohort).astype("float32")
        ix=np.arange(len(y));roles=np.full(len(y),"query",dtype=object);roles[train]="build";roles[tune]="tune_model"
        kind="side_mlp" if method=="blocknet.disco.global" else "transformer"
        arr=tokens if kind=="transformer" else np.empty((len(y),0,1),np.float32)
        mu,sd=feature_scaler(arr,train) if kind=="transformer" else (np.empty((0,1)),np.empty((0,1)))
        opt=dict(cfg["direct_blocknet"],_binary=binary,binary_objective="brier",initialize_prevalence=True)
        devc=torch.device(cfg["block_deep"]["device"]);torch.set_num_threads(threads);torch.manual_seed(opt["seed"])
        opt["amp"]="bf16" if devc.type=="cuda" and torch.cuda.is_bf16_supported() else "off"
        mdl,history,info=fit_net(kind,arr,ix,side,y,roles,tokens.tokens,opt,devc,mu,sd)
        raw,_=predict_model(mdl,kind,arr,ix,side,ix,mu,sd,opt,devc,info["ymean"],info["ysd"])
        if binary:raw=expit(raw)
        atomic_table(history,out/"training.tsv")
        store.update(design=design,columns=cols,config=opt,state={k:v.detach().cpu() for k,v in mdl.state_dict().items()},
                     info=info,tokens=tokens.tokens,mu=mu,sd=sd)
    # Calibration alone selects a convex blend including unchanged Disco. No TEST labels enter.
    yy=pd.to_numeric(aligned(dev,cohort.id,allow_missing=True).y).to_numpy(float)
    cal=cohort.role.eq("calibration").to_numpy()&np.isfinite(yy)
    if not cal.any():raise ValueError("Missing CALIBRATION labels for the Disco blend")
    disco=scores["disco.prediction"].to_numpy(float)
    prior=float(cfg.get("disco_blocks",{}).get("anchor_prior_n",200))
    if not np.isfinite(raw).all() or not np.isfinite(disco).all():raise ValueError("Nonfinite hybrid predictions")
    final,blend,anchors=blend_with_disco(raw,disco,gates,yy,cal,prior)
    cols_out={}
    if method=="block.glmnet.disco":
        matched,_,_=blend_with_disco(base,disco,gates,yy,cal,prior)
        cols_out[method+".offset.prediction"]=matched
    atomic_table(pd.DataFrame({"id":cohort.id,**cols_out,method+".prediction":final,
                   method+".global_blend.prediction":disco+blend*(raw-disco),method+".unblended":raw}),out/"predictions.tsv.gz")
    store.update(blend=blend,anchor_blend=anchors);joblib.dump(store,out/"model.joblib",compress=2)
    atomic_json(dict(method=method,trait=trait,train=len(train),tune=len(tune),calibration=int(cal.sum()),
                     global_blend=blend,anchor_blend=dict(zip(POPS,anchors.tolist())),anchor_prior_n=prior,regional_tokens=0 if method=="blocknet.disco.global" else tokens.shape[1],channels=tokens.shape[2],test_labels_used=False),out/"model.json")
