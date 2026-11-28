"""CPU controls on the identical consensus-block representation used by GRID.

Actual R glmnet and actual CPU LightGBM estimate corrections to the same frozen
CSx baseline. TRAIN uses family-OOF baseline residuals/offsets. TUNE selects
regularization, early stopping and a blend containing zero. No held-out Y enters.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import tempfile
import joblib
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from core2 import identified, aligned, atomic_table, atomic_json, check_roles
from abm2 import check_oof

PRIMARY = "csx.population"


class BlockTransform:
    """TRAIN-only centering and one scale, retaining relative block amplitude."""
    def fit(self, features, rows, chunk=512):
        self.width = int(features.shape[1])
        mean = np.zeros(self.width); m2 = mean.copy(); n = 0
        for start in range(0, len(rows), chunk):
            x = np.asarray(features[rows[start:start+chunk]], dtype=np.float64)
            if not np.isfinite(x).all():
                raise ValueError("Consensus scores must be finite; structural empty blocks are zero")
            count = len(x); batch_mean = x.mean(0); delta = batch_mean-mean
            m2 += np.square(x-batch_mean).sum(0)+delta**2*n*count/(n+count)
            mean += delta*count/(n+count); n += count
        if n < 20:
            raise ValueError("Block models need at least 20 observed TRAIN outcomes")
        variance = m2/n
        self.active = np.flatnonzero(variance > 1e-20)
        if not len(self.active):
            raise ValueError("All consensus block features are constant in TRAIN")
        self.mean = mean[self.active]
        self.scale = float(np.sqrt(variance[self.active].sum()))
        self.n_train = int(n)
        return self

    def transform(self, values):
        x = np.asarray(values, dtype=np.float32)
        if x.ndim != 2 or x.shape[1] != self.width or not np.isfinite(x).all():
            raise ValueError("Consensus feature width/nonfinite input differs from fitted model")
        return np.asarray((x[:, self.active]-self.mean)/self.scale, dtype=np.float32)

    def rows(self, features, rows, chunk=512):
        out = np.empty((len(rows), len(self.active)), dtype=np.float32)
        for start in range(0, len(rows), chunk):
            ix = rows[start:start+chunk]
            out[start:start+len(ix)] = self.transform(features[ix])
        return out


def combine(base, correction, binary, blend, link=False):
    base = np.asarray(base, float)
    if float(blend) == 0:
        return base.copy()
    if binary and link:
        return expit(logit(np.clip(base, 1e-7, 1-1e-7))+float(blend)*correction)
    raw = base+float(blend)*correction
    return np.clip(raw, 1e-7, 1-1e-7) if binary else raw


def inputs(cohort, outcomes, predictions, trait):
    cohort = identified(cohort); check_roles(cohort)
    p = aligned(predictions, cohort.id, "block model baseline")
    truth = aligned(outcomes, cohort.id, "development labels", allow_missing=True)
    if truth.loc[cohort.role.eq("test"),"y"].notna().any():
        raise ValueError("Development labels contain forbidden TEST outcomes")
    allowed = cohort.role.isin(["train", "tune"]).to_numpy()
    labels = np.full(len(cohort),np.nan)
    labels[allowed] = pd.to_numeric(truth.loc[allowed,"y"], errors="raise").to_numpy(float)
    if np.isinf(labels).any():
        raise ValueError("Nonfinite development label")
    train = np.flatnonzero(cohort.role.eq("train").to_numpy() & np.isfinite(labels))
    tune = np.flatnonzero(cohort.role.eq("tune").to_numpy() & np.isfinite(labels))
    if len(train) < 20 or len(tune) < 10:
        raise ValueError("Block models need observed TRAIN and separate TUNE labels")
    binary = trait == "t2dm"
    if binary and any(not np.isin(labels[r], [0, 1]).all() or len(np.unique(labels[r])) != 2 for r in (train,tune)):
        raise ValueError("Binary block models need cases and controls in TRAIN and TUNE")
    base = pd.to_numeric(p[PRIMARY+".prediction"], errors="raise").to_numpy(float)
    observed_train = np.zeros(len(cohort), dtype=bool)
    observed_train[train] = True
    oof = check_oof(cohort, p, observed_train)
    if not np.isfinite(base).all() or not np.isfinite(oof[train]).all():
        raise ValueError("Real family-held-out CSx predictions are mandatory")
    if binary and ((base<0).any() or (base>1).any() or (oof[train]<0).any() or (oof[train]>1).any()):
        raise ValueError("Binary baseline is not a probability")
    return cohort, labels, base, oof, train, tune, binary


def fit_glmnet(features, transform, y, base, oof, train, tune, binary, output, config, seed, rscript, threads=1):
    """R sees TRAIN/TUNE feature rows only, with no CAL/TEST phenotype file."""
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    rows = np.r_[train,tune]
    with tempfile.TemporaryDirectory(prefix="input-", dir=output) as tmp:
        folder = Path(tmp); width = len(transform.active)
        matrix = np.memmap(folder/"X.bin", mode="w+", dtype="<f8", shape=(len(rows), width))
        for start in range(0,len(rows),512):
            matrix[start:start+512] = transform.transform(features[rows[start:start+512]])
        matrix.flush(); del matrix
        offset = np.r_[oof[train],base[tune]]
        data = pd.DataFrame({"role":["train"]*len(train)+["tune"]*len(tune),
                             "y":y[rows], "base":offset})
        atomic_table(data,folder/"labels.tsv")
        specification = {"n":len(rows),"p":width,"binary":binary,"seed":int(seed),
            "alpha_grid":config.get("alpha_grid",[0,.5,1]),"nlambda":int(config.get("nlambda",40)),
            "maxit":int(config.get("maxit",100000)),"blend_grid":[0,.25,.5,1]}
        from parallel_runtime import worker_count
        matrix_gib=len(rows)*width*8/2**30
        specification["workers"]=worker_count(threads,len(specification["alpha_grid"]),
            memory_per_worker_gib=.5+4*matrix_gib,reserve_gib=2*matrix_gib)
        atomic_json(specification,folder/"config.json")
        command=[str(rscript),str(Path(__file__).with_name("block_glmnet.R")),str(folder),str(output)]
        with (output/"command.log").open("w") as log:
            subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,
                env={**os.environ,**{name:"1" for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS")}})
    coef = pd.read_csv(output/"coefficients.tsv",sep="\t")
    if (not {"feature","coefficient"} <= set(coef) or len(coef)!=width+1
            or not np.array_equal(pd.to_numeric(coef.feature,errors="raise"),np.arange(width+1))
            or not np.isfinite(coef.coefficient).all()):
        raise ValueError("Invalid glmnet coefficient export")
    metadata=json.loads((output/"fit.json").read_text())
    if metadata.get("test_outcomes_used") is not False or metadata.get("train_n") != len(train):
        raise ValueError("glmnet fit did not retain the common label boundary")
    return {"kind":"glmnet","coef":coef.coefficient.to_numpy(float),"blend":float(metadata["blend"]),
            "link":bool(binary),"fit":metadata}


def fit_lightgbm(xtrain, xtune, y, base, oof, train, tune, binary, config, seed, threads):
    import lightgbm as lgb
    target=y[train]-oof[train]; validation=y[tune]-base[tune]
    best_loss=float(np.square(y[tune]-base[tune]).mean())
    selected={"kind":"lightgbm","model":None,"blend":0.,"link":False}
    records=[]
    for leaves in config.get("lightgbm_leaves",[7,15]):
        params=dict(objective="regression",metric="l2",num_leaves=int(leaves),
            learning_rate=.03,min_data_in_leaf=max(20,min(200,len(train)//50)),
            lambda_l2=10.,feature_fraction=.8,bagging_fraction=.8,bagging_freq=1,
            max_bin=63,num_threads=int(threads),seed=int(seed),verbosity=-1,
            deterministic=True,force_col_wise=True,device_type="cpu")
        tr=lgb.Dataset(xtrain,label=target,free_raw_data=False)
        va=lgb.Dataset(xtune,label=validation,reference=tr,free_raw_data=False)
        model=lgb.train(params,tr,num_boost_round=int(config.get("lightgbm_rounds",400)),
            valid_sets=[va],callbacks=[lgb.early_stopping(int(config.get("lightgbm_patience",30)),verbose=False)])
        correction=model.predict(xtune,num_iteration=model.best_iteration)
        for blend in (.25,.5,1.):
            p=combine(base[tune],correction,binary,blend)
            value=float(np.square(y[tune]-p).mean())
            records.append(dict(num_leaves=int(leaves),blend=blend,iterations=int(model.best_iteration),tune_loss=value))
            if value < best_loss:
                best_loss=value;selected.update(model=model,blend=blend)
    selected["fit"]={"implementation":"actual CPU LightGBM","version":lgb.__version__,
        "training_target":"Y minus family-OOF CSx baseline; squared-error residual learner",
        "selection_metric":"Brier" if binary else "MSE","tune_loss":best_loss,"candidates":records,
        "train_n":len(train),"tune_n":len(tune),"test_outcomes_used":False}
    return selected


def corrections(model, values):
    if model["kind"]=="glmnet":
        return model["coef"][0]+values@model["coef"][1:]
    if model["model"] is None:
        return np.zeros(len(values))
    return model["model"].predict(values,num_iteration=model["model"].best_iteration)


def feature_schema(manifest):
    manifest=pd.DataFrame(manifest)
    columns=[c for c in ("block_id","name","map","population","chrom38","start38","end38","chrom37","start37","end37") if c in manifest]
    if not columns:
        raise ValueError("Block feature catalog needs identities or GRCh38 coordinates")
    result=manifest[columns].fillna("").astype(str).to_dict("records")
    if len({json.dumps(x,sort_keys=True) for x in result})!=len(result):
        raise ValueError("Duplicate block feature identities")
    return result


def predict_frozen(bundle, features, base, ids, batch_size=512, feature_manifest=None):
    if isinstance(bundle,(str,Path)):bundle=joblib.load(bundle)
    if len(features)!=len(ids) or len(base)!=len(ids):
        raise ValueError("Block model row count/ID alignment mismatch")
    if feature_manifest is None or feature_schema(feature_manifest)!=bundle["feature_schema"]:
        raise ValueError("Block query features need the same ordered feature catalog as training")
    if batch_size<1 or len(set(map(str,ids)))!=len(ids):
        raise ValueError("Prediction batch size must be positive and IDs unique")
    base=np.asarray(base,float)
    if not np.isfinite(base).all() or (bundle["binary"] and ((base<0).any() or (base>1).any())):
        raise ValueError("Invalid frozen baseline query predictions")
    out=pd.DataFrame({"id":np.asarray(ids,str)})
    for name,model in bundle["models"].items():
        pred=np.empty(len(ids))
        for start in range(0,len(ids),batch_size):
            stop=min(len(ids),start+batch_size)
            x=bundle["transform"].transform(features[np.arange(start,stop)])
            pred[start:stop]=combine(np.asarray(base)[start:stop],corrections(model,x),
                bundle["binary"],model["blend"],model["link"])
        out[name+".prediction"]=pred
    return out


def fit_block_models(cohort,outcomes,predictions,features,feature_manifest,output_dir,*,trait,
                     config=None,seed=20261115,threads=1,rscript="Rscript"):
    config=config or {}
    cohort,y,base,oof,train,tune,binary=inputs(cohort,outcomes,predictions,trait)
    if len(features)!=len(cohort) or features.shape[1]!=len(feature_manifest):
        raise ValueError("Consensus feature roster/catalog mismatch")
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    transform=BlockTransform().fit(features,train)
    models={}
    if config.get("glmnet",True):
        models["block.glmnet"]=fit_glmnet(features,transform,y,base,oof,train,tune,binary,out/"glmnet",config,seed,rscript,threads=threads)
    if config.get("lightgbm",True):
        xtrain=transform.rows(features,train);xtune=transform.rows(features,tune)
        models["block.lightgbm"]=fit_lightgbm(xtrain,xtune,y,base,oof,train,tune,binary,config,seed,threads)
        del xtrain,xtune
    if not models:raise ValueError("No block control model enabled")
    bundle={"transform":transform,"models":models,"binary":binary,
            "feature_catalog":feature_manifest.to_dict("records"),
            "feature_schema":feature_schema(feature_manifest),"baseline":PRIMARY}
    joblib.dump(bundle,out/"model.joblib")
    prediction=predict_frozen(bundle,features,base,cohort.id,feature_manifest=feature_manifest)
    atomic_table(prediction,out/"predictions.tsv.gz")
    metadata={"trait":trait,"n_catalog":features.shape[1],"n_active":len(transform.active),
        "train_n":len(train),"tune_n":len(tune),"common_baseline":PRIMARY,
        "normalization":"TRAIN centered, single scale preserving block amplitude",
        "methods":{k:{"blend":v["blend"],**v["fit"]} for k,v in models.items()},
        "calibration_outcomes_used":False,"test_outcomes_used":False}
    atomic_json(metadata,out/"model.json")
    return prediction,metadata
