"""Fair target-population calibration, genuine family-held-out residuals, and Disco stacking."""
from __future__ import annotations
from output_layout import project_path, trait_path, module_root
import argparse
import json
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.metrics import log_loss, mean_squared_error
from core2 import (POPS, VERSION, NumericDesign, read_table, identified, aligned,
                   atomic_table, atomic_json, check_roles, group_folds)

from parallel_runtime import ordered_map, worker_count

class CalibratedModel:
    def __init__(self, columns, task, penalty=1.):
        self.columns = list(columns); self.task = task; self.penalty = float(penalty)
    def fit(self, frame, y):
        y = np.asarray(y, float); good = np.isfinite(y)
        if good.sum() < max(10, len(self.columns)+2): raise ValueError("Insufficient labels for target calibration")
        if self.task == "binary" and (not np.isin(y[good], [0,1]).all() or len(np.unique(y[good])) != 2):
            raise ValueError("Binary calibration needs cases and controls")
        self.transformer = NumericDesign(self.columns)
        x = self.transformer.fit_transform(frame.loc[good])
        if self.task == "binary":
            self.model = LogisticRegression(C=1/max(self.penalty,1e-12), solver="lbfgs", max_iter=2000, tol=1e-8)
        else: self.model = Ridge(alpha=self.penalty)
        self.model.fit(x, y[good]); self.n_fit = int(good.sum())
        return self
    def predict(self, frame):
        x = self.transformer.transform(frame)
        return self.model.predict_proba(x)[:,1] if self.task == "binary" else self.model.predict(x)
    def coefficients(self):
        coef = np.asarray(self.model.coef_).reshape(-1)
        return pd.DataFrame({"feature":self.columns,"coefficient_standardized":coef,
            "coefficient_original":coef/self.transformer.scale,"train_mean":self.transformer.mean,
            "train_sd":self.transformer.scale,"train_median":self.transformer.median})

def loss(y, prediction, task):
    y, prediction = np.asarray(y,float),np.asarray(prediction,float)
    good = np.isfinite(y)&np.isfinite(prediction)
    if good.sum() < 10: return np.inf
    if task == "binary":
        if len(np.unique(y[good])) != 2: return np.inf
        return log_loss(y[good],np.clip(prediction[good],1e-8,1-1e-8),labels=[0,1])
    return mean_squared_error(y[good],prediction[good])

def choose_model(frame, y, fit, tune, candidates, covariates, task, penalties):
    records=[]; best=None
    for label, scores in candidates:
        for penalty in penalties:
            model=CalibratedModel([*covariates,*scores],task,penalty).fit(frame.loc[fit],y[fit])
            value=loss(y[tune],model.predict(frame.loc[tune]),task)
            records.append({"candidate":label,"penalty":penalty,"tune_loss":value,
                            "n_fit":model.n_fit,"n_tune":int(np.isfinite(y[tune]).sum())})
            if best is None or value < best[0]: best=(value,model,label)
    if best is None or not np.isfinite(best[0]): raise ValueError("No valid model on independent tune set")
    return best[1],best[2],records

def method_candidates(manifest):
    d=manifest.copy(); d["phi"]=d.phi.astype(str)
    d["family"]=d.family.replace({"prscs":"cs"})
    methods={"covariate":[("covariates",[])]}
    joint=[]
    for phi,rows in d[d.family.eq("csx")].groupby("phi",sort=False):
        mapping=dict(zip(rows.population,rows.score_name))
        if all(p in mapping for p in POPS): joint.append((phi,[mapping[p] for p in POPS]))
    if not joint: raise ValueError("No complete four-ancestry CSx score set")
    methods["csx.population"]=joint
    for method,phi_auto in (("csx.auto",True),("csx.meta",False)):
        rows=d[d.family.eq("csx_meta") & (d.phi.eq("auto") if phi_auto else ~d.phi.eq("auto"))]
        if not rows.empty: methods[method]=[(r.phi,[r.score_name]) for r in rows.itertuples()]
    singles=[]
    for pop in POPS:
        rows=d[d.family.eq("cs")&d.population.eq(pop)]
        if rows.empty: raise ValueError(f"Missing independent PRS-CS {pop} baseline")
        choices=[(f"{pop}:{r.phi}",[r.score_name]) for r in rows.itertuples()]
        methods[f"cs.{pop}"]=choices;singles+=choices
    methods["cs.best"]=singles
    independent=[]
    for phi,rows in d[d.family.eq("cs")].groupby("phi",sort=False):
        mapping=dict(zip(rows.population,rows.score_name))
        if all(p in mapping for p in POPS): independent.append((phi,[mapping[p] for p in POPS]))
    if independent: methods["cs.linear_combination"]=independent
    return methods

def inputs(cfg, root, trait):
    root=Path(root)
    cohort=identified(read_table(project_path(root, "tmp/cohort/participants.tsv.gz")));check_roles(cohort)
    scores=aligned(read_table(project_path(root, f"tmp/{trait}/scores.tsv.gz")),cohort.id,"raw scores")
    manifest=pd.read_csv(project_path(root, f"tmp/{trait}/inference/weights/manifest.tsv"),sep="\t",dtype=str)
    names=manifest.score_name.tolist()
    if set(names)-set(scores):raise ValueError("Scoring output missing weight-manifest score names")
    for col in names:
        scores[col]=pd.to_numeric(scores[col],errors="coerce")
        if not np.isfinite(scores[col]).all(): raise ValueError(f"Missing/nonfinite raw PRS: {col}")
    frame=cohort.merge(scores[["id",*names]],on="id",validate="one_to_one",sort=False)
    for col in cfg["benchmark"]["covariates"]:
        if col not in frame: raise ValueError(f"Missing common covariate {col}")
    # The fitting entry has no reason to read the separately stored test outcomes.
    outcomes=aligned(read_table(project_path(root, f"tmp/cohort/{trait}.development.tsv.gz")),frame.id,"development outcomes",allow_missing=True)
    y=pd.to_numeric(outcomes.y,errors="coerce").to_numpy(float)
    if np.isfinite(y[frame.role.eq("test")]).any():raise ValueError("Development outcome file contains test labels")
    task="binary" if trait=="t2dm" else "continuous"
    return frame,y,manifest,task

def fixed_input_set(candidates, requested=None):
    """A declared discovery phi, independent of target labels and grid ordering."""
    if requested is None:
        fixed=[choice for choice in candidates if choice[0] != "auto"]
        return (fixed or candidates)[0]
    for label,columns in candidates:
        if str(requested)=="auto" and label=="auto": return label,columns
        if str(requested)!="auto" and label!="auto" and float(label)==float(requested):
            return label,columns
    raise ValueError(f"benchmark.input_phi={requested} is absent from inference.phi_grid; include it explicitly")

def fit_models(cfg, root, trait):
    root=Path(root)
    frame,y,manifest,task=inputs(cfg,root,trait)
    fit=frame.role.eq("train").to_numpy();tune=frame.role.eq("tune").to_numpy()
    opts=cfg["benchmark"];cov=opts["covariates"];penalties=opts.get("penalties",[.01,1.,100.])
    methods=method_candidates(manifest)
    geometry_label,geometry_cols=fixed_input_set(methods["csx.population"],opts.get("input_phi"))
    methods["csx.input"]=[(geometry_label,geometry_cols)]
    prediction=frame[["id"]].copy();models={};tuning=[];coefficients=[];selected={}
    budget=int(cfg.get("runtime",{}).get("threads",1))
    frame_gib=frame.memory_usage(deep=True).sum()/2**30
    count=worker_count(budget,len(methods),memory_per_worker_gib=.5+4*frame_gib,reserve_gib=frame_gib)
    def fit_method(name,candidates):
        model,label,records=choose_model(frame,y,fit,tune,candidates,cov,task,penalties)
        return name,model,label,records,model.predict(frame),model.coefficients().assign(method=name,selected=label)
    print(f"[fit] independent model processes={count}; native threads/process={max(1,budget//count)}",flush=True)
    for name,model,label,records,values,coef in ordered_map(fit_method,methods.items(),count,inner_threads=max(1,budget//count)):
        prediction[f"{name}.prediction"]=values
        models[name]=model;selected[name]=label
        tuning += [{"method":name,**r} for r in records]
        coefficients.append(coef)
    # The geometry uses a fixed discovery-only weight set. It is NOT chosen by Y.
    for pop,col in zip(POPS,geometry_cols):prediction[f"csx.{pop}.raw_score"]=frame[col].to_numpy(float)
    nfold=int(opts.get("oof_folds",5));folds=np.full(len(frame),-1,dtype=int)
    folds[fit]=group_folds(frame.loc[fit,"group_id"],nfold,int(cfg["cohort"]["seed"]))
    prediction["oof_fold"]=folds;fold_records=[]
    # Retune within each OOF fit; tuning the full fit once would let a held-out
    # donor's own label influence the chosen model/phi indirectly.
    def fit_fold(name,k):
        train=fit&(folds!=k);held=fit&(folds==k)
        if set(frame.loc[train,"group_id"])&set(frame.loc[held,"group_id"]):raise ValueError("OOF family leakage")
        model,label,records=choose_model(frame,y,train,tune,methods[name],cov,task,penalties)
        record={"method":name,"fold":k,"selected":label,"penalty":model.penalty,
                "n_train":int(train.sum()),"n_held":int(held.sum()),"family_overlap":0}
        return name,np.flatnonzero(held),model.predict(frame.loc[held]),record
    tasks=[(name,k) for name in ("csx.population","covariate") for k in range(nfold)]
    count=worker_count(budget,len(tasks),memory_per_worker_gib=.5+4*frame_gib,reserve_gib=frame_gib)
    oofs={name:np.full(len(frame),np.nan) for name in ("csx.population","covariate")}
    for name,rows,values,record in ordered_map(fit_fold,tasks,count,inner_threads=max(1,budget//count)):
        oofs[name][rows]=values;fold_records.append(record)
    for name,values in oofs.items():prediction[f"{name}.oof_prediction"]=values
    if not np.isfinite(prediction.loc[fit,"csx.population.oof_prediction"]).all():raise ValueError("Incomplete OOF baseline")
    atomic_table(prediction,trait_path(root, trait, "baseline.predictions.tsv.gz"))
    atomic_table(frame[["id","group_id","role"]].assign(oof_fold=folds),trait_path(root, trait, "oof_folds.tsv.gz"))
    atomic_table(pd.DataFrame(tuning),project_path(root, f"qc/{trait}.benchmark_tuning.tsv"))
    atomic_table(pd.concat(coefficients,ignore_index=True),project_path(root, f"qc/{trait}.coefficients.tsv"))
    atomic_json({"schema":VERSION,"task":task,"primary":"csx.population","selected":selected,
        "geometry_phi":geometry_label,"geometry_columns":geometry_cols,"covariates":cov,
        "normalization":"TRAIN-only medians, means and SDs; each OOF fold refits transforms",
        "folds":fold_records,"tuning_loss":"log_loss" if task=="binary" else "MSE",
        "training_budget_shared_with_abm":True,"test_outcomes_read":False},trait_path(root, trait, "benchmark.json"))
    joblib.dump(models,trait_path(root, trait, "benchmark_models.joblib"))
    return prediction

def fit_disco(cfg, root, trait):
    from disco2 import DiscoDivas,RawScoreInterpolation,DiscoFitError
    root=Path(root);out=trait_path(root, trait, "disco");out.mkdir(parents=True,exist_ok=True)
    frame,y,manifest,task=inputs(cfg,root,trait)
    predictions=aligned(read_table(trait_path(root, trait, "baseline.predictions.tsv.gz")),frame.id,"benchmark predictions")
    meta=json.loads((trait_path(root, trait, "benchmark.json")).read_text());cov=cfg["benchmark"]["covariates"]
    score_cols=meta["geometry_columns"] # fixed predeclared input-PRS phi, as for ABM geometry
    opts=cfg["disco"];fit=frame.role.eq("train").to_numpy();tune=frame.role.eq("tune").to_numpy()
    folds=pd.to_numeric(predictions.oof_fold).to_numpy(int);records=[];model_records=[]
    args=dict(score_cols=score_cols,task=task,pc_cols=[f"PC{k}" for k in range(1,int(opts.get("distance_pcs",10))+1)],
              covariate_cols=cov,min_anchor_n=int(opts.get("min_anchor_n",100)),min_cases=int(opts.get("min_cases",20)),
              penalty=float(opts.get("penalty",1e-6)),source_method="PRS-CSx")
    oof=np.full(len(frame),np.nan);fitted=None;status="complete";reason=""
    try:
        def fit_fold(k):
            tr=fit&(folds!=k)&np.isfinite(y);held=fit&(folds==k)
            model=DiscoDivas(**args)
            try:
                model.fit(frame.loc[tr],y[tr],covariates=NumericDesign(cov).fit_transform(frame.loc[tr]))
                return k,np.flatnonzero(held),model.predict(frame.loc[held]),model.anchor_table().assign(fold=str(k)),model.summary(),None
            except DiscoFitError as exc:
                return k,None,None,model.anchor_table().assign(fold=str(k)),None,str(exc)
        nfold=int(cfg["benchmark"].get("oof_folds",5));budget=int(cfg.get("runtime",{}).get("threads",1))
        frame_gib=frame.memory_usage(deep=True).sum()/2**30
        count=worker_count(budget,nfold,memory_per_worker_gib=.5+4*frame_gib,reserve_gib=frame_gib)
        print(f"[disco] OOF processes={count}; native threads/process={max(1,budget//count)}",flush=True)
        pieces=list(ordered_map(fit_fold,((k,) for k in range(nfold)),count,inner_threads=max(1,budget//count)))
        for k,held,values,anchors,summary,error in pieces:
            records.append(anchors)
            if error is not None:raise DiscoFitError(error)
            oof[held]=values;model_records.append({"fold":k,**summary})
        if not np.isfinite(oof[fit]).all():raise DiscoFitError("Nonfinite cross-fitted genetic scores")
        tr=fit&np.isfinite(y);fitted=DiscoDivas(**args)
        try:fitted.fit(frame.loc[tr],y[tr],covariates=NumericDesign(cov).fit_transform(frame.loc[tr]))
        finally:records.append(fitted.anchor_table().assign(fold="full"))
        genetic=fitted.predict(frame)
        if not np.isfinite(genetic).all():raise DiscoFitError("Unstable interpolation geometry for some participants; see details")
        cal_frame=frame.copy();cal_frame["disco_genetic"]=genetic
        cal_frame.loc[fit,"disco_genetic"]=oof[fit]
        model,_,selection=choose_model(cal_frame,y,fit,tune,[("crossfit_genetic",["disco_genetic"])],cov,task,cfg["benchmark"].get("penalties",[.01,1.,100.]))
        predictions["disco.prediction"]=model.predict(cal_frame)
        predictions["disco.genetic_score"]=cal_frame.disco_genetic.to_numpy()
        atomic_table(fitted.predict_details(frame).assign(id=frame.id.to_numpy()),out/"interpolation.tsv.gz")
        atomic_table(fitted.coefficient_table(),project_path(root, f"qc/{trait}.disco_coefficients.tsv"))
        atomic_table(pd.DataFrame(selection),project_path(root, f"qc/{trait}.disco_tuning.tsv"))
        joblib.dump({"disco":fitted,"calibration":model},out/"model.joblib")
    except DiscoFitError as exc:
        status="unavailable";reason=str(exc)
        # Omission is explicit. Never substitute a global model under this name.
        predictions=predictions.drop(columns=[c for c in ("disco.prediction","disco.genetic_score") if c in predictions])
        if fitted is not None and getattr(fitted,"is_fitted_",False):
            atomic_table(fitted.predict_details(frame).assign(id=frame.id.to_numpy()),out/"interpolation.tsv.gz")
    if records:atomic_table(pd.concat(records,ignore_index=True),project_path(root, f"qc/{trait}.disco_anchors.tsv"))
    if opts.get("raw_interpolation_ablation",True):
        raw=RawScoreInterpolation(dict(zip(POPS,score_cols)),pc_cols=args["pc_cols"],min_anchor_n=int(opts.get("min_anchor_n",100)))
        try:
            raw.fit(frame.loc[fit]);raw_frame=frame.copy();raw_frame["raw_disco"]=raw.predict(frame)
            if not np.isfinite(raw_frame.raw_disco).all():raise DiscoFitError("Nonfinite raw interpolation")
            model,_,_=choose_model(raw_frame,y,fit,tune,[("raw_input",["raw_disco"])],cov,task,cfg["benchmark"].get("penalties",[.01,1.,100.]))
            predictions["disco.raw_interpolation.prediction"]=model.predict(raw_frame)
        except DiscoFitError as exc:model_records.append({"raw_ablation_unavailable":str(exc)})
    atomic_table(predictions,trait_path(root, trait, "predictions.tsv.gz"))
    atomic_json({"schema":VERSION,"status":status,"reason":reason,"input_phi":meta["geometry_phi"],
        "source":"PRS-CSx-input DiscoDivas; full fine-tuning input layer, not a replication of all original paper datasets",
        "downstream_calibration":"TRAIN family-cross-fitted genetic score + common covariates; TUNE selects penalty",
        "fold_summaries":model_records,"full":fitted.summary() if status=="complete" else None,"test_outcomes_read":False},out/"model.json")
    return predictions

if __name__=="__main__":
    from core2 import load_config
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("mode",choices=["fit","disco"]);p.add_argument("--config",required=True);p.add_argument("--root",required=True);p.add_argument("--trait",required=True)
    a=p.parse_args();(fit_models if a.mode=="fit" else fit_disco)(load_config(a.config),a.root,a.trait)
