"""Day-365 common-test comparison and a blinded baseline molecular profile of early cases."""
import argparse
from pathlib import Path
import json
import os
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from threadpoolctl import threadpool_limits

from abm import CensoringKM
from c1_protocol import _read_group_file, ensure_roster, roster_record
from dx_layout import scope_for_read
from file_lock import file_lock

HERE = Path(__file__).resolve().parent
HORIZON = 10.
SEED = 2026


def read_table(path):
    return pd.read_csv(path,sep="\t",dtype={"eid":str,"group":str})


def frozen_scale(x, train):
    """Estimate means/SDs on eligible training rows, in bounded column blocks."""
    keep, mu, sd = [], [], []
    for i in range(0,x.shape[1],32):
        a=np.asarray(x[:,i:i+32][train],float)
        a[~np.isfinite(a)] = np.nan
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            m,s=np.nanmean(a,axis=0),np.nanstd(a,axis=0)
        ok=(np.isfinite(a).sum(axis=0)>=3) & np.isfinite(s) & (s>1e-8)
        keep.extend(i+np.flatnonzero(ok)); mu.extend(m[ok]); sd.extend(s[ok])
    mat=np.empty((len(x),len(keep)),dtype=np.float32,order="C")
    for j, (k,m,s) in enumerate(zip(keep,mu,sd)):
        col=(x[:,k]-m)/s
        mat[:,j]=np.where(np.isfinite(col),np.clip(col,-8,8),0)
    return mat,np.asarray(keep,dtype=int)


def fitted(x, rows, y, weights, target):
    # Match mean-loss penalty across training populations of different sizes.
    weight=weights[rows]
    model=LogisticRegression(C=100./weight.sum(),solver="lbfgs",max_iter=600,tol=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter("error",ConvergenceWarning)
        model.fit(x[rows],y[rows],sample_weight=weight)
    return model.predict_proba(x[target])[:,1]


def auc_layout(scores):
    order=np.argsort(scores,kind="stable")
    first=np.r_[0,1+np.flatnonzero(np.diff(scores[order]))]
    return order,first


def metrics(y,w,p,order,starts):
    # Count concordant pairs, assigning half credit to tied scores.
    cases=np.add.reduceat((w*y)[order],starts)
    controls=np.add.reduceat((w*(1-y))[order],starts)
    mass=cases.sum()*controls.sum()
    auc=np.sum(cases*(controls.cumsum()-.5*controls))/mass if mass>0 else np.nan
    # The upper cap rounds to 1 in float32. Promote before clipping for the log terms.
    log_p=np.clip(np.asarray(p,dtype=np.float64),1.e-8,1.-1.e-8)
    loss=np.where(y==1,-np.log(log_p),-np.log1p(-log_p))
    p=np.clip(p,1.e-8,1-1.e-8)
    return np.array([auc,np.average((p-y)**2,weights=w),
                     np.average(loss,weights=w),
                     np.average(p-y,weights=w)])


def compare(y,w,g,broad,clean,reps,label,reference,model,seed=SEED):
    orders=[auc_layout(p) for p in (broad,clean)]
    _,idx=np.unique(g,return_inverse=True)
    blocks=idx.max()+1
    def score(weights):
        return np.array([metrics(y,weights,p,*order) for p,order in zip((broad,clean),orders)])
    before,after=score(w)
    rng=np.random.default_rng(seed)
    diffs=[]
    for _ in range(reps):
        freq=np.bincount(rng.integers(blocks,size=blocks),minlength=blocks)
        a,b=score(w*freq[idx])
        diffs.append(b-a)
    rows=[]
    for k,metric in enumerate(("AUC_IPCW","Brier_IPCW","LogLoss_IPCW","predicted_minus_observed_risk")):
        d=np.asarray(diffs)[:,k]
        keep=np.isfinite(d)
        lo,hi=np.quantile(d[keep],[.025,.975]) if keep.sum()>=max(20,reps*.8) else (np.nan,np.nan)
        rows.append(dict(analysis=label,reference=reference,model=model,metric=metric,reference_value=before[k],
                    model_value=after[k],difference=after[k]-before[k],lower=lo,upper=hi,
                    test_known_N=len(y),test_cases=int(y.sum()),test_families=blocks,valid_bootstraps=int(keep.sum())))
    return rows


def assign_profile_groups(d):
    early=d.any1_state.eq("excluded_year1")
    prev=d.any1_state.eq("prevalent") & d.icd_state.eq("prevalent") & d.ref_state.eq("prevalent")
    late=d.any1_state.eq("incident") & d.icd_state.eq("incident") & d.ref_state.eq("incident") & d.any_first_days.le(365.25*HORIZON)
    return np.select([prev,late,early & d.any_prior,early & d.any_undated,
                early & d.any_early_days.le(90),early & d.any_early_days.le(180),early,
                d.any1_state.eq("uncertain"),d.any1_state.eq("control")],
               ["prevalent_anchor","later_incident_anchor","early_with_prior","early_with_undated",
                "early_day_1_90","early_day_91_180","early_day_181_365","uncertain_timing","control"],default="other")


def profiles(d,raw,names,train,test,private,groups=None):
    groups=assign_profile_groups(d) if groups is None else groups
    pos=groups=="prevalent_anchor"; neg=groups=="later_incident_anchor"
    fit=train & (pos|neg)
    counts=pd.crosstab(pd.Series(groups,name="profile_group"),pd.Series(np.where(test,"test","training"),name="split")).reset_index()
    out={"timing_groups":counts}
    query=test & (groups!="other")
    if min((fit&pos).sum(),(fit&neg).sum())<10 or not query.any():
        return out,"insufficient_training_anchors_or_test"
    x,sel=frozen_scale(raw,fit)
    if not x.shape[1]:
        return out,"no_estimable_biomarkers"
    weights=np.ones(len(d))
    for mask in (pos,neg):
        weights[fit&mask]=fit.sum()/(2*(fit&mask).sum())
    y=pos.astype(int)
    models={"omics":x}
    # Compare with simple clinical descriptors on the identical anchor sample.
    cv=[];used=[]
    for key in ("age","sex","drug.lipid"):
        if key=="drug.lipid" and os.getenv("LE8_TREATMENT_BASELINE_CONFIRMED","").lower() not in ("true","1","yes"): continue
        if key not in d: continue
        if key=="sex":
            s=d[key].astype("string").str.lower().map({"m":1,"male":1,"f":0,"female":0,"1":1,"0":0,"1.0":1,"0.0":0})
        else: s=pd.to_numeric(d[key],errors="coerce")
        cv.append(s.to_numpy(dtype=float,na_value=np.nan));used.append(key)
    if cv:
        z,available=frozen_scale(np.column_stack(cv),fit)
        if z.shape[1]:
            models["clinical_descriptors"]=z
            models["omics_plus_descriptors"]=np.column_stack((x,z))
    descriptors=",".join(used[i] for i in available) if cv else ""
    scores=[]; feature_rows=[]; validation=[]
    individual=d.loc[query,["eid"]].copy()
    individual["profile_group"]=groups[query]
    for label,matrix in models.items():
        p=fitted(matrix,fit,y,weights,query)
        individual[label]=p
        qg=groups[query]; anchors=(pos|neg)[query]
        if len(np.unique(y[query][anchors]))==2:
            true=y[query][anchors]; pred=p[anchors]
            validation.append(dict(model=label,clinical_variables=descriptors if label!="omics" else "",prevalent_N=int(true.sum()),later_incident_N=int((1-true).sum()),
                heldout_anchor_AUC=roc_auc_score(true,pred),prevalent_sensitivity_at_half=float((pred[true==1]>=.5).mean()),
                later_incident_specificity_at_half=float((pred[true==0]<.5).mean())))
        for g in np.unique(qg):
            ps=p[qg==g]
            scores.append(dict(model=label,profile_group=g,N=len(ps),q25=np.quantile(ps,.25),median=np.median(ps),
                               q75=np.quantile(ps,.75),fraction_above_half=float((ps>=.5).mean())))
    if "baseline_lipid_i0" in d:
        med = pd.to_numeric(d.loc[query,"baseline_lipid_i0"],errors="coerce").to_numpy()
        qg = groups[query]
        score = individual["omics"].to_numpy()
        rows = []
        for label,mask in (("nonuser",med==0),("user",med==1),("unknown",~np.isin(med,[0,1]))):
            anchors = mask & np.isin(qg,["prevalent_anchor","later_incident_anchor"])
            is_prevalent = qg[anchors]=="prevalent_anchor"
            auc = roc_auc_score(is_prevalent,score[anchors]) if len(np.unique(is_prevalent))==2 else np.nan
            for group in np.unique(qg[mask]):
                values = score[mask & (qg==group)]
                rows.append(dict(baseline_lipid_treatment=label,profile_group=group,N=len(values),
                    q25=np.quantile(values,.25),median=np.median(values),q75=np.quantile(values,.75),
                    heldout_anchor_AUC=auc,prevalent_anchor_N=int(is_prevalent.sum()),
                    later_incident_anchor_N=int((~is_prevalent).sum())))
        out["timing_medication"] = pd.DataFrame(rows)
    for j in sel:
        col=np.asarray(raw[:,j],float)
        ref=col[fit & np.isfinite(col)]; center=ref.mean(); scale=ref.std()
        for g in np.unique(groups[query]):
            values=col[query & (groups==g) & np.isfinite(col)]
            if len(values):
                feature_rows.append(dict(feature=names[j],profile_group=g,N=len(values),mean=values.mean(),median=np.median(values),
                    mean_z=(values.mean()-center)/scale))
    individual.to_csv(private,sep="\t",index=False,compression="gzip")
    out.update(timing_scores=pd.DataFrame(scores),timing_anchor_validation=pd.DataFrame(validation),
               timing_biomarkers=pd.DataFrame(feature_rows))
    return out,"estimated"


def native_predictions(layer,date,method,hashes):
    scope=scope_for_read(layer,date)
    def obj(f): return json.loads(f.read_text()) if f.is_file() else {}
    dx=obj(scope/"qc/dx.json")
    if dx.get("endpoint_sha256")!=hashes.get(date) or not dx.get("endpoint_sha256"):
        return None
    directory=scope/"c1_correlate"
    final=obj(directory/"qc/C1_FINAL_INDEX.json")
    if not str(final.get("status","")).startswith("completed"): return None
    phases=[obj(scope/"logs"/(k+".status.json")) for k in ("analysis","prepare","primary","evaluation")]
    last=max((r for r in phases if r),key=lambda r:r.get("started",""),default={})
    if last.get("status") in ("RUN", "RUNNING", "EXIT", "FAILED", "INTERRUPTED", "SKIPPED", "BLOCKED"): return None
    info=final.get("methods",{}).get(method,{})
    if not info.get("path"): return None
    version=directory/"tmp/models"/method/Path(info["path"]).name
    sig=info.get("stored_signature") or info.get("signature")
    if not sig or obj(version/"COMPLETE.json").get("signature")!=sig: return None
    pred=version/"test_predictions_LOCAL_ONLY.tsv.gz"
    if not pred.is_file(): return None
    return pd.read_csv(pred,sep="\t",usecols=["eid","probability"],dtype={"eid":str}).set_index("eid",verify_integrity=True).probability


def prediction(d,raw,train,test,reps,layer,hashes,dates):
    # Training arms differ only by early first diagnoses without a prior or
    # undated history. All controls and later cases are shared.
    clean=d.any1_state.isin(["control","incident"]).to_numpy()
    early=(d.any1_state.eq("excluded_year1") & ~d.any_prior & ~d.any_undated).to_numpy()
    base=train & (clean|early)
    time=np.where(d.any1_state.eq("incident")|early,d.any_first_days,d.any_followup_days)/365.25
    event=(d.any1_state.eq("incident").to_numpy()|early).astype(int)
    status=[]; results=[]
    if not base.any() or not (train & clean).any() or not (test & clean).any():
        return pd.DataFrame(),[dict(analysis="controlled_training",status="insufficient_eligible_participants")]
    try:
        km=CensoringKM().fit(time[base],event[base],HORIZON,.01)
        yy,ww=km.labels_weights(time[base],event[base])
        y=np.zeros(len(d),int); w=np.zeros(len(d)); y[base]=yy; w[base]=ww
        # Condition the evaluation's censoring distribution on eligibility at the landmark.
        evalu=test & clean
        km_test=CensoringKM().fit(time[train&clean],event[train&clean],HORIZON,.01)
        y[evalu],w[evalu]=km_test.labels_weights(time[evalu],event[evalu])
    except ValueError as exc:
        if "Insufficient training censoring support" not in str(exc): raise
        return pd.DataFrame(),[dict(analysis="controlled_training",status=str(exc))]
    full=base & (w>0); filtered=full & clean; q=evalu & (w>0)
    def enough(mask): return min((y[mask]==0).sum(),(y[mask]==1).sum())>=2
    if not all(enough(s) for s in (full,filtered,q)):
        return pd.DataFrame(),[dict(analysis="controlled_training",status="insufficient_known_cases_or_controls")]
    x,_=frozen_scale(raw,train & clean)
    if x.shape[1] and (full&early).any():
        p0=fitted(x,full,y,w,q);p1=fitted(x,filtered,y,w,q)
        results.extend(compare(y[q],w[q],d.group[q].to_numpy(),p0,p1,reps,"controlled_training","any_all_first_dates","any1"))
        status.append(dict(analysis="controlled_training",status="estimated",training_full_N=int(full.sum()),
                  training_filtered_N=int(filtered.sum()),early_training_cases_removed=int((full&early).sum()),
                  test_eligible_N=int(evalu.sum()),test_known_N=int(q.sum()),horizon_from_baseline=HORIZON))
    else:
        status.append(dict(analysis="controlled_training",status="no_exclusive_early_cases_or_estimable_biomarkers"))
    if layer is not None:
        for model in ("glmnet","abm","abm_clinical","lightgbm"):
            pc=native_predictions(layer,"any1",model,hashes)
            for date in dict.fromkeys(["icd","ref",*[a for a in dates if a!="any1"]]):
                pr=native_predictions(layer,date,model,hashes)
                info=dict(analysis="saved_C1_common_test",method=model,reference=date,reference_endpoint_sha256=hashes.get(date))
                if pc is None or pr is None:
                    status.append(dict(**info,status="compatible_frozen_test_predictions_unavailable")); continue
                first=d.eid.map(pr).to_numpy(float); second=d.eid.map(pc).to_numpy(float)
                sub=q & np.isfinite(first) & np.isfinite(second)
                if not enough(sub):
                    status.append(dict(**info,status="insufficient_common_cases_or_controls")); continue
                paired = compare(y[sub],w[sub],d.group[sub].to_numpy(),first[sub],second[sub],reps,
                                 "saved_C1_common_test",date+"."+model,"any1."+model)
                results.extend(dict(row,reference_endpoint_sha256=hashes[date]) for row in paired)
                status.append(dict(**info,status="estimated",test_eligible_N=int(evalu.sum()),test_known_N=int(sub.sum())))
    return pd.DataFrame(results),status


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("layer-dir","snapshot","ygb-phe","roster","r-bin","trait","layer","dates"):
        p.add_argument("--"+key,required=True)
    p.add_argument("--group-file");p.add_argument("--group-col")
    p.add_argument("--cores",type=int,default=2)
    p.add_argument("--bootstrap",type=int,default=300)
    a=p.parse_args()
    if a.bootstrap < 20: p.error("--bootstrap requires at least 20 family draws")
    layer=Path(a.layer_dir); tmp=layer/"tmp/any1/c1_correlate"; qc=layer/"qc/any1/c1_correlate"
    dest=scope_for_read(layer,"any1")/"c1_correlate/timing.xlsx"
    for folder in (tmp,qc,dest.parent): folder.mkdir(parents=True,exist_ok=True)
    snapshot=Path(a.snapshot)
    meta=json.loads(snapshot.with_name("manifest.json").read_text())
    hashes={r["dx"]:r["endpoint_sha256"] for r in meta["counts"]}
    phe=Path(a.ygb_phe); roster=Path(a.roster)
    with file_lock(qc/"timing.lock"), tempfile.TemporaryDirectory(dir=tmp) as work, threadpool_limits(limits=max(1,a.cores)):
        work=Path(work)
        valid_input=next(r for r in meta["counts"] if r["dx"]=="any1")
        tables={}
        if valid_input["unavailable"]==valid_input["cohort_N"]:
            status=[dict(analysis="timing",status="endpoint_unavailable")]
        else:
            family=a.group_file or (None if a.group_col else str(phe/"rap/raw/ygb_kinship_groups.tsv"))
            if not roster.is_file():
                ensure_roster(ygb_phe=phe,output=roster,seed=SEED,group_file=family or "",
                              group_col=a.group_col or "",r_bin=a.r_bin)
            rec=roster_record(roster)
            roles=read_table(roster)[["eid","role"]]
            if rec.get("group_file"):
                families=_read_group_file(rec["group_file"],a.r_bin)[["eid","group"]]
            elif rec.get("group_col"):
                from c1_protocol import _cohort_projection, sha
                allrds=phe/"Rdata/all.rds"
                families=_cohort_projection(allrds,sha(allrds),a.r_bin,rec["group_col"])
            else:
                raise ValueError("The frozen C1 roster lacks family groups")
            roles=roles.merge(families,on="eid",validate="one_to_one",how="left")
            if roles.group.isna().any() or roles.groupby("group").role.nunique().max()>1:
                raise ValueError("Invalid family assignment in the frozen outer roster")
            cmd=[a.r_bin,str(HERE/"timing_matrix.R"),str(phe/"Rdata"/(a.layer+".rds")),str(work)]
            subprocess.run(cmd,check=True)
            ids=pd.read_csv(work/"ids.txt",header=None,dtype=str).iloc[:,0]
            names=(work/"features.txt").read_text().splitlines()
            x=np.memmap(work/"baseline.f32",dtype="<f4",mode="r",shape=(len(ids),len(names)),order="F")
            needed={"eid","icd_state","ref_state","any1_state","any_prior","any_undated",
                    "any_first_days","any_followup_days","any_early_days","age","sex","drug.lipid","baseline_lipid_i0"}
            phen=pd.read_csv(snapshot,sep="\t",usecols=lambda c:c in needed,dtype={"eid":str})
            aligned=pd.DataFrame({"eid":ids,"assay_row":np.arange(len(ids))}).merge(phen,on="eid",validate="one_to_one")
            d=aligned.merge(roles,on="eid",validate="one_to_one")
            if not len(d): raise ValueError("No assay participants matched the frozen cohort")
            raw=x if len(d)==len(x) and np.array_equal(d.assay_row,np.arange(len(x))) else x[d.assay_row]
            observed=np.zeros(len(d),bool)
            for k in range(0,len(names),32): observed |= np.isfinite(raw[:,k:k+32]).any(axis=1)
            if not observed.all():
                raw=raw[observed];d=d.loc[observed].reset_index(drop=True)
            if not d.role.isin(["training","test"]).all(): raise ValueError("Unknown outer role")
            train=d.role.eq("training").to_numpy();test=~train
            tables,msg=profiles(d,raw,names,train,test,work/"scores.tsv.gz")
            pair,status=prediction(d,raw,train,test,a.bootstrap,layer,hashes,a.dates.split(","))
            tables["timing_paired"]=pair
            status.append(dict(analysis="blinded_profiles",status=msg,assay_N=len(ids),overlap_N=len(d),
                        rows_without_assay_values=int((~observed).sum()),family_roster="frozen C1 outer roles"))
        tables["timing_status"]=pd.DataFrame(status)
        with pd.ExcelWriter(work/"timing.xlsx",engine="openpyxl") as book:
            for name,df in tables.items():
                df.assign(trait=a.trait,layer=a.layer,dx="any1",endpoint_sha256=hashes["any1"]).to_excel(book,sheet_name=name,index=False)
        saved=tmp/"timing_scores_LOCAL_ONLY.tsv.gz"
        if (work/"scores.tsv.gz").exists():
            (work/"scores.tsv.gz").chmod(0o600);os.replace(work/"scores.tsv.gz",saved)
        else: saved.unlink(missing_ok=True)
        os.replace(work/"timing.xlsx",dest)
        print(dest)


if __name__=="__main__":
    main()
