"""Fit metabolic timing anchors and quantify the early-case contrast in the selected trait."""
import numpy as np
import pandas as pd
from scipy.stats import t
from timing_study import profiles

EARLY=("early_day_1_90","early_day_91_180","early_day_181_365")


def evidence(d,raw,names,train,test,path,reps,seed):
    tables,msg=profiles(d,raw,names,train,test,path)
    keep=("timing_scores","timing_anchor_validation","timing_medication","timing_groups","timing_biomarkers")
    result={"original_"+k.removeprefix("timing_"):tables[k] for k in keep if k in tables}
    if "original_anchor_validation" in result:result["original_validation"]=result.pop("original_anchor_validation")
    result["original_status"]=pd.DataFrame([dict(analysis="original_any1_anchor_protocol",status=msg)])
    if not path.is_file():return result
    s=pd.read_csv(path,sep="\t",dtype={"eid":str}).merge(d,on="eid",validate="one_to_one")
    s=s[s.profile_group.isin(["early_with_prior",*EARLY])].copy()
    s["prior"]=s.profile_group.eq("early_with_prior")
    meds=pd.to_numeric(s.get("baseline_lipid_i0",pd.Series(np.nan,index=s.index)),errors="coerce")
    rows=[]
    for stratum,mask in (("all",s.prior.notna()),("prior_also_dated",~s.any_undated),
                         ("baseline_nonuser",meds.eq(0)),("baseline_user",meds.eq(1))):
        a=s[mask]
        if a.prior.nunique()<2:continue
        g,inverse=np.unique(a.group,return_inverse=True); z=a.omics.to_numpy(); prior=a.prior.to_numpy()
        idx0=np.argsort(z[~prior]); idx1=np.argsort(z[prior])
        z0,z1=z[~prior][idx0],z[prior][idx1];g0,g1=inverse[~prior][idx0],inverse[prior][idx1]
        def median(vals,weight):
            cs=np.cumsum(weight); n=cs[-1]
            if n==0:return np.nan
            return (vals[np.searchsorted(cs,(n+1)//2)]+vals[np.searchsorted(cs,n//2+1)])/2
        rng=np.random.default_rng(seed);dif=[]
        for _ in range(reps):
            freq=np.bincount(rng.integers(len(g),size=len(g)),minlength=len(g))
            dif.append(median(z1,freq[g1])-median(z0,freq[g0]))
        valid=np.asarray(dif);valid=valid[np.isfinite(valid)]
        ci=np.quantile(valid,[.025,.975]) if len(valid)>=max(20,.8*reps) else [np.nan,np.nan]
        rows.append(dict(stratum=stratum,prior_N=len(z1),no_prior_N=len(z0),families=len(g),
                    prior_median=np.median(z1),no_prior_median=np.median(z0),
                    median_difference=np.median(z1)-np.median(z0),lower=ci[0],upper=ci[1],valid_bootstraps=len(valid)))
    result["early_median_contrast"]=pd.DataFrame(rows)
    # The molecular score is frozen. Model its mean among early cases with and without a prior record.
    s["age_numeric"]=pd.to_numeric(s.get("age"),errors="coerce")
    s["sex_numeric"]=s.get("sex",pd.Series(index=s.index,dtype=float)).astype(str).str.lower().map(
        {"male":1,"m":1,"1":1,"1.0":1,"female":0,"f":0,"0":0,"0.0":0})
    s["lipid_i0"]=meds.where(meds.isin([0,1]))
    rows=[]
    for cov in ([],["age_numeric","sex_numeric"],["age_numeric","sex_numeric","lipid_i0"]):
        cols=["prior",*cov]; a=s.dropna(subset=cols)
        x=np.column_stack([np.ones(len(a)),a[cols].to_numpy(float)])
        n,k=x.shape;groups,g=np.unique(a.group,return_inverse=True);m=len(groups)
        if n<=k or m<=k or np.linalg.matrix_rank(x)<k:continue
        z=a.omics.to_numpy(); beta=np.linalg.lstsq(x,z,rcond=None)[0]
        meat=np.zeros((m,k));np.add.at(meat,g,x*(z-x@beta)[:,None])
        bread=np.linalg.inv(x.T@x);var=bread@(meat.T@meat)@bread*m/(m-1)*(n-1)/(n-k)
        se=np.sqrt(max(0,var[1,1]));radius=t.ppf(.975,m-1)*se
        rows.append(dict(adjustment=",".join(cov) or "none",N=n,families=m,prior_N=int(a.prior.sum()),
                         mean_score_difference=beta[1],lower=beta[1]-radius,upper=beta[1]+radius))
    result["early_adjusted_contrast"]=pd.DataFrame(rows)
    return result
