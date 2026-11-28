"""Baseline metabolite scores in people with first recorded ICD-10 diagnosis on days 1–365."""
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from c1_protocol import stable_hash
from timing_study import fitted, frozen_scale

OTHER = ("gp","srd","icd9","procedure","ref","death")
NO_PRIOR = "no_prior_record_found"


def other_history(d):
    dated=d[["src_"+s+"_days" for s in OTHER]]
    raw_gp=d.src_raw_gp_days.lt(0)
    pre=dated.lt(0).any(axis=1)
    same_day=dated.eq(0).any(axis=1)
    group=np.select([raw_gp,d.src_gp_days.lt(0),d.src_icd9_days.lt(0),d.src_procedure_days.lt(0),
                     d.src_srd_days.lt(0),d.src_ref_days.lt(0),pre,same_day,d.any_undated],
                    ["raw_GP_prior","GP_cache_or_other_definition_prior","ICD9_prior","procedure_prior",
                     "self_report_prior","REF_prior","other_prior","day_zero_only","unresolved_timing"],default=NO_PRIOR)
    return group


def contrast(score,label,group,cut,reps,seed,method):
    """Conditional family bootstrap on disjoint test groups; threshold and fit stay frozen."""
    group=np.asarray(group);label=np.asarray(label,bool);score=np.asarray(score,float)
    n1=int(label.sum());n0=len(label)-n1
    unique,inv=np.unique(group,return_inverse=True); n=len(unique)
    row=dict(estimator=method,positive_N=n1,reference_N=n0,families=n)
    if min(n0,n1)<2 or n<2:return dict(row,status="insufficient_test_groups")
    high=(score>=cut).astype(float)
    order=np.argsort(score,kind="stable")
    start=np.r_[0,np.flatnonzero(np.diff(score[order]))+1]
    def est(w):
        a=w*label;b=w*(~label)
        if min(a.sum(),b.sum())==0:return np.nan
        if method=="fraction_above_cut_difference":
            return a@high/a.sum()-b@high/b.sum()
        pos=np.add.reduceat(a[order],start);neg=np.add.reduceat(b[order],start)
        return pos@(np.cumsum(neg)-.5*neg)/(a.sum()*b.sum())
    rng=np.random.default_rng(seed);boot=[]
    for _ in range(reps):
        mult=np.bincount(rng.integers(n,size=n),minlength=n)
        boot.append(est(mult[inv]))
    boot=np.asarray(boot);boot=boot[np.isfinite(boot)]
    lo,hi=np.quantile(boot,[.025,.975]) if len(boot)>=max(20,reps*.8) else (np.nan,np.nan)
    return dict(row,status="estimated",estimate=est(np.ones(len(score))),lower=lo,upper=hi,valid_bootstraps=len(boot))


def analyze_icd_first(d,raw,train,test,folder,reps,seed):
    required={"src_icd10_days","src_icd10_undated","src_raw_gp_days","any_early_sources",
              *("src_"+s+"_days" for s in OTHER)}
    if not required<=set(d):raise ValueError("Source dates missing. Run the updated AOD preparation with the C1 workflow.")
    canonical=d.icd_state.eq("incident") & d.icd_days.between(1,365)
    conflict=d.src_icd10_days.le(0)
    invalid=d.src_icd10_undated
    early=(canonical & ~conflict & ~invalid).to_numpy()
    # Anchors use ICD-10 alone. Other sources remain concealed until scores and cut are fixed.
    prev=(d.icd_state.eq("prevalent") & d.icd_days.lt(-365) & ~invalid).to_numpy()
    late=(d.icd_state.eq("incident") & d.icd_days.between(366,3652) & d.src_icd10_days.gt(365) & ~invalid).to_numpy()
    cal_family=d.group.map({g:stable_hash(g,seed+711)%5==0 for g in d.group.unique()}).to_numpy()
    build=train & ~cal_family & (prev|late); cal=train & cal_family & late
    history=other_history(d)
    old=(test & d.any1_state.eq("excluded_year1") & d.any_prior).to_numpy()
    old_audit=pd.DataFrame({"saved_ICD_first":np.select([canonical,d.icd_state.eq("prevalent"),d.icd_state.eq("incident")],
                             ["days_1_365","baseline_or_earlier","after_year_one"],default="other"),
                            "ICD_any_version_prior_or_day0":conflict,"first_year_sources":d.any_early_sources.fillna(""),
                            "other_history":history}).loc[old]
    result={"icd_any1_overlap":old_audit.value_counts().rename("N").reset_index()}
    flows=[]
    for name,mask in (("analysis_cohort",np.ones(len(d),bool)),("heldout_test",test)):
        for stage,rows in (("all_met",mask),("canonical_first_ICD_days_1_365",mask & canonical),
                           ("earlier_or_day0_ICD_version",mask & canonical & conflict),
                           ("unresolved_ICD_no_dated_conflict",mask & canonical & ~conflict & invalid),
                           ("eligible_first_ICD",mask & early)):
            flows.append(dict(split=name,stage=stage,N=int(rows.sum())))
    result["icd_entry_flow"]=pd.DataFrame(flows)
    query=test & (prev|late|early)
    primary=test & early
    signatures=pd.DataFrame({s+"_pre":d["src_"+s+"_days"].lt(0) for s in ("raw_gp",*OTHER)})
    signatures["any_day0"]=d[["src_"+s+"_days" for s in OTHER]].eq(0).any(axis=1)
    signatures["any_undated"]=d.any_undated
    result["icd_prior_sources"]=signatures.loc[primary].value_counts().rename("N").reset_index()
    summary=dict(build_prevalent_N=int((build & prev).sum()),build_later_N=int((build & late).sum()),
                 threshold_later_N=int(cal.sum()),eligible_test_N=int(primary.sum()))
    def status(msg,**fields):
        result["icd_score_status"]=pd.DataFrame([dict(summary,status=msg,**fields)])
        return result
    if min((build & prev).sum(),(build & late).sum(),cal.sum())<10 or not query.any():
        return status("insufficient_independent_ICD_anchors_or_test")
    z,features=frozen_scale(raw,build)
    if not z.shape[1]:return status("no_estimable_metabolites")
    weights=np.ones(len(d))
    for cls in (prev,late):
        weights[build & cls]=build.sum()/(2*(build & cls).sum())
    score=np.full(len(d),np.nan);predict=cal|query
    score[predict]=fitted(z,build,prev.astype(int),weights,predict)
    cut=float(np.quantile(score[cal],.90))
    # No outcome-record fields enter this classifier or the training-family cutoff.
    period=np.select([prev,early,late],["prevalent_anchor","first_year_ICD","later_ICD"],default="other")
    s=d.loc[query,["eid","group","icd_days"]].assign(period=period[query],history=history[query],
                  score=score[query],above_cut=score[query]>=cut)
    meds=pd.to_numeric(d.get("baseline_lipid_i0",pd.Series(np.nan,index=d.index)),errors="coerce")
    s["lipid_treatment"]=np.select([meds[query].eq(0),meds[query].eq(1)],["nonuser","user"],default="unknown")
    s.to_csv(folder/"icd_first.tsv.gz",sep="\t",index=False,compression="gzip")
    anchors=s.period.isin(["prevalent_anchor","later_ICD"])
    y=s.period[anchors].eq("prevalent_anchor")
    anchor_auc=roc_auc_score(y,s.score[anchors]) if y.nunique()==2 else np.nan
    rows=[];tests=[]
    for med,data in [("all",s),*list(s.groupby("lipid_treatment"))]:
        for (period_name,hist),v in data.groupby(["period","history"]):
            rows.append(dict(period=period_name,history=hist,baseline_lipid=med,N=len(v),high_N=int(v.above_cut.sum()),
                     fraction_high=v.above_cut.mean(),q25=v.score.quantile(.25),median=v.score.median(),q75=v.score.quantile(.75)))
        e=data.period.eq("first_year_ICD"); no=data.history.eq(NO_PRIOR)
        positive=e & data.history.eq("raw_GP_prior")
        for label,yes,ref,method in (
             ("early_no_record_vs_later_no_record",e & no,data.period.eq("later_ICD") & no,"fraction_above_cut_difference"),
             ("raw_GP_prior_vs_no_record_within_early_ICD",positive,e & no,"record_evidence_AUC")):
            use=yes|ref;part=data[use]
            tests.append(dict(comparison=label,baseline_lipid=med,
                 **contrast(part.score,yes[use],part.group,cut,reps,seed,method)))
    result["icd_first_scores"]=pd.DataFrame(rows)
    result["icd_first_tests"]=pd.DataFrame(tests)
    return status("estimated",cut=cut,threshold_percentile=90,features=len(features),heldout_anchor_AUC=anchor_auc)
