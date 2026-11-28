"""ICD-10 first-year cases classified by other sources of pre-baseline records."""
import numpy as np
import pandas as pd
from scipy.stats import t
from sklearn.exceptions import ConvergenceWarning
from cad_icd_first import contrast
from timing_study import fitted, frozen_scale

GROUPS = ("no_prior_any_source","REF_prior","GP_prior","SRD_prior")


def first_year_groups(d):
    canonical = d.icd_state.eq("incident") & d.icd_days.between(1,365)
    earlier_icd = d.src_icd10_days.le(0)
    eligible = canonical & ~earlier_icd & ~d.src_icd10_undated
    flags = pd.DataFrame({
        GROUPS[0]: d.any_first_days.gt(0) & ~d.any_undated,
        GROUPS[1]: d.ref_days.lt(0),
        GROUPS[2]: d.src_gp_days.lt(0),
        GROUPS[3]: d.src_srd_days.lt(0)
    },index=d.index) & eligible.to_numpy()[:,None]
    return eligible, flags


def feature_profiles(d,raw,names,flags):
    """Use measured values. Cluster the difference in group means by frozen family."""
    sel=flags.any(axis=1).to_numpy()
    if not sel.any():return pd.DataFrame()
    x=raw[sel]
    masks=flags.loc[sel].to_numpy()
    fam=pd.factorize(d.group[sel],sort=False)[0]; nf=fam.max()+1
    rows=[]
    for j,name in enumerate(names):
        v=np.asarray(x[:,j],float)
        summaries=[];influence=[];families=[]
        for g,col in enumerate(GROUPS):
            mask=masks[:,g] & np.isfinite(v)
            obs=v[mask]; n=len(obs); ids=fam[mask]
            mean=obs.mean() if n else np.nan
            sd=obs.std(ddof=1) if n>1 else np.nan
            quant=np.quantile(obs,[.25,.5,.75]) if n else (np.nan,)*3
            summaries.append(dict(group=col,feature=name,N=n,mean=mean,SD=sd,q25=quant[0],median=quant[1],q75=quant[2]))
            influence.append(np.bincount(ids,weights=(obs-mean)/max(n,1),minlength=nf))
            families.append(np.bincount(ids,minlength=nf)>0)
        ref=summaries[0]
        for g,row in enumerate(summaries):
            row["mean_z_vs_no_prior"]=(row["mean"]-ref["mean"])/ref["SD"] if ref["SD"]>0 else np.nan
            if g and min(row["N"],ref["N"])>=2:
                k=(families[0]|families[g]).sum()
                if k>=2:
                    diff=row["mean"]-ref["mean"]
                    se=np.sqrt(k/(k-1)*np.square(influence[g]-influence[0]).sum())
                    row.update(mean_difference=diff,families=int(k))
                    if se>0:
                        ci=t.ppf(.975,k-1)*se
                        row.update(lower=diff-ci,upper=diff+ci,p=2*t.sf(abs(diff/se),k-1))
            rows.append(row)
    df=pd.DataFrame(rows)
    df["FDR"]=np.nan
    if "p" in df:
        p=df.p.dropna().sort_values(); n=len(p)
        df.loc[p.index,"FDR"]=np.minimum(1,np.minimum.accumulate((p.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1])
    return df


def four_groups(d,raw,names,train,test,folder,reps,seed):
    eligible,flags=first_year_groups(d)
    # Use strict pre-baseline evidence. Day-zero and unresolved-only records cannot enter group one.
    entered=eligible.to_numpy()
    combos=flags.astype(int).assign(eligible=entered)
    coverage=[]
    for split,subset in (("all_assayed",np.ones(len(d),bool)),("outer_test",test)):
        group_rows=(d.icd_state.eq("incident") & d.icd_days.between(1,365)).to_numpy() & subset
        for k,mask in (("canonical_ICD_days_1_365",group_rows),("same_source_earlier_or_day_zero",group_rows & d.src_icd10_days.le(0)),
                       ("undated_ICD_only_exclusion",group_rows & d.src_icd10_days.gt(0) & d.src_icd10_undated),
                       ("eligible",entered & subset),("distinct_in_four_groups",subset & flags.any(axis=1)),
                       ("overlap_two_or_three_prior_sources",subset & flags.sum(axis=1).gt(1)),
                       ("other_only_day_zero_or_unresolved",subset & entered & ~flags.any(axis=1))):
            coverage.append(dict(split=split,stage=k,N=int(mask.sum())))
    overlap=combos.loc[eligible,list(GROUPS)].value_counts().rename("N").reset_index()
    groups=pd.DataFrame(dict(group=GROUPS,N_all_assayed=flags.sum().to_numpy(),
                             N_outer_test=flags.loc[test].sum().to_numpy()))
    result={"four_groups":groups,"four_flow":pd.DataFrame(coverage),"four_overlap":overlap,
            "four_features":feature_profiles(d,raw,names,flags)}
    # Keep the target cases out of both sides of the training contrast.
    prev=(d.icd_state.eq("prevalent") & d.icd_days.lt(0) & ~d.src_icd10_undated).to_numpy()
    late=(d.icd_state.eq("incident") & d.icd_days.between(366,3652) &
          d.src_icd10_days.gt(365) & ~d.src_icd10_undated).to_numpy()
    fit=train & (prev|late); query=test & (entered|prev|late)
    counts=dict(train_prevalent=int((train & prev).sum()),train_later=int((train & late).sum()),
                test_prevalent=int((test & prev).sum()),test_later=int((test & late).sum()))
    def finish(status,**extra):
        result["four_anchors"]=pd.DataFrame([dict(status=status,**counts,**extra)])
        return result
    if min(counts["train_prevalent"],counts["train_later"])<10 or not query.any():
        return finish("insufficient_anchors_or_test")
    # Scaling and fitting need only anchor training and requested test rows.
    subset=np.flatnonzero(fit|query)
    x,keep=frozen_scale(raw[subset],fit[subset])
    if not len(keep):return finish("no_estimable_features")
    weights=np.ones(len(d))
    for label in (prev,late):weights[fit & label]=fit.sum()/(2*(fit & label).sum())
    score=np.full(len(d),np.nan)
    try:
        score[query]=fitted(x,fit[subset],prev[subset].astype(int),weights[subset],query[subset])
    except ConvergenceWarning:
        return finish("anchor_ridge_did_not_converge",features=len(keep))
    local=d.loc[query,["eid","group"]].rename(columns={"group":"family"}).join(flags.loc[query])
    local["anchor_or_query"]=np.select([prev[query],late[query]],["prevalent_anchor","later_anchor"],default="first_year")
    local["similarity_score"]=score[query]
    local.to_csv(folder/"four_scores.tsv.gz",sep="\t",index=False,compression="gzip")
    stats=[]
    for j,name in enumerate(GROUPS):
        mask=flags[name].to_numpy() & test; p=score[mask]
        row=dict(group=name,score_q25=np.nan,score_median=np.nan,score_q75=np.nan)
        if len(p):row["score_q25"],row["score_median"],row["score_q75"]=np.quantile(p,[.25,.5,.75])
        if j:
            ix=mask|(test & flags[GROUPS[0]].to_numpy())
            est=contrast(score[ix],mask[ix],d.group[ix],0.,reps,seed,"AUC")
            row.update(AUC_vs_no_prior=est.get("estimate"),AUC_lower=est.get("lower"),AUC_upper=est.get("upper"),
                       AUC_status=est["status"])
        stats.append(row)
    result["four_groups"]=groups.merge(pd.DataFrame(stats),on="group",validate="one_to_one")
    ix=test & (prev|late)
    test_auc=contrast(score[ix],prev[ix],d.group[ix],0.,reps,seed,"AUC")
    anchor_scores={}
    for label,mask in (("prevalent",prev),("later",late)):
        values=score[test & mask]
        if len(values):
            for key,value in zip(("q25","median","q75"),np.quantile(values,[.25,.5,.75])):
                anchor_scores[label+"_"+key]=value
    return finish("estimated",**anchor_scores,features=len(keep),anchor_AUC=test_auc.get("estimate"),
                  anchor_lower=test_auc.get("lower"),anchor_upper=test_auc.get("upper"))
