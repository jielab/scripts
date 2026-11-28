"""Outcome-blind prediction of 10-year risk and expected Brier loss for each at-risk participant."""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.exceptions import ConvergenceWarning
from abm import CensoringKM
from c1_protocol import stable_hash
from timing_study import fitted, HORIZON

PRIOR = 10.
MAX_DONORS = 32


def calibrated_errors(d, query, donors, risk, y, w, geo, sex, valid, cores):
    """Learn a risk-only expected error and a locally weighted error on independent donor families."""
    edges=np.unique(np.quantile(risk[donors],np.arange(1,10)/10))
    bins=np.searchsorted(edges,risk,side="right")
    known=donors[w[donors]>0]
    loss=(y[known]-risk[known])**2
    prior=np.average(loss,weights=w[known])
    size=len(edges)+1
    mass=np.bincount(bins[known],weights=w[known],minlength=size)
    sums=np.bincount(bins[known],weights=w[known]*loss,minlength=size)
    reference=(sums+PRIOR*prior)/(mass+PRIOR)
    background=reference[bins[query]]
    expected=background.copy()
    used=np.zeros(len(query),int); ess=np.zeros(len(query))
    # Eligible outcome donors come from other families than all queries. Tie choices use eid only.
    loc=known[valid[known]]
    order=np.argsort(d.eid.to_numpy()[loc],kind="stable")
    loc=loc[order];loc=loc[~d.group.iloc[loc].duplicated().to_numpy()]
    for b in np.unique(bins[query]):
        for code in (0,1):
            bank=loc[(bins[loc]==b) & (sex[loc]==code)]
            target=np.flatnonzero(valid[query] & (bins[query]==b) & (sex[query]==code))
            if not len(bank) or not len(target):continue
            tree=cKDTree(geo[bank])
            ranks=list(range(1,min(MAX_DONORS,len(bank))+1))
            for start in range(0,len(target),4096):
                ix=target[start:start+4096]
                dist,j=tree.query(geo[query[ix]],k=ranks,p=np.inf,workers=cores,
                                  distance_upper_bound=np.nextafter(1.,np.inf))
                ok=np.isfinite(dist); near=bank[np.minimum(j,len(bank)-1)]; weight=w[near]*ok
                total=weight.sum(axis=1);used[ix]=ok.sum(axis=1)
                eff=np.divide(total**2,(weight**2).sum(axis=1),out=np.zeros(len(ix)),where=total>0)
                local=np.divide((weight*(y[near]-risk[near])**2).sum(axis=1),total,
                                out=np.zeros(len(ix)),where=total>0)
                expected[ix]=(eff*local+PRIOR*background[ix])/(eff+PRIOR)
                ess[ix]=eff
    return expected,background,used,ess,bins[query]


def evaluation(rows,y,w,groups,reps,seed):
    test=rows.role.eq("test") & rows.score_status.eq("scored")
    known=test & (w>0)
    q=np.flatnonzero(known)
    calibrated=[]
    for band in ("all","lower","middle","higher"):
        elig=test if band=="all" else test & rows.reliability_band.eq(band)
        ix=np.flatnonzero(elig & (w>0))
        row=dict(band=band,eligible_N=int(elig.sum()),known_N=len(ix),cases=int(y[ix].sum()))
        if len(ix):
            risk=rows.risk_10y.to_numpy()[ix]; loss=(risk-y[ix])**2
            for column,key in (("risk_10y","mean_risk"),("predictivity","mean_predictivity"),
                               ("expected_Brier","mean_expected_Brier"),("risk_only_Brier","mean_risk_only_Brier")):
                row[key]=np.average(rows[column].to_numpy()[ix],weights=w[ix])
            row.update(observed_Brier=np.average(loss,weights=w[ix]),neighbor_coverage=float((rows.donor_families.to_numpy()[ix]>0).mean()))
        calibrated.append(row)
    if not len(q):return pd.DataFrame(calibrated),pd.DataFrame([dict(status="no_known_test_outcomes")])
    z=rows.iloc[q]; weight=w[q];err=(z.risk_10y.to_numpy()-y[q])**2
    # Test whether the local estimator predicts realized loss better than an estimator using risk alone.
    gain=(err-z.expected_Brier.to_numpy())**2-(err-z.risk_only_Brier.to_numpy())**2
    bins=z.risk_bin.to_numpy(int); levels=z.within_risk_level.to_numpy()
    hi=levels=="higher";lo=levels=="lower";size=int(bins.max()+1)
    mh=np.bincount(bins[hi],weights=weight[hi],minlength=size)
    ml=np.bincount(bins[lo],weights=weight[lo],minlength=size)
    shared=(mh>0)&(ml>0)
    target=(mh+ml)[shared]
    if len(target):target/=target.sum()
    def estimates(wt):
        added=np.average(gain,weights=wt)
        means=[]
        for m in (hi,lo):
            mass=np.bincount(bins[m],weights=wt[m],minlength=size)[shared]
            if not len(mass) or (mass<=0).any():return added,np.nan
            losses=np.bincount(bins[m],weights=wt[m]*err[m],minlength=size)[shared]
            means.append(target@(losses/mass))
        return added,means[0]-means[1]
    point=estimates(weight)
    _,fam=np.unique(groups[q],return_inverse=True);n=fam.max()+1;rng=np.random.default_rng(seed)
    boot=[]
    for _ in range(reps):
        mult=np.bincount(rng.integers(n,size=n),minlength=n)
        boot.append(estimates(weight*mult[fam]))
    out=[]
    for j,name in enumerate(("local_minus_risk_only_error_prediction_MSE","higher_minus_lower_observed_Brier_within_risk")):
        samp=np.asarray(boot)[:,j];samp=samp[np.isfinite(samp)]
        ci=np.quantile(samp,[.025,.975]) if len(samp)>=max(20,.8*reps) else [np.nan,np.nan]
        selected=np.ones(len(q),bool) if j==0 else (hi|lo)&shared[bins]
        out.append(dict(contrast=name,estimate=point[j],lower=ci[0],upper=ci[1],
                        N=int(selected.sum()),cases=int(y[q][selected].sum()),
                        families=len(np.unique(groups[q][selected])),shared_risk_bins=int(shared.sum()),
                        valid_bootstraps=len(samp),status="estimated" if np.isfinite(point[j]) else "no_within_risk_separation"))
    return pd.DataFrame(calibrated),pd.DataFrame(out)


def score_roster(d, eligible, count, valid, unavailable_reason):
    """Retain every assayed person with an explicit reason when a risk score is unavailable."""
    rows=d[["eid","role"]].copy()
    rows["score_status"]=np.where(eligible,unavailable_reason,"outside_ICD_incident_risk_set")
    for name in ("risk_10y","expected_Brier","predictivity","risk_only_Brier","family_ESS","donor_families","risk_bin"):
        rows[name]=np.nan
    rows["reference_twin_families"]=count;rows["geometry_complete"]=valid
    rows["model_rotation"]=-1;rows["reliability_band"]="unavailable";rows["within_risk_level"]="unavailable"
    return rows


def unavailable_predictivity(d, eligible, folder, reason):
    rows=score_roster(d,eligible,np.nan,False,reason)
    rows.to_csv(folder/"predictivity.tsv.gz",sep="\t",index=False,compression="gzip")
    validation,tests=evaluation(rows,np.zeros(len(d)),np.zeros(len(d)),d.group.to_numpy(),0,0)
    return {"predictivity":validation,"predictivity_test":tests,
            "predictivity_fit":pd.DataFrame([dict(status=reason)])}


def personal_predictivity(d,full_x,geo,sex,valid,count,eligible,folder,cores,seed,reps):
    """In each rotation two development folds fit risk, a third supplies honest losses, and the fourth is scored."""
    train=d.role.eq("training").to_numpy();test=d.role.eq("test").to_numpy()
    fold=d.group.map({g:stable_hash(g,seed+3109)%4 for g in d.group.unique()}).to_numpy()
    time=d.icd_t2e.to_numpy(float);event=d.icd_state.eq("incident").to_numpy(int)
    rows=score_roster(d,eligible,count,valid,"insufficient_development_support")
    known_y=np.zeros(len(d),int);known_w=np.zeros(len(d));audit=[]
    for rotation in range(4):
        donor_fold=(rotation+1)%4
        fit=train & eligible & (fold!=rotation)&(fold!=donor_fold)
        donor=np.flatnonzero(train & eligible & (fold==donor_fold))
        check=train & eligible & (fold==rotation)
        query=np.flatnonzero(check | (test & eligible & (rotation==0)))
        info=dict(rotation=rotation,fit_eligible_N=int(fit.sum()),donor_eligible_N=len(donor),score_N=len(query))
        if not fit.any() or not len(donor) or not len(query):
            audit.append(dict(info,status="empty_partition"));continue
        km=CensoringKM().fit(time[fit],event[fit],HORIZON,.01)
        y=np.zeros(len(d),int);w=np.zeros(len(d))
        try:
            y[fit],w[fit]=km.labels_weights(time[fit],event[fit])
            y[donor],w[donor]=km.labels_weights(time[donor],event[donor])
        except ValueError as exc:
            if "Insufficient training censoring support" not in str(exc):raise
            audit.append(dict(info,status="insufficient_censoring_support"));continue
        build=fit & (w>0);known=donor[w[donor]>0]
        for label,idx in (("fit",build),("donor",known)):
            info[label+"_cases"]=int(y[idx].sum());info[label+"_controls"]=int((y[idx]==0).sum())
        if min(info["fit_cases"],info["fit_controls"],info["donor_cases"],info["donor_controls"])<5:
            audit.append(dict(info,status="insufficient_known_cases_or_controls"));continue
        ix=np.r_[donor,query];risk=np.full(len(d),np.nan)
        try:
            risk[ix]=fitted(full_x,build,y,w,ix)
        except ConvergenceWarning:
            audit.append(dict(info,status="risk_ridge_did_not_converge"))
            continue
        e,b,n,ess,bins=calibrated_errors(d,query,donor,risk,y,w,geo,sex,valid,cores)
        reliability=1-e
        rows.loc[query,["risk_10y","expected_Brier","predictivity","risk_only_Brier",
                        "donor_families","family_ESS","risk_bin"]]=np.column_stack([risk[query],e,reliability,b,n,ess,bins])
        rows.loc[query,"model_rotation"]=rotation
        rows.loc[query,"score_status"]="scored"
        # Define validation bands on held-out development features, independently of all test features/outcomes.
        calibration=check[query]
        if calibration.any():
            cut=np.quantile(reliability[calibration],[1/3,2/3])
            rows.loc[query,"reliability_band"]=np.select([reliability<cut[0],reliability>cut[1]],
                                                         ["lower","higher"],default="middle")
            info.update(band_cut_low=cut[0],band_cut_high=cut[1])
            for bnum in np.unique(bins[calibration]):
                sub=bins==bnum;med=np.median(reliability[sub & calibration])
                rows.loc[query[sub],"within_risk_level"]=np.select([reliability[sub]>med,reliability[sub]<med],
                                                                   ["higher","lower"],default="tied")
        audit.append(dict(info,status="estimated",donors_with_measured_geometry=int(valid[known].sum())))
        if rotation==0:
            outer_km=km
    # Outcome read for test evaluation occurs after every risk and predictivity score is frozen.
    rows.to_csv(folder/"predictivity.tsv.gz",sep="\t",index=False,compression="gzip",float_format="%.17g")
    qt=np.flatnonzero(test & rows.score_status.eq("scored").to_numpy())
    if len(qt):
        try:known_y[qt],known_w[qt]=outer_km.labels_weights(time[qt],event[qt])
        except ValueError as exc:
            if "Insufficient training censoring support" not in str(exc):raise
            audit.append(dict(rotation=0,status="test_IPCW_unavailable"))
    validation,contrasts=evaluation(rows,known_y,known_w,d.group.to_numpy(),reps,seed)
    return {"predictivity":validation,"predictivity_test":contrasts,"predictivity_fit":pd.DataFrame(audit)}
