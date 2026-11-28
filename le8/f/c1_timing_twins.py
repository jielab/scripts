"""Family-held-out, ten-year timing and molecular-twin experiments on a single ICD outcome."""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits

from abm import CensoringKM
from c1_protocol import ensure_roster, roster_record, stable_hash, _cohort_projection, sha, _read_group_file
from file_lock import file_lock
from timing_study import HORIZON, HERE, read_table, frozen_scale, fitted, profiles, compare, metrics, auc_layout

ROUTINE = ("bmi", "sbp", "dbp", "bb_TC", "bb_HDL", "bb_HBA1C")
SEX = {"male":1, "m":1, "1":1, "1.0":1, "female":0, "f":0, "0":0, "0.0":0}
SHEETS = ("study_status","exclusions","fit_test","contrasts","selection","twin_geometry",
          "twin_counts","timing_groups","timing_scores","timing_anchor_validation","timing_medication","timing_biomarkers")
METRICS = ("AUC_IPCW","Brier_IPCW","LogLoss_IPCW","predicted_minus_observed_risk")


def timing(d):
    """Apply exclusions without changing ICD times or events for retained people."""
    broad = d.icd_state.isin(["incident","control"]) & d.icd_t2e.gt(0)
    early = d.any_early_days.between(1,365)
    near_before = d.any_near_before_days.between(-365,0)
    observed = d.any_first_days.le(d.any_followup_days)
    cross = d.any_prior & d.icd_state.isin(["incident","control"])
    other_case = (d.icd_state.eq("control") | d.icd_t2e.gt(HORIZON)) & observed & d.any_first_days.between(1,365.25*HORIZON)
    gap = d.source_first_span_days.gt(365)
    flags = {"first_year_any_source":early, "year_before_or_day_zero":near_before,
             "prior_vs_ICD_at_risk":cross, "other_case_vs_ICD_10y_none":other_case,
             "source_first_span_gt365":gap, "unresolved_date":d.any_undated,
             "short_followup":d.any_followup_days.le(365), "sources_unavailable":d.any1_state.eq("unavailable")}
    clean = broad & ~pd.DataFrame(flags).any(axis=1)
    prev = d.icd_days.lt(-365) & d.ref_days.lt(-365) & ~d.any_undated & ~early & ~near_before & ~gap
    later = clean & d.icd_days.gt(365) & d.ref_days.gt(365) & d.icd_t2e.le(HORIZON)
    groups = np.select([prev,later,early & d.any_prior,near_before,early & d.any_undated,
                        early & d.any_early_days.le(90),early & d.any_early_days.le(180),early,
                        cross,other_case,gap,d.any_undated,clean & d.icd_state.eq("control")],
                ["prevalent_anchor","later_incident_anchor","early_with_prior","year_before_or_day_zero",
                 "early_with_undated","early_day_1_90","early_day_91_180","early_day_181_365",
                 "other_source_prior","other_source_case_vs_ICD_10y_none","source_gap_gt365",
                 "uncertain_timing","control"],default="other")
    return broad.to_numpy(), clean.to_numpy(), flags, groups


def representatives(d,mask):
    loc = np.flatnonzero(mask)
    order = np.argsort(d.eid.to_numpy()[loc],kind="stable")
    loc = loc[order]
    return loc[~d.group.iloc[loc].duplicated().to_numpy()]


def coordinates(d, mats, reference, seed):
    """Learn each assay projection in reference families only. One unit is one matching caliper."""
    missing=pd.Series(np.nan,index=d.index)
    age=pd.to_numeric(d.get("age",missing),errors="coerce").to_numpy(float)
    sex=d.get("sex",missing).astype("string").str.strip().str.lower().map(SEX).to_numpy(dtype=float,na_value=np.nan)
    clinical=[age,sex]
    axes=[age/5.]
    details=["exact sex","age within 5 years"]
    valid=np.isfinite(age) & np.isfinite(sex)
    for col in ROUTINE:
        if col not in d: continue
        v=pd.to_numeric(d[col],errors="coerce").to_numpy(float)
        observed=v[reference]; observed=observed[np.isfinite(observed)]
        if len(observed)<3 or observed.std()<1e-8: continue
        clinical.append(v)
        axes.append((v-observed.mean())/(1.5*observed.std()))
        valid &= np.isfinite(v)
        details.append(col)
    if len(clinical)<4:valid[:]=False
    clin_x,_=frozen_scale(np.column_stack(clinical),reference)
    if not clin_x.shape[1]:clin_x=np.zeros((len(d),1))
    standardized=[]
    for layer,raw in mats.items():
        # At least half of each requested assay must actually be measured for a twin match.
        nobs=np.zeros(len(d),int)
        for start in range(0,raw.shape[1],32):
            nobs+=np.isfinite(raw[:,start:start+32]).sum(axis=1)
        valid &= nobs >= raw.shape[1]/2
        x,_=frozen_scale(raw,reference)
        rank=min(3,x.shape[1],len(reference)-1)
        if rank<1: return None
        pca=PCA(n_components=rank,svd_solver="randomized",random_state=seed)
        pca.fit(x[reference])
        useful=pca.explained_variance_>1.e-8
        if not useful.any(): return None
        z=pca.transform(x)[:,useful]/np.sqrt(pca.explained_variance_[useful])
        axes.extend((z/1.5).T)
        details.append(layer+"_PCs="+str(int(useful.sum())))
        standardized.append(x)
    return np.column_stack(axes),sex,valid,clin_x,np.column_stack([clin_x,*standardized]),details


def twin_counts(d,geo,sex,valid,ref,cores):
    count=np.zeros(len(d),int)
    for code in (0,1):
        bank=ref[valid[ref] & (sex[ref]==code)]
        query=np.flatnonzero(valid & (sex==code))
        if len(bank) and len(query):
            tree=cKDTree(geo[bank])
            count[query]=tree.query_ball_point(geo[query],1.,p=np.inf,workers=cores,return_length=True)
    # Reference relatives can see their own family's representative. Remove that match.
    family_bank=pd.Series(ref,index=d.group.iloc[ref])
    own=d.group.map(family_bank)
    loc=np.flatnonzero(own.notna() & valid)
    donor=own.iloc[loc].to_numpy(int)
    same=valid[donor] & (sex[loc]==sex[donor]) & (np.max(np.abs(geo[loc]-geo[donor]),axis=1)<=1.)
    count[loc[same]]-=1
    return count


def neighbor_risk(d,geo,sex,valid,train,query,y,w,background,cores,residual=False):
    # Restrict outcome donors to one per family. Match in baseline coordinates only.
    donor=representatives(d,train & valid)
    pred=background.copy()
    for code in (0,1):
        bank=donor[sex[donor]==code]
        ix=np.flatnonzero(valid[query] & (sex[query]==code))
        if not len(bank) or not len(ix): continue
        tree=cKDTree(geo[bank]); k=min(32,len(bank))
        for start in range(0,len(ix),4096):
            sel=ix[start:start+4096]
            dist,idx=tree.query(geo[query[sel]],k=list(range(1,k+1)),p=np.inf,
                                distance_upper_bound=np.nextafter(1.,np.inf),workers=cores)
            ok=np.isfinite(dist); matched=bank[np.minimum(idx,len(bank)-1)]
            wm=w[matched]*ok
            mass=wm.sum(axis=1)
            local=np.divide((wm*y[matched]).sum(axis=1),mass,out=np.zeros(len(sel)),where=mass>0)
            eff=np.divide(mass**2,(wm**2).sum(axis=1),out=np.zeros(len(sel)),where=mass>0)
            if residual:
                pred[sel]=np.clip(background[sel]+eff*local/(eff+10),1.e-8,1-1.e-8)
            else:
                pred[sel]=(eff*local+10*background[sel])/(eff+10)
    return pred


def study(d,mats,names,folder,seed,cores,reps,trait=None):
    d=d.reset_index(drop=True)
    dev=d.role.eq("training").to_numpy(); test=d.role.eq("test").to_numpy()
    u=d.group.map({g:stable_hash(g,seed+2027)%100 for g in d.group.unique()}).to_numpy()
    ref=representatives(d,dev & (u<20))
    gate=dev & (u>=20) & (u<35)
    fit=dev & (u>=35)
    broad,clean,flags,groups=timing(d)
    rows={key:[] for key in SHEETS}
    profiles_out={}
    def status(part,reason,**kw):rows["study_status"].append(dict(analysis=part,status=reason,**kw))
    def finish():
        return {**{k:pd.DataFrame(v) for k,v in rows.items() if v},**profiles_out}
    for split,mask in (("train",dev),("test",test)):
        for flag,values in {"broad_eligible":broad,"clean_eligible":clean,**flags}.items():
            rows["exclusions"].append(dict(split=split,reason=flag,N=int((mask & values).sum()),
                     among_broad=int((mask & broad & values).sum()),denominator=int(mask.sum())))
    if len(mats)==1:
        layer=next(iter(mats))
        profiles_out,msg=profiles(d,mats[layer],names[layer],dev,test,folder/"scores.tsv.gz",groups=groups)
        status("profile",msg)
        from c1_four_groups import four_groups
        profiles_out.update(four_groups(d,mats[layer],names[layer],dev,test,folder,reps,seed))
        if layer=="met":
            from cad_timing_evidence import evidence
            profiles_out.update(evidence(d,mats[layer],names[layer],dev,test,
                                folder/"original_scores.tsv.gz",reps,seed))
            from cad_icd_first import analyze_icd_first
            profiles_out.update(analyze_icd_first(d,mats[layer],dev,test,folder,reps,seed))
    from c1_personal_predictivity import personal_predictivity, unavailable_predictivity
    if len(ref)<3 or not broad.any():
        profiles_out.update(unavailable_predictivity(d,broad,folder,"insufficient_reference_or_risk_set"))
        status("risk","insufficient_partition"); return finish()
    coord=coordinates(d,mats,ref,seed)
    if coord is None:
        profiles_out.update(unavailable_predictivity(d,broad,folder,"insufficient_baseline_dimensions"))
        status("risk","insufficient_baseline_dimensions"); return finish()
    geo,sex,valid,clin_x,full_x,details=coord
    twin_available=bool(valid[ref].any())
    if not twin_available:status("twins","no_reference_family_with_complete_matching_variables")
    count=twin_counts(d,geo,sex,valid,ref,cores)
    profiles_out.update(personal_predictivity(d,full_x,geo,sex,valid,count,broad,folder,cores,seed,reps))
    if not (gate & broad).any() or not (fit & broad).any() or not (test & broad).any():
        status("risk","insufficient_partition"); return finish()
    # Threshold is frozen before looking at risk outcomes. Ties are retained, zero never qualifies.
    # Use the entire baseline gate, so the rule never requires a future timing classification.
    cut=max(1,int(np.ceil(np.median(count[gate]))))
    high=count>=cut
    members=d[["eid","role"]].assign(reference_fold=np.select([dev & (u<20),gate,fit],
                   ["reference","gate","fit"],default="test"),twin_families=count,
                   geometry_complete=valid,high_support=high,rule_threshold=cut)
    members.to_csv(folder/"twins.tsv.gz",sep="\t",index=False,compression="gzip")
    rows["twin_geometry"].append(dict(variables="; ".join(details),reference_families=len(ref),
        valid_reference_families=int(valid[ref].sum()),gate_N=int(gate.sum()),threshold=cut,
        gate_coverage=float(high[gate].mean()),max_risk_donors=32,prior_strength=10,
        caliper_in_reference_SD=1.5,assay_observed_fraction=.5))
    time=d.icd_t2e.to_numpy(float)
    event=d.icd_state.eq("incident").to_numpy(int)
    y=np.zeros(len(d),int); w=np.zeros(len(d))
    try:
        km=CensoringKM().fit(time[fit & broad],event[fit & broad],HORIZON,.01)
        y[broad],w[broad]=km.labels_weights(time[broad],event[broad])
    except ValueError as exc:
        if "Insufficient training censoring support" not in str(exc):raise
        status("risk",str(exc));return finish()
    for tag,mask in (("fit",fit),("gate",gate),("test",test)):
        for risk,cohort in (("broad",broad),("clean",clean)):
            for name,level in (("all",np.ones(len(d),bool)),("high",high),("low",~high)):
                m=mask & cohort & level
                rows["twin_counts"].append(dict(split=tag,cohort=risk,support=name,N=int(m.sum()),
                    denominator=int((mask & cohort).sum()),geometry_complete=int((m & valid).sum()),
                    known_10y_N=int((m & (w>0)).sum()),cases_10y=int((m & (w>0) & (y==1)).sum()),
                    p25=np.quantile(count[m],.25) if m.any() else np.nan,
                    median=np.median(count[m]) if m.any() else np.nan,p75=np.quantile(count[m],.75) if m.any() else np.nan))
    query=np.flatnonzero(test & broad & (w>0))
    preds={}; train_ns={}; residual_ns={}
    def enough(m):return min((y[m]==0).sum(),(y[m]==1).sum())>=5
    for arm,include in (("broad",broad),("clean",clean),("clean_high",clean & high)):
        use=fit & include & (w>0)
        if not enough(use) or not len(query):
            status("risk_"+arm,"insufficient_known_cases_or_controls",train_N=int(use.sum()));continue
        train_ns[arm]=int(use.sum())
        base=fitted(clin_x,use,y,w,query)
        preds[arm,"clinical_ridge"]=base
        donor_mask=gate & include & (w>0) & valid
        donor_rows=np.flatnonzero(donor_mask)
        out_sample=fitted(full_x,use,y,w,np.r_[query,donor_rows])
        molecular=out_sample[:len(query)]
        preds[arm,"clinical_omics_ridge"]=molecular
        # Gate outcomes are opened only after the support cut and the ridge model are frozen.
        # Estimate the molecular model's out-of-family errors, then borrow those errors locally.
        residual=np.zeros(len(d))
        residual[donor_rows]=y[donor_rows]-out_sample[len(query):]
        residual_ns[arm]=len(representatives(d,donor_mask))
        if not twin_available:continue
        preds[arm,"omics_ridge_residual_twins"]=neighbor_risk(
            d,geo,sex,valid,donor_mask,query,residual,w,molecular,cores,residual=True)
        preds[arm,"twins"]=neighbor_risk(d,geo,sex,valid,use,query,y,w,base,cores)
    scopes={"broad_all":np.ones(len(query),bool),"clean_all":clean[query],
            "clean_high":(clean & high)[query],"clean_low":(clean & ~high)[query]}
    for scope,m in scopes.items():
        q=query[m]
        if not enough(q):
            status(scope,"insufficient_test_cases_or_controls",N=len(q)); continue
        for (arm,model),scores in preds.items():
            p=scores[m]; val=metrics(y[q],w[q],p,*auc_layout(p))
            rows["fit_test"].append(dict(test_set=scope,train_set=arm,model=model,train_N=train_ns[arm],
                residual_donor_families=residual_ns[arm] if model=="omics_ridge_residual_twins" else 0,
                test_known_N=len(q),test_cases=int(y[q].sum()),**dict(zip(METRICS,val))))
        def pair(label,a,b):
            if a not in preds or b not in preds:return
            results=compare(y[q],w[q],d.group.to_numpy()[q],preds[a][m],preds[b][m],reps,label,
                            "/".join(a),"/".join(b))
            rows["contrasts"].extend(dict(item,test_set=scope) for item in results)
        for model in ("clinical_ridge","clinical_omics_ridge","twins","omics_ridge_residual_twins"):
            pair("H1_remove_ambiguous",("broad",model),("clean",model))
            if scope.startswith("clean"):
                pair("H2_remove_low_support_train",("clean",model),("clean_high",model))
        if scope.startswith("clean"):
            pair("strong_model_neighbor_increment",("clean","clinical_omics_ridge"),("clean","omics_ridge_residual_twins"))
            pair("ABM_vs_molecular_ridge",("clean","clinical_omics_ridge"),("clean","twins"))
    if all(("clean",m) in preds for m in ("clinical_ridge","clinical_omics_ridge","twins","omics_ridge_residual_twins")):
        for target,pop in (("all_clean",clean),("complete_twin_coordinates",clean & valid)):
            results=risk_standardized(d,query,pop,high,y,w,preds,reps,seed)
            rows["selection"].extend(dict(row,eligible_set=target) for row in results)
    if ("clean","twins") in preds:
        pd.DataFrame(dict(eid=d.eid.to_numpy()[query],group=d.group.to_numpy()[query],y=y[query],w=w[query],
                          clean=clean[query],high=high[query],twins=preds["clean","twins"])).to_csv(
                folder/"risk.tsv.gz",sep="\t",compression="gzip",index=False,float_format="%.17g")
    status("risk","evaluated" if rows["fit_test"] else "no_evaluable_test_set",horizon_years=HORIZON,
           bootstrap_draws=reps)
    return finish()


def risk_standardized(d,q,clean,high,y,w,preds,reps,seed):
    """Compare high vs low on the fixed clean-trained models within baseline risk deciles."""
    keep=clean[q]
    if not keep.any():return []
    risk=preds["clean","clinical_ridge"][keep]
    edges=np.unique(np.quantile(risk,np.arange(1,10)/10))
    bins=np.searchsorted(edges,risk,side="right")
    ix=q[keep]; nbin=len(edges)+1
    h=high[ix]; weight=w[ix]
    mass_h=np.bincount(bins[h],weights=weight[h],minlength=nbin)
    mass_l=np.bincount(bins[~h],weights=weight[~h],minlength=nbin)
    overlap=(mass_h>0)&(mass_l>0)
    if not overlap.any():return []
    # Target the pooled distribution on common support. Keep this target fixed in resamples.
    target=(mass_h+mass_l)[overlap]; target/=target.sum()
    level=bins+nbin*h.astype(int)
    _,family=np.unique(d.group.to_numpy()[ix],return_inverse=True)
    out=[]
    for label in ("twins","clinical_omics_ridge","omics_ridge_residual_twins"):
        err=(preds["clean",label][keep]-y[ix])**2
        def delta(weights):
            total=np.bincount(level,weights=weights,minlength=2*nbin).reshape(2,nbin)[:,overlap]
            if (total<=0).any():return np.nan
            sums=np.bincount(level,weights=weights*err,minlength=2*nbin).reshape(2,nbin)[:,overlap]
            return np.diff((sums/total)@target)[0]
        blocks=family.max()+1; rng=np.random.default_rng(seed); boot=[]
        for _ in range(reps):
            freq=np.bincount(rng.integers(blocks,size=blocks),minlength=blocks)
            boot.append(delta(weight*freq[family]))
        vals=np.asarray(boot); vals=vals[np.isfinite(vals)]
        lo,hi=np.quantile(vals,[.025,.975]) if len(vals)>=max(20,.8*reps) else (np.nan,np.nan)
        ok=overlap[bins]
        out.append(dict(model=label,contrast="high_minus_low_risk_standardized_Brier",difference=delta(weight),
               lower=lo,upper=hi,valid_bootstraps=len(vals),overlap_strata=int(overlap.sum()),
               high_known_N=int((h & ok).sum()),low_known_N=int((~h & ok).sum()),
               high_cases=int(y[ix][h & ok].sum()),low_cases=int(y[ix][~h & ok].sum()),
               clean_known_N=len(ix),standardization="pooled clinical-risk deciles on shared support"))
    return out


def load_roles(a):
    phe=Path(a.ygb_phe)
    group_file=a.group_file or (None if a.group_col else str(phe/"rap/raw/ygb_kinship_groups.tsv"))
    ensure_roster(ygb_phe=phe,output=a.roster,seed=a.seed,group_file=group_file or "",
                  group_col=a.group_col or "",r_bin=a.r_bin)
    rec=roster_record(a.roster)
    roles=read_table(a.roster)[["eid","role"]]
    if rec.get("group_file"):
        families=_read_group_file(rec["group_file"],a.r_bin)[["eid","group"]]
    else:
        allrds=phe/"Rdata/all.rds"
        families=_cohort_projection(allrds,sha(allrds),a.r_bin,rec["group_col"])
    roles=roles.merge(families,on="eid",validate="one_to_one",how="left")
    if (roles.group.isna().any() or not roles.role.isin(["training","test"]).all()
            or roles.groupby("group").role.nunique().max()>1):
        raise ValueError("Invalid family assignments for the frozen C1 split")
    return roles,rec


def load_assay(phe,layer,folder,r_bin):
    folder.mkdir()
    subprocess.run([r_bin,str(HERE/"timing_matrix.R"),str(phe/"Rdata"/(layer+".rds")),str(folder)],check=True)
    ids=pd.read_csv(folder/"ids.txt",header=None,dtype=str).iloc[:,0]
    names=(folder/"features.txt").read_text().splitlines()
    raw=np.memmap(folder/"baseline.f32",dtype="<f4",mode="r",shape=(len(ids),len(names)),order="F")
    observed=np.zeros(len(ids),bool)
    for j in range(0,len(names),32):observed|=np.isfinite(raw[:,j:j+32]).any(axis=1)
    index=pd.Series(np.arange(len(ids))[observed],index=ids[observed],name=layer)
    if not index.index.is_unique:raise ValueError("Duplicate assay IDs")
    return index,names,raw


def run(a):
    base=Path(a.trait_dir); dest=base/"main/results/c1_timing_twins"
    tmp=base/"tmp/c1_timing_twins"; qc=base/"qc/c1_timing_twins"
    for p in (dest,tmp,qc):p.mkdir(parents=True,exist_ok=True)
    meta=json.loads(Path(a.snapshot).with_name("manifest.json").read_text())
    hashes={r["dx"]:r["endpoint_sha256"] for r in meta["counts"]}
    layers=list(dict.fromkeys(a.layers.split(",")))
    if not set(layers)<={"prot","met"}:raise ValueError("Unknown molecular layer")
    views=[*layers,*(["prot_met"] if len(layers)==2 else [])]
    needed={"eid","age","sex","drug.lipid","baseline_lipid_i0","icd_days","ref_days","icd_state",
            "ref_state","icd_t2e","any1_state","any_early_days","any_first_days","any_followup_days",
            "any_prior","any_undated","source_first_span_days","any_near_before_days",*ROUTINE}
    with file_lock(qc/"study.lock"),tempfile.TemporaryDirectory(dir=tmp) as td,threadpool_limits(limits=a.cores):
        work=Path(td); assays={}; roles,rec=load_roles(a)
        d=pd.read_csv(a.snapshot,sep="\t",usecols=lambda c:c in needed or c.startswith("src_") or c=="any_early_sources",dtype={"eid":str})
        required_sources={"src_icd10_days","src_icd10_undated","src_gp_days","src_srd_days"}
        if required_sources.difference(d.columns):
            raise ValueError("The AOD snapshot predates the four-group study. Run le8.sh c1_correlate to regenerate it.")
        d=d.merge(roles,on="eid",validate="one_to_one")
        code={f:sha(HERE/f) for f in ("c1_timing_twins.py","timing_study.py","timing_matrix.R",
                                        "cad_timing_evidence.py","cad_icd_first.py","abm.py","c1_protocol.py",
                                        "c1_four_groups.py","c1_personal_predictivity.py")}
        for view in views:
            parts=view.split("_"); folder=work/("study_"+view);folder.mkdir()
            assay_files={}
            for k in parts:
                source=Path(a.ygb_phe)/"Rdata"/(k+".rds"); stat=source.stat()
                assay_files[k]=[str(source.resolve()),stat.st_size,stat.st_mtime_ns]
            spec=[base.name,view,meta["endpoint_sha256"],rec["roster_sha256"],a.seed,a.bootstrap,code,assay_files,
                  os.getenv("LE8_TREATMENT_BASELINE_CONFIRMED","")]
            signature=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
            if a.replace=="FALSE" and (dest/(view+".xlsx")).is_file():
                last=pd.read_excel(dest/(view+".xlsx"),sheet_name="provenance").iloc[0]
                if last.get("fit_signature")==signature:
                    private=json.loads(last["private_hashes"])
                    if all((tmp/name).is_file() and sha(tmp/name)==value for name,value in private.items()):
                        print(f"[C1] reuse {base.name} {view}",flush=True);continue
            for k in parts:
                if k not in assays:assays[k]=load_assay(Path(a.ygb_phe),k,work/k,a.r_bin)
            frame=d
            for k in parts:frame=frame.merge(assays[k][0],left_on="eid",right_index=True,validate="one_to_one")
            frame=frame.reset_index(drop=True)
            print(f"[C1] {base.name} {view} assay/phenotype/family overlap N={len(frame)}",flush=True)
            if len(frame) and not frame.icd_state.eq("unavailable").all():
                mats={k:assays[k][2][frame[k].to_numpy()] for k in parts}
                tables=study(frame,mats,{k:assays[k][1] for k in parts},folder,a.seed,a.cores,a.bootstrap,trait=base.name)
                del mats
            else:tables={"study_status":pd.DataFrame([dict(analysis="study",status="no_overlap_or_ICD_unavailable")])}
            local=[("scores.tsv.gz",view+"_profile_scores_LOCAL_ONLY.tsv.gz"),
                   ("four_scores.tsv.gz",view+"_four_group_scores_LOCAL_ONLY.tsv.gz"),
                   ("predictivity.tsv.gz",view+"_predictivity_LOCAL_ONLY.tsv.gz"),
                   ("twins.tsv.gz",view+"_twins_LOCAL_ONLY.tsv.gz"),("risk.tsv.gz",view+"_risk_LOCAL_ONLY.tsv.gz")]
            if view=="met":
                local.extend([("original_scores.tsv.gz","met_original_profiles_LOCAL_ONLY.tsv.gz"),
                              ("icd_first.tsv.gz","met_icd_first_LOCAL_ONLY.tsv.gz")])
            private={target:sha(folder/source) for source,target in local if (folder/source).is_file()}
            tables["provenance"]=pd.DataFrame([dict(trait=base.name,layer=view,dx="icd",seed=a.seed,
                    snapshot_sha256=meta["endpoint_sha256"],roster_sha256=rec["roster_sha256"],icd_sha256=hashes["icd"],
                    any1_sha256=hashes["any1"],cohort_N=len(d),assay_overlap_N=len(frame),protocol="timing_twins_v1",
                    fit_signature=signature,private_hashes=json.dumps(private,sort_keys=True))])
            staged=folder/"study.xlsx"
            with pd.ExcelWriter(staged,engine="openpyxl") as book:
                for key,tbl in tables.items():tbl.to_excel(book,sheet_name=key,index=False,freeze_panes=(1,0))
            for source,target in local:
                if target in private:
                    (folder/source).chmod(0o600);os.replace(folder/source,tmp/target)
                else:(tmp/target).unlink(missing_ok=True)
            os.replace(staged,dest/(view+".xlsx"))
            print(f"[C1] published {dest/(view+'.xlsx')}",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("trait-dir","snapshot","ygb-phe","roster","layers","r-bin"):p.add_argument("--"+key,required=True)
    p.add_argument("--group-file");p.add_argument("--group-col")
    p.add_argument("--replace",choices=["TRUE","FALSE"],default="TRUE")
    p.add_argument("--cores",type=int,default=1);p.add_argument("--seed",type=int,default=2026)
    p.add_argument("--bootstrap",type=int,default=300)
    args=p.parse_args()
    if args.bootstrap<20 or args.cores<1:p.error("Need >=20 family draws and >=1 CPU")
    run(args)
