"""One development-frozen selective-prediction rule evaluated across every model.

The selector targets excess Brier error of the declared glmnet anchor. It is
model-anchored, not an intrinsic probability that a person's future is knowable.
Gate outcomes fit error heads; independent calibration features freeze coverage;
neither audit nor test outcomes are accepted by the fitting API.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.special import logit
from sklearn.ensemble import HistGradientBoostingRegressor

import abm

VERSION = "LE8-C1-common-selector-2026-11-16-v1"
PRIMARY = "P1_incident_support"
RISK_ONLY = "P0_risk_only"
CONTROLS = (RISK_ONLY, "raw_count_support", "low_risk", "confidence", "random")
COVERAGES = (.1, .2, .3, .4, .5, .6, .8, 1.)
RISK_BINS = 10


def _area(y,x):
    return (np.trapezoid if hasattr(np,"trapezoid") else np.trapz)(y,x)


def _align(ids, pred):
    ids=np.asarray(ids).astype(str);other=np.asarray(pred["eid"]).astype(str)
    if len(set(other)) != len(other) or set(ids) != set(other):
        raise ValueError("A selector projection differs from the frozen shared roster")
    index={key:i for i,key in enumerate(other)};order=np.asarray([index[key] for key in ids])
    out={}
    for key,value in pred.items():
        if key=="diagnostics":
            out[key]={k:np.asarray(v)[order] for k,v in value.items() if np.asarray(v).shape==other.shape}
        elif isinstance(value,np.ndarray) and value.shape==other.shape:out[key]=value[order]
        else:out[key]=value
    return out


def risk_features(p):
    p=np.clip(np.asarray(p,float),1e-7,1-1e-7)
    return np.c_[p,p*(1-p),logit(p)]


def molecular_features(p, projection):
    """Build-donor evidence, never query outcomes or query-batch summaries."""
    n=len(p);d=projection.get("diagnostics",{})
    def get(key,default=0.):
        x=np.asarray(d.get(key,np.full(n,default)),float)
        if x.shape!=(n,):raise ValueError("Molecular selector feature is misaligned: "+key)
        return np.nan_to_num(x,nan=default,posinf=default,neginf=default)
    return np.c_[risk_features(p),np.log1p(np.maximum(get("eligible_match_count"),0)),
        np.log1p(np.maximum(get("family_ess"),0)),get("twin_fraction"),
        np.minimum(get("radius_ratio",100),100),np.minimum(get("nearest_radius_ratio",100),100),
        get("supported"),get("local_event_rate")-p,
        get("local_clinical_residual"),get("local_clinical_oof_brier_excess"),
        get("local_clinical_oof_logloss_excess")]


@dataclass
class ErrorHead:
    name:str
    seed:int=2026
    bins:int=10
    min_events:int=20
    iterations:int=120
    residual_model:object=None
    error_model:object=None
    status:str="not_fitted"
    risk_edges:object=None

    def fit(self,X,p,y,w,qc):
        X,p,y,w,qc=map(np.asarray,(X,p,y,w,qc))
        use=qc.astype(bool)&(w>0)&np.isfinite(p)&np.all(np.isfinite(X),axis=1)
        if (len(X)!=len(p) or y.shape!=p.shape or w.shape!=p.shape or
                not np.isin(y[use],[0,1]).all() or not np.isfinite(w).all() or np.any(w<0)):
            raise ValueError("Invalid common-selector gate arrays")
        cases=int(np.sum(y[use]));controls=int(np.sum(use))-cases
        if min(cases,controls)<self.min_events or use.sum()<100:
            self.status=f"insufficient_tune_gate_information:N={use.sum()};cases={cases};controls={controls}"
            return self
        self.risk_edges=np.unique(np.quantile(p[use],np.linspace(0,1,self.bins+1)))[1:-1]
        labels=np.searchsorted(self.risk_edges,p,side="right")
        loss=(y-p)**2;residual=loss.copy()
        for label in np.unique(labels[use]):
            ix=use&(labels==label)
            residual[ix]-=np.average(loss[ix],weights=w[ix])
        params=dict(max_iter=self.iterations,max_leaf_nodes=3,learning_rate=.04,
            min_samples_leaf=max(30,min(200,int(use.sum()/20))),l2_regularization=20,
            early_stopping=False,random_state=self.seed)
        self.residual_model=HistGradientBoostingRegressor(**params).fit(X[use],residual[use],sample_weight=w[use]/w[use].mean())
        self.error_model=HistGradientBoostingRegressor(**params).fit(X[use],loss[use],sample_weight=w[use]/w[use].mean())
        self.status=f"fitted_tune_gate:N={use.sum()};cases={cases};controls={controls}"
        return self

    def score(self,X):
        if self.residual_model is None:return np.full(len(X),np.nan),np.full(len(X),np.nan)
        return -self.residual_model.predict(X),np.clip(self.error_model.predict(X),0,1)


def _control_scores(p,raw_count,ids,seed):
    d=raw_count.get("diagnostics",{}) if raw_count is not None else {}
    n=len(p)
    count=np.asarray(d.get("eligible_match_count",np.full(n,np.nan)),float)
    ess=np.asarray(d.get("family_ess",np.full(n,np.nan)),float)
    return dict(raw_count_support=np.log1p(count)+np.log1p(ess),low_risk=-np.asarray(p),
        confidence=np.abs(np.asarray(p)-.5),
        random=np.asarray([abm.tie_hash(i,seed+808)/2**64 for i in ids]))


def fit_common_selector(gate,calibration,gate_y,gate_w,*,gate_molecular,calibration_molecular,
                        gate_count=None,calibration_count=None,seed=2026,bins=10,
                        min_events=20,iterations=120,coverages=COVERAGES):
    """Fit without an audit/test truth parameter; selection calibration is feature-only."""
    p=np.asarray(gate["probability"],float);cp=np.asarray(calibration["probability"],float)
    qc=np.asarray(gate["technical_qc"],bool);cq=np.asarray(calibration["technical_qc"],bool)
    names=(RISK_ONLY,PRIMARY)
    gx={RISK_ONLY:risk_features(p),PRIMARY:molecular_features(p,gate_molecular)}
    cx={RISK_ONLY:risk_features(cp),PRIMARY:molecular_features(cp,calibration_molecular)}
    heads={name:ErrorHead(name,seed+101*j,bins,min_events,iterations).fit(gx[name],p,gate_y,gate_w,qc)
           for j,name in enumerate(names)}
    scores=_control_scores(cp,calibration_count,calibration["eid"],seed)
    expected_edges={}
    for name,head in heads.items():
        scores[name],error=head.score(cx[name])
        valid=cq&np.isfinite(error)
        expected_edges[name]=np.unique(np.quantile(error[valid],np.linspace(0,1,11)))[1:-1] if valid.any() else np.array([])
    rules={name:{float(q):abm.FrozenRule().fit(score,cq&np.isfinite(score),
        np.asarray(calibration["eid"]).astype(str),float(q),seed+1601) for q in coverages}
        for name,score in scores.items() if np.isfinite(score).any()}
    edges={k:np.unique(np.quantile(cp[cq&np.isfinite(cp)],np.linspace(0,1,k+1)))[1:-1]
           for k in (10,50)}
    return dict(version=VERSION,heads=heads,rules=rules,risk_edges=edges,expected_error_edges=expected_edges,seed=seed,
        coverages=tuple(map(float,coverages)),primary_selector=PRIMARY,primary_coverage=.5,
        fit_role="tune_gate",threshold_role="calibration_fit features only",
        gate_id_hash=abm.fingerprint(gate["eid"]),calibration_id_hash=abm.fingerprint(calibration["eid"]),
        selector_target="risk-residual Brier error of frozen anchor model",
        intrinsic_individual_correctness_probability=False)


def project_common_selector(fitted,anchor,molecular,raw_count=None):
    p=np.asarray(anchor["probability"],float);qc=np.asarray(anchor["technical_qc"],bool)
    ids=np.asarray(anchor["eid"]).astype(str)
    scores=_control_scores(p,raw_count,ids,fitted["seed"]);expected={}
    for name,head in fitted["heads"].items():
        X=risk_features(p) if name==RISK_ONLY else molecular_features(p,molecular)
        scores[name],expected[name]=head.score(X)
    result={}
    for name,rules in fitted["rules"].items():
        score=scores[name];valid=qc&np.isfinite(score)
        result[name]=dict(score=score,expected_error=expected.get(name,np.full(len(p),np.nan)),
            masks={q:rule.apply(score,valid,ids) for q,rule in rules.items()},valid=valid)
    return result


def freeze_common_selector(c,method_paths,out,protocol_path=None):
    """Read only development labels and frozen feature/prediction exports."""
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    paths={k:Path(v) for k,v in method_paths.items()}
    anchor="glmnet" if "glmnet" in paths else "clinical" if "clinical" in paths else next(iter(paths))
    molecular="abm_clinical" if "abm_clinical" in paths else "abm" if "abm" in paths else anchor
    protocol=json.loads(Path(protocol_path).read_text()) if protocol_path else {}
    declared_anchor=protocol.get("primary_risk_model","glmnet")
    if declared_anchor!="glmnet":
        raise ValueError("The current CVD protocol declares glmnet as its primary selector target")
    bins=int(protocol.get("risk_bins",10));coverages=tuple(protocol.get("coverage_grid",COVERAGES))
    if bins!=10 or float(protocol.get("primary_coverage",.5))!=.5:
        raise ValueError("The CVD shared-selector protocol fixes10 primary bins and50% primary coverage")
    loaded={}
    for role in ("tune_gate","calibration_fit","calibration_audit","test"):
        ids=c.ids[c.parts[role]]
        current={}
        for method in dict.fromkeys((anchor,molecular,"abm")):
            if method not in paths:continue
            file=paths[method]/f"{role}_predictions_LOCAL_ONLY.joblib"
            if not file.is_file():raise ValueError(f"{method}: refit required to export independent {role} predictions")
            pred=_align(ids,joblib.load(file))
            if not np.array_equal(np.asarray(pred["group"]).astype(str),c.groups[c.parts[role]]):
                raise ValueError("Selector and shared family assignments differ")
            current[method]=pred
        loaded[role]=current
    g,cal=loaded["tune_gate"],loaded["calibration_fit"];ix=c.parts["tune_gate"]
    fitted=fit_common_selector(g[anchor],cal[anchor],c.y[ix],c.w[ix],gate_molecular=g[molecular],
        calibration_molecular=cal[molecular],gate_count=g.get("abm"),calibration_count=cal.get("abm"),
        seed=c.cfg.seed,bins=bins,min_events=c.cfg.min_events,iterations=c.cfg.gate_iterations,coverages=coverages)
    fitted.update(anchor_model=anchor,molecular_evidence_model=molecular,protocol=protocol)
    projected={role:project_common_selector(fitted,items[anchor],items[molecular],items.get("abm"))
        for role,items in loaded.items() if role in ("calibration_audit","test")}
    bundle=dict(fitted=fitted,projected=projected,eid={role:c.ids[c.parts[role]] for role in projected},
        group={role:c.groups[c.parts[role]] for role in projected},
        anchor_probability={role:loaded[role][anchor]["probability"] for role in projected},
        technical_qc={role:loaded[role][anchor]["technical_qc"] for role in projected})
    required_evidence=("local_clinical_residual","local_clinical_oof_brier_excess","local_clinical_oof_logloss_excess")
    missing_evidence=[key for key in required_evidence if key not in loaded["test"][molecular].get("diagnostics",{})]
    primary_status=("unavailable_primary_glmnet_secondary_"+anchor if anchor!=declared_anchor else
        "completed_reduced_evidence_secondary" if missing_evidence else "completed_primary_P1_vs_P0")
    manifest=dict(version=VERSION,anchor_model=anchor,declared_anchor_model=declared_anchor,
        primary_anchor_changed=anchor!=declared_anchor,molecular_evidence_model=molecular,
        selector_target=fitted["selector_target"],selector_fit_role="tune_gate",
        selector_threshold_role="calibration_fit features only",audit_outcomes_used=False,test_outcomes_used=False,
        status=primary_status if PRIMARY in projected["test"] else "insufficient_primary_reliability_support",
        missing_primary_evidence_fields=missing_evidence,
        primary_risk_bins=10,sensitivity_risk_bins=50,primary_coverage=.5,
        primary_selector=PRIMARY,control_selectors=list(CONTROLS),
        same_frozen_masks_for_all_models=True,model_neutral_selector=False,
        test_provenance="previously_inspected_or_unverified_exploratory",
        discovery_or_replication="CVD_exploratory_rerun",protocol=protocol,
        head_status={k:v.status for k,v in fitted["heads"].items()},
        selected_roster_hashes={name:{str(q):abm.fingerprint(bundle["eid"]["test"][mask])
            for q,mask in value["masks"].items()} for name,value in projected["test"].items()})
    joblib.dump(bundle,out/"common_selector_LOCAL_ONLY.joblib",compress=3)
    for role,current in projected.items():
        rows=dict(eid=bundle["eid"][role],group=bundle["group"][role],
            anchor_probability=bundle["anchor_probability"][role])
        for name,projection in current.items():
            rows[name+"_score"]=projection["score"]
            rows[name+"_expected_Brier_error"]=projection["expected_error"]
            for q,mask in projection["masks"].items():rows[f"{name}_selected_{q:g}"]=mask
        pd.DataFrame(rows).to_csv(out/f"{role}_common_selector_LOCAL_ONLY.tsv.gz",sep="\t",index=False)
    for path in out.glob("*LOCAL_ONLY*"):path.chmod(0o600)
    abm.json_write(out/"FROZEN_BEFORE_TEST.json",manifest)
    bundle["manifest"]=manifest
    return bundle


class FamilyCells:
    """Sparse family sums: bootstrap does not rescan all people for every stratum."""
    def __init__(self,groups,masks,w,losses,labels=None,bins=1):
        levels,inverse=np.unique(np.asarray(groups).astype(str),return_inverse=True)
        self.families=len(levels);self.names=list(masks);self.bins=bins
        labels=np.zeros(len(groups),int) if labels is None else np.asarray(labels,int)
        rr=[];cc=[];ww=[];people=[]
        for j,name in enumerate(self.names):
            use=np.flatnonzero(np.asarray(masks[name],bool)&(w>0))
            rr.append(j*bins+labels[use]);cc.append(inverse[use]);ww.append(w[use]);people.append(use)
        row=np.concatenate(rr);col=np.concatenate(cc);weight=np.concatenate(ww);person=np.concatenate(people)
        shape=(len(masks)*bins,len(levels))
        self.den=csr_matrix((weight,(row,col)),shape=shape)
        self.num={key:csr_matrix((weight*np.asarray(value)[person],(row,col)),shape=shape)
                  for key,value in losses.items()}

    def sums(self,mult):
        den=np.asarray(self.den@mult).reshape(len(self.names),self.bins)
        num={key:np.asarray(matrix@mult).reshape(den.shape) for key,matrix in self.num.items()}
        return den,num


def _interval(values,point,B,family_n=None):
    a=np.asarray(values,float);a=a[np.isfinite(a)];minimum=max(20,int(np.ceil(.8*B)))
    ok=len(a)>=minimum and np.isfinite(point) and (family_n is None or family_n>=2)
    return dict(lower=float(np.quantile(a,.025)) if ok else np.nan,
        upper=float(np.quantile(a,.975)) if ok else np.nan,bootstrap_valid=len(a),
        inference_status="evaluated" if ok else "insufficient_independent_families" if family_n is not None and family_n<2 else "insufficient_valid_family_bootstrap")


def _std_pair(den,num,all_index,a,b,fixed=None):
    common=(den[a]>0)&(den[b]>0)&(den[all_index]>0) if fixed is None else fixed
    if not np.any(common) or np.any(den[a,common]<=0) or np.any(den[b,common]<=0):
        return np.nan,np.nan,common,0.
    mass=den[all_index,common];total=den[all_index].sum()
    va=np.average(num[a,common]/den[a,common],weights=mass)
    vb=np.average(num[b,common]/den[b,common],weights=mass)
    return float(va),float(vb),common,float(mass.sum()/total) if total else 0.


def evaluate_common_selector(bundle,predictions,truth,km,*,bootstrap=1000,seed=2026):
    """All tests are conditional on the already saved common selector bundle."""
    from c1_evaluation import point_metrics, _risk_balance
    ids=np.asarray(truth["eid"]).astype(str)
    if not np.array_equal(ids,np.asarray(bundle["eid"]["test"]).astype(str)):
        raise ValueError("Common-selector test truth order differs from frozen IDs")
    aligned={name:_align(ids,pred) for name,pred in predictions.items()}
    groups=np.asarray(bundle["group"]["test"]).astype(str)
    for pred in aligned.values():
        if not np.array_equal(groups,np.asarray(pred["group"]).astype(str)):
            raise ValueError("Method families differ from the frozen selector")
    qc=np.asarray(bundle["technical_qc"]["test"],bool).copy()
    qc &= np.logical_and.reduce([np.isfinite(v["probability"]) for v in aligned.values()])
    y,w=km.labels_weights(truth["time"],truth["event"]);w=np.where(qc,w,0.)
    projected=bundle["projected"]["test"];anchor=bundle["fitted"]["anchor_model"]
    p=np.asarray(bundle["anchor_probability"]["test"],float)
    qgrid=bundle["fitted"]["coverages"];B=int(bootstrap)
    curves=[];pop=[];losses={}
    for model,pred in aligned.items():
        prob=np.asarray(pred["probability"],float)
        losses[(model,"Brier")]=(y-prob)**2
        losses[(model,"LogLoss")]=abm.logloss_rows(y,prob)
        pop.append(dict(model=model,population="same frozen test participants",**point_metrics(y,w,prob,qc)))
        for head,projection in projected.items():
            for q,mask in projection["masks"].items():
                mask=mask&qc
                curves.append(dict(risk_model=model,head=head,coverage_target=q,
                    actual_coverage=float(mask.sum()/qc.sum()) if qc.any() else np.nan,
                    primary_risk_model=model==anchor,same_people_across_models=True,
                    selected_roster_hash=abm.fingerprint(ids[mask]),**point_metrics(y,w,prob,mask)))
    status=[dict(head=name,status=head.status,selector_target_model=anchor,fit_role="tune_gate",
        threshold_role="calibration_fit features only",primary_status=bundle["manifest"]["status"])
        for name,head in bundle["fitted"]["heads"].items()]
    status.append(dict(head="P2_incident_plus_prevalent",status="unavailable_no_prevalent_selector_in_CVD_primary_protocol",
        historical_outcomes_mixed_with_incident=False))
    tables={"reliability_status":pd.DataFrame(status),"reliability_population_performance":pd.DataFrame(pop),
            "risk_coverage":pd.DataFrame(curves)}
    error_calibration=[]
    for head,edges in bundle["fitted"].get("expected_error_edges",{}).items():
        if head not in projected:continue
        error=projected[head]["expected_error"];labels=np.searchsorted(edges,error,side="right")
        for label in range(len(edges)+1):
            use=qc&np.isfinite(error)&(w>0)&(labels==label)
            if not use.any():continue
            expected=float(np.average(error[use],weights=w[use]))
            observed=float(np.average((y[use]-p[use])**2,weights=w[use]))
            error_calibration.append(dict(head=head,risk_model=anchor,bin=label,known_N=int(use.sum()),
                cases=int(y[use].sum()),mean_expected_Brier_error=expected,observed_Brier_error=observed,
                calibration_difference=observed-expected,bin_edges_source="calibration_fit_features_only",
                correctness_probability=False))
    tables["common_error_calibration"]=pd.DataFrame(error_calibration)
    if PRIMARY not in projected:
        return tables
    controls=[name for name in CONTROLS if name in projected]
    masks={name:projection["masks"][.5]&qc for name,projection in projected.items()}
    masks["all"]=qc
    rng=np.random.default_rng(seed+101)
    increments=[];matched_rows=[];plans=[]
    for bins,edges in bundle["fitted"]["risk_edges"].items():
        labels=np.searchsorted(edges,p,side="right");nb=len(edges)+1
        current=dict(masks)
        for name in (PRIMARY,"raw_count_support"):
            if name in masks:
                current[name+"_riskmatched"]=abm.risk_matched_mask(masks[name],qc,p,ids,edges,seed+412)
        plan=FamilyCells(groups,current,w,losses,labels,nb)
        den,num=plan.sums(np.ones(plan.families));index={name:j for j,name in enumerate(plan.names)}
        specs=[]
        pairs=[(PRIMARY,control,"head_increment") for control in controls]
        pairs += [(name,name+"_riskmatched","selection_vs_riskmatched")
                  for name in (PRIMARY,"raw_count_support") if name in masks]
        for model,reference,kind in pairs:
            a,b=index[model],index[reference]
            ma,mb,smd,ks=_risk_balance(p,current[model],current[reference],w)
            for (risk_model,metric),value in num.items():
                va,vb,common,mass=_std_pair(den,value,index["all"],a,b)
                row=dict(model=model,reference=reference,risk_model=risk_model,
                    primary_risk_model=risk_model==anchor,contrast=kind,metric=metric,coverage_target=.5,
                    risk_bins=int(bins),primary_risk_bins=bins==10,model_value=va,reference_value=vb,difference=va-vb,
                    bins_used=int(common.sum()),reference_weight_fraction=mass,
                    standardized_population="full risk distribution" if mass>=.999 else "common-support risk strata only",
                    full_risk_distribution_supported=mass>=.999,overlap_warning=mass<.8,
                    model_selected_N=int(current[model].sum()),reference_selected_N=int(current[reference].sum()),
                    model_cases=int(np.sum(y[current[model]])),reference_cases=int(np.sum(y[current[reference]])),
                    model_mean_risk=ma,reference_mean_risk=mb,risk_SMD=smd,risk_KS=ks,
                    selected_roster_hash=abm.fingerprint(ids[current[model]]),
                    reference_roster_hash=abm.fingerprint(ids[current[reference]]),
                    same_frozen_masks_across_models=True,paired_family_draws=True,bootstrap_unit="family",
                    standardization="same populated frozen development risk strata; shared all-person IPCW mass")
                target=increments if kind=="head_increment" else matched_rows
                target.append(row);specs.append((row,(risk_model,metric),a,b,common,[]))
        plans.append((plan,index,specs))
    # Identical family multiplicities for all methods, heads and both stratum grids.
    family_n=plans[0][0].families
    for _ in range(B):
        mult=np.bincount(rng.integers(0,family_n,size=family_n),minlength=family_n)
        for plan,index,specs in plans:
            den,num=plan.sums(mult)
            for row,key,a,b,common,draws in specs:
                va,vb,_,_=_std_pair(den,num[key],index["all"],a,b,common)
                draws.append(va-vb)
    for plan,index,specs in plans:
        for row,key,a,b,common,draws in specs:row.update(_interval(draws,row["difference"],B,plan.families))
    tables["reliability_increment_50"]=pd.DataFrame(increments)
    tables["common_riskmatched_selection"]=pd.DataFrame(matched_rows)
    # AURC is descriptive for every model; predeclared interval inference concerns
    # the glmnet anchor, while all-model paired intervals above use common masks.
    curve_masks={(name,q):proj["masks"][q]&qc for name,proj in projected.items() for q in qgrid}
    anchor_losses={k:v for k,v in losses.items() if k[0]==anchor}
    aurc_plan=FamilyCells(groups,curve_masks,w,anchor_losses)
    index={key:j for j,key in enumerate(aurc_plan.names)}
    def aurc_from_sums(den,num,name):
        ix=[index[(name,q)] for q in qgrid]
        values=np.divide(num[ix,0],den[ix,0],out=np.full(len(ix),np.nan),where=den[ix,0]>0)
        return float(_area(values,qgrid)/(qgrid[-1]-qgrid[0])) if np.isfinite(values).all() else np.nan
    aurc=[];spec=[];frame=pd.DataFrame(curves)
    for risk_model in aligned:
        for reference in controls:
            for metric in ("Brier","LogLoss"):
                values=[]
                for head in (PRIMARY,reference):
                    z=frame[frame.risk_model.eq(risk_model)&frame["head"].eq(head)].sort_values("coverage_target")
                    values.append(float(_area(z[metric],z.coverage_target)/(qgrid[-1]-qgrid[0])))
                row=dict(model=PRIMARY,reference=reference,risk_model=risk_model,primary_risk_model=risk_model==anchor,
                    metric="AURC_"+metric,difference=values[0]-values[1],model_aurc=values[0],reference_aurc=values[1],
                    coverage_axis="predeclared_target",coverage_min=qgrid[0],coverage_max=qgrid[-1],
                    lower=np.nan,upper=np.nan,bootstrap_valid=0,inference_status="descriptive_secondary_model",
                    bootstrap_unit="family",selector_target_model=anchor)
                aurc.append(row)
                if risk_model==anchor:spec.append((row,reference,metric,[]))
    rng=np.random.default_rng(seed+313)
    for _ in range(B):
        mult=np.bincount(rng.integers(0,aurc_plan.families,size=aurc_plan.families),minlength=aurc_plan.families)
        den,num=aurc_plan.sums(mult)
        for row,reference,metric,draws in spec:
            values=num[(anchor,metric)]
            draws.append(aurc_from_sums(den,values,PRIMARY)-aurc_from_sums(den,values,reference))
    for row,reference,metric,draws in spec:row.update(_interval(draws,row["difference"],B,aurc_plan.families))
    tables["aurc_contrasts"]=pd.DataFrame(aurc)
    tables["support_trend"]=pd.DataFrame(support_error_trends(
        aligned.get("abm"),groups,qc,w,p,y,bundle["fitted"]["risk_edges"][10],B,seed))
    return tables


def support_error_trends(raw_count,groups,qc,w,p,y,edges,B,seed):
    """Does raw twin support predict less anchor-model error within risk strata?"""
    diagnostics=raw_count.get("diagnostics",{}) if raw_count is not None else {}
    rows=[];plans=[];labels=np.searchsorted(edges,p,side="right");bins=len(edges)+1
    def slope(den,num,metric):
        mass=den[0];use=mass>0
        if not use.any():return np.nan
        wx=num["x"][0,use];wy=num[metric][0,use];ww=mass[use]
        numerator=np.sum(num["x_"+metric][0,use]-wx*wy/ww)
        denominator=np.sum(num["x2"][0,use]-wx**2/ww)
        return float(numerator/denominator) if denominator>1e-12 else np.nan
    for key in ("eligible_match_count","family_ess"):
        if key not in diagnostics:
            rows.extend(dict(support_metric=key,error_metric=metric,slope=np.nan,lower=np.nan,upper=np.nan,
                        bootstrap_valid=0,N_used=0,inference_status="unavailable_raw_count_diagnostic")
                        for metric in ("Brier","LogLoss"))
            continue
        support=np.asarray(diagnostics[key],float)
        valid=qc&np.isfinite(support)&(support>=0)
        x=np.log1p(np.maximum(np.nan_to_num(support),0))
        losses=dict(Brier=(y-p)**2,LogLoss=abm.logloss_rows(y,p))
        values=dict(x=x,x2=x*x,**losses,**{"x_"+name:x*loss for name,loss in losses.items()})
        plan=FamilyCells(groups,{"all":valid},w,values,labels,bins)
        den,num=plan.sums(np.ones(plan.families));spec=[]
        for metric in losses:
            row=dict(support_metric=key,transformation="log1p",risk_adjustment="frozen development10-risk-bin fixed effects",
                error_metric=metric,slope=slope(den,num,metric),N_used=int(np.sum(valid&(w>0))),
                risk_bins=10,source="pure abm count geometry",bootstrap_unit="family")
            rows.append(row);spec.append((row,metric,[]))
        plans.append((plan,spec))
    if not plans:return rows
    rng=np.random.default_rng(seed+444);family_n=plans[0][0].families
    for _ in range(B):
        mult=np.bincount(rng.integers(0,family_n,size=family_n),minlength=family_n)
        for plan,spec in plans:
            den,num=plan.sums(mult)
            for row,metric,draws in spec:draws.append(slope(den,num,metric))
    for plan,spec in plans:
        for row,metric,draws in spec:row.update(_interval(draws,row["slope"],B,family_n))
    return rows


def write_common_results(directory,bundle,tables):
    """Keep aggregate report compatibility; participant arrays remain private."""
    import shutil,tempfile
    from le8_common import write_tables
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    # A fresh atomically published directory prevents older packed workbooks or
    # old P1/P2 tables from surviving beside the new selector's current results.
    with tempfile.TemporaryDirectory(prefix=".common-reliability-",dir=directory.parent) as tmp:
        stage=Path(tmp)/"result";stage.mkdir()
        for path in directory.iterdir():
            if path.is_file() and ("LOCAL_ONLY" in path.name or path.name=="FROZEN_BEFORE_TEST.json"):
                shutil.copy2(path,stage/path.name)
        for name,frame in tables.items():frame.to_csv(stage/(name+".tsv"),sep="\t",index=False)
        scope=directory.parent.parent
        write_tables(stage,tables,"c1_correlate",scope.parent.name,scope.parent.parent.name,
            VERSION,workbook="common_reliability.results.xlsx")
        manifest=dict(bundle["manifest"],test_outcomes_used=True,
            inference="paired family bootstrap conditional on frozen risk models, clinical OOF fits and common selector")
        abm.json_write(stage/"reliability_contract.json",manifest)
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in stage.iterdir() if p.is_file()}
        abm.json_write(stage/"COMPLETE.json",dict(version=VERSION,status=manifest["status"],files=files))
        previous=Path(tmp)/"previous";directory.rename(previous)
        try:stage.rename(directory)
        except BaseException:
            previous.rename(directory)
            raise


def run_existing_scope(scope,bootstrap,seed,*,method_paths=None,prepared_source=None,protocol_path=None):
    """Standalone reliability uses the same new path, without refitting risk models."""
    from types import SimpleNamespace
    from c1_evaluation import locate_method,locate_prepared_source,load_km
    from c1_correlate import completion
    scope=Path(scope)
    protocol=json.loads(Path(protocol_path).read_text()) if protocol_path else {}
    if protocol.get("version",VERSION)!=VERSION:
        raise ValueError("Standalone reliability requires the current exploratory CVD protocol")
    if method_paths is None:
        method_paths={name:locate_method(scope,name) for name in
            ("clinical","glmnet","lightgbm","abm","abm_clinical","abm_deep")}
    paths={name:Path(path) for name,path in method_paths.items() if path is not None}
    if not paths:return dict(status="skipped_no_frozen_models",scope=str(scope))
    for name,path in paths.items():
        if not completion(path):raise ValueError(name+": frozen model integrity check failed")
    anchor="glmnet" if "glmnet" in paths else "clinical" if "clinical" in paths else next(iter(paths))
    predictions={name:joblib.load(path/"test_predictions_LOCAL_ONLY.joblib") for name,path in paths.items()}
    ids=np.asarray(predictions[anchor]["eid"]).astype(str)
    source=Path(prepared_source) if prepared_source is not None else locate_prepared_source(scope,ids)
    manifest=json.loads((paths[anchor]/"MODEL_FROZEN.json").read_text())
    cfg=abm.Config(**manifest["config"])
    if cfg.grouping_unit!="family":raise ValueError("Verified family grouping is required")
    km=load_km(paths[anchor])
    with np.load(source/"development_LOCAL_ONLY.npz",allow_pickle=False) as d:
        did=d["eid"].astype(str);dg=d["group"].astype(str);roles=d["role"].astype(str)
        gate=roles=="tune_gate"
        gy,gw=km.labels_weights(d["time"][gate],d["event"][gate])
    with np.load(source/"test_X_LOCAL_ONLY.npz",allow_pickle=False) as d:
        tid=d["eid"].astype(str);tg=d["group"].astype(str)
    if not np.array_equal(ids,tid):raise ValueError("Prepared test roster and frozen methods differ")
    all_ids=np.r_[did,tid];groups=np.r_[dg,tg];all_roles=np.r_[roles,np.repeat("test",len(tid))]
    check=pd.DataFrame(dict(group=groups,role=all_roles)).groupby("group").role.nunique()
    if (check>1).any():raise ValueError("A family crosses selector development/evaluation roles")
    parts={role:np.flatnonzero(all_roles==role) for role in (*abm.ROLES,"test")}
    y=np.zeros(len(all_ids),int);w=np.zeros(len(all_ids))
    gate_projection=joblib.load(paths[anchor]/"tune_gate_predictions_LOCAL_ONLY.joblib")
    gate_projection=_align(all_ids[parts["tune_gate"]],gate_projection)
    y[parts["tune_gate"]]=gy;w[parts["tune_gate"]]=gw*np.asarray(gate_projection["technical_qc"],bool)
    c=SimpleNamespace(ids=all_ids,groups=groups,parts=parts,y=y,w=w,cfg=cfg)
    out=scope/"qc/reliability"
    bundle=freeze_common_selector(c,paths,out,protocol_path)
    # First semantic test-outcome read is strictly after FROZEN_BEFORE_TEST.
    truth=abm.read_npz(source/"test_truth_LOCAL_ONLY.npz")
    tables=evaluate_common_selector(bundle,predictions,truth,km,bootstrap=bootstrap,seed=seed)
    write_common_results(out,bundle,tables)
    return dict(status=bundle["manifest"]["status"],scope=str(scope),output=str(out),
        version=VERSION,risk_models_refitted=False)
