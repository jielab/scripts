#!/usr/bin/env python3
from __future__ import annotations
import hashlib, numpy as np, pandas as pd
from c1_lock_roster import COVERAGES, PRIMARY_SELECTOR, eid_hash

ROUTES=("risk_resid_error","logloss_resid","local_stability","model_consensus","hybrid_predictability",
        "support_quality","absolute_error","rank_error","support_only","lowrisk","confidence","random")

def _stable_hash(value,salt):
    return int.from_bytes(hashlib.sha256(f"{salt}:{value}".encode()).digest()[:8],"big")

def _metrics(time,event,score,km):
    time=np.asarray(time,float);event=np.asarray(event,int);score=np.asarray(score,float)
    y,w=km.labels_weights(time,event);ok=(w>0)&np.isfinite(score)
    out=dict(N=int(ok.sum()),cases=int(y[ok].sum()),Brier=np.nan,LogLoss=np.nan,AUC=np.nan)
    if not ok.any(): return out
    p=np.clip(score[ok],1e-7,1-1e-7)
    out["Brier"]=float(np.average((y[ok]-p)**2,weights=w[ok]))
    out["LogLoss"]=float(np.average(-(y[ok]*np.log(p)+(1-y[ok])*np.log(1-p)),weights=w[ok]))
    if len(np.unique(y[ok]))==2:
        from sklearn.metrics import roc_auc_score
        out["AUC"]=float(roc_auc_score(y[ok],p,sample_weight=w[ok]))
    return out

def _risk_match(selected,qc,risk,ids,salt,bins=10):
    selected=np.asarray(selected,bool);qc=np.asarray(qc,bool);risk=np.asarray(risk,float);ids=np.asarray(ids,str)
    ok=qc&np.isfinite(risk);out=np.zeros(len(ids),bool)
    selected=selected&ok
    if not ok.any():return out
    q=np.unique(np.quantile(risk[ok],np.linspace(0,1,bins+1)))
    if len(q)<3: labels=np.zeros(len(ids),int)
    else:
        q[0],q[-1]=-np.inf,np.inf;labels=np.searchsorted(q[1:-1],risk,side="right")
    for b in np.unique(labels[ok]):
        pool=np.flatnonzero(ok&(labels==b));n=int(np.sum(selected&(labels==b)))
        order=sorted(pool,key=lambda j:_stable_hash(ids[j],salt))
        out[np.asarray(order[:n],int)]=True
    return out

def _aurc(cov,val):
    x=np.asarray(cov,float);y=np.asarray(val,float);ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<2 or not ok.all():return np.nan
    x,y=x[ok],y[ok];o=np.argsort(x);x,y=x[o],y[o]
    return float(np.trapezoid(y,x)/(x[-1]-x[0])) if x[-1]>x[0] else np.nan

def _finite_mean(values):
    values=np.asarray(values,float)
    finite=values[np.isfinite(values)]
    return float(finite.mean()) if len(finite) else np.nan

def infer_selector(prediction,time,event,groups,km,coverages=(.2,.4,.5,.6,.8,1.0),routes=ROUTES,B=1000,seed=2026,
                   grouping_unit="individual_unverified"):
    if grouping_unit not in {"family", "individual_unverified"}:
        raise ValueError("Invalid grouping unit")
    if tuple(coverages)!=COVERAGES or PRIMARY_SELECTOR not in routes:
        raise ValueError("C1 LOCKPAPER requires the locked coverage grid and primary selector")
    if int(B)!=B or B<1: raise ValueError("Bootstrap replicates must be positive")
    if pd.isna(np.asarray(groups)).any(): raise ValueError("Missing family/group identifiers")
    ids=np.asarray(prediction["eid"],str);groups=np.asarray(groups,str)
    time=np.asarray(time,float);event=np.asarray(event,int)
    qc=np.asarray(prediction["diagnostics"]["technical_qc"],bool)
    risk=np.asarray(prediction["probability"]["elasticnet_full"],float)
    if not len(ids) or len(set(ids))!=len(ids) or any(len(x)!=len(ids) for x in (time,event,groups,qc,risk)):
        raise ValueError("Invalid selector evaluation roster")
    if pd.isna(groups).any() or np.any(np.char.strip(groups)==""):
        raise ValueError("Missing family/group identifiers")
    qc=qc&np.isfinite(risk)
    masks={k:np.asarray(v,bool)&qc for k,v in prediction.get("masks",{}).items()}
    if any(f"{PRIMARY_SELECTOR}:{c:g}" not in masks for c in coverages):
        raise ValueError("Frozen primary selector must provide all six coverage masks")
    rows=[];controls={}
    for route in routes:
        for cov in coverages:
            key=f"{route}:{cov:g}"
            if key not in masks:continue
            sel=masks[key];rm=_risk_match(sel,qc,risk,ids,f"lockpaper:{route}:{cov:g}")
            controls[(route,float(cov))]=rm
            a=_metrics(np.asarray(time)[sel],np.asarray(event)[sel],risk[sel],km)
            b=_metrics(np.asarray(time)[rm],np.asarray(event)[rm],risk[rm],km)
            rows.append(dict(route=route,coverage_target=float(cov),coverage=float(sel.mean()),N=a["N"],cases=a["cases"],
              selected_N=int(sel.sum()),riskmatched_total_N=int(rm.sum()),
              selected_eid_hash=eid_hash(ids[sel]),riskmatched_eid_hash=eid_hash(ids[rm]),
              riskmatched_N=b["N"],riskmatched_cases=b["cases"],Brier=a["Brier"],riskmatched_Brier=b["Brier"],
              delta_Brier=a["Brier"]-b["Brier"],LogLoss=a["LogLoss"],riskmatched_LogLoss=b["LogLoss"],
              delta_LogLoss=a["LogLoss"]-b["LogLoss"],AUC=a["AUC"],riskmatched_AUC=b["AUC"],delta_AUC=a["AUC"]-b["AUC"]))
    curve=pd.DataFrame(rows)
    point=pd.DataFrame([dict(route=r,AURC_delta_Brier=_aurc(g.coverage_target,g.delta_Brier),
                             AURC_delta_LogLoss=_aurc(g.coverage_target,g.delta_LogLoss),
                             mean_delta_AUC=_finite_mean(g.delta_AUC))
                        for r,g in curve.groupby("route",sort=False)])
    u,inv=np.unique(groups,return_inverse=True);order=np.argsort(inv,kind="stable")
    pos=np.split(order,np.cumsum(np.bincount(inv))[:-1]);rng=np.random.default_rng(seed+9181)
    boots={r:[] for r in point.route}
    for _ in range(int(B)):
        ix=np.concatenate([pos[j] for j in rng.integers(0,len(u),len(u))])
        for route in boots:
            rr=[]
            for cov in coverages:
                key=f"{route}:{cov:g}"
                if key not in masks or (route,float(cov)) not in controls:continue
                sm=masks[key][ix];cm=controls[(route,float(cov))][ix]
                a=_metrics(np.asarray(time)[ix][sm],np.asarray(event)[ix][sm],risk[ix][sm],km)
                b=_metrics(np.asarray(time)[ix][cm],np.asarray(event)[ix][cm],risk[ix][cm],km)
                rr.append((float(cov),a["Brier"]-b["Brier"],a["LogLoss"]-b["LogLoss"],a["AUC"]-b["AUC"]))
            if len(rr)>=2:
                z=np.asarray(rr,float);boots[route].append((_aurc(z[:,0],z[:,1]),_aurc(z[:,0],z[:,2]),_finite_mean(z[:,3])))
    infer=[]
    minimum_valid=max(20,int(np.ceil(.8*B)))
    for _,p in point.iterrows():
        arr=np.asarray(boots[p.route],float)
        def ci(j):
            v=arr[:,j] if arr.ndim==2 else np.array([]);v=v[np.isfinite(v)]
            return (float(np.quantile(v,.025)),float(np.quantile(v,.975)),len(v)) if len(v)>=minimum_valid else (np.nan,np.nan,len(v))
        c0,c1,c2=ci(0),ci(1),ci(2)
        def inference_status(value,interval):
            return ("insufficient_information" if not np.isfinite(value) else
                    "evaluated" if interval[2]>=minimum_valid else "insufficient_bootstrap")
        s0=inference_status(p.AURC_delta_Brier,c0)
        s1=inference_status(p.AURC_delta_LogLoss,c1)
        s2=inference_status(p.mean_delta_AUC,c2)
        infer.append(dict(route=p.route,AURC_delta_Brier=p.AURC_delta_Brier,AURC_delta_Brier_lo=c0[0],AURC_delta_Brier_hi=c0[1],
                          AURC_delta_LogLoss=p.AURC_delta_LogLoss,AURC_delta_LogLoss_lo=c1[0],AURC_delta_LogLoss_hi=c1[1],
                          mean_delta_AUC=p.mean_delta_AUC,mean_delta_AUC_lo=c2[0],mean_delta_AUC_hi=c2[1],
                          valid_boot=min(c0[2],c1[2],c2[2]),requested_boot=int(B),minimum_valid_boot=minimum_valid,
                          Brier_bootstrap_valid=c0[2],LogLoss_bootstrap_valid=c1[2],AUC_bootstrap_valid=c2[2],
                          Brier_inference_status=s0,LogLoss_inference_status=s1,AUC_inference_status=s2,
                          status=s0 if s0!="evaluated" else s1 if s1!="evaluated" else s2,
                          bootstrap_unit=grouping_unit,
                          independent_family_support_verified=grouping_unit == "family",
                          uncertainty="conditional on frozen models, selector, risk-matched controls and censoring estimator"))
    return curve,pd.DataFrame(infer)
