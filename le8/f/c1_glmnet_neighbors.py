"""Add a held-out residual bank to the frozen native molecular LASSO."""
import numpy as np
import pandas as pd
import abm


def freeze(c,probability):
    """Called before opening test outcomes. Gate families are excluded from the native LASSO and calibrator."""
    donor=c.parts["tune_gate"]; donor=donor[c.w[donor]>0]
    tune=c.parts["tune_model"]; tune=tune[c.qc[tune]]
    q=c.parts["test"]
    audit=dict(donor_role="tune_gate",radius_role="tune_model_features",donor_N=len(donor),
               donor_families=len(np.unique(c.groups[donor])),donor_cases=int(c.y[donor].sum()),
               fallback="frozen calibrated glmnet",test_outcomes_used=False)
    if len(np.unique(c.y[donor]))<2 or len(donor)<c.cfg.k or len(tune)<10:
        return None,dict(audit,status="insufficient_development_support"),None
    bank=abm.TwinBank().fit(c.z[donor],c.y[donor],c.w[donor],c.ids[donor],c.groups[donor],
               c.y[donor]-probability[donor],c.cfg,c.common_geometry)
    try:
        bank.calibrate_radius(c.z[tune],c.ids[tune],c.groups[tune])
    except ValueError as error:
        if str(error) not in ("Not enough independent queries to calibrate matching radius","Collapsed matching geometry"):raise
        return None,dict(audit,status=str(error)),None
    _,diag=bank.predict(c.z[q],c.ids[q],c.groups[q],probability[q])
    use=diag["supported"] & c.qc[q]
    corrected=probability[q].copy()
    corrected[use]=np.clip(corrected[use]+diag["borrow_fraction"][use]*diag["donor_loss"][use],1e-8,1-1e-8)
    return dict(probability=corrected,supported=use,families=diag["family_count"]),dict(audit,status="frozen",
                    radius=bank.radius,prior_strength=c.cfg.prior_strength,test_supported_N=int(use.sum())),bank


def evaluate(predictions,truth,km,reps,seed):
    from c1_evaluation import _align
    from timing_study import compare, metrics, auc_layout
    native=predictions.get("glmnet",{}); arm=native.get("frozen_residual_arm")
    if arm is None:return pd.DataFrame()
    ids,g,time,event,models,qc,finite=_align(predictions,truth)
    order=pd.Index(np.asarray(native["eid"],str)).get_indexer(ids)
    if np.any(order<0):raise ValueError("Mismatched glmnet roster")
    corrected=np.asarray(arm["probability"])[order]; high=np.asarray(arm["supported"])[order]
    base=models["glmnet"]["probability"]; y,w=km.labels_weights(time,event)
    common=qc & finite & np.isfinite(corrected)
    out=[]
    for name,mask in (("all_common",common),("supported",common & high),("unsupported",common & ~high)):
        sub=mask & (w>0)
        if min((y[sub]==1).sum(),(y[sub]==0).sum())<2:continue
        if reps>0:
            rows=compare(y[sub],w[sub],g[sub],base[sub],corrected[sub],reps,"frozen_glmnet_neighbor_increment",
                         "glmnet","glmnet_plus_neighbors",seed=seed)
        else:
            a,b=[metrics(y[sub],w[sub],p[sub],*auc_layout(p[sub])) for p in (base,corrected)]
            rows=[dict(metric=k,reference_value=x,model_value=z,difference=z-x,lower=np.nan,upper=np.nan,
                       test_known_N=int(sub.sum()),test_cases=int(y[sub].sum())) for k,x,z in zip(
                 ("AUC_IPCW","Brier_IPCW","LogLoss_IPCW","predicted_minus_observed_risk"),a,b)]
        out.extend(dict(row,subset=name,eligible_N=int(mask.sum()),common_N=int(common.sum()),
             coverage=float(mask.sum()/common.sum()),residual_donor_families=native["metadata"]["residual_bank"]["donor_families"])
             for row in rows)
    return pd.DataFrame(out)
