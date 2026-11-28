"""Validation-only probability calibration and held-out top-risk enrichment."""
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logit

from common import outdir
from data import load_meta
from metrics import CensoringKM, censor_for, ipcw_weights
from train import risk_cache


def calibration_slope(p, y, weight):
    x = logit(np.clip(p, 1e-7, 1-1e-7))
    use = np.isfinite(x) & (weight > 0)
    x,y,w = x[use], y[use].astype(float), weight[use]
    if y.sum() < 10 or len(y)-y.sum() < 10 or np.ptp(x) == 0:
        return None
    w = w/w.sum()
    def objective(par):
        eta = par[0] + par[1]*x
        d = w*(expit(eta)-y)
        return w@(np.logaddexp(0,eta)-y*eta), np.array([d.sum(),d@x])
    fit = minimize(objective,[0.,1.],jac=True,method='L-BFGS-B')
    return fit.x if fit.success and np.isfinite(fit.x).all() else None


def profiles(cfg):
    folder = outdir(cfg,'main/results')
    meta = load_meta(cfg)
    train, _ = risk_cache(cfg,'train')
    valid, _ = risk_cache(cfg,'validation')
    test, _ = risk_cache(cfg,'test')
    vpred = np.load(folder/'flux.validation.npy',mmap_mode='r')
    tpred = np.load(folder/'flux.test.npy',mmap_mode='r')
    horizon = float(cfg['readout']['horizon_years'])
    rows, groups = [], []
    for j,name in enumerate(meta['endpoint_names']):
        def data(a):
            t=a['event_time'][:,j]
            c=censor_for(a,name,meta)
            m=a['risk_mask'][:,j] & (c>0)
            return m, np.minimum(t,c)[m], (t<=c)[m], t[m]<=horizon
        mt, time, event, _ = data(train)
        if not len(time):
            continue
        km=CensoringKM.fit(time,event)
        mv,tv,ev,_ = data(valid)
        mv_cases,_,vw,vs=ipcw_weights(tv,ev,horizon,km,weight_cap=cfg['evaluation']['weight_cap'])
        m,t,e,y = data(test)
        p=np.asarray(tpred[:,j])[m]
        if not len(p) or not np.isfinite(p).all():
            continue
        cases,_,w,ts=ipcw_weights(t,e,horizon,km,weight_cap=cfg['evaluation']['weight_cap'])
        par=calibration_slope(vpred[:,j][mv],mv_cases,vw) if vs['status']=='ok' else None
        row=dict(endpoint=name,n=len(p),validation_cases=int(mv_cases.sum()),status='insufficient_support')
        if par is not None and ts['status']=='ok':
            new=expit(par[0]+par[1]*logit(np.clip(p,1e-7,1-1e-7)))
            row.update(status='ok',intercept=par[0],slope=par[1],
                       raw_ipcw_brier=np.mean(w*(cases-p)**2), calibrated_ipcw_brier=np.mean(w*(cases-new)**2))
        rows.append(row)
        # Stable eid order breaks score ties; observed incidence is recorded-event incidence.
        order=np.argsort(-p,kind='stable')
        for frac in (0.01,0.1):
            k=int(np.ceil(frac*len(p)))
            a,b=y[order[:k]],y[order[k:]]
            other=b.mean() if len(b) else np.nan
            groups.append(dict(endpoint=name,top_fraction=frac,n=len(p),cases=int(y.sum()),top_n=k,
                top_cases=int(a.sum()),top_incidence=a.mean(),overall_incidence=y.mean(),remaining_incidence=other,
                risk_ratio=a.mean()/other if other>0 else np.nan,case_capture=a.sum()/y.sum() if y.any() else np.nan))
    pd.DataFrame(rows,columns=['endpoint','n','validation_cases','status','intercept','slope','raw_ipcw_brier','calibrated_ipcw_brier']).to_csv(folder/'flux.validation_calibration.tsv',sep='\t',index=False)
    pd.DataFrame(groups,columns=['endpoint','top_fraction','n','cases','top_n','top_cases','top_incidence','overall_incidence','remaining_incidence','risk_ratio','case_capture']).to_csv(folder/'flux.enrichment.tsv',sep='\t',index=False)
