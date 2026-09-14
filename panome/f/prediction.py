"""Locked train/validation/test survival benchmarks and marginal PWAS baseline."""
import json
import joblib
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from sksurv.util import Surv
from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored, concordance_index_ipcw, cumulative_dynamic_auc, brier_score
from sksurv.nonparametric import kaplan_meier_estimator
from data import design,csv_list


def cindex(y,score):
    return float(concordance_index_censored(y['event'],y['time'],score)[0])


def marginal_cox(x,y,features):
    """Vectorized univariate Breslow Cox MLEs, training only; no clinical adjustment."""
    order=np.argsort(y['time']);t=y['time'][order];e=y['event'][order]
    starts=np.r_[0,np.flatnonzero(np.diff(t))+1];deaths=np.add.reduceat(e.astype(float),starts)
    rows=[]
    for start in range(0,x.shape[1],64):
        z=x[order,start:start+64].astype(float);beta=np.zeros(z.shape[1]);event_sum=z[e].sum(axis=0)
        for it in range(40):
            eta=z*beta;offset=eta.max(axis=0);w=np.exp(eta-offset)
            s0=np.cumsum(w[::-1],axis=0)[::-1][starts]
            s1=np.cumsum((w*z)[::-1],axis=0)[::-1][starts]
            s2=np.cumsum((w*z*z)[::-1],axis=0)[::-1][starts]
            mu=s1/s0
            gradient=event_sum-(deaths[:,None]*mu).sum(axis=0)
            information=(deaths[:,None]*np.maximum(s2/s0-mu*mu,0)).sum(axis=0)
            step=np.clip(gradient/np.maximum(information,1e-9),-.5,.5)
            beta+=step
            if np.max(np.abs(step))<1e-7:break
        se=1/np.sqrt(np.maximum(information,1e-9));zs=beta/se
        for j in range(z.shape[1]):
            rows.append(dict(feature=features[start+j],beta=beta[j],se=se[j],z=zs[j],p=2*norm.sf(abs(zs[j])),
                             converged=bool(abs(step[j])<1e-5)))
    frame=pd.DataFrame(rows);order=np.argsort(frame.p.values);q=np.empty(len(frame))
    q[order]=np.minimum.accumulate((frame.p.values[order]*len(frame)/np.arange(1,len(frame)+1))[::-1])[::-1]
    frame['fdr']=np.minimum(q,1)
    return frame


def predict_all(x,p,z,pca,g,features,a,out):
    part=p.split.values;tr=part=='train';va=part=='validation';te=part=='test'
    y=Surv.from_arrays(p.event.astype(bool),p.time)
    covs=csv_list(a.covariates)
    transformer=design(covs,csv_list(a.categorical))
    c=transformer.fit_transform(p.loc[tr]);c_all=transformer.transform(p)
    nonconstant=np.std(c,axis=0)>1e-8;c_all=c_all[:,nonconstant]
    joblib.dump(dict(transformer=transformer,nonconstant=nonconstant),out/'clinical_preprocessor.joblib')
    pwas=marginal_cox(x[tr],y[tr],features)
    pwas.to_csv(out/'training_marginal_pwas.csv',index=False)
    selected=pwas.loc[pwas.converged].sort_values('p').head(a.pwas_top)
    if selected.empty: raise ValueError('No converged marginal Cox coefficients')
    lookup={s:i for i,s in enumerate(features)}
    selected_ix=np.array([lookup[s] for s in selected.feature]);score=x[:,selected_ix]@selected.beta.values
    joblib.dump(dict(indices=selected_ix,beta=selected.beta.values),out/'pwas_score.joblib')
    # Largest discovery state is the neutral numerical reference; no outcome-selected labeling.
    hard=np.eye(g['soft'].shape[1])[g['state']][:,1:]
    blocks={'clinical':c_all,'clinical_pwas':np.c_[c_all,score],
            'clinical_pca':np.c_[c_all,pca], 'clinical_ae':np.c_[c_all,z],
            'clinical_state':np.c_[c_all,hard], 'clinical_soft_state':np.c_[c_all,g['soft'][:,1:]],
            'clinical_diffusion':np.c_[c_all,g['diffusion']],
            'clinical_panome':np.c_[c_all,z,g['soft'][:,1:]],
            'clinical_elasticnet':np.c_[c_all,x]}
    metrics=[];tuning=[];test=p.loc[te,[a.id_col,'time','event']].copy();test['state']=g['state'][te]+1
    test['novel']=g['novelty'][te]
    predicted={};calibration=[];errors=[]
    for name,block in blocks.items():
        print('Survival benchmark:',name,flush=True)
        scaler=StandardScaler().fit(block[tr]);b=scaler.transform(block)
        active=np.std(b[tr],axis=0)>1e-8;b=b[:,active]
        if name=='clinical_elasticnet':
            model=CoxnetSurvivalAnalysis(l1_ratio=.5,n_alphas=a.coxnet_alphas,alpha_min_ratio=.01,
                    max_iter=100000,fit_baseline_model=True).fit(b[tr],y[tr])
            scores=[cindex(y[va],model.predict(b[va],alpha=alpha)) for alpha in model.alphas_]
            best=int(np.argmax(scores));alpha=float(model.alphas_[best])
            for aa,ss in zip(model.alphas_,scores):tuning.append(dict(model=name,alpha=aa,validation_c=ss))
            pred=model.predict(b[te],alpha=alpha)
            survival=model.predict_survival_function(b[te],alpha=alpha)
        else:
            candidates=[]
            for alpha in [0.1,1.,10.,100.]:
                m=CoxPHSurvivalAnalysis(alpha=alpha,ties='breslow',n_iter=200).fit(b[tr],y[tr])
                value=cindex(y[va],m.predict(b[va]));candidates.append((value,alpha,m))
                tuning.append(dict(model=name,alpha=alpha,validation_c=value))
            value,alpha,model=max(candidates,key=lambda v:v[0])
            pred=model.predict(b[te]);survival=model.predict_survival_function(b[te])
        joblib.dump(dict(scaler=scaler,active=active,model=model,alpha=alpha),out/f'{name}.joblib')
        predicted[name]=pred;test[name+'_log_hazard']=pred
        metrics.append(dict(model=name,metric='Harrell_C',horizon=0,value=cindex(y[te],pred),n=int(te.sum()),events=int(y[te]['event'].sum())))
        # Truncate held-out follow-up inside training support for IPCW evaluation.
        tau=min(float(np.quantile(y[tr]['time'],.99)),float(y[te]['time'].max()))-1e-6
        yt=y[te].copy();yt['event']&=yt['time']<=tau;yt['time']=np.minimum(yt['time'],tau)
        for horizon in a.horizons:
            if not (yt['time'].min()<horizon<tau):
                errors.append(dict(model=name,horizon=horizon,reason=f'Outside evaluable follow-up support (tau={tau})'));continue
            try:
                surv=np.array([fn(horizon) for fn in survival]);risk=1-surv
                auc=float(cumulative_dynamic_auc(y[tr],yt,pred,[horizon])[0][0])
                bs=float(brier_score(y[tr],yt,surv[:,None],[horizon])[1][0])
                uno=float(concordance_index_ipcw(y[tr],yt,pred,tau=horizon)[0])
                test[f'{name}_net_risk_{horizon:g}y']=risk
                for metric,value in [('AUC_IPCW',auc),('Brier_IPCW',bs),('Uno_C',uno)]:
                    metrics.append(dict(model=name,metric=metric,horizon=horizon,value=value,n=len(yt),events=int((yt['event']&(yt['time']<=horizon)).sum())))
                groups=pd.qcut(pd.Series(risk),5,labels=False,duplicates='drop').to_numpy()
                for group in np.unique(groups[np.isfinite(groups)]):
                    ix=groups==group;tt,ss=kaplan_meier_estimator(yt['event'][ix],yt['time'][ix])
                    supported=yt['time'][ix].max()>=horizon
                    observed=1-ss[np.searchsorted(tt,horizon,side='right')-1] if (tt<=horizon).any() else 0.
                    calibration.append(dict(model=name,horizon=horizon,quintile=int(group)+1,n=int(ix.sum()),
                         predicted=float(risk[ix].mean()),observed_net_risk=float(observed) if supported else np.nan,
                         at_risk=int((yt['time'][ix]>=horizon).sum())))
            except ValueError as exc: errors.append(dict(model=name,horizon=horizon,reason=str(exc)))
    pd.DataFrame(metrics).to_csv(out/'test_metrics.csv',index=False)
    pd.DataFrame(tuning).to_csv(out/'validation_tuning.csv',index=False)
    test.to_csv(out/'test_predictions.csv',index=False)
    pd.DataFrame(calibration).to_csv(out/'test_calibration.csv',index=False)
    (out/'metric_limitations.json').write_text(json.dumps(errors,indent=2))
    # Paired test bootstrap: uncertainty conditional on these trained models.
    rng=np.random.default_rng(a.seed);boot=[];yt=y[te]
    for rep in range(a.bootstrap):
        ix=rng.integers(0,len(yt),len(yt))
        if not yt['event'][ix].any():continue
        try:
            base=cindex(yt[ix],predicted['clinical'][ix])
            for name,pred in predicted.items():
                value=cindex(yt[ix],pred[ix]);boot.append(dict(replicate=rep,model=name,c=value,delta_vs_clinical=value-base))
        except ValueError:continue
    pd.DataFrame(boot).to_csv(out/'paired_bootstrap.csv',index=False)
    if boot:
        ci=pd.DataFrame(boot).groupby('model')[['c','delta_vs_clinical']].quantile([.025,.975]).unstack()
        ci.columns=['_'.join(map(str,c)) for c in ci.columns];ci.to_csv(out/'test_c_intervals.csv')
    # Fixed 2y landmark sensitivity, using locked baseline scores; no retraining or post-baseline X.
    lag=a.lag_years;eligible=yt['time']>lag
    rows=[]
    if eligible.sum()>10 and yt['event'][eligible].sum()>1:
        yl=yt[eligible].copy();yl['time']-=lag
        for name,pred in predicted.items():rows.append(dict(model=name,lag_years=lag,n=int(eligible.sum()),events=int(yl['event'].sum()),Harrell_C=cindex(yl,pred[eligible])))
    pd.DataFrame(rows).to_csv(out/'landmark_sensitivity.csv',index=False)
