"""Locked GRID evaluation. All methods use identical frozen subset masks.

No fitting, calibration, threshold selection, or test-driven model selection occurs
here. Regression R2 is 1-SSE/SST, never squared residual correlation.
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from twin_core import hash_number, state_hash, id_hash

METRIC_DIRECTIONS={"MSE":"lower","RMSE":"lower","MAE":"lower","prediction_R2":"higher","MSE_over_train_variance":"lower","AUC":"higher","AUPRC":"higher","Brier":"lower","LogLoss":"lower"}


def metrics(y,p,binary,training_variance,weights=None):
    y=np.asarray(y,float);p=np.asarray(p,float);w=np.ones(len(y)) if weights is None else np.asarray(weights,float)
    good=w>0;y,p,w=y[good],p[good],w[good];n=float(w.sum())
    if not n:return {k:np.nan for k in (('AUC','AUPRC','Brier','LogLoss') if binary else ('MSE','RMSE','MAE','prediction_R2','MSE_over_train_variance'))}
    if not np.isfinite(y).all() or not np.isfinite(p).all():raise ValueError('Nonfinite prediction inside a locked common evaluation mask')
    error=y-p;mse=float(np.average(error**2,weights=w))
    if binary:
        if np.any((p<0)|(p>1)):raise ValueError('Binary benchmarks require probabilities, not PRS values/logits')
        pp=np.clip(p,1e-7,1-1e-7);two=len(np.unique(y))==2
        return dict(AUC=float(roc_auc_score(y,p,sample_weight=w)) if two else np.nan,
                    AUPRC=float(average_precision_score(y,p,sample_weight=w)) if two else np.nan,
                    Brier=mse,LogLoss=float(np.average(-(y*np.log(pp)+(1-y)*np.log1p(-pp)),weights=w)))
    var=float(np.average((y-np.average(y,weights=w))**2,weights=w))
    return dict(MSE=mse,RMSE=math.sqrt(mse),MAE=float(np.average(np.abs(error),weights=w)),
                prediction_R2=1-mse/var if var>0 else np.nan,
                MSE_over_train_variance=mse/training_variance if training_variance>0 else np.nan)


class FixedMetricCache:
    """Pre-sort frozen prediction ties once; bootstrap changes weights, not scores.
    Weighted AUC/AP then use O(n) cumulative sums instead of O(n log n) sorting
    in every bootstrap draw. Losses and outcome moments are also precomputed.
    """
    def __init__(self,y,p,binary,training_variance):
        self.y=np.asarray(y,float);self.p=np.asarray(p,float);self.binary=binary;self.tv=training_variance
        if not np.isfinite(self.y).all() or not np.isfinite(self.p).all():raise ValueError('Invalid common-mask metric input')
        self.err2=(self.y-self.p)**2;self.abserr=np.abs(self.y-self.p)
        if binary:
            if not np.isin(self.y,[0,1]).all() or np.any((self.p<0)|(self.p>1)):raise ValueError('Invalid binary probability metric')
            values,self.bin=np.unique(self.p,return_inverse=True);self.nb=len(values)
            pp=np.clip(self.p,1e-7,1-1e-7);self.logloss=-(self.y*np.log(pp)+(1-self.y)*np.log1p(-pp))
    def __call__(self,w):
        w=np.asarray(w,float);total=w.sum()
        keys=('AUC','AUPRC','Brier','LogLoss') if self.binary else ('MSE','RMSE','MAE','prediction_R2','MSE_over_train_variance')
        if not total:return {k:np.nan for k in keys}
        mse=float(w@self.err2/total)
        if self.binary:
            pos=np.bincount(self.bin,weights=w*self.y,minlength=self.nb)
            neg=np.bincount(self.bin,weights=w*(1-self.y),minlength=self.nb)
            np_,nn=pos.sum(),neg.sum();auc=ap=np.nan
            if np_>0 and nn>0:
                auc=float(pos@(np.cumsum(neg)-.5*neg)/(np_*nn))
                tp=np.cumsum(pos[::-1]);fp=np.cumsum(neg[::-1]);den=tp+fp
                precision=np.divide(tp,den,out=np.zeros_like(tp),where=den>0)
                ap=float(pos[::-1]@precision/np_)
            return dict(AUC=auc,AUPRC=ap,Brier=mse,LogLoss=float(w@self.logloss/total))
        mean=float(w@self.y/total);var=float(w@((self.y-mean)**2)/total)
        return dict(MSE=mse,RMSE=math.sqrt(mse),MAE=float(w@self.abserr/total),prediction_R2=1-mse/var if var>0 else np.nan,MSE_over_train_variance=mse/self.tv if self.tv>0 else np.nan)


def composition_control(selected,eligible,risk,ancestry,ids,edges,seed):
    """Outcome-blind descriptive control; batch-composition dependent, NOT deployment."""
    bins=np.searchsorted(edges,risk);mask=np.zeros(len(ids),bool)
    for pop in np.unique(ancestry):
        for b in np.unique(bins[ancestry==pop]):
            scope=(ancestry==pop)&(bins==b);ix=np.flatnonzero(eligible&scope)
            order=sorted(ix,key=lambda i:hash_number(ids[i],seed));count=int(np.sum(selected&scope))
            mask[np.array(order[:count],int)]=True
    return mask


def evaluate(model,prediction,truth,*,external=None,bootstrap=None):
    """truth: eid,y and optionally ancestry. All test subjects remain in the flow."""
    cfg=model.config
    if prediction['frozen_hash']!=state_hash(model):raise ValueError('Frozen model identity changed')
    truth=pd.DataFrame(truth).copy();truth['eid']=truth.eid.astype(str)
    ids=prediction['eid'];groups=prediction['group'];n=len(ids)
    if truth.eid.duplicated().any() or set(truth.eid)!=set(ids):raise ValueError('Truth must cover exactly the frozen query roster')
    truth=truth.set_index('eid').loc[ids]
    y=pd.to_numeric(truth.y,errors='raise').to_numpy(float);known=np.isfinite(y)
    if cfg.trait_type=='binary' and not np.isin(y[known],[0,1]).all():raise ValueError('Binary endpoint must be coded 0/1')
    pop=truth.ancestry.fillna('UNASSIGNED').astype(str).to_numpy() if 'ancestry' in truth else np.full(n,'UNASSIGNED')
    qc=prediction['diagnostics']['technical_qc'];base_scope=known&qc
    if base_scope.sum()<2:raise ValueError('Insufficient evaluable technical-QC passing outcomes')
    p={k:np.asarray(v,float) for k,v in prediction['predictions'].items()};availability=[]
    for name,values in (external or {}).items():
        if name in p:raise ValueError('External method collides with an internal method name')
        values=np.asarray(values,float)
        if len(values)!=n or not np.isfinite(values[base_scope]).all():
            availability.append(dict(method=name,status='unavailable_incomplete_frozen_roster'));continue
        if cfg.trait_type=='binary' and np.any((values[base_scope]<0)|(values[base_scope]>1)):
            raise ValueError(f'{name}: expected probabilities, not uncalibrated scores')
        p[name]=values;availability.append(dict(method=name,status='same_roster_prediction_available'))
    for name,v in p.items():
        if not np.isfinite(v[base_scope]).all():raise ValueError(f'{name}: incomplete common-mask predictions; do not silently use pairwise intersections')
    selected=prediction['selected'];masks={'all':base_scope,'selected':base_scope&selected,'rejected':base_scope&~selected}
    for name,v in prediction['masks'].items():masks[name]=base_scope&v
    # Count matching uses predictor values only, before looking at y. Missing-y
    # proportions are then reported rather than balancing using test outcomes.
    mask_random=composition_control(selected,qc,np.nan_to_num(p[model.fallback]),pop,ids,model.risk_edges,cfg.seed+471)
    masks['risk_ancestry_matched_random']=base_scope&mask_random
    pool=np.flatnonzero(qc);risk=np.nan_to_num(p[model.fallback]);center=.5 if cfg.trait_type=='binary' else model.prediction_center
    order=sorted(pool,key=lambda i:(-abs(risk[i]-center),hash_number(ids[i],cfg.seed+472)))
    cc=np.zeros(n,bool);cc[np.asarray(order[:int(selected.sum())],int)]=True
    masks['confidence_count_matched']=base_scope&cc
    scopes=['ALL',*sorted(set(pop))];rows=[];lookup={};binary=cfg.trait_type=='binary'
    for ancestry in scopes:
        am=np.ones(n,bool) if ancestry=='ALL' else pop==ancestry
        total=int(am.sum());cases=int(np.sum(am&known&(y==1))) if binary else None
        for subset,mask in masks.items():
            use=mask&am;count=int(use.sum())
            for name,values in p.items():
                met=metrics(y[use],values[use],binary,model.training_variance)
                row=dict(trait=cfg.trait,ancestry=ancestry,subset=subset,method=name,N=count,original_N=total,
                         actual_coverage=count/total if total else np.nan,selected_policy_coverage=float(np.sum(am&selected)/total) if total else np.nan,
                         cases=int(np.sum(y[use]==1)) if binary else np.nan,controls=int(np.sum(y[use]==0)) if binary else np.nan,
                         case_coverage=float(np.sum(y[use]==1)/cases) if binary and cases else np.nan,
                         outcome_variance=float(np.var(y[use])) if count else np.nan,
                         status='estimated' if count else 'empty',evaluation_scope=prediction['evaluation_scope'],**met)
                rows.append(row);lookup[(ancestry,subset,name)]=met
    primary=cfg.primary;full='ABM_full';comparisons=[
        ('selection_effect',full,'selected',full,'all'),
        ('training_effect_same_subset',primary,'selected',full,'selected'),
        ('combined_effect',primary,'selected',full,'all'),
        ('training_effect_all',primary,'all',full,'all'),
        ('pruning_specificity',primary,'selected','ABM_random','selected'),
        ('beyond_baseline_confidence',full,'selected',full,'confidence_count_matched'),
        ('beyond_risk_ancestry_composition',full,'selected',full,'risk_ancestry_matched_random'),
        ('fallback_policy_all','ABM_policy','all',model.fallback,'all')]
    benchmark_names=[name for name in p if name not in {primary,full,'ABM_soft','ABM_random','Clinical','ABM_policy'}]
    for name in benchmark_names:
        comparisons.extend([(f'benchmark_selection:{name}',name,'selected',name,'all'),
                            (f'benchmark_same_subset:{name}',primary,'selected',name,'selected'),
                            (f'benchmark_combined:{name}',primary,'selected',name,'all')])
    B=int(cfg.bootstrap if bootstrap is None else bootstrap)
    if B<0:raise ValueError('Bootstrap count must be nonnegative')
    metric_names=list(next(iter(lookup.values())))
    contrasts=[]
    # Cluster resampling keeps each selected group nested inside each original draw.
    for ancestry in scopes:
        am=np.ones(n,bool) if ancestry=='ALL' else pop==ancestry
        scope_indices=np.flatnonzero(am);ug,gi=np.unique(groups[am],return_inverse=True)
        if not len(ug):continue
        needed=set((method,subset) for _,a,b,c,d in comparisons for method,subset in [(a,b),(c,d)])
        scoped_masks={subset:masks[subset][am] for _,subset in needed}
        yy=y[am];vv={k:v[am] for k,v in p.items()}
        cached={(method,subset):FixedMetricCache(yy[scoped_masks[subset]],vv[method][scoped_masks[subset]],binary,model.training_variance)
                for method,subset in needed}
        collected={(label,metric):[] for label,*_ in comparisons for metric in metric_names}
        rng=np.random.default_rng(cfg.seed+702)
        for _ in range(B):
            gw=np.bincount(rng.integers(0,len(ug),len(ug)),minlength=len(ug));w=gw[gi]
            bm={}
            for method,subset in needed:
                use=scoped_masks[subset];bm[(method,subset)]=cached[(method,subset)](w[use])
            for label,a,s,b,t in comparisons:
                for metric in metric_names:
                    collected[(label,metric)].append(bm[(a,s)][metric]-bm[(b,t)][metric])
        for label,a,s,b,t in comparisons:
            for metric in metric_names:
                samples=np.asarray(collected[(label,metric)]);samples=samples[np.isfinite(samples)]
                enough=len(samples)>=max(20,math.ceil(.8*B));bounds=np.quantile(samples,[.025,.975]) if enough else [np.nan,np.nan]
                contrasts.append(dict(trait=cfg.trait,ancestry=ancestry,contrast=label,method=a,subset=s,reference=b,reference_subset=t,metric=metric,
                                      delta=lookup[(ancestry,s,a)][metric]-lookup[(ancestry,t,b)][metric],lower=bounds[0],upper=bounds[1],valid_boot=len(samples),requested_boot=B,
                                      direction=METRIC_DIRECTIONS[metric],grouping_status=cfg.grouping_status,
                                      inference='conditional_on_frozen_models_and_thresholds',multiplicity='primary contrast prespecified; others exploratory',
                                      status='interval_estimated' if enough else 'point_only_or_insufficient_bootstrap'))
    confidence_rows=[];confidence=prediction['diagnostics']['pre_prediction_confidence']
    for ancestry in scopes:
        am=np.ones(n,bool) if ancestry=='ALL' else pop==ancestry
        for lo in range(0,100,10):
            mask=base_scope&am&(confidence>=lo)&((confidence<lo+10) if lo<90 else confidence<=100)
            if mask.sum()<10:continue
            for name in (primary,full,model.fallback):
                row=dict(trait=cfg.trait,ancestry=ancestry,method=name,confidence_lo=lo,confidence_hi=lo+10,N=int(mask.sum()),
                         mean_confidence=float(confidence[mask].mean()),matching_supported_fraction=float(prediction['diagnostics']['matching_supported'][mask].mean()),
                         **metrics(y[mask],p[name][mask],binary,model.training_variance))
                confidence_rows.append(row)
    state_metric='Brier' if binary else 'MSE'
    tab=pd.DataFrame(contrasts);check=tab[(tab.ancestry=='ALL')&(tab.metric==state_metric)]
    selected_known=masks['selected'];enough_info=selected_known.sum()>=30
    if binary:enough_info &= min(np.sum(y[selected_known]==v) for v in (0,1))>=10
    def supported(label):
        d=check[check.contrast==label]
        return bool(enough_info and len(d)==1 and d.iloc[0].valid_boot>=200 and d.iloc[0].upper<0)
    status=dict(trait=cfg.trait,task=cfg.trait_type,evaluation_scope=prediction['evaluation_scope'],
                original_N=n,known_outcome_N=int(known.sum()),technical_QC_excluded=int((~qc).sum()),
                selected_N=int(selected.sum()),selected_with_observed_outcome_N=int(selected_known.sum()),actual_policy_coverage=float(selected.mean()),
                selection_loss_supported=supported('selection_effect'),training_loss_supported=supported('training_effect_same_subset'),
                beyond_confidence_loss_supported=supported('beyond_baseline_confidence'),primary_loss_metric=state_metric,
                grouping_status=cfg.grouping_status,clinical_release=False,external_validation=False,
                dataset_has_been_inspected=True,roster_hash=id_hash(ids),frozen_hash=prediction['frozen_hash'])
    flow=pd.DataFrame(dict(reason=prediction['reason'],outcome_observed=known)).groupby(['reason','outcome_observed'],dropna=False).size().reset_index(name='N')
    tables=dict(metrics=pd.DataFrame(rows),contrasts=tab,confidence_validation=pd.DataFrame(confidence_rows),
                cohort_flow=flow,status=pd.DataFrame([status]),comparator_availability=pd.DataFrame(availability))
    for d in tables.values():
        d['run_id']=prediction['frozen_hash'];d['trait']=cfg.trait;d['grouping_status']=cfg.grouping_status
    return tables,status
