"""Locked evaluation for TwinABM; never updates models, masks or thresholds."""
from __future__ import annotations
import math
from typing import Any
import numpy as np
from sklearn.metrics import roc_auc_score
from abm_twin_core import logloss_rows, pairwise_error_target, model_hash


def uno_c(time, event, probability, km):
    """Horizon-truncated Uno C with strict earlier event times and half risk ties.

    Marginal censoring assumption matches the fixed-horizon IPCW implementation.
    Equal event times are not ordered. Returns comparable-pair diagnostics.
    """
    time, event, prob = map(np.asarray, (time, event, probability))
    levels = np.unique(prob); rank = np.searchsorted(levels, prob) + 1
    tree = np.zeros(len(levels) + 1)
    def add(k):
        while k < len(tree): tree[k] += 1; k += k & -k
    def prefix(k):
        total = 0.
        while k > 0: total += tree[k]; k -= k & -k
        return total
    order = np.argsort(-time, kind="stable")
    num = den = 0.; later = pairs = start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and time[order[stop]] == time[order[start]]: stop += 1
        for i in order[start:stop]:
            if event[i] == 1 and time[i] <= km.horizon and later:
                g = float(km.at(np.array([time[i]]), left=True)[0])
                if g < km.minimum: return np.nan, pairs
                weight = 1/g**2
                low = prefix(int(rank[i]-1)); equal = prefix(int(rank[i])) - low
                num += weight*(low + .5*equal); den += weight*later; pairs += later
        for i in order[start:stop]: add(int(rank[i])); later += 1
        start = stop
    return (num/den if den else np.nan), pairs


def metrics(time, event, probability, km):
    n = len(time)
    if not n: return dict(N=0, AUC_IPCW=np.nan, Uno_C_horizon=np.nan, Brier_IPCW=np.nan,
                         LogLoss_IPCW=np.nan, cases=0, controls=0, status="empty_subset")
    time, event, probability = map(np.asarray, (time, event, probability))
    if not np.isfinite(probability).all(): raise ValueError("Nonfinite probabilities inside evaluation mask")
    y, w = km.labels_weights(time, event); ok = w > 0
    row = dict(N=n, known_N=int(ok.sum()), cases=int(y.sum()), controls=int(np.sum(ok & (y == 0))),
        status="completed", Brier_normalization="Hajek: sum(IPCW*loss)/sum(IPCW)")
    row["AUC_IPCW"] = float(roc_auc_score(y[ok], probability[ok], sample_weight=w[ok])) if len(set(y[ok])) == 2 else np.nan
    row["Uno_C_horizon"], row["comparable_pairs"] = uno_c(time, event, probability, km)
    row["Brier_IPCW"] = float(np.average((y[ok]-probability[ok])**2, weights=w[ok])) if ok.any() else np.nan
    row["LogLoss_IPCW"] = float(np.average(logloss_rows(y[ok], probability[ok]), weights=w[ok])) if ok.any() else np.nan
    if not np.isfinite(row["AUC_IPCW"]): row["status"] = "insufficient_case_control_information"
    return row


def risk_matched_mask(selected, eligible, risk, ids, risk_edges, seed):
    """X/prediction-only negative control: no outcomes enter its selection.

    This is a descriptive, coverage-composition-matched control on current query X,
    not the primary deployable policy. Its membership can depend on batch makeup.
    """
    from abm_twin_core import tie_hash
    bins = np.searchsorted(risk_edges, risk, side="right")
    mask = np.zeros(len(ids), bool)
    for label in np.unique(bins):
        pool = np.flatnonzero(eligible & (bins == label))
        order = sorted(pool, key=lambda j: tie_hash(str(ids[j]), seed))
        count = int(np.sum(selected & (bins == label)))
        mask[np.asarray(order[:count], int)] = True
    return mask


def evaluate(model, prediction, time, event, truth_ids, bootstrap=None, full_curves=True):
    """AUC(selected) vs AUC(all) uses joint FAMILY resampling of the full cohort.

    Evaluates selection effects and same-subset training effects separately.
    No claim of individual certainty, external validation or guaranteed coverage.
    """
    a = model.config
    ids, groups = prediction["eid"], prediction["group"]
    if not np.array_equal(ids.astype(str), np.asarray(truth_ids).astype(str)):
        raise ValueError("Truth IDs/order must exactly match frozen prediction roster")
    time, event = np.asarray(time, float), np.asarray(event)
    if len(time) != len(ids) or len(event) != len(ids) or not np.isfinite(time).all() or np.any(time <= 0) or not np.isin(event, [0,1]).all():
        raise ValueError("Invalid evaluation outcomes")
    if prediction["frozen_model_hash"] != model_hash(model):
        raise ValueError("Predictions do not match the supplied frozen model")
    p = prediction["probability"]
    qc = prediction["diagnostics"]["technical_qc"].astype(bool)
    selected = prediction["selected"].astype(bool)
    masks = dict(all=qc, selected=selected, rejected=qc & ~selected)
    if full_curves: masks.update(prediction["masks"])
    masks["risk_matched_random"] = risk_matched_mask(selected,qc,p["elasticnet_full"],ids,model.risk_edges,a.seed+412)
    # Exact-coverage confidence control is descriptive and outcome-blind; it is
    # NOT allowed to alter the deployable, already-frozen primary selection.
    from abm_twin_core import tie_hash
    pool=np.flatnonzero(qc)
    conf=np.abs(p['elasticnet_full']-.5)
    order=sorted(pool,key=lambda j:(-conf[j],tie_hash(str(ids[j]),a.seed+413)))
    masks['confidence_count_matched']=np.zeros(len(ids),bool)
    masks['confidence_count_matched'][np.asarray(order[:int(selected.sum())],int)]=True
    y, w = model.km.labels_weights(time, event)
    rows = []
    for subset, mask in masks.items():
        for name, prob in p.items():
            valid = mask & np.isfinite(prob)
            row = metrics(time[valid],event[valid],prob[valid],model.km)
            row.update(model=name,subset=subset,requested_N=int(mask.sum()),coverage=float(valid.mean()),
                case_coverage=float(y[valid].sum()/y.sum()) if y.sum() else np.nan,
                control_coverage=float(np.sum((w[valid]>0)&(y[valid]==0))/np.sum((w>0)&(y==0))) if np.sum((w>0)&(y==0)) else np.nan,
                IPCW_case_coverage=float(np.sum(w[valid]*y[valid])/np.sum(w*y)) if np.sum(w*y)>0 else np.nan,
                family_N=len(set(groups[valid])),horizon=a.horizon)
            rows.append(row)
    lookup = {(r['model'],r['subset']):r for r in rows}
    primary, baseline = a.primary, 'elasticnet_full'
    pairs = [('selection_effect',baseline,'selected',baseline,'all'),
             ('training_effect_same_subset',primary,'selected',baseline,'selected'),
             ('combined_vs_bulk',primary,'selected',baseline,'all'),
             ('all_cohort_training_effect',primary,'all',baseline,'all'),
             ('fallback_policy', 'fallback_policy','all',baseline,'all'),
             ('twin_library_cleaning','twin_clean','selected','twin_full','selected'),
             ('beyond_risk_composition',baseline,'selected',baseline,'risk_matched_random'),
             ('beyond_confidence_selection',baseline,'selected',baseline,'confidence_count_matched')]
    if 'elasticnet_random' in p:
        pairs.append(('pruning_specificity','elasticnet_clean','selected','elasticnet_random','selected'))
    if 'lightgbm_full' in p:
        pairs.extend([('selected_vs_full_lightgbm',primary,'selected','lightgbm_full','all'),
                      ('same_subset_vs_lightgbm',primary,'selected','lightgbm_full','selected')])
    metric_names = ['AUC_IPCW','Uno_C_horizon','Brier_IPCW','LogLoss_IPCW']
    B = a.bootstrap if bootstrap is None else int(bootstrap)
    if B < 1: raise ValueError("bootstrap must be positive")
    u, inverse = np.unique(groups,return_inverse=True)
    order = np.argsort(inverse, kind='stable')
    positions = np.split(order, np.cumsum(np.bincount(inverse))[:-1])
    rng = np.random.default_rng(a.seed+301)
    draws = {(name,m):[] for name,*_ in pairs for m in metric_names}
    needed = set((mn,sn) for _,m1,s1,m2,s2 in pairs for mn,sn in [(m1,s1),(m2,s2)])
    for b in range(B):
        ix = np.concatenate([positions[j] for j in rng.integers(0,len(u),len(u))])
        batch = {}
        for mn,sn in needed:
            jj = ix[masks[sn][ix] & np.isfinite(p[mn][ix])]
            batch[(mn,sn)] = metrics(time[jj],event[jj],p[mn][jj],model.km)
        for label,m1,s1,m2,s2 in pairs:
            for metric in metric_names:
                draws[(label,metric)].append(batch[(m1,s1)][metric]-batch[(m2,s2)][metric])
    contrasts = []
    for label,m1,s1,m2,s2 in pairs:
        for metric in metric_names:
            values = np.asarray(draws[(label,metric)]); values = values[np.isfinite(values)]
            enough = len(values) >= max(20,math.ceil(.8*B))
            bounds = np.quantile(values,[.025,.975]) if enough else [np.nan,np.nan]
            contrasts.append(dict(contrast=label,metric=metric,model=m1,subset=s1,reference=m2,reference_subset=s2,
              delta=lookup[(m1,s1)][metric]-lookup[(m2,s2)][metric],lower=float(bounds[0]),upper=float(bounds[1]),
              valid_boot=len(values),requested_boot=B,bootstrap_unit=a.grouping_unit,
              uncertainty='conditional on frozen models, selection rules and build-only censoring model',
              multiplicity='predeclared primary contrast only; other intervals exploratory'))
    c = {(r['contrast'],r['metric']):r for r in contrasts}
    sel = lookup[(primary,'selected')]
    info = sel['cases'] >= a.audit_min_cases and sel['controls'] >= a.audit_min_controls
    def supported(label,metric='AUC_IPCW',lower=True):
        rr = c[(label,metric)]
        return bool(info and rr['valid_boot'] >= max(200,math.ceil(.8*B)) and (rr['lower']>0 if lower else rr['upper']<0))
    status = dict(information_adequate=bool(info),selected_cases=sel['cases'],selected_controls=sel['controls'],
      selective_discrimination_supported=supported('combined_vs_bulk'),
      selection_only_supported=supported('selection_effect'),
      training_gain_supported=supported('training_effect_same_subset'),
      absolute_error_improvement_supported=supported('combined_vs_bulk','Brier_IPCW',False),
      beyond_risk_composition_supported=supported('beyond_risk_composition'),
      beyond_confidence_selection_supported=supported('beyond_confidence_selection'),
      clinical_release=False,external_validation=False,
      interpretation='Independent audit can validate a subset policy even without same-subset algorithm superiority; no individual guarantee')
    # Transparent test calibration diagnostics only; never used to update predictions.
    calibration = []
    for subset in ('all','selected','rejected'):
        for name in (primary,baseline):
            ok = masks[subset] & (w>0)
            yy,pp,ww = y[ok],p[name][ok],w[ok]
            row = dict(model=name,subset=subset,N=int(ok.sum()),evaluation_only=True)
            if len(set(yy)) == 2:
                from abm_twin_core import Calibrator
                cc = Calibrator().fit(pp,yy,ww)
                row.update(intercept=float(cc.model.intercept_[0]),slope=float(cc.model.coef_[0,0]),
                           observed_expected=float(np.sum(ww*yy)/np.sum(ww*pp)))
            else: row['status']='insufficient_information'
            calibration.append(row)
    strata_rows=[]
    strata={'baseline_risk_decile':np.searchsorted(model.risk_edges,p[baseline],side='right'),
            'frozen_quality_group':np.searchsorted(model.quality_edges,prediction['selector_scores'][a.selector],side='right')}
    for kind,labels in strata.items():
        for label in np.unique(labels):
            for subset,mask in [('all',qc),('selected',selected)]:
                use=mask & (labels==label)
                if not use.any(): continue
                pr=metrics(time[use],event[use],p[primary][use],model.km)
                br=metrics(time[use],event[use],p[baseline][use],model.km)
                known=use & (w>0)
                strata_rows.append(dict(stratum_type=kind,stratum=int(label),subset=subset,
                  N=int(use.sum()),cases=int(y[use].sum()),cutpoints_source='calibration_fit',
                  candidate_brier=pr['Brier_IPCW'],baseline_brier=br['Brier_IPCW'],
                  brier_gain=br['Brier_IPCW']-pr['Brier_IPCW'],
                  candidate_AUC=pr['AUC_IPCW'],baseline_AUC=br['AUC_IPCW'],
                  mean_predicted_squared_error=float(np.average(prediction['estimated_error']['absolute_error'][known],weights=w[known])) if known.any() else np.nan))
    confidence_rows=[]
    confidence=prediction['diagnostics']['pre_prediction_confidence']
    # Fixed score intervals: no outcome-driven binning and no score calibration on test.
    for lo in range(0,100,10):
        mask=qc&(confidence>=lo)&((confidence<lo+10) if lo<90 else confidence<=100)
        if mask.sum()<10:continue
        for name in (primary,baseline):
            rr=metrics(time[mask],event[mask],p[name][mask],model.km)
            rr.update(model=name,confidence_lo=lo,confidence_hi=lo+10,
                      mean_confidence=float(np.mean(confidence[mask])),
                      supported_fraction=float(np.mean(prediction['diagnostics']['supported'][mask])))
            confidence_rows.append(rr)
    return dict(metrics=rows,contrasts=contrasts,status=status,calibration=calibration,grouping_unit=a.grouping_unit,
       confidence_validation=confidence_rows,
       strata_diagnostics=strata_rows,
       roster_hash=__import__('abm_twin_core').fingerprint(ids),version=prediction['version'],
       layer=a.layer,trait=a.trait,run_id=prediction['frozen_model_hash'],evaluation_scope=prediction.get('evaluation_scope','unspecified'),
       note='Censoring weights are marginal, and inference is internal/conditional; refit-bootstrap and external validation are additional analyses')
