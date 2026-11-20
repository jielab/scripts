#!/usr/bin/env python3
"""Evaluate frozen GRID predictions and support masks on an unchanged test roster.

Only this evaluator uses test outcomes. It never refits a prediction model or
changes thresholds. Calibration slopes are diagnostic estimates, never applied
back to predictions. All public outputs are aggregate.
"""
from __future__ import annotations
import argparse
import html
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize, brentq
from scipy.special import expit
from sklearn.metrics import roc_auc_score, average_precision_score, log_loss

from parallel_runtime import ordered_map, worker_count

PRIMARY = "csx.population"

def boolean_mask(series):
    """Accept saved booleans without treating the string 'False' as true."""
    if series.isna().any():
        raise ValueError("A frozen selection mask contains missing values")
    if pd.api.types.is_bool_dtype(series):
        return series.to_numpy(bool)
    mapped = series.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False,
         "1.0": True, "0.0": False})
    if mapped.isna().any():
        raise ValueError("A frozen selection mask is not Boolean")
    return mapped.to_numpy(bool)



def interval(values, repetitions):
    values = np.asarray(values, float)
    good = values[np.isfinite(values)]
    minimum = max(2, int(np.ceil(.8 * repetitions)))
    status = ("not_requested" if repetitions == 0 else
              "available" if len(good) >= minimum else "insufficient_valid_replicates")
    low, high = (np.quantile(good, [.025, .975]) if status == "available"
                 else (np.nan, np.nan))
    return float(low), float(high), len(good), status



def paired_cluster_contrasts(y, predictions, masks, groups, repetitions, seed, workers=1):
    """One group resample for nested selected/all and paired method/reference losses.

predictions[0] is the prespecified reference model. All models must predict all rows; missing outcomes are
removed identically for every model. Masks stay fixed in every replicate.
"""
    y = np.asarray(y, float)
    predictions = np.asarray(predictions, float)
    observed = np.isfinite(y)
    m, k = len(predictions), len(masks)
    point = np.full((m, k), np.nan)
    boot_selection = np.full((repetitions, m, k), np.nan)
    boot_count = np.full_like(boot_selection, np.nan)
    if not observed.any():
        return point, boot_selection, boot_count, 0
    _, group = np.unique(np.asarray(groups, str)[observed], return_inverse=True)
    ng = int(group.max()) + 1
    loss = np.square(predictions[:, observed] - y[observed])
    counts = np.zeros((k, ng))
    sums = np.zeros((m, k, ng))
    for j, selected in enumerate(masks):
        chosen = np.asarray(selected, bool)[observed]
        counts[j] = np.bincount(group[chosen], minlength=ng)
        for a in range(m):
            sums[a, j] = np.bincount(group[chosen], weights=loss[a, chosen], minlength=ng)
    with np.errstate(divide="ignore", invalid="ignore"):
        point = sums.sum(2) / counts.sum(1)[None, :]
    if ng < 2:
        return point, boot_selection, boot_count, ng
    rng = np.random.default_rng(seed)
    def draw(weights):
        denominator = np.einsum("g,kg->k", weights, counts)
        numerator = np.einsum("g,mkg->mk", weights, sums)
        with np.errstate(divide="ignore", invalid="ignore"):
            means = numerator / denominator[None, :]
        return means - means[:, :1], means - means[:1, :]
    arguments=((np.bincount(rng.integers(0,ng,size=ng),minlength=ng).astype(float),) for _ in range(repetitions))
    count=worker_count(workers,repetitions)
    for b,(selection,reference) in enumerate(ordered_map(draw,arguments,count)):
        boot_selection[b]=selection;boot_count[b]=reference
    return point, boot_selection, boot_count, ng



def risk_standardization(y, predictions, selected, groups, cells, repetitions, seed, workers=1):
    """Fixed overlap estimand; cluster bootstrap includes both selection arms.

Risk bins are frozen in development elsewhere. Overlap cell weights are
min(n_selected,n_rejected), calculated on the observed evaluation roster using
only ancestry, predictions, selection, and outcome availability, never Y values.
Weights/cells stay fixed during bootstrap. Replicates missing either arm in a
required cell are reported as invalid, not silently substituted or reweighted.
"""
    m = len(predictions)
    empty = dict(risk_standardized_selected_loss=np.nan,
                 risk_standardized_rejected_loss=np.nan,
                 risk_standardized_difference=np.nan,
                 risk_difference_low=np.nan, risk_difference_high=np.nan,
                 risk_valid_bootstrap_replicates=0, risk_ci_status="not_estimable",
                 overlap_selected_n=0, overlap_rejected_n=0,
                 overlap_selected_fraction=0., overlap_rejected_fraction=0.,
                 overlap_cells=0, risk_status="no_common_ancestry_risk_cells")
    observed = np.isfinite(y)
    if cells is None:
        return [dict(empty, risk_status="calibration_reference_unavailable") for _ in range(m)]
    if not observed.any():
        return [dict(empty, risk_status="no_observed_outcomes") for _ in range(m)]
    y = np.asarray(y, float)[observed]
    pred = np.asarray(predictions, float)[:, observed]
    sel = np.asarray(selected, bool)[observed]
    _, g = np.unique(np.asarray(groups, str)[observed], return_inverse=True)
    _, c = np.unique(np.asarray(cells, str)[observed], return_inverse=True)
    ng, nc = int(g.max()) + 1, int(c.max()) + 1
    # 0 = rejected, 1 = selected. Collapse repeated family/cell/arm entries.
    bucket = (g * nc + c) * 2 + sel.astype(int)
    key, inverse = np.unique(bucket, return_inverse=True)
    bucket_group = key // (2 * nc)
    bucket_cell_arm = key % (2 * nc)
    bucket_n = np.bincount(inverse).astype(float)
    bucket_loss = np.asarray([np.bincount(inverse, weights=np.square(p-y)) for p in pred])
    n = np.bincount(bucket_cell_arm, weights=bucket_n, minlength=2*nc).reshape(nc, 2)
    common = (n[:, 0] > 0) & (n[:, 1] > 0)
    if not common.any():
        return [empty.copy() for _ in range(m)]
    weights = np.minimum(n[common, 0], n[common, 1])
    weights = weights / weights.sum()
    means = []
    for a in range(m):
        sums = np.bincount(bucket_cell_arm, weights=bucket_loss[a], minlength=2*nc).reshape(nc, 2)
        means.append(weights @ (sums[common] / n[common]))
    means = np.asarray(means)
    boot = np.full((repetitions, m), np.nan)
    rng = np.random.default_rng(seed)
    if ng >= 2:
        def draw(gw):
            values=np.full(m,np.nan)
            bw = gw[bucket_group]
            bn = np.bincount(bucket_cell_arm, weights=bucket_n*bw, minlength=2*nc).reshape(nc, 2)
            if np.any(bn[common] <= 0):return values
            for a in range(m):
                bs = np.bincount(bucket_cell_arm, weights=bucket_loss[a]*bw, minlength=2*nc).reshape(nc, 2)
                arm_means = weights @ (bs[common]/bn[common])
                values[a] = arm_means[1] - arm_means[0]
            return values
        arguments=((np.bincount(rng.integers(0,ng,size=ng),minlength=ng),) for _ in range(repetitions))
        for b,values in enumerate(ordered_map(draw,arguments,worker_count(workers,repetitions))):boot[b]=values
    result = []
    for a in range(m):
        low, high, valid, status = interval(boot[:, a], repetitions)
        result.append(dict(risk_standardized_selected_loss=means[a, 1],
            risk_standardized_rejected_loss=means[a, 0],
            risk_standardized_difference=means[a, 1]-means[a, 0],
            risk_difference_low=low, risk_difference_high=high,
            risk_valid_bootstrap_replicates=valid, risk_ci_status=status,
            overlap_selected_n=int(n[common, 1].sum()),
            overlap_rejected_n=int(n[common, 0].sum()),
            overlap_selected_fraction=float(n[common, 1].sum()/sel.sum()) if sel.any() else 0.,
            overlap_rejected_fraction=float(n[common, 0].sum()/(~sel).sum()) if (~sel).any() else 0.,
            overlap_cells=int(common.sum()), risk_status="descriptive_fixed_overlap_standardization"))
    return result


def read_table(path):
    return pd.read_csv(path,sep="\t",dtype={"id":str,"group_id":str,"ancestry":str})


def aligned(frame, ids, label):
    if "id" not in frame or frame.id.isna().any() or frame.id.duplicated().any():
        raise ValueError(label+": unique nonmissing IDs required")
    if not set(ids)<=set(frame.id):
        raise ValueError(label+": missing roster IDs; no intersection filtering")
    return frame.set_index("id").loc[list(ids)].reset_index()


def binary_calibration(y,p):
    """Diagnostic intercept/slope; never used to modify frozen predictions."""
    if len(np.unique(y))<2:
        return dict(calibration_intercept=np.nan,calibration_slope=np.nan,
                    calibration_in_the_large=np.nan,calibration_status="one_class")
    z=np.log(np.clip(p,1e-8,1-1e-8)/np.clip(1-p,1e-8,1-1e-8))
    citl=float(brentq(lambda a:expit(a+z).mean()-y.mean(),-60.,60.))
    if np.std(z)<1e-12:
        return dict(calibration_intercept=np.nan,calibration_slope=np.nan,
                    calibration_in_the_large=citl,calibration_status="constant_prediction")
    x=np.column_stack([np.ones(len(z)),z])
    def objective(beta):
        eta=x@beta
        return float(np.mean(np.logaddexp(0,eta)-y*eta)),x.T@(expit(eta)-y)/len(y)
    fit=minimize(objective,[0.,1.],jac=True,method="BFGS",options={"maxiter":200,"gtol":1e-6})
    valid=bool(fit.success and np.isfinite(fit.x).all())
    return dict(calibration_intercept=float(fit.x[0]) if valid else np.nan,
        calibration_slope=float(fit.x[1]) if valid else np.nan,
        calibration_in_the_large=citl,calibration_status="diagnostic_fit" if valid else "fit_unavailable")


def metrics(y,p,binary,fit_calibration=False,covariate=None):
    y,p=np.asarray(y,float),np.asarray(p,float)
    valid=np.isfinite(y)&np.isfinite(p)
    reference=None if covariate is None else np.asarray(covariate,float)
    if reference is not None:
        if reference.shape!=y.shape or not np.isfinite(reference).all():
            raise ValueError("Covariate baseline must predict the same people completely")
        reference=reference[valid]
    y,p=y[valid],p[valid]
    out=dict(n_evaluable=len(y),primary_loss=np.nan,mae=np.nan,mean_prediction=np.nan,
             outcome_mean=np.nan,mean_prediction_minus_outcome=np.nan,
             covariate_reference_status="unavailable" if reference is None else "no_observed_outcomes")
    if not len(y):
        return out
    mse=float(np.mean((y-p)**2))
    out.update(primary_loss=mse,mae=float(np.mean(abs(y-p))),mean_prediction=float(p.mean()),
               outcome_mean=float(y.mean()),mean_prediction_minus_outcome=float(p.mean()-y.mean()),
               outcome_variance=float(np.var(y)),prediction_variance=float(np.var(p)),
               prediction_to_outcome_sd_ratio=float(np.std(p)/np.std(y)) if np.std(y)>0 else np.nan)
    if reference is not None:
        reference_loss=float(np.mean((y-reference)**2))
        out.update(covariate_reference_status="available_frozen_train_model",
                   covariate_primary_loss=reference_loss)
        if binary:
            reference_logloss=float(log_loss(y,np.clip(reference,1e-12,1-1e-12),labels=[0,1]))
            out.update(covariate_logloss=reference_logloss,
                delta_logloss_vs_covariate=float(log_loss(y,np.clip(p,1e-12,1-1e-12),labels=[0,1]))-reference_logloss,
                brier_difference_vs_covariate=mse-reference_loss,
                brier_skill_vs_covariate=1-mse/reference_loss if reference_loss>0 else np.nan)
        else:
            variance=float(np.var(y))
            out.update(partial_r2_vs_covariate=1-mse/reference_loss if reference_loss>0 else np.nan,
                       incremental_r2_vs_covariate=(reference_loss-mse)/variance if variance>0 else np.nan)
    if binary:
        both=len(np.unique(y))==2
        reference=float(y.mean()*(1-y.mean()))
        out.update(n_cases=int(y.sum()),event_fraction=float(y.mean()),brier=mse,
            auroc=float(roc_auc_score(y,p)) if both else np.nan,
            auprc=float(average_precision_score(y,p)) if both else np.nan,
            logloss=float(log_loss(y,np.clip(p,1e-12,1-1e-12),labels=[0,1])),
            brier_skill_vs_subset_prevalence=1-mse/reference if reference else np.nan)
        if fit_calibration:
            out.update(binary_calibration(y,p))
        else:
            out.update(calibration_intercept=np.nan,calibration_slope=np.nan,
                       calibration_in_the_large=np.nan,calibration_status="slope_diagnostic_full_cohort_only")
    else:
        variance=float(np.var(y));pv=float(np.var(p))
        slope=float(np.mean((p-p.mean())*(y-y.mean()))/pv) if pv else np.nan
        out.update(rmse=float(np.sqrt(mse)),r2=1-mse/variance if variance else np.nan,
            correlation=float(np.corrcoef(y,p)[0,1]) if variance and pv else np.nan,
            calibration_slope=slope,calibration_intercept=float(y.mean()-slope*p.mean()) if pv else np.nan,
            calibration_status="diagnostic_linear_fit" if pv else "constant_prediction")
    return out


class FrozenDiscrimination:
    """Cache score/tie ordering, then compute weighted AUROC and average precision.

    Family bootstrap weights change; prediction order and tie definitions do not.
    This avoids sorting N predictions again in every bootstrap replicate.
    """
    def __init__(self,y,prediction):
        y,prediction=np.asarray(y,float),np.asarray(prediction,float)
        self.order=np.argsort(prediction,kind="stable")
        ordered=prediction[self.order]
        self.starts=np.r_[0,1+np.flatnonzero(ordered[1:]!=ordered[:-1])] if len(y) else np.array([],int)
        self.positive=y[self.order]

    def scores(self,weights):
        if not len(self.order):
            return np.array([np.nan,np.nan])
        weight=np.asarray(weights,float)[self.order]
        positive=np.add.reduceat(weight*self.positive,self.starts)
        negative=np.add.reduceat(weight*(1-self.positive),self.starts)
        npos,nneg=positive.sum(),negative.sum()
        if npos<=0 or nneg<=0:
            return np.array([np.nan,np.nan])
        below=np.cumsum(negative)-negative
        auc=np.sum(positive*(below+.5*negative))/(npos*nneg)
        tp=np.cumsum(positive[::-1]);total=np.cumsum((positive+negative)[::-1])
        precision=np.divide(tp,total,out=np.zeros_like(tp),where=total>0)
        ap=np.sum(positive[::-1]*precision)/npos
        return np.asarray([auc,ap])


def paired_discrimination(y,predictions,selections,groups,ancestry,repetitions,seed,workers=1):
    """Predeclared binary contrasts on all/support20 and all/ancestry strata.

    The selected mask and all predictors are frozen. Each stratum resamples its
    declared family groups, applying identical multiplicities to every method.
    Average precision is the noninterpolated PR summary used by sklearn.
    """
    y=np.asarray(y,float);groups=np.asarray(groups,str);ancestry=np.asarray(ancestry,str)
    comparisons=("disco","abm","prscs.best")+tuple(
        name for name in ("csx.input","block.glmnet","block.lightgbm","abm.deep") if name in predictions)
    masks=("all","selected_support_0p20")
    strata=[("ALL",np.ones(len(y),bool))]+[(a,ancestry==a) for a in sorted(set(ancestry))]
    rows=[]
    for stratum_index,(stratum,in_stratum) in enumerate(strata):
        observed=in_stratum&np.isfinite(y)
        _,group=np.unique(groups[observed],return_inverse=True)
        ng=int(group.max())+1 if len(group) else 0
        for label in masks:
            missing_mask=label not in selections
            selected=(np.zeros(len(y),bool) if missing_mask else np.asarray(selections[label],bool))
            take=selected[observed]
            yy=y[observed][take];gg=group[take]
            ncase=int(yy.sum());ncontrol=len(yy)-ncase
            selected_groups=len(np.unique(gg))
            active=[PRIMARY]+[name for name in comparisons if name in predictions]
            scorers={name:FrozenDiscrimination(yy,np.asarray(predictions[name])[observed][take]) for name in active}
            points={name:scorer.scores(np.ones(len(yy))) for name,scorer in scorers.items()}
            draws={name:np.full((repetitions,2),np.nan) for name in active if name!=PRIMARY}
            if ncase and ncontrol and selected_groups>=2 and ng>=2:
                # Same seed within a stratum also pairs all/support20 resampling.
                rng=np.random.default_rng(seed+stratum_index)
                def draw(frequency):
                    weight=frequency[gg]
                    reference=scorers[PRIMARY].scores(weight)
                    return {name:scorers[name].scores(weight)-reference for name in draws}
                arguments=((np.bincount(rng.integers(0,ng,size=ng),minlength=ng),) for _ in range(repetitions))
                for b,values in enumerate(ordered_map(draw,arguments,worker_count(workers,repetitions))):
                    for name in draws:draws[name][b]=values[name]
            for name in comparisons:
                available=name in scorers
                point=points[name] if available else np.array([np.nan,np.nan])
                delta=point-points[PRIMARY]
                result=dict(method=name,reference=PRIMARY,ancestry=stratum,selection=label,
                    n_roster=int(in_stratum.sum()),n_evaluable=len(yy),n_cases=ncase,n_controls=ncontrol,
                    n_selected_groups=selected_groups,n_bootstrap_groups=ng,
                    method_auroc=point[0],reference_auroc=points[PRIMARY][0],delta_auroc=delta[0],
                    method_average_precision=point[1],reference_average_precision=points[PRIMARY][1],
                    delta_average_precision=delta[1],bootstrap_replicates=repetitions,
                    bootstrap_unit="declared_family_group",
                    contrast_status="mask_unavailable" if missing_mask else "method_unavailable" if not available
                        else "no_observed_outcomes" if not len(yy) else "one_class" if not ncase or not ncontrol
                        else "fewer_than_two_selected_groups" if selected_groups<2
                        else "low_case_or_control_count_lt20" if min(ncase,ncontrol)<20 else "available",
                    inference="Pointwise; conditional on frozen models and masks; no multiplicity correction")
                for j,metric in enumerate(("auroc","average_precision")):
                    values=draws[name][:,j] if available else np.full(repetitions,np.nan)
                    low,high,valid,status=interval(values,repetitions)
                    result.update({f"delta_{metric}_low":low,f"delta_{metric}_high":high,
                                   f"{metric}_valid_replicates":valid,f"{metric}_ci_status":status})
                rows.append(result)
    return pd.DataFrame(rows)


def evaluate_predictions(cohort,outcomes,benchmark,abm,output_dir,*,trait,bootstrap=200,seed=20261115,workers=1):
    required={"id","group_id","ancestry","role"}
    if not required<=set(cohort) or cohort[list(required)].isna().any().any() or cohort.id.duplicated().any():
        raise ValueError("Invalid evaluation roster")
    if cohort.groupby("group_id").role.nunique().max()!=1:
        raise ValueError("A group crosses development and test roles")
    test=cohort.loc[cohort.role.eq("test")].reset_index(drop=True)
    cal=cohort.loc[cohort.role.eq("calibration")].reset_index(drop=True)
    if not len(test) or not len(cal):
        raise ValueError("Separate test and calibration rosters required")
    bt=aligned(benchmark,test.id,"benchmark test predictions")
    at=aligned(abm,test.id,"ABM test predictions")
    if "y" in benchmark or "y" in abm:
        raise ValueError("Frozen prediction snapshots must not include outcomes")
    for frame,label in [(bt,"benchmark"),(at,"ABM")]:
        if "role" in frame and not frame.role.eq("test").all():
            raise ValueError(label+" test roles differ")
        if "group_id" in frame and not np.array_equal(frame.group_id.astype(str),test.group_id.astype(str)):
            raise ValueError(label+" test family mapping differs")
    # Test labels enter only here, after frozen models and selections are loaded.
    truth=aligned(outcomes,test.id,"test outcomes")
    y=pd.to_numeric(truth.y,errors="raise").to_numpy(float)
    binary=trait=="t2dm"
    if np.isinf(y).any() or (binary and not np.isin(y[np.isfinite(y)],[0.,1.]).all()):
        raise ValueError("Invalid test outcomes")
    columns={}
    for source in (bt,at):
        for column in source:
            if not column.endswith(".prediction"):
                continue
            name=column[:-11]
            if name in columns:
                raise ValueError("Duplicate prediction method: "+name)
            p=pd.to_numeric(source[column],errors="raise").to_numpy(float)
            if not np.isfinite(p).all() or (binary and ((p<0)|(p>1)).any()):
                raise ValueError("Incomplete/out-of-range frozen predictions: "+name)
            columns[name]=p
    if PRIMARY not in columns or "abm" not in columns:
        raise ValueError("Primary CSx and ABM predictions are required")
    names=[PRIMARY]+sorted(x for x in columns if x!=PRIMARY)
    predictions=np.asarray([columns[n] for n in names])
    covariate=columns.get("covariate")
    mask_columns=sorted(c for c in at if c.startswith("selected_"))
    if not mask_columns:
        raise ValueError("No calibration-frozen selection masks")
    selections=[("all",np.ones(len(test),bool))]+[(c,boolean_mask(at[c])) for c in mask_columns]
    masks=[x[1] for x in selections]
    means,selection_boot,reference_boot,ng=paired_cluster_contrasts(y,predictions,masks,test.group_id,bootstrap,seed,workers=workers)
    global_index=names.index("abm.global_residual") if "abm.global_residual" in names else None
    # A common reference is prespecified, never the test-best model.
    bc=aligned(benchmark,cal.id,"benchmark calibration predictions")
    cp=pd.to_numeric(bc[PRIMARY+".prediction"],errors="raise").to_numpy(float)
    if not np.isfinite(cp).all():
        raise ValueError("Incomplete calibration risk reference; no test-derived replacement")
    edges=np.unique(np.quantile(cp,np.arange(.1,1.,.1)))
    bins=np.searchsorted(edges,columns[PRIMARY],side="right")
    cells=np.asarray([json.dumps([str(a),int(b)],separators=(",",":")) for a,b in zip(test.ancestry,bins)])
    rows,ancestry_rows=[],[]
    observed=np.isfinite(y)
    total_cases=float(y[observed].sum()) if binary else np.nan
    total_controls=float(observed.sum()-total_cases) if binary else np.nan
    error_head_audit={}
    error_cols=[c for c in at if c.endswith(".expected_loss")]
    reference_error=np.square(y-columns[PRIMARY])
    if {"risk_error.expected_loss","support_error.expected_loss"}<=set(error_cols):
        # These evaluate accuracy of the expected-error estimates, not accuracy
        # of disease predictions. Targets and both heads are fixed in training.
        pe=np.asarray([at["risk_error.expected_loss"],at["support_error.expected_loss"]],float)
        error_mean,_,error_boot,_=paired_cluster_contrasts(reference_error,pe,[np.ones(len(test),bool)],test.group_id,bootstrap,seed,workers=workers)
        low,high,valid,status=interval(error_boot[:,1,0],bootstrap)
        error_head_audit=dict(risk_error_head_mse=float(error_mean[0,0]),support_error_head_mse=float(error_mean[1,0]),
            support_minus_risk_error_head_mse=float(error_mean[1,0]-error_mean[0,0]),
            difference_low=low,difference_high=high,valid_bootstrap_replicates=valid,ci_status=status)
    for j,(label,selected) in enumerate(selections):
        comparison=risk_standardization(y,predictions,selected,test.group_id,cells,bootstrap,seed,workers=workers)
        if label=="all":
            for value in comparison:
                value["risk_status"]="not_applicable_all_has_no_rejected_group"
        selected_cases=float(y[selected&observed].sum()) if binary else np.nan
        selected_controls=float((selected&observed).sum()-selected_cases) if binary else np.nan
        for a,name in enumerate(names):
            sl,sh,sn,ss=interval(selection_boot[:,a,j],bootstrap)
            pl,ph,pn,ps=interval(reference_boot[:,a,j],bootstrap)
            global_contrast={}
            if global_index is not None:
                gl,gh,gn,gs=interval(reference_boot[:,a,j]-reference_boot[:,global_index,j],bootstrap)
                global_contrast=dict(method_minus_global_residual_loss=means[a,j]-means[global_index,j],
                    global_difference_low=gl,global_difference_high=gh,global_valid_replicates=gn,global_ci_status=gs)
            perf=metrics(y[selected],predictions[a,selected],binary,fit_calibration=label=="all",
                         covariate=None if covariate is None else covariate[selected])
            rows.append(dict(trait=trait,method=name,selection=label,primary_reference=PRIMARY,
                selector_source="supervised_deep_embedding_exploratory" if label.startswith("selected_deep_")
                    else "CPU_frozen_block_support" if label.startswith("selected_support_")
                    else "OOF_expected_error_model" if label.startswith("selected_") else "all_test",
                n_roster=len(test),n_selected_roster=int(selected.sum()),actual_coverage=float(selected.mean()),
                n_missing_y=int((selected&~observed).sum()),loss_name="brier" if binary else "mse",
                case_retention=selected_cases/total_cases if total_cases>0 else np.nan,
                control_retention=selected_controls/total_controls if total_controls>0 else np.nan,
                selected_minus_all_loss=means[a,j]-means[a,0],selection_low=sl,selection_high=sh,
                selection_valid_replicates=sn,selection_ci_status=ss,
                method_minus_primary_loss=means[a,j]-means[0,j],paired_low=pl,paired_high=ph,
                paired_valid_replicates=pn,paired_ci_status=ps,n_bootstrap_groups=ng,
                bootstrap_unit="declared_family_group",bootstrap_replicates=bootstrap,
                risk_reference="calibration_frozen_csx_population_x_ancestry",**perf,**comparison[a],**global_contrast))
            for ancestry,index in test.groupby("ancestry",sort=True).groups.items():
                group=np.zeros(len(test),bool);group[np.asarray(list(index),int)]=True
                take=group&selected
                ancestry_rows.append(dict(trait=trait,method=name,selection=label,ancestry=str(ancestry),
                    n_ancestry_roster=int(group.sum()),n_selected_roster=int(take.sum()),
                    ancestry_coverage=float(take.sum()/group.sum()),n_missing_y=int((take&~observed).sum()),
                    **metrics(y[take],predictions[a,take],binary,
                              covariate=None if covariate is None else covariate[take])))
    performance,ancestry=pd.DataFrame(rows),pd.DataFrame(ancestry_rows)
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
    performance.to_csv(out/f"{trait}.performance.tsv",sep="\t",index=False)
    ancestry.to_csv(out/f"{trait}.ancestry.tsv",sep="\t",index=False)
    details=dict(trait=trait,primary=PRIMARY,test_N=len(test),risk_calibration_N=len(cal),
        risk_bin_edges=edges.tolist(),error_head_audit=error_head_audit,
        inference="Exploratory pointwise family-group bootstrap conditional on frozen models and masks",
        risk_estimand="Fixed overlap weights min(selected_N,rejected_N) in ancestry x frozen primary-risk bins",
        calibration="Test intercepts/slopes are diagnostics only; no recalibration applied",
        covariate_reference="Frozen TRAIN-fitted covariate model evaluated on the same people; no test covariate regression",
        continuous_increment_definition="partial_R2=1-SSE_method/SSE_covariate; incremental_R2=(SSE_covariate-SSE_method)/SST; both may be negative")
    sex_diagnostics=[]
    if "sex" in test:
        # Grouping numeric NaNs with dropna=False can make pandas construct an
        # invalid null Categorical. Keep missing sex explicit in diagnostics;
        # model fitting already uses its frozen TRAIN imputation policy.
        sex_labels=test["sex"].astype("string").fillna("missing")
        for sex,index in test.groupby(sex_labels,sort=True).groups.items():
            rows=np.asarray(list(index),int)
            for a,name in enumerate(names):
                sex_diagnostics.append(dict(method=name,sex=str(sex),n_roster=len(rows),
                                            **metrics(y[rows],predictions[a,rows],binary,
                                                      covariate=None if covariate is None else covariate[rows])))
    details["sex_diagnostics"]=sex_diagnostics
    discrimination=None
    if binary:
        discrimination=paired_discrimination(y,columns,dict(selections),test.group_id,test.ancestry,bootstrap,seed,workers=workers)
        discrimination.insert(0,"trait",trait)
        discrimination.to_csv(out/f"{trait}.paired_discrimination.tsv",sep="\t",index=False)
        details["paired_discrimination"]=dict(file=f"{trait}.paired_discrimination.tsv",
            contrasts=["disco", "abm", "prscs.best"],reference=PRIMARY,
            selections=["all", "selected_support_0p20"],strata="ALL and each ancestry",
            ordering="Fixed prediction/tie ordering; family multiplicity weights",
            contrast_direction="Positive AUROC/average-precision differences favor the comparison method",
            average_precision="Noninterpolated average precision; not trapezoidal PR area",
            inference="Pointwise conditional on frozen models and masks; no multiplicity adjustment; flag cells with fewer than20 cases or controls")
    (out/f"{trait}.evaluation.json").write_text(json.dumps(details,indent=2,allow_nan=True)+"\n")
    explanatory=("Masks, models and the primary CSx comparison are frozen before test evaluation. "
        "Selected-minus-all is an association between selection and error; the all group includes selected people. "
        "Method-minus-primary compares methods on identical people. Negative loss differences favor the first term. "
        "The common-risk contrast uses calibration-only prediction deciles crossed with ancestry; read the overlap "
        "fractions and valid bootstrap counts. It is descriptive and does not remove all confounding. "
        "Support counts describe a specified genetic-score space, not biological twins or guaranteed individual accuracy. "
        "Family ESS summarizes weight/count concentration, not verified statistical independence. "
        "Confidence intervals are pointwise, conditional on frozen models; they exclude GWAS and training uncertainty. "
        "Raw neighbour mean is explicitly an ablation. Nominal full coverage can retain only supported people. "
        "Missing outcomes are excluded identically for all methods. Binary slope diagnostics are fitted only for full-test rows; "
        "subset prevalence Brier skill is descriptive. Covariate-relative diagnostics instead compare with the frozen "
        "TRAIN-fitted covariate model on exactly the same people. Partial R² is 1−SSE(method)/SSE(covariate); "
        "incremental R² divides that SSE improvement by the outcome total sum of squares. Both can be negative. "
        "For binary outcomes, paired discrimination comparisons are predeclared; confidence intervals are pointwise, "
        "uncorrected for multiplicity, and conditional on the fixed models. Positive AUROC/average-precision differences "
        "favor the comparison method. Low case/control counts are flagged.")
    body=f"<h1>{html.escape(trait)}: frozen GRID evaluation</h1><p>{html.escape(explanatory)}</p>"
    body+="<h2>Expected-error model audit</h2><pre>"+html.escape(json.dumps(error_head_audit,indent=2))+"</pre>"
    if sex_diagnostics:
        body+="<h2>Full-test mean error and variance by sex</h2>"+pd.DataFrame(sex_diagnostics).to_html(index=False,border=0)
    body+="<h2>Same-person comparisons</h2>"+performance.to_html(index=False,border=0)
    body+="<h2>Ancestry coverage and performance</h2>"+ancestry.to_html(index=False,border=0)
    if discrimination is not None:
        body+="<h2>Predeclared paired discrimination contrasts</h2>"+discrimination.to_html(index=False,border=0)
    (out/f"{trait}.report.html").write_text("<!doctype html><html><head><meta charset='utf-8'>"
        "<style>body{font:14px system-ui;margin:28px}p{max-width:1100px;line-height:1.5}"
        "table{border-collapse:collapse;font-size:11px}td,th{padding:5px;border:1px solid #ddd}</style>"
        "</head><body>"+body+"</body></html>")
    return performance,ancestry,details


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("cohort","outcomes","predictions","abm-predictions","trait","output-dir"):
        parser.add_argument("--"+name,required=True)
    parser.add_argument("--bootstrap",type=int,default=200);parser.add_argument("--seed",type=int,default=20261115)
    args=parser.parse_args(argv)
    if args.bootstrap<0:
        parser.error("bootstrap must be nonnegative")
    evaluate_predictions(read_table(args.cohort),read_table(args.outcomes),read_table(args.predictions),
        read_table(args.abm_predictions),args.output_dir,trait=args.trait,bootstrap=args.bootstrap,seed=args.seed)


if __name__=="__main__":
    main()


