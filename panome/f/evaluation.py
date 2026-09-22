"""Censor-aware metrics and paired uncertainty on one locked test set."""
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored, cumulative_dynamic_auc, brier_score, concordance_index_ipcw
from sksurv.nonparametric import kaplan_meier_estimator, CensoringDistributionEstimator
from common import dump

def survival_y(p):
    return Surv.from_arrays(p.event.to_numpy(bool), p.time.to_numpy(float))

def concordance(p, score):
    return float(concordance_index_censored(p.event.astype(bool), p.time, score)[0])

def primary_metric(p, score, kind):
    return concordance(p, score) if kind == "survival" else float(r2_score(p.target, score))

def comparable_y(train, test, horizon):
    yt = survival_y(test)
    ytr = survival_y(train)
    support = float(np.max(ytr["time"]))
    tau = min(float(np.quantile(ytr["time"], .99)), support) - 1e-6
    if not yt["time"].min() < horizon < min(tau, yt["time"].max()):
        raise ValueError("Horizon outside training/test follow-up support")
    g = CensoringDistributionEstimator().fit(ytr)
    if g.predict_proba([horizon])[0] <= .01:
        raise ValueError("Training censoring survival at horizon <= 0.01")
    yt["event"] &= yt["time"] <= tau
    yt["time"] = np.minimum(yt["time"], tau)
    return ytr, yt

def evaluate(train, test, score, survival_at, kind, horizons, name):
    rows, calibration, limitations = [], [], []
    if kind == "quantitative":
        for metric, value in [("R2", r2_score(test.target, score)),
                               ("RMSE", np.sqrt(mean_squared_error(test.target, score))),
                               ("MAE", mean_absolute_error(test.target, score))]:
            rows.append(dict(model=name, metric=metric, horizon=0, value=float(value), n=len(test)))
        groups = pd.qcut(pd.Series(score), 5, labels=False, duplicates="drop").to_numpy()
        for group in np.unique(groups[np.isfinite(groups)]):
            ix = groups == group
            calibration.append(dict(model=name, horizon=0, group=int(group)+1,
                n=int(ix.sum()), predicted=float(score[ix].mean()),
                observed=float(test.target.to_numpy()[ix].mean())))
        return rows, calibration, limitations, {}
    rows.append(dict(model=name, metric="Harrell_C", horizon=0,
                     value=concordance(test, score), n=len(test), events=int(test.event.sum())))
    risks = {}
    for horizon in horizons:
        try:
            ytr, yt = comparable_y(train, test, horizon)
            survival = np.asarray(survival_at(horizon), float)
            if not np.isfinite(survival).all() or np.any((survival<0) | (survival>1)):
                raise ValueError("Invalid survival probabilities")
            risks[horizon] = 1-survival
            auc = float(cumulative_dynamic_auc(ytr, yt, score, [horizon])[0][0])
            bs = float(brier_score(ytr, yt, survival[:, None], [horizon])[1][0])
            uno = float(concordance_index_ipcw(ytr, yt, score, tau=horizon)[0])
            for metric, value in [("AUC_IPCW", auc), ("Brier_IPCW", bs), ("Uno_C", uno)]:
                rows.append(dict(model=name, metric=metric, horizon=horizon, value=value,
                                 n=len(test), events=int((yt["event"] & (yt["time"]<=horizon)).sum())))
            groups = pd.qcut(pd.Series(1-survival), 5, labels=False, duplicates="drop").to_numpy()
            for group in np.unique(groups[np.isfinite(groups)]):
                ix = groups == group
                t, s = kaplan_meier_estimator(yt["event"][ix], yt["time"][ix])
                hit = np.searchsorted(t, horizon, side="right")-1
                obs = 1-s[hit] if hit >= 0 else 0
                if np.max(yt["time"][ix]) < horizon:
                    obs = np.nan
                calibration.append(dict(model=name, horizon=horizon, group=int(group)+1,
                    n=int(ix.sum()), predicted=float(np.mean(1-survival[ix])),
                    observed=obs, at_risk=int((yt["time"][ix]>=horizon).sum())))
        except (ValueError, ArithmeticError) as exc:
            limitations.append(dict(model=name, horizon=horizon, reason=str(exc)))
    return rows, calibration, limitations, risks

def paired_uncertainty(p, predictions, kind, groups, count, seed, out):
    rng = np.random.default_rng(seed)
    comparisons = [(m, "clinical") for m in predictions if m != "clinical"]
    comparisons += [("panome_varying", "clinical_ae"),
                    ("panome_neighborhood", "clinical_ae"),
                    ("panome_transformer", "transformer_column"),
                    ("panome_transformer", "transformer_random")]
    comparisons = list(dict.fromkeys((m, b) for m,b in comparisons if m in predictions and b in predictions))
    metrics = {m: primary_metric(p, values, kind) for m,values in predictions.items()}
    boot, contrasts = [], []
    group_rows = None
    if groups is not None:
        group_rows = [np.flatnonzero(np.asarray(groups)==g) for g in np.unique(groups)]
    for repeat in range(count):
        if group_rows is None:
            ix = rng.integers(0, len(p), len(p))
        else:
            draw = rng.integers(0, len(group_rows), len(group_rows))
            ix = np.concatenate([group_rows[j] for j in draw])
        values = {}
        for model, pred in predictions.items():
            try:
                value = primary_metric(p.iloc[ix], pred[ix], kind)
                if not np.isfinite(value):
                    continue
                values[model] = value
                boot.append(dict(replicate=repeat+1, model=model, value=values[model]))
            except ValueError:
                continue
        for model, reference in comparisons:
            if model in values and reference in values:
                contrasts.append(dict(replicate=repeat+1, model=model, reference=reference,
                    delta=values[model]-values[reference]))
    b = pd.DataFrame(boot, columns=["replicate","model","value"])
    d = pd.DataFrame(contrasts, columns=["replicate","model","reference","delta"])
    b.to_csv(out/"paired_bootstrap.csv", index=False)
    d.to_csv(out/"paired_bootstrap_contrasts.csv", index=False)
    result = []
    for model, reference in comparisons:
        values = d.loc[(d.model==model)&(d.reference==reference), "delta"].to_numpy()
        low, high = np.quantile(values, [.025,.975]) if len(values)>=10 else (np.nan,np.nan)
        result.append(dict(model=model, reference=reference, metric="Harrell_C" if kind=="survival" else "R2",
            delta=metrics[model]-metrics[reference], lower=low, upper=high,
            successful_replicates=len(values), uncertainty="conditional_on_fitted_models"))
    pd.DataFrame(result, columns=["model","reference","metric","delta","lower","upper",
                                  "successful_replicates","uncertainty"]).to_csv(out/"paired_contrasts.csv", index=False)

def landmarks(train, test, predictions, a, out):
    rows = []
    if a.outcome_type == "survival":
        for lag in a.landmarks:
            keep = test.time.to_numpy() > lag
            sub = test.loc[keep].copy()
            sub["time"] -= lag
            if len(sub)<20 or sub.event.sum()<5:
                continue
            for model, pred in predictions.items():
                try:
                    value = concordance(sub, pred[keep])
                    rows.append(dict(model=model, landmark=lag, n=len(sub),
                        events=int(sub.event.sum()), Harrell_C=value))
                except ValueError:
                    continue
    pd.DataFrame(rows, columns=["model","landmark","n","events","Harrell_C"]).to_csv(
        out/"landmark_sensitivity.csv", index=False)
