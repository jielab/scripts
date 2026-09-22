"""Fixed-horizon IPCW learning. Death is censored: this estimates net risk.

Censoring distribution is estimated in the reference-building sample only.
Weights require independent censoring; no clipping conceals lack of support.
"""
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score


class CensoringKM:
    def fit(self, p, horizon, min_g=.05):
        t = p.time.to_numpy(float)
        e = p.event.to_numpy(bool)
        self.horizon = float(horizon)
        self.times, inverse, counts = np.unique(t, return_inverse=True, return_counts=True)
        at_risk = len(t) - np.r_[0, np.cumsum(counts)[:-1]]
        # With tied event/censor times, observed failures leave before censoring.
        failures = np.bincount(inverse, weights=e, minlength=len(counts))
        censored = counts - failures
        denom = at_risk - failures
        step = np.divide(censored, denom, out=np.zeros_like(failures), where=denom > 0)
        self.survival = np.cumprod(1 - step)
        self.min_g = min_g
        if np.sum(t >= horizon) < 10 or self.at([horizon])[0] < min_g:
            raise ValueError("Horizon lacks censoring support in build sample; shorten --horizon")
        return self

    def at(self, t, left=False):
        ix = np.searchsorted(self.times, np.asarray(t), side="left" if left else "right") - 1
        return np.where(ix >= 0, self.survival[np.maximum(ix, 0)], 1.)

    def labels_weights(self, p):
        t, e = p.time.to_numpy(float), p.event.to_numpy(bool)
        case = e & (t <= self.horizon)
        control = t > self.horizon
        g = np.ones(len(t))
        g[case] = self.at(t[case], left=True)
        g[control] = self.at(np.full(control.sum(), self.horizon))
        known = case | control
        if np.any(g[known] < self.min_g):
            raise ValueError("IPCW unsupported: censor survival below --min-censor-survival")
        w = np.zeros(len(t))
        w[known] = 1 / g[known]
        return case.astype(int), w


def loss(y, p):
    p = np.clip(p, 1e-6, 1-1e-6)
    return -(y*np.log(p)+(1-y)*np.log1p(-p))


def metrics(y, w, p):
    y, w, p = np.asarray(y), np.asarray(w), np.asarray(p)
    if not len(y) or not np.isfinite(p).all() or w.sum() <= 0:
        return dict(AUC_IPCW=np.nan, AUPRC_IPCW=np.nan, Brier_IPCW=np.nan,
                    Brier_Hajek=np.nan, LogLoss_IPCW=np.nan, observed_IPCW=np.nan)
    auc = ap = np.nan
    if len(np.unique(y[w > 0])) == 2:
        auc = roc_auc_score(y, p, sample_weight=w)
        ap = average_precision_score(y, p, sample_weight=w)
    return dict(AUC_IPCW=float(auc), AUPRC_IPCW=float(ap),
                Brier_IPCW=float(np.mean(w*(y-p)**2)),
                Brier_Hajek=float(np.average((y-p)**2, weights=w)),
                LogLoss_IPCW=float(np.mean(w*loss(y, p))),
                observed_IPCW=float(np.average(y, weights=w)))


def bootstrap_contrasts(y, w, predictions, masks, groups, count, seed):
    import pandas as pd
    rng = np.random.default_rng(seed)
    pairs = [(m, "elasticnet") for m in predictions if m != "elasticnet"]
    pairs += [("panome", b) for b in ["all_reference", "random_panel", "diversity_panel", "copy1_topfit"]
              if b in predictions]
    pairs += [("panome_clinical", "clinical")]
    rows = []
    for subset, mask in masks.items():
        idx = np.flatnonzero(mask)
        if len(idx) < 30 or len(np.unique(y[idx][w[idx] > 0])) < 2:
            continue
        units = ([idx[groups[idx] == g] for g in np.unique(groups[idx])]
                 if groups is not None else None)
        point = {m: metrics(y[idx], w[idx], p[idx]) for m, p in predictions.items()}
        draws = {(m, b, k): [] for m, b in pairs for k in ["AUC_IPCW", "Brier_IPCW"]}
        for _ in range(count):
            draw = (np.concatenate([units[j] for j in rng.integers(0, len(units), len(units))])
                    if units is not None else rng.choice(idx, len(idx), replace=True))
            scores = {m: metrics(y[draw], w[draw], p[draw]) for m, p in predictions.items()}
            for key in draws:
                m, b, k = key
                delta = scores[m][k] - scores[b][k]
                if np.isfinite(delta):
                    draws[key].append(delta)
        for (m, b, k), values in draws.items():
            ci = np.quantile(values, [.025, .975]) if len(values) >= 20 else [np.nan]*2
            rows.append(dict(subset=subset, model=m, reference=b, metric=k,
                             delta=point[m][k]-point[b][k], lower=ci[0], upper=ci[1],
                             replicates=len(values), uncertainty="conditional_on_frozen_pipeline"))
    return pd.DataFrame(rows)
