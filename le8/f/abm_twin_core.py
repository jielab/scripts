"""Research implementation of training-donor pruning + outcome-blind twin selection.

This is a new, independently runnable numeric core, NOT a replacement for all of
cohort_inputs.py. The host keeps phenotype definitions, assay mapping, family rosters,
and its existing output exporter. No test labels are accepted by predict().

Estimand: fixed-horizon net risk under independent censoring. Death is treated as
censoring when represented as event=0. This is NOT a competing-risk estimator.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any
import warnings

import joblib
import numpy as np
from scipy.special import expit, logit
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import StratifiedGroupKFold

VERSION = "ABM-TWIN-2026-10-14"
ROLES = ("build", "tune_model", "tune_gate", "calibration_fit", "calibration_audit")


@dataclass(frozen=True)
class Config:
    trait: str = "cvd_cad"
    layer: str = "prot"
    seed: int = 12345
    horizon: float = 10.0
    folds: int = 5
    repeats: int = 2
    keep_fraction: float = 0.95
    pruning: str = "profile_concordance"  # literal model--outcome agreement; class/risk control retained
    panel_size: int = 10
    always_include_indices: tuple[int, ...] = ()
    profile_kind: str = "molecular_predictors"
    risk_bins: int = 5
    soft_floor: float = 0.1
    primary: str = "elasticnet_clean"
    selector: str = "support_quality"
    primary_coverage: float = 0.5
    coverages: tuple[float, ...] = (0.2, 0.4, 0.5, 0.6, 0.8, 1.0)
    c_grid: tuple[float, ...] = (0.01, 0.1, 1.0)
    pilot_c: float = 0.1
    max_iter: int = 2000
    lightgbm: bool = False
    tree_iterations: int = 150
    cores: int = 4
    feature_missing: float = 0.3
    sample_missing: float = 0.2
    dimensions: int = 24
    k: int = 64
    min_twins: int = 20
    min_family_ess: float = 15.0
    radius_quantile: float = 0.9
    prior_strength: float = 20.0
    gate_iterations: int = 120
    min_censor_survival: float = 0.05
    min_events: int = 20
    audit_min_cases: int = 30
    audit_min_controls: int = 30
    bootstrap: int = 1000
    grouping_unit: str = 'family'

    def validate(self) -> None:
        if self.grouping_unit not in {'family','individual_unverified'}:
            raise ValueError('Invalid grouping unit')
        if not (0 < self.keep_fraction <= 1 and self.horizon > 0):
            raise ValueError("Invalid retention fraction or horizon")
        if self.pruning not in {"profile_concordance", "raw_loss", "class_risk_stratified"}:
            raise ValueError("Unknown pruning mode")
        if self.selector not in {"support_quality", "absolute_error", "rank_error", "support_only"}:
            raise ValueError("Unknown outcome-blind selector")
        if self.primary_coverage not in self.coverages or any(not 0 < q <= 1 for q in self.coverages):
            raise ValueError("Primary coverage must be predeclared in coverages")
        if self.panel_size < 2 or any(j < 0 for j in self.always_include_indices):
            raise ValueError("Panel size >=2 and valid fixed covariate indices required")
        if self.folds < 2 or self.repeats < 1 or self.k < self.min_twins:
            raise ValueError("Invalid cross-fitting/matching parameters")
        if not (0 < self.radius_quantile < 1 and 0 < self.soft_floor <= 1):
            raise ValueError("Invalid radius/soft-weight setting")
        if (self.risk_bins < 1 or self.dimensions < 1 or self.cores < 1 or self.gate_iterations < 1
            or self.min_family_ess <= 0 or self.prior_strength < 0 or self.min_events < 2
            or not 0 < self.feature_missing <= 1 or not 0 <= self.sample_missing < 1
            or not 0 < self.min_censor_survival <= 1 or any(c <= 0 for c in self.c_grid)):
            raise ValueError("Invalid fitting/quality thresholds")
        if self.audit_min_cases < 2 or self.audit_min_controls < 2 or self.bootstrap < 1:
            raise ValueError("Invalid audit information requirements")


def tie_hash(value: str, seed: int) -> int:
    return int.from_bytes(hashlib.sha256(f"{seed}:{value}".encode()).digest()[:8], "big")


def fingerprint(ids: np.ndarray) -> str:
    return hashlib.sha256("\n".join(map(str, ids)).encode()).hexdigest()


def model_hash(value: Any) -> str:
    """Canonical content hash, invariant to pickle's pointer/alias bookkeeping."""
    def digest(x):
        h = hashlib.sha256()
        def put(v): h.update(v if isinstance(v, bytes) else str(v).encode())
        if x is None: put("none")
        elif isinstance(x, np.generic): return digest(x.item())
        elif isinstance(x, np.ndarray):
            put("array:"); put(x.dtype.str); put(repr(x.shape))
            if x.dtype.hasobject:
                for item in x.ravel(order="C"): put(digest(item))
            else: put(np.ascontiguousarray(x).tobytes())
        elif isinstance(x, (str, bytes, int, float, bool)):
            put(type(x).__name__); put(repr(x))
        elif isinstance(x, dict):
            put("dict")
            for k in sorted(x, key=lambda v: (type(v).__name__,repr(v))): put(digest(k)); put(digest(x[k]))
        elif isinstance(x, (list, tuple)):
            put(type(x).__name__)
            for item in x: put(digest(item))
        elif isinstance(x, (set, frozenset)):
            put(type(x).__name__)
            for item in sorted(digest(v) for v in x): put(item)
        elif type(x).__module__.startswith("lightgbm") and hasattr(x, "model_to_string"):
            put("lightgbm_booster"); put(x.model_to_string())
        elif hasattr(x, "__dict__"):
            put(type(x).__module__+"."+type(x).__qualname__); put(digest(vars(x)))
        elif hasattr(x, "__getstate__"):
            put(type(x).__module__+"."+type(x).__qualname__); put(digest(x.__getstate__()))
        else:
            raise TypeError("Unsupported fitted-state hash type: " + repr(type(x)))
        return h.digest()
    return digest(value).hex()


def json_write(path: Path, value: Any) -> None:
    def clean(x):
        if isinstance(x, dict): return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (tuple, list)): return [clean(v) for v in x]
        if isinstance(x, np.ndarray): return clean(x.tolist())
        if isinstance(x, np.generic): return clean(x.item())
        if isinstance(x, float) and not np.isfinite(x): return None
        return x
    path.write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


class CensoringKM:
    """Training-only reverse KM; fail explicitly when censoring positivity fails."""
    def fit(self, time, event, horizon: float, minimum: float):
        t, e = np.asarray(time, float), np.asarray(event, int)
        self.horizon, self.minimum = horizon, minimum
        self.times, inverse, totals = np.unique(t, return_inverse=True, return_counts=True)
        censored = np.bincount(inverse, weights=(e == 0), minlength=len(totals))
        risk = np.cumsum(totals[::-1])[::-1]
        self.survival = np.cumprod(1 - censored / risk)
        return self

    def at(self, time, left=False):
        pos = np.searchsorted(self.times, time, side="left" if left else "right") - 1
        return np.where(pos < 0, 1., self.survival[np.maximum(pos, 0)])

    def labels_weights(self, time, event):
        t, e = np.asarray(time, float), np.asarray(event, int)
        case = (e == 1) & (t <= self.horizon)
        control = t > self.horizon
        # A censoring exactly at horizon is conservatively not treated as a control.
        y = case.astype(int)
        w = np.zeros(len(t))
        gc = self.at(t[case], left=True)
        gh = float(self.at(np.array([self.horizon]))[0])
        if (len(gc) and np.min(gc) < self.minimum) or (control.any() and gh < self.minimum):
            raise ValueError("Insufficient training censoring support at the chosen horizon")
        w[case] = 1 / gc
        w[control] = 1 / gh
        return y, w


class NumericPrep:
    """Medians, retained columns and scale learned exclusively in this fit's data."""
    def fit(self, X, config: Config):
        x = np.asarray(X, float)
        self.input_columns = x.shape[1]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            variance = np.nanvar(x, axis=0)
        self.keep = np.flatnonzero((np.mean(np.isnan(x), axis=0) < config.feature_missing) & (variance > 1e-12))
        if len(self.keep) < 2: raise ValueError("Fewer than two usable numeric features")
        z = x[:, self.keep]
        self.median = np.nanmedian(z, axis=0)
        z = np.where(np.isfinite(z), z, self.median)
        self.mean, self.scale = z.mean(0), np.maximum(z.std(0), 1e-8)
        return self

    def transform(self, X):
        x = np.asarray(X, float)
        if x.ndim != 2 or x.shape[1] != self.input_columns or np.isinf(x).any():
            raise ValueError("Feature schema mismatch or infinite input")
        z = x[:, self.keep]
        return (np.where(np.isfinite(z), z, self.median) - self.mean) / self.scale

    def missing(self, X):
        return np.mean(np.isnan(np.asarray(X)[:, self.keep]), axis=1)


def check_classes(y, w, minimum=2):
    if min(np.sum((w > 0) & (y == k)) for k in (0, 1)) < minimum:
        raise ValueError("Insufficient known cases/controls; do not treat censored rows as controls")


def fit_classifier(X, y, w, config: Config, kind="elasticnet", parameter=0.1):
    ok = np.asarray(w) > 0
    check_classes(y, w, config.min_events)
    if kind == "elasticnet":
        model = LogisticRegression(C=float(parameter), solver="saga", penalty="elasticnet",
                                   l1_ratio=0.5, max_iter=config.max_iter, random_state=config.seed)
    elif kind == "lightgbm":
        from lightgbm import LGBMClassifier
        model = LGBMClassifier(n_estimators=config.tree_iterations, num_leaves=int(parameter),
                    learning_rate=.03, min_child_samples=40, reg_lambda=10, n_jobs=config.cores,
                    random_state=config.seed, verbosity=-1, deterministic=True, force_col_wise=True)
    else:
        raise ValueError(kind)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        model.fit(np.asarray(X)[ok], np.asarray(y)[ok], sample_weight=np.asarray(w)[ok] / np.mean(np.asarray(w)[ok]))
    return model


def logloss_rows(y, p):
    p = np.clip(p, 1e-7, 1 - 1e-7)
    return -(y * np.log(p) + (1-y) * np.log1p(-p))


def choose_model_panel(prep, X, y, w, config):
    """Selection uses only this training role. Never import full-cohort ProtWAS hits.

    Pilot elastic-net coefficients rank omics within an a-priori assay budget.
    A constant fitted pilot is a failure, not permission to select arbitrary hits.
    Clinical columns named by explicit indices are retained outside the omic budget.
    """
    if not set(config.always_include_indices)<=set(range(np.asarray(X).shape[1])):raise ValueError("Mandatory clinical column index outside input schema")
    fit = fit_classifier(prep.transform(X), y, w, config, parameter=config.pilot_c)
    importance = np.abs(fit.coef_[0])
    mandatory = np.array([j for j,k in enumerate(prep.keep) if k in config.always_include_indices], int)
    candidate = np.array([j for j,k in enumerate(prep.keep) if k not in config.always_include_indices], int)
    if len(candidate) and not np.any(importance[candidate] > 1e-10):
        raise ValueError("Panel pilot has no nonzero molecular coefficients; change the development-only specification")
    order = candidate[np.lexsort((prep.keep[candidate], -importance[candidate]))]
    chosen = np.sort(np.r_[mandatory, order[:config.panel_size]])
    if len(chosen) < 2: raise ValueError("Too few selected predictors")
    original = prep.keep[chosen].copy()
    for name in ('keep','median','mean','scale'):
        setattr(prep, name, getattr(prep, name)[chosen])
    prep.panel_source = "this fit's training data only"
    prep.panel_pilot_coefficients = fit.coef_[0][chosen].copy()
    prep.panel_pilot_zero_count = int(np.sum(np.abs(prep.panel_pilot_coefficients)<1e-10))
    return original


class ModelProfileGeometry:
    """A vector of selected standardized predictors, NOT a 1-D risk score.

    Opposing predictor changes cannot cancel as they could in beta'X.
    Every selected predictor retains a positive distance weight.
    """
    def fit(self, X, coefficients):
        beta = np.abs(np.asarray(coefficients, float))
        self.weights = np.sqrt(.25 + .75 * beta / max(beta.mean(), 1e-8))
        self.n_components_ = X.shape[1]
        return self
    def transform(self, X):
        return np.asarray(X) * self.weights


def crossfit_quality(X, time, event, ids, groups, config: Config):
    """Own outcomes are NEVER in the teacher fit. They ARE used to score donors.

    This use of build labels is deliberate supervised data selection. It is not
    evidence that discarded labels are erroneous. Fold preprocessing and KM are
    refitted; hyperparameters are fixed pilot settings, not tuned on recipient labels.
    """
    n = len(ids)
    prob = np.full((config.repeats, n), np.nan)
    losses = np.full_like(prob, np.nan)
    logs = []
    strat = ((np.asarray(event) == 1) & (np.asarray(time) <= config.horizon)).astype(int)
    for repeat in range(config.repeats):
        splitter = StratifiedGroupKFold(config.folds, shuffle=True, random_state=config.seed + repeat)
        for fold, (tr, va) in enumerate(splitter.split(X, strat, groups)):
            if set(groups[tr]) & set(groups[va]): raise AssertionError("Family leakage in OOF")
            prep = NumericPrep().fit(X[tr], config)
            km = CensoringKM().fit(time[tr], event[tr], config.horizon, config.min_censor_survival)
            y, w = km.labels_weights(time[tr], event[tr])
            panel = choose_model_panel(prep, X[tr], y, w, config)
            model = fit_classifier(prep.transform(X[tr]), y, w, config, parameter=config.pilot_c)
            p = model.predict_proba(prep.transform(X[va]))[:, 1]
            yy, ww = km.labels_weights(time[va], event[va])
            prob[repeat, va] = p
            losses[repeat, va] = np.where(ww > 0, logloss_rows(yy, p), np.nan)
            logs.append(dict(repeat=repeat, fold=fold, train_n=len(tr), recipient_n=len(va),
                             recipient_hash=fingerprint(ids[va]), family_disjoint=True, panel_input_indices=panel.tolist(),
                             selection_in_teacher_only=True))
    if not np.isfinite(prob).all(): raise AssertionError("Incomplete OOF predictions")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        quality = np.nanmedian(losses, axis=0)
    return dict(oof_probability=prob.mean(0), oof_sd=prob.std(0),
                oof_loss=quality, model_concordance=np.exp(-quality), repeat_losses=losses, audit=logs)


def prune_donors(y, w, quality, ids, config: Config):
    """Never deletes rows from the evaluation population. Returns training weights."""
    valid = (w > 0) & np.isfinite(quality["oof_loss"])
    p, err = quality["oof_probability"], quality["oof_loss"]
    edges = np.unique(np.quantile(p[valid], np.linspace(0, 1, config.risk_bins + 1)))[1:-1]
    bins = np.searchsorted(edges, p, side="right")
    # Raw-loss pruning is provided as the literal user hypothesis, not silently discarded.
    strata = y * (config.risk_bins + 1) + bins if config.pruning == "class_risk_stratified" else np.zeros(len(y), int)
    keep, soft, control = np.zeros(len(y), bool), np.zeros(len(y)), np.zeros(len(y), bool)
    ties = np.asarray([tie_hash(s, config.seed) for s in ids], np.uint64)
    audit = []
    for label in np.unique(strata[valid]):
        idx = np.flatnonzero(valid & (strata == label))
        order = idx[np.lexsort((ties[idx], err[idx]))]
        count = max(1, math.ceil(config.keep_fraction * len(idx)))
        keep[order[:count]] = True
        rank = (np.arange(len(order)) + .5) / len(order)
        soft[order] = config.soft_floor + (1-config.soft_floor) * (1-rank)
        audit.append(dict(stratum=int(label), total=len(idx), retained=count))
    # Control training subset preserves CASE x OOF-RISK-bin composition of pruning.
    for label in np.unique((y * (config.risk_bins + 1) + bins)[valid]):
        idx = np.flatnonzero(valid & ((y * (config.risk_bins + 1) + bins) == label))
        r = sorted(idx, key=lambda j: tie_hash(ids[j], config.seed + 701))
        control[np.asarray(r[:int(keep[idx].sum())], int)] = True
    return dict(keep=keep, soft=soft, random_keep=control, eligible=valid, strata=strata, audit=audit)


class Calibrator:
    def fit(self, p, y, w):
        check_classes(y, w)
        ok = w > 0
        z = logit(np.clip(p, 1e-7, 1-1e-7)).reshape(-1, 1)
        self.model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
        self.model.fit(z[ok], y[ok], sample_weight=w[ok]/w[ok].mean())
        return self
    def predict(self, p):
        return self.model.predict_proba(logit(np.clip(p, 1e-7, 1-1e-7)).reshape(-1, 1))[:, 1]


class TwinBank:
    """All-radius candidates counted exactly; only <=k are used for borrowing.

    The supplied frozen model-panel geometry defines proximity, not query outcomes. Radius is learned
    from independent tune_model queries; family and self matches are excluded.
    """
    def fit(self, X, y, w, ids, groups, losses, config: Config, geometry=None):
        self.config = config
        self.pca = geometry or PCA(n_components=min(config.dimensions, X.shape[1], len(X)-1),
                                  svd_solver="randomized", random_state=config.seed).fit(X)
        self.Z = self.pca.transform(X) / np.sqrt(self.pca.n_components_)
        self.y, self.w = np.asarray(y), np.asarray(w)
        self.ids, self.groups, self.losses = ids.copy(), groups.copy(), losses.copy()
        from sklearn.neighbors import KDTree
        self.index = KDTree(self.Z,leaf_size=40) if self.Z.shape[1] <= 24 else None
        self.radius = None
        return self

    def calibrate_radius(self, X, ids, groups):
        z = self.pca.transform(X) / np.sqrt(self.pca.n_components_)
        kth = []
        for start in range(0, len(z), 128):
            query=z[start:start+128]
            if self.index is not None:
                distance,near=self.index.query(query,k=min(self.config.k+64,len(self.Z)))
                for j,(ds,nn) in enumerate(zip(distance,near)):
                    ok=(self.groups[nn]!=groups[start+j])&(self.ids[nn]!=ids[start+j])
                    good=ds[ok]
                    if len(good)>=self.config.k:kth.append(good[self.config.k-1]);continue
                    # Unusually large families: exact full-radius fallback for this query.
                    dd=euclidean_distances(query[j:j+1],self.Z)[0]
                    dd[(self.groups==groups[start+j])|(self.ids==ids[start+j])]=np.inf
                    kth.append(np.partition(dd,min(self.config.k,len(dd))-1)[min(self.config.k,len(dd))-1])
            else:
                d = euclidean_distances(query, self.Z)
                valid = (groups[start:start+128, None] != self.groups[None, :]) & (ids[start:start+128, None] != self.ids[None, :])
                d[~valid] = np.inf
                kk = min(self.config.k, d.shape[1]) - 1
                kth.extend(np.partition(d, kk, axis=1)[:, kk])
        finite = np.asarray(kth)[np.isfinite(kth)]
        if len(finite) < 10: raise ValueError("Not enough independent queries to calibrate matching radius")
        self.radius = float(np.quantile(finite, self.config.radius_quantile))
        if self.radius <= 1e-12: raise ValueError("Collapsed matching geometry")
        return self

    def predict(self, X, ids, groups, fallback):
        if self.radius is None: raise RuntimeError("Radius is not frozen")
        z = self.pca.transform(X) / np.sqrt(self.pca.n_components_)
        n = len(z)
        out = {k: np.zeros(n) for k in ("eligible_match_count", "matched_count", "family_count", "family_ess",
            "donor_ess", "donor_loss", "radius_ratio", "nearest_radius_ratio", "label_disagreement", "borrow_fraction")}
        out["donor_loss"].fill(10.)
        out["radius_ratio"].fill(100.)
        out["nearest_radius_ratio"].fill(100.)
        prob = np.asarray(fallback, float).copy()
        for start in range(0, n, 128):
            query=z[start:start+128]
            if self.index is not None:
                near,distances=self.index.query_radius(query,r=self.radius,return_distance=True)
            else:
                dense=euclidean_distances(query,self.Z)
                near=[np.flatnonzero(dd<=self.radius) for dd in dense]
                distances=[dd[ii] for dd,ii in zip(dense,near)]
            for row,(nn,ds) in enumerate(zip(near,distances)):
                i = start + row
                allowed=(self.groups[nn]!=groups[i])&(self.ids[nn]!=ids[i])&(self.w[nn]>0)
                idx=nn[allowed];ds=ds[allowed]
                out["eligible_match_count"][i] = len(idx)
                if not len(idx): continue
                # Stable tie order is independent of donor outcomes.
                order=np.lexsort((self.ids[idx],ds))[:self.config.k]
                idx=idx[order];ds=ds[order]
                weight = np.exp(-.5*(ds/self.radius)**2) * self.w[idx]
                weight /= weight.sum()
                _, inv = np.unique(self.groups[idx], return_inverse=True)
                fw = np.bincount(inv, weights=weight)
                ess, fess = 1/np.sum(weight**2), 1/np.sum(fw**2)
                local = float(weight @ self.y[idx])
                borrowing = fess / (fess + self.config.prior_strength)
                prob[i] = borrowing * local + (1-borrowing) * fallback[i]
                for name, value in dict(matched_count=len(idx), family_count=len(fw), family_ess=fess, donor_ess=ess,
                    donor_loss=weight @ self.losses[idx], radius_ratio=ds[-1]/self.radius,
                    nearest_radius_ratio=ds[0]/self.radius,
                    label_disagreement=local*(1-local), borrow_fraction=borrowing).items(): out[name][i] = value
        out["supported"] = (out["matched_count"] >= self.config.min_twins) & (out["family_ess"] >= self.config.min_family_ess)
        return np.clip(prob, 1e-7, 1-1e-7), out


def pairwise_error_target(y, w, probability):
    """Per-person, weighted pair-discordance on independent gate-training data.

    This is an AUC-oriented surrogate, not a guarantee about a jointly selected
    subpopulation's AUC. O(n log n), no n-by-n pair matrix.
    """
    p, y, w = np.asarray(probability), np.asarray(y), np.asarray(w)
    target = np.full(len(p), .5)
    check_classes(y, w)
    for label in (0, 1):
        other = (y != label) & (w > 0)
        order = np.argsort(p[other], kind="stable")
        po, wo = p[other][order], w[other][order]
        cumulative = np.r_[0., np.cumsum(wo)]
        i = np.flatnonzero((y == label) & (w > 0))
        left = cumulative[np.searchsorted(po, p[i], side="left")]
        right = cumulative[np.searchsorted(po, p[i], side="right")]
        total = cumulative[-1]
        target[i] = ((total-right)+.5*(right-left))/total if label == 1 else (left+.5*(right-left))/total
    return target


class FrozenRule:
    def fit(self, score, eligible, ids, coverage, seed):
        self.seed, self.coverage = seed, coverage
        idx = np.flatnonzero(eligible & np.isfinite(score))
        self.support_limited = len(idx) < math.ceil(coverage * len(ids))
        # Coverage is relative to ALL calibration queries, never only supported ones.
        count = min(len(idx), math.ceil(coverage * len(ids)))
        self.threshold, self.tie = np.inf, 0
        if count:
            ties = np.array([tie_hash(v, seed) for v in ids], np.uint64)
            order = idx[np.lexsort((ties[idx], -np.asarray(score)[idx]))]
            boundary = order[count-1]
            self.threshold, self.tie = float(score[boundary]), int(ties[boundary])
            if count == len(idx): self.threshold = -np.inf
        return self
    def apply(self, score, eligible, ids):
        ties = np.array([tie_hash(v, self.seed) for v in ids], np.uint64)
        score = np.asarray(score)
        return eligible & np.isfinite(score) & ((score > self.threshold) | ((score == self.threshold) & (ties <= self.tie)))


class TwinABM:
    def __init__(self, config: Config = Config()):
        config.validate()
        self.config = config

    @staticmethod
    def validate_features(X, ids, groups):
        x = np.asarray(X, float); ids, groups = np.asarray(ids).astype(str), np.asarray(groups).astype(str)
        if x.ndim != 2 or len(x) != len(ids) or len(ids) != len(groups) or np.isinf(x).any():
            raise ValueError("Invalid numeric feature input")
        invalid = {"", "nan", "none", "na", "<na>"}
        if len(set(ids)) != len(ids) or any(v.strip().lower() in invalid for v in ids) or any(v.strip().lower() in invalid for v in groups):
            raise ValueError("IDs must be unique and all family groups explicitly defined")
        return x, ids, groups

    def fit(self, X, time, event, ids, groups, roles, feature_names):
        a = self.config
        X, ids, groups = self.validate_features(X, ids, groups)
        time, event, roles = np.asarray(time, float), np.asarray(event), np.asarray(roles).astype(str)
        if len(time) != len(X) or len(event) != len(X) or len(roles) != len(X) or not np.isfinite(time).all() or (time <= 0).any() or not np.isin(event, [0, 1]).all():
            raise ValueError("Invalid survival outcomes")
        if set(roles) != set(ROLES): raise ValueError("Development must contain exactly the five declared roles, never test")
        group_roles = {}
        for g, r in zip(groups, roles):
            if g in group_roles and group_roles[g] != r: raise ValueError("Family crosses development roles")
            group_roles[g] = r
        self.feature_names = list(map(str, feature_names))
        if len(self.feature_names) != X.shape[1] or len(set(self.feature_names)) != len(self.feature_names):
            raise ValueError("Invalid named numeric schema")
        self.development_ids, self.development_groups = ids.copy(), groups.copy()
        self.audit_ids = ids[roles == "calibration_audit"].copy()
        self.audit_groups = groups[roles == "calibration_audit"].copy()
        parts = {r: np.flatnonzero(roles == r) for r in ROLES}
        b, t, gate, cal = (parts[r] for r in ROLES[:4])
        self.prep = NumericPrep().fit(X[b], a)
        self.km = CensoringKM().fit(time[b], event[b], a.horizon, a.min_censor_survival)
        # Audit outcomes are not even passed to label/weight computation during fit.
        non_audit = roles != "calibration_audit"
        y, w = np.zeros(len(X), int), np.zeros(len(X))
        y[non_audit], w[non_audit] = self.km.labels_weights(time[non_audit], event[non_audit])
        for r in ROLES[:4]: check_classes(y[parts[r]], w[parts[r]], a.min_events)
        self.panel_indices = choose_model_panel(self.prep, X[b], y[b], w[b], a)
        self.panel_names = [self.feature_names[j] for j in self.panel_indices]
        z = self.prep.transform(X)
        quality = crossfit_quality(X[b], time[b], event[b], ids[b], groups[b], a)
        pruning = prune_donors(y[b], w[b], quality, ids[b], a)
        multipliers = {"full":np.ones(len(b)), "clean":pruning["keep"].astype(float),
                       "soft":pruning["soft"], "random":pruning["random_keep"].astype(float)}
        self.models, self.model_tuning = {}, []
        for kind in (["elasticnet", "lightgbm"] if a.lightgbm else ["elasticnet"]):
            for arm, mult in multipliers.items():
                best = None
                for par in (a.c_grid if kind == "elasticnet" else (7, 15, 31)):
                    model = fit_classifier(z[b], y[b], w[b]*mult, a, kind, par)
                    p = model.predict_proba(z[t])[:, 1]
                    score = float(np.average(logloss_rows(y[t], p), weights=w[t]))
                    self.model_tuning.append(dict(model=kind+"_"+arm, parameter=par, logloss=score))
                    if best is None or score < best[0]: best = score, model
                self.models[kind+"_"+arm] = best[1]
        self.banks = {}
        eligible = pruning["eligible"]
        common_pca = ModelProfileGeometry().fit(z[b], self.models["elasticnet_full"].coef_[0])
        for arm, select in (("full", eligible), ("clean", pruning["keep"])):
            bb = b[select]
            bank = TwinBank().fit(z[bb], y[bb], w[bb], ids[bb], groups[bb], quality["oof_loss"][select], a, common_pca)
            bank.calibrate_radius(z[t], ids[t], groups[t])
            self.banks[arm] = bank
        self.calibrators = {}
        pred, _ = self._raw(X[cal], ids[cal], groups[cal])
        self.calibrators = {name:Calibrator().fit(p, y[cal], w[cal]) for name,p in pred.items()}
        if a.primary not in self.calibrators: raise ValueError("Requested primary model was not fitted")
        # Gate is trained on independent gate-role outcomes, never build residual labels.
        gp, gd = self._calibrated(X[gate], ids[gate], groups[gate])
        gx = self._gate_features(gp, gd)
        self.gate_models = {}
        targets = {"absolute_error":(y[gate]-gp[a.primary])**2,
                   "rank_error":pairwise_error_target(y[gate], w[gate], gp[a.primary])}
        known = w[gate] > 0
        for name, target in targets.items():
            model = HistGradientBoostingRegressor(max_iter=a.gate_iterations, max_leaf_nodes=5, learning_rate=.04,
                     min_samples_leaf=max(20, min(100, int(known.sum()/10))), l2_regularization=10,
                     early_stopping=False, random_state=a.seed+903)
            model.fit(gx[known], target[known], sample_weight=w[gate][known]/w[gate][known].mean())
            self.gate_models[name] = model
        cp, cd = self._calibrated(X[cal], ids[cal], groups[cal])
        score, _ = self._scores(cp, cd, ids[cal])
        self.rules = {}
        technical = cd["technical_qc"]
        for name, values in score.items():
            allowed = technical if name in ("confidence", "lowrisk", "random") else technical & cd["supported"]
            for q in a.coverages:
                self.rules[(name,q)] = FrozenRule().fit(values, allowed, ids[cal], q, a.seed+91)
        self.risk_edges = np.unique(np.quantile(cp["elasticnet_full"], np.linspace(0, 1, 11)))[1:-1]
        self.quality_edges = np.unique(np.quantile(score[a.selector], np.linspace(0,1,6)))[1:-1]
        self.training_audit = dict(oof=quality["audit"], pruning=pruning["audit"], known_donors=int(eligible.sum()),
            retained_donors=int(pruning["keep"].sum()), build_n=len(b),
            selected_panel=self.panel_names, panel_budget=a.panel_size,
            selected_zero_pilot_coefficients=self.prep.panel_pilot_zero_count,
            final_nonzero_coefficients={n:int(np.sum(np.abs(m.coef_[0])>1e-10)) for n,m in self.models.items() if hasattr(m,"coef_")},
            matching_space="selected predictor vector; coefficient-weighted distances",
            confidence_definition="100*geometric_mean(proximity,independent_support,donor_concordance,profile_radius); zero without support",
            retained_cases=int(np.sum(y[b][pruning["keep"]])), known_cases=int(np.sum(y[b])),
            literal_pruning_is_supervised_selection_not_proof_of_noise=True,
            roles={r:dict(N=len(parts[r]),cases=int(y[parts[r]].sum()) if r != "calibration_audit" else None,
                         id_hash=fingerprint(ids[parts[r]])) for r in ROLES})
        self.donor_table = dict(eid=ids[b], group=groups[b], oof_probability=quality["oof_probability"],
            oof_sd=quality["oof_sd"], oof_loss=quality["oof_loss"],
            model_concordance=quality["model_concordance"], retained=pruning["keep"],
            soft_weight=pruning["soft"], random_retained=pruning["random_keep"], horizon_known=eligible)
        self.frozen = True
        # calibration_audit labels have not been used for fitting or threshold choice.
        return self

    def _raw(self, X, ids, groups):
        z = self.prep.transform(X)
        pred = {name:model.predict_proba(z)[:, 1] for name,model in self.models.items()}
        diag = None
        for arm, bank in self.banks.items():
            pred["twin_"+arm], d = bank.predict(z, ids, groups, pred["elasticnet_full"])
            if arm == "clean": diag = d
        diag["technical_qc"] = self.prep.missing(X) <= self.config.sample_missing
        return pred, diag

    def _calibrated(self, X, ids, groups):
        pred, diag = self._raw(X, ids, groups)
        return {name:self.calibrators[name].predict(p) for name,p in pred.items()}, diag

    def _gate_features(self, p, d):
        return np.c_[p[self.config.primary], p["elasticnet_full"], p["twin_clean"],
              np.log1p(d["eligible_match_count"]), np.log1p(d["family_ess"]),
              d["radius_ratio"], d["donor_loss"], d["label_disagreement"], d["borrow_fraction"]]

    def _scores(self, p, d, ids):
        gx = self._gate_features(p, d)
        est = {n:m.predict(gx) for n,m in self.gate_models.items()}
        proximity = np.exp(-np.minimum(d["nearest_radius_ratio"],100)**2)
        support = np.clip(d["family_ess"] / self.config.k, 0, 1)
        concordance = np.exp(-np.maximum(d["donor_loss"],0))
        profile_radius = np.exp(-np.minimum(d["radius_ratio"],100)**2)
        confidence = np.where(d["supported"] & d["technical_qc"],
                              (proximity * support * concordance * profile_radius)**.25, 0.)
        d.update(proximity_score=proximity, independent_support_score=support,
                 donor_concordance_score=concordance, profile_radius_score=profile_radius,
                 pre_prediction_confidence=100*confidence)
        score = {"absolute_error":-est["absolute_error"], "rank_error":-est["rank_error"],
          "support_quality":confidence,
          "support_only":-d["radius_ratio"] + .01*np.log1p(d["family_ess"]),
          "confidence":np.abs(p["elasticnet_full"]-.5), "lowrisk":-p["elasticnet_full"],
          "random":np.array([tie_hash(v,self.config.seed+99)/2**64 for v in ids])}
        return score, est

    def predict(self, X, ids, groups, feature_names, *, allow_development_audit=False):
        """Outcome-blind API. There is intentionally no y/time/event argument."""
        X, ids, groups = self.validate_features(X, ids, groups)
        if list(map(str, feature_names)) != self.feature_names: raise ValueError("Feature names/order mismatch")
        if allow_development_audit and (not set(ids) <= set(self.audit_ids) or not set(groups) <= set(self.audit_groups)):
            raise ValueError("Audit override is restricted to the predeclared audit roster")
        if not allow_development_audit and (set(ids) & set(self.development_ids) or set(groups) & set(self.development_groups)):
            raise ValueError("Projection overlaps development participants/families")
        p, d = self._calibrated(X, ids, groups)
        scores, est = self._scores(p, d, ids)
        masks = {}
        for (name,q), rule in self.rules.items():
            eligible = d["technical_qc"] if name in ("confidence","lowrisk","random") else d["technical_qc"] & d["supported"]
            masks[f"{name}:{q:g}"] = rule.apply(scores[name], eligible, ids)
        selected = masks[f"{self.config.selector}:{self.config.primary_coverage:g}"]
        p["fallback_policy"] = np.where(d["technical_qc"], np.where(selected, p[self.config.primary], p["elasticnet_full"]), np.nan)
        reason = np.where(~d["technical_qc"],"technical_qc",np.where(~d["supported"],"no_adequate_twins",
                    np.where(selected,"selected","below_frozen_quality_cutoff")))
        return dict(eid=ids, group=groups, probability=p, diagnostics=d, estimated_error=est,
                    masks=masks, selected=selected, reason=reason, selector_scores=scores,
                    frozen_model_hash=model_hash(self), version=VERSION,
                    confidence_semantics="reference support, not probability of correct outcome",
                    evaluation_scope="calibration_audit" if allow_development_audit else "heldout_projection",
                    profile_kind=self.config.profile_kind, panel_features=self.panel_names)
