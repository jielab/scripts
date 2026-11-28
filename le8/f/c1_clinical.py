"""Clinical background shared by the clinical baseline and molecular residual twins.

All fitted inputs are incident development data. Historical prevalent outcomes
never enter this module. The query API accepts features, identifiers and family
groups only. The estimator is explicitly an L2 logistic clinical model, not glmnet.
"""
from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression

import abm

VERSION = "LE8-C1-clinical-background-2026-11-16-v1"


def family_folds(groups, seed, folds):
    """Outcome-blind balanced assignment of complete build families to folds."""
    groups = np.asarray(groups).astype(str)
    levels, inverse = np.unique(groups, return_inverse=True)
    if len(levels) < folds:
        raise ValueError("Fewer independent build families than clinical OOF folds")
    order = sorted(range(len(levels)), key=lambda j: abm.tie_hash(levels[j], seed))
    assigned = np.empty(len(levels), int)
    assigned[np.asarray(order)] = np.arange(len(levels)) % folds
    return assigned[inverse]


def _fit(x, y, w, parameter, seed, max_iter):
    abm.check_classes(y, w, minimum=2)
    model = LogisticRegression(C=float(parameter), solver="lbfgs",
        max_iter=max(2000, int(max_iter)), random_state=int(seed))
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        model.fit(x, y, sample_weight=w / np.mean(w))
    return model


def _choose(x, y, w, tune_x, tune_y, tune_w, cfg, seed):
    best, trials = None, []
    for parameter in cfg.c_grid:
        model = _fit(x, y, w, parameter, seed, cfg.max_iter)
        p = model.predict_proba(tune_x)[:, 1]
        loss = float(np.average(abm.logloss_rows(tune_y, p), weights=tune_w))
        trials.append(dict(C=float(parameter), tune_logloss=loss))
        # The declared order resolves exact ties without using evaluation data.
        if best is None or loss < best[0]:
            best = loss, model, float(parameter)
    return best[1], best[2], trials


@dataclass
class ClinicalBackground:
    columns: np.ndarray
    model: object
    calibrator: object
    raw_probability: np.ndarray
    probability: np.ndarray
    donors: np.ndarray
    selected_C: float
    tuning: list
    oof_probability: np.ndarray | None = None
    oof_raw_probability: np.ndarray | None = None
    oof_fold: np.ndarray | None = None
    oof_audit: list | None = None

    def predict_raw(self, z):
        return self.model.predict_proba(np.asarray(z)[:, self.columns])[:, 1]

    def predict(self, z):
        return self.calibrator.predict(self.predict_raw(z))


def clinical_background(c, *, require_oof=False):
    """Fit once per shared context; add honest family OOF residuals when needed.

    Each OOF model selects its penalty using the independent tune_model role
    after excluding the entire donor fold. In particular it does not reuse a
    penalty chosen by a model trained on that donor's outcome. Build clinical
    medians/scales are unsupervised and shared with all comparison methods.
    """
    cache = getattr(c, "clinical_cache", None)
    if cache is None:
        cache = c.clinical_cache = {}
    columns = np.flatnonzero(np.isin(c.prep.keep, c.cfg.always_include_indices))
    if not len(columns):
        raise ValueError("clinical/abm_clinical require explicit shared clinical columns")
    key = (tuple(columns), tuple(c.cfg.c_grid), c.cfg.seed, c.cfg.folds, c.cfg.max_iter)
    b = c.parts["build"]
    b = b[c.w[b] > 0]
    tune = c.parts["tune_model"]
    tune = tune[c.w[tune] > 0]
    cal = c.parts["calibration_fit"]
    if key not in cache:
        model, parameter, tuning = _choose(c.z[b][:, columns], c.y[b], c.w[b],
            c.z[tune][:, columns], c.y[tune], c.w[tune], c.cfg, c.cfg.seed)
        raw = model.predict_proba(c.z[:, columns])[:, 1]
        calibrator = abm.Calibrator().fit(raw[cal], c.y[cal], c.w[cal])
        cache[key] = ClinicalBackground(columns, model, calibrator, raw,
            calibrator.predict(raw), b, parameter, tuning)
    fitted = cache[key]
    if require_oof and fitted.oof_probability is None:
        folds = min(c.cfg.folds, len(np.unique(c.groups[b])))
        if folds < 2:
            raise ValueError("Clinical donor residuals need at least two independent family folds")
        labels = family_folds(c.groups[b], c.cfg.seed + 1907, folds)
        oof = np.full(len(b), np.nan)
        oof_raw = np.full(len(b), np.nan)
        audit = []
        for fold in range(folds):
            train, valid = labels != fold, labels == fold
            if set(c.groups[b[train]]) & set(c.groups[b[valid]]):
                raise AssertionError("A clinical OOF donor family leaked into its fitted model")
            # The censoring nuisance also excludes held-out donor families.
            # Otherwise their event/censoring records could affect fit weights
            # even though their rows were excluded from the clinical regression.
            train_build=c.parts["build"]
            train_build=train_build[c.qc[train_build]&~np.isin(c.groups[train_build],c.groups[b[valid]])]
            fold_km=abm.CensoringKM().fit(c.development["time"][train_build],
                c.development["event"][train_build],c.cfg.horizon,c.cfg.min_censor_survival)
            train_y,train_w=fold_km.labels_weights(c.development["time"][b[train]],c.development["event"][b[train]])
            tune_y,tune_w=fold_km.labels_weights(c.development["time"][tune],c.development["event"][tune])
            model, parameter, tuning = _choose(c.z[b[train]][:, columns], train_y,
                train_w, c.z[tune][:, columns], tune_y, tune_w,
                c.cfg, c.cfg.seed + fold + 1)
            oof_raw[valid] = model.predict_proba(c.z[b[valid]][:, columns])[:, 1]
            # Donor and query residuals must share the deployed calibrated scale.
            # A full-build calibration map would indirectly reuse this donor's
            # outcome through its clinical model, so the map is fold-specific.
            cal_y,cal_w=fold_km.labels_weights(c.development["time"][cal],c.development["event"][cal])
            cal_w=cal_w*c.qc[cal]
            fold_calibrator=abm.Calibrator().fit(model.predict_proba(c.z[cal][:,columns])[:,1],cal_y,cal_w)
            oof[valid]=fold_calibrator.predict(oof_raw[valid])
            audit.append(dict(fold=fold, train_N=int(train.sum()), heldout_N=int(valid.sum()),
                train_cases=int(c.y[b[train]].sum()), heldout_cases=int(c.y[b[valid]].sum()),
                train_family_hash=abm.fingerprint(np.unique(c.groups[b[train]])),
                heldout_family_hash=abm.fingerprint(np.unique(c.groups[b[valid]])),
                selected_C=parameter, tuning=tuning, family_overlap_N=0,
                censoring_fit="fold-training build families only",
                calibration_fit="independent calibration_fit role; fold-specific clinical model and probability map"))
        if not np.isfinite(oof).all():
            raise AssertionError("An eligible clinical donor lacks an OOF prediction")
        fitted.oof_probability, fitted.oof_fold, fitted.oof_audit = oof, labels, audit
        fitted.oof_raw_probability=oof_raw
    return fitted


def fit_clinical(c, out):
    from c1_models import finish_projection, shared_support
    fitted = clinical_background(c)
    return finish_projection(c, fitted.probability, shared_support(c), out,
        fitted=dict(prep=c.prep, km=c.km, clinical=fitted),
        metadata=dict(method="clinical", engine="sklearn L2 logistic clinical background",
            actual_device="cpu", clinical_columns=[c.names[c.prep.keep[j]] for j in fitted.columns],
            input_information="shared clinical variables only", molecular_risk_features=0,
            clinical_columns_shared_with_comparators=True, selected_C=fitted.selected_C,
            tuning=fitted.tuning, probability_calibration_role="calibration_fit",
            fit_build_N=len(fitted.donors), donor_case_N=int(c.y[fitted.donors].sum()),
            case_loss_pruning=False, version=VERSION))
