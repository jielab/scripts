#!/usr/bin/env python3
"""Training-frozen DiscoDivas with genuine fine-tuning-cohort input models.

References
----------
Ruan et al., AJHG 2026;113:1448-1460, doi:10.1016/j.ajhg.2026.05.006.
https://github.com/YunfengRuan/DiscoDivas/blob/main/DiscoDivas.R
https://github.com/YunfengRuan/DiscoDivas/blob/main/supplemental/generate.metrics.R

Every anchor model uses ALL supplied discovery-population scores and covariates,
fitted only in that anchor's fine-tuning cohort. Only the genetic linear term is
interpolated. In particular, four uncombined population scores are not relabelled
as four fine-tuned models. ``RawScoreInterpolation`` exposes that older procedure
as the explicitly named ``disco.raw_interpolation`` ablation.

All PC coordinates must already be in the SAME reference-projected space and
units. Distances use the supplied coordinates, with no per-axis or per-ancestry
rescaling. Anchor centers are medians of their actual training individuals. The
official interpolation is w_ik proportional to d_k r_k / D_ik, d = G^-1 1, with
r_k = 1 by default. r is an optional quality factor, not a distance statistic.

The official R script residualizes and standardizes scores in the supplied target
batch without using phenotypes. Here that transformation is learned once in TRAIN
and frozen, so predicting one person does not depend on other target people. The
default tiny fixed L2 penalty stabilizes cohort regressions; set penalty=0 for the
unpenalized fit. Binary anchor fits use logistic regression. These choices are
recorded and are not described as a numerical reproduction of the original paper.

When a downstream phenotype/calibration model is needed, obtain training scores
by cross-fitting this ENTIRE estimator, respecting family groups. Fit downstream
calibration on those out-of-fold training scores; refit this estimator on all
training data for target prediction. Fitting downstream calibration to in-sample
anchor scores would add another source of optimism.

Numerical boundaries are explicit: insufficient/constant anchors are omitted;
coincident or ill-conditioned anchor geometry is rejected, without pseudoinverse
or a relabelled global fallback. At an isolated nonzero-weight anchor center, the
continuous inverse-distance limit is a one-hot weight. Signed d/weights are kept
as in the published formula; a cancelling denominator gives an unavailable score.
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from scipy.special import expit
from sklearn.base import BaseEstimator


DEFAULT_PCS = tuple(f"PC{i}" for i in range(1, 11))
DEFAULT_COVARIATES = ("age", "sex") + tuple(f"PC{i}" for i in range(1, 21))
DEFAULT_ANCHORS = ("EUR", "AFR", "EAS", "SAS")


class DiscoFitError(ValueError):
    """The requested Disco model cannot be fitted with the supplied training data."""


def _matrix(frame, columns):
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError("Missing feature columns: " + ", ".join(missing))
    return frame.loc[:, list(columns)].apply(pd.to_numeric, errors="raise").to_numpy(float)


def _training_ids(frame, id_col, role_col):
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("X must be a pandas DataFrame with named score and PC columns")
    if role_col in frame and not frame[role_col].astype(str).eq("train").all():
        raise ValueError("fit accepts TRAIN rows only; split target roles before calling fit")
    ids = frame[id_col].astype(str) if id_col in frame else frame.index.astype(str)
    if pd.Index(ids).has_duplicates:
        raise ValueError("Training IDs must be unique")
    return tuple(ids)


class DistanceInterpolator:
    """The published distance formula, with documented limits and failure states."""

    def __init__(self, centers, quality=None, max_condition=1e12, tolerance=1e-12):
        self.centers = np.asarray(centers, dtype=float)
        if (self.centers.ndim != 2 or len(self.centers) < 2
                or not np.isfinite(self.centers).all()):
            raise DiscoFitError("At least two finite PC anchor centers are required")
        self.quality = (np.ones(len(self.centers)) if quality is None
                        else np.asarray(quality, dtype=float))
        if (self.quality.shape != (len(self.centers),)
                or not np.isfinite(self.quality).all() or np.any(self.quality <= 0)):
            raise ValueError("Active-anchor quality factors must be finite and positive")
        self.tolerance = tolerance
        self.distance_matrix = cdist(self.centers, self.centers)
        off_diagonal = self.distance_matrix[np.triu_indices(len(self.centers), 1)]
        if np.any(off_diagonal <= 0):
            raise DiscoFitError("Coincident anchor PC centers: G is singular")
        self.distance_scale = float(np.median(off_diagonal))
        scaled = self.distance_matrix / self.distance_scale
        self.condition_number = float(np.linalg.cond(scaled))
        if not np.isfinite(self.condition_number) or self.condition_number > max_condition:
            raise DiscoFitError("Ill-conditioned anchor PC geometry: no stable G inverse")
        try:
            self.d = np.linalg.solve(scaled, np.ones(len(self.centers))) / self.distance_scale
        except np.linalg.LinAlgError as exc:
            raise DiscoFitError("Singular anchor PC geometry") from exc
        if not np.isfinite(self.d).all():
            raise DiscoFitError("Nonfinite distance shrinkage coefficients")

    def weights(self, pcs, return_status=False):
        pcs = np.asarray(pcs, dtype=float)
        if pcs.ndim != 2 or pcs.shape[1] != self.centers.shape[1]:
            raise ValueError("PC dimensions do not match the fitted reference space")
        weights = np.full((len(pcs), len(self.centers)), np.nan)
        status = np.full(len(pcs), "missing_pc", dtype=object)
        good = np.isfinite(pcs).all(axis=1)
        rows = np.flatnonzero(good)
        if not len(rows):
            return (weights, status) if return_status else weights
        distances = cdist(pcs[rows], self.centers) / self.distance_scale
        strength = self.d * self.quality
        # This common normalization changes neither the formula nor its signs.
        strength = strength / np.max(np.abs(strength))
        active = np.abs(strength) > np.finfo(float).eps
        exact = (distances <= self.tolerance) & active[None, :]
        at_center = exact.sum(axis=1) == 1
        weights[rows[at_center]] = exact[at_center].astype(float)
        status[rows[at_center]] = "anchor_limit"
        ordinary = ~at_center
        if ordinary.any():
            dist = distances[ordinary]
            # d=0 contributes zero even when its distance is zero.
            raw = np.zeros_like(dist)
            np.divide(strength, dist, out=raw,
                      where=(dist > self.tolerance) & active[None, :])
            denominator = raw.sum(axis=1)
            total_absolute = np.abs(raw).sum(axis=1)
            stable = (np.isfinite(raw).all(axis=1) & (total_absolute > 0)
                      & (np.abs(denominator) > self.tolerance * total_absolute)
                      & (exact[ordinary].sum(axis=1) == 0))
            idx = rows[ordinary]
            status[idx] = "unstable_weight_denominator"
            weights[idx[stable]] = raw[stable] / denominator[stable, None]
            status[idx[stable]] = "ok"
        return (weights, status) if return_status else weights


class FrozenPCNormalizer:
    """OLS score~PC residuals with TRAIN-only coefficients, mean and sample SD."""

    def fit(self, pcs, scores):
        pcs, scores = np.asarray(pcs, float), np.asarray(scores, float)
        design = np.column_stack((np.ones(len(pcs)), pcs))
        self.coef_ = np.linalg.lstsq(design, scores, rcond=None)[0]
        residuals = scores - design @ self.coef_
        self.mean_ = np.mean(residuals, axis=0)
        self.scale_ = np.std(residuals, axis=0, ddof=1)
        raw_scale = np.maximum(np.std(scores, axis=0, ddof=1), 1e-12)
        self.valid_ = (np.isfinite(self.scale_)
                       & (self.scale_ > np.maximum(1e-12, raw_scale * 1e-10)))
        return self

    def transform(self, pcs, scores):
        design = np.column_stack((np.ones(len(pcs)), np.asarray(pcs, float)))
        residuals = np.asarray(scores, float) - design @ self.coef_
        return (residuals - self.mean_) / self.scale_


def _fit_regression(design, y, task, penalty, max_iter):
    """Intercept is unpenalized; all supplied score/covariate terms are retained."""
    x = np.column_stack((np.ones(len(design)), design))
    if task == "continuous":
        if penalty == 0:
            beta = np.linalg.lstsq(x, y, rcond=None)[0]
        else:
            reg = np.sqrt(penalty) * np.eye(x.shape[1])
            reg[0, 0] = 0
            beta = np.linalg.lstsq(np.vstack((x, reg)),
                                  np.concatenate((y, np.zeros(x.shape[1]))), rcond=None)[0]
        return beta, "linear", 0
    # Mean loss/gradient avoid sample-size-dependent optimizer stopping behavior.
    def objective(beta):
        eta = x @ beta
        loss = np.mean(np.logaddexp(0.0, eta) - y * eta)
        loss += penalty * np.dot(beta[1:], beta[1:]) / (2 * len(y))
        gradient = x.T @ (expit(eta) - y) / len(y)
        gradient[1:] += penalty * beta[1:] / len(y)
        return loss, gradient

    start = np.zeros(x.shape[1])
    prevalence = np.clip(np.mean(y), 1e-8, 1 - 1e-8)
    start[0] = np.log(prevalence / (1 - prevalence))
    fit = minimize(objective, start, jac=True, method="L-BFGS-B",
                   options={"maxiter": max_iter, "gtol": 1e-7, "ftol": 1e-12})
    if not fit.success or not np.isfinite(fit.x).all():
        raise DiscoFitError("Anchor logistic fit did not converge: " + str(fit.message))
    return fit.x, "logistic", int(fit.nit)


class DiscoDivas(BaseEstimator):
    """Fit cohort-specific multivariate PRS inputs, then interpolate genetic terms.

    ``X`` contains id, ancestry, optional role, named scores and reference PCs.
    ``y`` contains labels for exactly the passed TRAIN rows. ``covariates`` can be
    an aligned finite numeric matrix from the caller's common covariate encoder;
    otherwise ``covariate_cols`` are read as numeric. It is never required for
    prediction, because covariate predictions are excluded from the genetic score.
    The resulting score is NOT an outcome prediction or a disease probability.
    """

    def __init__(self, score_cols, task="continuous", anchors=DEFAULT_ANCHORS,
                 pc_cols=DEFAULT_PCS, covariate_cols=DEFAULT_COVARIATES,
                 ancestry_col="ancestry", id_col="id", role_col="role",
                 min_anchor_n=100, min_cases=20, penalty=1e-6, max_iter=2000,
                 quality=None, source_method="PRS-CSx"):
        self.score_cols = score_cols
        self.task = task
        self.anchors = anchors
        self.pc_cols = pc_cols
        self.covariate_cols = covariate_cols
        self.ancestry_col = ancestry_col
        self.id_col = id_col
        self.role_col = role_col
        self.min_anchor_n = min_anchor_n
        self.min_cases = min_cases
        self.penalty = penalty
        self.max_iter = max_iter
        self.quality = quality
        self.source_method = source_method

    def fit(self, X, y, *, covariates=None):
        self.is_fitted_ = False
        for name in ("geometry_", "normalizer_"):
            self.__dict__.pop(name, None)
        self.anchor_records_, self.anchor_models_, self.anchor_fit_ids_ = [], {}, {}
        self.fit_ids_ = _training_ids(X, self.id_col, self.role_col)
        if self.task not in ("continuous", "binary"):
            raise ValueError("task must be continuous or binary")
        if self.penalty < 0 or not np.isfinite(self.penalty):
            raise ValueError("penalty must be finite and nonnegative")
        if self.min_anchor_n < 2 or self.min_cases < 1:
            raise ValueError("min_anchor_n >= 2 and min_cases >= 1 are required")
        if self.ancestry_col not in X:
            raise ValueError("Missing fine-tuning cohort column: " + self.ancestry_col)
        self.score_cols_, self.pc_cols_ = tuple(self.score_cols), tuple(self.pc_cols)
        if not self.score_cols_ or not self.pc_cols_:
            raise ValueError("At least one discovery score and one reference PC are required")
        if len(set(self.score_cols_)) != len(self.score_cols_):
            raise ValueError("Duplicate discovery score names")
        if len(set(self.anchors)) != len(self.anchors):
            raise ValueError("Duplicate anchor names")
        scores, pcs = _matrix(X, self.score_cols_), _matrix(X, self.pc_cols_)
        if isinstance(y, pd.Series) and not y.index.equals(X.index):
            raise ValueError("A Series outcome must have the same row index as X")
        y = np.asarray(y, dtype=float).reshape(-1)
        if len(y) != len(X) or not np.isfinite(y).all():
            raise ValueError("y must contain one finite outcome per TRAIN row")
        if not np.isfinite(scores).all() or not np.isfinite(pcs).all():
            raise ValueError("Training discovery scores and reference PCs must be complete")
        if self.task == "binary" and not np.isin(y, [0, 1]).all():
            raise ValueError("Binary outcomes must be coded 0/1")
        if covariates is None:
            cov = _matrix(X, tuple(self.covariate_cols))
            self.covariate_names_ = tuple(self.covariate_cols)
        else:
            if isinstance(covariates, pd.DataFrame) and not covariates.index.equals(X.index):
                raise ValueError("A covariate DataFrame must have the same row index as X")
            cov = np.asarray(covariates, dtype=float)
            if cov.ndim == 1:
                cov = cov[:, None]
            self.covariate_names_ = (tuple(map(str, covariates.columns))
                                     if isinstance(covariates, pd.DataFrame)
                                     else tuple(f"covariate_{i+1}" for i in range(cov.shape[1])))
        if cov.ndim != 2 or len(cov) != len(X) or not np.isfinite(cov).all():
            raise ValueError("Covariates must be a finite numeric matrix aligned with TRAIN")
        self.score_mean_ = scores.mean(axis=0)
        self.score_scale_ = scores.std(axis=0)
        self.score_scale_[self.score_scale_ <= 1e-12] = 1.0
        self.covariate_mean_ = cov.mean(axis=0)
        self.covariate_scale_ = cov.std(axis=0)
        self.covariate_scale_[self.covariate_scale_ <= 1e-12] = 1.0
        scaled_scores = (scores - self.score_mean_) / self.score_scale_
        scaled_cov = (cov - self.covariate_mean_) / self.covariate_scale_
        design = np.column_stack((scaled_scores, scaled_cov))
        ancestry = X[self.ancestry_col].astype(str).str.upper().to_numpy()
        self.anchor_models_, self.anchor_fit_ids_ = {}, {}
        records, centers, genetic_terms, included = [], [], [], []
        self.method_ = "disco"
        self.normalization_ids_ = self.fit_ids_
        for anchor in self.anchors:
            anchor = str(anchor).upper()
            selected = ancestry == anchor
            n = int(selected.sum())
            cases = int(y[selected].sum()) if self.task == "binary" else None
            r = 1.0 if self.quality is None else float(self.quality.get(anchor, 1.0))
            row = {"anchor": anchor, "n_train": n, "n_cases": cases,
                   "n_controls": n - cases if cases is not None else None,
                   "status": "omitted", "reason": "", "quality": r,
                   "n_discovery_scores": len(self.score_cols_)}
            records.append(row)
            if not np.isfinite(r) or r < 0:
                raise ValueError("Anchor quality must be finite and nonnegative")
            if r == 0:
                row["reason"] = "zero_quality_weight"
                continue
            if n < self.min_anchor_n:
                row["reason"] = "insufficient_training_n"
                continue
            if self.task == "binary" and min(cases, n - cases) < self.min_cases:
                row["reason"] = "insufficient_training_cases_or_controls"
                continue
            if np.std(y[selected]) <= 1e-12:
                row["reason"] = "constant_training_outcome"
                continue
            try:
                beta, model_type, n_iter = _fit_regression(
                    design[selected], y[selected], self.task, self.penalty, self.max_iter)
            except DiscoFitError as exc:
                row["reason"] = str(exc)
                continue
            genetic = scaled_scores @ beta[1:1 + len(self.score_cols_)]
            pc_normalizer = FrozenPCNormalizer().fit(pcs, genetic[:, None])
            if not bool(pc_normalizer.valid_[0]):
                row["reason"] = "no_genetic_variance_after_pc_adjustment"
                continue
            self.anchor_models_[anchor] = {
                "intercept": float(beta[0]), "score_coef": beta[1:1 + len(self.score_cols_)],
                "covariate_coef": beta[1 + len(self.score_cols_):],
                "regression": model_type, "iterations": n_iter}
            self.anchor_fit_ids_[anchor] = tuple(np.asarray(self.fit_ids_)[selected])
            centers.append(np.median(pcs[selected], axis=0))
            genetic_terms.append(genetic)
            included.append(anchor)
            row.update(status="included", reason="", regression=model_type,
                       iterations=n_iter, n_covariates=cov.shape[1])
        self.anchor_records_ = records
        self.anchor_names_ = tuple(included)
        if len(included) < 2:
            raise DiscoFitError("Fewer than two usable fine-tuning cohorts; inspect anchor_table()")
        self.centers_ = np.asarray(centers)
        quality = [next(row["quality"] for row in records if row["anchor"] == a) for a in included]
        self.geometry_ = DistanceInterpolator(self.centers_, quality)
        self.normalizer_ = FrozenPCNormalizer().fit(pcs, np.column_stack(genetic_terms))
        self.n_features_in_ = len(self.score_cols_)
        self.is_fitted_ = True
        return self

    def _check_fit(self):
        if not getattr(self, "is_fitted_", False):
            raise DiscoFitError("DiscoDivas has not been successfully fitted")

    def predict_components(self, X):
        self._check_fit()
        scores, pcs = _matrix(X, self.score_cols_), _matrix(X, self.pc_cols_)
        scaled = (scores - self.score_mean_) / self.score_scale_
        betas = np.column_stack([self.anchor_models_[a]["score_coef"] for a in self.anchor_names_])
        return self.normalizer_.transform(pcs, scaled @ betas)

    def interpolation_weights(self, X):
        self._check_fit()
        return self.geometry_.weights(_matrix(X, self.pc_cols_))

    def predict(self, X):
        components = self.predict_components(X)
        weights = self.interpolation_weights(X)
        return np.sum(components * weights, axis=1)

    def predict_details(self, X):
        self._check_fit()
        weights, status = self.geometry_.weights(_matrix(X, self.pc_cols_), return_status=True)
        components = self.predict_components(X)
        predictions = np.sum(components * weights, axis=1)
        missing_scores = ~np.isfinite(_matrix(X, self.score_cols_)).all(axis=1)
        status[missing_scores] = "missing_discovery_score"
        ids = X[self.id_col].astype(str).to_numpy() if self.id_col in X else X.index.astype(str)
        out = pd.DataFrame({self.id_col: ids, self.method_: predictions, "status": status})
        out["has_signed_weight"] = np.any(weights < 0, axis=1)
        for j, a in enumerate(self.anchor_names_):
            out[f"weight_{a}"] = weights[:, j]
            out[f"component_{a}"] = components[:, j]
        return out

    def anchor_table(self):
        out = pd.DataFrame(getattr(self, "anchor_records_", []))
        if not out.empty and hasattr(self, "geometry_"):
            for j, a in enumerate(self.anchor_names_):
                mask = out.anchor.eq(a)
                out.loc[mask, "distance_shrinkage_d"] = self.geometry_.d[j]
                for k, name in enumerate(self.pc_cols_):
                    out.loc[mask, "median_" + name] = self.centers_[j, k]
        return out

    def coefficient_table(self):
        rows = []
        for anchor, model in getattr(self, "anchor_models_", {}).items():
            rows.append({"anchor": anchor, "kind": "intercept", "term": "intercept",
                         "coefficient": model["intercept"], "used_in_interpolation": False})
            for name, value in zip(self.score_cols_, model["score_coef"]):
                rows.append({"anchor": anchor, "kind": "discovery_score", "term": name,
                             "coefficient": float(value), "used_in_interpolation": True})
            for name, value in zip(self.covariate_names_, model["covariate_coef"]):
                rows.append({"anchor": anchor, "kind": "covariate", "term": name,
                             "coefficient": float(value), "used_in_interpolation": False})
        return pd.DataFrame(rows)

    def summary(self):
        self._check_fit()
        return {"method": self.method_, "input_method": self.source_method,
                "description": f"DiscoDivas with {self.source_method}-input fine-tuned cohort models",
                "task": self.task, "n_train": len(self.fit_ids_),
                "score_columns": list(self.score_cols_), "pc_columns": list(self.pc_cols_),
                "active_anchors": list(self.anchor_names_), "penalty": self.penalty,
                "pc_geometry": "unscaled common reference-projected PC coordinates",
                "anchor_centers": "medians of actual fine-tuning TRAIN cohorts",
                "normalization": "TRAIN-frozen PC residualization and sample-SD scaling",
                "output": "genetic score, not calibrated phenotype/risk",
                "distance_condition_number": self.geometry_.condition_number,
                "signed_distance_shrinkage": bool(np.any(self.geometry_.d < 0)),
                "fit_uses_test_labels": False, "fit_roles": ["train"],
                "anchors": self.anchor_table().replace({np.nan: None}).to_dict("records")}


class RawScoreInterpolation(DiscoDivas):
    """Named ablation: interpolate raw ancestry scores, without phenotype fine-tuning.

    The explicit ``anchor_score_cols`` mapping prevents accidental assumptions
    about score/ancestry order. It must not be presented as the main Disco model.
    """

    def __init__(self, anchor_score_cols, pc_cols=DEFAULT_PCS, ancestry_col="ancestry",
                 id_col="id", role_col="role", min_anchor_n=100, quality=None,
                 source_method="PRS-CSx"):
        self.anchor_score_cols = anchor_score_cols
        self.pc_cols = pc_cols
        self.ancestry_col = ancestry_col
        self.id_col = id_col
        self.role_col = role_col
        self.min_anchor_n = min_anchor_n
        self.quality = quality
        self.source_method = source_method

    def fit(self, X, y=None, *, covariates=None):
        self.is_fitted_ = False
        for name in ("geometry_", "normalizer_"):
            self.__dict__.pop(name, None)
        self.anchor_records_, self.anchor_models_, self.anchor_fit_ids_ = [], {}, {}
        self.fit_ids_ = _training_ids(X, self.id_col, self.role_col)
        self.normalization_ids_ = self.fit_ids_
        mapping = {str(a).upper(): col for a, col in self.anchor_score_cols.items()}
        self.score_cols_, self.pc_cols_ = tuple(dict.fromkeys(mapping.values())), tuple(self.pc_cols)
        scores, pcs = _matrix(X, self.score_cols_), _matrix(X, self.pc_cols_)
        if not np.isfinite(scores).all() or not np.isfinite(pcs).all():
            raise ValueError("Ablation training scores and PCs must be complete")
        ancestry = X[self.ancestry_col].astype(str).str.upper().to_numpy()
        self.score_mean_, self.score_scale_ = np.zeros(scores.shape[1]), np.ones(scores.shape[1])
        self.covariate_names_, self.anchor_models_, self.anchor_fit_ids_ = (), {}, {}
        centers, terms, included, records = [], [], [], []
        for a, col in mapping.items():
            selected = ancestry == a
            n = int(selected.sum())
            r = 1.0 if self.quality is None else float(self.quality.get(a, 1.0))
            row = {"anchor": a, "n_train": n, "quality": r,
                   "status": "omitted", "reason": "", "raw_score": col}
            records.append(row)
            if not np.isfinite(r) or r < 0:
                raise ValueError("Anchor quality must be finite and nonnegative")
            if n < self.min_anchor_n or r == 0:
                row["reason"] = "insufficient_training_n" if n < self.min_anchor_n else "zero_quality_weight"
                continue
            j = self.score_cols_.index(col)
            if not FrozenPCNormalizer().fit(pcs, scores[:, j:j+1]).valid_[0]:
                row["reason"] = "no_genetic_variance_after_pc_adjustment"
                continue
            beta = np.zeros(len(self.score_cols_))
            beta[j] = 1.0
            self.anchor_models_[a] = {"intercept": 0.0, "score_coef": beta,
                                      "covariate_coef": np.zeros(0), "regression": "none"}
            # These IDs determine only the center; no outcome fitting is performed.
            self.anchor_fit_ids_[a] = ()
            centers.append(np.median(pcs[selected], axis=0))
            terms.append(scores[:, j])
            included.append(a)
            row.update(status="included", reason="raw_interpolation_ablation")
        self.anchor_records_, self.anchor_names_ = records, tuple(included)
        if len(included) < 2:
            raise DiscoFitError("Fewer than two available raw-score interpolation anchors")
        self.centers_ = np.asarray(centers)
        quality = [next(row["quality"] for row in records if row["anchor"] == a) for a in included]
        self.geometry_ = DistanceInterpolator(self.centers_, quality)
        self.normalizer_ = FrozenPCNormalizer().fit(pcs, np.column_stack(terms))
        self.method_, self.task, self.penalty = "disco.raw_interpolation", "none", 0.0
        self.n_features_in_, self.is_fitted_ = len(self.score_cols_), True
        return self

    def summary(self):
        out = super().summary()
        out["description"] = "Raw ancestral-score interpolation ablation; no fine-tuned cohort PRS models"
        return out


def _read_table(path, id_col="id"):
    separator = "," if str(path).endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(path, sep=separator, dtype={id_col: str})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    fit = commands.add_parser("fit", help="Fit from TRAIN rows of a roster/score table")
    fit.add_argument("--features", required=True, help="id,group_id,role,ancestry,covariates,PCs,scores TSV")
    fit.add_argument("--outcomes", required=True, help="Separate id,y TSV; only TRAIN IDs are selected")
    fit.add_argument("--outcome-col", default="y")
    fit.add_argument("--score-cols", required=True, help="Comma-separated discovery score column names")
    fit.add_argument("--task", choices=("continuous", "binary"), default="continuous")
    fit.add_argument("--pc-cols", default=",".join(DEFAULT_PCS))
    fit.add_argument("--covariate-cols", default=",".join(DEFAULT_COVARIATES))
    fit.add_argument("--min-anchor-n", type=int, default=100)
    fit.add_argument("--min-cases", type=int, default=20)
    fit.add_argument("--penalty", type=float, default=1e-6)
    fit.add_argument("--source-method", default="PRS-CSx")
    fit.add_argument("--model", required=True, help="Local-only fitted model pickle output")
    predict = commands.add_parser("predict", help="Predict genetic scores; no outcomes are read")
    predict.add_argument("--features", required=True)
    predict.add_argument("--model", required=True)
    predict.add_argument("--out", required=True, help="Local-only TSV/TSV.GZ including weights and components")
    args = parser.parse_args(argv)
    if args.command == "fit":
        frame = _read_table(args.features)
        if "role" not in frame:
            raise ValueError("CLI training requires an explicit role column")
        if frame.id.duplicated().any():
            raise ValueError("Duplicate roster IDs")
        if "group_id" in frame and frame.groupby("group_id", dropna=False).role.nunique().max() > 1:
            raise ValueError("A family group appears in multiple data roles")
        train = frame.loc[frame.role.eq("train")].copy()
        outcomes = _read_table(args.outcomes)
        outcomes = outcomes.loc[outcomes.id.isin(train.id), ["id", args.outcome_col]]
        if outcomes.id.duplicated().any():
            raise ValueError("Duplicate TRAIN outcome IDs")
        y = outcomes.set_index("id").reindex(train.id)[args.outcome_col].to_numpy(float)
        model = DiscoDivas(tuple(args.score_cols.split(",")), task=args.task,
                           pc_cols=tuple(args.pc_cols.split(",")),
                           covariate_cols=tuple(c for c in args.covariate_cols.split(",") if c),
                           min_anchor_n=args.min_anchor_n, min_cases=args.min_cases,
                           penalty=args.penalty, source_method=args.source_method)
        model.fit(train, y)
        target = Path(args.model)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("wb") as handle:
            pickle.dump(model, handle, protocol=pickle.HIGHEST_PROTOCOL)
        model.anchor_table().to_csv(target.with_suffix(".anchors.tsv"), sep="\t", index=False)
        model.coefficient_table().to_csv(target.with_suffix(".coefficients.tsv"), sep="\t", index=False)
        target.with_suffix(".json").write_text(json.dumps(model.summary(), indent=2) + "\n")
    else:
        with Path(args.model).open("rb") as handle:
            model = pickle.load(handle)
        frame = _read_table(args.features)
        target = Path(args.out)
        target.parent.mkdir(parents=True, exist_ok=True)
        model.predict_details(frame).to_csv(target, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    # Keep pickle classes importable as disco2.* even when fitting via this CLI.
    from disco2 import main as imported_main
    raise SystemExit(imported_main())
