"""Same-input tabular baselines, trained without test-set model selection.

The caller supplies features derived only from the prediction-time prefix.
This module never constructs features from full follow-up records. Training
statistics fit the imputer/scaler. Validation labels choose prespecified
hyperparameters; test labels are neither accepted nor used by this API.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import joblib
import numpy as np
from scipy.optimize import minimize
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from threadpoolctl import threadpool_limits


@dataclass
class TabularPreprocessor:
    median_: np.ndarray | None = None
    scale_: np.ndarray | None = None

    def fit(self, features):
        x = np.asarray(features)
        if x.ndim != 2 or not len(x):
            raise ValueError("Training features must be a nonempty 2D matrix")
        median = np.empty(x.shape[1], dtype=np.float32)
        scale = np.empty_like(median)
        # Column-wise work avoids a second dense float64 cohort matrix.
        for k in range(x.shape[1]):
            column = np.array(x[:, k], dtype=float, copy=True)
            observed = np.isfinite(column)
            center = float(np.median(column[observed])) if observed.any() else 0.0
            column[~observed] = center
            spread = float(np.std(column))
            median[k], scale[k] = center, spread if np.isfinite(spread) and spread > 1e-8 else 1.0
        self.median_, self.scale_ = median, scale
        return self

    def transform(self, features):
        if self.median_ is None:
            raise RuntimeError("Preprocessor has not been fitted")
        x = np.asarray(features)
        if x.ndim != 2 or x.shape[1] != len(self.median_):
            raise ValueError("Feature matrix does not match fitted preprocessor")
        out = np.empty(x.shape, dtype=np.float32)
        for start in range(0, len(x), 4096):
            block = np.array(x[start:start + 4096], dtype=np.float32, copy=True)
            block = np.where(np.isfinite(block), block, self.median_)
            out[start:start + len(block)] = (block - self.median_) / self.scale_
        return out


class RidgeCoxPH:
    """Breslow-tie Cox proportional-hazards model with ridge regularization.

    The likelihood gradient uses a scalar reverse cumulative sum followed by
    one matrix product; it does not allocate a cumulative n-by-p matrix. The
    fit is an ordinary Cox model, not a claim to reproduce DeepSurv.
    """
    def __init__(self, alpha=1.0, max_iter=150, tolerance=1e-7):
        self.alpha, self.max_iter, self.tolerance = float(alpha), int(max_iter), float(tolerance)
        if not np.isfinite(self.alpha) or self.alpha < 0 or self.max_iter < 1 or self.tolerance <= 0:
            raise ValueError("Invalid Cox regularization or optimization settings")

    @staticmethod
    def _objective(beta, x, event, ends, event_counts, alpha):
        eta = x @ beta
        shift = float(np.max(eta))
        risk = np.exp(eta - shift)
        denominators = np.cumsum(risk)[ends]
        if np.any(denominators <= 0):
            return np.inf, np.full_like(beta, np.nan)
        count = max(int(event.sum()), 1)
        value = (np.dot(event_counts, np.log(denominators) + shift) - np.dot(event, eta)) / count
        increments = np.zeros(len(event), dtype=float)
        increments[ends] = event_counts / denominators
        multiplier = np.cumsum(increments[::-1])[::-1]
        gradient = x.T @ (risk * multiplier - event) / count
        value += 0.5 * alpha * float(np.dot(beta, beta))
        gradient += alpha * beta
        return float(value), np.asarray(gradient)

    def fit(self, features, observed_time, event):
        x = np.asarray(features, dtype=np.float32)
        t, e = np.asarray(observed_time, dtype=float), np.asarray(event, dtype=bool)
        if x.ndim != 2 or len(x) != len(t) or t.shape != e.shape or t.ndim != 1:
            raise ValueError("Cox feature/time/event arrays have inconsistent shapes")
        valid = np.isfinite(t) & (t > 0) & np.all(np.isfinite(x), axis=1)
        x, t, e = x[valid], t[valid], e[valid]
        if not len(x) or e.sum() == 0:
            raise ValueError("Cox fitting requires at least one event")
        order = np.argsort(-t, kind="stable")
        x, t, e = x[order], t[order], e[order]
        starts = np.r_[0, np.flatnonzero(t[1:] != t[:-1]) + 1]
        ends = np.r_[starts[1:] - 1, len(t) - 1]
        event_counts = np.add.reduceat(e.astype(float), starts)
        result = minimize(self._objective, np.zeros(x.shape[1]),
                          args=(x, e.astype(float), ends, event_counts, self.alpha),
                          method="L-BFGS-B", jac=True,
                          options={"maxiter": self.max_iter, "ftol": self.tolerance})
        if not np.isfinite(result.fun) or not np.all(np.isfinite(result.x)):
            raise ValueError("Cox optimization produced nonfinite parameters")
        self.coef_ = result.x
        self.converged_ = bool(result.success)
        self.optimization_message_ = str(result.message)
        self.n_iter_ = int(result.nit)
        eta = x @ self.coef_
        shift = float(np.max(eta))
        denominators = np.cumsum(np.exp(eta - shift))[ends]
        have_event = event_counts > 0
        event_times = t[ends][have_event]
        # Keep the baseline in the shifted risk scale for numerical stability.
        increments = event_counts[have_event] / denominators[have_event]
        ascending = np.argsort(event_times)
        self.event_times_ = event_times[ascending]
        self.baseline_cumulative_ = np.cumsum(increments[ascending])
        self.risk_shift_ = shift
        return self

    def predict(self, features):
        return np.asarray(features) @ self.coef_

    def predict_proba(self, features, horizon_years):
        idx = np.searchsorted(self.event_times_, float(horizon_years), side="right") - 1
        if idx < 0:
            return np.zeros(len(features), dtype=float)
        log_ch = np.log(self.baseline_cumulative_[idx]) + self.predict(features) - self.risk_shift_
        return -np.expm1(-np.exp(np.clip(log_ch, -745, 30)))


def _matrix(value, name, n_rows=None):
    x = np.asarray(value)
    if x.ndim == 1:
        x = x[:, None]
    if x.ndim != 2 or (n_rows is not None and len(x) != n_rows):
        raise ValueError(f"{name} must be a compatible 2D matrix")
    return x


def _auc(labels, predictions, mask):
    keep = np.asarray(mask, dtype=bool) & np.isfinite(predictions) & np.isfinite(labels)
    y = np.asarray(labels)[keep].astype(bool)
    return float(roc_auc_score(y, predictions[keep])) if y.any() and (~y).any() else float("nan")


def _safe_name(value):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip(".")
    return name or "endpoint"


def _prediction_array(shape, output_dir, name):
    if output_dir is None:
        result = np.empty(shape, dtype=np.float32)
    else:
        result = np.lib.format.open_memmap(Path(output_dir) / name, mode="w+", dtype=np.float32, shape=shape)
    result[:] = np.nan
    return result


def fit_predict_baselines(X_train, labels_train, masks_train,
                          X_val, labels_val, masks_val, X_test, *,
                          methods=("logistic", "histgb"), seed=2026,
                          endpoint_names=None, survival_train=None,
                          survival_val=None, horizon_years=5.0,
                          min_cases=10, max_iter=200, n_jobs=1,
                          output_dir=None, options=None, return_models=False):
    """Fit per-endpoint classifiers or Cox models and return held-out scores.

    Parameters
    ----------
    labels_*, masks_* : (people, endpoints) arrays
        Endpoint-specific labels and eligibility masks prepared by the caller.
        For the recorded task, short-follow-up controls remain controls. A
        complete-follow-up sensitivity uses different masks and a named run.
    survival_train, survival_val : tuple(time, event), optional
        Each array has shape (people, endpoints). Required for method='cox'.
        The Cox fit administratively truncates these observations at the
        requested prediction horizon.
    options : dict, optional
        logistic: {'C': [0.1, 1.0]}; histgb: {'max_leaf_nodes': [15, 31],
        'learning_rate': 0.05, 'l2_regularization': 1.0};
        lightgbm: {'num_leaves': [15,31], 'learning_rate': 0.05};
        cox: {'alpha': [0.01, 0.1, 1.0]}.
    output_dir : path, optional
        Store the fitted training preprocessor and endpoint estimators using
        joblib. Only load these trusted locally generated files. Predicted
        score matrices become memory-mapped NPY files, avoiding E-fold copies.

    No test outcomes are accepted. No model is refitted on validation data.
    All output probabilities refer only to horizon_years.
    """
    from metrics import harrell_c_index
    methods = tuple(methods)
    allowed = {"logistic", "histgb", "lightgbm", "cox"}
    if not methods or len(set(methods)) != len(methods) or set(methods) - allowed:
        raise ValueError("methods must be unique names from logistic,histgb,lightgbm,cox")
    if not np.isfinite(horizon_years) or horizon_years <= 0 or min_cases < 1 or max_iter < 1 or n_jobs < 1:
        raise ValueError("Invalid baseline horizon or resource settings")
    x_train, x_val, x_test = map(np.asarray, (X_train, X_val, X_test))
    if any(x.ndim != 2 for x in (x_train, x_val, x_test)) or len({x.shape[1] for x in (x_train, x_val, x_test)}) != 1:
        raise ValueError("All baseline feature matrices must have identical columns")
    y_train = _matrix(labels_train, "labels_train", len(x_train))
    y_val = _matrix(labels_val, "labels_val", len(x_val))
    m_train = _matrix(masks_train, "masks_train", len(x_train)).astype(bool)
    m_val = _matrix(masks_val, "masks_val", len(x_val)).astype(bool)
    if y_train.shape != m_train.shape or y_val.shape != m_val.shape or y_train.shape[1] != y_val.shape[1]:
        raise ValueError("Baseline label and mask arrays do not align")
    endpoints = list(endpoint_names) if endpoint_names is not None else [str(i) for i in range(y_train.shape[1])]
    if len(endpoints) != y_train.shape[1] or len(set(map(str, endpoints))) != len(endpoints):
        raise ValueError("endpoint_names must be unique and match the target columns")
    for labels, mask in ((y_train, m_train), (y_val, m_val)):
        selected = labels[mask & np.isfinite(labels)]
        if np.any((selected != 0) & (selected != 1)):
            raise ValueError("Classifier labels must be 0 or 1 on eligible rows")
    if "cox" in methods and (survival_train is None or survival_val is None):
        raise ValueError("Cox requires survival_train and survival_val")
    if "lightgbm" in methods:
        try:
            from lightgbm import LGBMClassifier
        except ImportError as exc:
            raise RuntimeError("Requested LightGBM is not installed. Run flux/install.sh or omit lightgbm.") from exc
    else:
        LGBMClassifier = None
    st, sv = None, None
    if "cox" in methods:
        st = tuple(_matrix(v, "survival_train", len(x_train)) for v in survival_train)
        sv = tuple(_matrix(v, "survival_val", len(x_val)) for v in survival_val)
        if len(st) != 2 or len(sv) != 2 or any(v.shape != y_train.shape for v in st) or any(v.shape != y_val.shape for v in sv):
            raise ValueError("Cox survival arrays must match endpoint target matrices")
    directory = Path(output_dir) if output_dir is not None else None
    if directory is not None:
        directory.mkdir(parents=True, exist_ok=True)
    preprocessor = TabularPreprocessor().fit(x_train)
    x_train = preprocessor.transform(x_train)
    x_val = preprocessor.transform(x_val)
    x_test = preprocessor.transform(x_test)
    if directory is not None:
        joblib.dump(preprocessor, directory / "preprocess.joblib", compress=3)
    options = options or {}
    result = {"predictions": {}, "status": [], "preprocessor": preprocessor}
    if return_models:
        result["models"] = {}
    for method in methods:
        pv = _prediction_array((len(x_val), len(endpoints)), directory, f"{method}.validation.npy")
        pt = _prediction_array((len(x_test), len(endpoints)), directory, f"{method}.test.npy")
        result["predictions"][method] = {"validation": pv, "test": pt,
                                           "score_kind": "probability", "horizon_years": float(horizon_years)}
        if return_models:
            result["models"][method] = {}
        settings = options.get(method, {})
        for k, endpoint in enumerate(endpoints):
            train_mask = m_train[:, k] & np.isfinite(y_train[:, k])
            val_mask = m_val[:, k] & np.isfinite(y_val[:, k])
            if method == "cox":
                train_time, train_event = st[0][:, k].astype(float), st[1][:, k].astype(bool)
                val_time, val_event = sv[0][:, k].astype(float), sv[1][:, k].astype(bool)
                train_mask &= np.isfinite(train_time) & (train_time > 0)
                val_mask &= np.isfinite(val_time) & (val_time > 0)
                train_event = train_event & (train_time <= horizon_years)
                train_time = np.minimum(train_time, horizon_years)
                val_event = val_event & (val_time <= horizon_years)
                val_time = np.minimum(val_time, horizon_years)
                cases = int(train_event[train_mask].sum())
            else:
                cases = int(y_train[train_mask, k].sum())
            controls = int(train_mask.sum()) - cases
            status = {"method": method, "endpoint": str(endpoint), "train_n": int(train_mask.sum()),
                      "train_cases": cases, "train_controls": controls, "validation_n": int(val_mask.sum()),
                      "status": "ok", "selected_parameter": None, "validation_metric": float("nan"),
                      "selection": "validation", "model_path": ""}
            if cases < min_cases or controls < min_cases:
                status["status"] = "insufficient_training_cases_or_controls"
                result["status"].append(status)
                continue
            xx, yy = x_train[train_mask], y_train[train_mask, k].astype(np.int32)
            if method == "logistic":
                parameter_name, candidates = "C", settings.get("C", [0.1, 1.0])
            elif method == "histgb":
                parameter_name, candidates = "max_leaf_nodes", settings.get("max_leaf_nodes", [15, 31])
            elif method == "lightgbm":
                parameter_name, candidates = "num_leaves", settings.get("num_leaves", [15, 31])
            else:
                parameter_name, candidates = "alpha", settings.get("alpha", [0.01, 0.1, 1.0])
            if np.isscalar(candidates):
                candidates = [candidates]
            if not len(candidates):
                raise ValueError(f"No baseline candidate values for {method}")
            best_model, best_prediction, best_value, best_parameter = None, None, -np.inf, None
            try:
                with threadpool_limits(limits=n_jobs):
                    for candidate in candidates:
                        if method == "logistic":
                            estimator = LogisticRegression(C=float(candidate), solver="lbfgs", max_iter=max_iter,
                                                           random_state=seed)
                            estimator.fit(xx, yy)
                            prediction = estimator.predict_proba(x_val)[:, 1] if len(x_val) else np.array([])
                            metric = _auc(y_val[:, k], prediction, val_mask)
                        elif method == "histgb":
                            estimator = HistGradientBoostingClassifier(max_leaf_nodes=int(candidate), max_iter=max_iter,
                                        learning_rate=float(settings.get("learning_rate", 0.05)),
                                        l2_regularization=float(settings.get("l2_regularization", 1.0)),
                                        early_stopping=False, random_state=seed)
                            estimator.fit(xx, yy)
                            prediction = estimator.predict_proba(x_val)[:, 1] if len(x_val) else np.array([])
                            metric = _auc(y_val[:, k], prediction, val_mask)
                        elif method == "lightgbm":
                            estimator = LGBMClassifier(num_leaves=int(candidate), n_estimators=max_iter,
                                        learning_rate=float(settings.get("learning_rate", 0.05)),
                                        reg_lambda=float(settings.get("reg_lambda", 1.0)), n_jobs=n_jobs,
                                        random_state=seed, verbosity=-1, deterministic=True, force_col_wise=True)
                            estimator.fit(xx, yy)
                            prediction = estimator.predict_proba(x_val)[:, 1] if len(x_val) else np.array([])
                            metric = _auc(y_val[:, k], prediction, val_mask)
                        else:
                            estimator = RidgeCoxPH(alpha=float(candidate), max_iter=max_iter)
                            estimator.fit(xx, train_time[train_mask], train_event[train_mask])
                            prediction = estimator.predict_proba(x_val, horizon_years)
                            metric = harrell_c_index(val_time[val_mask], val_event[val_mask], estimator.predict(x_val[val_mask]))
                        selection_value = metric if np.isfinite(metric) else -np.inf
                        if best_model is None or selection_value > best_value:
                            best_model, best_prediction, best_value, best_parameter = estimator, prediction, selection_value, candidate
                    pv[:, k] = best_prediction
                    if len(x_test):
                        pt[:, k] = (best_model.predict_proba(x_test, horizon_years) if method == "cox"
                                    else best_model.predict_proba(x_test)[:, 1])
                status["selected_parameter"] = f"{parameter_name}={best_parameter}"
                status["validation_metric"] = float(best_value) if np.isfinite(best_value) else float("nan")
                if not np.isfinite(best_value):
                    status["selection"] = "first_prespecified_candidate_no_evaluable_validation"
                if method == "cox" and not best_model.converged_:
                    status["status"] = "optimizer_not_converged"
                elif method == "logistic" and np.max(best_model.n_iter_) >= max_iter:
                    status["status"] = "optimizer_iteration_limit"
                if directory is not None:
                    model_path = directory / f"{method}.{k:04d}.{_safe_name(endpoint)}.joblib"
                    joblib.dump({"estimator": best_model, "endpoint": str(endpoint), "horizon_years": float(horizon_years),
                                 "method": method, "n_features": x_train.shape[1]}, model_path, compress=3)
                    status["model_path"] = str(model_path)
                if return_models:
                    result["models"][method][str(endpoint)] = best_model
            except (ValueError, FloatingPointError) as exc:
                # Keep unsupported endpoints explicit rather than fabricating a
                # constant-risk fit or excluding them without an audit row.
                pv[:, k], pt[:, k] = np.nan, np.nan
                status["status"] = "fit_failed"
                status["error"] = str(exc)
            result["status"].append(status)
            del xx, yy
        for array in (pv, pt):
            if hasattr(array, "flush"):
                array.flush()
    return result
