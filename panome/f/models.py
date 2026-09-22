"""Classical comparators and strictly out-of-fold prototype quality."""
import inspect
import warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.exceptions import ConvergenceWarning
from sklearn.ensemble import HistGradientBoostingClassifier
from scipy.special import logit, expit
from scipy.optimize import minimize
from preprocess import MolecularPreprocessor
from survival import CensoringKM, loss, metrics
from common import words, log


def fit_logistic(x, y, w, c, ratio, seed, max_iter=2000):
    keep = w > 0
    if np.bincount(y[keep], minlength=2).min() < 5:
        raise ValueError("Too few observed cases/controls for a logistic fit")
    kw = dict(C=float(c), solver="saga" if ratio else "lbfgs", max_iter=max_iter,
              tol=1e-4, random_state=seed)
    # sklearn 1.8+ encodes the penalty by l1_ratio; earlier supported versions
    # require an explicit penalty. Do not emit thousands of deprecation warnings.
    if inspect.signature(LogisticRegression).parameters["penalty"].default == "deprecated":
        kw["l1_ratio"] = ratio
    else:
        kw["penalty"] = "elasticnet" if ratio else "l2"
        if ratio:
            kw["l1_ratio"] = ratio
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model = LogisticRegression(**kw).fit(
            np.asarray(x[keep], dtype=np.float64), y[keep],
            sample_weight=w[keep]/np.mean(w[keep]))
    if any(issubclass(v.category, ConvergenceWarning) for v in caught):
        raise ValueError(f"Logistic model did not converge in {max_iter} iterations")
    return model


def tune_logistic(x, y, w, build, tune, name, a, ratio):
    best = None
    rows = []
    for c in a.c_grid:
        try:
            model = fit_logistic(x[build], y[build], w[build], c, ratio, a.seed, a.max_iter)
            pred = model.predict_proba(x[tune])[:, 1]
            value = np.average(loss(y[tune], pred), weights=w[tune])
            rows.append(dict(model=name, C=c, validation_logloss=value, status="completed"))
            if best is None or value < best[0]:
                best = (value, model)
        except ValueError as exc:
            rows.append(dict(model=name, C=c, validation_logloss=np.nan, status=str(exc)))
    if best is None:
        raise ValueError(f"All candidates failed: {name}; see --max-iter and event counts")
    return best[1], rows


class RiskCalibrator:
    """Monotone logistic recalibration; optional additive clinical log-odds.

    Coefficients are nonnegative. Constant raw scores remain intercept-only.
    Fit only on the dedicated calibration subset; never on test outcomes.
    """
    def fit(self, raw, y, w, clinical=None):
        x = logit(np.clip(raw, 1e-5, 1-1e-5))[:, None]
        if clinical is not None:
            x = np.c_[x, logit(np.clip(clinical, 1e-5, 1-1e-5))]
        self.mean, self.sd = x.mean(0), x.std(0)
        self.sd[self.sd < 1e-8] = 1
        z = np.c_[np.ones(len(x)), (x-self.mean)/self.sd]
        wn = w/w.sum()
        def objective(b):
            eta = z @ b
            val = np.sum(wn*(np.logaddexp(0, eta)-y*eta)) + .0001*np.sum(b[1:]**2)
            grad = z.T @ (wn*(expit(eta)-y)) + np.r_[0, .0002*b[1:]]
            return val, grad
        prior = np.clip(np.average(y, weights=w), 1e-5, 1-1e-5)
        fit = minimize(objective, np.r_[logit(prior), np.ones(x.shape[1])], jac=True,
                       method="L-BFGS-B", bounds=[(None, None)]+[(0, None)]*x.shape[1])
        if not fit.success:
            raise ValueError("Risk recalibration failed: "+fit.message)
        self.coef = fit.x
        return self

    def predict(self, raw, clinical=None):
        x = logit(np.clip(raw, 1e-5, 1-1e-5))[:, None]
        if clinical is not None:
            x = np.c_[x, logit(np.clip(clinical, 1e-5, 1-1e-5))]
        return expit(self.coef[0]+((x-self.mean)/self.sd) @ self.coef[1:])


def oof_quality(raw, people, features, a, out):
    n, nf = raw.shape
    repeat_predictions = np.full((a.oof_repeats, n), np.nan)
    repeat_gains = np.full_like(repeat_predictions, np.nan)
    linear_predictions = np.full_like(repeat_predictions, np.nan)
    tree_predictions = np.full_like(repeat_predictions, np.nan)
    coefs, fit_rows = [], []
    ids = people[a.id_col].to_numpy(str)
    groups = people[a.group_col].to_numpy(str) if a.group_col else None
    # Unknown horizon outcomes form their own fold stratum.
    y_all = (people.event.eq(1) & people.time.le(a.horizon)).to_numpy(int)
    known = people.time.gt(a.horizon).to_numpy() | y_all.astype(bool)
    strata = np.where(known, y_all, 2)
    for repeat in range(a.oof_repeats):
        cv = (StratifiedGroupKFold(a.folds, shuffle=True, random_state=a.seed+repeat)
              if groups is not None else StratifiedKFold(a.folds, shuffle=True, random_state=a.seed+repeat))
        for fold, (tr, va) in enumerate(cv.split(raw, strata, groups)):
            prep = MolecularPreprocessor(a.feature_missing, words(a.residualize),
                                         words(a.categorical), a.transform).fit(raw[tr], people.iloc[tr])
            xt, _ = prep.transform(raw[tr], people.iloc[tr])
            xv, mv = prep.transform(raw[va], people.iloc[va])
            train_ok = np.mean(~np.isfinite(raw[tr][:, prep.keep]), axis=1) <= a.sample_missing
            km = CensoringKM().fit(people.iloc[tr[train_ok]], a.horizon, a.min_censor_survival)
            yt, wt = km.labels_weights(people.iloc[tr[train_ok]])
            yv, wv = km.labels_weights(people.iloc[va])
            model = fit_logistic(xt[train_ok], yt, wt, a.teacher_c, .5,
                                 a.seed+repeat*a.folds+fold, a.max_iter)
            linear_pred = model.predict_proba(xv)[:, 1]
            linear_predictions[repeat, va] = linear_pred
            pred = linear_pred
            if a.quality_teacher == "ensemble":
                tree = HistGradientBoostingClassifier(max_iter=a.quality_trees, max_leaf_nodes=7,
                    learning_rate=.05, min_samples_leaf=30, l2_regularization=10,
                    early_stopping=False, random_state=a.seed+repeat*a.folds+fold)
                known_train = wt > 0
                tree.fit(xt[train_ok][known_train], yt[known_train],
                         sample_weight=wt[known_train]/wt[known_train].mean())
                tree_pred = tree.predict_proba(xv)[:,1]
                tree_predictions[repeat,va] = tree_pred
                pred = .5*(linear_pred+tree_pred)
            prior = np.average(yt, weights=wt)
            gain = loss(yv, prior)-loss(yv, pred)
            valid = (wv > 0) & (1-mv.mean(1) <= a.sample_missing)
            repeat_predictions[repeat, va] = pred
            repeat_gains[repeat, va[valid]] = gain[valid]
            beta = np.zeros(nf)
            beta[prep.keep] = model.coef_[0]
            coefs.append(beta)
            fit_rows.append(dict(repeat=repeat+1, fold=fold+1, train_n=len(tr), validation_n=len(va),
                                 retained_features=len(prep.keep), iterations=int(model.n_iter_[0])))
        log("DONE", "OOF quality", f"repeat={repeat+1}/{a.oof_repeats}")
    finite = np.all(np.isfinite(repeat_gains), axis=0)
    mean_gain = np.full(n, np.nan)
    sd_gain = np.full(n, np.nan)
    mean_gain[finite] = repeat_gains[:, finite].mean(0)
    sd_gain[finite] = repeat_gains[:, finite].std(0)
    fraction = (repeat_gains > 0).mean(0)
    reliable = finite & (mean_gain > a.min_gain) & (fraction >= a.fit_stability)
    table = pd.DataFrame({a.id_col: ids, "horizon_label": y_all, "known_label": known,
                          "OOF_probability": repeat_predictions.mean(0),
                          "OOF_probability_SD": repeat_predictions.std(0),
                          "OOF_elasticnet_probability": linear_predictions.mean(0),
                          "OOF_tree_probability": tree_predictions.mean(0),
                          "OOF_logloss": loss(y_all, repeat_predictions.mean(0)),
                          "fit_gain": mean_gain, "fit_gain_SD": sd_gain,
                          "positive_gain_fraction": fraction, "reliable_candidate": reliable})
    table.to_csv(out/"reference_candidates.csv", index=False)
    pd.DataFrame(fit_rows).to_csv(out/"oof_fits.csv", index=False)
    beta = np.asarray(coefs)
    magnitude = np.mean(np.abs(beta), axis=0)
    stability = np.abs(np.mean(np.sign(beta), axis=0))
    importance = magnitude*stability
    pd.DataFrame(dict(feature=features, mean_abs_beta=magnitude,
                      sign_stability=stability, metric_importance=importance)).to_csv(
                          out/"metric_features.csv", index=False)
    return table, importance
