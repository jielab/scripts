"""Shared incident-risk fitting, individual error assessment and historical transfer.

No function in this module accepts test outcomes. Development audit outcomes
are not used by fit either. New models share the input roster and a build-only
censoring estimator; source-domain historical labels never become target labels.
"""
from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace
import json
import subprocess
import shutil
import warnings

import joblib
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import log_loss

import abm

VERSION = "LE8-C1-unified-models-2026-10-22-v1"
HERE = Path(__file__).resolve().parent


def context(development, test_x, config):
    """One numerical contract. The only outcomes accepted are development labels."""
    config.validate()
    d, t = development, test_x
    names = list(map(str, d["feature_names"]))
    if names != list(map(str, t["feature_names"])):
        raise ValueError("Development/test feature schema mismatch")
    if set(d["role"].astype(str)) != set(abm.ROLES):
        raise ValueError("Exactly five development roles are required")
    X = np.concatenate([d["X"], t["X"]])
    ids = np.concatenate([d["eid"], t["eid"]]).astype(str)
    groups = np.concatenate([d["group"], t["group"]]).astype(str)
    roles = np.concatenate([d["role"], np.repeat("test", len(t["eid"]))]).astype(str)
    abm.TwinABM.validate_features(X, ids, groups)
    rg = pd.DataFrame({"g": groups, "r": roles}).groupby("g").r.nunique()
    if (rg > 1).any():
        raise ValueError("A family crosses fitting/tuning/calibration/test roles")
    forbidden = {"time", "event", "outcome", "target", "prevalent_status", "role", "split", "diagnosis_date"}
    if forbidden & {n.lower() for n in names} or any(n.lower().endswith((".t2e", ".yt2e", ".bi2e")) for n in names):
        raise ValueError("Outcome metadata in predictors")
    if len(names) != X.shape[1] or len(set(names)) != len(names):
        raise ValueError("Feature names must be unique and aligned")
    parts = {r: np.flatnonzero(roles == r) for r in (*abm.ROLES, "test")}
    fit_rows = np.concatenate([parts[r] for r in abm.ROLES[:4]])
    dt, de = np.asarray(d["time"], float), np.asarray(d["event"])
    if dt.shape != (len(d["eid"]),) or de.shape != dt.shape:
        raise ValueError("Development outcomes must align with the feature roster")
    if (not np.isfinite(dt[fit_rows]).all() or (dt[fit_rows] <= 0).any()
            or not np.isin(de[fit_rows], [0, 1]).all()):
        raise ValueError("Invalid fitting-role survival outcomes")
    b = parts["build"]
    prep = abm.NumericPrep().fit(X[b], config)
    qc = prep.technical_qc(X)
    km = abm.CensoringKM().fit(d["time"][b][qc[b]], d["event"][b][qc[b]],
                             config.horizon, config.min_censor_survival)
    y, w = np.zeros(len(X), int), np.zeros(len(X))
    y[fit_rows], w[fit_rows] = km.labels_weights(d["time"][fit_rows], d["event"][fit_rows])
    w *= qc
    for r in ("build", "tune_model", "calibration_fit"):
        abm.check_classes(y[parts[r]], w[parts[r]], config.min_events)
    z = prep.transform(X)
    return SimpleNamespace(development=d, test_x=t, X=X, ids=ids, groups=groups, roles=roles,
        parts=parts, names=names, cfg=config, prep=prep, km=km, y=y, w=w, qc=qc, z=z,
        common_support=None, common_geometry=None, support_cache={})


def shared_support(c):
    """Outcome-blind geometry; only labeled, technically valid build references.

    The binary outcome is used by TwinBank's risk estimator but the returned
    global-model support features depend only on geometry, weights and families.
    The returned label_disagreement and donor_loss are deliberately omitted.
    """
    # Shallow method contexts share this cache, while changing scientifically
    # relevant support parameters obtains a fresh geometry/radius/diagnostic set.
    cache_key = tuple(getattr(c.cfg, name, 0) for name in (
        "k", "radius_k_ref", "dimensions", "seed", "min_twins", "min_family_ess",
        "radius_quantile", "prior_strength")) + (tuple(c.cfg.always_include_indices),)
    cache = getattr(c, "support_cache", None)
    if cache is None:
        cache = c.support_cache = {}
    if cache_key in cache:
        c.common_support, c.common_geometry, c.common_bank = cache[cache_key]
        return c.common_support
    b, tune = c.parts["build"], c.parts["tune_model"]
    b = b[c.w[b] > 0]
    columns = np.flatnonzero(~np.isin(c.prep.keep, c.cfg.always_include_indices))
    cfg = replace(c.cfg, geometry="all_molecular_pca")
    geometry = abm.FullProfileGeometry().fit(c.z[c.parts["build"]][c.qc[c.parts["build"]]], columns, cfg)
    bank = abm.TwinBank().fit(c.z[b], c.y[b], c.w[b], c.ids[b], c.groups[b], np.zeros(len(b)), cfg, geometry)
    bank.calibrate_radius(c.z[tune][c.qc[tune]], c.ids[tune][c.qc[tune]], c.groups[tune][c.qc[tune]])
    ix = np.concatenate([c.parts[r] for r in ("tune_gate", "calibration_fit", "calibration_audit", "test")])
    _, local = bank.predict(c.z[ix], c.ids[ix], c.groups[ix], np.full(len(ix), .5))
    keys = ("eligible_match_count", "matched_count", "eligible_reference_count", "twin_fraction",
            "family_ess", "donor_ess", "matching_radius", "radius_ratio", "nearest_radius_ratio", "supported")
    diagnostics = {}
    for key in keys:
        diagnostics[key] = np.zeros(len(c.ids), dtype=bool if key == "supported" else float)
        diagnostics[key][ix] = local[key]
    diagnostics["technical_qc"] = c.qc.copy()
    diagnostics["supported"] &= c.qc
    c.common_support = diagnostics
    c.common_geometry = geometry
    c.common_bank = bank
    cache[cache_key] = (diagnostics, geometry, bank)
    return diagnostics


def _subset(values, ix):
    return {k: np.asarray(v)[ix] for k, v in values.items() if np.asarray(v).ndim == 1}


class IndividualAssessment:
    """A model-specific expected-error estimate, not probability of being correct.

    Error heads use tune_gate outcomes. Risk-bin mean loss is subtracted before
    fitting the primary ranking head, so low predicted risk alone is not the
    definition of reliability. Selection thresholds use calibration_fit features
    and frozen predictions. This class does not read calibration_fit labels,
    but upstream probability calibrators used those labels. Thus calibration_fit
    is a reused development role, not an independent validation of the error
    head or a distribution-free coverage guarantee. calibration_audit/test are
    untouched until all predictions and selection rules have been frozen.
    """
    def __init__(self, config):
        self.config = config

    @staticmethod
    def features(p, diagnostics):
        p = np.clip(np.asarray(p, float), 1e-7, 1-1e-7)
        n = len(p)
        get = lambda k, default=0: np.asarray(diagnostics.get(k, np.full(n, default)), float)
        return np.c_[p, p * (1-p), np.log1p(get("eligible_match_count")),
                     np.log1p(get("family_ess")), get("twin_fraction"),
                     np.minimum(get("radius_ratio", 100), 100),
                     np.minimum(get("nearest_radius_ratio", 100), 100)]

    def fit(self, gate_probability, gate_diagnostics, gate_y, gate_w,
            calibration_probability, calibration_diagnostics, calibration_ids):
        p, y, w = map(np.asarray, (gate_probability, gate_y, gate_w))
        if p.ndim != 1 or y.shape != p.shape or w.shape != p.shape:
            raise ValueError("Error-assessment gate inputs must be aligned vectors")
        if (not np.isfinite(p).all() or np.any((p < 0) | (p > 1))
                or not np.isfinite(w).all() or np.any(w < 0)
                or not np.isin(y[w > 0], [0, 1]).all()):
            raise ValueError("Invalid gate probabilities, labels or weights")
        qc = np.asarray(gate_diagnostics["technical_qc"], bool)
        known = qc & (w > 0) & np.isfinite(p)
        self.risk_edges = np.unique(np.quantile(p[known], np.linspace(0, 1, 11)))[1:-1] if known.any() else np.arange(.1, 1, .1)
        self.models = {}
        self.status = "support_only_insufficient_gate_cases"
        self.gate_N = int(known.sum())
        self.gate_cases = int(np.sum(y[known]))
        if known.sum() >= 50 and min(np.sum(y[known] == 1), np.sum(y[known] == 0)) >= self.config.min_events:
            loss = (y-p)**2
            bins = np.searchsorted(self.risk_edges, p)
            residual = loss.copy()
            self.risk_means = {}
            for label in np.unique(bins[known]):
                use = known & (bins == label)
                avg = float(np.average(loss[use], weights=w[use]))
                self.risk_means[int(label)] = avg
                residual[use] -= avg
            gx = self.features(p, gate_diagnostics)
            for target, name in ((loss, "error"), (residual, "residual")):
                head = HistGradientBoostingRegressor(max_iter=self.config.gate_iterations,
                    max_leaf_nodes=5, learning_rate=.04, min_samples_leaf=max(15, min(100, int(known.sum()/12))),
                    l2_regularization=10, early_stopping=False, random_state=self.config.seed + 503)
                head.fit(gx[known], target[known], sample_weight=w[known]/w[known].mean())
                self.models[name] = head
            self.status = "fitted_model_specific_risk_residual_error"
        cp = np.asarray(calibration_probability)
        if (cp.ndim != 1 or len(cp) != len(calibration_ids) or not np.isfinite(cp).all()
                or np.any((cp < 0) | (cp > 1))):
            raise ValueError("Invalid calibration-feature probability roster")
        _, score = self.raw_scores(cp, calibration_diagnostics)
        cq = np.asarray(calibration_diagnostics["technical_qc"], bool)
        supported = np.asarray(calibration_diagnostics.get("supported", cq), bool)
        allowed = cq & supported & np.isfinite(score)
        self.rules = {float(q): abm.FrozenRule().fit(score, allowed, np.asarray(calibration_ids).astype(str), q,
                        self.config.seed+1601) for q in self.config.coverages}
        self.score_reference = np.sort(score[cq & np.isfinite(score)])
        return self

    def raw_scores(self, p, d):
        x = self.features(p, d)
        if self.models:
            error = np.clip(self.models["error"].predict(x), 0, 1)
            score = -self.models["residual"].predict(x)
        else:
            error = np.full(len(p), np.nan)
            score = np.log1p(np.asarray(d.get("eligible_match_count", np.zeros(len(p)))))
        return error, score

    def project(self, p, diagnostics, ids, groups):
        error, score = self.raw_scores(p, diagnostics)
        qc = np.asarray(diagnostics["technical_qc"], bool)
        supported = qc & np.asarray(diagnostics.get("supported", qc), bool)
        masks = {q: rule.apply(score, supported, np.asarray(ids).astype(str)) for q, rule in self.rules.items()}
        return dict(eid=np.asarray(ids).astype(str), group=np.asarray(groups).astype(str),
            probability=np.asarray(p), technical_qc=qc, supported=supported,
            masks=masks, predicted_error=error, predictability=1-error,
            residual_predictability=score, predictability_scale="one_minus_expected_Brier_error",
            predictability_is_correctness_probability=False,
            diagnostics=diagnostics, risk_edges=self.risk_edges,
            risk_edges_source="tune_gate_frozen_before_evaluation", assessment_status=self.status)


def _fit_assessment(c, probability, diagnostics):
    g, cal = c.parts["tune_gate"], c.parts["calibration_fit"]
    return IndividualAssessment(c.cfg).fit(probability[g], _subset(diagnostics, g), c.y[g], c.w[g],
        probability[cal], _subset(diagnostics, cal), c.ids[cal])


def finish_projection(c, probability, diagnostics, out, *, fitted, metadata, extra_probabilities=None):
    """Freeze everything before the caller is allowed to load test truth."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    probability = np.asarray(probability, float)
    if probability.shape != (len(c.ids),) or not np.isfinite(probability).all():
        raise ValueError("Every person must receive a finite probability (fallback required)")
    if np.any((probability < 0) | (probability > 1)):
        raise ValueError("Invalid probability")
    assessment = _fit_assessment(c, probability, diagnostics)
    fitted.update(assessment=assessment, feature_names=c.names, config=asdict(c.cfg), version=VERSION)
    joblib.dump(fitted, out / "model_LOCAL_ONLY.joblib", compress=3)
    manifest = {**metadata, "version": VERSION, "config": asdict(c.cfg),
        "roles": {k: {"N": len(v), "id_hash": abm.fingerprint(c.ids[v])} for k,v in c.parts.items()},
        "test_outcomes_used": False, "assessment_status": assessment.status,
        "assessment_gate_N": assessment.gate_N, "assessment_gate_cases": assessment.gate_cases,
        "predictability_definition": "one minus predicted conditional squared probability error; not correctness probability",
        "assessment_calibration_role": "calibration_fit reused for upstream probability calibration and feature-only selection thresholds; validation is calibration_audit/test",
        "selection": "frozen outcome-blind test queries; unsupported persons retain risk predictions",
        "test_provenance": "previously_inspected_or_unverified_exploratory"}
    abm.json_write(out / "MODEL_FROZEN.json", manifest)
    saved = {}
    for role in ("calibration_audit", "test"):
        ix = c.parts[role]
        p = assessment.project(probability[ix], _subset(diagnostics, ix), c.ids[ix], c.groups[ix])
        p["metadata"] = metadata
        if extra_probabilities:
            p["probability_arms"] = {name: np.asarray(values)[ix] for name, values in extra_probabilities.items()}
        joblib.dump(p, out / f"{role}_predictions_LOCAL_ONLY.joblib", compress=3)
        frame = {**p["diagnostics"], "eid": p["eid"], "group": p["group"],
            "probability": p["probability"], "technical_qc": p["technical_qc"],
            "supported": p["supported"], "expected_Brier_error": p["predicted_error"],
            "predictability": p["predictability"],
            "risk_residual_predictability": p["residual_predictability"],
            **{f"selected_{q:g}": mask for q,mask in p["masks"].items()}}
        pd.DataFrame(frame).to_csv(out / f"{role}_predictions_LOCAL_ONLY.tsv.gz", sep="\t", index=False)
        saved[role] = p
    return saved


def fit_glmnet(c, out, r_bin="Rscript"):
    """Actual R glmnet binomial LASSO; never relabel a sklearn surrogate."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    work = out / "glmnet_input_LOCAL_ONLY"; work.mkdir()
    np.asarray(c.z, dtype="<f8").tofile(work / "X.bin")
    pd.DataFrame({"role": c.roles, "y": c.y, "weight": c.w}).to_csv(work / "roles.tsv", sep="\t", index=False)
    abm.json_write(work / "config.json", dict(n=len(c.ids), p=c.z.shape[1], seed=c.cfg.seed,
        maxit=max(100000, c.cfg.max_iter), feature_names=[c.names[j] for j in c.prep.keep]))
    subprocess.run([r_bin, str(HERE / "c1_unified_glmnet.R"), str(work), str(out)], check=True)
    raw = np.fromfile(out / "raw_probability_LOCAL_ONLY.bin", dtype="<f8")
    if len(raw) != len(c.ids):
        raise ValueError("R glmnet prediction row count mismatch")
    # These files duplicate the prepared inputs and frozen model; discard only
    # this private bridge's intermediates after successful R projection.
    shutil.rmtree(work)
    (out / "raw_probability_LOCAL_ONLY.bin").unlink()
    cal = c.parts["calibration_fit"]
    calibrator = abm.Calibrator().fit(raw[cal], c.y[cal], c.w[cal])
    p = calibrator.predict(raw)
    tuning = json.loads((out / "glmnet_fit.json").read_text())
    return finish_projection(c, p, shared_support(c), out, fitted=dict(prep=c.prep, km=c.km,
        support_bank=c.common_bank, geometry=c.common_geometry,
        calibrator=calibrator, r_model_file="glmnet_model_LOCAL_ONLY.rds"), metadata=dict(
        method="glmnet", engine="R glmnet", target="IPCW fixed-horizon incident logistic LASSO",
        scientific_status="matched-protocol implementation inspired by Pradeep; not exact publication reproduction",
        actual_device="cpu", tuning=tuning, fit_build_N=int(np.sum(c.w[c.parts['build']] > 0)),
        molecular_budget="all build-QC-passed features", case_loss_pruning=False))


def fit_lightgbm(c, out):
    from lightgbm import LGBMClassifier
    b, t, cal = (c.parts[r] for r in ("build", "tune_model", "calibration_fit"))
    b, t = b[c.w[b] > 0], t[c.w[t] > 0]
    best, trials = None, []
    for leaves in (7, 15, 31):
        fit = LGBMClassifier(n_estimators=c.cfg.tree_iterations, num_leaves=leaves,
            learning_rate=.03, min_child_samples=40, reg_lambda=10, n_jobs=c.cfg.cores,
            random_state=c.cfg.seed, verbosity=-1, deterministic=True, force_col_wise=True)
        fit.fit(c.z[b], c.y[b], sample_weight=c.w[b]/c.w[b].mean())
        pred = fit.predict_proba(c.z[t])[:,1]
        score = float(np.average(abm.logloss_rows(c.y[t], pred), weights=c.w[t]))
        trials.append(dict(num_leaves=leaves, tune_logloss=score))
        if best is None or score < best[0]: best = (score, fit)
    raw = best[1].predict_proba(c.z)[:,1]
    calibrator = abm.Calibrator().fit(raw[cal], c.y[cal], c.w[cal])
    return finish_projection(c, calibrator.predict(raw), shared_support(c), out,
        fitted=dict(prep=c.prep, km=c.km, support_bank=c.common_bank, geometry=c.common_geometry,
        model=best[1], calibrator=calibrator), metadata=dict(
        method="lightgbm", engine="LightGBM", actual_device="cpu", tuning=trials,
        fit_build_N=len(b), molecular_budget="all build-QC-passed features", case_loss_pruning=False))


def fit_twin(c, config, out, method):
    d = c.development
    # Even validation inside a backend must not depend on audit outcomes.
    time, event = np.asarray(d["time"], float).copy(), np.asarray(d["event"]).copy()
    audit = np.asarray(d["role"]).astype(str) == "calibration_audit"
    time[audit], event[audit] = config.horizon + 1., 0
    model = abm.TwinABM(config).fit(d["X"], time, event, d["eid"], d["group"], d["role"], c.names)
    # These private projections use named development roles solely for meta-
    # assessment. The frozen model's public test projection remains strict.
    p, diagnostics = model._calibrated(c.X, c.ids, c.groups)
    probability = p[config.primary]
    # Global fallback remains available for unsupported molecular neighbourhoods.
    unsupported = ~np.asarray(diagnostics["supported"], bool)
    probability = np.where(unsupported, p["elasticnet_full"], probability)
    metadata = dict(method=method, engine="local molecular twins", actual_device="cpu",
        primary_route=config.primary, training_audit=model.training_audit,
        fallback_route="calibrated elasticnet_full when twin support unavailable",
        case_loss_pruning=method == "abm_legacy", molecular_budget=config.panel_selection,
        model_feature_count=len(model.panel_indices), matching_feature_count=len(model.geometry_molecular_features))
    saved = finish_projection(c, probability, diagnostics, out,
        fitted=dict(twin=model, km=c.km), metadata=metadata, extra_probabilities=p)
    return saved


def fit_prevalent(c, history, out, history_config=None):
    """History transfer and its architecture-matched no-history control.

    Both nonlinear heads are trained on exactly the same incident build rows.
    Only the additional historical features differ. alpha in [0,.25] is chosen
    on tune_model before a single final mixture calibration on calibration_fit.
    The .25 bound applies to raw probability borrowing; recalibration can change
    the final probability difference, which is separately recorded.
    """
    from c1_prevalent import PrevalentConfig, PrevalentPrototypeTransformer
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    hcfg = PrevalentConfig(**{"seed": c.cfg.seed, **(history_config or {})})
    mol = np.asarray([j for j in range(len(c.names)) if j not in c.cfg.always_include_indices], int)
    if history is None:
        raise ValueError("Historical source was not prepared. Rebuild C1 inputs with include_prevalent; legacy incident NPZ cannot supply it.")
    hnames = list(map(str, history["feature_names"]))
    if hnames != [c.names[j] for j in mol]:
        raise ValueError("Historical molecular feature order differs from incident input")
    transformer = PrevalentPrototypeTransformer(hcfg).fit(c.X[:,mol], c.ids, c.groups, c.roles,
        history["prevalent_X"], history["prevalent_ids"], history["prevalent_groups"],
        history["prevalent_roles"], incident_metadata=history["incident_metadata"],
        prevalent_metadata=history["prevalent_metadata"], feature_names=hnames)
    H, hd = transformer.transform_with_diagnostics(c.X[:,mol], history["incident_metadata"], c.ids, c.groups)
    d = {**shared_support(c), **hd}
    # Same representation and model capacity in the no-history ablation.
    embedding = c.common_geometry.transform(c.z) / np.sqrt(c.common_geometry.n_components_)
    clinical = c.z[:, np.isin(c.prep.keep, c.cfg.always_include_indices)]
    Z = np.c_[embedding, clinical]
    b, tune, cal = (c.parts[r] for r in ("build", "tune_model", "calibration_fit"))
    b, tune = b[c.w[b] > 0], tune[c.w[tune] > 0]
    def fit_head(features):
        model = HistGradientBoostingClassifier(max_iter=c.cfg.tree_iterations, max_leaf_nodes=7,
            learning_rate=.04, min_samples_leaf=max(20, min(100, len(b)//30)), l2_regularization=10,
            early_stopping=False, random_state=c.cfg.seed)
        model.fit(features[b], c.y[b], sample_weight=c.w[b]/c.w[b].mean())
        return model, model.predict_proba(features)[:,1]
    baseline, p0 = fit_head(Z)
    source_available = bool(transformer.available_)
    if source_available:
        auxiliary, ph = fit_head(np.c_[Z,H])
    else:
        auxiliary, ph = None, p0.copy()
    trials = []
    for alpha in ((0., .1, .25) if source_available else (0.,)):
        probability = (1-alpha)*p0 + alpha*ph
        score = float(np.average(abm.logloss_rows(c.y[tune], probability[tune]), weights=c.w[tune]))
        trials.append(dict(alpha=alpha, tune_logloss=score))
    choice = min(trials, key=lambda row: (row["tune_logloss"], row["alpha"]))
    alpha = choice["alpha"]
    raw = (1-alpha)*p0 + alpha*ph
    raw_delta = float(np.max(np.abs(raw-p0)))
    if raw_delta > .25 + 1e-12:
        raise AssertionError("Historical borrowing cap violated")
    calibrator0 = abm.Calibrator().fit(p0[cal], c.y[cal], c.w[cal])
    control = calibrator0.predict(p0)
    # alpha=0 is the exact no-history control, including when the source is
    # unavailable. Reuse calibration, avoiding an unnecessary duplicate fit.
    calibrator = calibrator0 if alpha == 0 else abm.Calibrator().fit(raw[cal], c.y[cal], c.w[cal])
    p = control.copy() if alpha == 0 else calibrator.predict(raw)
    metadata = dict(method="abm_prevalent", engine="historical prototypes with incident-only nonlinear risk heads",
        actual_device="cpu", historical_source=history.get("audit", {}), source_transform=transformer.audit_,
        alpha=alpha, alpha_trials=trials, raw_probability_borrow_cap=.25,
        max_raw_probability_change=raw_delta, max_calibrated_probability_change=float(np.max(np.abs(p-control))),
        fit_build_N=len(b), fit_build_horizon_cases=int(np.sum(c.y[b])),
        fit_build_all_followup_events=int(np.sum(np.asarray(c.development["event"])[b] == 1)),
        historical_rows_in_incident_loss=0, case_loss_pruning=False,
        history_use_status="no_history_selected_alpha_zero" if alpha == 0 else "history_blend_selected_on_tune_model",
        comparator="abm_prevalent_control: same PCA, same HGB, same incident build/tune/calibration roles, no history features",
        support_definition="incident twins and prevalent-source support counted separately",
        reverse_causality="age/sex alignment and disclosed source weights do not identify a causal prediagnosis profile")
    fitted = dict(prep=c.prep, km=c.km, geometry=c.common_geometry, support_bank=c.common_bank,
        transformer=transformer, baseline=baseline, auxiliary=auxiliary, alpha=alpha,
        calibrator=calibrator, calibrator0=calibrator0, molecular_columns=mol)
    saved = finish_projection(c, p, d, out, fitted=fitted, metadata=metadata)
    ctrl = finish_projection(c, control, shared_support(c), out / "nohistory_control",
        fitted=dict(prep=c.prep, km=c.km, geometry=c.common_geometry, support_bank=c.common_bank,
                    model=baseline, calibrator=calibrator0),
        metadata=dict(method="abm_prevalent_control", engine="architecture-matched no-history control",
            actual_device="cpu", fit_build_N=len(b), case_loss_pruning=False,
            historical_rows_in_incident_loss=0))
    return saved, ctrl


def fit_deep(c, out, network_config, device="auto"):
    import c1_abm_gpu as deep
    # A literal outcome-free placeholder for test/audit prevents accidental
    # future additions to the deep backend from seeing their real labels.
    time = np.full(len(c.ids), c.cfg.horizon+1.)
    event = np.zeros(len(c.ids), int)
    allowed = np.concatenate([c.parts[r] for r in abm.ROLES[:4]])
    time[allowed], event[allowed] = c.development["time"][allowed], c.development["event"][allowed]
    data = dict(X=c.X, eid=c.ids, family_group=c.groups, role=c.roles, time=time, event=event,
        feature_names=c.names, clinical_feature_names=[c.names[j] for j in c.cfg.always_include_indices])
    fitted = deep.train(data, c.cfg, network_config, Path(out) / "transformer", device=device, report=False)
    diagnostics = fitted["diagnostics"]
    diagnostics["technical_qc"] = fitted["technical_qc"]
    return finish_projection(c, fitted["probability"], diagnostics, out,
        fitted=dict(backend_artifact="transformer/preprocessing_and_bank_LOCAL_ONLY.joblib", km=c.km),
        metadata=dict(method="abm_deep", engine="feature-token Transformer", actual_device=fitted["device"],
            case_loss_pruning=False, source_roles="same incident build/tune/calibration roles",
            model_support_geometry="learned molecular embedding", cuda_parity_claim=False))
