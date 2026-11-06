"""Predictor adapters with one shared side contract and disjoint calibration."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.special import expit, logit
from sklearn.linear_model import Ridge, LogisticRegression, LinearRegression
from sklearn.preprocessing import StandardScaler

from .losses import row_loss
from .contracts import ordered_ids_sha256


def require_labels(y, role, binary, *, min_n=30, min_class=10):
    yy = np.asarray(y, float)
    if len(yy) < min_n or not np.isfinite(yy).all():
        raise ValueError(f"{role}: insufficient observed labels; no holdout borrowing")
    if binary and (not np.isin(yy, [0, 1]).all() or
                   min(np.sum(yy == v) for v in (0, 1)) < min_class):
        raise ValueError(f"{role}: insufficient cases/controls; unavailable")


@dataclass
class PredictionCalibration:
    binary: bool
    model: object

    @classmethod
    def fit(cls, prediction, y, binary):
        require_labels(y, "prediction_calibration", binary)
        x = np.asarray(prediction, float)
        if not np.isfinite(x).all():
            raise ValueError("Incomplete predictions in prediction_calibration")
        x = logit(np.clip(x, 1e-7, 1 - 1e-7)) if binary else x
        m = (LogisticRegression(C=1e6, max_iter=3000) if binary else LinearRegression())
        m.fit(x[:, None], y)
        return cls(binary, m)

    def predict(self, prediction):
        p = np.asarray(prediction, float)
        out = np.full(len(p), np.nan)
        ok = np.isfinite(p)
        if ok.any():
            x = logit(np.clip(p[ok], 1e-7, 1 - 1e-7)) if self.binary else p[ok]
            out[ok] = (self.model.predict_proba(x[:, None])[:, 1] if self.binary
                       else self.model.predict(x[:, None]))
        return out


def fit_linear(x, y, roles, binary, config):
    build = roles == "build"
    tune = roles == "tune_model"
    grid = config.get("logistic_c_grid", [.01, .1, 1.]) if binary else config.get("ridge_grid", [1., 10., 100., 1000.])
    best, history = None, []
    for parameter in grid:
        m = (LogisticRegression(C=parameter, max_iter=3000, solver="lbfgs")
             if binary else Ridge(alpha=parameter))
        m.fit(x[build], y[build])
        p = m.predict_proba(x[tune])[:, 1] if binary else m.predict(x[tune])
        score = float(np.mean(row_loss(y[tune], p, "binary" if binary else "continuous")))
        history.append(dict(parameter=parameter, tune_loss=score))
        if best is None or score < best[0]:
            best = score, m
    return best[1], history


def fit_moe(x, y, roles, binary, csx_columns, config):
    from abm_routes import MixtureExperts
    best, history = None, []
    build, tune = roles == "build", roles == "tune_model"
    for k in config.get("moe_clusters", [4, 8, 16]):
        m = MixtureExperts().fit(
            x[build][:, csx_columns], x[build], y[build], np.ones(build.sum(), bool),
            k, int(config["seed"]), int(config.get("moe_min_n", 100)),
            int(config.get("moe_min_class", 10)))
        # Explicit trait semantics; MixtureExperts is reused only after these checks.
        if m.binary != binary:
            raise ValueError("MoE inferred trait type conflicts with declared outcome")
        p, _ = m.predict(x[tune][:, csx_columns], x[tune])
        score = float(np.mean(row_loss(y[tune], p, "binary" if binary else "continuous")))
        history.append(dict(clusters=k, tune_loss=score, fallback_clusters=m.fallback_clusters))
        if best is None or score < best[0]:
            best = score, m
    return best[1], history


def lineage(frame, predictor):
    record = {"predictor": predictor, "query_outcomes_used": False, "roles": {}}
    for role in ("build", "tune_model", "prediction_calibration"):
        d = frame.loc[frame.role.eq(role)]
        record["roles"][role] = dict(n=len(d), ids_sha256=ordered_ids_sha256(d.id),
                                    groups_sha256=ordered_ids_sha256(d.group_id))
    record["training_role"] = "build"
    record["selection_role"] = "tune_model"
    record["calibrator_role"] = "prediction_calibration"
    record["head_labels_out_of_fit"] = ["tune_gate", "head_calibration", "policy_calibration"]
    return record


def fit_predictors(frame, features, csx_columns, candidate, config, *,
                   blocks=None, block_rows=None, tokens=None, device="cpu"):
    """No test/audit Y is consumed. Caller passes NaN for all holdout outcomes."""
    binary = config["_binary"]
    y = frame.y.to_numpy(float)
    roles = frame.role.to_numpy(str)
    for role in ("build", "tune_model", "prediction_calibration"):
        require_labels(y[roles == role], role, binary)
    x = frame[features].to_numpy(float)
    if not np.isfinite(x).all():
        raise ValueError("Predictors require all declared side covariates; no silent imputation")
    scale = StandardScaler().fit(x[roles == "build"])
    z = scale.transform(x)
    base, base_history = fit_linear(z, y, roles, binary, config)
    raw_base = base.predict_proba(z)[:, 1] if binary else base.predict(z)
    base_cal = PredictionCalibration.fit(raw_base[roles == "prediction_calibration"],
                                         y[roles == "prediction_calibration"], binary)
    bundle = dict(features=features, csx_columns=csx_columns, scale=scale, base=base,
                  base_calibration=base_cal, candidate_kind=candidate, binary=binary,
                  config=config, histories={"linear": base_history})
    bundle["training_objectives"] = dict(linear="logloss" if binary else "mse",
        candidate=("brier" if candidate == "blocknet" else "logloss") if binary else "mse",
        selection="brier" if binary else "mse",
        candidate_target="build_standardized_y" if candidate == "blocknet" and not binary else "native_y",
        assessment="stored_separately_in_assessor")
    candidate_roles = roles
    raw_candidate = np.full(len(frame), np.nan)
    if candidate == "moe":
        m, hist = fit_moe(z, y, roles, binary, csx_columns, config)
        raw_candidate, _ = m.predict(z[:, csx_columns], z)
        bundle.update(candidate=m)
        bundle["histories"]["moe"] = hist
    elif candidate == "blocknet":
        import torch
        from blocknet import feature_scaler, train, predict_model, stage_block_matrix
        original_blocks = blocks
        blocks = stage_block_matrix(blocks, config)
        bundle["block_staging"] = dict(strategy="ram" if blocks is not original_blocks else "mmap",
                                      requested_gib=config.get("block_ram_cache_gib", 16),
                                      actual_bytes=blocks.nbytes if blocks is not original_blocks else 0)
        cfg = dict(config, amp=config.get("amp", "bf16") if device == "cuda" else "off",
                   binary_objective="brier", initialize_prevalence=True)
        torch.manual_seed(int(config["seed"]))
        valid = np.zeros(len(frame), bool)
        for start in range(0, len(frame), 512):
            b = np.asarray(blocks[block_rows[start:start + 512]])
            valid[start:start + len(b)] = np.isfinite(b).reshape(len(b), -1).all(1)
        net_roles = np.where(valid, roles, "genotype_qc_unavailable")
        candidate_roles = net_roles
        for role in ("build", "tune_model"):
            require_labels(y[net_roles == role], role, binary)
        mu, sd = feature_scaler(blocks, block_rows[net_roles == "build"])
        net, hist, info = train("transformer", blocks, block_rows, z.astype("float32"),
                                y, net_roles, tokens, cfg, torch.device(device), mu, sd)
        raw, _ = predict_model(net, "transformer", blocks, block_rows, z.astype("float32"),
                               np.flatnonzero(valid), mu, sd, cfg, torch.device(device),
                               info["ymean"], info["ysd"])
        raw_candidate[valid] = expit(raw) if binary else raw
        bundle.update(candidate_state={k: v.detach().cpu() for k, v in net.state_dict().items()},
                      block_mean=mu, block_sd=sd, tokens=tokens, block_info=info, config=cfg)
        bundle["histories"]["blocknet"] = hist.to_dict("records")
        del net
        if device == "cuda":
            torch.cuda.empty_cache()
    else:
        raise ValueError("candidate must be moe or blocknet")
    cal = roles == "prediction_calibration"
    candidate_cal = PredictionCalibration.fit(raw_candidate[cal], y[cal], binary)
    bundle["candidate_calibration"] = candidate_cal
    bundle["lineage"] = lineage(frame, "csx-" + candidate)
    bundle["lineage"].update(feature_columns=list(features), training_objectives=bundle["training_objectives"],
        supervision={role: dict(base_n=int(np.sum(roles == role)),
                                candidate_n=int(np.sum(candidate_roles == role)))
                     for role in ("build", "tune_model", "prediction_calibration")},
        calibration=dict(base_type=type(base_cal.model).__name__,
            candidate_type=type(candidate_cal.model).__name__,
            base_coef=base_cal.model.coef_.tolist(), base_intercept=np.asarray(base_cal.model.intercept_).tolist(),
            candidate_coef=candidate_cal.model.coef_.tolist(),
            candidate_intercept=np.asarray(candidate_cal.model.intercept_).tolist()))
    predictions = frame[["id", "group_id", "role"]].copy()
    predictions["p_base"] = base_cal.predict(raw_base)
    predictions["p_candidate"] = candidate_cal.predict(raw_candidate)
    return bundle, predictions


def predict(bundle, frame, *, blocks=None, block_rows=None, device="cpu"):
    x = frame[bundle["features"]].to_numpy(float)
    if not np.isfinite(x).all():
        raise ValueError("Missing declared predictor covariates; baseline cannot be fabricated")
    z = bundle["scale"].transform(x)
    base = bundle["base"]
    p0 = base.predict_proba(z)[:, 1] if bundle["binary"] else base.predict(z)
    if bundle["candidate_kind"] == "moe":
        p1, _ = bundle["candidate"].predict(z[:, bundle["csx_columns"]], z)
    else:
        import torch
        from blocknet import BlockTransformer, predict_model
        cfg = dict(bundle["config"])
        if device == "cpu":
            cfg["amp"] = "off"
        net = BlockTransformer(bundle["tokens"], blocks.shape[2], z.shape[1], cfg).to(device)
        net.load_state_dict(bundle["candidate_state"])
        info = bundle["block_info"]
        p1, _ = predict_model(net, "transformer", blocks, block_rows, z.astype("float32"),
                              np.arange(len(frame)), bundle["block_mean"], bundle["block_sd"],
                              cfg, torch.device(device), info["ymean"], info["ysd"])
        if bundle["binary"]:
            p1 = expit(p1)
    return bundle["base_calibration"].predict(p0), bundle["candidate_calibration"].predict(p1)
