#!/usr/bin/env python3
"""The C1 molecular Transformer, with the same architecture on CPU or CUDA."""
from __future__ import annotations
import importlib.util
import json
import os
from pathlib import Path
import random
import time
import numpy as np
import pandas as pd
import torch
from torch import nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
import joblib
import abm
VERSION = "LE8-C1-deep-current"
HERE = Path(__file__).resolve().parent


def host_module():
    spec = importlib.util.spec_from_file_location("le8_deep_host", HERE / "0.host.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def setup_device(seed, cores=2, device="cpu"):
    """One scientific network on CPU/CUDA; explicit CUDA never silently falls back."""
    from c1_protocol import resolve_device
    policy = resolve_device("abm_deep", device)
    actual = policy["actual"]
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(max(1, int(cores)))
    if actual == "cuda":
        torch.cuda.manual_seed_all(seed)
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.backends.cudnn.benchmark = False
        torch.backends.cuda.enable_flash_sdp(False)
        torch.backends.cuda.enable_mem_efficient_sdp(False)
        torch.backends.cuda.enable_math_sdp(True)
    torch.use_deterministic_algorithms(True)
    label = torch.cuda.get_device_name(0) if actual == "cuda" else "CPU"
    print(f"[LE8 abm_deep] device={actual} ({label}); model=float32; seed={seed}", flush=True)
    return actual


def synchronize(device):
    if str(device).startswith("cuda"):
        torch.cuda.synchronize()


def validate_config(c):
    if c["architecture"] != "feature_token_latent_transformer":
        raise ValueError("Unknown deep architecture")
    for key in ("d_model", "heads", "latent_tokens", "layers", "batch_size", "epochs", "patience", "neighbors"):
        if type(c[key]) is not int or c[key] < 1:
            raise ValueError(f"Invalid {key}")
    if c["d_model"] % c["heads"] or not 0 <= c["dropout"] < 1:
        raise ValueError("Invalid attention dimensions/dropout")
    if c["coverages"] != [0.2,0.4,0.5,0.6,0.8,1.0] or any(c[k] <= 0 for k in ("learning_rate", "gradient_clip")):
        raise ValueError("Invalid optimization/coverage settings")
    if c["weight_decay"] < 0 or c["min_delta"] < 0:
        raise ValueError("Invalid regularization/early stopping")

    if c["transformer_matching_space"] != "molecular":
        raise ValueError("Transformer matching must exclude clinical predictors")


class MolecularTransformer(nn.Module):
    """Feature-specific numeric tokens + missingness, latent cross-attention,
    then self-attention among latent tokens. Latent cross-attention scales linearly
    in selected molecular tokens; this controlled experiment reuses the frozen
    CPU development-selected panel and can later extend to larger panels.
    Attention is within each participant; other participants are never tokens.
    """
    def __init__(self, n_features, config, clinical_positions=()):
        super().__init__()
        clinical_positions = list(map(int, clinical_positions))
        molecular_positions = [j for j in range(n_features) if j not in clinical_positions]
        if not molecular_positions:
            raise ValueError("No molecular tokens")
        self.register_buffer("molecular_positions", torch.tensor(molecular_positions, dtype=torch.long))
        self.register_buffer("clinical_positions", torch.tensor(clinical_positions, dtype=torch.long))
        n_features = len(molecular_positions)
        d = config["d_model"]
        self.value_weight = nn.Parameter(torch.empty(n_features, d))
        self.feature_bias = nn.Parameter(torch.empty(n_features, d))
        self.missing_embedding = nn.Parameter(torch.empty(n_features, d))
        self.latents = nn.Parameter(torch.empty(config["latent_tokens"], d))
        for p in (self.value_weight, self.feature_bias, self.missing_embedding, self.latents):
            nn.init.normal_(p, std=.02)
        self.input_norm = nn.LayerNorm(d)
        self.cross = nn.MultiheadAttention(d, config["heads"], dropout=config["dropout"], batch_first=True)
        # Construct each block separately; do not clone identical initialized layers.
        self.blocks = nn.ModuleList([
            nn.TransformerEncoderLayer(d, config["heads"], dim_feedforward=d*4,
                                       dropout=config["dropout"], activation="gelu",
                                       batch_first=True, norm_first=True)
            for _ in range(config["layers"])
        ])
        self.output_norm = nn.LayerNorm(d)
        self.head = nn.Linear(d + len(clinical_positions), 1)

    def forward(self, x, missing):
        clinical = x[:, self.clinical_positions]
        x, missing = x[:, self.molecular_positions], missing[:, self.molecular_positions]
        tokens = (x[..., None] * self.value_weight + self.feature_bias
                  + missing[..., None] * self.missing_embedding)
        query = self.latents[None].expand(len(x), -1, -1)
        tokens = self.input_norm(tokens)
        z = query + self.cross(query, tokens, tokens, need_weights=False)[0]
        for block in self.blocks:
            z = block(z)
        embedding = self.output_norm(z.mean(dim=1))
        return self.head(torch.cat([embedding, clinical], dim=1)).squeeze(-1), embedding


def transform(prep, X):
    missing = ~np.isfinite(np.asarray(X)[:, prep.keep])
    return np.asarray(prep.transform(X), dtype=np.float32), missing.astype(np.float32)


def predict(model, x, missing, indices, batch_size):
    model.eval()
    device = next(model.parameters()).device
    scores, embeddings = [], []
    with torch.inference_mode():
        for start in range(0, len(indices), batch_size):
            ix = indices[start:start+batch_size]
            logits, z = model(torch.as_tensor(x[ix], device=device),
                              torch.as_tensor(missing[ix], device=device))
            scores.append(logits.cpu().numpy())
            embeddings.append(z.cpu().numpy())
    return np.concatenate(scores), np.concatenate(embeddings)


def metrics(y, w, p):
    ok = (w > 0) & np.isfinite(p)
    yy, ww, pp = y[ok], w[ok], np.clip(p[ok], 1e-7, 1-1e-7)
    if len(np.unique(yy)) != 2:
        return dict(n=int(ok.sum()), status="insufficient_classes")
    return dict(n=int(ok.sum()), events=int(yy.sum()), status="ok",
                auc=float(roc_auc_score(yy, pp, sample_weight=ww)),
                average_precision=float(average_precision_score(yy, pp, sample_weight=ww)),
                brier=float(np.average((yy-pp)**2, weights=ww)),
                logloss=float(np.average(abm.logloss_rows(yy, pp), weights=ww)))


def shared_knn():
    spec = importlib.util.spec_from_file_location("le8_deep_kernels", HERE.parents[1] / "0f/abm.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.BlockedKNN


def distance_selection(distances, qc, gate, coverages):
    valid = qc[gate] & np.isfinite(distances[gate])
    if not valid.any():
        raise ValueError("No eligible tune_gate distances")
    thresholds = {float(q): float(np.quantile(distances[gate][valid], q)) for q in coverages}
    masks = {q: qc & np.isfinite(distances) & (distances <= t) for q, t in thresholds.items()}
    return thresholds, masks


def train(data, cfg, network, out, *, device="cpu"):
    device = setup_device(cfg.seed, cfg.cores, device)
    validate_config(network)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    X = np.asarray(data["X"])
    ids, groups, roles = (np.asarray(data[k]).astype(str) for k in ("eid", "family_group", "role"))
    host_module().check_group_roles(groups, roles)
    if set(roles) != set(abm.ROLES) | {"test"} or len(set(ids)) != len(ids):
        raise ValueError("Invalid shared cohort roles/IDs")
    names = list(map(str, data["feature_names"]))
    forbidden = {"time", "event", "outcome", "target", "prevalent_status", "role", "split", "diagnosis_date", "t2e", "Yt2e"}
    if (X.ndim != 2 or len(X) != len(ids) or len(names) != X.shape[1]
            or len(set(names)) != len(names) or forbidden & set(names)
            or any(n.lower().endswith((".t2e", ".yt2e", ".bi2e")) for n in names)):
        raise ValueError("Invalid feature schema or outcome in predictors")
    if (not np.isfinite(data["time"]).all() or np.any(np.asarray(data["time"]) <= 0)
            or not np.isin(data["event"], [0,1]).all()):
        raise ValueError("Invalid survival outcomes")
    parts = {r: np.flatnonzero(roles == r) for r in (*abm.ROLES, "test")}
    b = parts["build"]
    from dataclasses import replace
    cfg = replace(cfg, always_include_indices=tuple(names.index(n) for n in data["clinical_feature_names"]))
    prep = abm.NumericPrep().fit(X[b], cfg)
    qc = prep.technical_qc(X)
    km = abm.CensoringKM().fit(data["time"][b][qc[b]], data["event"][b][qc[b]], cfg.horizon, cfg.min_censor_survival)
    y, w = km.labels_weights(data["time"], data["event"])
    w = w * qc
    x, missing = transform(prep, X)
    train_idx = b[w[b] > 0]
    tune = parts["tune_model"]
    abm.check_classes(y[train_idx], w[train_idx], cfg.min_events)
    abm.check_classes(y[tune], w[tune], cfg.min_events)
    clinical_pos = np.flatnonzero(np.isin(prep.keep, prep.clinical_columns))
    model = MolecularTransformer(x.shape[1], network, clinical_pos).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=network["learning_rate"], weight_decay=network["weight_decay"])
    rng = np.random.default_rng(cfg.seed)
    best, best_epoch, best_state, history = np.inf, 0, None, []
    normalized_w = w / w[train_idx].mean()
    synchronize(device)
    training_started = time.perf_counter()
    for epoch in range(1, network["epochs"]+1):
        started = time.monotonic()
        model.train()
        indices = rng.permutation(train_idx)
        for start in range(0, len(indices), network["batch_size"]):
            ix = indices[start:start+network["batch_size"]]
            optimizer.zero_grad(set_to_none=True)
            logits, _ = model(torch.as_tensor(x[ix], device=device), torch.as_tensor(missing[ix], device=device))
            loss = (nn.functional.binary_cross_entropy_with_logits(
                logits, torch.as_tensor(y[ix], dtype=torch.float32, device=device), reduction="none")
                * torch.as_tensor(normalized_w[ix], dtype=torch.float32, device=device)).mean()
            if not torch.isfinite(loss):
                raise RuntimeError("Non-finite abm_deep training loss")
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), network["gradient_clip"], error_if_nonfinite=True)
            optimizer.step()
        score, _ = predict(model, x, missing, tune, network["batch_size"])
        probability = 1 / (1 + np.exp(-np.clip(score, -50, 50)))
        valid = metrics(y[tune], w[tune], probability)["logloss"]
        improved = valid < best - network["min_delta"]
        if improved:
            best, best_epoch = valid, epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            torch.save(dict(version=VERSION, network=network, state_dict=best_state,
                            features=np.asarray(data["feature_names"])[prep.keep].tolist(),
                            epoch=epoch, validation_loss=valid), out / "best_checkpoint_LOCAL_ONLY.pt")
        history.append(dict(epoch=epoch, tune_logloss=valid, seconds=time.monotonic()-started, best_epoch=best_epoch))
        pd.DataFrame(history).to_csv(out / "training_history.tsv", sep="\t", index=False)
        print(f"[LE8 abm_deep] epoch={epoch}/{network['epochs']} tune_logloss={valid:.6f} "
              f"best={best_epoch} seconds={history[-1]['seconds']:.1f}", flush=True)
        if epoch-best_epoch >= network["patience"]:
            break
    synchronize(device)
    training_seconds = time.perf_counter() - training_started
    model.load_state_dict(best_state)
    score, z = predict(model, x, missing, np.arange(len(x)), network["batch_size"])
    cal = parts["calibration_fit"]
    cal = cal[w[cal] > 0]
    abm.check_classes(y[cal], w[cal], cfg.min_events)
    calibrator = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000).fit(
        score[cal, None], y[cal], sample_weight=w[cal])
    if calibrator.n_iter_.max() >= calibrator.max_iter:
        raise RuntimeError("Deep model probability calibration did not converge")
    p = calibrator.predict_proba(score[:, None])[:, 1]
    # Exact retrieval on the selected device with molecular-only embeddings and family/self exclusion.
    bank_idx = train_idx
    bank = shared_knn()(z[bank_idx], ids[bank_idx], groups[bank_idx],
                       [abm.tie_hash(v, cfg.seed) for v in ids[bank_idx]], device=device, dtype="float64")
    # Set support geometry using only the independent tune_gate role. No
    # holdout-derived threshold, quantile or label is used by this procedure.
    neighbors = min(network["neighbors"], len(bank_idx))
    radius_reference_k = min(int(getattr(cfg, "radius_k_ref", 0)) or neighbors, len(bank_idx))
    distances = np.full(len(ids), np.nan)
    gate = parts["tune_gate"]
    gate_qc = gate[qc[gate]]
    _, gate_ds = bank.query(z[gate_qc], ids[gate_qc], groups[gate_qc], max(neighbors, radius_reference_k))
    distances[gate_qc] = gate_ds[:, neighbors - 1]
    valid_gate = qc[gate] & np.isfinite(distances[gate])
    if not valid_gate.any():
        raise ValueError("No independent donors for coverage calibration")
    reference_distance = gate_ds[:, radius_reference_k - 1]
    if not np.isfinite(reference_distance).any():
        raise ValueError("No independent references for support-radius calibration")
    support_radius = float(np.quantile(reference_distance[np.isfinite(reference_distance)], cfg.radius_quantile))
    names = ("eligible_match_count", "matched_count", "eligible_reference_count", "twin_fraction",
             "family_count", "family_ess", "donor_ess", "matching_radius", "radius_ratio", "nearest_radius_ratio")
    diagnostics = {name: np.full(len(ids), np.nan) for name in names}
    diagnostics["technical_qc"] = qc.copy()
    diagnostics["supported"] = np.zeros(len(ids), bool)
    bank_groups = groups[bank_idx]
    for role in ("tune_model", "tune_gate", "calibration_fit", "calibration_audit", "test"):
        ix = parts[role]
        good = ix[qc[ix]]
        if not len(good):
            continue
        near, ds, count = bank.query(z[good], ids[good], groups[good], neighbors,
                                     radius=support_radius, return_counts=True)
        distances[good] = ds[:, -1]
        diagnostics["eligible_match_count"][good] = count
        # check_group_roles established disjoint build/non-build families.
        diagnostics["eligible_reference_count"][good] = len(bank_idx)
        diagnostics["twin_fraction"][good] = count / len(bank_idx)
        diagnostics["matching_radius"][good] = support_radius
        denominator = max(support_radius, np.finfo(float).eps)
        diagnostics["radius_ratio"][good] = ds[:, -1] / denominator
        diagnostics["nearest_radius_ratio"][good] = ds[:, 0] / denominator
        for local, row in enumerate(good):
            valid = (near[local] >= 0) & np.isfinite(ds[local]) & (ds[local] <= support_radius)
            chosen = near[local][valid]
            diagnostics["matched_count"][row] = len(chosen)
            if not len(chosen):
                diagnostics["donor_ess"][row] = diagnostics["family_ess"][row] = 0.
                diagnostics["family_count"][row] = 0.
                continue
            # Support weights are outcome-blind distance kernels. The network
            # does not borrow neighbor labels, so do not label these risk weights.
            mass = np.exp(-.5 * (ds[local][valid] / denominator)**2)
            mass /= mass.sum()
            diagnostics["donor_ess"][row] = 1. / np.sum(mass**2)
            levels, inverse = np.unique(bank_groups[chosen], return_inverse=True)
            family_mass = np.bincount(inverse, weights=mass)
            diagnostics["family_ess"][row] = 1. / np.sum(family_mass**2)
            diagnostics["family_count"][row] = len(levels)
    diagnostics["supported"] = (qc & (diagnostics["eligible_match_count"] >= cfg.min_twins) &
                                  (diagnostics["family_ess"] >= cfg.min_family_ess))
    eligible_match_count = diagnostics["eligible_match_count"]
    eligible_reference_count = diagnostics["eligible_reference_count"]
    thresholds, masks = distance_selection(distances, qc, gate, network["coverages"])
    selected = masks[0.5]
    risk_edges = np.unique(np.quantile(p[gate][valid_gate], np.linspace(0, 1, 11)[1:-1]))
    # Freeze all learned parameters before any outcome-based report is made.
    joblib.dump(dict(prep=prep, km=km, calibrator=calibrator, thresholds=thresholds, risk_edges=risk_edges,
                     support_radius=support_radius, clinical_positions=clinical_pos,
                     bank_embeddings=z[bank_idx], bank_ids=ids[bank_idx], bank_groups=groups[bank_idx],
                     feature_names=data["feature_names"], network=network, seed=cfg.seed),
                out / "preprocessing_and_bank_LOCAL_ONLY.joblib", compress=3)
    abm.json_write(out / "MODEL_FROZEN.json", dict(version=VERSION, best_epoch=best_epoch,
        training_ids_hash=abm.fingerprint(ids[train_idx]), tuning_ids_hash=abm.fingerprint(ids[tune]),
        roles={k: len(v) for k,v in parts.items()}, preprocessing="build only", censoring="build only",
        selection="tune_model loss only", calibration="calibration_fit only", horizon=cfg.horizon,
        confidence="kth independent-family donor distance; tune_gate threshold; not clinical certainty",
        matching_radius=support_radius, radius_reference_k=radius_reference_k,
        raw_count_definition="All eligible incident build references within frozen radius; no K cap",
        support_ess_definition="Gaussian distance-kernel mass among at most K within-radius neighbors; no outcome weights",
        inference_geometry="molecular-only learned embedding",
        model="abm_deep", algorithm_equivalent_to_twin=False, same_architecture_across_devices=True,
        cross_device_numerical_parity_verified=False, exploratory=True, confirmatory=False,
        clinical_predictors_excluded_from_matching=True,
        risk_input="molecular embedding + baseline clinical predictors",
        selector="embedding_neighbor_distance", selector_status="exploratory", outcome_blind=True, config=network, seed=cfg.seed,
        torch_version=torch.__version__, device=device,
        device_name=torch.cuda.get_device_name(0) if device == "cuda" else "CPU", clinical_release=False))
    test = parts["test"]
    private = dict(eid=ids[test], family=groups[test], probability=p[test], selected=selected[test],
                   technical_qc=qc[test], donor_distance=distances[test],
                   eligible_match_count=eligible_match_count[test],
                   eligible_reference_count=eligible_reference_count[test],
                   matching_radius=np.repeat(support_radius, len(test)),
                   **{f"selected_{q:g}": mask[test] for q, mask in masks.items()})
    pd.DataFrame(private).to_csv(out / "test_predictions_LOCAL_ONLY.tsv.gz", sep="\t", index=False)
    return dict(probability=p, diagnostics=diagnostics, technical_qc=qc, selected=selected, masks=masks,
                donor_distance=distances, eligible_match_count=eligible_match_count,
                eligible_reference_count=eligible_reference_count, matching_radius=support_radius,
                thresholds=thresholds, parts=parts, freeze_path=str(out / "MODEL_FROZEN.json"),
                device=device, model=model, prep=prep, km=km, calibrator=calibrator,
                training_case_N=int(np.sum(np.asarray(data["event"])[b][qc[b]] == 1)),
                incident_cases_removed_by_loss=0)
