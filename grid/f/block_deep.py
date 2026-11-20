#!/usr/bin/env python3
"""Optional consensus-block Transformer with a conservative CSx residual head.

This follows LE8 c1's participant-wise feature tokens and latent attention, but
uses regional PRS rather than proteins and no incident-survival/IPCW machinery.
There is no attention between people. The learned embedding provides optional
support diagnostics only: it never borrows labels from embedding neighbours.
PyTorch is imported only by explicitly requested training/prediction.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import time

import joblib
import numpy as np
import pandas as pd
from scipy.special import expit, logit
from scipy.spatial import cKDTree

from core2 import (aligned, atomic_json, atomic_table, check_roles, identified,
                   read_table)
from abm2 import check_oof, development_outcomes, eligible_nearest, batches

VERSION = "grid-block-deep-20261116"
PRIMARY = "csx.population"
METHOD = "abm.deep"
COVERAGES = (.2, .4, .5, .6, .8, 1.)
DEFAULTS = dict(device="cpu", d_model=32, heads=4, latent_tokens=8, layers=2,
                dropout=.1, batch_size=128, epochs=40, patience=6,
                learning_rate=.001, weight_decay=.01, gradient_clip=1.,
                min_delta=1e-5, gate_se=1., alpha_grid=[0., .25, .5, 1.],
                support=True, support_dimensions=8, radius_k=64,
                radius_quantile=.9, support_neighbors=64,
                min_support=20, min_group_ess=15.)


def torch_api():
    try:
        import torch
    except ImportError as exc:
        raise RuntimeError("The optional abm_deep stage requires PyTorch; install the CPU or CUDA wheel matching this machine") from exc
    return torch


def settings(config=None):
    result = {**deepcopy(DEFAULTS), **(config or {})}
    for key in ("d_model", "heads", "latent_tokens", "layers", "batch_size", "epochs",
                "patience", "support_dimensions", "radius_k", "support_neighbors", "min_support"):
        if isinstance(result[key], bool) or int(result[key]) != result[key] or result[key] < 1:
            raise ValueError(f"abm_deep.{key} must be a positive integer")
        result[key] = int(result[key])
    if result["d_model"] % result["heads"]:
        raise ValueError("abm_deep.d_model must be divisible by heads")
    if not 0 <= result["dropout"] < 1 or not 0 < result["radius_quantile"] < 1:
        raise ValueError("Invalid dropout or support-radius quantile")
    for key in ("learning_rate", "gradient_clip", "min_group_ess"):
        if not np.isfinite(result[key]) or result[key] <= 0:
            raise ValueError(f"Invalid positive abm_deep.{key}")
    for key in ("weight_decay", "min_delta", "gate_se"):
        if not np.isfinite(result[key]) or result[key] < 0:
            raise ValueError(f"Invalid nonnegative abm_deep.{key}")
    alpha = sorted(set(map(float, result["alpha_grid"])))
    if not alpha or alpha[0] != 0 or not np.isfinite(alpha).all() or alpha[-1] > 1:
        raise ValueError("abm_deep.alpha_grid must include 0 and lie in [0,1]")
    result["alpha_grid"] = alpha
    if result["device"] not in {"cpu", "gpu", "cuda", "auto"}:
        raise ValueError("abm_deep.device must be cpu, gpu, cuda, or auto")
    return result


def setup_device(requested, seed, threads):
    torch = torch_api()
    requested = "cuda" if requested == "gpu" else requested
    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was explicitly requested for abm_deep but is unavailable; no CPU fallback was performed")
    # Match the shared LE8 policy: auto remains CPU unless GPU is explicit.
    device = "cpu" if requested == "auto" else requested
    if device not in {"cpu", "cuda"}:
        raise ValueError("Unknown deep device")
    if threads < 1:
        raise ValueError("Deep worker threads must be positive")
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.manual_seed(int(seed))
    torch.set_num_threads(int(threads))
    if device == "cuda":
        torch.cuda.manual_seed_all(int(seed))
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.backends.cudnn.benchmark = False
        torch.backends.cuda.enable_flash_sdp(False)
        torch.backends.cuda.enable_mem_efficient_sdp(False)
        torch.backends.cuda.enable_math_sdp(True)
    torch.use_deterministic_algorithms(True)
    return torch, device


def feature_schema(manifest, width):
    manifest = pd.DataFrame(manifest).copy()
    if len(manifest) != width or width < 1:
        raise ValueError("Regional feature manifest does not match the matrix width")
    columns = [c for c in ("block_id", "name", "map", "population", "chrom38", "start38", "end38") if c in manifest]
    if not columns:
        raise ValueError("Regional feature manifest needs block IDs, names, or GRCh38 coordinates")
    schema = manifest[columns].fillna("").astype(str).to_dict("records")
    identities = [json.dumps(row, sort_keys=True) for row in schema]
    if len(set(identities)) != width:
        raise ValueError("Duplicate regional feature identities")
    return schema


def fit_transformer_input(features, train, batch_size):
    """TRAIN-only mean and one RMS scale; preserve relative block amplitudes."""
    width = features.shape[1]
    mean, m2, seen = np.zeros(width), np.zeros(width), 0
    for rows in batches(np.asarray(train), batch_size):
        values = np.asarray(features[rows], dtype=np.float64)
        if not np.isfinite(values).all():
            raise ValueError("Regional PRS must be finite before deep fitting")
        n = len(values)
        local_mean = values.mean(0)
        delta = local_mean - mean
        m2 += ((values-local_mean)**2).sum(0) + delta**2 * seen*n/(seen+n)
        mean += delta*n/(seen+n)
        seen += n
    if seen < 3:
        raise ValueError("Deep fitting needs at least three observed TRAIN participants")
    variance = m2/seen
    keep = np.flatnonzero(variance > 1e-14)
    if not len(keep):
        raise ValueError("No variable consensus-block PRS in TRAIN")
    return dict(mean=mean[keep], scale=float(np.sqrt(variance[keep].mean())),
                keep=keep, input_width=width, training_rows=seen,
                normalization="TRAIN mean per block; one common RMS scale preserves relative block variances")


def transform(features, prep, rows):
    values = np.asarray(features[rows], dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != prep["input_width"] or not np.isfinite(values).all():
        raise ValueError("Query regional features have the wrong width or nonfinite values")
    return np.asarray((values[:, prep["keep"]]-prep["mean"])/prep["scale"], dtype=np.float32)


def build_model(n_features, config):
    torch = torch_api()
    nn = torch.nn

    class BlockTransformer(nn.Module):
        def __init__(self):
            super().__init__()
            d = config["d_model"]
            self.value_weight = nn.Parameter(torch.empty(n_features, d))
            self.block_bias = nn.Parameter(torch.empty(n_features, d))
            self.latents = nn.Parameter(torch.empty(config["latent_tokens"], d))
            for value in (self.value_weight, self.block_bias, self.latents):
                nn.init.normal_(value, std=.02)
            self.input_norm = nn.LayerNorm(d)
            self.cross = nn.MultiheadAttention(d, config["heads"], dropout=config["dropout"], batch_first=True)
            self.blocks = nn.ModuleList([
                nn.TransformerEncoderLayer(d, config["heads"], dim_feedforward=4*d,
                    dropout=config["dropout"], activation="gelu", batch_first=True, norm_first=True)
                for _ in range(config["layers"])
            ])
            self.output_norm = nn.LayerNorm(d)
            self.head = nn.Linear(d, 1)
            # Start at the externally supplied strong baseline, not at mean Y.
            nn.init.zeros_(self.head.weight)
            nn.init.zeros_(self.head.bias)

        def forward(self, x):
            tokens = self.input_norm(x[..., None]*self.value_weight+self.block_bias)
            query = self.latents[None].expand(len(x), -1, -1)
            z = query + self.cross(query, tokens, tokens, need_weights=False)[0]
            for layer in self.blocks:
                z = layer(z)
            embedding = self.output_norm(z.mean(1))
            return self.head(embedding).squeeze(-1), embedding

    return BlockTransformer()


def predict_network(model, features, prep, rows, batch_size):
    torch = torch_api()
    device = next(model.parameters()).device
    scores, embeddings = [], []
    model.eval()
    with torch.inference_mode():
        for batch in batches(np.asarray(rows), batch_size):
            x = torch.as_tensor(transform(features, prep, batch), device=device)
            score, embedding = model(x)
            scores.append(score.cpu().numpy())
            embeddings.append(embedding.cpu().numpy())
    if not scores:
        return np.empty(0), np.empty((0, model.head.in_features), dtype=np.float32)
    result, latent = np.concatenate(scores).astype(float), np.concatenate(embeddings)
    if not np.isfinite(result).all() or not np.isfinite(latent).all():
        raise RuntimeError("Nonfinite deep inference")
    return result, latent


def combine(base, delta, task, scale, alpha):
    base = np.asarray(base, float)
    if alpha == 0:
        return base.copy()
    return expit(logit(np.clip(base, 1e-8, 1-1e-8))+alpha*delta) if task == "binary" else base+alpha*scale*delta


def primary_rows(y, prediction):
    return np.square(np.asarray(y, float)-np.asarray(prediction, float))


def cluster_mean_se(differences, groups):
    differences = np.asarray(differences, float)
    _, codes = np.unique(np.asarray(groups, str), return_inverse=True)
    count = len(np.unique(codes))
    if count < 2:
        return np.inf
    score = np.bincount(codes, weights=differences-differences.mean())
    return float(np.sqrt(count/(count-1)*np.square(score).sum())/len(differences))


def choose_alpha(y, base, delta, groups, task, scale, options):
    baseline = primary_rows(y, base)
    records, best_alpha, best_loss = [], 0., float(baseline.mean())
    for alpha in options["alpha_grid"]:
        candidate = combine(base, delta, task, scale, alpha)
        loss = primary_rows(y, candidate)
        difference = baseline-loss
        se = cluster_mean_se(difference, groups)
        accepted = alpha == 0 or difference.mean() > options["gate_se"]*se
        records.append(dict(alpha=alpha, tune_loss=float(loss.mean()),
                            mean_improvement=float(difference.mean()), family_cluster_se=se,
                            passes_gate=bool(accepted)))
        if accepted and loss.mean() < best_loss:
            best_alpha, best_loss = alpha, float(loss.mean())
    return best_alpha, records


def embedding_projection(embeddings, train, dimensions):
    values = np.asarray(embeddings[train], float)
    center = values.mean(0)
    covariance = (values-center).T@(values-center)/len(values)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    total = float(eigenvalues[eigenvalues > 0].sum())
    rank = min(dimensions, int((eigenvalues > 1e-12).sum()))
    if rank < 1 or total <= 1e-12:
        return None
    chosen = order[:rank]
    return dict(center=center, components=eigenvectors[:, chosen],
                dimensions=rank, captured_train_embedding_variance=float(eigenvalues[chosen].sum()/total))


def query_support(bank, z, groups, options, workers=1):
    """Exact projected-embedding counts without materializing all radius members.

Counts include all eligible donors. Family ESS intentionally uses only the at
most K matched neighbours; the two quantities are separately named/reported.
"""
    donor_z, donor_groups = bank["z"], bank["groups"]
    tree = cKDTree(donor_z)
    group_rows = pd.Series(np.arange(len(donor_groups))).groupby(donor_groups).apply(lambda s:s.to_numpy()).to_dict()
    k = min(options["support_neighbors"], bank["effective_k"])
    count, ess, matched, nearest, distance = (np.zeros(len(z), dtype=float) for _ in range(5))
    for rows in batches(np.arange(len(z)), options["batch_size"]):
        ii, dd = eligible_nearest(tree, donor_groups, z[rows], groups[rows], k, workers)
        counts = tree.query_ball_point(z[rows], bank["radius"], workers=workers, return_length=True)
        for local, row in enumerate(rows):
            own = group_rows.get(groups[row], np.empty(0, int))
            own_n = int((np.linalg.norm(donor_z[own]-z[row], axis=1) <= bank["radius"]).sum()) if len(own) else 0
            count[row] = max(0, int(counts[local])-own_n)
            valid = (ii[local] >= 0) & np.isfinite(dd[local]) & (dd[local] <= bank["radius"])
            chosen = ii[local][valid]
            matched[row] = len(chosen)
            if len(chosen):
                _, sizes = np.unique(donor_groups[chosen], return_counts=True)
                ess[row] = len(chosen)**2/np.square(sizes.astype(float)).sum()
            nearest[row], distance[row] = dd[local, 0], dd[local, -1]
    supported = (count >= options["min_support"]) & (ess >= options["min_group_ess"]) & np.isfinite(distance)
    return dict(n_support=count.astype(int), matched_count=matched.astype(int),
                matched_family_ess=ess, nearest_distance=nearest, kth_distance=distance,
                supported=supported)


def fit_support(embeddings, cohort, train, options, workers):
    projection = embedding_projection(embeddings, train, options["support_dimensions"])
    if projection is None:
        return None, {"status":"degenerate_embedding", "uses_query_outcomes":False}
    z = (np.asarray(embeddings, float)-projection["center"])@projection["components"]
    groups = cohort.group_id.to_numpy(str)
    donor_groups = groups[train]
    maximum_family = int(pd.Series(donor_groups).value_counts().max())
    effective_k = min(options["radius_k"], len(train)-maximum_family)
    if effective_k < 1:
        return None, {"status":"insufficient_independent_donors", "uses_query_outcomes":False}
    tree = cKDTree(z[train])
    calibration = np.flatnonzero(cohort.role.eq("calibration"))
    distances = []
    for rows in batches(calibration, options["batch_size"]):
        _, d = eligible_nearest(tree, donor_groups, z[rows], groups[rows], effective_k, workers)
        distances.extend(d[:, -1])
    if not np.isfinite(distances).all():
        raise ValueError("Calibration participants lack enough independent TRAIN donors")
    radius = max(float(np.quantile(distances, options["radius_quantile"])), 1e-10)
    bank = dict(z=z[train], groups=donor_groups, ids=cohort.id.to_numpy(str)[train],
                projection=projection, radius=radius, effective_k=effective_k)
    diagnostics = query_support(bank, z, groups, options, workers)
    reference = diagnostics["kth_distance"][calibration]
    bank["thresholds"] = {str(q):float(np.quantile(reference, q)) if q < 1 else float("inf") for q in COVERAGES}
    return bank, dict(status="complete", radius=radius, requested_radius_k=options["radius_k"],
                      effective_radius_k=effective_k, dimensions=projection["dimensions"],
                      captured_train_embedding_variance=projection["captured_train_embedding_variance"],
                      fit="TRAIN supervised embeddings; CALIBRATION X only for radius and thresholds",
                      uses_query_outcomes=False, borrows_neighbor_outcomes=False,
                      count="All eligible TRAIN participants in frozen projected-embedding radius; no K cap",
                      family_ess="Concentration of family counts among at most K within-radius matched neighbours",
                      interpretation="Exploratory support in a supervised embedding, not the outcome-independent primary ABM geometry")


def support_columns(frame, embeddings, cohort, bank, options, workers):
    if bank is None:
        return frame
    projection = bank["projection"]
    z = (np.asarray(embeddings, float)-projection["center"])@projection["components"]
    diagnostics = query_support(bank, z, cohort.group_id.to_numpy(str), options, workers)
    for key, value in diagnostics.items():
        frame[f"{METHOD}.{key}"] = value
    for q in COVERAGES:
        suffix = f"{q:.2f}".replace(".", "p")
        frame[f"{METHOD}.selected_support_{suffix}"] = diagnostics["supported"] & (diagnostics["kth_distance"] <= bank["thresholds"][str(q)])
    return frame


def atomic_torch(value, path):
    temporary = Path(path).with_name(Path(path).name+f".{os.getpid()}.partial")
    torch_api().save(value, temporary)
    temporary.replace(path)


def fit_block_deep(cohort, outcomes, predictions, features, feature_manifest,
                   output_dir, *, trait, config=None, seed=20261116, threads=1):
    options = settings(config)
    cohort = identified(cohort, "Deep roster")
    check_roles(cohort)
    predictions = aligned(predictions, cohort.id, "Deep baseline predictions")
    if len(features.shape) != 2 or features.shape[0] != len(cohort):
        raise ValueError("Regional matrix and roster rows differ")
    schema = feature_schema(feature_manifest, features.shape[1])
    outcomes = identified(outcomes, "Deep development outcomes")
    # Keep exactly the same information boundary as CPU ABM. In particular,
    # calibration and test values are not converted, validated or fitted.
    y = development_outcomes(outcomes, cohort)
    task = "binary" if trait == "t2dm" else "continuous"
    train = np.flatnonzero(cohort.role.eq("train").to_numpy() & np.isfinite(y))
    tune = np.flatnonzero(cohort.role.eq("tune").to_numpy() & np.isfinite(y))
    if len(train) < 10 or len(tune) < 10:
        raise ValueError("Deep fitting needs at least ten observed TRAIN and TUNE outcomes")
    if task == "binary":
        for rows in (train, tune):
            if not np.isin(y[rows], [0,1]).all() or len(np.unique(y[rows])) != 2:
                raise ValueError("Binary deep fitting needs cases and controls in TRAIN and TUNE")
    base = pd.to_numeric(predictions[PRIMARY+".prediction"], errors="raise").to_numpy(float)
    if not np.isfinite(base).all() or (task == "binary" and not np.all((base >= 0)&(base <= 1))):
        raise ValueError("Invalid frozen baseline prediction")
    observed_train = np.zeros(len(cohort), bool)
    observed_train[train] = True
    oof = check_oof(cohort, predictions, observed_train)
    if task == "binary" and not np.all((oof[train] >= 0)&(oof[train] <= 1)):
        raise ValueError("Invalid baseline OOF probabilities")
    prep = fit_transformer_input(features, train, options["batch_size"])
    scale = float(np.std(y[train]-oof[train])) if task == "continuous" else 1.
    scale = max(scale, 1e-8)
    torch, device = setup_device(options["device"], seed, threads)
    model = build_model(len(prep["keep"]), options).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=options["learning_rate"], weight_decay=options["weight_decay"])
    rng = np.random.default_rng(seed)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    best, best_state, best_epoch, history = np.inf, None, 0, []
    started = time.monotonic()
    for epoch in range(1, options["epochs"]+1):
        model.train()
        total, count = 0., 0
        for rows in batches(rng.permutation(train), options["batch_size"]):
            x = torch.as_tensor(transform(features, prep, rows), device=device)
            target = torch.as_tensor(y[rows], dtype=torch.float32, device=device)
            optimizer.zero_grad(set_to_none=True)
            delta, _ = model(x)
            if task == "binary":
                offset = torch.as_tensor(logit(np.clip(oof[rows],1e-8,1-1e-8)), dtype=torch.float32, device=device)
                loss = torch.nn.functional.binary_cross_entropy_with_logits(offset+delta, target)
            else:
                residual = torch.as_tensor((y[rows]-oof[rows])/scale, dtype=torch.float32, device=device)
                loss = torch.nn.functional.mse_loss(delta, residual)
            if not torch.isfinite(loss):
                raise RuntimeError("Nonfinite abm_deep training loss")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), options["gradient_clip"], error_if_nonfinite=True)
            optimizer.step()
            total += float(loss.detach())*len(rows)
            count += len(rows)
        delta, _ = predict_network(model, features, prep, tune, options["batch_size"])
        raw = combine(base[tune], delta, task, scale, 1.)
        value = float(primary_rows(y[tune], raw).mean())
        if value < best-options["min_delta"]:
            best, best_epoch = value, epoch
            best_state = {k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
            atomic_torch(dict(version=VERSION, state_dict=best_state), out/"model.pt")
        history.append(dict(epoch=epoch, train_loss=total/count, tune_loss=value, best_epoch=best_epoch))
        atomic_table(pd.DataFrame(history), out/"training.tsv")
        print(f"[abm.deep] {trait}: epoch={epoch}; tune_loss={value:.6g}; best={best_epoch}; device={device}", flush=True)
        if epoch-best_epoch >= options["patience"]:
            break
    if best_state is None:
        raise RuntimeError("Deep training produced no finite checkpoint")
    model.load_state_dict(best_state)
    delta, embeddings = predict_network(model, features, prep, np.arange(len(cohort)), options["batch_size"])
    alpha, choices = choose_alpha(y[tune], base[tune], delta[tune], cohort.group_id.to_numpy(str)[tune], task, scale, options)
    atomic_table(pd.DataFrame(choices), out/"tuning.tsv")
    bank, support_meta = fit_support(embeddings, cohort, train, options, threads) if options["support"] else (None,{"status":"disabled"})
    final = pd.DataFrame({"id":cohort.id, f"{METHOD}.prediction":combine(base, delta, task, scale, alpha),
                          f"{METHOD}.raw_prediction":combine(base, delta, task, scale, 1.)})
    final = support_columns(final, embeddings, cohort, bank, options, threads)
    bundle = dict(version=VERSION, model_file="model.pt", options=options, task=task, trait=trait,
                  prep=prep, feature_schema=schema, alpha=alpha, residual_scale=scale,
                  bank=bank, seed=seed, primary=PRIMARY)
    temporary = out/f"bundle.joblib.{os.getpid()}.partial"
    joblib.dump(bundle, temporary, compress=3)
    temporary.replace(out/"bundle.joblib")
    atomic_table(final, out/"predictions.tsv.gz")
    metadata = dict(version=VERSION, status="complete", trait=trait, task=task, method=METHOD,
        architecture="participant-wise consensus-block tokens; latent cross-attention; latent self-attention; CSx residual/offset head",
        baseline=PRIMARY, baseline_training_residual="family OOF only", alpha=alpha,
        gate="TUNE MSE/Brier improvement greater than gate_se times family-cluster SE; alpha0 exact baseline fallback",
        learned_embedding_borrows_labels=False, test_outcomes_read=False, calibration_outcomes_used=False,
        n_train=len(train), n_tune=len(tune), n_input_blocks=features.shape[1], n_variable_blocks=len(prep["keep"]),
        feature_normalization=prep["normalization"], residual_scale=scale, best_epoch=best_epoch,
        train_loss="baseline-offset binary cross-entropy" if task=="binary" else "normalized OOF residual MSE",
        tune_loss="Brier" if task=="binary" else "MSE", support=support_meta,
        requested_device=options["device"], actual_device=device, torch_version=str(torch.__version__),
        cuda_tested_in_this_run=device=="cuda", same_architecture_cpu_cuda=True, cross_device_numerical_parity_verified=False,
        parameter_count=sum(p.numel() for p in model.parameters()), elapsed_seconds=time.monotonic()-started,
        config=options)
    atomic_json(metadata, out/"model.json")
    return final, metadata


def predict_frozen(bundle_dir, cohort, predictions, features, feature_manifest, *,
                   device="cpu", batch_size=128, threads=1):
    location = Path(bundle_dir)
    if location.is_file():
        location = location.parent
    bundle = joblib.load(location/"bundle.joblib")
    if bundle.get("version") != VERSION:
        raise ValueError("Unknown deep bundle version")
    cohort = identified(cohort, "Deep query cohort")
    if "group_id" not in cohort or cohort.group_id.isna().any():
        raise ValueError("Query cohort needs known family/group IDs for support exclusion")
    predictions = aligned(predictions, cohort.id, "Deep query baseline")
    if features.shape[0] != len(cohort) or feature_schema(feature_manifest,features.shape[1]) != bundle["feature_schema"]:
        raise ValueError("Deep query features do not match the ordered training schema")
    options = {**bundle["options"], "batch_size":int(batch_size)}
    if batch_size < 1:
        raise ValueError("Prediction batch_size must be positive")
    torch, actual = setup_device(device, bundle["seed"], threads)
    model = build_model(len(bundle["prep"]["keep"]), options).to(actual)
    state = torch.load(location/bundle["model_file"], map_location="cpu", weights_only=True)
    if state["version"] != VERSION:
        raise ValueError("Unknown deep weights version")
    model.load_state_dict(state["state_dict"])
    base = pd.to_numeric(predictions[bundle["primary"]+".prediction"], errors="raise").to_numpy(float)
    if not np.isfinite(base).all() or (bundle["task"]=="binary" and not np.all((base >= 0)&(base <= 1))):
        raise ValueError("Invalid query baseline prediction")
    delta, embeddings = predict_network(model,features,bundle["prep"],np.arange(len(cohort)),batch_size)
    result = pd.DataFrame({"id":cohort.id,
        f"{METHOD}.prediction":combine(base,delta,bundle["task"],bundle["residual_scale"],bundle["alpha"]),
        f"{METHOD}.raw_prediction":combine(base,delta,bundle["task"],bundle["residual_scale"],1.)})
    return support_columns(result,embeddings,cohort,bundle["bank"],options,threads)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle",required=True)
    parser.add_argument("--cohort",required=True)
    parser.add_argument("--predictions",required=True)
    parser.add_argument("--features-npy",required=True)
    parser.add_argument("--features-ids",required=True)
    parser.add_argument("--features-manifest",required=True)
    parser.add_argument("--output",required=True)
    parser.add_argument("--device",choices=["cpu","gpu","cuda","auto"],default="cpu")
    parser.add_argument("--batch-size",type=int,default=128)
    parser.add_argument("--threads",type=int,default=1)
    args = parser.parse_args(argv)
    cohort = identified(read_table(args.cohort))
    ids = identified(read_table(args.features_ids))
    matrix = np.load(args.features_npy,mmap_mode="r")
    if len(matrix)!=len(ids) or set(cohort.id)-set(ids.id):
        raise ValueError("Query matrix ID roster mismatch")
    from abm2 import IndexedFeatures
    rows = pd.Series(np.arange(len(ids)),index=ids.id).loc[cohort.id].to_numpy()
    features = IndexedFeatures(matrix,rows)
    result = predict_frozen(args.bundle,cohort,read_table(args.predictions),features,
        read_table(args.features_manifest),device=args.device,batch_size=args.batch_size,threads=args.threads)
    atomic_table(result,args.output)


if __name__ == "__main__":
    main()
