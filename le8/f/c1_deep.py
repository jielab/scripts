#!/usr/bin/env python3
"""C1 Transformer training and legacy frozen-model GPU parity."""
from __future__ import annotations


# Transformer and frozen GPU parity
# Shared PNG byte limit; PDF/SVG export settings are preserved.
import sys as _png_sys
from pathlib import Path as _PngPath
_png_sys.path.append(str(_PngPath(__file__).resolve().parents[2] / "0f"))
from png_output import save_png_figure

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import time
import tempfile
import fcntl

import numpy as np
import pandas as pd
import torch
from torch import nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
import joblib

import abm

VERSION = "LE8-C1-deep-2026-10-22-v1"
HERE = Path(__file__).resolve().parent


def host_module():
    spec = importlib.util.spec_from_file_location("le8_gpu_host", HERE / "0.host.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def setup_device(seed, cores=2, device="auto"):
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


def cuda_setup(seed, cores=2):
    """Compatibility entry for historical frozen-CUDA parity experiments."""
    return setup_device(seed, cores, "cuda")


def synchronize(device):
    if str(device).startswith("cuda"):
        torch.cuda.synchronize()


def validate_config(c):
    if c["architecture"] != "feature_token_latent_transformer":
        raise ValueError("Unknown GPU architecture")
    for key in ("d_model", "heads", "latent_tokens", "layers", "batch_size", "epochs", "patience", "neighbors", "bootstrap"):
        if type(c[key]) is not int or c[key] < 1:
            raise ValueError(f"Invalid {key}")
    if c["d_model"] % c["heads"] or not 0 <= c["dropout"] < 1:
        raise ValueError("Invalid attention dimensions/dropout")
    if c["coverages"] != [0.2,0.4,0.5,0.6,0.8,1.0] or any(c[k] <= 0 for k in ("learning_rate", "gradient_clip")):
        raise ValueError("Invalid optimization/coverage settings")
    if c["weight_decay"] < 0 or c["min_delta"] < 0:
        raise ValueError("Invalid regularization/early stopping")

    for key in ("run_parity", "run_transformer"):
        if type(c[key]) is not bool:
            raise ValueError(f"Invalid {key}")
    if not (c["run_parity"] or c["run_transformer"]):
        raise ValueError("Enable at least one GPU branch")
    for key in ("parity_query_block", "parity_donor_block", "parity_benchmark_n", "parity_timing_repeats"):
        if type(c[key]) is not int or c[key] < 1:
            raise ValueError(f"Invalid {key}")
    for key, maximum in (("parity_probability_tolerance", 1e-8), ("parity_continuous_tolerance", 1e-7)):
        if not 0 < c[key] <= maximum:
            raise ValueError(f"Invalid {key}")
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
    spec = importlib.util.spec_from_file_location("le8_gpu_kernels", HERE.parents[1] / "0f/abm.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.BlockedKNN


CPU_FILES = ("development_LOCAL_ONLY.npz", "test_X_LOCAL_ONLY.npz",
             "test_truth_LOCAL_ONLY.npz", "fit/model_LOCAL_ONLY.joblib",
             "fit/MODEL_FROZEN.json", "host_contract.json", "COMPLETE.json",
             "test_predictions_LOCAL_ONLY.confidence.tsv.gz", "test_predictions_LOCAL_ONLY.joblib")


def source_receipt(source):
    files = {}
    for name in CPU_FILES:
        path = source / name
        if not path.is_file():
            raise ValueError(f"Missing CPU prerequisite {path}; run c1_abm first")
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(8*1024*1024), b""):
                h.update(block)
        files[name] = h.hexdigest()
    return dict(path=str(source.resolve()), files=files)


def load_cpu_inputs(source, args, cfg):
    """Read-only adapter for the existing CPU publication format; no raw inputs."""
    host = host_module()
    if not host.completed(source):
        raise ValueError(f"CPU ABM is incomplete: {source}; run c1_abm before c1_abm_GPU")
    host.check_input_reuse(source, host.host_config(args))
    receipt = source_receipt(source)
    cpu = abm.load_model(source / "fit/model_LOCAL_ONLY.joblib")
    from dataclasses import asdict
    ignored = {"cores", "parallel_workers", "always_include_indices", "grouping_unit"}
    changed = [k for k,v in asdict(cfg).items() if k not in ignored and getattr(cpu.config,k,None) != v]
    if changed:
        raise ValueError("CPU/GPU configuration mismatch: " + ",".join(changed))
    dev = abm.read_npz(source / "development_LOCAL_ONLY.npz")
    test = abm.read_npz(source / "test_X_LOCAL_ONLY.npz")
    truth = abm.read_npz(source / "test_truth_LOCAL_ONLY.npz")
    names = list(map(str,dev["feature_names"]))
    if (names != list(map(str,test["feature_names"])) or names != cpu.feature_names
        or not np.array_equal(test["eid"],truth["eid"])
        or not np.array_equal(dev["eid"].astype(str),cpu.development_ids)
        or not np.array_equal(dev["group"].astype(str),cpu.development_groups)):
        raise ValueError("CPU prepared schema/ID alignment mismatch")
    data = dict(X=np.concatenate([dev["X"],test["X"]]),
        time=np.concatenate([dev["time"],truth["time"]]),
        event=np.concatenate([dev["event"],truth["event"]]),
        eid=np.concatenate([dev["eid"],test["eid"]]),
        family_group=np.concatenate([dev["group"],test["group"]]),
        role=np.concatenate([dev["role"],np.repeat("test",len(test["eid"]))]),
        feature_names=names,clinical_feature_names=[names[i] for i in cpu.config.always_include_indices])
    host.check_group_roles(data["family_group"],data["role"])
    if len(set(data["eid"])) != len(data["eid"]) or set(dev["role"]) != set(abm.ROLES):
        raise ValueError("CPU development/test roles or IDs are invalid")
    audit = json.loads((source / "host_contract.json").read_text())
    validate_cpu_identity(cpu, dev, test, truth, audit, source)
    audit.update(cpu_source=receipt, reuse=["prepared inputs", "family roles", "feature panel",
        "fitted preprocessing", "fitted censoring model", "CPU predictions"],
        not_recomputed=["raw data preparation", "CPU crossfit", "CPU model fitting",
                        "CPU donor pruning", "molecular deviation", "other LE8 modules"])
    return data,cpu,audit,receipt


def lexicographic_tie_keys(ids):
    ids = np.asarray(ids, str)
    order = np.argsort(ids, kind="stable")
    rank = np.empty(len(ids), dtype=np.int64)
    rank[order] = np.arange(len(ids), dtype=np.int64)
    return rank

class GpuFrozenTwinBank:
    """CUDA executor for one already-fitted CPU TwinBank.
    Fits nothing and changes no scientific parameter.
    """

    def __init__(self, bank, query_block=256, donor_block=8192):
        if bank.radius is None:
            raise ValueError("Frozen CPU TwinBank has no calibrated radius")

        self.bank = bank
        self.radius = float(bank.radius)
        self.k = int(bank.config.k)
        self.prior_strength = float(bank.config.prior_strength)

        # CPU ultimately excludes donors with w<=0.
        self.original_donor_index = np.flatnonzero(np.asarray(bank.w) > 0)
        if not len(self.original_donor_index):
            raise ValueError("No eligible frozen Twin donors")

        donor_Z = np.asarray(bank.Z, float)[self.original_donor_index]
        donor_ids = np.asarray(bank.ids, str)[self.original_donor_index]
        donor_groups = np.asarray(bank.groups, str)[self.original_donor_index]

        tie = lexicographic_tie_keys(donor_ids)

        t0 = time.perf_counter()
        self.index = shared_knn()(
            donor_Z,
            donor_ids,
            donor_groups,
            tie,
            device="cuda",
            dtype="float64",
            query_block=int(query_block),
            donor_block=int(donor_block),
            tie_tolerance=0.0,
        )
        torch.cuda.synchronize()
        self.upload_seconds = time.perf_counter() - t0

    def predict(self, z, ids, groups, fallback):
        ids = np.asarray(ids, str)
        groups = np.asarray(groups, str)
        fallback = np.asarray(fallback, float)

        q = np.asarray(self.bank.pca.transform(z), float) / np.sqrt(
            self.bank.pca.n_components_
        )

        torch.cuda.synchronize()
        t0 = time.perf_counter()
        near, distances, counts = self.index.query(
            q,
            ids,
            groups,
            k=min(self.k, len(self.original_donor_index)),
            radius=self.radius,
            return_counts=True,
        )
        torch.cuda.synchronize()
        retrieval_seconds = time.perf_counter() - t0

        names = (
            "eligible_match_count","matched_count","family_count","family_ess",
            "donor_ess","donor_loss","radius_ratio","nearest_radius_ratio",
            "label_disagreement","borrow_fraction"
        )
        diag = {name: np.zeros(len(q), dtype=float) for name in names}
        diag["donor_loss"].fill(10.0)
        diag["radius_ratio"].fill(100.0)
        diag["nearest_radius_ratio"].fill(100.0)

        prob = fallback.copy()
        bank_y = np.asarray(self.bank.y, float)
        bank_w = np.asarray(self.bank.w, float)
        bank_loss = np.asarray(self.bank.losses, float)
        bank_group = np.asarray(self.bank.groups, str)

        post0 = time.perf_counter()

        for i in range(len(q)):
            ok = (
                (near[i] >= 0)
                & np.isfinite(distances[i])
                & (distances[i] <= self.radius)
            )
            rel = near[i][ok]
            ds = distances[i][ok]

            diag["eligible_match_count"][i] = counts[i]
            if not len(rel):
                continue

            idx = self.original_donor_index[rel]
            if len(idx) > self.k:
                idx = idx[:self.k]
                ds = ds[:self.k]

            weight = np.exp(-0.5 * np.square(ds / self.radius)) * bank_w[idx]
            den = weight.sum()
            if not np.isfinite(den) or den <= 0:
                continue
            weight /= den

            _, inv = np.unique(bank_group[idx], return_inverse=True)
            fw = np.bincount(inv, weights=weight)
            donor_ess = 1.0 / np.sum(np.square(weight))
            family_ess = 1.0 / np.sum(np.square(fw))
            local = float(weight @ bank_y[idx])
            borrowing = family_ess / (family_ess + self.prior_strength)

            prob[i] = borrowing * local + (1.0 - borrowing) * fallback[i]

            diag["matched_count"][i] = len(idx)
            diag["family_count"][i] = len(fw)
            diag["family_ess"][i] = family_ess
            diag["donor_ess"][i] = donor_ess
            diag["donor_loss"][i] = float(weight @ bank_loss[idx])
            diag["radius_ratio"][i] = float(ds[-1] / self.radius)
            diag["nearest_radius_ratio"][i] = float(ds[0] / self.radius)
            diag["label_disagreement"][i] = local * (1.0 - local)
            diag["borrow_fraction"][i] = borrowing

        post_seconds = time.perf_counter() - post0

        diag["supported"] = (
            (diag["matched_count"] >= self.bank.config.min_twins)
            & (diag["family_ess"] >= self.bank.config.min_family_ess)
        )

        return (
            np.clip(prob, 1e-7, 1-1e-7),
            diag,
            {
                "retrieval_seconds": retrieval_seconds,
                "postprocess_seconds": post_seconds,
                "n_queries": len(q),
            },
        )

def parity_selector(cpu_model, cpu_prediction, gpu_probability, gpu_diag):
    ext = cpu_model.final_predictability
    if ext is None:
        raise ValueError("CPU model lacks final predictability extension")
    if ext.primary != "risk_resid_error" or ext.audit.get("predictability_target_model") != "elasticnet_full":
        raise ValueError("LOCKPAPER primary selector changed")

    p = {
        k: np.asarray(v).copy()
        for k,v in cpu_prediction["probability"].items()
    }
    p.update({k: np.asarray(v) for k,v in gpu_probability.items()})

    scores = ext.scores(cpu_model, p, gpu_diag)
    eligible = (
        np.asarray(gpu_diag["technical_qc"], bool)
        & np.asarray(gpu_diag["supported"], bool)
    )

    masks = {}
    for q in cpu_model.config.coverages:
        key = f"{ext.primary}:{q:g}"
        masks[key] = ext.rules[(ext.primary,float(q))].apply(
            scores[ext.primary],
            eligible,
            np.asarray(cpu_prediction["eid"], str),
        )

    selected = masks[f"{ext.primary}:{cpu_model.config.primary_coverage:g}"]
    return scores, masks, selected

def validate_cpu_identity(cpu, dev, test, truth, audit, source):
    required = dict(grouping_unit="family", grouping_source="UKB KING kinship connected components",
                    family_disjoint_outer_split=True, family_disjoint_ABM_roles=True)
    if any(audit.get(k) != v for k, v in required.items()):
        raise ValueError("GPU requires verified KING family-disjoint CPU artifacts")
    from c1_protocol import sha, roster_record, grouping_audit
    record = roster_record(audit["outer_roster"])
    if record["roster_sha256"] != audit.get("outer_roster_sha256"):
        raise ValueError("CPU outer roster SHA changed")
    current = grouping_audit(record)
    if not audit.get("group_map_sha256") or current.get("group_map_sha256") != audit["group_map_sha256"]:
        raise ValueError("KING group-map SHA mismatch")
    group_file = audit["input_config"].get("group_file")
    if not group_file or sha(group_file) != audit["group_map_sha256"]:
        raise ValueError("KING group-map source changed")
    for role in abm.ROLES:
        if abm.fingerprint(dev["eid"][dev["role"].astype(str) == role].astype(str)) != cpu.training_audit["roles"][role]["id_hash"]:
            raise ValueError("CPU role identity mismatch")
    prediction = abm.load_model(source / "test_predictions_LOCAL_ONLY.joblib")
    if not np.array_equal(np.asarray(prediction["eid"], str), test["eid"].astype(str)) or not np.array_equal(np.asarray(prediction["group"], str), test["group"].astype(str)):
        raise ValueError("CPU prediction roster/groups mismatch")
    confidence = pd.read_csv(source / "test_predictions_LOCAL_ONLY.confidence.tsv.gz", sep="\t", dtype={"eid": str})
    if not np.array_equal(confidence.eid.to_numpy(), test["eid"].astype(str)):
        raise ValueError("CPU confidence roster mismatch")
    if cpu.config.seed != 2026 or cpu.config.bootstrap != 1000 or cpu.config.primary_coverage != .5 or list(cpu.config.coverages) != [.2,.4,.5,.6,.8,1.]:
        raise ValueError("Frozen CPU LOCKPAPER protocol mismatch")


def distance_selection(distances, qc, gate, coverages):
    valid = qc[gate] & np.isfinite(distances[gate])
    if not valid.any():
        raise ValueError("No eligible tune_gate distances")
    thresholds = {float(q): float(np.quantile(distances[gate][valid], q)) for q in coverages}
    masks = {q: qc & np.isfinite(distances) & (distances <= t) for q, t in thresholds.items()}
    return thresholds, masks


def risk_matched_evaluation(parts, masks, qc, p, ids, y, w, edges, seed):
    rows = []
    for role in ("calibration_audit", "test"):
        ix = parts[role]
        for q, mask in masks.items():
            selected = mask[ix]
            control = abm.risk_matched_mask(selected, qc[ix], p[ix], ids[ix], edges, seed)
            a = metrics(y[ix][selected], w[ix][selected], p[ix][selected])
            b = metrics(y[ix][control], w[ix][control], p[ix][control])
            rows.append(dict(role=role, population="selected_vs_risk_matched", target_coverage=q,
                coverage=float(selected.mean()), **a,
                **{f"risk_matched_{k}": b.get(k, np.nan) for k in ("brier", "logloss", "auc")},
                **{f"delta_{k}": a.get(k, np.nan)-b.get(k, np.nan) for k in ("brier", "logloss", "auc")},
                selector="embedding_neighbor_distance", exploratory=True, outcome_blind=True))
    return rows


def parity_comparison(reference, probability, diag, masks, selected, network):
    rows = []
    exact = {"eligible_match_count", "matched_count", "family_count", "supported", "technical_qc"}
    for kind, actual, expected in (("probability", probability, reference["probability"]),
                                    ("diagnostic", diag, reference["diagnostics"])):
        for name, value in actual.items():
            a, b = np.asarray(value), np.asarray(expected[name])
            delta = np.abs(a.astype(float)-b.astype(float))
            finite = bool(np.isfinite(delta).all())
            maximum = float(delta.max()) if len(delta) and finite else (0. if not len(delta) else float("inf"))
            tolerance = network["parity_probability_tolerance"] if kind == "probability" else network["parity_continuous_tolerance"]
            passed = np.array_equal(a,b) if name in exact else finite and maximum <= tolerance
            corr = float(np.corrcoef(a,b)[0,1]) if len(a)>1 and np.std(a)>0 and np.std(b)>0 else np.nan
            rows.append(dict(kind=kind, field=name, N=len(a), max_abs_difference=maximum,
                max_abs_probability_difference=maximum if kind == "probability" else np.nan,
                mean_abs_probability_difference=float(delta.mean()) if kind == "probability" else np.nan,
                correlation=corr, tolerance=0 if name in exact else tolerance, passed=bool(passed)))
    selector = [dict(mask=k, N=len(v), mismatches=int(np.count_nonzero(v != reference["masks"][k])),
                     passed=bool(np.array_equal(v,reference["masks"][k]))) for k,v in masks.items()]
    selector.append(dict(mask="selected", N=len(selected), mismatches=int(np.count_nonzero(selected != reference["selected"])),
                         passed=bool(np.array_equal(selected,reference["selected"]))))
    return rows, selector


def benchmark_bank(bank, executor, z, ids, groups, fallback, network, route):
    from threadpoolctl import threadpool_limits
    n = min(network["parity_benchmark_n"], len(z))
    if not n:
        raise ValueError("No QC queries for parity benchmark")
    args = (z[:n], ids[:n], groups[:n], fallback[:n])
    cpu_times, gpu_times, retrieval, post = [], [], [], []
    with threadpool_limits(limits=2):
        bank.predict(*args)
        executor.predict(*args)
        for _ in range(network["parity_timing_repeats"]):
            start = time.perf_counter()
            bank.predict(*args)
            cpu_times.append(time.perf_counter()-start)
            torch.cuda.synchronize()
            start = time.perf_counter()
            _, _, timing = executor.predict(*args)
            torch.cuda.synchronize()
            gpu_times.append(time.perf_counter()-start)
            retrieval.append(timing["retrieval_seconds"])
            post.append(timing["postprocess_seconds"])
    cpu_seconds, gpu_seconds = float(np.median(cpu_times)), float(np.median(gpu_times))
    return dict(route=route, query_n=n, query_ids_hash=abm.fingerprint(ids[:n]),
        bank_size=len(bank.Z), eligible_bank_size=len(executor.original_donor_index), dimension=bank.Z.shape[1],
        k=bank.config.k, radius=bank.radius, gpu_upload_seconds=executor.upload_seconds,
        cpu_end_to_end_seconds=cpu_seconds, gpu_end_to_end_seconds=gpu_seconds,
        gpu_retrieval_seconds=float(np.median(retrieval)), gpu_postprocess_seconds=float(np.median(post)),
        cpu_queries_per_second=n/cpu_seconds, gpu_queries_per_second=n/gpu_seconds, speedup=cpu_seconds/gpu_seconds)


def run_parity(cpu, data, source, out, network):
    out.mkdir(parents=True, exist_ok=True)
    test = np.asarray(data["role"]) == "test"
    X, ids, groups = (np.asarray(data[k])[test] for k in ("X", "eid", "family_group"))
    reference = abm.load_model(source / "test_predictions_LOCAL_ONLY.joblib")
    z, qc = cpu.prep.transform(X), cpu.prep.technical_qc(X)
    fallback = cpu.models["elasticnet_full"].predict_proba(z)[:,1]
    probability, timing_rows, diag = {}, [], None
    for prefix, banks in (("twin_", cpu.banks), ("twin_integrated_", cpu.banks_integrated)):
        for arm in ("full", "clean"):
            route, bank = prefix+arm, banks[arm]
            executor = GpuFrozenTwinBank(bank, network["parity_query_block"], network["parity_donor_block"])
            raw = fallback.copy()
            raw[qc], d, _ = executor.predict(z[qc], ids[qc], groups[qc], fallback[qc])
            probability[route] = cpu.calibrators[route].predict(raw)
            if prefix == "twin_":
                alias = "twin_molecular_"+arm
                probability[alias] = cpu.calibrators[alias].predict(raw)
                timing_rows.append(benchmark_bank(bank, executor, z[qc], ids[qc], groups[qc], fallback[qc], network, route))
                if arm == "clean":
                    diag = {k: np.zeros(len(X), dtype=v.dtype) for k,v in d.items()}
                    for k,v in d.items():
                        diag[k][qc] = v
                    for k,v in (("donor_loss",10.),("radius_ratio",100.),("nearest_radius_ratio",100.)):
                        diag[k][~qc] = v
            del executor
    diag.update(technical_qc=qc, molecular_missing_fraction=cpu.prep.missing(X),
                frozen_molecular_missing_fraction=cpu.prep.missing(X,panel=False))
    scores, masks, selected = parity_selector(cpu, reference, probability, diag)
    rows, selector = parity_comparison(reference, probability, diag, masks, selected, network)
    pd.DataFrame(rows).to_csv(out / "numerical_parity.tsv", sep="\t", index=False)
    pd.DataFrame(selector).to_csv(out / "selector_parity.tsv", sep="\t", index=False)
    pd.DataFrame(timing_rows).to_csv(out / "timing.tsv", sep="\t", index=False)
    abm.json_write(out / "timing_summary.json", dict(label="Twin matching GPU speedup", rows=timing_rows,
        repeats=network["parity_timing_repeats"], aggregation="median", warmup=True, cpu_threads=2,
        scope="frozen geometry transform, query transfer, exact search, return transfer, donor borrowing and diagnostics; bank upload reported separately"))
    if not all(r["passed"] for r in rows+selector):
        raise RuntimeError("Frozen CPU/GPU parity failed; see numerical_parity.tsv and selector_parity.tsv")
    joblib.dump(dict(eid=ids, group=groups, probability=probability, diagnostics=diag,
                     selector_scores=scores, masks=masks, selected=selected), out / "parity_predictions_LOCAL_ONLY.joblib", compress=3)
    report = dict(algorithm_equivalent=True, numerical_parity="passed",
        matching_space="frozen CPU molecular-only geometry; integrated sensitivity also checked",
        selector="frozen CPU risk_resid_error", test_outcomes_used_for_selection=False)
    report.update(parity_summary(rows, selector, timing_rows))
    abm.json_write(out / "PARITY_COMPLETE.json", dict(version=VERSION, **report))
    return report


def train(data, cfg, network, out, cpu_model=None, *, device="auto", report=True):
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
    prep = cpu_model.prep if cpu_model is not None else abm.NumericPrep().fit(X[b], cfg)
    qc = prep.technical_qc(X)
    km = cpu_model.km if cpu_model is not None else abm.CensoringKM().fit(data["time"][b][qc[b]], data["event"][b][qc[b]], cfg.horizon, cfg.min_censor_survival)
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
        raise RuntimeError("GPU model probability calibration did not converge")
    p = calibrator.predict_proba(score[:, None])[:, 1]
    # Exact CUDA retrieval with molecular-only embeddings and family/self exclusion.
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
        roles={k: len(v) for k,v in parts.items()}, preprocessing="reused CPU build fit" if cpu_model is not None else "build only", censoring="reused CPU build fit" if cpu_model is not None else "build only",
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
    if not report:
        pd.DataFrame(private).to_csv(out / "test_predictions_LOCAL_ONLY.tsv.gz", sep="\t", index=False)
        return dict(probability=p, diagnostics=diagnostics, technical_qc=qc, selected=selected, masks=masks,
                    donor_distance=distances, eligible_match_count=eligible_match_count,
                    eligible_reference_count=eligible_reference_count, matching_radius=support_radius,
                    thresholds=thresholds, parts=parts, freeze_path=str(out / "MODEL_FROZEN.json"),
                    device=device, model=model, prep=prep, km=km, calibrator=calibrator,
                    training_case_N=int(np.sum(np.asarray(data["event"])[b][qc[b]] == 1)),
                    incident_cases_removed_by_loss=0)
    rows = []
    for role in ("calibration_audit", "test"):
        ix = parts[role]
        rows.append(dict(role=role, population="all_QC", **metrics(y[ix], w[ix], p[ix])))
        use = ix[selected[ix]]
        rows.append(dict(role=role, population="selected", coverage=float(selected[ix].mean()),
                         **metrics(y[use], w[use], p[use])))
    rows.extend(risk_matched_evaluation(parts, masks, qc, p, ids, y, w, risk_edges, cfg.seed))
    perf = out / "performance"
    perf.mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(perf / "metrics.tsv", sep="\t", index=False)
    test = parts["test"]
    pd.DataFrame(dict(private, time=np.asarray(data["time"])[test],
                      event=np.asarray(data["event"])[test], y=y[test], weight=w[test])).to_csv(
            out / "test_predictions_LOCAL_ONLY.tsv.gz", sep="\t", index=False)
    plot_performance(out)
    # Compare core prediction on identical rows/panel, excluding preparation,
    # calibration, neighbor search and disk writes. These are different models.
    bench = parts["test"][:min(2048,len(parts["test"]))]
    predict(model,x,missing,bench,network["batch_size"])
    synchronize(device)
    started = time.perf_counter()
    predict(model,x,missing,bench,network["batch_size"])
    synchronize(device)
    gpu_seconds = time.perf_counter()-started
    timing = [dict(step="Transformer training and validation",seconds=training_seconds,n=len(train_idx)),
              dict(step=f"Transformer {device} core inference",seconds=gpu_seconds,n=len(bench))]
    if cpu_model is not None:
        from threadpoolctl import threadpool_limits
        # CPU models may have per-estimator n_jobs; do not mutate frozen CPU artifacts.
        cpu_z = cpu_model.prep.transform(X[bench])
        for name,estimator in cpu_model.models.items():
            with threadpool_limits(limits=2):
                estimator.predict_proba(cpu_z)
                started = time.perf_counter()
                estimator.predict_proba(cpu_z)
                elapsed = time.perf_counter()-started
            timing.append(dict(step=name+" CPU core inference",seconds=elapsed,n=len(bench)))
    pd.DataFrame(timing).to_csv(perf / "step_timings.tsv",sep="\t",index=False)
    abm.json_write(perf / "timing_scope.json",dict(
        comparison="Different models on identical frozen input panel and test rows; not equivalent algorithm speedup",
        inference_rows_hash=abm.fingerprint(ids[bench]),inference_n=len(bench),
        warmup=True,actual_device=device,CUDA_synchronized=device == "cuda",GPU_transfer_included=device == "cuda",
        CPU_threads=2,
        GPU_host_threads=torch.get_num_threads(),
        CPU_training_seconds=None,CPU_training_note="CPU training reused; no retrospective training speedup claimed"))
    return rows


def compare_cpu(target, cpu, seed, bootstrap):
    """Same participant/outcome intersection, same IPCW, paired family bootstrap.
    Merely reading completed CPU predictions never changes either model.
    """
    perf = target / "performance"
    perf.mkdir(exist_ok=True)
    source = cpu / "test_predictions_LOCAL_ONLY.confidence.tsv.gz"
    truth = cpu / "test_truth_LOCAL_ONLY.npz"
    if not source.exists() or not truth.exists():
        abm.json_write(perf / "cpu_comparison_status.json", dict(status="pending_CPU_predictions",
            instruction="Rerun c1_abm_GPU after CPU completes; completed GPU training is reused"))
        return
    g = pd.read_csv(target / "test_predictions_LOCAL_ONLY.tsv.gz", sep="\t", dtype={"eid": str, "family": str})
    c = pd.read_csv(source, sep="\t", dtype={"eid": str})
    t = np.load(truth)
    ct = pd.DataFrame(dict(eid=t["eid"].astype(str), cpu_time=t["time"], cpu_event=t["event"]))
    risk_columns = [k for k in ("fallback_risk", "candidate_risk", "policy_risk") if k in c]
    if "fallback_risk" not in risk_columns:
        raise ValueError("CPU reference predictions lack elasticnet_full risk")
    d = g.merge(c[["eid", *risk_columns]], on="eid", validate="one_to_one").merge(ct, on="eid", validate="one_to_one")
    if len(d) != len(g) or len(d) != len(c) or not len(d) or not np.array_equal(d.event, d.cpu_event) or not np.allclose(d.time, d.cpu_time, rtol=0, atol=1e-8):
        raise ValueError("CPU/GPU endpoint mismatch or no paired subjects")
    y,w = d.y.to_numpy(),d.weight.to_numpy()
    gpu,base = d.probability.to_numpy(),d.fallback_risk.to_numpy()
    gm,cm = metrics(y,w,gpu),metrics(y,w,base)
    labels = dict(fallback_risk="CPU_elasticnet_full", candidate_risk="CPU_primary_candidate", policy_risk="CPU_fallback_policy")
    rows=[dict(model="transformer_GPU",**gm)] + [dict(model=labels[k], **metrics(y,w,d[k].to_numpy())) for k in risk_columns]
    pd.DataFrame(rows).to_csv(perf / "paired_CPU_GPU.tsv", sep="\t", index=False)
    _, codes = np.unique(d.family, return_inverse=True)
    n = codes.max()+1
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(bootstrap):
        counts = np.bincount(rng.integers(n, size=n), minlength=n)
        ww = w*counts[codes]
        a,b = metrics(y,ww,gpu),metrics(y,ww,base)
        if a["status"] == b["status"] == "ok":
            deltas.append({k:a[k]-b[k] for k in ("auc","brier","logloss")})
    intervals = {k:np.quantile([v[k] for v in deltas],[.025,.975]).tolist() for k in deltas[0]} if deltas else {}
    abm.json_write(perf / "cpu_comparison_status.json", dict(status="paired_exploratory", n=len(d),
        gpu_test_n=len(g), cpu_test_n=len(c), family_bootstrap_valid=len(deltas),
        difference="GPU minus CPU", intervals_95=intervals, weights="same frozen CPU build-derived IPCW censoring model reused for both CPU and GPU predictions",
        CPU_model="elasticnet_full", confirmatory=False))

    return {"delta_" + k: finite_value(gm.get(k, np.nan) - cm.get(k, np.nan))
            for k in ("auc", "brier", "logloss")}


def plot_performance(target):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve
    history = pd.read_csv(target / "training_history.tsv", sep="\t")
    test = pd.read_csv(target / "test_predictions_LOCAL_ONLY.tsv.gz", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(history.epoch, history.tune_logloss)
    axes[0].set(xlabel="Epoch", ylabel="Tuning IPCW log loss", title="Development only")
    ok = test.weight > 0
    if test.loc[ok, "y"].nunique() == 2:
        fpr, tpr, _ = roc_curve(test.loc[ok, "y"], test.loc[ok, "probability"], sample_weight=test.loc[ok, "weight"])
        axes[1].plot(fpr, tpr, label="abm_deep")
        axes[1].legend()
    axes[1].plot([0,1], [0,1], "--", color="gray")
    axes[1].set(xlabel="False positive rate", ylabel="True positive rate", title="Exploratory holdout ROC")
    fig.tight_layout()
    save_png_figure(fig, target / "performance/training_and_ROC.png", dpi=160)
    plt.close(fig)


# Explicit aggregate-file allowlist: never discover review files recursively.
SUMMARY_SOURCES = {
    "parity_numerical": "parity/numerical_parity.tsv",
    "parity_selector": "parity/selector_parity.tsv",
    "parity_timing": "parity/timing.tsv",
    "transformer_metrics": "transformer/performance/metrics.tsv",
    "transformer_cpu_comparison": "transformer/performance/paired_CPU_GPU.tsv",
    "transformer_training": "transformer/training_history.tsv",
}
REVIEW_SOURCES = (
    "GPU_REPORT.json", "GPU_results.xlsx", "host_contract.json",
    *SUMMARY_SOURCES.values(), "parity/timing_summary.json", "parity/PARITY_COMPLETE.json",
    "transformer/performance/cpu_comparison_status.json",
    "transformer/performance/step_timings.tsv", "transformer/performance/timing_scope.json",
    "transformer/TRANSFORMER_COMPLETE.json",
)


def finite_value(value):
    return float(value) if value is not None and np.isfinite(value) else None


def parity_summary(rows, selectors, timings):
    probabilities = [r for r in rows if r["kind"] == "probability"]
    result = dict(probability_parity_passed=all(r["passed"] for r in probabilities) if probabilities else None,
                  selector_parity_passed=all(r["passed"] for r in selectors) if selectors else None,
                  device=torch.cuda.get_device_name(0) if torch.cuda.is_available() else None)
    for row in timings:
        route = row["route"]
        if route not in ("twin_full", "twin_clean"):
            continue
        for label, source in (("speedup", "speedup"), ("cpu_seconds", "cpu_end_to_end_seconds"),
                              ("gpu_seconds", "gpu_end_to_end_seconds")):
            result[route + "_" + label] = finite_value(row.get(source))
    counts = {r["query_n"] for r in timings}
    result["benchmark_query_N"] = int(next(iter(counts))) if len(counts) == 1 else None
    return result


def transformer_summary(rows, comparison):
    """Copy existing TEST aggregates; do not refit or recompute inference."""
    result = {}
    for row in rows:
        if row.get("role") != "test":
            continue
        population = row.get("population")
        prefix = {"all_QC": "all_QC", "selected": "selected_50"}.get(population)
        if prefix:
            for key, label in (("auc", "AUC"), ("brier", "Brier"), ("logloss", "LogLoss")):
                result[prefix + "_" + label] = finite_value(row.get(key))
        if population == "selected":
            result["coverage_50"] = finite_value(row.get("coverage"))
        if population == "selected_vs_risk_matched" and row.get("target_coverage") == .5:
            for key, label in (("auc", "AUC"), ("brier", "Brier"), ("logloss", "LogLoss")):
                result["delta_" + label + "_riskmatched_50"] = finite_value(row.get("delta_" + key))
    for key, label in (("auc", "AUC"), ("brier", "Brier"), ("logloss", "LogLoss")):
        result["paired_delta_" + label + "_vs_CPU_elasticnet_full"] = finite_value((comparison or {}).get("delta_" + key))
    return result


def write_gpu_summary(target):
    target = Path(target)
    tables = {}
    for name, relative in SUMMARY_SOURCES.items():
        path = target / relative
        if path.is_file():
            table = pd.read_csv(path, sep="\t")
            tables[name] = table
    report = json.loads((target / "GPU_REPORT.json").read_text())
    tables["summary"] = pd.DataFrame([
        dict(branch=branch, **report.get(branch, {})) for branch in ("parity", "transformer")])
    if "timing" in report:
        tables["module_timing"] = pd.DataFrame([report["timing"]])
    host_module()  # Registers the existing safe inline-string OOXML writer.
    from le8_common import write_xlsx
    write_xlsx(target / "GPU_results.xlsx", tables)


def write_review_manifest(target):
    target = Path(target)
    files = [dict(path=name, sha256=hashlib.sha256((target / name).read_bytes()).hexdigest())
             for name in REVIEW_SOURCES if (target / name).is_file()]
    abm.json_write(target / "REVIEW_FILES.json", dict(
        purpose="Small aggregate files sufficient for CPU/GPU scientific and timing review", files=files))


def finish_gpu(host, stage, audit):
    host.finish(stage, audit)
    write_gpu_summary(stage)
    write_review_manifest(stage)
    files = {}
    for path in stage.rglob("*"):
        if path.is_file() and path.name != "COMPLETE.json":
            h = hashlib.sha256()
            with path.open("rb") as stream:
                for block in iter(lambda: stream.read(8*1024*1024), b""):
                    h.update(block)
            files[str(path.relative_to(stage))] = h.hexdigest()
    abm.json_write(stage / "COMPLETE.json", dict(status="completed", version=VERSION, files=files))


def main():
    module_started = time.perf_counter()
    os.umask(0o077)
    host = host_module()
    parser = host.parser()
    parser.add_argument("--gpu-config", default=str(HERE.parent / "config/abm_GPU.json"))
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    args = parser.parse_args()
    cfg = abm.load_config(args.config)
    from dataclasses import replace
    cfg = replace(cfg, trait=args.Y, layer=args.biom, seed=args.seed, cores=args.cores)
    cfg.validate()
    network = json.loads(Path(args.gpu_config).read_text())
    validate_config(network)
    device = setup_device(args.seed, args.cores, args.device)
    if device == "cuda" and not args.preflight and os.getenv("LE8_GPU_LEASE_HELD") != "1":
        raise RuntimeError("CUDA lease required; launch c1_correlate --methods abm_deep --device cuda")
    if network["run_parity"] and device != "cuda":
        network = dict(network, run_parity=False)
        print("[LE8 abm_deep] Legacy CUDA parity skipped on CPU; scientific Transformer is unchanged", flush=True)
    if args.preflight:
        # CPU prerequisites may be produced earlier in this same pipeline.
        print("[LE8 GPU] Preflight passed; requires completed c1_abm artifacts at execution", flush=True)
        return
    target = args.analysis_root / args.Y / args.biom / "c1_abm_GPU"
    target.parent.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(str(target.resolve()).encode()).hexdigest()[:20]
    with open(f"/tmp/le8-gpu-{key}.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        source = target.parent / "c1_abm/abm_twin"
        data, cpu_model, audit, receipt = load_cpu_inputs(source,args,cfg)
        report_path = target / "GPU_REPORT.json"
        compatible = report_path.is_file() and json.loads(report_path.read_text()).get("version") == VERSION
        if host.completed(target) and compatible and args.replace != "TRUE":
            host.check_input_reuse(target, host.host_config(args))
            saved = json.loads((target / "host_contract.json").read_text())
            if {k:v for k,v in saved.get("cpu_source", {}).items() if k != "code_files"} != receipt or saved.get("gpu_config") != network:
                raise ValueError("CPU source or GPU configuration changed; use --replace TRUE")
            print("[LE8 GPU] SKIP verified completed branches", flush=True)
            return
        work = Path(tempfile.mkdtemp(prefix="_abm_GPU_staging_", dir=target.parent))
        try:
            print("[LE8 GPU] REUSE completed CPU inputs, preprocessing and censoring model",flush=True)
            report = dict(version=VERSION, cpu_reference="frozen c1_abm/abm_twin",
                          parity={"status": "disabled"}, transformer={"status": "disabled"})
            parity_seconds = transformer_seconds = None
            from threadpoolctl import threadpool_limits
            with threadpool_limits(limits=2):
                if network["run_parity"]:
                    started = time.perf_counter()
                    report["parity"] = run_parity(cpu_model, data, source, work / "parity", network)
                    parity_seconds = time.perf_counter() - started
                if network["run_transformer"]:
                    started = time.perf_counter()
                    rows = train(data, cfg, network, work / "transformer", cpu_model=cpu_model, device=device)
                    comparison = compare_cpu(work / "transformer", source, args.seed, network["bootstrap"])
                    report["transformer"] = dict(model="gpu_native_transformer", algorithm_equivalent=False,
                        algorithm_equivalent_to_CPU=False, exploratory=True, confirmatory=False,
                        matching_space="molecular-only learned embedding",
                        risk_model="molecular embedding + baseline clinical predictors",
                        selector="embedding_neighbor_distance", test_outcomes_used_for_selection=False)
                    report["transformer"].update(status="complete", **transformer_summary(rows, comparison))
                    transformer_seconds = time.perf_counter() - started
                    abm.json_write(work / "transformer/TRANSFORMER_COMPLETE.json", dict(version=VERSION, **report["transformer"]))
            report["timing"] = dict(module_wall_seconds=time.perf_counter() - module_started,
                parity_wall_seconds=parity_seconds, transformer_wall_seconds=transformer_seconds,
                scope="Module start through branch completion; excludes final report packaging and publication. Workflow runtime, not algorithm speedup.")
            abm.json_write(work / "GPU_REPORT.json", report)
            audit.update(version=VERSION, gpu_config=network)
            if source_receipt(source) != receipt:
                raise RuntimeError("CPU artifacts changed during GPU experiment; refusing publication")
            finish_gpu(host, work, audit)
            host.publish(work, target)
            print(f"[LE8 GPU] DONE {target}", flush=True)
        except BaseException:
            print(f"[LE8 GPU] Failed; diagnostics retained at {work}", flush=True)
            raise


if __name__ == "__main__":
    import sys
    sys.exit(main())
