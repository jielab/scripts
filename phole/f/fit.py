#!/usr/bin/env python3
"""Fit the PHOLE cis + alpha*beta + sparse-gamma graph decomposition."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    atomic_dataframe,
    atomic_json,
    chromosome_order,
    cis_rows,
    log,
    metric_dict,
    normalize_chr,
    parse_csv_numbers,
    position_kernel,
    randomized_svd,
    robust_center_scale,
    sample_trans_pairs,
    str2bool,
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--train-frac", type=float, default=0.80)
    parser.add_argument("--n-train", type=int, default=None)
    parser.add_argument("--split", choices=["random", "chromosome"], default="random")
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--rank", default="auto")
    parser.add_argument("--rank-grid", default="2,4,8,12,16,24")
    parser.add_argument("--inner-frac", type=float, default=0.20)
    parser.add_argument("--cis-flank", type=int, default=100_000)
    parser.add_argument("--gamma-z", type=float, default=4.0)
    parser.add_argument("--gamma-max-edges", type=int, default=1_000_000)
    parser.add_argument("--bandwidths", default="100000,1000000,10000000,100000000")
    parser.add_argument("--global-weight", type=float, default=0.02)
    parser.add_argument("--metric-pairs", type=int, default=1_000_000)
    parser.add_argument("--replace", default="FALSE")
    return parser.parse_args()


def resolve_n_train(n: int, n_train: int | None, train_frac: float) -> int:
    if n < 3:
        raise SystemExit(f"ERROR: at least three retained proteins are required, got {n}")
    if n_train is not None:
        if not 2 <= n_train < n:
            raise SystemExit(f"ERROR: --n-train must be in [2,{n - 1}], got {n_train}")
        return int(n_train)
    if not 0 < train_frac < 1:
        raise SystemExit(f"ERROR: --train-frac must be in (0,1), got {train_frac}")
    resolved = int(np.floor(n * train_frac))
    if not 2 <= resolved < n:
        raise SystemExit(
            f"ERROR: --train-frac {train_frac} gives {resolved} training proteins out of {n}; "
            "use a compatible fraction or --n-train"
        )
    return resolved


def split_targets(target: pd.DataFrame, n_train: int, method: str, seed: int) -> tuple[np.ndarray, np.ndarray]:
    n = len(target)
    if not 2 <= n_train < n:
        raise SystemExit(f"ERROR: --n-train must be in [2,{n - 1}], got {n_train}")
    rng = np.random.default_rng(seed)
    if method == "random":
        order = rng.permutation(n)
        train = np.sort(order[:n_train])
        test = np.sort(order[n_train:])
        return train, test

    # Whole-chromosome cold-start split. Dynamic programming picks a chromosome
    # subset with a target count closest to n - n_train.
    desired_test = n - n_train
    chrom_groups = {
        chrom: np.asarray(index, dtype=int)
        for chrom, index in target.groupby(target["CHR"].map(normalize_chr)).groups.items()
    }
    states: dict[int, list[str]] = {0: []}
    for chrom in sorted(chrom_groups, key=chromosome_order):
        count = len(chrom_groups[chrom])
        additions = {total + count: picked + [chrom] for total, picked in states.items()}
        for total, picked in additions.items():
            if total not in states:
                states[total] = picked
    best = min((x for x in states if 0 < x < n), key=lambda x: abs(x - desired_test))
    test_chrom = set(states[best])
    test = np.flatnonzero(target["CHR"].map(normalize_chr).isin(test_chrom).to_numpy())
    train = np.setdiff1d(np.arange(n), test)
    return train, test


def fit_components(
    w: np.ndarray,
    source: pd.DataFrame,
    target: pd.DataFrame,
    rank: int,
    cis_flank: int,
    seed: int,
) -> dict[str, np.ndarray]:
    x = np.asarray(w, dtype=np.float32)
    finite_count = np.sum(np.isfinite(x), axis=1)
    row_base = np.divide(
        np.nansum(x, axis=1, dtype=np.float64),
        finite_count,
        out=np.zeros(x.shape[0], dtype=np.float64),
        where=finite_count > 0,
    ).astype(np.float32)
    residual = x - row_base[:, None]
    residual = np.where(np.isfinite(residual), residual, 0.0).astype(np.float32)
    col_base = np.mean(residual, axis=0).astype(np.float32)
    residual -= col_base[None, :]

    src_chr = source["CHR"].map(normalize_chr).to_numpy(object)
    src_pos = pd.to_numeric(source["POS"], errors="coerce").to_numpy(float)
    cis_alpha = np.zeros(x.shape[1], dtype=np.float32)
    for j in range(x.shape[1]):
        mask = cis_rows(src_chr, src_pos, target.iloc[j]["CHR"], target.iloc[j]["POS"], cis_flank)
        vals = residual[mask, j]
        vals = vals[np.isfinite(vals)]
        if vals.size:
            cis_alpha[j] = float(np.median(vals))
            residual[mask, j] -= cis_alpha[j]

    u, singular, vt = randomized_svd(residual, rank=rank, seed=seed, n_iter=2, oversample=8)
    alpha = (u * singular[None, :]).astype(np.float32)
    beta = vt.astype(np.float32)
    return {
        "row_base": row_base,
        "col_base": col_base,
        "cis_alpha": cis_alpha,
        "alpha": alpha,
        "beta": beta,
        "singular": singular.astype(np.float32),
        "residual_before_factor": residual,
    }


def predict_components(
    model: dict[str, np.ndarray],
    source: pd.DataFrame,
    train_meta: pd.DataFrame,
    query_meta: pd.DataFrame,
    bandwidths: list[float],
    global_weight: float,
    cis_flank: int,
    rank: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    k = model["alpha"].shape[1] if rank is None else int(rank)
    kernel = position_kernel(train_meta, query_meta, bandwidths, global_weight)
    col = kernel @ model["col_base"].astype(float)
    cis = kernel @ model["cis_alpha"].astype(float)
    beta_query = model["beta"][:k, :].astype(float) @ kernel.T
    predicted = model["row_base"].astype(float)[:, None] + col[None, :]
    predicted += model["alpha"][:, :k].astype(float) @ beta_query

    src_chr = source["CHR"].map(normalize_chr).to_numpy(object)
    src_pos = pd.to_numeric(source["POS"], errors="coerce").to_numpy(float)
    for j in range(len(query_meta)):
        mask = cis_rows(src_chr, src_pos, query_meta.iloc[j]["CHR"], query_meta.iloc[j]["POS"], cis_flank)
        predicted[mask, j] += cis[j]
    return predicted.astype(np.float32), beta_query.astype(np.float32), col.astype(np.float32), cis.astype(np.float32)


def sampled_prediction(
    model: dict[str, np.ndarray],
    beta_query: np.ndarray,
    col_query: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    rank: int,
) -> np.ndarray:
    pred = model["row_base"][rows].astype(float) + col_query[cols].astype(float)
    pred += np.sum(
        model["alpha"][rows, :rank].astype(float) * beta_query[:rank, cols].T.astype(float), axis=1
    )
    return pred


def choose_rank(
    matrix: np.ndarray,
    source: pd.DataFrame,
    target: pd.DataFrame,
    train_index: np.ndarray,
    rank_grid: list[int],
    inner_frac: float,
    bandwidths: list[float],
    global_weight: float,
    cis_flank: int,
    metric_pairs: int,
    seed: int,
) -> tuple[int, pd.DataFrame]:
    if not 0.05 <= inner_frac <= 0.5:
        raise SystemExit("ERROR: --inner-frac must be between 0.05 and 0.5")
    rng = np.random.default_rng(seed + 11)
    shuffled = rng.permutation(train_index)
    n_valid = max(2, int(round(len(shuffled) * inner_frac)))
    inner_valid = np.sort(shuffled[:n_valid])
    inner_train = np.sort(shuffled[n_valid:])
    max_rank = min(max(rank_grid), len(inner_train) - 1, len(source) - 1)
    rank_grid = sorted({x for x in rank_grid if 1 <= x <= max_rank})
    if not rank_grid:
        raise SystemExit("ERROR: no usable ranks in --rank-grid")
    log(f"Internal cold-start rank selection: train={len(inner_train)} valid={len(inner_valid)} ranks={rank_grid}")
    inner_model = fit_components(
        matrix[:, inner_train], source, target.iloc[inner_train].reset_index(drop=True), max(rank_grid), cis_flank, seed + 17
    )
    _, beta_valid, col_valid, _ = predict_components(
        inner_model, source, target.iloc[inner_train].reset_index(drop=True),
        target.iloc[inner_valid].reset_index(drop=True), bandwidths, global_weight, cis_flank,
    )
    row_pair, col_pair = sample_trans_pairs(
        source, target.iloc[inner_valid].reset_index(drop=True), metric_pairs, cis_flank, seed + 19
    )
    observed = np.asarray(matrix[:, inner_valid])[row_pair, col_pair]
    rows = []
    for rank in rank_grid:
        predicted = sampled_prediction(inner_model, beta_valid, col_valid, row_pair, col_pair, rank)
        metrics = metric_dict(observed, predicted)
        rows.append({"RANK": rank, **metrics})
    results = pd.DataFrame(rows)
    finite = results[np.isfinite(results["R2"])]
    if finite.empty:
        best = int(results.loc[results["RMSE"].idxmin(), "RANK"])
    else:
        best = int(finite.loc[finite["R2"].idxmax(), "RANK"])
    log(f"Selected rank={best} by internal trans-only R2")
    return best, results


def save_npy_atomic(path: Path, array: np.ndarray) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.npy")
    np.save(tmp, array)
    os.replace(tmp, path)


def main() -> None:
    args = arguments()
    matrix_dir = Path(args.matrix_dir)
    out = Path(args.out)
    replace = str2bool(args.replace)
    if (out / "fit.done").exists() and (out / "model.npz").exists() and not replace:
        log(f"PHOLE fit exists: {out / 'model.npz'}")
        return
    for required in ["score.npy", "source.tsv.gz", "target.tsv.gz", "manifest.json"]:
        if not (matrix_dir / required).exists():
            raise SystemExit(f"ERROR: missing prepared input: {matrix_dir / required}")
    out.mkdir(parents=True, exist_ok=True)
    matrix = np.load(matrix_dir / "score.npy", mmap_mode="r")
    source = pd.read_csv(matrix_dir / "source.tsv.gz", sep="\t", dtype={"SOURCE": str, "CHR": str})
    target = pd.read_csv(
        matrix_dir / "target.tsv.gz", sep="\t", dtype={"PROTEIN": str, "GENE": str, "CHR": str}
    )
    if matrix.shape != (len(source), len(target)):
        raise SystemExit(f"ERROR: matrix/metadata shape mismatch: {matrix.shape} vs {(len(source), len(target))}")
    bandwidths = parse_csv_numbers(args.bandwidths, float)
    if any(x <= 0 for x in bandwidths):
        raise SystemExit("ERROR: all --bandwidths must be positive")
    if args.global_weight < 0:
        raise SystemExit("ERROR: --global-weight must be non-negative")

    n_train = resolve_n_train(len(target), args.n_train, args.train_frac)
    train_index, test_index = split_targets(target, n_train, args.split, args.seed)
    split = target[["TARGET_INDEX", "PROTEIN", "GENE", "CHR", "POS"]].copy()
    split["SET"] = "test"
    split.loc[train_index, "SET"] = "train"
    atomic_dataframe(split, out / "split.tsv.gz")
    log(f"Cold-start split: train={len(train_index)} test={len(test_index)} method={args.split}")

    max_possible = min(len(source) - 1, len(train_index) - 1)
    if args.rank.lower() == "auto":
        grid = [int(x) for x in parse_csv_numbers(args.rank_grid, int)]
        rank, rank_cv = choose_rank(
            matrix, source, target, train_index, grid, args.inner_frac, bandwidths,
            args.global_weight, args.cis_flank, min(args.metric_pairs, 500_000), args.seed,
        )
    else:
        rank = int(args.rank)
        if not 1 <= rank <= max_possible:
            raise SystemExit(f"ERROR: rank must be in [1,{max_possible}]")
        rank_cv = pd.DataFrame([{"RANK": rank, "N": np.nan, "PEARSON": np.nan, "SPEARMAN": np.nan,
                                 "RMSE": np.nan, "R2": np.nan, "TOP1_RECALL": np.nan,
                                 "TOP01_RECALL": np.nan}])
    atomic_dataframe(rank_cv, out / "rank.cv.tsv")

    train_meta = target.iloc[train_index].reset_index(drop=True)
    test_meta = target.iloc[test_index].reset_index(drop=True)
    log(f"Fit rank-{rank} alpha*beta graph on {len(train_index)} protein children")
    model = fit_components(matrix[:, train_index], source, train_meta, rank, args.cis_flank, args.seed)
    prediction, beta_test, col_test, cis_test = predict_components(
        model, source, train_meta, test_meta, bandwidths, args.global_weight, args.cis_flank
    )
    save_npy_atomic(out / "prediction.test.npy", prediction)

    row_pair, col_pair = sample_trans_pairs(
        source, test_meta, args.metric_pairs, args.cis_flank, args.seed + 101
    )
    observed = np.asarray(matrix[:, test_index])[row_pair, col_pair]
    predicted = prediction[row_pair, col_pair]
    finite_pair = np.isfinite(observed) & np.isfinite(predicted)
    if int(finite_pair.sum()) < 3:
        raise SystemExit("ERROR: fewer than three finite held-out trans pairs")
    if not np.all(finite_pair):
        log(f"Drop {int((~finite_pair).sum()):,} sampled held-out trans pair(s) with missing scores")
    row_pair = row_pair[finite_pair]
    col_pair = col_pair[finite_pair]
    observed = observed[finite_pair]
    predicted = predicted[finite_pair]
    metrics = metric_dict(observed, predicted)
    metric_table = pd.DataFrame([{"SCOPE": "cold_start_test_trans", **metrics}])
    atomic_dataframe(metric_table, out / "metrics.tsv")

    # Incremental held-out gain per factor; this is predictive evidence, not a
    # causal mediation test.
    incremental_rows = []
    baseline_pred = model["row_base"][row_pair].astype(float) + col_test[col_pair].astype(float)
    baseline_sse = float(np.sum((observed - baseline_pred) ** 2))
    current = baseline_pred.copy()
    previous_sse = baseline_sse
    total_var = float(np.sum((observed - np.mean(observed)) ** 2))
    for h in range(rank):
        current += model["alpha"][row_pair, h].astype(float) * beta_test[h, col_pair].astype(float)
        sse = float(np.sum((observed - current) ** 2))
        incremental_rows.append({
            "FACTOR": h + 1,
            "CV_DELTA_R2": (previous_sse - sse) / total_var if total_var > 0 else np.nan,
            "CV_CUMULATIVE_R2": 1.0 - sse / total_var if total_var > 0 else np.nan,
        })
        previous_sse = sse
    incremental = pd.DataFrame(incremental_rows)

    # Sparse gamma contains training edges not explained by cis or alpha*beta.
    residual = model.pop("residual_before_factor") - model["alpha"] @ model["beta"]
    center, sigma = robust_center_scale(residual.ravel())
    gamma_cut = abs(args.gamma_z) * sigma
    gamma_mask = np.abs(residual - center) >= gamma_cut
    gr, gc = np.where(gamma_mask)
    if len(gr) > args.gamma_max_edges:
        strength = np.abs(residual[gr, gc] - center)
        keep = np.argpartition(strength, -args.gamma_max_edges)[-args.gamma_max_edges:]
        gr, gc = gr[keep], gc[keep]
    gamma = pd.DataFrame({
        "SOURCE_INDEX": gr,
        "TARGET_LOCAL_INDEX": gc,
        "GAMMA": residual[gr, gc],
        "GAMMA_Z": (residual[gr, gc] - center) / sigma,
    })
    gamma["TARGET_INDEX"] = train_index[gamma["TARGET_LOCAL_INDEX"].to_numpy(int)]
    gamma = gamma.merge(source[["SOURCE_INDEX", "SOURCE", "CHR", "POS"]], on="SOURCE_INDEX", how="left")
    gamma = gamma.merge(target[["TARGET_INDEX", "PROTEIN", "GENE"]], on="TARGET_INDEX", how="left")
    gamma = gamma.sort_values("GAMMA_Z", key=lambda x: np.abs(x), ascending=False)
    atomic_dataframe(gamma, out / "gamma.edges.tsv.gz")

    # Export alpha and beta using explicit long tables.
    factor_names = np.arange(1, rank + 1, dtype=int)
    alpha_long = pd.DataFrame({
        "SOURCE_INDEX": np.repeat(source["SOURCE_INDEX"].to_numpy(int), rank),
        "FACTOR": np.tile(factor_names, len(source)),
        "ALPHA": model["alpha"].reshape(-1),
    }).merge(source[["SOURCE_INDEX", "SOURCE", "CHR", "POS"]], on="SOURCE_INDEX", how="left")
    atomic_dataframe(alpha_long, out / "alpha.source.tsv.gz")
    beta_all = np.full((rank, len(target)), np.nan, dtype=np.float32)
    beta_all[:, train_index] = model["beta"]
    beta_all[:, test_index] = beta_test
    beta_long = pd.DataFrame({
        "TARGET_INDEX": np.repeat(target["TARGET_INDEX"].to_numpy(int), rank),
        "FACTOR": np.tile(factor_names, len(target)),
        "BETA": beta_all.T.reshape(-1),
        "BETA_KIND": np.repeat(np.where(np.isin(np.arange(len(target)), train_index), "fitted", "position_predicted"), rank),
    }).merge(target[["TARGET_INDEX", "PROTEIN", "GENE", "CHR", "POS"]], on="TARGET_INDEX", how="left")
    atomic_dataframe(beta_long, out / "beta.target.tsv.gz")

    cis_table = target[["TARGET_INDEX", "PROTEIN", "GENE", "CHR", "POS"]].copy()
    cis_all = np.full(len(target), np.nan, dtype=float)
    cis_kind = np.full(len(target), "position_predicted", dtype=object)
    cis_all[train_index] = model["cis_alpha"]
    cis_all[test_index] = cis_test
    cis_kind[train_index] = "fitted"
    cis_table["CIS_ALPHA"] = cis_all
    cis_table["ALPHA_KIND"] = cis_kind
    atomic_dataframe(cis_table, out / "cis.alpha.tsv.gz")

    # A black-hole candidate is a stable predictive factor whose source loading
    # is localized, broad on the target side, and not already anchored by a
    # measured coding-gene position.
    src_chr = source["CHR"].map(normalize_chr).to_numpy(object)
    src_pos = pd.to_numeric(source["POS"], errors="coerce").to_numpy(float)
    measured_chr = target["CHR"].map(normalize_chr).to_numpy(object)
    measured_pos = pd.to_numeric(target["POS"], errors="coerce").to_numpy(float)
    black_rows = []
    top_rows = []
    factor_variance = model["singular"].astype(float) ** 2
    factor_variance /= max(float(np.sum(factor_variance)), np.finfo(float).eps)
    for h in range(rank):
        source_weight = model["alpha"][:, h].astype(float) ** 2
        source_weight /= max(float(source_weight.sum()), np.finfo(float).eps)
        target_weight = beta_all[h, :].astype(float) ** 2
        target_weight = np.nan_to_num(target_weight, nan=0.0)
        target_weight /= max(float(target_weight.sum()), np.finfo(float).eps)
        peak = int(np.argmax(source_weight))
        local = (src_chr == src_chr[peak]) & np.isfinite(src_pos) & (np.abs(src_pos - src_pos[peak]) <= 1_000_000)
        local_fraction = float(source_weight[local].sum())
        source_effective = float(1.0 / max(np.sum(source_weight ** 2), np.finfo(float).eps))
        target_effective = float(1.0 / max(np.sum(target_weight ** 2), np.finfo(float).eps))
        measured_peak = bool(np.any(
            (measured_chr == src_chr[peak]) & np.isfinite(measured_pos) &
            (np.abs(measured_pos - src_pos[peak]) <= args.cis_flank)
        ))
        cv_gain = float(incremental.loc[incremental["FACTOR"] == h + 1, "CV_DELTA_R2"].iloc[0])
        breadth = min(1.0, target_effective / max(25.0, 0.1 * len(target)))
        novelty = 0.25 if measured_peak else 1.0
        score = max(0.0, cv_gain) * local_fraction * breadth * novelty
        black_rows.append({
            "FACTOR": h + 1,
            "PEAK_SOURCE": source.iloc[peak]["SOURCE"],
            "PEAK_CHR": source.iloc[peak]["CHR"],
            "PEAK_POS": source.iloc[peak]["POS"],
            "MEASURED_PROTEIN_ANCHOR": measured_peak,
            "SOURCE_LOCAL_1MB_FRACTION": local_fraction,
            "SOURCE_EFFECTIVE_GENES": source_effective,
            "TARGET_EFFECTIVE_PROTEINS": target_effective,
            "FACTOR_VARIANCE_FRACTION": factor_variance[h],
            "CV_DELTA_R2": cv_gain,
            "BLACK_HOLE_SCORE": score,
            "INTERPRETATION": "putative latent proteome-regulatory factor; not a proven unknown protein",
        })
        top = np.argsort(source_weight)[-20:][::-1]
        for order, source_index in enumerate(top, start=1):
            top_rows.append({
                "FACTOR": h + 1,
                "RANK": order,
                "SOURCE_INDEX": source_index,
                "SOURCE": source.iloc[source_index]["SOURCE"],
                "CHR": source.iloc[source_index]["CHR"],
                "POS": source.iloc[source_index]["POS"],
                "ALPHA": model["alpha"][source_index, h],
                "WEIGHT": source_weight[source_index],
            })
    black = pd.DataFrame(black_rows).merge(incremental, on=["FACTOR", "CV_DELTA_R2"], how="left")
    black = black.sort_values("BLACK_HOLE_SCORE", ascending=False)
    atomic_dataframe(black, out / "black_holes.tsv")
    atomic_dataframe(pd.DataFrame(top_rows), out / "factor.top_sources.tsv.gz")

    model_path = out / "model.npz"
    tmp_model = out / f".model.{os.getpid()}.npz"
    np.savez_compressed(
        tmp_model,
        train_index=train_index,
        test_index=test_index,
        row_base=model["row_base"],
        col_base=model["col_base"],
        cis_alpha=model["cis_alpha"],
        alpha=model["alpha"],
        beta=model["beta"],
        beta_test=beta_test,
        col_test=col_test,
        cis_test=cis_test,
        singular=model["singular"],
        metric_row=row_pair,
        metric_col=col_pair,
    )
    os.replace(tmp_model, model_path)
    atomic_json(out / "fit.manifest.json", {
        "format": "phole-model-v1",
        "matrix_dir": str(matrix_dir.resolve()),
        "rank": rank,
        "n_train": len(train_index),
        "n_test": len(test_index),
        "train_fraction_requested": args.train_frac,
        "n_train_override": args.n_train,
        "train_fraction_actual": len(train_index) / len(target),
        "split": args.split,
        "seed": args.seed,
        "cis_flank": args.cis_flank,
        "gamma_z": args.gamma_z,
        "gamma_center": center,
        "gamma_sigma": sigma,
        "bandwidths": bandwidths,
        "global_weight": args.global_weight,
        "model_equation": "W = row + child + cis*C + alpha*beta + gamma + epsilon",
        "causal_claim": False,
    })
    (out / "fit.done").write_text("done\n", encoding="utf-8")
    log(f"PHOLE fit done: {model_path}")


if __name__ == "__main__":
    main()
