#!/usr/bin/env python3
"""Permutation validation for PHOLE cold-start prediction and latent factors."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from common import atomic_dataframe, bh_fdr, log, metric_dict, position_kernel, read_json, str2bool


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-dir", required=True)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--perm", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--replace", default="FALSE")
    return parser.parse_args()


def prediction_on_pairs(
    row_base: np.ndarray,
    alpha: np.ndarray,
    beta_query: np.ndarray,
    col_query: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    rank: int,
) -> np.ndarray:
    pred = row_base[rows].astype(float) + col_query[cols].astype(float)
    if rank:
        pred += np.sum(alpha[rows, :rank].astype(float) * beta_query[:rank, cols].T.astype(float), axis=1)
    return pred


def main() -> None:
    args = arguments()
    matrix_dir = Path(args.matrix_dir)
    model_dir = Path(args.model_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    replace = str2bool(args.replace)
    done = out / "validate.done"
    if done.exists() and (out / "validation.tsv").exists() and not replace:
        log(f"PHOLE validation exists: {out / 'validation.tsv'}")
        return
    if args.perm < 1:
        raise SystemExit("ERROR: --perm must be positive")

    matrix = np.load(matrix_dir / "score.npy", mmap_mode="r")
    target = pd.read_csv(
        matrix_dir / "target.tsv.gz", sep="\t", dtype={"PROTEIN": str, "GENE": str, "CHR": str}
    )
    model = np.load(model_dir / "model.npz")
    manifest = read_json(model_dir / "fit.manifest.json")
    train_index = model["train_index"].astype(int)
    test_index = model["test_index"].astype(int)
    train_meta = target.iloc[train_index].reset_index(drop=True)
    test_meta = target.iloc[test_index].reset_index(drop=True)
    rows = model["metric_row"].astype(int)
    cols = model["metric_col"].astype(int)
    observed = np.asarray(matrix[:, test_index])[rows, cols].astype(float)
    alpha = model["alpha"]
    beta = model["beta"]
    row_base = model["row_base"]
    col_base = model["col_base"]
    rank = alpha.shape[1]
    bandwidths = [float(x) for x in manifest["bandwidths"]]
    global_weight = float(manifest["global_weight"])

    true_kernel = position_kernel(train_meta, test_meta, bandwidths, global_weight)
    true_beta = beta.astype(float) @ true_kernel.T
    true_col = true_kernel @ col_base.astype(float)
    true_pred = prediction_on_pairs(row_base, alpha, true_beta, true_col, rows, cols, rank)
    finite_pair = np.isfinite(observed) & np.isfinite(true_pred)
    if int(finite_pair.sum()) < 3:
        raise SystemExit("ERROR: fewer than three finite saved held-out trans pairs")
    if not np.all(finite_pair):
        log(f"Drop {int((~finite_pair).sum()):,} saved held-out trans pair(s) with missing scores")
    rows = rows[finite_pair]
    cols = cols[finite_pair]
    observed = observed[finite_pair]
    true_pred = true_pred[finite_pair]
    true_metrics = metric_dict(observed, true_pred)

    true_sse = []
    previous = prediction_on_pairs(row_base, alpha, true_beta, true_col, rows, cols, 0)
    previous_sse = float(np.sum((observed - previous) ** 2))
    total_var = float(np.sum((observed - np.mean(observed)) ** 2))
    for h in range(1, rank + 1):
        pred = prediction_on_pairs(row_base, alpha, true_beta, true_col, rows, cols, h)
        sse = float(np.sum((observed - pred) ** 2))
        true_sse.append((previous_sse - sse) / total_var if total_var > 0 else np.nan)
        previous_sse = sse

    rng = np.random.default_rng(args.seed + 701)
    null_rows = []
    factor_null = np.full((args.perm, rank), np.nan, dtype=float)
    log(f"Position-label permutation validation: B={args.perm} sampled trans edges={len(rows):,}")
    for b in range(args.perm):
        # Shuffle the coding-gene positions assigned to held-out proteins while
        # leaving their observed effect profiles untouched.
        permuted_meta = test_meta.copy()
        order = rng.permutation(len(test_meta))
        for column in ["CHR", "START", "END", "POS", "GENE"]:
            if column in permuted_meta:
                permuted_meta[column] = test_meta.iloc[order][column].to_numpy()
        kernel = position_kernel(train_meta, permuted_meta, bandwidths, global_weight)
        beta_query = beta.astype(float) @ kernel.T
        col_query = kernel @ col_base.astype(float)
        pred = prediction_on_pairs(row_base, alpha, beta_query, col_query, rows, cols, rank)
        metrics = metric_dict(observed, pred)
        null_rows.append({"PERM": b + 1, **metrics})

        previous = prediction_on_pairs(row_base, alpha, beta_query, col_query, rows, cols, 0)
        previous_sse = float(np.sum((observed - previous) ** 2))
        for h in range(1, rank + 1):
            current = prediction_on_pairs(row_base, alpha, beta_query, col_query, rows, cols, h)
            sse = float(np.sum((observed - current) ** 2))
            factor_null[b, h - 1] = (previous_sse - sse) / total_var if total_var > 0 else np.nan
            previous_sse = sse
        if (b + 1) % 10 == 0 or b + 1 == args.perm:
            log(f"permutations: {b + 1}/{args.perm}")

    null = pd.DataFrame(null_rows)
    metric_direction = {
        "PEARSON": "greater",
        "SPEARMAN": "greater",
        "R2": "greater",
        "TOP1_RECALL": "greater",
        "TOP01_RECALL": "greater",
        "RMSE": "less",
    }
    validation_rows = []
    for metric, direction in metric_direction.items():
        observed_metric = float(true_metrics[metric])
        null_metric = pd.to_numeric(null[metric], errors="coerce").to_numpy(float)
        null_metric = null_metric[np.isfinite(null_metric)]
        if direction == "greater":
            extreme = int(np.sum(null_metric >= observed_metric))
        else:
            extreme = int(np.sum(null_metric <= observed_metric))
        p = (extreme + 1.0) / (len(null_metric) + 1.0)
        validation_rows.append({
            "METRIC": metric,
            "OBSERVED": observed_metric,
            "NULL_MEAN": float(np.mean(null_metric)),
            "NULL_SD": float(np.std(null_metric, ddof=1)) if len(null_metric) > 1 else np.nan,
            "EMPIRICAL_P": p,
            "ALTERNATIVE": direction,
        })
    validation = pd.DataFrame(validation_rows)
    atomic_dataframe(validation, out / "validation.tsv")
    atomic_dataframe(null, out / "permutation.metrics.tsv.gz")

    black = pd.read_csv(model_dir / "black_holes.tsv", sep="\t")
    p_factor = []
    for h in range(rank):
        null_h = factor_null[:, h]
        null_h = null_h[np.isfinite(null_h)]
        p_factor.append((1.0 + np.sum(null_h >= true_sse[h])) / (1.0 + len(null_h)))
    factor_validation = pd.DataFrame({
        "FACTOR": np.arange(1, rank + 1),
        "CV_DELTA_R2_RECALCULATED": true_sse,
        "FACTOR_EMPIRICAL_P": p_factor,
    })
    factor_validation["FACTOR_FDR"] = bh_fdr(factor_validation["FACTOR_EMPIRICAL_P"])
    black = black.merge(factor_validation, on="FACTOR", how="left")
    black["STATISTICALLY_SUPPORTED"] = (
        (black["FACTOR_FDR"] <= 0.05)
        & (black["CV_DELTA_R2_RECALCULATED"] > 0)
        & (black["SOURCE_LOCAL_1MB_FRACTION"] >= 0.25)
    )
    black = black.sort_values(["STATISTICALLY_SUPPORTED", "BLACK_HOLE_SCORE"], ascending=[False, False])
    atomic_dataframe(black, out / "black_holes.validated.tsv")
    factor_null_long = pd.DataFrame(
        {
            "PERM": np.repeat(np.arange(1, args.perm + 1), rank),
            "FACTOR": np.tile(np.arange(1, rank + 1), args.perm),
            "NULL_DELTA_R2": factor_null.reshape(-1),
        }
    )
    atomic_dataframe(factor_null_long, out / "factor.permutation.tsv.gz")
    done.write_text("done\n", encoding="utf-8")
    log(f"PHOLE validation done: {out / 'validation.tsv'}")


if __name__ == "__main__":
    main()
