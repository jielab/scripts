#!/usr/bin/env python3
"""Permutation validation for PHOLE cold-start prediction and latent factors."""

from __future__ import annotations

import argparse
import os
import time
import gpu
from pathlib import Path

import numpy as np
import pandas as pd


# 🚩 Project utilities
import importlib.util as _module_util
import sys as _module_sys
from pathlib import Path as _ModulePath

_common_name = "phole_common"
if _common_name not in _module_sys.modules:
	_spec = _module_util.spec_from_file_location(_common_name, _ModulePath(__file__).with_name("0.common.py"))
	_common = _module_util.module_from_spec(_spec)
	_module_sys.modules[_common_name] = _common
	_spec.loader.exec_module(_common)
else:
	_common = _module_sys.modules[_common_name]
atomic_dataframe = _common.atomic_dataframe
bh_fdr = _common.bh_fdr
log = _common.log
metric_dict = _common.metric_dict
position_kernel = _common.position_kernel
read_json = _common.read_json
str2bool = _common.str2bool


def arguments() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--matrix-dir", required=True)
	parser.add_argument("--model-dir", required=True)
	parser.add_argument("--out", required=True)
	parser.add_argument("--perm", type=int, default=100)
	parser.add_argument("--seed", type=int, default=20260828)
	parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default=os.getenv("PHOLE_DEVICE", "auto"))
	parser.add_argument("--gpu-batch", type=int, default=250_000)
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
	started = time.perf_counter()
	device_info = gpu.preflight(args.device)
	args.device = device_info["backend"]
	log(f"PHOLE permutation backend={args.device}")
	matrix_dir = Path(args.matrix_dir)
	model_dir = Path(args.model_dir)
	out = Path(args.out)
	out.mkdir(parents=True, exist_ok=True)
	replace = str2bool(args.replace)
	done = out / "validate.done"
	if done.exists() and (out / "validation.tsv").exists() and not replace:
		from json import loads
		timing = out / "validation.timing.json"
		saved = loads(timing.read_text()) if timing.exists() else {"backend": "cpu"}
		if saved["backend"] != args.device:
			raise SystemExit("ERROR: existing validation uses another backend; use a separate --out or --replace TRUE")
		log(f"PHOLE validation exists: {out / 'validation.tsv'}")
		return
	if args.perm < 1:
		raise SystemExit("ERROR: --perm must be positive")

	matrix = np.load(matrix_dir / "score.npy", mmap_mode="r")
	target = pd.read_csv(matrix_dir / "target.tsv.gz", sep="\t", dtype={"PROTEIN": str, "GENE": str, "CHR": str})
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

	predictor = gpu.PairPredictor(row_base, alpha, beta, col_base, args.gpu_batch) if args.device == "cuda" else None
	prediction_seconds = 0.0
	kernel_seconds = 0.0
	def evaluate(kernel, row_indices, col_indices, y):
		nonlocal prediction_seconds
		t0 = time.perf_counter()
		if predictor is not None:
			answer = predictor.evaluate(kernel, row_indices, col_indices, y)
		else:
			bq = beta.astype(float) @ kernel.T
			cq = kernel @ col_base.astype(float)
			pred = prediction_on_pairs(row_base, alpha, bq, cq, row_indices, col_indices, rank)
			sse = np.array([np.sum((y-prediction_on_pairs(row_base, alpha, bq, cq, row_indices, col_indices, h))**2) for h in range(rank+1)])
			answer = pred, sse
		prediction_seconds += time.perf_counter()-t0
		return answer

	true_kernel = position_kernel(train_meta, test_meta, bandwidths, global_weight)
	true_pred, _ = evaluate(true_kernel, rows, cols, observed)
	finite_pair = np.isfinite(observed) & np.isfinite(true_pred)
	if int(finite_pair.sum()) < 3:
		raise SystemExit("ERROR: fewer than three finite saved held-out trans pairs")
	rows, cols, observed = rows[finite_pair], cols[finite_pair], observed[finite_pair]
	true_pred, true_sse_array = evaluate(true_kernel, rows, cols, observed)
	true_metrics = metric_dict(observed, true_pred)
	total_var = float(np.sum((observed - np.mean(observed)) ** 2))
	true_sse = -np.diff(true_sse_array)/total_var if total_var > 0 else np.full(rank, np.nan)

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
		t0 = time.perf_counter()
		kernel = position_kernel(train_meta, permuted_meta, bandwidths, global_weight)
		kernel_seconds += time.perf_counter()-t0
		pred, sse = evaluate(kernel, rows, cols, observed)
		metrics = metric_dict(observed, pred)
		null_rows.append({"PERM": b + 1, **metrics})
		factor_null[b] = -np.diff(sse)/total_var if total_var > 0 else np.nan
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
		validation_rows.append(
			{
				"METRIC": metric,
				"OBSERVED": observed_metric,
				"NULL_MEAN": float(np.mean(null_metric)),
				"NULL_SD": float(np.std(null_metric, ddof=1)) if len(null_metric) > 1 else np.nan,
				"EMPIRICAL_P": p,
				"ALTERNATIVE": direction,
			}
		)
	validation = pd.DataFrame(validation_rows)
	atomic_dataframe(validation, out / "validation.tsv")
	atomic_dataframe(null, out / "permutation.metrics.tsv.gz")

	black = pd.read_csv(model_dir / "black_holes.tsv", sep="\t")
	p_factor = []
	for h in range(rank):
		null_h = factor_null[:, h]
		null_h = null_h[np.isfinite(null_h)]
		p_factor.append((1.0 + np.sum(null_h >= true_sse[h])) / (1.0 + len(null_h)))
	factor_validation = pd.DataFrame(
		{
			"FACTOR": np.arange(1, rank + 1),
			"CV_DELTA_R2_RECALCULATED": true_sse,
			"FACTOR_EMPIRICAL_P": p_factor,
		}
	)
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
	_common.atomic_json(out / "validation.timing.json", dict(module="phole-permutation",
		**{**device_info, "gpu_peak_memory_mb": gpu.memory_mb(args.device)},
		wall_seconds=time.perf_counter()-started, prediction_seconds=prediction_seconds,
		position_kernel_seconds=kernel_seconds, permutations=args.perm, pairs=len(rows), rank=rank, seed=args.seed))
	done.write_text("done\n", encoding="utf-8")
	log(f"PHOLE validation done: {out / 'validation.tsv'}")


if __name__ == "__main__":
	main()
