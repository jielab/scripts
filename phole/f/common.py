#!/usr/bin/env python3
"""Shared utilities for the PHOLE proteome-effect graph pipeline."""

from __future__ import annotations

import gzip
import json
import math
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


def log(message: str) -> None:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {message}", file=sys.stderr, flush=True)


def str2bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().upper()
    if text not in {"TRUE", "FALSE"}:
        raise ValueError(f"expected TRUE or FALSE, got: {value}")
    return text == "TRUE"


def normalize_chr(value: object) -> str:
    text = str(value).strip().replace("chr", "").replace("CHR", "")
    if text.endswith(".0"):
        text = text[:-2]
    text = text.upper()
    aliases = {"23": "X", "24": "Y", "25": "MT", "M": "MT"}
    return aliases.get(text, text)


def chromosome_order(value: object) -> int:
    chrom = normalize_chr(value)
    if chrom.isdigit():
        return int(chrom)
    return {"X": 23, "Y": 24, "MT": 25}.get(chrom, 10_000)


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def atomic_dataframe(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".gz" if path.name.endswith(".gz") else ""
    fd, tmp = tempfile.mkstemp(prefix=f".{path.stem}.", suffix=suffix, dir=path.parent)
    os.close(fd)
    try:
        df.to_csv(tmp, sep="\t", index=False, compression="gzip" if suffix else None)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_csv_numbers(value: str, cast=float) -> list:
    out = [cast(x.strip()) for x in str(value).split(",") if x.strip()]
    if not out:
        raise ValueError(f"empty numeric list: {value}")
    return out


def robust_center_scale(values: np.ndarray, mask: np.ndarray | None = None) -> tuple[float, float]:
    x = np.asarray(values, dtype=float)
    keep = np.isfinite(x)
    if mask is not None:
        keep &= np.asarray(mask, dtype=bool)
    x = x[keep]
    if x.size == 0:
        return 0.0, 1.0
    center = float(np.median(x))
    mad = float(np.median(np.abs(x - center))) * 1.4826
    if not np.isfinite(mad) or mad <= 1e-8:
        mad = float(np.std(x))
    if not np.isfinite(mad) or mad <= 1e-8:
        mad = 1.0
    return center, mad


def cis_rows(
    source_chr: np.ndarray,
    source_pos: np.ndarray,
    target_chr: object,
    target_pos: object,
    flank: int,
) -> np.ndarray:
    chrom = normalize_chr(target_chr)
    try:
        pos = float(target_pos)
    except (TypeError, ValueError):
        return np.zeros(len(source_chr), dtype=bool)
    if not np.isfinite(pos):
        return np.zeros(len(source_chr), dtype=bool)
    return (source_chr == chrom) & np.isfinite(source_pos) & (np.abs(source_pos - pos) <= flank)


def position_kernel(
    train_meta: pd.DataFrame,
    query_meta: pd.DataFrame,
    bandwidths: Sequence[float],
    global_weight: float,
) -> np.ndarray:
    """Query-by-train kernel using only coding-gene positions.

    Cross-chromosome pairs receive only ``global_weight``. Same-chromosome
    pairs also receive the mean of several exponential distance kernels.
    """
    tr_chr = train_meta["CHR"].map(normalize_chr).to_numpy(object)
    qu_chr = query_meta["CHR"].map(normalize_chr).to_numpy(object)
    tr_pos = pd.to_numeric(train_meta["POS"], errors="coerce").to_numpy(float)
    qu_pos = pd.to_numeric(query_meta["POS"], errors="coerce").to_numpy(float)
    same = qu_chr[:, None] == tr_chr[None, :]
    valid = np.isfinite(qu_pos)[:, None] & np.isfinite(tr_pos)[None, :]
    dist = np.abs(qu_pos[:, None] - tr_pos[None, :])
    kernel = np.full(dist.shape, float(global_weight), dtype=np.float64)
    local = np.zeros(dist.shape, dtype=np.float64)
    for bandwidth in bandwidths:
        local += np.exp(-dist / float(bandwidth))
    local /= max(len(bandwidths), 1)
    kernel += np.where(same & valid, local, 0.0)
    row_sum = kernel.sum(axis=1, keepdims=True)
    bad = ~np.isfinite(row_sum[:, 0]) | (row_sum[:, 0] <= 0)
    if np.any(bad):
        kernel[bad, :] = 1.0
        row_sum = kernel.sum(axis=1, keepdims=True)
    return kernel / row_sum


def randomized_svd(
    matrix: np.ndarray,
    rank: int,
    seed: int,
    n_iter: int = 2,
    oversample: int = 8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Dependency-light randomized truncated SVD."""
    a = np.asarray(matrix, dtype=np.float64, order="C")
    m, n = a.shape
    k = max(1, min(int(rank), m, n))
    ell = min(k + int(oversample), m, n)
    rng = np.random.default_rng(seed)
    omega = rng.normal(size=(n, ell))
    q, _ = np.linalg.qr(a @ omega, mode="reduced")
    for _ in range(max(0, int(n_iter))):
        z, _ = np.linalg.qr(a.T @ q, mode="reduced")
        q, _ = np.linalg.qr(a @ z, mode="reduced")
    small = q.T @ a
    ub, singular, vt = np.linalg.svd(small, full_matrices=False)
    u = q @ ub[:, :k]
    return u[:, :k], singular[:k], vt[:k, :]


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    a = np.asarray(x, dtype=float).ravel()
    b = np.asarray(y, dtype=float).ravel()
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() < 3:
        return float("nan")
    a = a[keep]
    b = b[keep]
    if np.std(a) <= 0 or np.std(b) <= 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def rankdata(values: np.ndarray) -> np.ndarray:
    """Average ranks, equivalent to scipy.stats.rankdata(method='average')."""
    x = np.asarray(values)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    sorted_x = x[order]
    starts = np.r_[0, np.flatnonzero(sorted_x[1:] != sorted_x[:-1]) + 1]
    ends = np.r_[starts[1:], len(x)]
    for start, end in zip(starts, ends):
        ranks[order[start:end]] = (start + end - 1) / 2.0 + 1.0
    return ranks


def metric_dict(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    obs = np.asarray(observed, dtype=float).ravel()
    pred = np.asarray(predicted, dtype=float).ravel()
    keep = np.isfinite(obs) & np.isfinite(pred)
    obs = obs[keep]
    pred = pred[keep]
    if obs.size < 3:
        return {"N": int(obs.size), "PEARSON": np.nan, "SPEARMAN": np.nan,
                "RMSE": np.nan, "R2": np.nan, "TOP1_RECALL": np.nan,
                "TOP01_RECALL": np.nan}
    pearson = correlation(obs, pred)
    spearman = correlation(rankdata(obs), rankdata(pred))
    rmse = float(np.sqrt(np.mean((obs - pred) ** 2)))
    denom = float(np.sum((obs - np.mean(obs)) ** 2))
    r2 = float(1.0 - np.sum((obs - pred) ** 2) / denom) if denom > 0 else np.nan

    def top_recall(frac: float) -> float:
        count = max(1, int(math.ceil(frac * obs.size)))
        truth = np.argpartition(obs, -count)[-count:]
        called = np.argpartition(pred, -count)[-count:]
        return float(np.intersect1d(truth, called, assume_unique=False).size / count)

    return {
        "N": int(obs.size),
        "PEARSON": pearson,
        "SPEARMAN": spearman,
        "RMSE": rmse,
        "R2": r2,
        "TOP1_RECALL": top_recall(0.01),
        "TOP01_RECALL": top_recall(0.001),
    }


def bh_fdr(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvalues), dtype=float)
    out = np.full(p.shape, np.nan)
    keep = np.isfinite(p)
    if not keep.any():
        return out
    vals = p[keep]
    order = np.argsort(vals)
    ranked = vals[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    out[np.flatnonzero(keep)] = q[inv]
    return out


def sample_trans_pairs(
    source: pd.DataFrame,
    target: pd.DataFrame,
    n_pairs: int,
    cis_flank: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    m, n = len(source), len(target)
    wanted = min(int(n_pairs), m * n)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    got = 0
    src_chr = source["CHR"].map(normalize_chr).to_numpy(object)
    src_pos = pd.to_numeric(source["POS"], errors="coerce").to_numpy(float)
    tgt_chr = target["CHR"].map(normalize_chr).to_numpy(object)
    tgt_pos = pd.to_numeric(target["POS"], errors="coerce").to_numpy(float)
    while got < wanted:
        batch = min(max(10_000, wanted - got), 1_000_000)
        r = rng.integers(0, m, size=batch)
        c = rng.integers(0, n, size=batch)
        is_cis = (src_chr[r] == tgt_chr[c]) & np.isfinite(src_pos[r]) & np.isfinite(tgt_pos[c]) & (
            np.abs(src_pos[r] - tgt_pos[c]) <= cis_flank
        )
        r, c = r[~is_cis], c[~is_cis]
        take = min(len(r), wanted - got)
        rows.append(r[:take])
        cols.append(c[:take])
        got += take
    return np.concatenate(rows), np.concatenate(cols)

