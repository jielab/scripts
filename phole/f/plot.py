#!/usr/bin/env python3
"""Create compact diagnostic figures for a PHOLE run."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import log, str2bool


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-dir", required=True)
    parser.add_argument("--model-dir", required=True)
    parser.add_argument("--validation-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--replace", default="FALSE")
    return parser.parse_args()


def main() -> None:
    args = arguments()
    matrix_dir = Path(args.matrix_dir)
    model_dir = Path(args.model_dir)
    validation_dir = Path(args.validation_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    replace = str2bool(args.replace)
    expected = [out / "prediction.png", out / "black_holes.png", out / "codenome_factors.png"]
    if all(x.exists() for x in expected) and not replace:
        log(f"PHOLE plots exist: {out}")
        return

    matrix = np.load(matrix_dir / "score.npy", mmap_mode="r")
    model = np.load(model_dir / "model.npz")
    prediction = np.load(model_dir / "prediction.test.npy", mmap_mode="r")
    test_index = model["test_index"].astype(int)
    rows = model["metric_row"].astype(int)
    cols = model["metric_col"].astype(int)
    obs = np.asarray(matrix[:, test_index])[rows, cols].astype(float)
    pred = prediction[rows, cols].astype(float)
    keep = np.isfinite(obs) & np.isfinite(pred)
    obs, pred = obs[keep], pred[keep]
    if len(obs) > 250_000:
        rng = np.random.default_rng(args.seed + 901)
        pick = rng.choice(len(obs), 250_000, replace=False)
        obs, pred = obs[pick], pred[pick]

    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    hb = ax.hexbin(obs, pred, gridsize=75, mincnt=1, bins="log", cmap="viridis")
    low = float(np.nanpercentile(np.r_[obs, pred], 0.5))
    high = float(np.nanpercentile(np.r_[obs, pred], 99.5))
    ax.plot([low, high], [low, high], linestyle="--", color="#d62728", linewidth=1.2)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_xlabel("Observed held-out trans score")
    ax.set_ylabel("Predicted held-out trans score")
    ax.set_title("PHOLE cold-start prediction")
    fig.colorbar(hb, ax=ax, label="log10 cell count")
    fig.tight_layout()
    fig.savefig(out / "prediction.png", dpi=240, facecolor="white")
    plt.close(fig)

    black_file = validation_dir / "black_holes.validated.tsv"
    if not black_file.exists():
        black_file = model_dir / "black_holes.tsv"
    black = pd.read_csv(black_file, sep="\t").sort_values("BLACK_HOLE_SCORE", ascending=False).head(20)
    labels = [f"F{int(f)}: {g} (chr{c})" for f, g, c in zip(black["FACTOR"], black["PEAK_SOURCE"], black["PEAK_CHR"])]
    supported = black.get("STATISTICALLY_SUPPORTED", pd.Series(False, index=black.index)).astype(bool)
    colors = np.where(supported, "#d62728", "#4c78a8")
    fig, ax = plt.subplots(figsize=(10.0, max(4.8, 0.36 * len(black) + 1.8)))
    y = np.arange(len(black))
    ax.barh(y, black["BLACK_HOLE_SCORE"], color=colors)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Exploratory black-hole score")
    ax.set_title("Localized latent proteome-effect factors")
    ax.grid(axis="x", color="0.9", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(out / "black_holes.png", dpi=240, facecolor="white")
    plt.close(fig)

    source = pd.read_csv(matrix_dir / "source.tsv.gz", sep="\t", dtype={"SOURCE": str, "CHR": str})
    alpha = model["alpha"].astype(float)
    top_factor = min(alpha.shape[1], 8)
    source_order = np.lexsort((pd.to_numeric(source["POS"], errors="coerce").fillna(0), source["CHR"].astype(str)))
    # Bin the long codenome to keep the image legible and lightweight.
    n_bin = min(500, len(source))
    groups = np.array_split(source_order, n_bin)
    binned = np.vstack([np.nanmean(alpha[index, :top_factor], axis=0) for index in groups]).T
    vmax = float(np.nanpercentile(np.abs(binned), 99))
    vmax = max(vmax, np.finfo(float).eps)
    fig, ax = plt.subplots(figsize=(13.0, 4.8))
    image = ax.imshow(binned, aspect="auto", interpolation="nearest", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_yticks(np.arange(top_factor), [f"Factor {x}" for x in range(1, top_factor + 1)])
    ax.set_xlabel("Codenome position (chromosome-ordered bins)")
    ax.set_title("Source loadings α across the coding genome")
    fig.colorbar(image, ax=ax, label="α loading")
    fig.tight_layout()
    fig.savefig(out / "codenome_factors.png", dpi=240, facecolor="white")
    plt.close(fig)
    log(f"PHOLE plots done: {out}")


if __name__ == "__main__":
    main()
