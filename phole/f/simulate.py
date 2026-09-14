#!/usr/bin/env python3
"""Generate a small MAGMA-shaped synthetic proteome for PHOLE self-testing."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

from common import log


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("--n-source", type=int, default=600)
    parser.add_argument("--n-protein", type=int, default=72)
    parser.add_argument("--rank", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260828)
    return parser.parse_args()


def main() -> None:
    args = arguments()
    out = Path(args.out)
    root = out / "gwas"
    root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    n_chr = 6
    source_chr = np.repeat(np.arange(1, n_chr + 1), int(np.ceil(args.n_source / n_chr)))[: args.n_source]
    source_pos = np.empty(args.n_source, dtype=int)
    for chrom in range(1, n_chr + 1):
        index = np.flatnonzero(source_chr == chrom)
        source_pos[index] = np.sort(rng.integers(1_000_000, 180_000_000, size=len(index)))
    source_id = np.asarray([f"GENE_{x:06d}" for x in range(args.n_source)])
    target_source = np.sort(rng.choice(args.n_source, size=args.n_protein, replace=False))
    target_chr = source_chr[target_source]
    target_pos = source_pos[target_source]
    proteins = np.asarray([f"PROT_{x:03d}" for x in range(args.n_protein)])

    alpha = np.zeros((args.n_source, args.rank), dtype=float)
    peaks = rng.choice(np.setdiff1d(np.arange(args.n_source), target_source), size=args.rank, replace=False)
    for h, peak in enumerate(peaks):
        same = source_chr == source_chr[peak]
        distance = np.abs(source_pos - source_pos[peak])
        alpha[:, h] = rng.normal(0, 0.02, args.n_source)
        alpha[same, h] += (1.8 + 0.25 * h) * np.exp(-distance[same] / 2_000_000)

    beta = np.zeros((args.rank, args.n_protein), dtype=float)
    for h in range(args.rank):
        phase = 0.4 + h * 0.9
        beta[h, :] = np.sin(target_pos / (14_000_000 + h * 5_000_000) + phase)
        beta[h, :] += 0.35 * np.cos(target_chr * (h + 1))
        beta[h, :] += rng.normal(0, 0.08, args.n_protein)
    row_base = rng.normal(0, 0.10, args.n_source)
    w = row_base[:, None] + alpha @ beta + rng.normal(0, 0.45, (args.n_source, args.n_protein))
    for j, mother in enumerate(target_source):
        w[mother, j] += rng.normal(4.5, 0.35)
    # A few sparse direct gamma edges.
    for _ in range(max(10, args.n_protein // 2)):
        w[rng.integers(args.n_source), rng.integers(args.n_protein)] += rng.choice([-1, 1]) * rng.uniform(3, 5)

    bed = pd.DataFrame({
        "CHR": target_chr,
        "START": target_pos - 2_000,
        "END": target_pos + 2_000,
        "PROTEIN": proteins,
        "GENE": [f"GENE_{x:06d}" for x in target_source],
    })
    bed.to_csv(out / "protein.b38.bed", sep="\t", header=False, index=False)
    gene_loc = out / "synthetic.NCBI38.gene.loc"
    pd.DataFrame({
        "GENE": source_id,
        "CHR": source_chr,
        "START": source_pos - 5_000,
        "STOP": source_pos + 5_000,
        "STRAND": "+",
        "SYMBOL": source_id,
    }).to_csv(gene_loc, sep="\t", header=False, index=False)

    for j, protein in enumerate(proteins):
        # Match the current gwas_post.sh layout:
        # <project>/<category>/<GWAS>/magma/<GWAS>.genes.out
        magma_dir = root / "common" / protein / "magma"
        magma_dir.mkdir(parents=True, exist_ok=True)
        z = w[:, j]
        p = norm.sf(z)
        dat = pd.DataFrame({
            "GENE": source_id,
            "CHR": source_chr,
            "START": source_pos - 5_000,
            "STOP": source_pos + 5_000,
            "NSNPS": rng.integers(2, 30, args.n_source),
            "NPARAM": rng.integers(1, 10, args.n_source),
            "N": 50_000,
            "ZSTAT": z,
            "P": np.clip(p, np.finfo(float).tiny, 1.0),
        })
        dat.to_csv(magma_dir / f"{protein}.genes.out", sep="\t", index=False, float_format="%.8g")
        (magma_dir / "magma.done").write_text("2026-08-28 00:00:00\n", encoding="utf-8")
        pd.DataFrame({
            "key": ["gwas", "grch", "gene_loc"],
            "value": [protein, "38", str(gene_loc)],
        }).to_csv(magma_dir / "magma.meta.tsv", sep="\t", index=False)
    pd.DataFrame({
        "FACTOR": np.arange(1, args.rank + 1),
        "TRUE_PEAK_SOURCE": source_id[peaks],
        "TRUE_CHR": source_chr[peaks],
        "TRUE_POS": source_pos[peaks],
    }).to_csv(out / "truth.tsv", sep="\t", index=False)
    log(f"Synthetic MAGMA proteome done: {root}")


if __name__ == "__main__":
    main()
