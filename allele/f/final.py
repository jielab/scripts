#!/usr/bin/env python3
"""Create compact final tables and figures for the four research questions."""
from __future__ import annotations
import argparse
import json
import os
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from parquet_io import batches


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=220, bbox_inches="tight")
    plt.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--work", required=True)
    p.add_argument("--results", required=True)
    p.add_argument("--out-dir", required=True)
    a = p.parse_args()
    for filename in [
        "geometry_candidates.parquet",
        "stage1_metrics.json",
        "stage2_metrics.json",
        "stage2_conditional_candidate.tsv",
        "stage2_conditional_all_candidates.tsv",
        "stage2_site_multiallelic.tsv",
        "stage2_triallelic_secondary.tsv",
        "stage2_tetra_sites.tsv",
    ]:
        path = os.path.join(a.results, filename)
        if not os.path.isfile(path):
            raise SystemExit(
                f"Missing result: {path}; finish geometry/stage1/stage2 before final"
            )
    os.makedirs(a.out_dir, exist_ok=True)
    for filename in [
        "Fig1_occupancy.png",
        "Fig2_observed_unobserved.png",
        "Fig3_triallelic_affinity.png",
        "Fig4_multiallelic_model.png",
    ]:
        path = os.path.join(a.out_dir, filename)
        if os.path.isfile(path):
            os.remove(path)
    tally = Counter()
    for part in batches(os.path.join(a.work, "sites.parquet"), ["occupancy"]):
        tally.update(part.occupancy.value_counts().to_dict())
    counts = pd.DataFrame(sorted(tally.items()), columns=["occupancy", "n_sites"])
    counts["allele_count"] = counts.occupancy + 1
    counts.to_csv(
        os.path.join(a.out_dir, "Table1_occupancy_counts.tsv"), sep="\t", index=False
    )
    plt.figure(figsize=(6, 4))
    plt.bar(counts.allele_count.astype(str), counts.n_sites)
    plt.xlabel("Alleles observed at locus")
    plt.ylabel("Number of loci")
    plt.title("Allelic occupancy")
    savefig(os.path.join(a.out_dir, "Fig1_occupancy.png"))
    gp = os.path.join(a.results, "geometry_candidates.parquet")
    if os.path.exists(gp):
        g = pd.read_parquet(gp)
        sub = g[g.occupancy.isin([1, 2])].copy()
        vals = [
            sub[(sub.occupancy == 2) & (sub.observed == 1)].effect_magnitude.dropna(),
            sub[(sub.occupancy == 2) & (sub.observed == 0)].effect_magnitude.dropna(),
        ]
        if all(len(x) for x in vals):
            plt.figure(figsize=(6, 4))
            plt.boxplot(
                vals, tick_labels=["observed ALT", "unobserved ALT"], showfliers=False
            )
            plt.ylabel("AlphaGenome effect magnitude")
            plt.title("Triallelic loci: observed vs unrealized substitution")
            savefig(os.path.join(a.out_dir, "Fig2_observed_unobserved.png"))
        tri = sub[(sub.occupancy == 2) & (sub.observed == 1) & (sub.is_primary == 0)]
        if len(tri):
            plt.figure(figsize=(6, 4))
            plt.hist(tri.affinity_to_primary.dropna(), bins=40)
            plt.axvline(0, ls="--")
            plt.xlabel("Affinity to primary ALT (>0 closer to primary ALT)")
            plt.ylabel("Secondary alleles")
            plt.title("Functional position of the third allele")
            savefig(os.path.join(a.out_dir, "Fig3_triallelic_affinity.png"))
    for src, dst in [
        ("stage2_conditional_candidate.tsv", "Table2_within_locus_retention.tsv"),
        ("stage2_conditional_all_candidates.tsv", "Table2b_all_candidates.tsv"),
        ("stage2_site_multiallelic.tsv", "Table3_multiallelic_model.tsv"),
        ("stage2_triallelic_secondary.tsv", "Table4_triallelic_geometry.tsv"),
        ("stage2_tetra_sites.tsv", "Table5_tetra_sites.tsv"),
    ]:
        q = os.path.join(a.results, src)
        if os.path.exists(q):
            pd.read_csv(q, sep="\t").to_csv(
                os.path.join(a.out_dir, dst), sep="\t", index=False
            )
    coefp = os.path.join(a.results, "stage2_site_multiallelic.tsv")
    if os.path.exists(coefp):
        c = pd.read_csv(coefp, sep="\t")
        if len(c):
            plt.figure(figsize=(7, 4))
            y = np.arange(len(c))
            plt.axvline(0, lw=1)
            plt.scatter(c.coef, y)
            plt.yticks(y, c.term)
            plt.xlabel("Log-odds coefficient per SD")
            plt.title("Predictors of multiallelic occupancy")
            savefig(os.path.join(a.out_dir, "Fig4_multiallelic_model.png"))
    with open(os.path.join(a.out_dir, "summary.txt"), "w") as f:
        f.write(
            f"Observed variable loci: {sum(tally.values()):,}\nTriallelic loci: {tally[2]:,}\nTetra-allelic loci: {tally[3]:,}\n"
        )
        f.write(
            "Interpretation: Stage 1 models relative mutation opportunity; Stage 2 asks whether AlphaGenome functional geometry adds retention information after that adjustment.\n"
        )
        f.write(
            "Table2 compares secondary vs unrealized ALT at triallelic sites, conditional on the observed primary anchor; Table2b compares all three ALTs using mutation score and magnitude only. Associations are exploratory, not causal retention estimates.\n"
        )
        for filename in ["stage1_metrics.json", "stage2_metrics.json"]:
            path = os.path.join(a.results, filename)
            if os.path.isfile(path):
                with open(path) as source:
                    metrics = json.load(source)
                with open(os.path.join(a.out_dir, filename), "w") as dest:
                    json.dump(metrics, dest, indent=2, allow_nan=False)
                f.write(filename + ": " + json.dumps(metrics, allow_nan=False) + "\n")
    print("Final outputs:", a.out_dir)


if __name__ == "__main__":
    main()
