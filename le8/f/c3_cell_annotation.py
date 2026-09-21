#!/usr/bin/env python3
"""Annotate frozen panels with externally published cell-enriched genes.

This is an annotation/enrichment analysis, not CIGMA, deconvolution or a claim
that a plasma protein was secreted by a particular cell. No atlas is fabricated.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from c3_cigma import bh


def annotate(universe, atlas, panels, outdir, prefix="cell"):
    u, a, p = (pd.read_csv(x, dtype=str) for x in (universe, atlas, panels))
    for frame, required in ((u, {"assay", "gene"}),
                            (a, {"gene", "cell_type", "source", "source_version"}),
                            (p, {"model", "feature"})):
        if not required <= set(frame.columns):
            raise ValueError(f"Missing columns: {required - set(frame.columns)}")
    if u.assay.isna().any() or u.assay.duplicated().any():
        raise ValueError("Universe must have unique nonmissing assays; leave unknown gene blank")
    if a[["gene", "cell_type", "source", "source_version"]].isna().any().any():
        raise ValueError("Atlas labels and provenance must be nonmissing")
    if len(a[["source", "source_version"]].drop_duplicates()) != 1:
        raise ValueError("Use one prespecified atlas/version per run; do not merge incompatible labels")
    a = a.drop_duplicates(["gene", "cell_type"])
    p = p[p.feature.notna()].drop_duplicates(["model", "feature"])
    if not set(p.feature) <= set(u.assay):
        raise ValueError("Every panel assay must be in the full measured assay universe")
    background = set(u.gene.dropna()) - {""}
    if not background:
        raise ValueError("No mapped background genes")
    annotated = p.merge(u, left_on="feature", right_on="assay", how="left").merge(a, on="gene", how="left")
    rows, coverage = [], []
    for model, group in p.groupby("model"):
        m = group.merge(u, left_on="feature", right_on="assay", how="left")
        genes = set(m.gene.dropna()) - {""}
        coverage.append(dict(model=model,assays=len(m),mapped_assays=int(m.gene.notna().sum()),
                             unique_genes=len(genes),cell_labelled_genes=len(genes & set(a.gene)),
                             background_genes=len(background)))
        for cell, table in a.groupby("cell_type"):
            cellgenes = set(table.gene) & background
            hit = len(genes & cellgenes)
            tab = [[hit, len(genes)-hit],
                   [len(cellgenes-genes), len(background-genes-cellgenes)]]
            odds, pv = fisher_exact(tab, alternative="greater") if genes else (np.nan, np.nan)
            rows.append(dict(model=model,cell_type=cell,hits=hit,panel_genes=len(genes),
                             atlas_genes_in_assay_background=len(cellgenes),background_genes=len(background),
                             odds_ratio=odds,p=pv))
    result = pd.DataFrame(rows)
    if len(result):
        result["FDR_all_panel_cell_tests"] = bh(result.p)
    outdir.mkdir(parents=True,exist_ok=True)
    annotated.to_csv(outdir/f"{prefix}.panel_annotation.csv",index=False)
    pd.DataFrame(coverage).to_csv(outdir/f"{prefix}.coverage.csv",index=False)
    result.to_csv(outdir/f"{prefix}.enrichment.csv",index=False)
    a[["source","source_version"]].drop_duplicates().assign(
        interpretation="External gene-expression annotation; putative cellular relevance, not plasma protein origin",
        denominator="Unique mapped genes across all measured assays; NPPB and NTPROBNP counted once for enrichment"
    ).to_csv(outdir/f"{prefix}.provenance.csv",index=False)
    if len(result):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        top = result.groupby("cell_type").hits.sum().nlargest(16).index
        z = result[result.cell_type.isin(top)].copy()
        z["evidence"] = -np.log10(z.FDR_all_panel_cell_tests.clip(lower=1e-300))
        # Average descriptive evidence across folds; never call this a combined P value.
        z["display_model"] = z.model.str.replace(r"\|[0-9]+$", "", regex=True)
        grid = z.pivot_table(index="cell_type", columns="display_model", values="evidence", aggfunc="mean")
        fig, ax = plt.subplots(figsize=(max(9, min(28, len(grid.columns)*.28)), 7))
        im = ax.imshow(grid.fillna(0), aspect="auto", cmap="Blues")
        ax.set_yticks(range(len(grid.index)), grid.index)
        ax.set_xticks(range(len(grid.columns)), grid.columns, rotation=90, fontsize=6)
        ax.set_title("External cell-expression annotation; descriptive mean across folds")
        fig.colorbar(im, ax=ax, label="Mean -log10 adjusted P (not a combined test)")
        fig.tight_layout()
        fig.savefig(outdir/f"{prefix}.enrichment.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--universe",type=Path,required=True)
    ap.add_argument("--atlas",type=Path,required=True)
    ap.add_argument("--panels",type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    ap.add_argument("--prefix",default="cell")
    arg = ap.parse_args()
    annotate(arg.universe,arg.atlas,arg.panels,arg.outdir,arg.prefix)
