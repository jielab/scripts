"""Publish same-person method comparisons from the aggregate TEST tables."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from core2 import atomic_table
from manuscript_stats import PRIMARY, SELECT

METHODS = {
    PRIMARY: "PRS-CSx", "disco": "DiscoDivas", "prsformer": "PRSformer",
    "abm": "GRID ABM", "abm.deep": "GRID ABM deep", "block.glmnet": "GRID glmnet",
    "block.lightgbm": "GRID LightGBM", "direct.blocknet": "GRID direct BlockNet",
    "direct.linear": "Global ridge", "direct.global_mlp":"Global MLP",
    "direct.regional_ridge":"Regional ridge",
}


def publish_comparison(out, tables, traits):
    full = tables["overall"].assign(ancestry="ALL")
    subgroup = tables["ancestry"]
    parts = [full, subgroup.loc[subgroup.ancestry.ne("ALL")], tables["screening"].assign(ancestry="ALL")]
    data = pd.concat(parts, ignore_index=True)
    data = data[data.method.isin(METHODS) & data.selection.isin(["all", SELECT])].copy()
    data["model"] = data.method.map(METHODS)
    keys = ["trait", "ancestry", "selection", "method"]
    if data.duplicated(keys).any():
        raise ValueError("Duplicate rows in held-out method comparison")
    counts = data.groupby(keys[:-1]).n_evaluable.nunique()
    if (counts > 1).any():
        raise ValueError("Method comparison sizes disagree within a TEST stratum")
    fields = keys + ["model","n_evaluable","n_cases","primary_loss","r2","partial_r2_vs_covariate","auroc","auprc"]
    data = data.reindex(columns=fields).sort_values(keys)
    atomic_table(data,out/"tables/table_method_comparison.tsv")
    paired = tables["paired"]
    subset = paired.loc[paired.method.isin(METHODS) & paired.reference.isin(METHODS)
                        & paired.comparison.eq("method_minus_reference")
                        & paired.metric.isin(["MSE","Brier"]) & paired.availability.eq("available")].copy()
    fig, axes = plt.subplots(1,len(traits),squeeze=False,figsize=(4.8*len(traits),4.5),layout="constrained")
    for ax,trait in zip(axes.flat,traits):
        rows = subset[(subset.trait==trait)&(subset.ancestry=="ALL")&(subset.selection=="all")&(subset.reference==PRIMARY)].dropna(subset=["difference"])
        ax.axvline(0, color="#84949b", linewidth=1)
        for i,row in enumerate(rows.itertuples()):
            if np.isfinite([row.ci_low,row.ci_high]).all():
                ax.plot([row.ci_low,row.ci_high],[i,i],color="#486773",lw=1.8)
            ax.scatter(row.difference,i,color="#075f61",zorder=3)
        ax.set_yticks(range(len(rows)),[METHODS[m] for m in rows.method])
        ax.invert_yaxis()
        ax.set_title(trait.upper() if trait!="height" else "Height")
        ax.set_xlabel(("Brier score" if trait=="t2dm" else "MSE")+" difference vs PRS-CSx")
        ax.spines[["top","right"]].set_visible(False)
    fig.savefig(out/"figures/fig5_method_comparison.png",dpi=300,bbox_inches="tight")
    plt.close(fig)
    status=tables["model_settings"][["trait","prsformer_status","disco_status"]]
    def records(frame):
        return frame.to_json(orient="records",double_precision=12).replace("<","\\u003c")
    # Embedded aggregate data only; the page also works from a local filesystem.
    html = Path(__file__).with_name("comparison.html").read_text(encoding="utf-8")
    payload = '{"results":'+records(data)+',"paired":'+records(subset)+',"status":'+records(status)+'}'
    (out/"comparison.html").write_text(html.replace("<!--DATA-->",payload),encoding="utf-8")
