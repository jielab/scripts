#!/usr/bin/env python3
"""Non-destructive, aggregate-only audit and compact figures for any 5C outcome.

This reads published tables, never selects features for prediction. Missing
modules and null results remain visible. Existing figures/tables are untouched.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm


def bh(p, family_size=None):
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    out = np.full(len(p), np.nan)
    if not ok.any():
        return out
    n = max(len(p), int(family_size or 0))
    ids = np.flatnonzero(ok)
    order = np.argsort(p[ok], kind="stable")
    adjusted = p[ids[order]] * n / np.arange(1, len(ids) + 1)
    out[ids[order]] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return out


def read(path):
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def aggregate(root):
    inventory, summary, concordance, windows, trajectories, pillars, contrasts = [], [], [], [], [], [], []
    for layer, prefix in (("prot", "pwas"), ("met", "mwas")):
        for module in ("c1_correlate", "c2_cause", "c3_coloc", "c4_connect", "c4_focus", "c5_consolidate"):
            files = sorted((root / layer / module).glob("*.csv"))
            if not files:
                inventory.append(dict(layer=layer, module=module, file="", rows=0, status="unavailable"))
            for f in files:
                d = read(f)
                inventory.append(dict(layer=layer, module=module, file=str(f.relative_to(root)),
                                      rows=len(d), status="available" if len(d) else "empty",
                                      sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
        c1 = root / layer / "c1_correlate"
        incident = read(c1 / f"{prefix}_incident_adj2.csv")
        pgs = read(c1 / f"{prefix}_pgs_incident_full_genetic.csv")
        for kind, d in (("measured", incident), ("biomarker_PGS", pgs)):
            if {"term", "p.value", "FDR"} <= set(d):
                summary.append(dict(layer=layer, analysis=kind, tested=len(d),
                                    finite_p=int(d["p.value"].notna().sum()),
                                    FDR05=int((d.FDR < .05).sum()),
                                    N_min=d.N_total.min(), N_max=d.N_total.max(),
                                    events_min=d.N_event.min(), events_max=d.N_event.max()))
        if len(incident) and len(pgs):
            cols = ["term", "beta", "std.error", "FDR", "N_total", "N_event"]
            z = incident[cols].merge(pgs[cols], on="term", how="outer", suffixes=("_measured", "_PGS"))
            z["layer"] = layer
            z["both_supported"] = (z.FDR_measured < .05) & (z.FDR_PGS < .05)
            z["same_direction"] = np.sign(z.beta_measured) == np.sign(z.beta_PGS)
            z["interpretation"] = "Parallel associations; different SD units/covariates/cohorts; not MR or a causal test"
            concordance.append(z)
        lm = read(c1 / f"{prefix}_incident_landmark_adj2.csv")
        if len(lm) and len(incident):
            # Legacy scans selected about 500 proteins using the same outcomes.
            # Retain original FDR and add the full assay x landmark family.
            lm["FDR_assay_landmark_family"] = bh(lm["p.value"], len(incident) * lm.landmark_years.nunique())
            lm["subset_selected_before_landmark"] = lm.term.nunique() < len(incident)
            lm["layer"] = layer
            trajectories.append(lm)
        rw = read(c1 / f"{prefix}_diagnosis_window_riskset_adj2.csv")
        if len(rw):
            z = rw.groupby(["side", "window_lo", "window_hi"], dropna=False).agg(
                events_min=("N_event", "min"), events_max=("N_event", "max"),
                N_min=("N_total", "min"), N_max=("N_total", "max")).reset_index()
            z["layer"] = layer
            z["warning"] = np.where((z.side == "Post-baseline incident") & (z.events_max == 0),
                                     "No observed events; verify registry coverage before interpreting lead time", "")
            windows.append(z)
        p = read(root / layer / "c4_focus/c4.focus.pillar_counts.csv")
        if len(p):
            p["layer"] = layer
            pillars.append(p)
        c = read(root / layer / "c4_focus/c4.focus.contrasts.csv")
        if len(c):
            c["layer"] = layer
            c["nominal_interval_excludes_zero"] = (c.delta_lo > 0) | (c.delta_hi < 0)
            c["interpretation"] = "Exploratory; all tested contrasts retained; no subgroup-heterogeneity claim from separate CIs"
            contrasts.append(c)
    combine = lambda xs: pd.concat(xs, ignore_index=True) if xs else pd.DataFrame()
    return dict(inventory=pd.DataFrame(inventory), c1_summary=pd.DataFrame(summary),
                measured_PGS=combine(concordance), event_windows=combine(windows),
                landmark_family=combine(trajectories), pillar_support=combine(pillars),
                all_c4_contrasts=combine(contrasts))


def cell_annotation(root, code_dir, outdir=None):
    """Use the full assayed universe, not a significant-only background."""
    import subprocess
    import sys
    output = outdir or root
    c1 = root / "prot/c1_correlate"
    d = read(c1 / "pwas_incident_adj2.csv")
    if d.empty or not {"term", "FDR"} <= set(d):
        return "C1 protein associations unavailable"
    universe = pd.DataFrame({"assay": d.term, "gene": d.term.str.upper().replace({"NTPROBNP": "NPPB"})})
    panels = [pd.DataFrame({"feature": d.loc[d.FDR < .05, "term"], "model": "C1_incident_FDR05"})]
    lm = read(c1 / "pwas_incident_landmark_adj2.csv")
    if len(lm):
        lm["q_full_family"] = bh(lm["p.value"], len(d) * lm.landmark_years.nunique())
        z = lm[(lm.landmark_years == 5) & (lm.q_full_family < .05)]
        panels.append(pd.DataFrame({"feature": z.term, "model": "C1_landmark5_full_family_FDR05"}))
    f = read(root / "prot/c4_focus/c4.focus.panel_members.csv")
    if {"feature", "model"} <= set(f):
        panels.append(f[["feature", "model"]].dropna())
    universe.to_csv(output / "c5.systematic.cell_universe.csv", index=False)
    pd.concat(panels, ignore_index=True).to_csv(output / "c5.systematic.cell_panels.csv", index=False)
    cmd = [sys.executable, str(code_dir / "c3_cell_annotation.py"), "--universe",
           str(output / "c5.systematic.cell_universe.csv"), "--panels", str(output / "c5.systematic.cell_panels.csv"),
           "--atlas", str(code_dir.parent / "data/cellage/atlas.csv"), "--outdir", str(output), "--prefix", "c5.systematic.cell"]
    subprocess.run(cmd, check=True)
    return "External CellAge cell labels; full-assay background; no inferred tissue of release"


def plots(t, root):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axs = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    for ax, title in zip(axs.flat, ["A  Tested associations", "B  Measured vs biomarker PGS", "C  Event ascertainment windows",
                                   "D  Supported LE8 domains", "E  Overall / subgroup comparisons", "F  Distal association support"]):
        ax.set_title(title, loc="left", fontweight="bold")
    d = t["c1_summary"]
    if len(d):
        labels = d.layer + " / " + d.analysis
        axs[0, 0].barh(labels, d.tested, color="#dce3eb", label="Tested")
        axs[0, 0].barh(labels, d.FDR05, color="#287a78", label="FDR < 0.05")
        axs[0, 0].legend(frameon=False)
    d = t["measured_PGS"]
    if len(d):
        for layer, color in (("prot", "#8d64ab"), ("met", "#df8d35")):
            z = d[d.layer == layer]
            axs[0, 1].scatter(z.beta_measured, z.beta_PGS, s=9, alpha=.4, color=color, label=layer)
        z = d[d.both_supported].sort_values("FDR_PGS").head(3)
        if root.name == "amr":
            z = pd.concat([z, d[d.term == "GlycA"]]).drop_duplicates("term")
        for _, row in z.iterrows():
            axs[0, 1].annotate(row.term, (row.beta_measured, row.beta_PGS), xytext=(4, 6), textcoords="offset points", fontsize=7)
        axs[0, 1].axhline(0, color="grey", lw=.7); axs[0, 1].axvline(0, color="grey", lw=.7)
        axs[0, 1].set(xlabel="Measured log HR / measured SD", ylabel="PGS log HR / PGS SD")
        axs[0, 1].legend(frameon=False)
    d = t["event_windows"]
    if len(d):
        for layer in d.layer.unique():
            z = d[(d.layer == layer) & (d.side == "Post-baseline incident")]
            axs[0, 2].plot(z.window_hi, z.events_max, "o-", label=layer)
        axs[0, 2].set(xlabel="Window upper bound (years)", ylabel="Events in window (feature max)")
        axs[0, 2].legend(frameon=False)
    d = t["pillar_support"]
    if len(d):
        z = d.pivot_table(index="component", columns=["layer", "cohort"], values="n", aggfunc="first").fillna(0)
        im = axs[1, 0].imshow(np.log1p(z), aspect="auto", cmap="Blues")
        axs[1, 0].set_yticks(range(len(z)), z.index)
        axs[1, 0].set_xticks(range(len(z.columns)), [" / ".join(x) for x in z.columns], rotation=25, ha="right")
        for i in range(len(z)):
            for j in range(len(z.columns)):
                axs[1, 0].text(j, i, str(int(z.iloc[i, j])), ha="center", va="center", fontsize=8,
                               color="white" if np.log1p(z.iloc[i, j]) > 4 else "black")
    d = t["all_c4_contrasts"]
    if len(d):
        z = d[(d.layer == "prot") & (d.model == "YSplus_YinYang_50") & (d.reference == "NS_50")]
        labels = z.stratum.str.replace(" baseline inflammation", "") + " / L" + z.landmark.astype(str)
        axs[1, 1].errorbar(z.delta_AUC, range(len(z)), xerr=[np.maximum(0,z.delta_AUC-z.delta_lo),np.maximum(0,z.delta_hi-z.delta_AUC)], fmt="o", color="#287a78")
        axs[1, 1].set_yticks(range(len(z)), labels, fontsize=7)
        axs[1, 1].axvline(0, color="grey", lw=.7)
        axs[1, 1].set_xlabel("YSplus YY vs NS; paired delta AUC (50 assays)")
    d = t["landmark_family"]
    if len(d):
        for layer in d.layer.unique():
            z = d[d.layer == layer].groupby("landmark_years")["FDR_assay_landmark_family"].apply(lambda x: (x < .05).sum())
            axs[1, 2].plot(z.index, z.values, "o-", label=layer)
        axs[1, 2].set(xlabel="Landmark (years)", ylabel="FDR < .05, full assay × landmark family")
        axs[1, 2].legend(frameon=False)
    # A C1-only trait gets meaningful C1 panels in the available space.
    if t["pillar_support"].empty and len(t["measured_PGS"]):
        for ax, layer in ((axs[1, 0], "prot"), (axs[1, 1], "met")):
            z = t["measured_PGS"].query("layer == @layer").nsmallest(8, "FDR_measured")
            pos = np.arange(len(z))
            ax.errorbar(z.beta_measured, pos, xerr=1.96*z["std.error_measured"], fmt="o", color="#287a78", label="Measured")
            ax.errorbar(z.beta_PGS, pos+.18, xerr=1.96*z["std.error_PGS"], fmt="s", color="#a276b4", label="Matched PGS")
            ax.set_yticks(pos, z.term, fontsize=8); ax.axvline(0, color="grey", lw=.7)
            ax.set_title(("D  " if layer=="prot" else "E  ")+layer+" leading measured signals", loc="left", fontweight="bold")
            ax.set_xlabel("Log HR per respective SD; different cohorts/covariates")
            ax.legend(frameon=False, fontsize=8)
    for ax in axs.flat:
        if not ax.has_data():
            ax.text(.5, .5, "Not available for this trait", transform=ax.transAxes, ha="center")
    fig.suptitle(f"{root.name}: systematic 5C evidence audit — descriptive, existing results", fontsize=14)
    for ext in ("png", "pdf"):
        fig.savefig(root / f"c5.systematic.Fig1.evidence_atlas.{ext}", dpi=220)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, help="Separate report location; source tables never modified")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--cell-annotation", action="store_true")
    args = ap.parse_args()
    root = args.root.resolve(); out = (args.outdir or root).resolve(); out.mkdir(parents=True, exist_ok=True)
    tables = aggregate(root)
    for name, d in tables.items():
        d.to_csv(out / f"c5.systematic.{name}.csv", index=False)
    if not args.no_plots:
        plots(tables, out)
    status = dict(source_root=str(root), interpretation="Descriptive consolidation; no independent validation or causal voting")
    if args.cell_annotation:
        status["cell_annotation"] = cell_annotation(root, Path(__file__).resolve().parent, out)
    (out / "c5.systematic.status.json").write_text(json.dumps(status, indent=2) + "\n")


if __name__ == "__main__":
    main()
