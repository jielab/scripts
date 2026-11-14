#!/usr/bin/env python3
"""Offline evidence report for the signed PHOLE graph pipeline.

Read-only reporting: no training, testing, imputation or MAGMA fallback.
Missing stages remain unrun; synthetic inputs are explicitly labelled.
"""
from __future__ import annotations
import argparse
import base64
import hashlib
import html
import json
import math
import os
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCHEMA = "phole-signed-report-v2"
MODELS = ["zero", "hotspot", "distance", "selected"]
COLORS = ["#c8d2d9", "#87a7b8", "#739c9d", "#286e88"]


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix-dir", type=Path, required=True, help="Primary signed v2 matrix directory")
    p.add_argument("--model-dir", type=Path, required=True)
    p.add_argument("--discovery-dir", type=Path)
    p.add_argument("--validation-dir", type=Path)
    p.add_argument("--cross-check-dir", type=Path, help="Separate summary cross-check snapshot")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--seed", type=int, default=20261105)
    p.add_argument("--max-points", type=int, default=100000)
    p.add_argument("--heatmap-loci", type=int, default=36)
    p.add_argument("--heatmap-assays", type=int, default=22)
    p.add_argument("--replace", choices=["TRUE", "FALSE"], default="TRUE")
    p.add_argument("--synthetic", action="store_true", help="Mark downstream model outputs as simulation")
    return p.parse_args()


def read_json(path):
    if path is None or not Path(path).is_file():
        return {}
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def read_table(path):
    if path is None or not Path(path).is_file():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[".", "nan", "NaN", ""])
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def records(frame):
    return plain(frame.to_dict("records")) if not frame.empty else []


def atomic_text(path, text):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def number(value, digits=3):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "—"
    if not math.isfinite(value):
        return "—"
    if value == 0:
        return "0"
    if abs(value) >= 10000 or abs(value) < .001:
        return f"{value:.{digits}g}"
    return f"{value:,.{digits}f}".rstrip("0").rstrip(".")


def integer(value):
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError, OverflowError):
        return "—"


def chr_order(value):
    value = str(value).removeprefix("chr").upper()
    try:
        return int(value)
    except ValueError:
        return {"X": 23, "Y": 24, "MT": 25, "M": 25}.get(value, 1000)


def evidence_origin(args, manifests):
    values = [m.get("synthetic") for m in manifests if "synthetic" in m]
    if args.synthetic or any(v is True for v in values):
        return "synthetic"
    if values and all(v is False for v in values):
        return "cohort"
    return "undeclared"


def load_inputs(args):
    mm = read_json(args.matrix_dir / "manifest.json")
    pm = read_json(args.matrix_dir.parent / "manifest.json")
    if pm.get("schema") != "phole-signed-panel-v2":
        pm = {}
    fm = read_json(args.model_dir / "fit.manifest.json")
    dm = read_json(args.discovery_dir / "discovery.manifest.json") if args.discovery_dir else {}
    vm = read_json(args.validation_dir / "validation.manifest.json") if args.validation_dir else {}
    vs = read_json(args.validation_dir / "summary.json") if args.validation_dir else {}
    sm = [read_json(folder / "simulation.manifest.json") for folder in
          (args.matrix_dir, args.matrix_dir.parent, args.model_dir, args.discovery_dir) if folder is not None]
    cs = read_json(args.cross_check_dir / "summary.json") if args.cross_check_dir else {}
    origin = evidence_origin(args, [mm, pm, fm, dm, vm, vs, *sm])
    warnings, arrays = [], {}
    source = target = pd.DataFrame()
    matrix_ready = False
    if mm:
        if mm.get("schema", mm.get("format")) != "phole-signed-matrix-v2":
            raise ValueError("Only phole-signed-matrix-v2 is accepted; MAGMA scores are not signed beta")
        if not mm.get("directional") or mm.get("input_scope") != "full":
            raise ValueError("Matrix must identify directional full-summary input")
        required = ["beta", "se", "observed", "cis"]
        missing = [name + ".npy" for name in required if not (args.matrix_dir / (name + ".npy")).is_file()]
        missing += [name for name in ("source.tsv", "target.tsv") if not (args.matrix_dir / name).is_file()]
        if missing:
            warnings.append("矩阵阶段未完成，缺少：" + ", ".join(missing))
        else:
            source, target = read_table(args.matrix_dir / "source.tsv"), read_table(args.matrix_dir / "target.tsv")
            for name in required:
                arrays[name] = np.load(args.matrix_dir / (name + ".npy"), mmap_mode="r", allow_pickle=False)
            if any(a.shape != (len(source), len(target)) for a in arrays.values()):
                raise ValueError("Matrix arrays and source/target metadata have different shapes")
            known_file = args.matrix_dir / "target_cis_known.npy"
            if known_file.is_file():
                arrays["known"] = np.load(known_file, mmap_mode="r", allow_pickle=False)
                if arrays["known"].shape != (len(target),):
                    raise ValueError("target_cis_known.npy shape mismatch")
            else:
                arrays["known"] = np.ones(len(target), dtype=bool)
                warnings.append("未提供 target_cis_known.npy；按已提供的 cis mask 显示。")
            matrix_ready = True
    else:
        warnings.append("尚无 signed matrix manifest；未将 .sig.txt 或 MAGMA 当作 β 矩阵。")
    metrics, folds = read_table(args.model_dir / "metrics.tsv"), read_table(args.model_dir / "fold-metrics.tsv")
    fit_ready = False
    if fm:
        if fm.get("schema") != "phole-signed-fit-v2":
            raise ValueError("Only phole-signed-fit-v2 model results are accepted")
        if matrix_ready and all((args.model_dir / f).is_file() for f in
                                ("prediction.oof.npy", "residual.oof.npy", "metrics.tsv")):
            for name, file in (("prediction", "prediction.oof.npy"), ("residual", "residual.oof.npy")):
                arrays[name] = np.load(args.model_dir / file, mmap_mode="r", allow_pickle=False)
                if arrays[name].shape != arrays["beta"].shape:
                    raise ValueError(f"{file} and primary matrix shapes differ")
            fit_ready = True
        else:
            warnings.append("fit manifest 存在，但矩阵、OOF 数组或 metrics 未齐备，不标为模型已完成。")
    elif not metrics.empty:
        warnings.append("metrics.tsv 缺少 v2 fit manifest，不作为已完成模型的性能展示。")
    if not fit_ready:
        metrics = folds = pd.DataFrame()
    candidates = read_table(args.discovery_dir / "candidates.tsv") if args.discovery_dir else pd.DataFrame()
    modules = read_table(args.discovery_dir / "latent_modules.tsv") if args.discovery_dir else pd.DataFrame()
    discovery_ready = bool(args.discovery_dir and (args.discovery_dir / "candidates.tsv").is_file() and fit_ready)
    if not discovery_ready:
        candidates = modules = pd.DataFrame()
    replication = read_table(args.validation_dir / "replication.tsv") if args.validation_dir else pd.DataFrame()
    validation_ready = bool(args.validation_dir and (args.validation_dir / "replication.tsv").is_file() and fit_ready)
    if not validation_ready:
        replication = pd.DataFrame()
    if origin == "undeclared" and matrix_ready:
        warnings.append("manifest 未显式声明 synthetic=true/false，不自动认定为真实队列结果。")
    return dict(matrix_meta=mm, panel_meta=pm, fit_meta=fm, discovery_meta=dm, validation_meta=vm,
                validation_summary=vs, cross_summary=cs, origin=origin, warnings=warnings, arrays=arrays,
                source=source, target=target, matrix_ready=matrix_ready, fit_ready=fit_ready,
                discovery_ready=discovery_ready, validation_ready=validation_ready, metrics=metrics,
                fold_metrics=folds, candidates=candidates, modules=modules, replication=replication)


def sample_effects(run, args):
    """Effect-independent uniform reservoir sampling; arrays stay memory mapped."""
    if not run["matrix_ready"]:
        return {}
    a = run["arrays"]
    nr, nc = a["beta"].shape
    row_counts, col_counts = np.zeros(nr, dtype=np.int64), np.zeros(nc, dtype=np.int64)
    rng = np.random.default_rng(args.seed)
    sample = {k: np.empty(0) for k in ("key", "obs", "pred", "residual", "se")}
    n_observed = n_trans = n_oof = 0
    for start in range(0, nr, 256):
        stop = min(nr, start + 256)
        beta, se = np.asarray(a["beta"][start:stop]), np.asarray(a["se"][start:stop])
        observed = np.asarray(a["observed"][start:stop], dtype=bool) & np.isfinite(beta) & np.isfinite(se) & (se > 0)
        n_observed += int(observed.sum())
        eligible = observed & ~np.asarray(a["cis"][start:stop], dtype=bool) & np.asarray(a["known"], dtype=bool)[None, :]
        row_counts[start:stop] = eligible.sum(axis=1)
        col_counts += eligible.sum(axis=0)
        n_trans += int(eligible.sum())
        if not run["fit_ready"]:
            continue
        pred, res = np.asarray(a["prediction"][start:stop]), np.asarray(a["residual"][start:stop])
        keep = eligible & np.isfinite(pred) & np.isfinite(res)
        n = int(keep.sum())
        n_oof += n
        if not n:
            continue
        obs, pred, res, std = beta[keep], pred[keep], res[keep], se[keep]
        if not np.allclose(obs - pred, res, rtol=1e-4, atol=1e-6):
            raise ValueError("OOF residuals disagree with observed beta minus OOF prediction")
        keys = rng.random(n)
        if n > args.max_points:
            take = np.argpartition(keys, args.max_points - 1)[:args.max_points]
            keys, obs, pred, res, std = [x[take] for x in (keys, obs, pred, res, std)]
        added = dict(key=keys, obs=obs, pred=pred, residual=res, se=std)
        sample = {k: np.concatenate([sample[k], added[k]]) for k in sample}
        if sample["key"].size > args.max_points:
            take = np.argpartition(sample["key"], args.max_points - 1)[:args.max_points]
            sample = {k: value[take] for k, value in sample.items()}
    run["counts"] = dict(n_source=nr, n_assay=nc,
                        n_protein=int(run["target"]["PROT"].nunique()) if "PROT" in run["target"] else nc,
                        n_observed=n_observed, n_trans=n_trans, n_oof=n_oof, n_missing=nr * nc - n_observed)
    expected = run["fit_meta"].get("n_observed_trans")
    if run["fit_ready"] and expected is not None and int(expected) != n_oof:
        raise ValueError("OOF edge count differs from fit manifest; do not combine stale runs")
    source, target = run["source"], run["target"]
    rows = sorted(np.flatnonzero(row_counts), key=lambda i: (chr_order(source.iloc[i].get("CHR", "")), float(source.iloc[i].get("POS", 0))))
    cols = sorted(np.flatnonzero(col_counts), key=lambda j: (str(target.iloc[j].get("PROT", "")), str(target.iloc[j].get("ASSAY", ""))))

    def spaced(values, limit):
        return np.asarray(values, dtype=int) if len(values) <= limit else np.asarray(values, dtype=int)[np.linspace(0, len(values) - 1, limit).round().astype(int)]

    ri, ci = spaced(rows, args.heatmap_loci), spaced(cols, args.heatmap_assays)
    heat = {"rows": ri, "cols": ci}
    if ri.size and ci.size:
        ix = np.ix_(ri, ci)
        beta, se = np.asarray(a["beta"][ix]), np.asarray(a["se"][ix])
        valid = np.asarray(a["observed"][ix], dtype=bool) & np.isfinite(beta) & np.isfinite(se) & (se > 0)
        valid &= ~np.asarray(a["cis"][ix], dtype=bool) & np.asarray(a["known"], dtype=bool)[ci][None, :]
        heat["observed"] = np.where(valid, beta, np.nan)
        if run["fit_ready"]:
            for key in ("prediction", "residual"):
                raw = np.asarray(a[key][ix])
                heat[key] = np.where(valid & np.isfinite(raw), raw, np.nan)
    sample["heat"] = heat
    return sample


def save_figure(fig, out, name, run):
    if run["origin"] == "synthetic":
        fig.text(.995, .995, "SYNTHETIC / PIPELINE TEST ONLY", ha="right", va="top", fontsize=10,
                 color="#a76a28", weight="bold", bbox=dict(facecolor="#fff5e5", edgecolor="none", pad=5))
    elif run["origin"] == "undeclared" and run["matrix_ready"]:
        fig.text(.995, .995, "INPUT ORIGIN NOT DECLARED", ha="right", va="top", fontsize=8, color="#8c7b67")
    temporary = out / (name + ".tmp.png")
    fig.savefig(temporary, dpi=165, facecolor="white", bbox_inches="tight")
    os.replace(temporary, out / (name + ".png"))
    plt.close(fig)
    return out / (name + ".png")


def no_data(ax, title, text):
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=12)
    ax.text(.5, .5, text, ha="center", va="center", transform=ax.transAxes, fontsize=11, color="#8295a0", wrap=True)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#dfe7eb")


def primary_metrics(frame):
    if frame.empty or "model" not in frame:
        return pd.DataFrame()
    if "scope" in frame:
        choice = frame[frame.scope.eq("grouped_out_of_fold_trans")]
        if choice.empty:
            choice = frame[frame.scope.eq(frame.iloc[0]["scope"])]
    else:
        choice = frame
    return choice.drop_duplicates("model").set_index("model").reindex(MODELS)


def metric_panel(ax, run):
    frame = primary_metrics(run["metrics"])
    if frame.empty or "rmse" not in frame:
        no_data(ax, "Out-of-fold model comparison", "MODEL NOT RUN\nNo performance values are inferred.")
        return
    vals = pd.to_numeric(frame.rmse, errors="coerce").to_numpy()
    xs = np.arange(4)
    ax.bar(xs, vals, color=COLORS, width=.62, zorder=2)
    for x, y in zip(xs, vals):
        if np.isfinite(y):
            ax.text(x, y, number(y), ha="center", va="bottom", fontsize=8)
    folds = run["fold_metrics"]
    if not folds.empty and {"model", "rmse"}.issubset(folds):
        for i, model in enumerate(MODELS):
            v = pd.to_numeric(folds.loc[folds.model.eq(model), "rmse"], errors="coerce").dropna().to_numpy()
            ax.scatter(np.full(len(v), i) + np.linspace(-.13, .13, len(v)), v, s=15, color="#294958",
                       edgecolor="white", linewidth=.4, alpha=.7, zorder=3)
    ax.set_xticks(xs, [m.capitalize() for m in MODELS], fontsize=9)
    ax.set_ylabel("RMSE (beta units; lower is better)")
    ax.set_title("Same held-out edges, four models", loc="left", fontsize=12, fontweight="bold", pad=12)
    ax.grid(axis="y", color="#e7edf0", linewidth=.7, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)


def prediction_panel(ax, run, sample):
    if not run["fit_ready"] or not sample.get("obs", np.empty(0)).size:
        no_data(ax, "Observed vs OOF prediction", "No eligible out-of-fold predictions.")
        return None
    obs, pred = sample["obs"], sample["pred"]
    lo, hi = np.quantile(np.r_[obs, pred], [.005, .995])
    if not np.isfinite(lo + hi) or hi <= lo:
        lo, hi = -1., 1.
    image = ax.hexbin(obs, pred, gridsize=46, mincnt=1, bins="log", cmap="Blues",
                      extent=(lo, hi, lo, hi), linewidths=0)
    ax.plot([lo, hi], [lo, hi], "--", color="#a76b52", linewidth=1)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("Observed signed beta"); ax.set_ylabel("OOF predicted beta")
    ax.set_title("Prediction on observed trans edges", loc="left", fontsize=12, fontweight="bold", pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    return image


def effect_limit(values):
    finite = np.abs(values[np.isfinite(values)])
    return max(float(np.quantile(finite, .98)), 1e-9) if finite.size else 1.


def heatmap_panel(ax, values, run, heat, title, labels=True, vmax=None):
    if values is None or not values.size:
        no_data(ax, title, "No eligible observed trans cells.")
        return None
    cmap = matplotlib.colormaps["RdBu_r"].copy(); cmap.set_bad("#e1e8ed")
    vmax = effect_limit(values) if vmax is None else vmax
    image = ax.imshow(np.ma.masked_invalid(values), aspect="auto", interpolation="nearest", cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=12)
    ri, ci = heat["rows"], heat["cols"]
    if labels:
        labs = [str(run["target"].iloc[j].get("ASSAY", run["target"].iloc[j].get("PROT", j))) for j in ci]
        ax.set_xticks(np.arange(len(ci)), labs, rotation=75, ha="right", fontsize=6.5)
        labs = [f"{run['source'].iloc[i].get('CHR', '?')}:{int(run['source'].iloc[i].get('POS', 0))}" for i in ri]
        every = max(1, math.ceil(len(labs) / 12))
        ax.set_yticks(np.arange(0, len(ri), every), labs[::every], fontsize=7)
    else:
        ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("Assays"); ax.set_ylabel("Frozen loci")
    return image


def candidate_panel(ax, run, limit=10):
    frame = run["candidates"].copy()
    if frame.empty or "RESIDUAL_RMS_Z" not in frame:
        text = "No residual-locus candidates in this run." if run["discovery_ready"] else "DISCOVERY NOT RUN\nNo candidate claims are inferred."
        no_data(ax, "Exploratory residual loci", text)
        return
    key = "rank" if "rank" in frame else "RESIDUAL_RMS_Z"
    frame["_order"] = pd.to_numeric(frame[key], errors="coerce")
    frame = frame.sort_values("_order", ascending=key == "rank").head(limit)
    vals = pd.to_numeric(frame.RESIDUAL_RMS_Z, errors="coerce").to_numpy()
    labs = [str(row.get("CANDIDATE", row.get("LOCUS", i))) for i, row in frame.iterrows()]
    labs = [s if len(s) <= 33 else s[:30] + "..." for s in labs]
    ax.barh(np.arange(len(vals)), vals, color="#668e9e", height=.65)
    ax.set_yticks(np.arange(len(vals)), labs, fontsize=7.5); ax.invert_yaxis()
    ax.set_xlabel("RMS of OOF residual / SE")
    ax.set_title("Exploratory residual loci", loc="left", fontsize=12, fontweight="bold", pad=12)
    ax.grid(axis="x", color="#e7edf0", linewidth=.7); ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)


def make_figures(run, args, sample):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.labelcolor": "#395c70",
                         "text.color": "#214454", "xtick.color": "#657d8b", "ytick.color": "#657d8b", "axes.edgecolor": "#cfdae1"})
    figures = {}
    fig, panels = plt.subplots(2, 2, figsize=(15.8, 10.8))
    metric_panel(panels[0, 0], run); prediction_panel(panels[0, 1], run, sample)
    heat = sample.get("heat", {})
    heatmap_panel(panels[1, 0], heat.get("observed"), run, heat, "Observed primary beta: systematic sample", labels=False)
    candidate_panel(panels[1, 1], run, 10)
    dataset = str(run["matrix_meta"].get("dataset", "primary")).upper()
    fig.suptitle("PHOLE / signed proteome-effect graph   |   " + dataset, x=.06, ha="left", fontsize=18, fontweight="bold", y=.98)
    c = run.get("counts", {})
    subtitle = f"{integer(c.get('n_source'))} loci   /   {integer(c.get('n_assay'))} assays   /   {integer(c.get('n_oof'))} observed OOF trans edges" if run["matrix_ready"] else "MODEL NOT RUN — input stages are reported explicitly"
    fig.text(.06, .943, subtitle, fontsize=11, color="#6c8796")
    fig.text(.06, .012, "Gray heatmap cells are masked, never filled with zero. Residual candidates are exploratory, not validated discoveries.", fontsize=9, color="#7d919d")
    fig.subplots_adjust(left=.10, right=.97, top=.885, bottom=.08, hspace=.42, wspace=.35)
    figures["overview"] = save_figure(fig, args.out, "overview", run)
    if run["matrix_ready"]:
        fig, panels = plt.subplots(1, 3, figsize=(17.5, 7.2))
        v = effect_limit(heat["observed"]) if "observed" in heat else 1.
        for ax, key, title in zip(panels, ("observed", "prediction", "residual"),
                                  ("Observed signed beta (trans)", "Out-of-fold prediction", "Out-of-fold residual")):
            im = heatmap_panel(ax, heat.get(key), run, heat, title, vmax=v if key != "residual" else None)
            if im is not None:
                fig.colorbar(im, ax=ax, shrink=.58, pad=.02, label="beta" if key != "residual" else "beta residual")
        fig.text(.05, .015, "Same loci and assays in all panels. Gray = unobserved / excluded cis / unavailable OOF. Colors clipped at 98% of shown absolute effects.", fontsize=9, color="#6d8493")
        fig.subplots_adjust(left=.075, right=.97, top=.90, bottom=.25, wspace=.48)
        figures["heatmap"] = save_figure(fig, args.out, "heatmap", run)
    if run["fit_ready"]:
        fig, panels = plt.subplots(1, 2, figsize=(12.7, 5.2))
        im = prediction_panel(panels[0], run, sample)
        if im is not None:
            fig.colorbar(im, ax=panels[0], shrink=.75, label="sampled cell count (log scale)")
        if sample.get("residual", np.empty(0)).size:
            z = sample["residual"] / sample["se"]; lim = max(float(np.quantile(np.abs(z), .99)), 1.)
            panels[1].hist(z, bins=65, range=(-lim, lim), color="#6e9caa", edgecolor="white", linewidth=.3)
            panels[1].axvline(0, color="#ae785c", linestyle="--", linewidth=1)
            panels[1].set_xlabel("(Observed beta - OOF prediction) / SE"); panels[1].set_ylabel("Sampled edges")
            panels[1].set_title("OOF residual distribution", loc="left", fontweight="bold")
            panels[1].spines[["top", "right"]].set_visible(False)
        else:
            no_data(panels[1], "OOF residual distribution", "No eligible predictions.")
        fig.text(.06, .015, f"Uniform effect-independent sample: {len(sample.get('obs', [])):,} / {run['counts']['n_oof']:,} OOF edges; axes trimmed for readability.", fontsize=8.5, color="#718a98")
        fig.subplots_adjust(left=.08, right=.97, top=.88, bottom=.18, wspace=.40)
        figures["prediction"] = save_figure(fig, args.out, "prediction", run)
    if run["discovery_ready"] and not run["candidates"].empty:
        fig, ax = plt.subplots(figsize=(10.6, max(5.2, min(len(run["candidates"]), 20) * .29 + 1.6)))
        candidate_panel(ax, run, 20)
        fig.text(.05, .014, "Ranking describes unexplained patterns; it is not a causal or replication claim.", fontsize=9, color="#778e9c")
        fig.subplots_adjust(left=.28, right=.96, top=.90, bottom=.13)
        figures["candidates"] = save_figure(fig, args.out, "candidates", run)
    return figures


def image_tag(path, alt):
    if path is None:
        return ""
    return '<img src="data:image/png;base64,' + base64.b64encode(path.read_bytes()).decode("ascii") + '" alt="' + html.escape(alt, quote=True) + '">'


def metric_table(run):
    if run["metrics"].empty:
        return '<div class="empty">模型尚未完成，不填入性能数值。</div>'
    columns = ["scope", "model", "n", "rmse", "r2", "pearson", "sign_agreement", "rmse_se"]
    labels = ["评估范围", "模型", "观测边数", "RMSE ↓", "R²", "Pearson r", "符号一致率", "RMSE / SE"]
    rows = []
    for row in records(run["metrics"]):
        cells = []
        for key in columns:
            value = row.get(key)
            if key in ("scope", "model"): text = str(value or "—")
            elif key == "n": text = integer(value)
            elif key == "sign_agreement" and value is not None: text = number(float(value) * 100, 1) + "%"
            else: text = number(value, 4)
            cells.append("<td>" + html.escape(text) + "</td>")
        rows.append('<tr class="' + ("selected" if row.get("model") == "selected" else "") + '">' + "".join(cells) + "</tr>")
    return '<div class="tablewrap"><table><thead><tr>' + "".join("<th>" + s + "</th>" for s in labels) + "</tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"


def model_conclusion(run):
    frame = primary_metrics(run["metrics"])
    if frame.empty or "rmse" not in frame:
        return "尚无可比较的 OOF 性能结果。"
    vals = pd.to_numeric(frame.rmse, errors="coerce")
    base, selected = vals.reindex(["zero", "hotspot", "distance"]).dropna(), vals.get("selected")
    if base.empty or selected is None or not np.isfinite(selected) or base.min() <= 0:
        return "性能以完整指标表为准；本次无法计算相对最佳基线的 RMSE 变化。"
    best = base.idxmin(); change = (float(base.min()) - float(selected)) / float(base.min()) * 100
    if change > 0:
        conclusion = f"所选模型的 OOF RMSE 比当前最好的 {best} 基线低 {number(change, 2)}%。"
    else:
        conclusion = f"所选模型尚未优于当前最好的 {best} 基线；其 OOF RMSE 高 {number(-change, 2)}%。"
    gain = run.get("fit_meta", {}).get("predictive_structure", {})
    required = [gain.get(key) for key in ("loss_reduction", "ci_low", "ci_high")]
    if all(isinstance(x, (int, float)) and np.isfinite(x) for x in required):
        conclusion += (f" 相对 hotspot 的标准化损失减少 {number(required[0], 4)}，"
                       f"双向 block/protein bootstrap 95% CI [{number(required[1], 4)}, {number(required[2], 4)}]。")
        if gain.get("predictive_structure_supported") is True:
            conclusion += " 本次内部矩阵补全通过预设的改进门槛；仍需独立队列检验。"
        else:
            conclusion += " 当前未通过预测结构支持门槛；点估计改善不能证明可推广的规律。"
    else:
        conclusion += " 尚无可核验的 bootstrap 改进区间，未据此声明统计显著或外部泛化成立。"
    return conclusion + " 残差结果按探索性线索解释。"


def stage(title, state, detail, done=False):
    return '<div class="stage ' + ("done" if done else "") + '"><span>' + html.escape(title) + '</span><b>' + html.escape(state) + '</b><small>' + html.escape(detail) + '</small></div>'


def make_report(run, args, figures, generated):
    origin = run["origin"]
    origin_label = {"synthetic": "模拟数据运行", "cohort": "真实队列输入", "undeclared": "输入来源未显式声明"}[origin]
    if not run["matrix_ready"]: origin_label = "后续模型尚未运行"
    banner = ("本页后续矩阵、模型、候选与验证均来自 synthetic 模拟数据，仅用于验证流程。它们不是本次真实蛋白组分析的发现。"
              if origin == "synthetic" else "本页逐阶段标记已有产物。残差位点与潜在模块属于探索性候选，不能直接称为新的生物学机制或已验证的“black hole”。")
    cs = run["cross_summary"]; has_cross = bool(cs.get("datasets")); csc = cs.get("counts", {})
    stages = [
        stage("摘要跨队列核查", "已提供独立快照" if has_cross else "本页未附快照",
              integer(csc.get("proteins")) + " 蛋白 · " + integer(csc.get("priority")) + " 优先候选" if has_cross else "请查看 cross-check.html 的真实摘要结果", has_cross),
        stage("Signed β 矩阵", "已构建" if run["matrix_ready"] else "未完成", origin_label, run["matrix_ready"]),
        stage("分组 OOF 模型", "已运行" if run["fit_ready"] else "未运行",
              str(run["fit_meta"].get("folds", "—")) + " 折 · rank " + str(run["fit_meta"].get("rank", "—")) if run["fit_ready"] else "不推断性能", run["fit_ready"]),
        stage("残差候选筛选", "已生成" if run["discovery_ready"] else "未运行",
              integer(len(run["candidates"])) + " 候选记录" if run["discovery_ready"] else "不推断发现", run["discovery_ready"]),
        stage("独立队列验证", "已有结果" if run["validation_ready"] else "未运行",
              integer(len(run["replication"])) + " 候选 × 队列记录" if run["validation_ready"] else "未验证不写作未复现", run["validation_ready"])]
    c = run.get("counts", {})
    stats = [("冻结位点", c.get("n_source")), ("检测 / assay", c.get("n_assay")),
             ("已观测 OOF trans 边", c.get("n_oof")), ("未观测矩阵单元", c.get("n_missing"))]
    cards = "".join('<div class="metric"><span>' + label + '</span><b>' + integer(value) + '</b></div>' for label, value in stats)
    warnings = '<details class="notes"><summary>输入状态与解释边界</summary><ul>' + "".join("<li>" + html.escape(s) + "</li>" for s in run["warnings"]) + "</ul></details>" if run["warnings"] else ""
    payload = plain(dict(origin=origin, candidates=records(run["candidates"]), modules=records(run["modules"]),
                         replication=records(run["replication"]), fold_metrics=records(run["fold_metrics"]), generated=generated))
    provenance = [
        ("数据集", run["matrix_meta"].get("dataset", "—")), ("输入来源", origin_label),
        ("效应定义", run["matrix_meta"].get("effect_semantics", "有符号 β；按 source effect allele 协调")),
        ("评估范围", run["fit_meta"].get("evaluation", "未运行")),
        ("面板选择边界", run["panel_meta"].get("internal_evaluation_scope", run["fit_meta"].get("panel_selection", "未提供"))),
        ("矩阵目录", str(args.matrix_dir)), ("模型目录", str(args.model_dir)),
        ("候选目录", str(args.discovery_dir or "未提供")), ("验证目录", str(args.validation_dir or "未提供"))]
    prov = "".join("<tr><th>" + html.escape(k) + "</th><td>" + html.escape(str(v)) + "</td></tr>" for k, v in provenance)
    repl = {
        "__GENERATED__": html.escape(generated), "__ORIGIN__": html.escape(origin_label),
        "__BANNER_CLASS__": "synthetic" if origin == "synthetic" else "", "__BANNER__": html.escape(banner),
        "__STAGES__": "".join(stages), "__STATS__": cards,
        "__OVERVIEW__": image_tag(figures["overview"], "signed graph分析概览"),
        "__METRICS__": metric_table(run), "__CONCLUSION__": html.escape(model_conclusion(run)),
        "__HEATMAP__": image_tag(figures.get("heatmap"), "观测beta、OOF预测及残差热图") or '<div class="empty">矩阵尚未构建；没有将未知单元填为零来绘图。</div>',
        "__PREDICTION__": image_tag(figures.get("prediction"), "OOF预测与残差分布") or '<div class="empty">模型未运行，没有 OOF 预测。</div>',
        "__CANDIDATE_FIGURE__": image_tag(figures.get("candidates"), "探索性残差位点"),
        "__WARNINGS__": warnings, "__PROVENANCE__": prov,
        "__PAYLOAD__": json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")}
    result = REPORT_TEMPLATE
    for marker, value in repl.items():
        result = result.replace(marker, value)
    return result


REPORT_TEMPLATE = r"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>PHOLE · Signed graph evidence</title>
<style>
:root{font-family:Inter,"Segoe UI","Microsoft YaHei",sans-serif;color:#234252;background:#f3f6f8;font-size:14px}*{box-sizing:border-box}body{margin:0;line-height:1.65}header{background:#123345;color:#fff;padding:30px max(24px,calc((100vw - 1440px)/2))}header small{letter-spacing:2px;color:#9dc0d0;font-size:11px}h1{font-size:28px;margin:9px 0 7px;line-height:1.3}header p{color:#c4d9e3;margin:0;font-size:13px}.head{display:flex;justify-content:space-between;align-items:center;gap:20px}.origin{font-size:12px;background:#ffffff13;border:1px solid #ffffff28;border-radius:6px;padding:8px 12px;white-space:nowrap}main{max-width:1488px;margin:auto;padding:24px}.banner{padding:14px 18px;border:1px solid #d1e2e8;background:#eaf3f7;color:#456e82;font-size:13px;border-radius:8px;margin-bottom:20px}.banner.synthetic{background:#fff5e5;border-color:#ebd4ad;color:#8c6027;font-weight:550}.stages{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-bottom:18px}.stage{background:#fff;border:1px solid #dce5eb;border-radius:8px;padding:15px;min-height:109px}.stage span{font-size:11px;color:#7b909e;display:block}.stage b{font-size:17px;display:block;margin:5px 0;color:#7d909b;font-weight:600}.stage.done b{color:#347865}.stage small{font-size:10px;color:#8195a0}.stats{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-bottom:22px}.metric{border:1px solid #dce5eb;border-radius:8px;background:#fff;padding:15px 18px}.metric span{font-size:11px;color:#718b9a;display:block}.metric b{font-size:27px;font-weight:600;letter-spacing:-.5px}.card{background:#fff;border:1px solid #dce5eb;border-radius:10px;overflow:hidden;margin-bottom:22px}.pad{padding:21px 23px}.kicker{font-size:11px;color:#3d7e9b;letter-spacing:1.3px;font-weight:650}h2{font-size:20px;margin:4px 0 7px;font-weight:650}.sub,.note{color:#78909f;font-size:12px}.sub{margin:3px 0 0}.note{padding:0 23px 19px}.callout{background:#f0f6f9;border-left:3px solid #6b9aaf;padding:12px 15px;font-size:13px;margin-top:17px}.card img{display:block;width:100%;height:auto}.figure{padding:0 8px}.tablewrap{overflow:auto;max-height:540px}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}th{background:#f1f5f7;color:#5c7b8d;font-size:11px;font-weight:600;padding:12px;text-align:left;white-space:nowrap;position:sticky;top:0;z-index:1}td{padding:10px 12px;border-bottom:1px solid #eaf0f3;font-size:11px;color:#456477;white-space:nowrap}.selected td{background:#edf6f3;font-weight:550}.toolbar{display:flex;justify-content:space-between;gap:12px;align-items:center;flex-wrap:wrap;margin-top:17px}input,button,select{font:inherit;border:1px solid #cddde6;border-radius:5px;padding:7px 10px;background:#fff;color:#456b80;font-size:12px}input{width:280px}button{cursor:pointer}button:hover:not(:disabled){background:#edf5f8}button:disabled{opacity:.35;cursor:default}input:focus-visible,button:focus-visible,summary:focus-visible{outline:3px solid #83b8d4;outline-offset:2px}.pager{padding:12px 20px;background:#fafcfd;display:flex;justify-content:space-between;gap:12px;align-items:center;color:#7c909d;font-size:11px;border-top:1px solid #e2ebf0}.pager button+button{margin-left:5px}.empty{padding:28px 23px;color:#8b9da8;font-size:13px}.tag{display:inline-block;font-size:10px;padding:2px 6px;background:#edf3f6;color:#5b8093;border-radius:4px}.tag.good{background:#e6f4eb;color:#337d59}.tag.warn{background:#fff4e5;color:#a17c40}.notes{margin:20px 0;color:#728b9b;font-size:12px}.notes summary{cursor:pointer;color:#416c84}.notes ul{padding-left:20px}.provenance td{white-space:normal;overflow-wrap:anywhere}.provenance th{position:static;min-width:130px;vertical-align:top;font-size:11px}.footer{font-size:10px;color:#8fa2ae;margin:20px 3px 0}.footnotes{font-size:12px;color:#79909e;margin-top:14px}.footnotes p{margin:9px 0}details.fold{padding:16px 23px;color:#698699;font-size:12px}details.fold summary{cursor:pointer}details.fold .tablewrap{margin-top:15px}@media(max-width:950px){.stages{grid-template-columns:repeat(3,1fr)}.stats{grid-template-columns:repeat(2,1fr)}main{padding:17px}.origin{display:none}}@media(max-width:600px){.stages{grid-template-columns:repeat(2,1fr)}header{padding:25px 20px}h1{font-size:23px}.pad{padding:18px}.toolbar input{width:100%}.stats{gap:8px}.metric b{font-size:24px}.tablewrap{max-height:440px}.pager{padding:11px 15px}.card img{min-width:640px}.figure{overflow:auto}}
</style></head><body>
<header><div class="head"><div><small>PHOLE / SIGNED PROTEOME-EFFECT GRAPH</small><h1>从有符号效应到可复核的残差候选</h1><p>区分已经运行的分析、已有的验证证据与尚未检验的解释。</p></div><div class="origin">__ORIGIN__</div></div></header>
<main><div class="banner __BANNER_CLASS__">__BANNER__</div><section class="stages" aria-label="各阶段状态">__STAGES__</section><section class="stats" aria-label="输入规模">__STATS__</section>
<section class="card"><div class="pad"><div class="kicker">01 / 分析概览</div><h2>当前这一次运行，支持什么结论？</h2><p class="sub">摘要跨队列核查与后续 signed graph 分开计数。模拟产物始终明确标记。</p></div><div class="figure">__OVERVIEW__</div></section>
<section class="card"><div class="pad"><div class="kicker">02 / 预测性能</div><h2>在同一批留出边上比较模型</h2><p class="sub">zero、hotspot、distance 与 selected 使用同一评估范围。这里评估分组矩阵补全，并非从未见过蛋白的 cold-start 预测。当前 fit 使用 primary 完整矩阵的可观测 trans 单元；edge-check 的 stable-cis/trans 文件作为证据输出，没有自动成为训练筛选器。</p><div class="callout">__CONCLUSION__</div></div>__METRICS__<div class="figure">__PREDICTION__</div><p class="note">符号一致率按模型产物的定义报告；不同 assay 的 β 单位可能不同。池化性能不能代替逐蛋白检验或独立队列验证。</p><details class="fold"><summary>查看逐折指标</summary><div id="foldTable"></div></details></section>
<section class="card"><div class="pad"><div class="kicker">03 / 直接看矩阵</div><h2>观测 β、OOF 预测与残差</h2><p class="sub">按基因组顺序与 assay 名称等距抽取小样本；不按效应强弱挑选。三张图使用相同的行和列。</p></div><div class="figure">__HEATMAP__</div><p class="note">灰色表示未观测、cis 未纳入当前 trans 比较，或没有 OOF 预测；没有把未知 β 填为 0。热图保留原始 β，跨 assay 的颜色绝对大小不应解释为相同的生物学效应单位。</p></section>
<section class="card"><div class="pad"><div class="kicker">04 / 探索性候选</div><h2>模型尚未解释的位点</h2><p class="sub">大残差意味着当前模型解释不足，原因也可能是遗漏变量、测量误差或模型形式。候选名称不构成因果发现。</p><div class="toolbar"><input type="search" id="candidateSearch" placeholder="搜索候选、位点或邻近已测基因" aria-label="搜索残差候选"><button id="candidateExport">导出筛选候选 TSV</button></div></div><div class="figure">__CANDIDATE_FIGURE__</div><div id="candidateTable"></div></section>
<section class="card"><div class="pad"><div class="kicker">05 / 潜在模块</div><h2>局部集中，还是分布较广？</h2><p class="sub">有效位点数、有效蛋白数与局部权重描述模块形态；它们不会自行证明存在新的实体或调控机制。</p><div class="toolbar"><input type="search" id="moduleSearch" placeholder="搜索模块或峰值位点" aria-label="搜索潜在模块"><button id="moduleExport">导出模块 TSV</button></div></div><div id="moduleTable"></div></section>
<section class="card"><div class="pad"><div class="kicker">06 / 外部证据</div><h2>独立队列能否支持候选的效应图谱？</h2><p class="sub">每行是一项“候选 × 验证队列”检验。未测、输入不足或无法匹配的记录保留为不可检验，不能改写成阴性。</p><div class="toolbar"><input type="search" id="replicationSearch" placeholder="搜索候选、队列或状态" aria-label="搜索外部复现"><button id="replicationExport">导出验证 TSV</button></div></div><div id="replicationTable"></div><div class="pad footnotes"><p>PROFILE_CORR、符号一致率与经验 P/Q 需结合可比蛋白数、null 数和 STATUS 一起解读。达到复现规则表示当前定义下的效应图谱支持，不能直接解释为因果机制成立。</p></div></section>
__WARNINGS__
<details class="notes"><summary>运行来源与可追溯信息</summary><div class="tablewrap" style="margin-top:12px;max-height:none"><table class="provenance"><tbody>__PROVENANCE__</tbody></table></div></details>
<div class="footer">生成时间：__GENERATED__ · PHOLE signed graph · 本报告和全部图像可离线查看。</div>
</main><script id="payload" type="application/json">__PAYLOAD__</script><script>
'use strict';
const DATA=JSON.parse(document.getElementById('payload').textContent),$=id=>document.getElementById(id),esc=v=>String(v==null?'':v).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function fmt(v){if(v===null||v===undefined||v==='')return '—';if(typeof v!=='number')return String(v);if(!Number.isFinite(v))return '—';if(Number.isInteger(v))return v.toLocaleString('en-US');return Math.abs(v)<.001&&v!==0?v.toExponential(2):Number(v.toPrecision(5)).toString()}
function save(name,rows){if(!rows.length)return;const keys=Object.keys(rows[0]),cell=v=>{let s=String(v==null?'':v).replace(/[\t\r\n]/g,' ');if(/^[=+@]/.test(s)||/^-(?!\d)/.test(s))s="'"+s;return s;},text='\uFEFF'+[keys,...rows.map(r=>keys.map(k=>r[k]))].map(row=>row.map(cell).join('\t')).join('\n')+'\n',u=URL.createObjectURL(new Blob([text],{type:'text/tab-separated-values;charset=utf-8'})),a=document.createElement('a');a.href=u;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(u),2000)}
function tableController(key,id,columns,labels,searchId,exportId){const rows=DATA[key]||[];let page=0,filtered=rows;function render(){const size=40,start=page*size,body=filtered.slice(start,start+size).map(r=>'<tr>'+columns.map(k=>{const value=r[k];if(k==='STATUS'){const s=String(value||''),cls=['replicated','profile_replicated','supported'].includes(s)?'good':/untestable|missing|insufficient|unavailable/.test(s)?'warn':'';return '<td><span class="tag '+cls+'">'+esc(fmt(value))+'</span></td>';}return '<td title="'+esc(fmt(value))+'">'+esc(fmt(value))+'</td>';}).join('')+'</tr>').join('');$(id).innerHTML='<div class="tablewrap"><table><thead><tr>'+labels.map(l=>'<th>'+esc(l)+'</th>').join('')+'</tr></thead><tbody>'+(body||'<tr><td colspan="'+columns.length+'" class="empty">当前没有可展示记录；未运行或未检验不等于阴性。</td></tr>')+'</tbody></table></div><div class="pager"><span>'+filtered.length.toLocaleString()+' 条记录'+(filtered.length?' · 第 '+(page+1)+' / '+Math.ceil(filtered.length/size)+' 页':'')+'</span><div><button data-prev '+(page===0?'disabled':'')+'>上一页</button><button data-next '+(start+size>=filtered.length?'disabled':'')+'>下一页</button></div></div>';$(id).querySelector('[data-prev]').onclick=()=>{page--;render()};$(id).querySelector('[data-next]').onclick=()=>{page++;render()};if(exportId)$(exportId).disabled=!filtered.length;}
if(searchId)$(searchId).oninput=()=>{const q=$(searchId).value.trim().toLowerCase();filtered=rows.filter(r=>!q||Object.values(r).some(v=>String(v||'').toLowerCase().includes(q)));page=0;render()};if(exportId)$(exportId).onclick=()=>save('phole.'+key+'.tsv',filtered);render();}
tableController('candidates','candidateTable',['rank','CANDIDATE','CHR','POS','N_PROTEINS','N_LARGE_RESIDUAL','RESIDUAL_RMS_Z','NEARBY_MEASURED_GENES','STATUS'],['排名','候选','CHR','POS','蛋白数','大残差数','残差 RMS Z','附近已测基因','状态'],'candidateSearch','candidateExport');
tableController('modules','moduleTable',['MODULE','PEAK_LOCUS','CHR','POS','EFFECTIVE_LOCI','EFFECTIVE_PROTEINS','LOCAL_WEIGHT','EXPLAINED_FRACTION','INTERPRETATION'],['模块','峰值位点','CHR','POS','有效位点数','有效蛋白数','局部权重','解释比例','描述'],'moduleSearch','moduleExport');
tableController('replication','replicationTable',['CANDIDATE','DATASET','N_PROTEINS','PROFILE_CORR','SIGN_AGREEMENT','EMPIRICAL_P','EMPIRICAL_Q','N_MATCHED_NULL','STATUS','REASON'],['候选','队列','可比蛋白数','图谱相关','符号一致率','经验 P','经验 Q','匹配 null 数','状态','原因'],'replicationSearch','replicationExport');
tableController('fold_metrics','foldTable',['fold','model','n','rmse','r2','pearson','sign_agreement','rmse_se'],['折','模型','边数','RMSE','R²','Pearson r','符号一致率','RMSE / SE']);
</script></body></html>"""


def fingerprint(args):
    paths = [Path(__file__), args.matrix_dir / "manifest.json", args.matrix_dir.parent / "manifest.json",
             args.matrix_dir / "source.tsv", args.matrix_dir / "target.tsv",
             args.model_dir / "fit.manifest.json", args.model_dir / "metrics.tsv", args.model_dir / "fold-metrics.tsv"]
    paths += [args.matrix_dir / (n + ".npy") for n in ("beta", "se", "observed", "cis", "target_cis_known")]
    paths += [args.model_dir / f for f in ("prediction.oof.npy", "residual.oof.npy")]
    if args.discovery_dir:
        paths += [args.discovery_dir / f for f in ("candidates.tsv", "latent_modules.tsv", "discovery.manifest.json")]
    if args.validation_dir:
        paths += [args.validation_dir / f for f in ("replication.tsv", "summary.json", "validation.manifest.json")]
    if args.cross_check_dir: paths += [args.cross_check_dir / "summary.json"]
    for folder in (args.matrix_dir, args.matrix_dir.parent, args.model_dir, args.discovery_dir):
        if folder is not None: paths.append(folder / "simulation.manifest.json")
    info = []
    for path in paths:
        if path.is_file():
            st = path.stat()
            item = dict(path=str(path.resolve()), size=st.st_size, mtime_ns=st.st_mtime_ns)
            if path.suffix != ".npy": item["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            info.append(item)
        else:
            info.append(dict(path=str(path.resolve()), missing=True))
    config = dict(seed=args.seed, max_points=args.max_points, heatmap_loci=args.heatmap_loci,
                  heatmap_assays=args.heatmap_assays, synthetic=args.synthetic)
    return hashlib.sha256(json.dumps([info, config], sort_keys=True).encode()).hexdigest(), info


def main():
    args = arguments()
    if args.max_points < 10 or args.heatmap_loci < 1 or args.heatmap_assays < 1:
        raise ValueError("max-points must be >=10 and heatmap dimensions positive")
    args.out.mkdir(parents=True, exist_ok=True)
    signature, inputs = fingerprint(args)
    previous = read_json(args.out / "report.manifest.json")
    if previous.get("signature") == signature and args.replace != "TRUE":
        artifacts = previous.get("artifacts", [])
        if artifacts and all((args.out / name).is_file() for name in artifacts):
            print(f"Reuse matching report: {args.out / 'report.html'}"); return
    run = load_inputs(args); sample = sample_effects(run, args); figures = make_figures(run, args, sample)
    for name in ("heatmap.png", "prediction.png", "candidates.png"):
        if name not in {p.name for p in figures.values()}: (args.out / name).unlink(missing_ok=True)
    generated = datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(timespec="seconds")
    atomic_text(args.out / "report.html", make_report(run, args, figures, generated))
    synthetic = True if run["origin"] == "synthetic" else False if run["origin"] == "cohort" else None
    summary = plain(dict(schema=SCHEMA, generated=generated, origin=run["origin"], synthetic=synthetic,
        stages={k: run[k] for k in ("matrix_ready", "fit_ready", "discovery_ready", "validation_ready")},
        cross_check_snapshot=bool(run["cross_summary"].get("datasets")), counts=run.get("counts", {}),
        candidates=len(run["candidates"]), latent_modules=len(run["modules"]), validation_records=len(run["replication"]),
        validation_summary=run["validation_summary"], conclusion=model_conclusion(run), warnings=run["warnings"],
        heatmap_sampling="systematic across chromosome/position ordered loci and alphabetically ordered assays; not selected by effect",
        scatter_sampling="uniform priorities among observed OOF trans edges, independent of effect",
        n_scatter_sample=len(sample.get("obs", [])), missing_values="masked; never filled with zero",
        evaluation=run["fit_meta"].get("evaluation")))
    atomic_text(args.out / "report-summary.json", json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    artifacts = ["report.html", "report-summary.json"] + [p.name for p in figures.values()]
    manifest = dict(schema=SCHEMA, signature=signature, generated=generated, origin=run["origin"],
                    synthetic=synthetic, artifacts=artifacts, inputs=inputs,
                    figure_semantics="signed beta, observed mask, OOF predictions; no causal discovery assertion")
    atomic_text(args.out / "report.manifest.json", json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    atomic_text(args.out / "plot.done", signature + "\n")
    print(f"PHOLE report: {args.out / 'report.html'}")
    print(json.dumps(summary["stages"], ensure_ascii=False))


if __name__ == "__main__":
    main()
