"""Figures and individual molecular contrasts; all outputs are descriptive."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def matched_pairs(ids, risk, profiles, names, caliper, max_pairs):
    order = np.argsort(risk, kind="stable")
    used, rows = set(), []
    # Greedy examples selected without outcomes. No discrete subtype condition.
    # Sample anchors across the risk distribution, not just its lowest tail.
    for a in np.random.default_rng(2026).permutation(order):
        if a in used:
            continue
        lo = np.searchsorted(risk[order], risk[a]-caliper)
        hi = np.searchsorted(risk[order], risk[a]+caliper, side="right")
        candidates = np.array([v for v in order[lo:hi] if v != a and v not in used], int)
        if not len(candidates):
            continue
        distance = np.linalg.norm(profiles[candidates]-profiles[a], axis=1)
        b = candidates[np.argmax(distance)]
        top = np.argsort(-np.abs(profiles[a]-profiles[b]))[:3]
        rows.append(dict(eid_A=ids[a], eid_B=ids[b], predicted_risk_A=risk[a],
                         predicted_risk_B=risk[b], profile_distance=float(distance.max()),
                         largest_contrasts=";".join(names[j] for j in top)))
        used.update([a, b])
        if len(rows) >= max_pairs:
            break
    return pd.DataFrame(rows, columns=["eid_A", "eid_B", "predicted_risk_A", "predicted_risk_B",
                                       "profile_distance", "largest_contrasts"])


def figures(out):
    import json
    synthetic = json.loads((out/"manifest.json").read_text())["config"].get("demo",False)
    metrics = pd.read_csv(out/"test_metrics.csv")
    curve = pd.read_csv(out/"coverage_curve.csv")
    fig, ax = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    if synthetic:
        fig.suptitle("SYNTHETIC VALIDATION — not UK Biobank results", fontsize=12)
    all_rows = metrics[metrics.subset == "all"].sort_values("AUC_IPCW")
    ax[0].barh(all_rows.model, all_rows.AUC_IPCW, color=["#B95146" if str(s).startswith("panome") else "#547F9A" for s in all_rows.model])
    ax[0].set(xlabel="Held-out IPCW AUC", xlim=(.45, 1), title="Prediction in the full test cohort")
    for name in ["panome", "panome_clinical", "elasticnet", "clinical"]:
        sub = curve[curve.model == name].sort_values("coverage")
        ax[1].plot(sub.coverage, sub.Brier_IPCW, "o-", label=name)
    ax[1].set(xlabel="Fraction with supported molecular matching", ylabel="IPCW Brier score (lower is better)",
              title="Identical participants compared at every threshold")
    ax[1].legend(frameon=False, fontsize=9)
    for ext in ["png", "pdf"]:
        fig.savefig(out/f"Fig_prediction_coverage.{ext}", dpi=180)
    plt.close(fig)
    pairs = pd.read_csv(out/"same_risk_pairs.csv", dtype={"eid_A": str, "eid_B": str})
    profiles = pd.read_csv(out/"molecular_profiles.csv", dtype={0: str})
    if len(pairs):
        idcol = profiles.columns[0]
        profiles[idcol] = profiles[idcol].astype(str)
        pmap = profiles.set_index(idcol)
        example = pairs.sort_values("profile_distance", ascending=False).iloc[0]
        chosen = pmap.loc[[str(example.eid_A), str(example.eid_B)]]
        fig, ax = plt.subplots(figsize=(12, 3.4), layout="constrained")
        lim = max(2., float(np.max(np.abs(chosen.to_numpy()))))
        im = ax.imshow(chosen.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
        ax.set_xticks(np.arange(chosen.shape[1]), chosen.columns, rotation=35, ha="right")
        ax.set_yticks([0,1], [f"A | risk {example.predicted_risk_A:.1%}", f"B | risk {example.predicted_risk_B:.1%}"])
        ax.set_title(("SYNTHETIC: " if synthetic else "")+"Similar elastic-net risk, different molecular profiles (descriptive)")
        fig.colorbar(im, ax=ax, label="Build-standardized block coordinate")
        for ext in ["png", "pdf"]:
            fig.savefig(out/f"Fig_individual_profiles.{ext}", dpi=180)
        plt.close(fig)
    return True


def write_report(out, audit, readiness, selected, horizon):
    m = pd.read_csv(out/"test_metrics.csv")
    allrows = m[m.subset == "all"].sort_values("AUC_IPCW", ascending=False)
    lines = ["# Panome reference matching", "", f"Fixed horizon: {horizon:g} years; death-censored net risk.",
             "", f"Reference selection: {selected}. Panel readiness is a development diagnostic, not proof of accuracy.",
             f"Readiness: {readiness['status']}. Reasons: {', '.join(readiness['reasons']) or 'none'}.",
             "", "| Model | Full-test AUC | IPCW Brier |", "|---|---:|---:|"]
    for row in allrows.itertuples():
        lines.append(f"| {row.model} | {row.AUC_IPCW:.4f} | {row.Brier_IPCW:.5f} |")
    lines += ["", "## Interpretation", "",
              "The 1NN COPY comparator returns a reference person's observed horizon outcome (0/1); it is not an individual risk estimate.",
              "Panome prototypes retain real person IDs. Their local risks borrow outcomes from unselected build-set donors, excluding self and family.",
              "Calibration uses its own development subset. The outer test set is never used for metric learning, panel selection, tuning, calibration or support thresholds.",
              "Coverage curves compare every model on exactly the same people. High performance on accepted people alone does not establish superiority in the population.",
              "Reference fit scores are model-dependent OOF likelihood gains. Unselected people are not invalid data and a failed panel does not show the cohort is unusable.",
              "Same-risk pairs and molecular blocks describe variation; they do not establish causal CAD mechanisms or intervention effects.",
              "IPCW assumes independent censoring. These risks are not competing-risk cumulative incidences. External cohorts, stronger clinical comparators and repeat training remain necessary.",
              "Intervals in paired_contrasts.csv condition on this fitted pipeline; they do not include model-retraining uncertainty.",
              "", "## Files", "", "reference_candidates.csv, reference_panel.csv, reference_readiness.json, panel_tuning.csv, test_individuals.csv, test_reference_matches.csv, same_risk_pairs.csv, molecular_profiles.csv, coverage_curve.csv, paired_contrasts.csv."]
    (out/"REPORT.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
