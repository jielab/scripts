"""Evaluation only after label-free memberships/predictions have been frozen."""
from __future__ import annotations

import html
from pathlib import Path

import numpy as np
import pandas as pd


def metrics(y, prediction, task):
    y, p = np.asarray(y, float), np.asarray(prediction, float)
    valid = np.isfinite(y) & np.isfinite(p)
    y, p = y[valid], p[valid]
    out = dict(n_evaluable=len(y), primary_loss=np.nan, mae=np.nan)
    if not len(y):
        return out
    out.update(primary_loss=float(np.mean((y-p)**2)), mae=float(np.mean(np.abs(y-p))))
    if task == "binary":
        from sklearn.metrics import roc_auc_score, average_precision_score
        both = len(np.unique(y)) == 2
        out.update(n_cases=int(y.sum()), event_fraction=float(y.mean()),
                   auroc=float(roc_auc_score(y, p)) if both else np.nan,
                   auprc=float(average_precision_score(y, p)) if both else np.nan)
    else:
        variance = float(np.mean((y-y.mean())**2))
        out.update(rmse=float(np.sqrt(out["primary_loss"])),
                   r2=1-out["primary_loss"]/variance if variance else np.nan,
                   correlation=float(np.corrcoef(y, p)[0, 1]) if np.std(y) and np.std(p) else np.nan)
    return out


def bootstrap_mean(values, groups, repetitions, seed):
    v = np.asarray(values, float)
    valid = np.isfinite(v)
    v, groups = v[valid], np.asarray(groups, str)[valid]
    if not len(v):
        return dict(estimate=np.nan, low=np.nan, high=np.nan, groups=0)
    _, inverse = np.unique(groups, return_inverse=True)
    n = np.bincount(inverse)
    s = np.bincount(inverse, weights=v)
    out = dict(estimate=float(v.mean()), low=np.nan, high=np.nan, groups=len(n))
    if len(n) < 2 or repetitions < 2:
        return out
    rng = np.random.default_rng(seed)
    boot = np.empty(repetitions)
    for i in range(repetitions):
        draws = rng.integers(0, len(n), len(n))
        boot[i] = s[draws].sum()/n[draws].sum()
    out.update(low=float(np.quantile(boot, .025)), high=float(np.quantile(boot, .975)))
    return out


def same_risk_comparison(frame, selected):
    """Overlap-standardized descriptive comparison in frozen prediction bins.

    No query outcome chooses matches or weights. Rejected controls sharing a
    group with selected queries are excluded. Coarse bins need not control
    every determinant of error; the result is a descriptive comparison.
    """
    good = frame.y.notna().to_numpy() & np.isfinite(frame.p_population)
    selected = np.asarray(selected, bool)
    selected_groups = set(frame.loc[selected, "group_id"])
    independent_controls = ~frame.group_id.isin(selected_groups).to_numpy()
    cells = []
    for risk_bin in sorted(frame.risk_bin.unique()):
        members = good & frame.risk_bin.eq(risk_bin).to_numpy()
        a, b = frame.loc[members & selected], frame.loc[members & ~selected & independent_controls]
        if not len(a) or not len(b):
            continue
        cells.append((min(len(a), len(b)), len(a), len(b),
                      float(np.mean((a.y-a.p_population)**2)),
                      float(np.mean((b.y-b.p_population)**2)),
                      a.group_id.nunique(), b.group_id.nunique()))
    n_selected = int(np.sum(good & selected)); n_other = int(np.sum(good & ~selected))
    if not cells:
        return dict(standardized_selected_loss=np.nan, standardized_other_loss=np.nan,
                    standardized_difference=np.nan, overlap_selected_fraction=0.,
                    overlap_other_fraction=0., overlap_bins=0, comparison_status="no_common_risk_bins")
    values = np.asarray(cells)
    weights = values[:, 0]/values[:, 0].sum()
    a, b = float(weights@values[:, 3]), float(weights@values[:, 4])
    enough = values[:, 0].sum() >= 20 and values[:, 5].sum() >= 10 and values[:, 6].sum() >= 10
    return dict(standardized_selected_loss=a, standardized_other_loss=b,
                standardized_difference=a-b,
                overlap_selected_fraction=float(values[:, 1].sum()/n_selected) if n_selected else 0.,
                overlap_other_fraction=float(values[:, 2].sum()/n_other) if n_other else 0.,
                overlap_bins=len(cells),
                comparison_status="descriptive_coarse_risk_standardization" if enough else "insufficient_overlap_counts")


def evaluate(frame, *, task, trait, role, repetitions=200, seed=20261102):
    if frame.id.duplicated().any():
        raise ValueError("Evaluation requires unique frozen query IDs")
    observed = frame.y.notna()
    if task == "binary" and not np.isin(frame.loc[observed, "y"], [0., 1.]).all():
        raise ValueError("Binary evaluation outcomes must be 0/1")
    rows = []
    subsets = [("all_with_explicit_fallback", np.ones(len(frame), bool), "p_population"),
               ("global_donor_prior_all", np.ones(len(frame), bool), "p_global_prior"),
               ("supported_count_only", frame.supported.to_numpy(bool), "p_abm"),
               ("unsupported_explicit_fallback", ~frame.supported.to_numpy(bool), "p_population")]
    subsets += [(state, frame.evidence_state.eq(state).to_numpy(), "p_population") for state in
                ("supported_lower_dispersion", "supported_heterogeneous")]
    for label, selected, prediction in subsets:
        d = frame.loc[selected]
        m = metrics(d.y, d[prediction], task)
        delta = (d.y-d[prediction])**2-(d.y-d.p_global_prior)**2
        ci = bootstrap_mean(delta, d.group_id, repetitions, seed)
        rows.append(dict(trait=trait, role=role, population=label, task=task,
                         loss_name="brier" if task == "binary" else "mse",
                         n_roster=len(frame), n_selected_roster=len(d),
                         actual_coverage=len(d)/len(frame), n_missing_y=int(d.y.isna().sum()),
                         delta_loss_vs_global_prior=ci["estimate"], delta_low=ci["low"], delta_high=ci["high"],
                         bootstrap_unit="declared_group", n_bootstrap_groups=ci["groups"], **m))
    risk = []
    total_cases = float(frame.loc[observed, "y"].sum()) if task == "binary" else np.nan
    total_controls = int(observed.sum())-total_cases if task == "binary" else np.nan
    for column in (c for c in frame if c.startswith("selected_")):
        selected = frame[column].to_numpy(bool)
        d = frame.loc[selected]
        nominal = float(column.removeprefix("selected_").replace("p", "."))
        m = metrics(d.y, d.p_population, task)
        ci = bootstrap_mean((d.y-d.p_population)**2, d.group_id, repetitions, seed)
        selected_cases = float(d.y.sum()) if task == "binary" else np.nan
        selected_controls = int(d.y.notna().sum())-selected_cases if task == "binary" else np.nan
        risk.append(dict(trait=trait, role=role, nominal_calibration_coverage=nominal,
                         actual_test_coverage=float(selected.mean()), n_selected_roster=int(selected.sum()),
                         primary_loss_low=ci["low"], primary_loss_high=ci["high"],
                         case_retention=selected_cases/total_cases if total_cases else np.nan,
                         control_retention=selected_controls/total_controls if total_controls else np.nan,
                         **m, **same_risk_comparison(frame, selected)))
    return pd.DataFrame(rows), pd.DataFrame(risk)


def subgroup_table(frame, *, task, trait, role):
    records = []
    for ancestry, d in frame.groupby("ancestry", dropna=False):
        records.append(dict(trait=trait, role=role, ancestry=str(ancestry), n_roster=len(d),
                            support_coverage=float(d.supported.mean()),
                            mean_twins=float(d.n_twins.mean()),
                            **metrics(d.y, d.p_population, task)))
    return pd.DataFrame(records)


def plot_risk(risk, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for role, d in risk.groupby("role"):
        axes[0].plot(d.actual_test_coverage, d.primary_loss, "o-", label=role)
        if d.case_retention.notna().any():
            axes[1].plot(d.actual_test_coverage, d.case_retention, "o-", label=role+" cases")
            axes[1].plot(d.actual_test_coverage, d.control_retention, "s--", label=role+" controls")
        else:
            axes[1].plot(d.actual_test_coverage, d.overlap_selected_fraction, "o-", label=role+" selected overlap")
            axes[1].plot(d.actual_test_coverage, d.overlap_other_fraction, "s--", label=role+" rejected overlap")
    axes[0].set(xlabel="Actual retained fraction of frozen roster", ylabel="Observed squared loss", title="Support-ranked retention")
    axes[1].set(xlabel="Actual retained fraction of frozen roster", ylabel="Fraction", ylim=(-.02, 1.02), title="Composition / common-risk overlap")
    for ax in axes:
        ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(output, dpi=160); plt.close(fig)


def write_trait_report(output, model, performance, risk, subgroups):
    caveats = (
        "The prediction is the arithmetic mean of every eligible training donor inside one frozen genetic-profile radius. "
        "Binary predictions are observed case counts divided by all donor counts. No fitted regression, CountHeads, "
        "risk-bin error head, or test-outcome donor pruning is used. "
        "Unsupported people receive no ABM prediction; an explicitly labelled global donor mean is retained for population evaluation. "
        "More support is not evidence of lower personal error. The displayed intervals and donor dispersion are descriptive; "
        "they do not guarantee individual or conditional coverage. Coarse same-predicted-risk standardization is descriptive, "
        "and its overlap fractions must be read with the difference. Previously inspected tests remain exploratory."
    )
    body = [f"<h1>GRID count ABM: {html.escape(model['trait'])}</h1>",
            "<p>Source: grid.count_abm.1; current completion receipt only. Historical blocknet and CountHeads reports are separate.</p>",
            f"<p>{html.escape(caveats)}</p>",
            f"<p>Group evidence: {html.escape(str(model['metadata'].get('grouping_status', 'unverified')))}; "
            f"GWAS overlap: {html.escape(str(model['metadata'].get('gwas_overlap_status', 'unverified')))}.</p>",
            "<h2>Population and supported prediction</h2>", performance.to_html(index=False, border=0),
            "<h2>Frozen selection thresholds, actual coverage</h2>", risk.to_html(index=False, border=0),
            '<img src="risk_coverage.png" alt="Risk coverage and composition" style="max-width:100%">',
            "<h2>Ancestry audit</h2>", subgroups.to_html(index=False, border=0),
            "<p>Complete donor memberships, reference outcomes, query identities, and examples are in the local-only directory and are not embedded in this report.</p>"]
    Path(output).write_text("<!doctype html><html><head><meta charset='utf-8'><title>GRID count ABM</title>"
                           "<style>body{font:15px system-ui;margin:32px;max-width:1500px}table{font-size:11px;border-collapse:collapse}"
                           "td,th{padding:5px;border-bottom:1px solid #ddd}h2{margin-top:32px}</style></head><body>"+
                           "\n".join(body)+"</body></html>")
