#!/usr/bin/env python3
"""Exploratory retention associations, conditional on sampled variable loci."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from parquet_io import selected_rows
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.discrete.conditional_models import ConditionalLogit

from model_utils import grouped_predictions


def fit_conditional(data, columns, output):
    """Report unidentifiable/nonconverged models explicitly, without fake estimates."""
    table = pd.DataFrame(columns=["term", "coef", "se", "z", "p", "OR_per_SD"])
    info = {"n_candidates": len(data), "n_loci": int(data.site_id.nunique())}
    try:
        if data.empty:
            raise ValueError("no informative matched loci")
        raw = data[columns]
        if not np.isfinite(raw.to_numpy(float)).all():
            raise ValueError("missing/nonfinite predictors")
        # A feature constant within every locus cannot enter conditional logit.
        centered = raw - raw.groupby(data.site_id).transform("mean")
        kept = [c for c in columns if np.linalg.norm(centered[c]) > 1e-10]
        info["excluded_constant_terms"] = [c for c in columns if c not in kept]
        if not kept:
            raise ValueError("no within-locus variation in predictors")
        design = StandardScaler().fit_transform(raw[kept])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fit = ConditionalLogit(
                data.observed.to_numpy(), design, groups=data.site_id.to_numpy()
            ).fit(disp=False, maxiter=500)
        if caught:
            raise ValueError("; ".join(str(w.message) for w in caught))
        if not np.isfinite(np.r_[fit.params, fit.bse]).all():
            raise ValueError("nonfinite estimates or standard errors")
        table = pd.DataFrame(
            {
                "term": kept,
                "coef": fit.params,
                "se": fit.bse,
                "z": fit.tvalues,
                "p": fit.pvalues,
                "OR_per_SD": np.exp(fit.params),
            }
        )
        info["status"] = "ok"
    except (ValueError, np.linalg.LinAlgError) as error:
        info.update(status="unavailable", reason=str(error))
        warnings.warn(f"{output.name}: {error}")
    table.to_csv(output, sep="\t", index=False)
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["geometry", "sites", "mutation", "out-dir"]:
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--seed", type=int, default=20261015)
    args = parser.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    geometry = pd.read_parquet(args.geometry)
    sites = pd.read_parquet(args.sites)
    mutation = selected_rows(args.mutation, geometry.site_id)
    if "mutation_prediction_source" not in mutation:
        raise SystemExit("Re-run stage1: held-out mutation scores are required")
    if (
        not mutation.loc[mutation.occupancy < 3, "mutation_prediction_source"]
        .eq("out_of_fold")
        .all()
    ):
        raise SystemExit("Stage 2 requires out-of-fold mutation scores for O<3")
    g = geometry.merge(
        mutation[["variant_id", "mutation_opportunity"]],
        on="variant_id",
        how="left",
        validate="one_to_one",
    )
    if not np.isfinite(g.mutation_opportunity).all():
        raise SystemExit(
            "Missing mutation scores for geometry candidates; re-run stage1"
        )

    # All three ALT candidates can enter a model with outcome-independent features.
    baseline = g[g.occupancy < 3].copy()
    metrics = {
        "all_candidates": fit_conditional(
            baseline,
            ["mutation_opportunity", "effect_magnitude"],
            out / "stage2_conditional_all_candidates.tsv",
        )
    }

    # The observed primary has affinity=cosine=1 and novelty=0 by construction.
    # Including it would encode the response in the predictors. Compare only the
    # remaining observed/unrealized pair at triallelic loci, conditional on anchor.
    pairs = g[(g.occupancy == 2) & (g.is_primary == 0)].copy()
    metrics["conditional"] = fit_conditional(
        pairs,
        [
            "mutation_opportunity",
            "effect_magnitude",
            "affinity_to_primary",
            "cosine_to_primary",
            "novelty_perpendicular",
        ],
        out / "stage2_conditional_candidate.tsv",
    )
    metrics["conditional"]["interpretation"] = (
        "Exploratory secondary-vs-unrealized association at triallelic loci, "
        "conditional on a present-day AF-selected observed primary anchor; not causal retention."
    )

    # Derive summaries from the supplied mutation file, never from stale neighbors.
    ms = (
        mutation.groupby("site_id")
        .agg(
            mutation_mean=("mutation_opportunity", "mean"),
            mutation_max=("mutation_opportunity", "max"),
        )
        .reset_index()
    )
    ss = sites.merge(ms, on="site_id", how="left", validate="one_to_one")
    ss["multi"] = (ss.occupancy >= 2).astype(int)
    columns = [
        "mutation_mean",
        "mutation_max",
        "potential_mean_magnitude",
        "potential_max_magnitude",
        "potential_pairwise_diversity",
    ]
    if not np.isfinite(ss[columns].to_numpy(float)).all():
        raise SystemExit("Missing/nonfinite site-level predictors")
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            max_iter=1000, class_weight="balanced", random_state=args.seed
        ),
    )
    # All transforms are fit within CV. Global geometry scaling and upstream
    # cross-fitting still make this exploratory, not a nested validation estimate.
    pred, auc, status = grouped_predictions(model, ss[columns], ss.multi, ss.chrom)
    ss["pred_multiallelic"] = pred
    coef = pd.DataFrame(columns=["term", "coef", "OR_per_SD"])
    fit_status = "unavailable: outcome has only one class"
    if ss.multi.nunique() == 2:
        model.fit(ss[columns], ss.multi)
        values = model.named_steps["logisticregression"].coef_[0]
        coef = pd.DataFrame(
            {"term": columns, "coef": values, "OR_per_SD": np.exp(values)}
        )
        fit_status = "ok"
    coef.to_csv(out / "stage2_site_multiallelic.tsv", sep="\t", index=False)
    ss.to_parquet(out / "stage2_sites.parquet", index=False)
    metrics.update(
        site_cv_auc=auc,
        site_cv_status=status,
        site_fit_status=fit_status,
        site_n=len(ss),
        site_multiallelic=int(ss.multi.sum()),
        validation_scope="exploratory sampled-cohort CV; not nested end-to-end or population-calibrated",
    )
    tri = g[(g.occupancy == 2) & (g.observed == 1) & (g.is_primary == 0)]
    tri[
        [
            "site_id",
            "variant_id",
            "af",
            "effect_magnitude",
            "distance_to_primary",
            "cosine_to_primary",
            "affinity_to_primary",
            "novelty_perpendicular",
        ]
    ].to_csv(out / "stage2_triallelic_secondary.tsv", sep="\t", index=False)
    sites[sites.occupancy == 3].to_csv(
        out / "stage2_tetra_sites.tsv", sep="\t", index=False
    )
    (out / "stage2_metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(metrics, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
