#!/usr/bin/env python3
"""Cross-fit mutation opportunity using exact counts of repeated feature patterns.

A weighted fit to aggregated patterns uses every candidate, without retaining
whole-genome candidates in RAM. These conditional, class-balanced scores are
relative opportunity, not calibrated mutation rates or selection coefficients.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from parquet_io import ParquetSink, batches

CAT = ["ref_alt", "context"]
NUM = ["transition", "cpg_transition"]
FEATURES = CAT + NUM


class CountModel:
    """Logistic regression equivalent to expanding integer frequency counts."""

    def __init__(self, seed):
        self.encoder = OneHotEncoder(handle_unknown="ignore")
        self.scaler = StandardScaler()
        self.lr = LogisticRegression(max_iter=1000, random_state=seed)

    def fit(self, data):
        count = data["count"].to_numpy(float)
        outcome = data.observed.to_numpy(int)
        totals = np.bincount(outcome, weights=count, minlength=2)
        if (totals == 0).any():
            raise ValueError("a training fold has only one class")
        self.encoder.fit(data[CAT])
        self.scaler.fit(data[NUM], sample_weight=count)
        balanced = count * count.sum() / (2 * totals[outcome])
        self.lr.fit(self.transform(data), outcome, sample_weight=balanced)
        return self

    def transform(self, data):
        return sparse.hstack(
            [self.encoder.transform(data[CAT]), self.scaler.transform(data[NUM])],
            format="csr",
        )

    def predict(self, data):
        return self.lr.predict_proba(self.transform(data))[:, 1]


def site_bucket(ids, seed):
    # All three alternatives and any repeated batches for one locus share a fold.
    return (
        pd.util.hash_pandas_object(ids.astype(str) + f":{seed}", index=False).to_numpy()
        % 5
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--seed", type=int, default=20261015)
    a = p.parse_args()
    columns = ["site_id", "chrom", "occupancy", "observed"] + FEATURES
    counts = Counter()
    for part in batches(a.candidates, columns):
        train = part[part.occupancy < 3].copy()
        train["bucket"] = site_bucket(train.site_id, a.seed)
        counts.update(
            train.groupby(["chrom", "bucket"] + FEATURES + ["observed"])
            .size()
            .to_dict()
        )
    if not counts:
        raise SystemExit(
            "Stage 1 needs biallelic/triallelic loci with observed and unrealized alternatives"
        )
    keys = ["chrom", "bucket"] + FEATURES + ["observed"]
    data = pd.DataFrame(
        [(*key, value) for key, value in counts.items()], columns=keys + ["count"]
    )
    chrom_counts = data.groupby("chrom")["count"].sum()
    if len(chrom_counts) > 1:
        group_column = "chrom"
        loads = np.zeros(min(5, len(chrom_counts)), dtype=np.int64)
        fold_map = {}
        for chrom, n in chrom_counts.sort_values(
            ascending=False, kind="stable"
        ).items():
            fold = int(np.argmin(loads))
            fold_map[chrom] = fold
            loads[fold] += n
        data["fold"] = data.chrom.map(fold_map)
    else:
        group_column = "site_id"
        fold_map = None
        data["fold"] = data.bucket
    data = data.groupby(["fold"] + FEATURES + ["observed"], as_index=False)[
        "count"
    ].sum()
    folds = sorted(data.fold.unique())
    if len(folds) < 2:
        raise SystemExit("Stage 1 cannot cross-fit: fewer than two validation groups")
    models = {}
    predictions = np.zeros(len(data))
    try:
        for fold in folds:
            model = CountModel(a.seed).fit(data[data.fold != fold])
            held = data.fold == fold
            predictions[held] = model.predict(data[held])
            models[fold] = model
        full = CountModel(a.seed).fit(data)
    except ValueError as error:
        raise SystemExit(f"Stage 1 cannot cross-fit: {error}")
    auc = float(roc_auc_score(data.observed, predictions, sample_weight=data["count"]))
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    candidate_sink = ParquetSink(out / "stage1_candidates.parquet")
    site_sink = ParquetSink(out / "stage1_sites.parquet")
    try:
        # prepare writes exactly three consecutive rows per locus; 60000 preserves boundaries.
        for part in batches(a.candidates):
            pred = full.predict(part)
            fold_ids = (
                part.chrom.map(fold_map).to_numpy()
                if fold_map
                else site_bucket(part.site_id, a.seed)
            )
            for fold, model in models.items():
                held = (fold_ids == fold) & (part.occupancy.to_numpy() < 3)
                if held.any():
                    pred[held] = model.predict(part.loc[held])
            saved = part[["site_id", "variant_id", "occupancy"]].copy()
            saved["mutation_opportunity"] = pred
            saved["mutation_prediction_source"] = np.where(
                part.occupancy < 3, "out_of_fold", "full_fit_tetra_only"
            )
            candidate_sink.write(saved)
            if not part.groupby("site_id").size().eq(3).all():
                raise ValueError(
                    "Stage 1 requires three consecutive candidates per locus; re-run prepare"
                )
            summary = (
                saved.groupby("site_id", sort=False)
                .agg(
                    mutation_mean=("mutation_opportunity", "mean"),
                    mutation_max=("mutation_opportunity", "max"),
                    mutation_min=("mutation_opportunity", "min"),
                    mutation_sd=("mutation_opportunity", "std"),
                )
                .reset_index()
            )
            site_sink.write(summary)
        candidate_sink.commit()
        site_sink.commit()
    finally:
        candidate_sink.close()
        site_sink.close()
    n_train = int(data["count"].sum())
    metrics = dict(
        cv_auc=auc,
        cv_status="ok",
        cv_group=group_column,
        n_train=n_train,
        n_sites=n_train // 3,
        fit="frequency-weighted exact feature counts; no training subsampling",
        folds=(
            "chromosome blocks"
            if fold_map
            else "stable locus hash into up to five folds"
        ),
        interpretation="relative alternative-specific mutation opportunity conditional on variable loci; O<3 uses held-out scores; class-weighted scores are not calibrated probabilities",
    )
    (out / "stage1_metrics.json").write_text(
        json.dumps(metrics, indent=2, allow_nan=False) + "\n"
    )
    print(
        f"Stage1 CV AUC={auc:.4f} candidates={n_train:,}; aggregated feature patterns={len(data):,}"
    )


if __name__ == "__main__":
    main()
