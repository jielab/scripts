#!/usr/bin/env python3
"""Train a lightweight no-expert PFI text model.

The shell script name is PubMedBERT-compatible, but this fallback uses a
transparent TF-IDF + ridge multi-output regressor so the pipeline can run without
GPU/BERT fine-tuning. It writes the same key output files used by downstream
steps. You can later replace this with the real PubMedBERT trainer without
changing the after-08 workflow.
"""
from __future__ import annotations

import argparse
import json
import math
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputRegressor
from sklearn.pipeline import Pipeline


def read_parquet(path: str) -> pd.DataFrame:
    return pq.read_table(path).to_pandas()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--training", required=True)
    ap.add_argument("--model_name", default="")
    ap.add_argument("--fallback_model_name", default="")
    ap.add_argument("--model_dir", required=True)
    ap.add_argument("--metrics_out", required=True)
    ap.add_argument("--valid_pred_out", required=True)
    ap.add_argument("--fit_summary_out", required=True)
    ap.add_argument("--epochs", type=int, default=3)
    ap.add_argument("--batch_size", type=int, default=8)
    ap.add_argument("--grad_accum", type=int, default=1)
    ap.add_argument("--max_length", type=int, default=512)
    ap.add_argument("--lr", type=float, default=2e-5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--fp16", action="store_true")
    args = ap.parse_args()

    df = read_parquet(args.training)
    if df.empty:
        raise SystemExit(f"ERROR: empty training dataset: {args.training}")
    text_col = "model_text" if "model_text" in df.columns else "blinded_text" if "blinded_text" in df.columns else None
    if text_col is None:
        raise SystemExit(f"ERROR: no model_text/blinded_text column in training data. Columns={list(df.columns)}")
    target_cols = [c for c in [
        "paper_score", "new_results", "new_methods", "new_data", "new_software",
        "clinical_public_health_value", "evidence_strength", "reproducibility", "claim_evidence_gap",
        "expected_journal_score_0_100",
    ] if c in df.columns]
    if "expected_journal_score_0_100" not in target_cols and "actual_journal_score_0_100" in df.columns:
        df["expected_journal_score_0_100"] = pd.to_numeric(df["actual_journal_score_0_100"], errors="coerce")
        target_cols.append("expected_journal_score_0_100")
    if not target_cols:
        raise SystemExit("ERROR: no numeric target columns found for training.")

    X = df[text_col].fillna("").astype(str).str.slice(0, max(1000, args.max_length * 20))
    Y = df[target_cols].apply(pd.to_numeric, errors="coerce").fillna(50).clip(0, 100)

    if len(df) >= 10:
        idx_train, idx_valid = train_test_split(np.arange(len(df)), test_size=max(0.2, min(0.4, 20/len(df))), random_state=args.seed)
    else:
        idx_train = np.arange(len(df)); idx_valid = np.arange(len(df))

    max_features = int(max(2000, min(50000, len(df) * 200)))
    model = Pipeline([
        ("tfidf", TfidfVectorizer(lowercase=True, max_features=max_features, ngram_range=(1, 2), min_df=1, max_df=0.98)),
        ("regressor", MultiOutputRegressor(Ridge(alpha=5.0, random_state=args.seed))),
    ])
    model.fit(X.iloc[idx_train], Y.iloc[idx_train])

    pred = np.clip(model.predict(X.iloc[idx_valid]), 0, 100)
    yv = Y.iloc[idx_valid].to_numpy()
    metrics = {
        "trainer": "tfidf_ridge_noexpert_fallback",
        "note": "This is a runnable weak-supervision fallback, not a final PubMedBERT model.",
        "n_rows": int(len(df)),
        "n_train": int(len(idx_train)),
        "n_valid": int(len(idx_valid)),
        "target_cols": target_cols,
        "text_col": text_col,
    }
    rows = []
    for i, c in enumerate(target_cols):
        rmse = float(math.sqrt(mean_squared_error(yv[:, i], pred[:, i]))) if len(yv) else float("nan")
        mae = float(mean_absolute_error(yv[:, i], pred[:, i])) if len(yv) else float("nan")
        r2 = float(r2_score(yv[:, i], pred[:, i])) if len(yv) >= 3 else float("nan")
        metrics[c] = {"rmse": rmse, "mae": mae, "r2": r2}
        rows.append({"target": c, "rmse": rmse, "mae": mae, "r2": r2})

    valid = df.iloc[idx_valid][[c for c in ["review_id", "article_id", "pmid", "pmcid", "year", "domain", "publication_type_proxy"] if c in df.columns]].copy()
    for i, c in enumerate(target_cols):
        valid[c] = yv[:, i]
        valid[f"pred_{c}"] = pred[:, i]
    if "expected_journal_score_0_100" in target_cols:
        valid["pred_expected_journal_score_0_100"] = valid["pred_expected_journal_score_0_100"] if "pred_expected_journal_score_0_100" in valid.columns else valid[f"pred_expected_journal_score_0_100"]
    if "paper_score" in target_cols:
        valid["pred_paper_score"] = valid[f"pred_paper_score"] if "pred_paper_score" in valid.columns else valid[f"pred_paper_score"]

    model_dir = Path(args.model_dir); model_dir.mkdir(parents=True, exist_ok=True)
    with open(model_dir / "model.pkl", "wb") as f:
        pickle.dump({"pipeline": model, "target_cols": target_cols, "text_col": text_col, "max_length": args.max_length}, f)
    with open(model_dir / "model_info.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)

    Path(args.metrics_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.metrics_out, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
    valid.to_csv(args.valid_pred_out, index=False)
    pd.DataFrame(rows).to_csv(args.fit_summary_out, index=False)
    print(f"OK: trained no-expert fallback model and wrote: {model_dir}")
    print(f"OK: wrote validation predictions: {args.valid_pred_out}")


if __name__ == "__main__":
    main()
