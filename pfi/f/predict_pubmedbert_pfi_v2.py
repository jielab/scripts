#!/usr/bin/env python3
"""Predict PFI scores for a full-text sample using the no-expert fallback model."""
from __future__ import annotations

import argparse
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


def norm_col(x: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(x).lower())).strip("_")


def first_existing(cols, candidates):
    lookup = {norm_col(c): c for c in cols}
    for c in candidates:
        if norm_col(c) in lookup:
            return lookup[norm_col(c)]
    return None


def dataset_for(path: str | Path) -> ds.Dataset:
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"ERROR: input does not exist: {p}")
    return ds.dataset(str(p), format="parquet", partitioning="hive") if p.is_dir() else ds.dataset(str(p), format="parquet")


def load_model(model_dir: str | Path):
    p = Path(model_dir) / "model.pkl"
    if not p.exists():
        raise SystemExit(f"ERROR: model file missing: {p}. Run 09b first.")
    with open(p, "rb") as f:
        obj = pickle.load(f)
    return obj


def write_parquet_batches(path: str | Path, frames):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    writer = None
    n = 0
    try:
        for df in frames:
            if df is None or df.empty:
                continue
            table = pa.Table.from_pandas(df, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(path, table.schema)
            writer.write_table(table)
            n += len(df)
    finally:
        if writer is not None:
            writer.close()
    return n


def iter_rows_from_dataset(path: str, full_text_only: bool, limit: int, sample_rows: int, batch_size: int, text_max_chars: int):
    d = dataset_for(path)
    cols_all = list(d.schema.names)
    id_cols = []
    for candidates in [["article_id", "id"], ["pmid"], ["pmcid"], ["doi"], ["year", "publication_year", "pub_year"], ["mesh_terms", "mesh"], ["has_pmc_full_text", "has_full_text"]]:
        c = first_existing(cols_all, candidates)
        if c and c not in id_cols:
            id_cols.append(c)
    text_col = first_existing(cols_all, ["model_text", "blinded_text", "full_text", "text", "article_text", "abstract", "title"])
    if text_col is None:
        raise SystemExit(f"ERROR: no text column found in prediction dataset. Available columns={cols_all}")
    scan_cols = list(dict.fromkeys(id_cols + [text_col]))
    full_col = first_existing(scan_cols, ["has_pmc_full_text", "has_full_text"])
    n_out = 0
    for batch in d.scanner(columns=scan_cols, batch_size=batch_size, use_threads=True).to_batches():
        df = batch.to_pandas()
        df.columns = [norm_col(c) for c in df.columns]
        tcol = norm_col(text_col)
        if full_text_only and full_col:
            fcol = norm_col(full_col)
            if fcol in df.columns:
                ok = df[fcol].fillna(0).astype(str).str.lower().isin(["1", "true", "t", "yes", "y"])
                df = df[ok].copy()
        if df.empty:
            continue
        df["model_text"] = df[tcol].fillna("").astype(str).str.slice(0, text_max_chars)
        df = df[df["model_text"].str.len() > 0].copy()
        if df.empty:
            continue
        if "article_id" not in df.columns:
            base = []
            for c in ["pmid", "pmcid", "doi"]:
                if c in df.columns:
                    base.append(c)
            df["article_id"] = df[base].fillna("").astype(str).agg("|".join, axis=1) if base else [f"row_{n_out+i}" for i in range(len(df))]
        need = [c for c in ["article_id", "pmid", "pmcid", "doi", "year", "mesh_terms", "has_pmc_full_text", "model_text"] if c in df.columns]
        df = df[need].copy()
        remain = None
        if limit and limit > 0:
            remain = limit - n_out
        elif sample_rows and sample_rows > 0:
            remain = sample_rows - n_out
        if remain is not None and remain <= 0:
            break
        if remain is not None and len(df) > remain:
            df = df.iloc[:remain].copy()
        n_out += len(df)
        yield df
        if (limit and n_out >= limit) or (sample_rows and sample_rows > 0 and n_out >= sample_rows):
            break


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blinded", required=True)
    ap.add_argument("--model_dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch_size", type=int, default=64)
    ap.add_argument("--max_length", type=int, default=512)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--progress_every", type=int, default=10000)
    ap.add_argument("--cpu_full_warning_rows", type=int, default=1000000)
    ap.add_argument("--sample_rows", type=int, default=0)
    ap.add_argument("--sample_seed", type=int, default=1)
    ap.add_argument("--sample_path", default="")
    ap.add_argument("--fresh_sample", action="store_true")
    ap.add_argument("--full_text_only", action="store_true")
    ap.add_argument("--fp16", action="store_true")
    ap.add_argument("--allow_cpu_full", action="store_true")
    args = ap.parse_args()

    obj = load_model(args.model_dir)
    model = obj["pipeline"]
    target_cols = obj["target_cols"]
    text_max_chars = max(1000, args.max_length * 20)

    sample_path = Path(args.sample_path) if args.sample_path else None
    if sample_path and sample_path.exists() and not args.fresh_sample:
        sample_df = pq.read_table(sample_path).to_pandas()
        batches = (sample_df.iloc[i:i+args.batch_size].copy() for i in range(0, len(sample_df), args.batch_size))
    else:
        # Build sample path while streaming. This is intentionally simple and deterministic: first eligible rows.
        rows = []
        for df in iter_rows_from_dataset(args.blinded, args.full_text_only, args.limit, args.sample_rows, max(args.batch_size*16, 4096), text_max_chars):
            rows.append(df)
        sample_df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        if sample_df.empty:
            raise SystemExit("ERROR: no rows available for prediction sample.")
        if sample_path:
            sample_path.parent.mkdir(parents=True, exist_ok=True)
            pq.write_table(pa.Table.from_pandas(sample_df, preserve_index=False), sample_path)
            print(f"OK: wrote prediction sample: {sample_path} ({len(sample_df):,} rows)")
        batches = (sample_df.iloc[i:i+args.batch_size].copy() for i in range(0, len(sample_df), args.batch_size))

    def pred_frames():
        done = 0
        for df in batches:
            texts = df["model_text"].fillna("").astype(str).str.slice(0, text_max_chars)
            pred = np.clip(model.predict(texts), 0, 100)
            out = df[[c for c in ["article_id", "pmid", "pmcid", "doi", "year", "mesh_terms", "has_pmc_full_text"] if c in df.columns]].copy()
            for i, c in enumerate(target_cols):
                out[f"pred_{c}"] = pred[:, i]
            if "paper_score" in target_cols:
                out["predicted_paper_score"] = out["pred_paper_score"]
            if "expected_journal_score_0_100" in target_cols:
                out["expected_journal_score_0_100"] = out["pred_expected_journal_score_0_100"]
                out["pred_expected_journal_score_0_100"] = out["pred_expected_journal_score_0_100"]
            done += len(out)
            if args.progress_every and done % args.progress_every < len(out):
                print(f"Predicted {done:,} rows")
            yield out

    n = write_parquet_batches(args.out, pred_frames())
    print(f"OK: wrote predictions: {args.out} ({n:,} rows)")


if __name__ == "__main__":
    main()
