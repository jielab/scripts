#!/usr/bin/env python3
"""Create the PFI training dataset in no-expert/weak-label mode.

The script is schema-tolerant and can recover journal metadata from
articles_master when the review key has only PMID/PMCID identifiers.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq
except Exception as e:
    raise SystemExit(f"ERROR: pyarrow is required: {e}")


def norm_col(x: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(x).lower())).strip("_")


def clean_title_value(x) -> str:
    s = str(x or "").lower().strip()
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[^a-z0-9 ]+", "", s)
    s = re.sub(r"\s+", "_", s).strip("_")
    return s


def clean_issn(x) -> str:
    s = re.sub(r"[^0-9Xx]", "", str(x or "")).upper()
    if len(s) == 8:
        return s[:4] + "-" + s[4:]
    return ""


def make_journal_key_row(row: pd.Series) -> str:
    for c in ["journal_key", "source_key", "journal_id"]:
        if c in row and str(row.get(c, "")).strip() not in {"", "nan", "None"}:
            return str(row[c]).strip().lower()
    for c, prefix in [("eissn", "eissn:"), ("issn", "issn:")]:
        if c in row:
            v = clean_issn(row[c])
            if v:
                return prefix + v
    for c in ["journal_title", "journal", "journal_name", "source_title", "medline_ta"]:
        if c in row:
            v = clean_title_value(row[c])
            if v:
                return "title:" + v
    return ""


def first_existing(cols, candidates):
    lookup = {norm_col(c): c for c in cols}
    for c in candidates:
        if norm_col(c) in lookup:
            return lookup[norm_col(c)]
    return None


def read_any(path: str | Path, columns: list[str] | None = None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    if p.is_dir() or p.suffix.lower() == ".parquet":
        d = ds.dataset(str(p), format="parquet", partitioning="hive")
        avail = list(d.schema.names)
        use = [c for c in (columns or avail) if c in avail]
        if not use:
            return pd.DataFrame()
        return d.to_table(columns=use).to_pandas()
    return pd.read_csv(p, dtype=str, low_memory=False)


def to_num(x, default=np.nan):
    return pd.to_numeric(x, errors="coerce").fillna(default)


def lookup_master(master_path: str, ids: pd.DataFrame) -> pd.DataFrame:
    if not master_path or not Path(master_path).exists() or ids.empty:
        return pd.DataFrame()
    d = ds.dataset(str(master_path), format="parquet", partitioning="hive")
    cols = list(d.schema.names)
    want_names = [
        "article_id", "pmid", "pmcid", "doi", "year", "publication_year", "pub_year",
        "journal_key", "source_key", "journal_id", "journal", "journal_title", "journal_name", "source_title", "medline_ta",
        "issn", "eissn", "mesh_terms", "mesh", "publication_type", "publication_type_proxy", "has_pmc_full_text",
        "authors", "author_names", "affiliations", "institutions", "first_author", "last_author",
    ]
    use = []
    for nm in want_names:
        hit = first_existing(cols, [nm])
        if hit and hit not in use:
            use.append(hit)
    join_cols = [c for c in [first_existing(use, ["article_id"]), first_existing(use, ["pmid"]), first_existing(use, ["pmcid"]), first_existing(use, ["doi"])] if c]
    if not join_cols:
        return pd.DataFrame()
    id_sets = {}
    for c in ["article_id", "pmid", "pmcid", "doi"]:
        if c in ids.columns:
            vals = set(ids[c].dropna().astype(str)) - {"", "nan", "None"}
            if vals:
                id_sets[c] = vals
    out = []
    for batch in d.scanner(columns=use, batch_size=65536, use_threads=True).to_batches():
        df = batch.to_pandas()
        df.columns = [norm_col(c) for c in df.columns]
        mask = pd.Series(False, index=df.index)
        for c, vals in id_sets.items():
            if c in df.columns:
                mask = mask | df[c].astype(str).isin(vals)
        if mask.any():
            out.append(df.loc[mask].copy())
    if not out:
        return pd.DataFrame()
    ans = pd.concat(out, ignore_index=True).drop_duplicates()
    return ans


def merge_by_any(left: pd.DataFrame, right: pd.DataFrame, suffix: str) -> pd.DataFrame:
    if right.empty:
        return left
    left = left.copy(); right = right.copy()
    left.columns = [norm_col(c) for c in left.columns]
    right.columns = [norm_col(c) for c in right.columns]
    for key in ["review_id", "article_id", "pmid", "pmcid", "doi"]:
        if key in left.columns and key in right.columns:
            r = right.drop_duplicates(key)
            rename = {}
            used = set(left.columns)
            keep = [key]
            for c in r.columns:
                if c == key:
                    continue
                new = c
                if new in used:
                    new = f"{c}{suffix}"
                    i = 2
                    while new in used or new in rename.values():
                        new = f"{c}{suffix}{i}"
                        i += 1
                rename[c] = new
                keep.append(c)
            if len(keep) == 1:
                return left
            return left.merge(r[keep].rename(columns=rename), on=key, how="left")
    return left


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--review_scores", required=True)
    ap.add_argument("--review_key", default="")
    ap.add_argument("--blinded", default="")
    ap.add_argument("--features", default="")
    ap.add_argument("--journal_rank", required=True)
    ap.add_argument("--master", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary_out", required=True)
    ap.add_argument("--fit_summary_out", required=True)
    ap.add_argument("--audit_csv_out", required=True)
    ap.add_argument("--max_rows", type=int, default=1000)
    ap.add_argument("--min_text_chars", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    review = read_any(args.review_scores)
    if review.empty:
        raise SystemExit(f"ERROR: review score file is empty or missing: {args.review_scores}")
    review.columns = [norm_col(c) for c in review.columns]

    key = read_any(args.review_key) if args.review_key else pd.DataFrame()
    if not key.empty:
        key.columns = [norm_col(c) for c in key.columns]
        review = merge_by_any(review, key, "_key")

    # If journal metadata are absent in key/review, recover them from master.
    id_cols = [c for c in ["review_id", "article_id", "pmid", "pmcid", "doi"] if c in review.columns]
    master_meta = lookup_master(args.master, review[id_cols].copy()) if args.master and id_cols else pd.DataFrame()
    if not master_meta.empty:
        review = merge_by_any(review, master_meta, "_master")

    # Coalesce duplicate columns created by merges.
    for base in ["year", "journal_key", "journal", "journal_title", "issn", "eissn", "mesh_terms", "domain", "publication_type_proxy", "has_pmc_full_text"]:
        alts = [c for c in review.columns if c == base or c.startswith(base + "_")]
        if alts:
            val = review[alts[0]].copy()
            for c in alts[1:]:
                val = val.where(val.astype(str).str.strip().ne("") & val.notna(), review[c])
            review[base] = val

    # Text column: prefer blinded/model text from review file; fallback to any text-like field.
    text_col = first_existing(review.columns, ["model_text", "blinded_text", "full_text", "text", "abstract", "title"])
    if text_col is None:
        raise SystemExit(f"ERROR: training data has no text column. Columns: {list(review.columns)}")
    review["model_text"] = review[text_col].fillna("").astype(str)
    review = review[review["model_text"].str.len() >= args.min_text_chars].copy()
    if review.empty:
        raise SystemExit("ERROR: no review rows passed min_text_chars after text-column detection.")

    if "article_id" not in review.columns:
        review["article_id"] = ""
    for c in ["article_id", "pmid", "pmcid", "doi"]:
        if c not in review.columns:
            review[c] = ""
    missing_aid = review["article_id"].astype(str).isin(["", "nan", "None"])
    review.loc[missing_aid, "article_id"] = review.loc[missing_aid, [c for c in ["pmid", "pmcid", "doi"] if c in review.columns]].fillna("").astype(str).agg("|".join, axis=1)

    if "year" not in review.columns:
        review["year"] = np.nan
    review["year"] = pd.to_numeric(review["year"], errors="coerce").fillna(2000).astype(int)
    if "journal_key" not in review.columns:
        review["journal_key"] = ""
    miss_jkey = review["journal_key"].fillna("").astype(str).str.len() == 0
    if miss_jkey.any():
        review.loc[miss_jkey, "journal_key"] = review.loc[miss_jkey].apply(make_journal_key_row, axis=1)

    rank = read_any(args.journal_rank)
    if rank.empty:
        raise SystemExit(f"ERROR: normalized journal rank is missing/empty: {args.journal_rank}")
    rank.columns = [norm_col(c) for c in rank.columns]
    if "jcr_year" not in rank.columns:
        rank["jcr_year"] = review["year"].median() if len(review) else 2000
    if "actual_journal_score_0_100" not in rank.columns:
        if "jif_percentile" in rank.columns:
            rank["actual_journal_score_0_100"] = rank["jif_percentile"]
        else:
            rank["actual_journal_score_0_100"] = 50.0
    rank["jcr_year"] = pd.to_numeric(rank["jcr_year"], errors="coerce").fillna(2000).astype(int)
    rank["actual_journal_score_0_100"] = to_num(rank["actual_journal_score_0_100"], 50).clip(0, 100)
    r1 = rank[[c for c in ["journal_key", "jcr_year", "jcr_category", "journal_domain", "actual_journal_score_0_100", "jif_quartile", "jif_percentile"] if c in rank.columns]].copy()
    r1 = r1.dropna(subset=["journal_key"]).drop_duplicates(["journal_key", "jcr_year"])
    review = review.merge(r1, left_on=["journal_key", "year"], right_on=["journal_key", "jcr_year"], how="left", suffixes=("", "_rank"))
    missing_score = review["actual_journal_score_0_100"].isna() if "actual_journal_score_0_100" in review.columns else pd.Series(True, index=review.index)
    if missing_score.any():
        r2 = rank.groupby("journal_key", as_index=False).agg(
            actual_journal_score_fallback=("actual_journal_score_0_100", "mean"),
            journal_domain_fallback=("journal_domain", lambda x: x.dropna().astype(str).iloc[0] if len(x.dropna()) else "unknown") if "journal_domain" in rank.columns else ("actual_journal_score_0_100", lambda x: "unknown"),
        )
        review = review.merge(r2, on="journal_key", how="left")
        if "actual_journal_score_0_100" not in review.columns:
            review["actual_journal_score_0_100"] = np.nan
        review.loc[review["actual_journal_score_0_100"].isna(), "actual_journal_score_0_100"] = review.loc[review["actual_journal_score_0_100"].isna(), "actual_journal_score_fallback"]
    review["actual_journal_score_0_100"] = to_num(review["actual_journal_score_0_100"], 50).clip(0, 100)
    review["expected_journal_score_0_100"] = review["actual_journal_score_0_100"]

    label_cols = ["new_results", "new_methods", "new_data", "new_software", "clinical_public_health_value", "evidence_strength", "reproducibility", "claim_evidence_gap", "paper_score"]
    for c in label_cols:
        if c not in review.columns:
            review[c] = 50.0
        review[c] = to_num(review[c], 50).clip(0, 100)

    if "domain" not in review.columns or review["domain"].fillna("").astype(str).eq("").all():
        jd = review["journal_domain"] if "journal_domain" in review.columns else review.get("journal_domain_fallback", "unknown")
        review["domain"] = pd.Series(jd).fillna("unknown").astype(str)
    if "publication_type_proxy" not in review.columns:
        review["publication_type_proxy"] = "research_article"

    keep = [
        "review_id", "article_id", "pmid", "pmcid", "doi", "year", "domain", "publication_type_proxy", "journal_key",
        "model_text", "actual_journal_score_0_100", "expected_journal_score_0_100",
    ] + label_cols
    for c in keep:
        if c not in review.columns:
            review[c] = ""
    out = review[keep].drop_duplicates("article_id").copy()
    if args.max_rows and len(out) > args.max_rows:
        out = out.sample(n=args.max_rows, random_state=args.seed).reset_index(drop=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), args.out)
    out.head(200).to_csv(args.audit_csv_out, index=False)

    def corr(a, b):
        a = pd.to_numeric(a, errors="coerce"); b = pd.to_numeric(b, errors="coerce")
        if a.notna().sum() < 3 or b.notna().sum() < 3:
            return np.nan
        return float(a.corr(b))
    summary = pd.DataFrame([{
        "n_training_rows": len(out),
        "mean_paper_score": float(out["paper_score"].mean()),
        "mean_actual_journal_score_0_100": float(out["actual_journal_score_0_100"].mean()),
        "cor_paper_score_journal_score": corr(out["paper_score"], out["actual_journal_score_0_100"]),
        "n_domains": int(out["domain"].nunique()),
        "label_source": "manual_or_pseudo_from_review_scores_file",
    }])
    summary.to_csv(args.summary_out, index=False)
    summary.to_csv(args.fit_summary_out, index=False)
    print(f"OK: wrote training dataset: {args.out} ({len(out):,} rows)")
    print(f"OK: wrote training summary: {args.summary_out}")


if __name__ == "__main__":
    main()
