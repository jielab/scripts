#!/usr/bin/env python3
"""Compute publication fairness results from PubMedBERT predictions."""


# 🚩 Imports, paths, and shared inputs
from __future__ import annotations

import argparse
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from make_pfi_training_dataset import normalize_rank_file

def sql_string(value: object) -> str:
    return str(value).replace("'", "''")

def parquet_source(path: str) -> str:
    p = Path(path)
    if p.is_dir() or p.suffix != ".parquet":
        return sql_string(p / "**" / "*.parquet")
    return sql_string(p)

FAMOUS_INSTITUTION_PATTERNS = [
    r"\bharvard\b", r"\bstanford\b", r"\boxford\b", r"\bcambridge\b", r"\bmit\b",
    r"\bmassachusetts general\b", r"\bmayo clinic\b", r"\bjohns hopkins\b", r"\bucsf\b",
    r"\byale\b", r"\bcolumbia university\b", r"\bduke university\b", r"\bcornell\b",
    r"\bnih\b", r"\bnational institutes of health\b", r"\bkarolinska\b",
    r"\bweizmann\b", r"\bmax planck\b", r"\bimperial college\b", r"\bucl\b",
]

def zscore_grouped(df: pd.DataFrame, value_col: str, group_cols: list[str], min_group_n: int = 30) -> pd.Series:
    z = pd.Series(np.nan, index=df.index, dtype="float64")
    grouped = df.groupby(group_cols, dropna=False)
    for _, idx in grouped.groups.items():
        idx = list(idx)
        if len(idx) < min_group_n:
            continue
        s = df.loc[idx, value_col].astype(float)
        sd = s.std(ddof=0)
        if sd and np.isfinite(sd):
            z.loc[idx] = (s - s.mean()) / sd
    missing = z.isna()
    if missing.any() and group_cols != ["domain"]:
        z2 = zscore_grouped(df.loc[missing].copy(), value_col, ["domain"], min_group_n=min_group_n)
        z.loc[missing] = z2
    missing = z.isna()
    if missing.any():
        s = df[value_col].astype(float)
        sd = s.std(ddof=0)
        z.loc[missing] = (s.loc[missing] - s.mean()) / sd if sd and np.isfinite(sd) else 0.0
    return z

def fairness_category(z: float) -> str:
    if not np.isfinite(z):
        return "unknown"
    if z >= 2:
        return "strongly_unfair_favored"
    if z >= 1:
        return "modestly_unfair_favored"
    if z <= -2:
        return "strongly_unfair_unfavored"
    if z <= -1:
        return "modestly_unfair_unfavored"
    return "fair"

def institution_proxy(affil: object, extra_patterns: list[str] | None = None) -> int:
    s = str(affil or "").lower()
    if not s or s == "nan":
        return 0
    patterns = FAMOUS_INSTITUTION_PATTERNS + (extra_patterns or [])
    return int(any(re.search(p, s) for p in patterns))

def load_extra_institution_patterns(path: str | None) -> list[str]:
    if not path:
        return []
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(p)
    col = None
    for c in ["institution", "institution_name", "name"]:
        if c in df.columns:
            col = c
            break
    if col is None:
        return []
    score_col = None
    for c in ["prestige_percentile", "rank_percentile", "score", "world_rank"]:
        if c in df.columns:
            score_col = c
            break
    if score_col and "rank" in score_col:
        df = df[pd.to_numeric(df[score_col], errors="coerce") <= 100]
    elif score_col:
        df = df[pd.to_numeric(df[score_col], errors="coerce") >= 90]
    names = df[col].dropna().astype(str).head(500).tolist()
    return [re.escape(x.lower()) for x in names if len(x) >= 4]

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--journal_rank", required=True)
    ap.add_argument("--rank_normalized", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--out_parquet", required=True)
    ap.add_argument("--summary_dir", required=True)
    ap.add_argument("--institution_rank", default=None)
    ap.add_argument("--max_csv_rows", type=int, default=0, help="0 writes all rows; set for a smaller CSV audit table")
    ap.add_argument("--full_text_only", action="store_true", help="Restrict fairness outputs to rows with PMC full text")
    args = ap.parse_args()

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_parquet).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_dir).mkdir(parents=True, exist_ok=True)
    rank_path = Path(args.rank_normalized)
    normalize_rank_file(args.journal_rank, rank_path)

    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    full_text_predicate = """
      AND CASE WHEN coalesce(m.full_text, m.intro_text, m.methods_text, m.results_text, m.discussion_text,
                             m.data_availability_text, m.ethics_text, m.trial_registration_text) IS NOT NULL
               THEN 1 ELSE 0 END = 1
    """ if args.full_text_only else ""
    q = f"""
    WITH rank AS (SELECT * FROM read_parquet('{rank_path.as_posix()}')),
    m AS (
      SELECT *,
        trim(regexp_replace(regexp_replace(lower(coalesce(journal, '')), '[^a-z0-9]+', ' ', 'g'), '\\s+', ' ', 'g')) AS journal_key
      FROM read_parquet('{parquet_source(args.master)}', union_by_name=true)
    ),
    p AS (
      SELECT *, coalesce(pmid, pmcid, doi) AS article_key
      FROM read_parquet('{parquet_source(args.predictions)}', union_by_name=true)
    )
    SELECT
      coalesce(m.pmid, p.pmid) AS pmid,
      coalesce(m.pmcid, p.pmcid) AS pmcid,
      coalesce(m.doi, p.doi) AS doi,
      coalesce(m.year, p.year) AS year,
      m.title, m.journal, rank.jcr_category, rank.domain, rank.journal_tier,
      rank.journal_score_0_100 AS actual_journal_score_0_100,
      p.expected_journal_score_0_100,
      p.predicted_paper_score,
      m.publication_types, m.mesh_terms, m.author_count, m.first_author, m.last_author,
      m.affiliations,
      CASE WHEN coalesce(m.full_text, m.intro_text, m.methods_text, m.results_text, m.discussion_text,
                         m.data_availability_text, m.ethics_text, m.trial_registration_text) IS NOT NULL
           THEN 1 ELSE 0 END AS has_pmc_full_text
    FROM p
    LEFT JOIN m ON coalesce(m.pmid, m.pmcid, m.doi) = coalesce(p.pmid, p.pmcid, p.doi)
    LEFT JOIN rank ON m.journal_key = rank.journal_key
    WHERE rank.journal_score_0_100 IS NOT NULL
    {full_text_predicate}
    """
    df = con.execute(q).fetchdf()
    con.close()
    print(f"Rows with predictions and journal ranks: {len(df):,}")
    if df.empty:
        raise RuntimeError("No prediction rows matched journal ranks.")

    df["field_proxy"] = df["mesh_terms"].fillna("Unknown").astype(str).str.split("|").str[0].replace("", "Unknown")
    df["publication_type_proxy"] = df["publication_types"].fillna("Unknown").astype(str).str.split("|").str[0].replace("", "Unknown")
    for c in ["actual_journal_score_0_100", "expected_journal_score_0_100", "predicted_paper_score"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").clip(0, 100)

    df["fairness_index"] = df["actual_journal_score_0_100"] - df["expected_journal_score_0_100"]
    df["paper_journal_gap"] = df["actual_journal_score_0_100"] - df["predicted_paper_score"]
    df["fairness_z"] = zscore_grouped(df, "fairness_index", ["year", "domain", "publication_type_proxy"])
    df["fairness_category"] = df["fairness_z"].map(fairness_category)
    df["favored_flag"] = df["fairness_category"].str.endswith("favored").astype(int)
    df["unfavored_flag"] = df["fairness_category"].str.endswith("unfavored").astype(int)
    df["strongly_favored_flag"] = df["fairness_category"].eq("strongly_unfair_favored").astype(int)
    df["strongly_unfavored_flag"] = df["fairness_category"].eq("strongly_unfair_unfavored").astype(int)

    extra_patterns = load_extra_institution_patterns(args.institution_rank)
    df["famous_institution_proxy"] = df["affiliations"].map(lambda x: institution_proxy(x, extra_patterns))

    keep = [
        "pmid", "pmcid", "doi", "year", "title", "journal", "jcr_category", "domain", "journal_tier",
        "actual_journal_score_0_100", "expected_journal_score_0_100", "predicted_paper_score",
        "fairness_index", "fairness_z", "paper_journal_gap", "fairness_category", "favored_flag", "unfavored_flag",
        "strongly_favored_flag", "strongly_unfavored_flag", "famous_institution_proxy",
        "author_count", "first_author", "last_author", "has_pmc_full_text", "field_proxy", "publication_type_proxy",
    ]
    keep = [c for c in keep if c in df.columns]
    out_df = df[keep].copy()
    out_df.to_parquet(args.out_parquet, index=False)

    csv_df = out_df
    if args.max_csv_rows and args.max_csv_rows > 0 and len(csv_df) > args.max_csv_rows:
        # Keep extreme papers plus a random background sample for manual audit.
        n_extreme = min(args.max_csv_rows // 2, len(csv_df))
        extreme = pd.concat([
            csv_df.nlargest(n_extreme // 2, "fairness_z"),
            csv_df.nsmallest(n_extreme // 2, "fairness_z"),
        ]).drop_duplicates()
        background = csv_df.drop(index=extreme.index, errors="ignore").sample(
            max(args.max_csv_rows - len(extreme), 0), random_state=1
        )
        csv_df = pd.concat([extreme, background]).drop_duplicates()
    csv_df.to_csv(args.out_csv, index=False)

    summary_dir = Path(args.summary_dir)
    out_df.groupby("fairness_category", dropna=False).size().reset_index(name="n").to_csv(summary_dir / "fairness_category_counts.csv", index=False)
    out_df.groupby(["year", "domain", "fairness_category"], dropna=False).size().reset_index(name="n").to_csv(summary_dir / "fairness_category_by_year_domain.csv", index=False)
    out_df.groupby(["domain"], dropna=False).agg(
        n=("fairness_index", "size"),
        mean_actual_journal_score=("actual_journal_score_0_100", "mean"),
        mean_expected_journal_score=("expected_journal_score_0_100", "mean"),
        mean_predicted_paper_score=("predicted_paper_score", "mean"),
        mean_fairness_index=("fairness_index", "mean"),
        sd_fairness_index=("fairness_index", "std"),
        strong_favored_n=("strongly_favored_flag", "sum"),
        strong_unfavored_n=("strongly_unfavored_flag", "sum"),
        famous_institution_rate=("famous_institution_proxy", "mean"),
    ).reset_index().to_csv(summary_dir / "fairness_summary_by_domain.csv", index=False)

    strong = out_df[out_df["fairness_category"].eq("strongly_unfair_favored")]
    if not strong.empty:
        strong.sort_values("fairness_z", ascending=False).head(500).to_csv(summary_dir / "top_strongly_favored_candidates.csv", index=False)
    weak = out_df[out_df["fairness_category"].eq("strongly_unfair_unfavored")]
    if not weak.empty:
        weak.sort_values("fairness_z", ascending=True).head(500).to_csv(summary_dir / "top_strongly_unfavored_candidates.csv", index=False)

    print(f"Wrote fairness parquet: {args.out_parquet}")
    print(f"Wrote fairness CSV: {args.out_csv}")
    print(f"Wrote summaries: {summary_dir}")

if __name__ == "__main__":
    main()
