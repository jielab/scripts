#!/usr/bin/env python3
"""Create provisional review inputs for bootstrapping the PFI training flow."""


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

def sql_string(value: object) -> str:
    return str(value).replace("'", "''")

def parquet_source(path: str) -> str:
    p = Path(path)
    if p.is_dir() or p.suffix != ".parquet":
        return sql_string(p / "**" / "*.parquet")
    return sql_string(p)

def norm_journal(x: object) -> str | None:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    s = re.sub(r"[^a-z0-9]+", " ", str(x).lower())
    s = re.sub(r"\s+", " ", s).strip()
    return s or None

def norm_id(x: object) -> str | None:
    if x is None or pd.isna(x):
        return None
    s = str(x).strip()
    if not s or s.lower() in {"nan", "none", "<na>"}:
        return None
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    return s

def article_key_from(df: pd.DataFrame) -> pd.Series:
    return (
        df["pmid"].map(norm_id)
        .combine_first(df["pmcid"].map(norm_id))
        .combine_first(df["doi"].map(norm_id))
    )

def infer_domain_from_journal(journal: object) -> str:
    s = str(journal or "").lower()
    if any(k in s for k in ["public health", "epidemiol", "environment", "occupational", "health policy"]):
        return "public_health"
    if any(k in s for k in ["clinical", "medicine", "surgery", "oncology", "cardiol", "neurol", "pediatric", "infectious"]):
        return "clinical_science"
    if any(k in s for k in ["cell", "molecular", "genetic", "genomic", "immunol", "microbiol", "neuroscience", "biochem"]):
        return "basic_medical_research"
    if any(k in s for k in ["method", "informatics", "statistics", "bioinformatics", "data"]):
        return "methods_resource"
    return "other_biomedicine"

def clip5(x: float) -> int:
    return int(max(0, min(5, round(float(x)))))

def add_provisional_scores(df: pd.DataFrame) -> pd.DataFrame:
    if "word_count" not in df.columns:
        text = df.get("review_text_preview", pd.Series("", index=df.index)).fillna("").astype(str)
        df["word_count"] = text.str.split().str.len()
    if "sample_size_max" not in df.columns:
        df["sample_size_max"] = 0
    text = df.get("review_text_preview", pd.Series("", index=df.index)).fillna("").astype(str).str.lower()

    text_flags = {
        "has_new_discovery_language": r"\b(?:novel|new|first|discover|identified|demonstrate|showed|revealed)\b",
        "has_replication": r"\b(?:replication|replicated|validated|validation|reproduced)\b",
        "has_method_innovation": r"\b(?:method|algorithm|model|pipeline|tool|protocol|assay|platform|sequencing|framework)\b",
        "has_external_validation": r"\b(?:external validation|independent cohort|validation cohort|held-out|test set)\b",
        "has_sensitivity_analysis": r"\b(?:sensitivity analysis|robustness|subgroup analysis|stratified analysis)\b",
        "has_trial_registration": r"\b(?:clinicaltrials\.gov|trial registration|registered trial)\b",
        "has_ethics": r"\b(?:ethics approval|institutional review board|informed consent)\b",
        "has_data_availability": r"\b(?:data availability|available from|deposited|repository|accession)\b",
        "has_randomized": r"\b(?:randomized|randomised|randomization|randomisation|clinical trial)\b",
        "has_meta_analysis": r"\b(?:meta-analysis|systematic review|pooled analysis)\b",
        "has_large_resource_language": r"\b(?:biobank|database|registry|atlas|genome-wide|nationwide|multi-center|multicenter|large-scale)\b",
        "has_open_science_language": r"\b(?:open source|code availability|github|software|publicly available)\b",
        "has_code_availability": r"\b(?:code availability|github|source code|software package)\b",
    }
    for c, pattern in text_flags.items():
        if c not in df.columns or pd.to_numeric(df[c], errors="coerce").fillna(0).sum() == 0:
            df[c] = text.str.contains(pattern, regex=True, na=False).astype(int)

    for c in df.columns:
        if c.startswith("has_"):
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).clip(0, 1)
    wc = pd.to_numeric(df.get("word_count"), errors="coerce").fillna(0)
    ss = pd.to_numeric(df.get("sample_size_max"), errors="coerce").fillna(0)
    ss_log = np.log1p(ss)
    ss_scaled = 5 * (ss_log - ss_log.min()) / (ss_log.max() - ss_log.min() + 1e-9)

    df["new_results"] = [
        clip5(1 + 2.5 * a + 1.0 * b + 0.5 * c)
        for a, b, c in zip(df.get("has_new_discovery_language", 0), df.get("has_replication", 0), wc.ge(1000).astype(int))
    ]
    df["new_methods"] = [
        clip5(1 + 2.2 * a + 1.2 * b + 0.8 * c)
        for a, b, c in zip(df.get("has_method_innovation", 0), df.get("has_external_validation", 0), df.get("has_sensitivity_analysis", 0))
    ]
    df["new_guidelines"] = [
        clip5(0.5 + 1.5 * a + 1.0 * b + 1.0 * c)
        for a, b, c in zip(df.get("has_trial_registration", 0), df.get("has_ethics", 0), df.get("has_data_availability", 0))
    ]
    df["big_impact"] = [
        clip5(1 + 1.5 * a + 1.0 * b + 0.5 * c)
        for a, b, c in zip(df.get("has_randomized", 0), df.get("has_meta_analysis", 0), ss_scaled)
    ]
    df["big_use"] = [
        clip5(1 + 1.6 * a + 0.9 * b + 0.8 * c)
        for a, b, c in zip(df.get("has_large_resource_language", 0), df.get("has_open_science_language", 0), df.get("has_code_availability", 0))
    ]
    df["claim_evidence_gap"] = [
        clip5(2.5 - 0.6 * a - 0.5 * b - 0.4 * c)
        for a, b, c in zip(df.get("has_replication", 0), df.get("has_external_validation", 0), df.get("has_sensitivity_analysis", 0))
    ]

    score = (
        0.25 * df["new_results"] / 5 * 100
        + 0.25 * df["new_methods"] / 5 * 100
        + 0.10 * df["new_guidelines"] / 5 * 100
        + 0.20 * df["big_impact"] / 5 * 100
        + 0.20 * df["big_use"] / 5 * 100
        - 0.15 * df["claim_evidence_gap"] / 5 * 100
    )
    df["paper_score"] = score.clip(0, 100).round(1)
    df["score_source_note"] = "provisional_auto_score_needs_manual_review"
    return df

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--blinded", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--paper_score_csv", required=True)
    ap.add_argument("--journal_rank_demo", required=True)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min_text_chars", type=int, default=300)
    args = ap.parse_args()

    Path(args.paper_score_csv).parent.mkdir(parents=True, exist_ok=True)
    Path(args.journal_rank_demo).parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    try:
        con.execute(f"SELECT setseed({float((args.seed % 10000) / 10000):.6f})")
    except duckdb.Error:
        pass
    sample_limit = max(int(args.n) * 4, int(args.n), 1)
    query = f"""
    WITH candidates AS (
      SELECT
        pmid, pmcid, doi, year, title, journal, author_count, publication_types, mesh_terms,
        CASE WHEN coalesce(full_text, intro_text, methods_text, results_text, discussion_text,
                           data_availability_text, ethics_text, trial_registration_text) IS NOT NULL
             THEN 1 ELSE 0 END AS has_pmc_full_text,
        coalesce(results_text, methods_text, discussion_text, intro_text, full_text, title, '') AS review_text_preview
      FROM read_parquet('{parquet_source(args.master)}', union_by_name=true)
      WHERE journal IS NOT NULL
        AND length(trim(journal)) > 0
        AND coalesce(full_text, intro_text, methods_text, results_text, discussion_text) IS NOT NULL
    ), sample AS (
      SELECT
        *,
        length(review_text_preview) AS text_chars
      FROM candidates
      WHERE length(review_text_preview) >= {int(args.min_text_chars)}
      ORDER BY random()
      LIMIT {sample_limit}
    )
    SELECT
      *
    FROM sample
    """
    df = con.execute(query).fetchdf()
    con.close()

    df["article_key"] = article_key_from(df)
    df = df.dropna(subset=["article_key"]).drop_duplicates("article_key")
    if len(df) > args.n:
        df = df.sample(args.n, random_state=args.seed)
    df = add_provisional_scores(df)
    df["domain_guess"] = df["journal"].map(infer_domain_from_journal)
    df["blinded_text_preview"] = (
        df["review_text_preview"]
        .fillna("")
        .astype(str)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
        .str.slice(0, 1800)
    )

    paper_cols = [
        "pmid", "pmcid", "doi", "article_key", "year", "title", "journal", "domain_guess",
        "has_pmc_full_text", "word_count", "sample_size_max",
        "new_results", "new_methods", "new_guidelines", "big_impact", "big_use", "claim_evidence_gap",
        "paper_score", "score_source_note", "blinded_text_preview",
    ]
    paper_cols = [c for c in paper_cols if c in df.columns]
    df[paper_cols].to_csv(args.paper_score_csv, index=False, na_rep="")

    journals = df["journal"].dropna().astype(str).value_counts().reset_index()
    journals.columns = ["journal", "n"]
    journals["journal_key"] = journals["journal"].map(norm_journal)
    journals = journals.dropna(subset=["journal_key"]).drop_duplicates("journal_key").reset_index(drop=True)
    journals["jcr_year"] = 2024
    journals["jcr_category"] = journals["journal"].map(lambda x: infer_domain_from_journal(x).replace("_", " ").title())
    journals["domain"] = journals["journal"].map(infer_domain_from_journal)
    # Demo-only synthetic percentile: frequent journals get a broad spread so the pipeline can run.
    ranks = journals["n"].rank(method="average", pct=True)
    journals["jcr_percentile"] = (20 + 75 * ranks).round(1)
    journals["impact_factor"] = (1 + 20 * ranks).round(2)
    journals["jcr_quartile"] = pd.cut(journals["jcr_percentile"], [-1, 25, 50, 75, 101], labels=["Q4", "Q3", "Q2", "Q1"]).astype(str)
    journals["journal_tier"] = journals["jcr_quartile"]
    journals["source"] = "DEMO_SYNTHETIC_REPLACE_WITH_REAL_JCR"
    journals["issn"] = ""
    journals["eissn"] = ""
    journals["journals_in_category"] = journals.groupby("jcr_category")["journal"].transform("count")
    journals["jif_rank"] = journals.groupby("jcr_category")["jcr_percentile"].rank(method="first", ascending=False).astype(int)
    journals[[
        "journal", "issn", "eissn", "jcr_year", "jcr_category", "domain",
        "impact_factor", "jif_rank", "journals_in_category", "jcr_percentile",
        "jcr_quartile", "journal_tier", "source"
    ]].to_csv(args.journal_rank_demo, index=False, na_rep="")

    print(f"Wrote paper-score review CSV: {args.paper_score_csv} ({len(df):,} rows)")
    print(f"Wrote demo journal-rank CSV: {args.journal_rank_demo} ({len(journals):,} journals)")

if __name__ == "__main__":
    main()
