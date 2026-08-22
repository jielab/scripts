#!/usr/bin/env python3
"""Build a supervised PubMedBERT training table for the PFI pipeline.

The script joins local PubMed/PMC text, transparent article features, and an
external journal-rank file. It creates two supervised targets:

1) journal_score_0_100: field-normalized journal prestige, usually from JCR
   percentile/quartile/impact factor within category.
2) paper_score: article-level quality/importance label. If a human/LLM
   score file is supplied, that label is used. Otherwise a deliberately weak
   transparent proxy is created from article features. The weak proxy is useful
   for bootstrapping, not for final claims.

The model is trained only on blinded text. Author names, affiliations, journal
names, funders, and publisher fields should already have been excluded by
make_blinded_text.py.
"""


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Imports, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Optional

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

def norm_journal(x: object) -> Optional[str]:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    x = str(x).lower()
    x = re.sub(r"[^a-z0-9]+", " ", x)
    x = re.sub(r"\s+", " ", x).strip()
    return x or None

def norm_id(x: object) -> Optional[str]:
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

def first_existing(df: pd.DataFrame, names: list[str]) -> Optional[str]:
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    return None

def infer_domain(category: object) -> str:
    s = str(category or "").lower()
    if not s or s == "nan":
        return "unknown"
    if any(k in s for k in [
        "public", "environment", "occupational", "epidemiology", "health care",
        "health policy", "health services", "population", "social sciences biomedical",
    ]):
        return "public_health"
    if any(k in s for k in [
        "clinical", "medicine", "surgery", "oncology", "cardiac", "cardiovascular",
        "neurology", "psychiatry", "pediatrics", "obstetrics", "radiology",
        "infectious", "gastro", "respiratory", "endocrin", "rheumatology",
        "urology", "dermatology", "hematology", "primary health",
    ]):
        return "clinical_science"
    if any(k in s for k in [
        "cell", "molecular", "genetics", "genomics", "biochemistry", "biophysics",
        "immunology", "microbiology", "neuroscience", "pharmacology", "physiology",
        "toxicology", "developmental biology", "biology", "chemistry medicinal",
    ]):
        return "basic_medical_research"
    if any(k in s for k in [
        "informatics", "statistics", "mathematical", "computer", "engineering biomedical",
        "methods", "multidisciplinary sciences",
    ]):
        return "methods_resource"
    return "other_biomedicine"

def normalize_rank_file(rank_file: str, out_path: Path) -> pd.DataFrame:
    rank = pd.read_csv(rank_file)
    rank.columns = [c.strip() for c in rank.columns]

    journal_col = first_existing(rank, ["journal", "journal_name", "full_journal_title", "title"])
    if journal_col is None:
        raise ValueError("journal rank file must contain a journal / journal_name column")
    rank["journal"] = rank[journal_col].astype(str)
    rank["journal_key"] = rank["journal"].map(norm_journal)

    cat_col = first_existing(rank, ["jcr_category", "category", "subject_category", "web_of_science_category", "field"])
    if cat_col:
        rank["jcr_category"] = rank[cat_col].fillna("unknown").astype(str)
    else:
        rank["jcr_category"] = "unknown"

    domain_col = first_existing(rank, ["domain", "broad_domain", "field_domain"])
    if domain_col:
        rank["domain"] = rank[domain_col].fillna("").astype(str).str.lower().str.replace(r"\s+", "_", regex=True)
        rank.loc[rank["domain"].isin(["", "nan", "none"]), "domain"] = rank["jcr_category"].map(infer_domain)
    else:
        rank["domain"] = rank["jcr_category"].map(infer_domain)

    score_col = first_existing(rank, [
        "journal_prestige_percentile", "jcr_percentile", "jif_percentile",
        "journal_percentile", "prestige_percentile", "percentile",
    ])
    if score_col:
        rank["journal_score_0_100"] = pd.to_numeric(rank[score_col], errors="coerce")
    else:
        jif_col = first_existing(rank, ["impact_factor", "jif", "journal_impact_factor", "if"])
        q_col = first_existing(rank, ["jcr_quartile", "quartile", "journal_tier", "tier"])
        if jif_col:
            rank["_jif"] = pd.to_numeric(rank[jif_col], errors="coerce")
            # Percentile within JCR category/domain is safer than raw IF because IF scale differs strongly by field.
            by = "jcr_category" if rank["jcr_category"].nunique(dropna=True) > 1 else "domain"
            rank["journal_score_0_100"] = rank.groupby(by)["_jif"].rank(pct=True) * 100.0
        elif q_col:
            q = rank[q_col].astype(str).str.upper().str.extract(r"Q?([1-4])")[0]
            approx = {"1": 87.5, "2": 62.5, "3": 37.5, "4": 12.5}
            rank["journal_score_0_100"] = q.map(approx).astype(float)
        else:
            raise ValueError(
                "journal rank file must contain one of: journal_prestige_percentile / JCR percentile / impact_factor / quartile"
            )

    tier_col = first_existing(rank, ["journal_tier", "jcr_quartile", "quartile", "tier"])
    if tier_col:
        rank["journal_tier"] = rank[tier_col].astype(str)
    else:
        rank["journal_tier"] = pd.cut(
            rank["journal_score_0_100"],
            bins=[-np.inf, 25, 50, 75, np.inf],
            labels=["Q4_like", "Q3_like", "Q2_like", "Q1_like"],
        ).astype(str)

    keep = ["journal", "journal_key", "jcr_category", "domain", "journal_tier", "journal_score_0_100"]
    rank = rank[keep].dropna(subset=["journal_key", "journal_score_0_100"])
    rank["journal_score_0_100"] = rank["journal_score_0_100"].clip(0, 100)

    # JCR may list the same journal in multiple categories. Use the mean
    # category percentile, mirroring JCR's "average JIF percentile" idea, while
    # preserving the category/domain labels for audit.
    rank = rank.groupby("journal_key", as_index=False).agg(
        journal=("journal", "first"),
        jcr_category=("jcr_category", lambda x: "|".join(sorted(set(map(str, x))))[:500]),
        domain=("domain", lambda x: "|".join(sorted(set(map(str, x))))[:200]),
        journal_tier=("journal_tier", lambda x: "|".join(sorted(set(map(str, x))))[:100]),
        journal_score_0_100=("journal_score_0_100", "mean"),
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rank.to_parquet(out_path, index=False)
    return rank

def read_paper_scores(path: Optional[str], keep_all: bool = False) -> Optional[pd.DataFrame]:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(path)
    if p.suffix.lower() in [".jsonl", ".ndjson"]:
        rows = [json.loads(line) for line in p.read_text(encoding="utf-8").splitlines() if line.strip()]
        lab = pd.DataFrame(rows)
    else:
        lab = pd.read_csv(path)
    if lab.empty:
        return None
    for c in ["pmid", "pmcid", "doi"]:
        if c not in lab.columns:
            lab[c] = None
    lab["article_key"] = article_key_from(lab)

    # Flexible label columns. Values may be 0-5, 0-10, or 0-100; rescaled below.
    components = {
        "new_results": ["new_results", "new_discovery", "discovery_strength", "novelty", "discovery"],
        "new_methods": ["new_methods", "methods_innovation", "innovation", "method_innovation"],
        "new_guidelines": ["new_guidelines", "guideline_change", "practice_guidelines"],
        "big_impact": ["big_impact", "societal_impact", "clinical_importance", "public_health_importance", "impact"],
        "big_use": ["big_use", "large_resource", "data_source", "sample_resource", "resource_scale"],
        "claim_gap": ["claim_evidence_gap", "overclaim", "claim_gap"],
    }
    vals = pd.DataFrame(index=lab.index)
    for name, aliases in components.items():
        col = first_existing(lab, aliases)
        vals[name] = pd.to_numeric(lab[col], errors="coerce") if col else np.nan

    # Infer scale for each component and map to 0-100. claim_gap is bad and subtracted.
    for c in vals.columns:
        m = vals[c].max(skipna=True)
        if pd.isna(m):
            continue
        if m <= 5:
            vals[c] = vals[c] / 5 * 100
        elif m <= 10:
            vals[c] = vals[c] / 10 * 100
        else:
            vals[c] = vals[c].clip(0, 100)

    def fill_component(name: str, default: float = 50.0) -> pd.Series:
        med = vals[name].median(skipna=True)
        if pd.isna(med):
            med = default
        return vals[name].fillna(med)

    direct_score_col = first_existing(lab, ["paper_score"])
    if direct_score_col:
        lab["paper_score_external"] = pd.to_numeric(lab[direct_score_col], errors="coerce").clip(0, 100)
    else:
        score = (
            0.25 * fill_component("new_results")
            + 0.25 * fill_component("new_methods")
            + 0.10 * fill_component("new_guidelines")
            + 0.20 * fill_component("big_impact")
            + 0.20 * fill_component("big_use")
            - 0.15 * vals["claim_gap"].fillna(0)
        )
        lab["paper_score_external"] = score.clip(0, 100)
    if keep_all:
        return lab.dropna(subset=["article_key"]).drop_duplicates("article_key")
    return lab[["article_key", "paper_score_external"]].dropna().drop_duplicates("article_key")

def weak_paper_score(df: pd.DataFrame) -> pd.Series:
    q = pd.to_numeric(df.get("quality_proxy_0_20"), errors="coerce").fillna(0) / 20 * 65
    wc = np.log1p(pd.to_numeric(df.get("word_count"), errors="coerce").fillna(0))
    wc = 10 * (wc - wc.min()) / (wc.max() - wc.min() + 1e-9)
    ss = np.log1p(pd.to_numeric(df.get("sample_size_max"), errors="coerce").fillna(0))
    ss = 15 * (ss - ss.min()) / (ss.max() - ss.min() + 1e-9)
    flags = pd.Series(0.0, index=df.index)
    for c, w in [
        ("has_randomized", 4), ("has_gwas", 3), ("has_mr", 3),
        ("has_external_validation", 4), ("has_sensitivity_analysis", 2),
        ("has_data_availability", 2), ("has_code_availability", 2),
        ("has_trial_registration", 2),
        ("has_method_innovation", 5), ("has_new_discovery_language", 5),
        ("has_large_resource_language", 5), ("has_replication", 3),
        ("has_open_science_language", 2),
    ]:
        if c in df.columns:
            flags = flags + w * pd.to_numeric(df[c], errors="coerce").fillna(0).clip(0, 1)
    return (q + wc + ss + flags).clip(0, 100)

def assign_split(df: pd.DataFrame, seed: int) -> pd.Series:
    rng = np.random.default_rng(seed)
    split = pd.Series("train", index=df.index)
    # Stratify coarsely by domain and journal-score quintile.
    try:
        bin_ = pd.qcut(df["journal_score_0_100"], 5, labels=False, duplicates="drop")
    except Exception:
        bin_ = pd.Series(0, index=df.index)
    strata = df["domain"].astype(str) + "_" + bin_.astype(str)
    for _, idx in df.groupby(strata).groups.items():
        idx = np.array(list(idx))
        r = rng.random(len(idx))
        split.iloc[idx[r < 0.10]] = "test"
        split.iloc[idx[(r >= 0.10) & (r < 0.20)]] = "valid"
    return split

def write_training_diagnostics(df: pd.DataFrame, summary_out: str | None, fit_summary_out: str | None, audit_csv_out: str | None) -> None:
    if summary_out:
        rows = []
        for cols in [["split"], ["domain"], ["journal_tier"], ["paper_score_source"], ["split", "domain"]]:
            cols = [c for c in cols if c in df.columns]
            if not cols:
                continue
            g = df.groupby(cols, dropna=False).agg(
                n=("article_key", "size"),
                full_text_n=("has_pmc_full_text", "sum"),
                mean_journal_score=("journal_score_0_100", "mean"),
                mean_paper_score=("paper_score", "mean"),
                mean_word_count=("word_count", "mean"),
            ).reset_index()
            g.insert(0, "summary_level", "+".join(cols))
            rows.append(g)
        if rows:
            out = pd.concat(rows, ignore_index=True, sort=False)
            Path(summary_out).parent.mkdir(parents=True, exist_ok=True)
            out.to_csv(summary_out, index=False)

    if fit_summary_out:
        rows = []
        for name, d in [("all", df)] + [(f"domain={k}", v) for k, v in df.groupby("domain", dropna=False)]:
            if len(d) < 3:
                continue
            rows.append({
                "group": name,
                "n": len(d),
                "journal_paper_pearson": d["journal_score_0_100"].corr(d["paper_score"], method="pearson"),
                "journal_paper_spearman": d["journal_score_0_100"].corr(d["paper_score"], method="spearman"),
                "mean_journal_score": d["journal_score_0_100"].mean(),
                "mean_paper_score": d["paper_score"].mean(),
                "full_text_rate": d["has_pmc_full_text"].mean(),
            })
        Path(fit_summary_out).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(fit_summary_out, index=False)

    if audit_csv_out:
        audit = df.drop(columns=["blinded_text"], errors="ignore").copy()
        audit["blinded_text_preview"] = df["blinded_text"].fillna("").astype(str).str.slice(0, 1500)
        Path(audit_csv_out).parent.mkdir(parents=True, exist_ok=True)
        audit.to_csv(audit_csv_out, index=False)

def maybe_write_fast_training_from_review_csv(
    labels_full: Optional[pd.DataFrame],
    rank: pd.DataFrame,
    out: Path,
    args: argparse.Namespace,
) -> bool:
    if labels_full is None or labels_full.empty:
        return False
    text_col = first_existing(labels_full, ["blinded_text", "blinded_text_preview", "review_text_preview", "text"])
    journal_col = first_existing(labels_full, ["journal", "journal_name", "full_journal_title"])
    if text_col is None or journal_col is None:
        return False

    df = labels_full.copy()
    df["journal"] = df[journal_col].astype(str)
    df["journal_key"] = df["journal"].map(norm_journal)
    df["blinded_text"] = df[text_col].fillna("").astype(str)
    if "has_pmc_full_text" in df.columns:
        df = df[pd.to_numeric(df["has_pmc_full_text"], errors="coerce").fillna(0).astype(int).eq(1)]
    df = df[df["blinded_text"].str.len() >= int(args.min_text_chars)]
    df = df.merge(
        rank[["journal_key", "jcr_category", "domain", "journal_tier", "journal_score_0_100"]],
        on="journal_key",
        how="inner",
    )
    if df.empty:
        return False
    if args.max_rows and args.max_rows > 0 and len(df) > args.max_rows:
        df = df.sample(int(args.max_rows), random_state=args.seed)

    for c in ["pmid", "pmcid", "doi", "year", "title", "has_pmc_full_text", "word_count", "sample_size_max"]:
        if c not in df.columns:
            df[c] = None
    df["paper_score"] = pd.to_numeric(df["paper_score_external"], errors="coerce").clip(0, 100)
    df["paper_score_weak"] = np.nan
    df["paper_score_source"] = "external"
    df["journal_score_0_100"] = pd.to_numeric(df["journal_score_0_100"], errors="coerce").clip(0, 100)
    df = df.dropna(subset=["journal_score_0_100", "paper_score", "blinded_text"])
    df["split"] = assign_split(df, args.seed)

    keep = [
        "article_key", "pmid", "pmcid", "doi", "year", "title", "journal", "jcr_category", "domain",
        "journal_tier", "has_pmc_full_text", "word_count", "sample_size_max", "paper_score_source",
        "journal_score_0_100", "paper_score", "paper_score_weak", "split", "blinded_text",
    ]
    keep = [c for c in keep if c in df.columns]
    out.parent.mkdir(parents=True, exist_ok=True)
    df[keep].to_parquet(out, index=False)

    summary = df.groupby(["split", "domain"], dropna=False).size().reset_index(name="n")
    summary_path = out.with_suffix(".summary.csv")
    summary.to_csv(summary_path, index=False)
    write_training_diagnostics(df, args.summary_out, args.fit_summary_out, args.audit_csv_out)
    print(f"Wrote training dataset from review CSV fast path: {out} ({len(df):,} rows)")
    print(f"Wrote summary: {summary_path}")
    if args.audit_csv_out:
        print(f"Wrote training audit CSV: {args.audit_csv_out}")
    if args.summary_out:
        print(f"Wrote training dataset summary: {args.summary_out}")
    if args.fit_summary_out:
        print(f"Wrote training label-fit summary: {args.fit_summary_out}")
    print("Label source counts:")
    print(df["paper_score_source"].value_counts(dropna=False).to_string())
    return True

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--blinded", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--journal_rank", required=True)
    ap.add_argument("--paper_scores", default=None, help="Optional CSV/JSONL with human/LLM paper-quality labels")
    ap.add_argument("--out", required=True)
    ap.add_argument("--rank_out", required=True)
    ap.add_argument("--audit_csv_out", default=None)
    ap.add_argument("--summary_out", default=None)
    ap.add_argument("--fit_summary_out", default=None)
    ap.add_argument("--max_rows", type=int, default=1000)
    ap.add_argument("--min_text_chars", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rank_path = Path(args.rank_out)
    rank = normalize_rank_file(args.journal_rank, rank_path)
    print(f"Normalized journal ranks: {len(rank):,} journals -> {rank_path}")
    labels_full = read_paper_scores(args.paper_scores, keep_all=True)
    if maybe_write_fast_training_from_review_csv(labels_full, rank, out, args):
        return
    labels = (
        labels_full[["article_key", "paper_score_external"]].dropna().drop_duplicates("article_key")
        if labels_full is not None and not labels_full.empty
        else None
    )

    limit = f"LIMIT {int(args.max_rows)}" if args.max_rows and args.max_rows > 0 else ""
    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    label_join = ""
    if labels is not None and not labels.empty:
        con.register("label_keys", labels[["article_key"]].dropna().drop_duplicates())
        label_join = "JOIN label_keys USING(article_key)"
    q = f"""
    WITH rank AS (
      SELECT * FROM read_parquet('{rank_path.as_posix()}')
    ), joined AS (
      SELECT
        coalesce(CAST(m.pmid AS VARCHAR), CAST(m.pmcid AS VARCHAR), CAST(m.doi AS VARCHAR)) AS article_key,
        coalesce(m.pmid, b.pmid, f.pmid) AS pmid,
        coalesce(m.pmcid, b.pmcid, f.pmcid) AS pmcid,
        coalesce(m.doi, b.doi, f.doi) AS doi,
        coalesce(m.year, b.year, f.year) AS year,
        m.title, m.journal, m.publication_types, m.mesh_terms,
        m.author_count,
        CASE WHEN coalesce(m.full_text, m.intro_text, m.methods_text, m.results_text, m.discussion_text,
                           m.data_availability_text, m.ethics_text, m.trial_registration_text) IS NOT NULL
             THEN 1 ELSE 0 END AS has_pmc_full_text,
        b.blinded_text,
        f.* EXCLUDE(pmid, pmcid, doi, year),
        trim(regexp_replace(regexp_replace(lower(coalesce(m.journal, '')), '[^a-z0-9]+', ' ', 'g'), '\\s+', ' ', 'g')) AS journal_key
      FROM read_parquet('{parquet_source(args.master)}', union_by_name=true, hive_partitioning=1) m
      JOIN read_parquet('{parquet_source(args.blinded)}', union_by_name=true, hive_partitioning=1) b
        ON coalesce(m.pmid, m.pmcid, m.doi) = coalesce(b.pmid, b.pmcid, b.doi)
      LEFT JOIN read_parquet('{parquet_source(args.features)}', union_by_name=true, hive_partitioning=1) f
        ON coalesce(m.pmid, m.pmcid, m.doi) = coalesce(f.pmid, f.pmcid, f.doi)
      WHERE b.blinded_text IS NOT NULL AND length(b.blinded_text) >= {int(args.min_text_chars)}
        AND coalesce(m.full_text, m.intro_text, m.methods_text, m.results_text, m.discussion_text,
                     m.data_availability_text, m.ethics_text, m.trial_registration_text) IS NOT NULL
    )
    SELECT joined.*, rank.jcr_category, rank.domain, rank.journal_tier, rank.journal_score_0_100
    FROM joined
    JOIN rank USING(journal_key)
    {label_join}
    ORDER BY random()
    {limit}
    """
    df = con.execute(q).fetchdf()
    con.close()
    print(f"Rows after journal-rank join and text filter: {len(df):,}")
    if df.empty:
        raise RuntimeError("No rows matched the journal-rank file. Check journal names / ISSN mapping.")

    for c in ["pmid", "pmcid", "doi"]:
        if c not in df.columns:
            df[c] = None
    df["article_key"] = df["article_key"].where(df["article_key"].notna(), article_key_from(df))
    df["paper_score_weak"] = weak_paper_score(df)

    if labels is not None and not labels.empty:
        df = df.merge(labels, on="article_key", how="left")
        df["paper_score"] = df["paper_score_external"].fillna(df["paper_score_weak"])
        df["paper_score_source"] = np.where(df["paper_score_external"].notna(), "external", "weak_proxy")
    else:
        df["paper_score"] = df["paper_score_weak"]
        df["paper_score_source"] = "weak_proxy"

    df["journal_score_0_100"] = pd.to_numeric(df["journal_score_0_100"], errors="coerce").clip(0, 100)
    df["paper_score"] = pd.to_numeric(df["paper_score"], errors="coerce").clip(0, 100)
    df = df.dropna(subset=["journal_score_0_100", "paper_score", "blinded_text"])
    df["split"] = assign_split(df, args.seed)

    keep = [
        "article_key", "pmid", "pmcid", "doi", "year", "title", "journal", "jcr_category", "domain",
        "journal_tier", "publication_types", "mesh_terms", "has_pmc_full_text", "author_count",
        "word_count", "sample_size_max", "quality_proxy_0_20", "paper_score_source",
        "journal_score_0_100", "paper_score", "paper_score_weak", "split", "blinded_text",
    ]
    keep = [c for c in keep if c in df.columns]
    df[keep].to_parquet(out, index=False)

    summary = df.groupby(["split", "domain"], dropna=False).size().reset_index(name="n")
    summary_path = out.with_suffix(".summary.csv")
    summary.to_csv(summary_path, index=False)
    write_training_diagnostics(df, args.summary_out, args.fit_summary_out, args.audit_csv_out)
    print(f"Wrote training dataset: {out} ({len(df):,} rows)")
    print(f"Wrote summary: {summary_path}")
    if args.audit_csv_out:
        print(f"Wrote training audit CSV: {args.audit_csv_out}")
    if args.summary_out:
        print(f"Wrote training dataset summary: {args.summary_out}")
    if args.fit_summary_out:
        print(f"Wrote training label-fit summary: {args.fit_summary_out}")
    print("Label source counts:")
    print(df["paper_score_source"].value_counts(dropna=False).to_string())

if __name__ == "__main__":
    main()
