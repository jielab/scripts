#!/usr/bin/env python3
"""Normalize JCR/SJR-style journal rank files for the PFI pipeline.

Output columns intentionally match the README contract:
  journal_key, journal_title, issn, eissn, jcr_year, jcr_category,
  journal_domain, impact_factor, journal_citation_indicator, jif_rank,
  jif_quartile, jif_percentile, total_cites, actual_journal_score_0_100

The input must contain real journal records. The schema CSV is documentation only.
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except Exception as e:  # pragma: no cover
    raise SystemExit("pyarrow is required: pip install pyarrow") from e


def norm_col(x: str) -> str:
    x = str(x).strip().lower()
    x = re.sub(r"[^a-z0-9]+", "_", x)
    return re.sub(r"_+", "_", x).strip("_")


def clean_title(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().lower()
    s = re.sub(r"&", " and ", s)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def clean_issn(x: object) -> str:
    if pd.isna(x):
        return ""
    s = re.sub(r"[^0-9Xx]", "", str(x))
    if len(s) == 8:
        return s[:4] + "-" + s[4:].upper()
    return str(x).strip().upper()


def first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s.astype(str).str.replace(",", "", regex=False).str.replace("%", "", regex=False), errors="coerce")


def parse_quartile(x: object) -> str:
    if pd.isna(x):
        return ""
    m = re.search(r"Q\s*([1-4])", str(x).upper())
    return f"Q{m.group(1)}" if m else str(x).strip().upper()


def percentile_from_rank(rank: pd.Series, n: pd.Series) -> pd.Series:
    rank = to_num(rank)
    n = to_num(n)
    out = 100.0 * (1.0 - (rank - 1.0) / n)
    out[(rank <= 0) | (n <= 0)] = np.nan
    return out.clip(0, 100)


def percentile_from_quartile(q: pd.Series) -> pd.Series:
    # Midpoint imputation: Q1=87.5, Q2=62.5, Q3=37.5, Q4=12.5.
    mp = {"Q1": 87.5, "Q2": 62.5, "Q3": 37.5, "Q4": 12.5}
    return q.map(mp).astype(float)


def infer_domain(category: object, title: object = "") -> str:
    txt = f"{category or ''} {title or ''}".lower()
    rules = [
        ("methods_data_software", r"bioinformatics|computational|statistics|mathematical|methods|software|database|data science|artificial intelligence|machine learning|medical informatics|imaging informatics"),
        ("public_health", r"public health|epidemiology|environmental|occupational|population|global health|health policy|health promotion|nutrition|social science|prevention"),
        ("health_services", r"health care|healthcare|health services|nursing|primary care|quality of care|implementation|medical education|health economics"),
        ("clinical_science", r"clinical|medicine|surgery|oncology|cardiology|neurology|psychiatry|pediatrics|obstetrics|radiology|diagnostic|therapeutic|trial|transplant|emergency"),
        ("basic_biomedical", r"biochemistry|molecular|cell|genetics|genomics|immunology|microbiology|neuroscience|physiology|pharmacology|toxicology|biology|biomedical"),
    ]
    for dom, pat in rules:
        if re.search(pat, txt):
            return dom
    return "other_biomedical"


def read_any_csv(path: Path) -> pd.DataFrame:
    # JCR exports can be comma-separated or tab-separated. Try permissive sniffing.
    try:
        return pd.read_csv(path, dtype=str, low_memory=False)
    except Exception:
        return pd.read_csv(path, dtype=str, sep="\t", low_memory=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", required=True)
    ap.add_argument("--schema", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--audit_csv", required=True)
    args = ap.parse_args()

    rank_path = Path(args.rank)
    df0 = read_any_csv(rank_path)
    raw_cols = list(df0.columns)
    df0.columns = [norm_col(c) for c in df0.columns]

    if set(df0.columns) == {"column", "required", "description"} or ("column" in df0.columns and "description" in df0.columns and len(df0) < 50):
        raise SystemExit(
            f"ERROR: {rank_path} looks like journal_rank_schema.csv, not a real journal-rank table. "
            "Use a real JCR/SJR/OpenAlex-derived file with journal records."
        )

    colmap = {
        "journal_key": ["journal_key", "source_id", "journal_id"],
        "journal_title": ["journal_title", "full_journal_title", "journal", "title", "source_title", "journal_name", "journals"],
        "issn": ["issn", "print_issn", "p_issn"],
        "eissn": ["eissn", "e_issn", "electronic_issn", "online_issn"],
        "jcr_year": ["jcr_year", "year", "rank_year", "sjr_year"],
        "jcr_category": ["jcr_category", "category", "categories", "subject_category", "web_of_science_categories", "asjc_category"],
        "journal_domain": ["journal_domain", "domain", "broad_domain", "field"],
        "impact_factor": ["impact_factor", "journal_impact_factor", "jif", "if", "cites_doc_2years", "cites_per_doc_2years", "sjr"],
        "journal_citation_indicator": ["journal_citation_indicator", "jci"],
        "jif_rank": ["jif_rank", "rank", "jif_rank_in_category", "sjr_rank", "rank_in_category"],
        "category_size": ["category_size", "journals_in_category", "n_journals", "total_journals", "category_count"],
        "jif_quartile": ["jif_quartile", "quartile", "jif_quartile_in_category", "sjr_quartile", "q"],
        "jif_percentile": ["jif_percentile", "percentile", "jif_percentile_in_category", "sjr_percentile", "cite_score_percentile", "citescore_percentile"],
        "total_cites": ["total_cites", "total_citations", "cites", "total_refs"],
    }

    out = pd.DataFrame(index=df0.index)
    for target, candidates in colmap.items():
        c = first_present(df0, candidates)
        out[target] = df0[c] if c else np.nan

    out["journal_title"] = out["journal_title"].fillna("").astype(str).str.strip()
    out["issn"] = out["issn"].map(clean_issn)
    out["eissn"] = out["eissn"].map(clean_issn)

    # Stable key preference: explicit key > eISSN > ISSN > normalized title.
    if out["journal_key"].isna().all() or (out["journal_key"].fillna("").astype(str).str.strip() == "").all():
        key = np.where(out["eissn"].astype(str).str.len() > 0, "eissn:" + out["eissn"].astype(str), "")
        key = np.where((key == "") & (out["issn"].astype(str).str.len() > 0), "issn:" + out["issn"].astype(str), key)
        title_key = out["journal_title"].map(clean_title).map(lambda x: "title:" + x if x else "")
        key = np.where(key == "", title_key, key)
        out["journal_key"] = key
    else:
        out["journal_key"] = out["journal_key"].astype(str).str.strip().str.lower()

    if out["jcr_year"].isna().all() or (out["jcr_year"].fillna("").astype(str).str.strip() == "").all():
        # Leave as missing rather than inventing a year. Downstream scripts can use nearest-year matching if implemented.
        out["jcr_year"] = np.nan
    else:
        out["jcr_year"] = pd.to_numeric(out["jcr_year"], errors="coerce").astype("Int64")

    out["jcr_category"] = out["jcr_category"].fillna("unknown").astype(str).str.strip()
    out.loc[out["jcr_category"].eq("") | out["jcr_category"].str.lower().eq("nan"), "jcr_category"] = "unknown"
    out["jif_quartile"] = out["jif_quartile"].map(parse_quartile)

    out["impact_factor"] = to_num(out["impact_factor"])
    out["journal_citation_indicator"] = to_num(out["journal_citation_indicator"])
    out["jif_rank"] = to_num(out["jif_rank"])
    out["category_size"] = to_num(out["category_size"])
    out["jif_percentile"] = to_num(out["jif_percentile"])
    out["total_cites"] = to_num(out["total_cites"])

    missing_percentile = out["jif_percentile"].isna()
    if missing_percentile.any() and out["jif_rank"].notna().any() and out["category_size"].notna().any():
        out.loc[missing_percentile, "jif_percentile"] = percentile_from_rank(out.loc[missing_percentile, "jif_rank"], out.loc[missing_percentile, "category_size"])

    missing_percentile = out["jif_percentile"].isna()
    if missing_percentile.any() and out["jif_quartile"].isin(["Q1", "Q2", "Q3", "Q4"]).any():
        out.loc[missing_percentile, "jif_percentile"] = percentile_from_quartile(out.loc[missing_percentile, "jif_quartile"])

    out["jif_percentile"] = out["jif_percentile"].clip(0, 100)
    out["actual_journal_score_0_100"] = out["jif_percentile"]

    if out["journal_domain"].isna().all() or (out["journal_domain"].fillna("").astype(str).str.strip() == "").all():
        out["journal_domain"] = [infer_domain(c, t) for c, t in zip(out["jcr_category"], out["journal_title"])]
    else:
        out["journal_domain"] = out["journal_domain"].fillna("").astype(str).str.strip().str.lower().str.replace(r"[^a-z0-9]+", "_", regex=True).str.strip("_")
        miss = out["journal_domain"].eq("") | out["journal_domain"].eq("nan")
        out.loc[miss, "journal_domain"] = [infer_domain(c, t) for c, t in zip(out.loc[miss, "jcr_category"], out.loc[miss, "journal_title"])]

    required = ["journal_key", "jcr_year", "jcr_category", "journal_domain", "jif_percentile", "jif_quartile"]
    audit = {
        "input_file": str(rank_path),
        "input_rows": int(len(df0)),
        "input_columns": ";".join(raw_cols),
        "output_rows_before_drop": int(len(out)),
        "missing_journal_key": int(out["journal_key"].isna().sum() + out["journal_key"].eq("").sum()),
        "missing_jif_percentile": int(out["jif_percentile"].isna().sum()),
        "missing_quartile": int(out["jif_quartile"].eq("").sum()),
        "missing_year": int(out["jcr_year"].isna().sum()),
    }

    out = out[out["journal_key"].fillna("").astype(str).str.len() > 0].copy()
    out = out[out["jif_percentile"].notna()].copy()
    if out.empty:
        raise SystemExit("ERROR: no usable journal rank records after normalization. Check rank/percentile/quartile columns.")

    keep = [
        "journal_key", "journal_title", "issn", "eissn", "jcr_year", "jcr_category", "journal_domain",
        "impact_factor", "journal_citation_indicator", "jif_rank", "jif_quartile", "jif_percentile",
        "actual_journal_score_0_100", "total_cites",
    ]
    out = out[keep].drop_duplicates()
    audit["output_rows"] = int(len(out))
    audit["n_unique_journal_keys"] = int(out["journal_key"].nunique())
    audit["n_domains"] = int(out["journal_domain"].nunique())

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), out_path)

    audit_path = Path(args.audit_csv)
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([audit]).to_csv(audit_path, index=False)
    print(f"OK: wrote {len(out):,} normalized journal-rank records to {out_path}")
    print(f"OK: wrote audit CSV to {audit_path}")


if __name__ == "__main__":
    main()
