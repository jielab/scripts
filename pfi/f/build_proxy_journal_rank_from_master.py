#!/usr/bin/env python3
"""Build a temporary proxy journal-rank CSV directly from the PFI master dataset.

This is only a run-through fallback when no real JCR/SJR journal rank file is
available. It estimates a within-domain/year journal placement proxy from the
number of full-text articles observed in the current PubMed/PMC-OA corpus. The
output intentionally mimics the PFI journal-rank schema so 08b can normalize it.

Do not interpret this proxy as a validated journal-prestige score.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import pyarrow.dataset as ds
except Exception as e:  # pragma: no cover
    raise SystemExit("pyarrow is required: pip install pyarrow") from e


def norm_col(x: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(x).lower())).strip("_")


def clean_title_value(x: object) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip().lower()
    s = s.replace("&", " and ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def clean_issn_series(s: pd.Series) -> pd.Series:
    def one(x: object) -> str:
        if pd.isna(x):
            return ""
        y = re.sub(r"[^0-9Xx]", "", str(x))
        if len(y) == 8:
            return y[:4] + "-" + y[4:].upper()
        return str(x).strip().upper()
    return s.map(one)


def infer_domain_text(category: object, title: object = "") -> str:
    txt = f"{category or ''} {title or ''}".lower()
    rules = [
        ("methods_data_software", r"bioinformatics|computational|statistics|mathematical|methods|software|database|data science|artificial intelligence|machine learning|medical informatics|imaging informatics"),
        ("public_health", r"public health|epidemiology|environmental|occupational|population|global health|health policy|health promotion|nutrition|social science|prevention"),
        ("health_services", r"health care|healthcare|health services|nursing|primary care|quality of care|implementation|medical education|health economics"),
        ("clinical_science", r"clinical|medicine|surgery|oncology|cardiology|neurology|psychiatry|pediatrics|obstetrics|radiology|diagnostic|therapeutic|trial|transplant|emergency|bmj|lancet|jama|new england"),
        ("basic_biomedical", r"biochemistry|molecular|cell|genetics|genomics|immunology|microbiology|neuroscience|physiology|pharmacology|toxicology|biology|biomedical|nature|science"),
    ]
    for dom, pat in rules:
        if re.search(pat, txt):
            return dom
    return "other_biomedical"


def choose_col(names: set[str], candidates: list[str]) -> str | None:
    low = {norm_col(n): n for n in names}
    for c in candidates:
        nc = norm_col(c)
        if nc in low:
            return low[nc]
    return None


def read_schema_names(path: Path) -> list[str]:
    dataset = ds.dataset(str(path), format="parquet")
    return list(dataset.schema.names)


def make_key(df: pd.DataFrame, c_jkey: str | None, c_eissn: str | None, c_issn: str | None, c_title: str | None) -> pd.Series:
    n = len(df)
    key = pd.Series([""] * n, index=df.index, dtype=object)
    if c_jkey and c_jkey in df.columns:
        key = df[c_jkey].fillna("").astype(str).str.strip().str.lower()
    if c_eissn and c_eissn in df.columns:
        e = clean_issn_series(df[c_eissn])
        key = np.where(pd.Series(key).astype(str).str.len() > 0, key, np.where(e.astype(str).str.len() > 0, "eissn:" + e.astype(str), ""))
    if c_issn and c_issn in df.columns:
        i = clean_issn_series(df[c_issn])
        key = np.where(pd.Series(key).astype(str).str.len() > 0, key, np.where(i.astype(str).str.len() > 0, "issn:" + i.astype(str), ""))
    if c_title and c_title in df.columns:
        t = df[c_title].map(clean_title_value)
        key = np.where(pd.Series(key).astype(str).str.len() > 0, key, np.where(t.astype(str).str.len() > 0, "title:" + t.astype(str), ""))
    return pd.Series(key, index=df.index).astype(str)


def percentile_from_rank(rank: pd.Series, n: pd.Series) -> pd.Series:
    rank = pd.to_numeric(rank, errors="coerce")
    n = pd.to_numeric(n, errors="coerce")
    denom = (n - 1).replace(0, np.nan)
    p = 100.0 * (1.0 - (rank - 1.0) / denom)
    p = p.where(n > 1, 50.0)
    return p.clip(0, 100).fillna(50.0)


def quartile_from_percentile(p: pd.Series) -> pd.Series:
    return pd.cut(p, bins=[-0.001, 25, 50, 75, 100.001], labels=["Q4", "Q3", "Q2", "Q1"]).astype(str)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True, help="PFI articles_master parquet dataset")
    ap.add_argument("--out", required=True, help="Output proxy rank CSV")
    ap.add_argument("--audit_csv", required=True)
    ap.add_argument("--min_year", type=int, default=2000)
    ap.add_argument("--full_text_only", action="store_true")
    ap.add_argument("--batch_size", type=int, default=200000)
    args = ap.parse_args()

    master = Path(args.master)
    if not master.exists():
        raise SystemExit(f"ERROR: master dataset is missing: {master}")

    names = read_schema_names(master)
    name_set = set(names)
    c_year = choose_col(name_set, ["year", "publication_year", "pub_year", "article_year", "jcr_year"])
    c_jkey = choose_col(name_set, ["journal_key", "source_key", "journal_id"])
    c_title = choose_col(name_set, ["journal_title", "journal", "journal_name", "full_journal_title", "source_title", "medline_ta"])
    c_issn = choose_col(name_set, ["issn", "print_issn", "journal_issn"])
    c_eissn = choose_col(name_set, ["eissn", "e_issn", "electronic_issn", "journal_eissn"])
    c_domain = choose_col(name_set, ["domain", "journal_domain", "article_domain", "field", "broad_domain"])
    c_cat = choose_col(name_set, ["jcr_category", "category", "journal_category", "subject_category", "mesh_heading_major", "mesh_terms"])
    c_full = choose_col(name_set, ["has_pmc_full_text", "full_text", "is_full_text", "pmc_full_text"])

    if not any([c_jkey, c_title, c_issn, c_eissn]):
        raise SystemExit("ERROR: could not find journal_key, journal title, ISSN, or eISSN in articles_master.")

    cols = [c for c in [c_year, c_jkey, c_title, c_issn, c_eissn, c_domain, c_cat, c_full] if c]
    dataset = ds.dataset(str(master), format="parquet")

    counts: Counter[tuple] = Counter()
    total_rows = 0
    kept_rows = 0
    for batch in dataset.to_batches(columns=cols, batch_size=args.batch_size):
        df = batch.to_pandas()
        total_rows += len(df)
        if df.empty:
            continue
        if c_full and args.full_text_only:
            ft = df[c_full].fillna(0).astype(str).str.lower().isin(["1", "true", "t", "yes", "y"])
            df = df[ft].copy()
            if df.empty:
                continue
        if c_year and c_year in df.columns:
            year = pd.to_numeric(df[c_year], errors="coerce").fillna(args.min_year).astype(int)
        else:
            year = pd.Series([args.min_year] * len(df), index=df.index)
        df = df[(year >= args.min_year) & (year <= 2100)].copy()
        year = year.loc[df.index]
        if df.empty:
            continue
        key = make_key(df, c_jkey, c_eissn, c_issn, c_title)
        title = df[c_title].fillna("").astype(str).str.strip() if c_title else pd.Series([""] * len(df), index=df.index)
        issn = clean_issn_series(df[c_issn]) if c_issn else pd.Series([""] * len(df), index=df.index)
        eissn = clean_issn_series(df[c_eissn]) if c_eissn else pd.Series([""] * len(df), index=df.index)
        cat_raw = df[c_cat].fillna("").astype(str).str.strip() if c_cat else pd.Series([""] * len(df), index=df.index)
        domain_raw = df[c_domain].fillna("").astype(str).str.strip().str.lower().str.replace(r"[^a-z0-9]+", "_", regex=True).str.strip("_") if c_domain else pd.Series([""] * len(df), index=df.index)
        domain = []
        cat = []
        for cr, dr, tt in zip(cat_raw, domain_raw, title):
            d = dr if dr and dr != "nan" else infer_domain_text(cr, tt)
            c = cr if cr and cr.lower() != "nan" else d
            # Very long MeSH/category strings create too many groups; use broad domain as category in proxy mode.
            if len(c) > 80 or ";" in c or "|" in c:
                c = d
            domain.append(d)
            cat.append(c)
        tmp = pd.DataFrame({
            "journal_key": key,
            "journal_title": title,
            "issn": issn,
            "eissn": eissn,
            "jcr_year": year,
            "jcr_category": cat,
            "journal_domain": domain,
        })
        tmp = tmp[tmp["journal_key"].astype(str).str.len() > 0]
        kept_rows += len(tmp)
        if tmp.empty:
            continue
        g = tmp.groupby(["journal_key", "journal_title", "issn", "eissn", "jcr_year", "jcr_category", "journal_domain"], dropna=False).size()
        for k, v in g.items():
            counts[tuple(k)] += int(v)

    if not counts:
        raise SystemExit("ERROR: no usable journal records found in articles_master for proxy rank.")

    rows = []
    for k, v in counts.items():
        rows.append((*k, v))
    out = pd.DataFrame(rows, columns=["journal_key", "journal_title", "issn", "eissn", "jcr_year", "jcr_category", "journal_domain", "n_articles_in_master"])

    # Add all-year count as a stabilizer; rank within year/category by all-year count first, then current-year count.
    all_year = out.groupby("journal_key", as_index=False)["n_articles_in_master"].sum().rename(columns={"n_articles_in_master": "n_articles_all_years"})
    out = out.merge(all_year, on="journal_key", how="left")
    out = out.sort_values(["jcr_year", "jcr_category", "n_articles_all_years", "n_articles_in_master", "journal_key"], ascending=[True, True, False, False, True])
    out["category_size"] = out.groupby(["jcr_year", "jcr_category"])["journal_key"].transform("nunique")
    out["jif_rank"] = out.groupby(["jcr_year", "jcr_category"])["n_articles_all_years"].rank(method="dense", ascending=False).astype(int)
    out["jif_percentile"] = percentile_from_rank(out["jif_rank"], out["category_size"])
    # Shrink extreme proxy percentiles because this is not real JCR/SJR.
    out["jif_percentile"] = 10 + 0.8 * out["jif_percentile"]
    out["jif_quartile"] = quartile_from_percentile(out["jif_percentile"])
    out["impact_factor"] = np.nan
    out["journal_citation_indicator"] = np.nan
    out["total_cites"] = np.nan
    out["proxy_source"] = "internal_article_volume_from_pmc_fulltext_master"

    keep = [
        "journal_key", "journal_title", "issn", "eissn", "jcr_year", "jcr_category", "journal_domain",
        "impact_factor", "journal_citation_indicator", "jif_rank", "category_size", "jif_quartile", "jif_percentile",
        "total_cites", "n_articles_in_master", "n_articles_all_years", "proxy_source",
    ]
    out = out[keep]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    audit = {
        "warning": "Proxy journal rank built from internal article volume only; use for pipeline testing, not validated final inference.",
        "master": str(master),
        "total_rows_scanned": int(total_rows),
        "rows_used": int(kept_rows),
        "output_rows": int(len(out)),
        "unique_journals": int(out["journal_key"].nunique()),
        "years": f"{int(out['jcr_year'].min())}-{int(out['jcr_year'].max())}",
        "columns_used": json.dumps({
            "year": c_year, "journal_key": c_jkey, "journal_title": c_title, "issn": c_issn, "eissn": c_eissn,
            "domain": c_domain, "category": c_cat, "full_text_flag": c_full,
        }, ensure_ascii=False),
    }
    Path(args.audit_csv).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([audit]).to_csv(args.audit_csv, index=False)
    print(f"OK: wrote proxy journal-rank CSV to {args.out}")
    print(f"OK: wrote proxy audit CSV to {args.audit_csv}")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
