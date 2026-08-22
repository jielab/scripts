#!/usr/bin/env python3
"""Compute PFI fairness index from predictions and journal-rank data."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


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


def read_any(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if p.is_dir() or p.suffix.lower() == ".parquet":
        return ds.dataset(str(p), format="parquet", partitioning="hive").to_table().to_pandas()
    return pd.read_csv(p, dtype=str, low_memory=False)


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
    use=[]
    for nm in want_names:
        c=first_existing(cols,[nm])
        if c and c not in use: use.append(c)
    id_sets={}
    for c in ["article_id","pmid","pmcid","doi"]:
        if c in ids.columns:
            vals=set(ids[c].dropna().astype(str))-{"","nan","None"}
            if vals: id_sets[c]=vals
    if not id_sets:
        return pd.DataFrame()
    out=[]
    for batch in d.scanner(columns=use, batch_size=65536, use_threads=True).to_batches():
        df=batch.to_pandas(); df.columns=[norm_col(c) for c in df.columns]
        mask=pd.Series(False,index=df.index)
        for c, vals in id_sets.items():
            if c in df.columns:
                mask = mask | df[c].astype(str).isin(vals)
        if mask.any(): out.append(df.loc[mask].copy())
    return pd.concat(out, ignore_index=True).drop_duplicates() if out else pd.DataFrame()


def merge_by_any(left, right, suffix):
    if right.empty: return left
    left=left.copy(); right=right.copy()
    left.columns=[norm_col(c) for c in left.columns]; right.columns=[norm_col(c) for c in right.columns]
    for key in ["article_id","pmid","pmcid","doi"]:
        if key in left.columns and key in right.columns:
            return left.merge(right.drop_duplicates(key), on=key, how="left", suffixes=("",suffix))
    return left


def coalesce_columns(df, base):
    alts=[c for c in df.columns if c==base or c.startswith(base+"_")]
    if not alts: return df
    val=df[alts[0]].copy()
    for c in alts[1:]:
        val=val.where(val.notna() & val.astype(str).str.strip().ne(""), df[c])
    df[base]=val
    return df


def infer_domain(mesh, journal=""):
    x=f"{mesh} {journal}".lower()
    if re.search(r"public health|epidemiol|population|environment|policy|global health", x): return "public_health"
    if re.search(r"clinical|trial|patient|surgery|oncology|cardiology|neurology|pediatr|psychiatry", x): return "clinical_science"
    if re.search(r"cell|molecular|gene|genom|protein|immunol|microbiol|biochem", x): return "basic_medical_research"
    if re.search(r"method|statistics|bioinformatics|database|software|algorithm|machine learning", x): return "methods_data_software"
    return "other_biomedical"


def robust_z(x):
    x=pd.to_numeric(x, errors="coerce")
    med=x.median()
    mad=(x-med).abs().median()
    if not np.isfinite(mad) or mad <= 1e-9:
        sd=x.std()
        return (x-x.mean())/(sd if np.isfinite(sd) and sd>1e-9 else 1.0)
    return 0.6745*(x-med)/mad


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--journal_rank", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--out_parquet", required=True)
    ap.add_argument("--summary_dir", required=True)
    ap.add_argument("--institution_rank", default="")
    ap.add_argument("--full_text_only", action="store_true")
    args=ap.parse_args()

    pred=read_any(args.predictions); pred.columns=[norm_col(c) for c in pred.columns]
    if pred.empty: raise SystemExit(f"ERROR: empty predictions: {args.predictions}")
    ids=pred[[c for c in ["article_id","pmid","pmcid","doi"] if c in pred.columns]].copy()
    meta=lookup_master(args.master, ids)
    dat=merge_by_any(pred, meta, "_master")
    for base in ["year","journal_key","journal","journal_title","issn","eissn","mesh_terms","domain","publication_type_proxy","has_pmc_full_text","authors","author_names","affiliations","institutions"]:
        dat=coalesce_columns(dat, base)
    if args.full_text_only and "has_pmc_full_text" in dat.columns:
        ft=dat["has_pmc_full_text"].fillna(0).astype(str).str.lower().isin(["1","true","t","yes","y"])
        dat=dat[ft].copy()
    if "year" not in dat.columns: dat["year"]=2000
    dat["year"]=pd.to_numeric(dat["year"], errors="coerce").fillna(2000).astype(int)
    if "journal_key" not in dat.columns: dat["journal_key"]=""
    miss=dat["journal_key"].fillna("").astype(str).str.len()==0
    if miss.any(): dat.loc[miss,"journal_key"]=dat.loc[miss].apply(make_journal_key_row, axis=1)

    rank=read_any(args.journal_rank); rank.columns=[norm_col(c) for c in rank.columns]
    if "jcr_year" not in rank.columns: rank["jcr_year"]=dat["year"].median() if len(dat) else 2000
    if "actual_journal_score_0_100" not in rank.columns:
        rank["actual_journal_score_0_100"]=pd.to_numeric(rank.get("jif_percentile",50), errors="coerce").fillna(50)
    rank["jcr_year"]=pd.to_numeric(rank["jcr_year"], errors="coerce").fillna(2000).astype(int)
    rank["actual_journal_score_0_100"]=pd.to_numeric(rank["actual_journal_score_0_100"], errors="coerce").fillna(50).clip(0,100)
    keep=[c for c in ["journal_key","jcr_year","jcr_category","journal_domain","actual_journal_score_0_100","jif_quartile","jif_percentile"] if c in rank.columns]
    r1=rank[keep].drop_duplicates(["journal_key","jcr_year"])
    dat=dat.merge(r1, left_on=["journal_key","year"], right_on=["journal_key","jcr_year"], how="left")
    if dat["actual_journal_score_0_100"].isna().any():
        r2=rank.groupby("journal_key", as_index=False).agg(actual_journal_score_fallback=("actual_journal_score_0_100","mean"))
        dat=dat.merge(r2,on="journal_key",how="left")
        dat.loc[dat["actual_journal_score_0_100"].isna(),"actual_journal_score_0_100"]=dat.loc[dat["actual_journal_score_0_100"].isna(),"actual_journal_score_fallback"]
    dat["actual_journal_score_0_100"]=pd.to_numeric(dat["actual_journal_score_0_100"], errors="coerce").fillna(50).clip(0,100)
    if "expected_journal_score_0_100" not in dat.columns:
        c = "pred_expected_journal_score_0_100" if "pred_expected_journal_score_0_100" in dat.columns else "predicted_paper_score" if "predicted_paper_score" in dat.columns else None
        dat["expected_journal_score_0_100"] = pd.to_numeric(dat[c], errors="coerce") if c else 50
    dat["expected_journal_score_0_100"]=pd.to_numeric(dat["expected_journal_score_0_100"], errors="coerce").fillna(50).clip(0,100)
    if "predicted_paper_score" not in dat.columns:
        dat["predicted_paper_score"] = pd.to_numeric(dat.get("pred_paper_score", dat["expected_journal_score_0_100"]), errors="coerce").fillna(50).clip(0,100)
    if "domain" not in dat.columns or dat["domain"].fillna("").astype(str).eq("").all():
        if "journal_domain" in dat.columns:
            dat["domain"]=dat["journal_domain"].fillna("unknown").astype(str)
        else:
            dat["domain"]=[infer_domain(m, j) for m,j in zip(dat.get("mesh_terms",pd.Series([""]*len(dat))), dat.get("journal",dat.get("journal_title",pd.Series([""]*len(dat)))))]
    if "publication_type_proxy" not in dat.columns: dat["publication_type_proxy"]="research_article"
    dat["fairness_index"]=dat["actual_journal_score_0_100"]-dat["expected_journal_score_0_100"]

    dat["fairness_z"]=np.nan
    group_cols=["year","domain","publication_type_proxy"]
    for _, idx in dat.groupby(group_cols, dropna=False).groups.items():
        idx=list(idx)
        if len(idx)>=30:
            dat.loc[idx,"fairness_z"]=robust_z(dat.loc[idx,"fairness_index"])
    if dat["fairness_z"].isna().any():
        for _, idx in dat[dat["fairness_z"].isna()].groupby("domain", dropna=False).groups.items():
            idx=list(idx)
            if len(idx)>=10:
                dat.loc[idx,"fairness_z"]=robust_z(dat.loc[idx,"fairness_index"])
    dat.loc[dat["fairness_z"].isna(),"fairness_z"]=robust_z(dat.loc[dat["fairness_z"].isna(),"fairness_index"])
    z=dat["fairness_z"]
    dat["fairness_category"]=np.select(
        [z>=2, z>=1, z<=-2, z<=-1],
        ["strongly_unfair_favored","modestly_unfair_favored","strongly_unfair_unfavored","modestly_unfair_unfavored"],
        default="fair"
    )

    aff = dat["affiliations"].fillna("").astype(str) if "affiliations" in dat.columns else pd.Series([""] * len(dat), index=dat.index)
    inst = dat["institutions"].fillna("").astype(str) if "institutions" in dat.columns else pd.Series([""] * len(dat), index=dat.index)
    inst_text = (aff + " " + inst).str.lower()
    famous_pat=r"harvard|stanford|mit|massachusetts institute|oxford|cambridge|johns hopkins|ucsf|columbia university|university of pennsylvania|yale|princeton|mayo clinic|nih|national institutes of health"
    dat["famous_institution_proxy"]=inst_text.str.contains(famous_pat, regex=True, na=False).astype(int)

    out_cols=[c for c in [
        "article_id","pmid","pmcid","doi","year","domain","publication_type_proxy","journal_key","journal","journal_title","issn","eissn",
        "actual_journal_score_0_100","expected_journal_score_0_100","predicted_paper_score","fairness_index","fairness_z","fairness_category",
        "famous_institution_proxy","authors","author_names","affiliations","institutions"
    ] if c in dat.columns]
    out=dat[out_cols].copy()
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv,index=False)
    pq.write_table(pa.Table.from_pandas(out, preserve_index=False), args.out_parquet)
    sdir=Path(args.summary_dir); sdir.mkdir(parents=True, exist_ok=True)
    out.groupby(["domain","fairness_category"], dropna=False).size().reset_index(name="n").to_csv(sdir/"fairness_category_by_domain_raw.csv", index=False)
    out["fairness_category"].value_counts().rename_axis("fairness_category").reset_index(name="n").to_csv(sdir/"fairness_category_counts_raw.csv", index=False)
    summary={"n_results": int(len(out)), "mean_fairness_index": float(out["fairness_index"].mean()), "sd_fairness_index": float(out["fairness_index"].std()), "n_strong_favored": int((out["fairness_category"]=="strongly_unfair_favored").sum()), "n_strong_unfavored": int((out["fairness_category"]=="strongly_unfair_unfavored").sum())}
    with open(sdir/"fairness_run_summary.json","w") as f: json.dump(summary,f,indent=2)
    print(f"OK: wrote fairness CSV: {args.out_csv} ({len(out):,} rows)")
    print(f"OK: wrote fairness parquet: {args.out_parquet}")

if __name__=="__main__":
    main()
