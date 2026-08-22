#!/usr/bin/env python3
"""Create simple audit tables for strongly favored/unfavored PFI outliers."""
from __future__ import annotations
import argparse
from pathlib import Path
import re
import pandas as pd
import pyarrow.parquet as pq


def read_any(path):
    p=Path(path)
    if not p.exists():
        raise SystemExit(f"ERROR: fairness file missing: {p}")
    if p.suffix.lower()=='.parquet':
        return pq.read_table(p).to_pandas()
    return pd.read_csv(p, low_memory=False)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--fairness', required=True)
    ap.add_argument('--out_dir', required=True)
    args=ap.parse_args()
    df=read_any(args.fairness)
    out=Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    if 'fairness_z' in df.columns:
        df['fairness_z']=pd.to_numeric(df['fairness_z'], errors='coerce')
    else:
        df['fairness_z']=0
    if 'fairness_category' not in df.columns:
        df['fairness_category']='fair'
    df.sort_values('fairness_z', ascending=False).head(1000).to_csv(out/'top_favored_outliers.csv', index=False)
    df.sort_values('fairness_z', ascending=True).head(1000).to_csv(out/'top_unfavored_outliers.csv', index=False)
    if 'famous_institution_proxy' in df.columns:
        inst=df.groupby('fairness_category', dropna=False).agg(n=('fairness_category','size'), famous_institution_rate=('famous_institution_proxy','mean')).reset_index()
        inst.to_csv(out/'famous_institution_proxy_by_fairness_category.csv', index=False)
    if 'journal_title' in df.columns:
        top=df[df['fairness_category'].astype(str).str.contains('strongly_unfair_favored', na=False)]
        top.groupby('journal_title', dropna=False).size().reset_index(name='n_strong_favored').sort_values('n_strong_favored', ascending=False).head(200).to_csv(out/'strong_favored_by_journal.csv', index=False)
    print(f"OK: wrote outlier audit tables to {out}")

if __name__=='__main__':
    main()
