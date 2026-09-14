#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Iterable, List, Optional

import duckdb
import numpy as np
import pandas as pd


def parquet_expr(path: str) -> str:
    p = Path(path)
    if p.is_dir() or (not p.suffix and not any(ch in path for ch in '*?[')):
        return str(p / '**' / '*.parquet')
    return str(p)


def read_table(path: str, columns: Optional[List[str]] = None) -> pd.DataFrame:
    p = str(path)
    if p.endswith('.csv') or p.endswith('.tsv'):
        sep = '\t' if p.endswith('.tsv') else ','
        return pd.read_csv(p, sep=sep, low_memory=False)
    return pd.read_parquet(p, columns=columns)


def normalize_text_key(x: object) -> str:
    if pd.isna(x):
        return ''
    s = str(x).strip().lower()
    s = re.sub(r'[^a-z0-9]+', ' ', s)
    s = re.sub(r'\s+', ' ', s).strip()
    return s


def normalize_issn(x: object) -> str:
    if pd.isna(x):
        return ''
    s = re.sub(r'[^0-9Xx]', '', str(x))
    if len(s) == 8:
        return s[:4] + '-' + s[4:].upper()
    return str(x).strip().upper()


def make_journal_key(df: pd.DataFrame) -> pd.Series:
    cols = {c.lower(): c for c in df.columns}
    if 'journal_key' in cols:
        return df[cols['journal_key']].map(normalize_text_key)
    keys = pd.Series([''] * len(df), index=df.index, dtype='object')
    for c in ['eissn', 'issn_electronic', 'issn', 'issn_print']:
        if c in cols:
            val = df[cols[c]].map(normalize_issn)
            keys = keys.mask(keys.eq('') & val.ne(''), val)
    for c in ['journal_title', 'journal', 'journal_name', 'journal_iso_abbrev']:
        if c in cols:
            val = df[cols[c]].map(normalize_text_key)
            keys = keys.mask(keys.eq('') & val.ne(''), val)
    return keys


def safe_cols(df: pd.DataFrame, cols: Iterable[str]) -> List[str]:
    return [c for c in cols if c in df.columns]


def require_cols(df: pd.DataFrame, cols: Iterable[str], name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f'{name} is missing required columns: {missing}')


def write_json(obj, path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def robust_z_by_group(df: pd.DataFrame, value: str, group_cols: List[str]) -> pd.Series:
    out = pd.Series(np.nan, index=df.index, dtype='float64')
    for _, idx in df.groupby(group_cols, dropna=False).groups.items():
        x = df.loc[idx, value].astype(float)
        med = np.nanmedian(x)
        mad = np.nanmedian(np.abs(x - med))
        if not np.isfinite(mad) or mad == 0:
            sd = np.nanstd(x)
            denom = sd if np.isfinite(sd) and sd > 0 else 1.0
            out.loc[idx] = (x - np.nanmean(x)) / denom
        else:
            out.loc[idx] = (x - med) / (1.4826 * mad)
    return out


def choose_first_existing(columns: Iterable[str], candidates: Iterable[str], default: Optional[str] = None) -> Optional[str]:
    cols = set(columns)
    for c in candidates:
        if c in cols:
            return c
    return default


def describe_parquet(path: str) -> List[str]:
    con = duckdb.connect()
    expr = parquet_expr(path)
    try:
        return [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM read_parquet('{expr}', hive_partitioning=1) LIMIT 0").fetchall()]
    finally:
        con.close()
