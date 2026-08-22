#!/usr/bin/env python3
"""Create weak-supervision pseudo scores for the 1,000-paper PFI review set.

This is a practical fallback when full expert review is unavailable. It is not a
claim that labels equal true paper quality. It creates transparent, auditable,
blinded-content labels using:
  1) blinded full-text section content;
  2) model-free article features, if available;
  3) optional article-level external metrics such as NIH iCite RCR or OpenAlex
     field/year citation percentiles;
  4) optional high/low anchor CSV files.

The script does not use actual journal score or journal name to construct the
paper-score dimensions, preventing direct journal-prestige leakage.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    import pyarrow.dataset as ds
except Exception:
    ds = None


def norm_col(x: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(x).lower())).strip("_")


def read_table(path: str | Path, columns: list[str] | None = None) -> pd.DataFrame:
    path = Path(path)
    if not str(path) or not path.exists():
        return pd.DataFrame()
    if path.is_dir() or path.suffix.lower() == ".parquet":
        if ds is None:
            return pd.DataFrame()
        dataset = ds.dataset(str(path), format="parquet")
        avail = set(dataset.schema.names)
        cols = [c for c in (columns or list(avail)) if c in avail]
        if not cols:
            return pd.DataFrame()
        return dataset.to_table(columns=cols).to_pandas()
    return pd.read_csv(path, dtype=str, low_memory=False)


def numeric(s: pd.Series, default: float = np.nan) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").fillna(default)


def clip01(x: pd.Series | np.ndarray | float) -> pd.Series | np.ndarray | float:
    return np.clip(x, 0.0, 1.0)


def rescale01(s: pd.Series, lo: float | None = None, hi: float | None = None) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    if lo is None:
        lo = float(np.nanpercentile(x, 5)) if x.notna().any() else 0.0
    if hi is None:
        hi = float(np.nanpercentile(x, 95)) if x.notna().any() else 1.0
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return pd.Series(np.zeros(len(s)), index=s.index)
    return pd.Series(np.clip((x - lo) / (hi - lo), 0, 1), index=s.index).fillna(0)


def any_col(df: pd.DataFrame, names: Iterable[str]) -> str | None:
    cols = {norm_col(c): c for c in df.columns}
    for n in names:
        if norm_col(n) in cols:
            return cols[norm_col(n)]
    return None


def text_cols(df: pd.DataFrame) -> list[str]:
    preferred = [
        "title", "abstract", "intro_text", "methods_text", "results_text", "discussion_text",
        "data_availability_text", "software_availability_text", "ethics_text", "trial_registration_text", "blinded_text", "text",
    ]
    out = [c for c in preferred if c in df.columns]
    if out:
        return out
    # Fallback: include short object columns that look like text-bearing fields.
    return [c for c in df.columns if df[c].dtype == object and re.search(r"text|abstract|title|method|result|discussion", c, re.I)]


def contains(series: pd.Series, pattern: str) -> pd.Series:
    return series.str.contains(pattern, case=False, regex=True, na=False).astype(float)


def count_terms(series: pd.Series, pattern: str) -> pd.Series:
    return series.str.count(pattern, flags=re.I).fillna(0).astype(float)


def load_anchor(path: str | Path, label: str) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    a = pd.read_csv(p, dtype=str, low_memory=False)
    a.columns = [norm_col(c) for c in a.columns]
    a["anchor_type"] = label
    return a


def merge_optional(left: pd.DataFrame, right: pd.DataFrame, suffix: str) -> pd.DataFrame:
    if right.empty:
        return left
    right.columns = [norm_col(c) for c in right.columns]
    left.columns = [norm_col(c) for c in left.columns]
    for key in ["pmid", "pmcid", "doi", "article_id"]:
        if key in left.columns and key in right.columns:
            keep = [key] + [c for c in right.columns if c != key]
            return left.merge(right[keep].drop_duplicates(key), on=key, how="left", suffixes=("", suffix))
    return left


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--review_csv", required=True)
    ap.add_argument("--review_key", default="")
    ap.add_argument("--features", default="")
    ap.add_argument("--article_metrics", default="")
    ap.add_argument("--high_anchor_csv", default="")
    ap.add_argument("--low_anchor_csv", default="")
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--audit_csv", required=True)
    ap.add_argument("--mode", default="weak_supervision")
    args = ap.parse_args()

    review = pd.read_csv(args.review_csv, dtype=str, low_memory=False)
    review.columns = [norm_col(c) for c in review.columns]

    key_cols = ["article_id", "pmid", "pmcid", "doi", "year", "journal_key", "journal_title", "publication_type_proxy", "domain"]
    key = read_table(args.review_key, columns=key_cols) if args.review_key else pd.DataFrame()
    if not key.empty:
        key.columns = [norm_col(c) for c in key.columns]
        review = merge_optional(review, key, "_key")

    # Features are optional; because they can be huge, only pull columns likely to be useful.
    feature_cols = [
        "article_id", "pmid", "pmcid", "doi", "n_chars", "text_chars", "abstract_chars", "methods_chars", "results_chars",
        "n_tables", "n_figures", "has_data_availability", "has_code_availability", "has_trial_registration",
        "has_ethics", "has_supplement", "has_retraction_word", "has_correction_word", "has_guideline_word",
    ]
    feat = read_table(args.features, columns=feature_cols) if args.features else pd.DataFrame()
    if not feat.empty:
        feat.columns = [norm_col(c) for c in feat.columns]
        review = merge_optional(review, feat, "_feat")

    metrics = read_table(args.article_metrics) if args.article_metrics else pd.DataFrame()
    if not metrics.empty:
        metrics.columns = [norm_col(c) for c in metrics.columns]
        review = merge_optional(review, metrics, "_metrics")

    high = load_anchor(args.high_anchor_csv, "high_anchor")
    low = load_anchor(args.low_anchor_csv, "low_anchor")
    if not high.empty or not low.empty:
        anchors = pd.concat([x for x in [high, low] if not x.empty], ignore_index=True)
        # Reduce to ids + anchor_type; multiple anchors concatenate.
        for key_id in ["pmid", "pmcid", "doi", "article_id"]:
            if key_id in review.columns and key_id in anchors.columns:
                tmp = anchors[[key_id, "anchor_type"]].dropna().drop_duplicates()
                tmp = tmp.groupby(key_id, as_index=False)["anchor_type"].agg(lambda x: ";".join(sorted(set(map(str, x)))))
                review = review.merge(tmp, on=key_id, how="left")
                break
    if "anchor_type" not in review.columns:
        review["anchor_type"] = ""
    review["anchor_type"] = review["anchor_type"].fillna("")

    tc = text_cols(review)
    if tc:
        text = review[tc].fillna("").agg("\n".join, axis=1).str.slice(0, 120000)
    else:
        text = pd.Series([""] * len(review), index=review.index)
    text_l = text.str.lower()

    text_len = text.str.len().astype(float)
    len_score = rescale01(text_len, lo=2000, hi=30000)

    # External article-level metrics: use only article-level signals, not journal rank.
    rcr_col = any_col(review, ["rcr", "relative_citation_ratio", "nih_rcr"])
    cite_pct_col = any_col(review, ["field_citation_percentile", "citation_percentile", "openalex_citation_percentile", "cited_by_percentile"])
    cited_col = any_col(review, ["cited_by_count", "openalex_cited_by_count", "citation_count", "citations"])
    rcr_score = rescale01(numeric(review[rcr_col], np.nan), lo=0, hi=5) if rcr_col else pd.Series(0, index=review.index)
    cite_pct_score = numeric(review[cite_pct_col], np.nan).clip(0, 100).fillna(0) / 100 if cite_pct_col else pd.Series(0, index=review.index)
    cited_score = rescale01(np.log1p(numeric(review[cited_col], 0)), lo=0, hi=math.log1p(500)) if cited_col else pd.Series(0, index=review.index)
    article_influence = pd.concat([rcr_score, cite_pct_score, cited_score], axis=1).max(axis=1).fillna(0)

    is_retracted = contains(text_l, r"\bretract(?:ed|ion)\b|expression of concern")
    for c in ["is_retracted", "retracted", "has_retraction_word", "expression_of_concern"]:
        if c in review.columns:
            is_retracted = np.maximum(is_retracted, numeric(review[c], 0).clip(0, 1))
    has_correction = contains(text_l, r"\bcorrection\b|corrigendum|erratum")
    for c in ["has_correction_word", "major_correction", "correction"]:
        if c in review.columns:
            has_correction = np.maximum(has_correction, numeric(review[c], 0).clip(0, 1))

    has_trial = contains(text_l, r"clinical trial|randomi[sz]ed|pragmatic trial|trial registration|nct\d{8}")
    has_cohort = contains(text_l, r"prospective cohort|population-based|longitudinal|biobank|registry")
    has_validation = contains(text_l, r"external validation|independent validation|replication cohort|validation cohort")
    has_meta = contains(text_l, r"systematic review|meta-analysis|meta analysis")
    has_guideline = contains(text_l, r"guideline|consensus statement|recommendation")
    has_data = contains(text_l, r"data availability|available at|deposited|repository|accession|zenodo|dryad|figshare|github|geo accession|dbgap|bioproject")
    has_code = contains(text_l, r"code availability|source code|github|gitlab|software package|r package|python package|cran|bioconductor|workflow|container|docker")
    has_software = contains(text_l, r"software|algorithm|toolkit|pipeline|package|database|web server|app|application")
    has_method = contains(text_l, r"novel method|new method|we developed|we propose|framework|assay|protocol|model architecture|statistical method|machine learning|deep learning")
    has_new_data = contains(text_l, r"new dataset|dataset|cohort|registry|trial|biobank|whole[- ]genome|single[- ]cell|proteom|metabolom|multi[- ]omics|database|resource")
    has_results = contains(text_l, r"we found|we observed|associated with|hazard ratio|odds ratio|risk ratio|confidence interval|p[ -]?value|significant|effect size")
    overclaim = contains(text_l, r"breakthrough|paradigm shift|revolutionary|definitive proof|cure|game[- ]changing|unprecedented")

    methods_chars = numeric(review["methods_chars"], np.nan) if "methods_chars" in review.columns else count_terms(text_l, r"method|statistical|analysis|model|protocol")
    results_chars = numeric(review["results_chars"], np.nan) if "results_chars" in review.columns else count_terms(text_l, r"result|figure|table|confidence interval|p[ -]?value")
    methods_score = rescale01(methods_chars, lo=100, hi=6000)
    results_score = rescale01(results_chars, lo=100, hi=6000)

    new_methods = 100 * clip01(0.25 + 0.20 * methods_score + 0.35 * has_method + 0.15 * has_software + 0.05 * article_influence)
    new_data = 100 * clip01(0.25 + 0.25 * len_score + 0.35 * has_new_data + 0.10 * has_data + 0.05 * article_influence)
    new_software = 100 * clip01(0.15 + 0.55 * has_code + 0.25 * has_software + 0.05 * article_influence)
    new_results = 100 * clip01(0.30 + 0.25 * results_score + 0.15 * has_results + 0.15 * article_influence + 0.10 * has_validation + 0.05 * has_trial)
    clinical_public_health_value = 100 * clip01(0.25 + 0.20 * has_trial + 0.20 * has_guideline + 0.10 * has_cohort + 0.15 * has_meta + 0.10 * article_influence)
    evidence_strength = 100 * clip01(0.35 + 0.15 * methods_score + 0.15 * results_score + 0.15 * has_trial + 0.10 * has_validation + 0.10 * has_cohort - 0.35 * is_retracted - 0.10 * has_correction)
    reproducibility = 100 * clip01(0.25 + 0.30 * has_data + 0.25 * has_code + 0.10 * has_trial + 0.10 * has_validation)
    claim_evidence_gap = 100 * clip01(0.10 + 0.35 * overclaim + 0.45 * is_retracted + 0.15 * has_correction - 0.15 * evidence_strength / 100 - 0.10 * reproducibility / 100)

    score = (
        0.18 * new_results
        + 0.16 * new_methods
        + 0.14 * new_data
        + 0.12 * new_software
        + 0.16 * clinical_public_health_value
        + 0.14 * evidence_strength
        + 0.10 * reproducibility
        - 0.10 * claim_evidence_gap
    )

    high_anchor = review["anchor_type"].str.contains("high", case=False, na=False).astype(float)
    low_anchor = review["anchor_type"].str.contains("low", case=False, na=False).astype(float)
    score = score + 8 * high_anchor - 12 * low_anchor - 25 * is_retracted
    score = np.clip(score, 0, 100)

    confidence = np.where((high_anchor == 1) | (low_anchor == 1) | (is_retracted > 0), "high_anchor", "")
    confidence = np.where((confidence == "") & ((article_influence >= 0.75) | (has_data + has_code + has_trial + has_validation >= 2)), "medium_proxy", confidence)
    confidence = np.where(confidence == "", "low_proxy", confidence)

    out = review.copy()
    out["new_results"] = np.round(new_results, 1)
    out["new_methods"] = np.round(new_methods, 1)
    out["new_data"] = np.round(new_data, 1)
    out["new_software"] = np.round(new_software, 1)
    out["clinical_public_health_value"] = np.round(clinical_public_health_value, 1)
    out["evidence_strength"] = np.round(evidence_strength, 1)
    out["reproducibility"] = np.round(reproducibility, 1)
    out["claim_evidence_gap"] = np.round(claim_evidence_gap, 1)
    out["paper_score_formula"] = np.round(score, 1)
    out["paper_score"] = np.round(score, 1)
    out["label_source"] = "weak_supervision_auto"
    out["label_confidence"] = confidence
    out["article_influence_proxy_0_1"] = np.round(article_influence, 3)
    out["anchor_type"] = out["anchor_type"].fillna("")

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)

    audit = {
        "mode": args.mode,
        "n": int(len(out)),
        "mean_paper_score": float(np.nanmean(score)),
        "sd_paper_score": float(np.nanstd(score)),
        "min_paper_score": float(np.nanmin(score)),
        "max_paper_score": float(np.nanmax(score)),
        "n_high_anchor": int(high_anchor.sum()),
        "n_low_anchor": int(low_anchor.sum()),
        "n_retracted_or_eoc_proxy": int((is_retracted > 0).sum()),
        "n_with_external_influence_signal": int((article_influence > 0).sum()),
        "label_confidence_counts": json.dumps(pd.Series(confidence).value_counts().to_dict(), ensure_ascii=False),
        "warning": "Weak-supervision labels are audit/proxy labels, not a substitute for independent expert truth labels.",
    }
    Path(args.audit_csv).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([audit]).to_csv(args.audit_csv, index=False)
    print(f"OK: wrote pseudo-scored review CSV to {args.out_csv}")
    print(f"OK: wrote audit CSV to {args.audit_csv}")
    print(json.dumps(audit, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
