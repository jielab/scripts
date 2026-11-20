#!/usr/bin/env python3
"""Create provisional review inputs for bootstrapping the PFI training flow."""


# 🚩 Imports, paths, and shared inputs
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
	return df["pmid"].map(norm_id).combine_first(df["pmcid"].map(norm_id)).combine_first(df["doi"].map(norm_id))


def infer_domain_from_journal(journal: object) -> str:
	s = str(journal or "").lower()
	if any(k in s for k in ["public health", "epidemiol", "environment", "occupational", "health policy"]):
		return "public_health"
	if any(
		k in s for k in ["clinical", "medicine", "surgery", "oncology", "cardiol", "neurol", "pediatric", "infectious"]
	):
		return "clinical_science"
	if any(
		k in s for k in ["cell", "molecular", "genetic", "genomic", "immunol", "microbiol", "neuroscience", "biochem"]
	):
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
		for a, b, c in zip(
			df.get("has_new_discovery_language", 0), df.get("has_replication", 0), wc.ge(1000).astype(int)
		)
	]
	df["new_methods"] = [
		clip5(1 + 2.2 * a + 1.2 * b + 0.8 * c)
		for a, b, c in zip(
			df.get("has_method_innovation", 0),
			df.get("has_external_validation", 0),
			df.get("has_sensitivity_analysis", 0),
		)
	]
	df["new_guidelines"] = [
		clip5(0.5 + 1.5 * a + 1.0 * b + 1.0 * c)
		for a, b, c in zip(
			df.get("has_trial_registration", 0), df.get("has_ethics", 0), df.get("has_data_availability", 0)
		)
	]
	df["big_impact"] = [
		clip5(1 + 1.5 * a + 1.0 * b + 0.5 * c)
		for a, b, c in zip(df.get("has_randomized", 0), df.get("has_meta_analysis", 0), ss_scaled)
	]
	df["big_use"] = [
		clip5(1 + 1.6 * a + 0.9 * b + 0.8 * c)
		for a, b, c in zip(
			df.get("has_large_resource_language", 0),
			df.get("has_open_science_language", 0),
			df.get("has_code_availability", 0),
		)
	]
	df["claim_evidence_gap"] = [
		clip5(2.5 - 0.6 * a - 0.5 * b - 0.4 * c)
		for a, b, c in zip(
			df.get("has_replication", 0), df.get("has_external_validation", 0), df.get("has_sensitivity_analysis", 0)
		)
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
    from review_inputs import main as export_reviews
    export_reviews()


if __name__ == "__main__":
    main()
