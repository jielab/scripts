#!/usr/bin/env python3
"""Build phenotype-free GRID evolutionary score modules from external annotations.

Input weights are additive effect-allele weights; age labels never change their
sign. Genome coordinates must already be normalized to GRCh37. Output is a
temporary exchange format for PLINK and the GRID fitting stage, not a table of
individual phenotypes or predictions.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import urllib.request

import numpy as np
import pandas as pd


# 🚩 Public data and fixed, explicitly non-demographic feature boundaries
POPULATIONS = ("EUR", "AFR", "EAS", "SAS")
AGE_BINS = ("young", "middle", "old", "uncertain", "unknown")
GEVA_URL = "https://human.genome.dating/bulk/atlas.chr{chromosome}.csv.gz"
GEVA_DOWNLOAD_PAGE = "https://human.genome.dating/download/index"
GEVA_FAQ = "https://human.genome.dating/info/faq"
GEVA_MD5 = {
	1: "573d944fd9c0820cda5e4f3e69f0436f", 2: "64928ae5acdab7f2112d77d72e4e05e7",
	3: "3b24b89f15f276761963dc0931a5ddc3", 4: "154db58a8ae87ca32ecbadd606323c8a",
	5: "4655f24184ca50cf9f41e9e0f41758ee", 6: "940f27c8ffc85bc19ea55c10ec435fdc",
	7: "012eb4f227853774ec462c1519f85342", 8: "e30a37a056cc175b944473114e1d2819",
	9: "afa2425bfc8f8d059e2a6d51f4f695f1", 10: "1e9020a64625962d63761a9ef0fbd32a",
	11: "5631f7fc07c8aecb264299fc66d29658", 12: "38ba8cae9809575fad66534ade7f60d5",
	13: "d0bdcf97613a784b6b7084be7ea1ab08", 14: "73477526e6661505c5e4cea99640620e",
	15: "69323fde9d17a3ba86fa1fd9c51cd504", 16: "f6ae0a54cc476ce5a9342f5eed45168d",
	17: "684dbce64241b16fffe97b5c44e1204a", 18: "9336c9832585688a9702fe2934ebc0e3",
	19: "8a7290197378a66c4befe0b38bb7b088", 20: "586c2bdc03742a4b1af602b412d0310c",
	21: "ec493b5518a6b68f29c8d658e2280703", 22: "0b97e098875a54adffcfabe2942358f2",
}
GEVA_AGE_COLUMNS = [
	"VariantID", "Chromosome", "Position", "AlleleRef", "AlleleAlt", "AlleleAnc",
	"DataSource", "AgeMode_Jnt", "AgeCI95Lower_Jnt", "AgeCI95Upper_Jnt", "QualScore_Jnt",
]
AGE_FIELDS = ["AGE_GEN", "AGE_LO", "AGE_HI", "AGE_QUAL", "AGE_SOURCE", "AGE_METHOD",
	"AGE_UNCERTAINTY", "AGE_REF", "AGE_ALT", "ANC"]
MISSING = {"", ".", "NA", "NAN", "NULL", "NONE"}


# 🚩 Input normalization; no strand or genome-build guessing
def chromosomes(value):
	answer = []
	for token in str(value).split(","):
		token = token.strip()
		if "-" in token:
			start, end = map(int, token.split("-", 1))
			if start > end:
				raise ValueError("Chromosome range must be ascending.")
			answer.extend(range(start, end + 1))
		elif token:
			answer.append(int(token))
	answer = sorted(set(answer))
	if not answer or any(x < 1 or x > 22 for x in answer):
		raise ValueError("Only human autosomes 1-22 are supported.")
	return answer


def population_list(value):
	answer = [x.strip().upper() for x in str(value).split(",") if x.strip()]
	if not answer or len(answer) != len(set(answer)) or any(x not in POPULATIONS for x in answer):
		raise ValueError("Populations must be unique members of EUR,AFR,EAS,SAS.")
	return answer


def require_columns(frame, names, label):
	missing = sorted(set(names) - set(frame.columns))
	if missing:
		raise ValueError(f"{label}: missing columns {','.join(missing)}")


def clean_text(series):
	result = series.fillna("").astype(str).str.strip()
	return result.mask(result.str.upper().isin(MISSING), "")


def numbers(series, label, allow_missing=True):
	clean = clean_text(series)
	result = pd.to_numeric(clean.replace("", np.nan), errors="coerce")
	bad = ((clean != "") & result.isna()) | np.isinf(result.to_numpy(dtype=float))
	if bad.any():
		raise ValueError(f"{label}: nonnumeric or infinite values, e.g. {clean.loc[bad].iloc[0]!r}")
	if not allow_missing and result.isna().any():
		raise ValueError(f"{label}: missing values are not permitted.")
	return result.astype(float)


def normalize_variants(frame, ref="REF", alt="ALT", label="variants"):
	frame = frame.copy()
	require_columns(frame, ["CHR", "BP", ref, alt], label)
	if "BUILD" in frame:
		build = clean_text(frame["BUILD"]).str.lower()
		if not build.isin(["grch37", "hg19", "37"]).all():
			raise ValueError(f"{label}: BUILD must explicitly be GRCh37, hg19, or 37.")
	chr_text = clean_text(frame["CHR"]).str.replace(r"^chr", "", regex=True)
	chrom = numbers(chr_text, f"{label}.CHR", False)
	bp = numbers(frame["BP"], f"{label}.BP", False)
	if ((chrom % 1 != 0) | (chrom < 1) | (chrom > 22)).any():
		raise ValueError(f"{label}: only chromosome integers 1-22 are accepted.")
	if ((bp % 1 != 0) | (bp < 1)).any():
		raise ValueError(f"{label}: BP must be a positive 1-based integer.")
	frame["CHR"] = chrom.astype(np.int16)
	frame["BP"] = bp.astype(np.int64)
	for column in [ref, alt]:
		frame[column] = clean_text(frame[column]).str.upper()
		if not frame[column].str.fullmatch("[ACGT]").all():
			raise ValueError(f"{label}.{column}: biallelic A/C/G/T SNPs only; normalize other variants upstream.")
	if (frame[ref] == frame[alt]).any():
		raise ValueError(f"{label}: the two alleles must differ.")
	ref_values, alt_values = frame[ref].to_numpy(dtype=str), frame[alt].to_numpy(dtype=str)
	pair_lo = np.where(ref_values < alt_values, ref_values, alt_values)
	pair_hi = np.where(ref_values < alt_values, alt_values, ref_values)
	frame["_key"] = frame["CHR"].astype(str) + ":" + frame["BP"].astype(str) + ":" + pair_lo + ":" + pair_hi
	return frame


def load_weights(path, pops, selected_chromosomes):
	frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
	frame.columns = frame.columns.str.strip()
	require_columns(frame, ["CHR", "BP", "SNP", "A1", "A2"], "weights")
	frame = normalize_variants(frame, "A1", "A2", "weights")
	frame = frame.loc[frame["CHR"].isin(selected_chromosomes)].copy()
	if frame.empty:
		raise ValueError("No weights remain on the selected chromosomes.")
	frame["SNP"] = clean_text(frame["SNP"])
	if (frame["SNP"] == "").any() or frame["SNP"].str.contains(r"\s", regex=True).any():
		raise ValueError("SNP identifiers must be nonempty and contain no whitespace.")
	if frame["SNP"].duplicated().any():
		raise ValueError("Duplicate SNP IDs in weights would make PLINK scoring ambiguous.")
	if frame["_key"].duplicated().any():
		raise ValueError("Duplicate position/allele pairs in weights.")
	present = [f"beta_{pop}" for pop in pops if f"beta_{pop}" in frame]
	if not present:
		raise ValueError("At least one requested beta_<POP> column is required.")
	for pop in pops:
		column = f"beta_{pop}"
		frame[column] = numbers(frame[column], column) if column in frame else np.nan
	if frame[[f"beta_{pop}" for pop in pops]].isna().all(axis=1).any():
		raise ValueError("Some SNP rows have no effect estimate in any requested population.")
	return frame.sort_values(["CHR", "BP", "SNP"]).reset_index(drop=True)


def read_chunks(path, size, sep="\t"):
	for frame in pd.read_csv(path, sep=sep, comment="#", dtype=str, keep_default_na=False,
		chunksize=size, skipinitialspace=True):
		frame.columns = frame.columns.str.strip()
		yield frame


def empty_annotation(index):
	frame = pd.DataFrame(index=index)
	for column in AGE_FIELDS:
		frame[column] = np.nan if column in ["AGE_GEN", "AGE_LO", "AGE_HI", "AGE_QUAL"] else ""
	return frame


# 🚩 Native GEVA reader: exact coordinate/allele matching and explicit provenance
def load_geva(weights, directory, chunk_size, selected_chromosomes, qc):
	annotation = empty_annotation(weights.index)
	rank = np.full(len(weights), 99, dtype=np.int16)
	lookup = dict(zip(weights["_key"], weights.index))
	source_rank = {"Combined": 0, "TGP": 1, "SGDP": 2}
	for chromosome in selected_chromosomes:
		if not (weights["CHR"] == chromosome).any():
			continue
		path = Path(directory) / f"atlas.chr{chromosome}.csv.gz"
		if not path.is_file():
			raise ValueError(f"Missing GEVA chromosome file: {path}. Use the explicit download command first.")
		seen = set()
		matched = 0
		for raw in read_chunks(path, chunk_size, ","):
			require_columns(raw, GEVA_AGE_COLUMNS, str(path))
			raw = raw.rename(columns={"Chromosome": "CHR", "Position": "BP",
				"AlleleRef": "REF", "AlleleAlt": "ALT"})
			frame = normalize_variants(raw, label=str(path))
			if not (frame["CHR"] == chromosome).all():
				raise ValueError(f"{path}: contains a different chromosome.")
			frame = frame.loc[frame["_key"].isin(lookup)].copy()
			if frame.empty:
				continue
			for column in ["AgeMode_Jnt", "AgeCI95Lower_Jnt", "AgeCI95Upper_Jnt", "QualScore_Jnt"]:
				frame[column] = numbers(frame[column], f"GEVA.{column}")
			frame["DataSource"] = frame["DataSource"].str.strip()
			if not frame["DataSource"].isin(source_rank).all():
				raise ValueError("Unrecognized GEVA DataSource; expected Combined,TGP,SGDP.")
			identities = list(zip(frame["_key"], frame["DataSource"]))
			if len(set(identities)) != len(identities) or seen.intersection(identities):
				raise ValueError("GEVA contains a duplicated SNP/source row.")
			seen.update(identities)
			matched += len(frame)
			frame["_rank"] = frame["DataSource"].map(source_rank)
			frame = frame.sort_values("_rank").drop_duplicates("_key", keep="first")
			frame.index = frame["_key"].map(lookup).to_numpy()
			frame = frame.loc[frame["_rank"].to_numpy() < rank[frame.index]].copy()
			if frame.empty:
				continue
			rank[frame.index] = frame["_rank"].to_numpy()
			values = {
				"AGE_GEN": frame["AgeMode_Jnt"], "AGE_LO": frame["AgeCI95Lower_Jnt"],
				"AGE_HI": frame["AgeCI95Upper_Jnt"], "AGE_QUAL": frame["QualScore_Jnt"],
				"AGE_SOURCE": "GEVA_Atlas:" + frame["DataSource"], "AGE_METHOD": "GEVA_joint_clock_mode",
				"AGE_UNCERTAINTY": "composite_posterior_95pct_bounds_may_be_overconfident",
				"AGE_REF": frame["REF"], "AGE_ALT": frame["ALT"], "ANC": frame["AlleleAnc"].str.strip().str.upper(),
			}
			for column, values_column in values.items():
				annotation.loc[frame.index, column] = values_column
		qc.append(("matching", f"geva_chr{chromosome}_matched_source_rows", "", matched,
			"Combined>TGP>SGDP; exact position and alleles, no strand complement."))
	return annotation


# 🚩 Canonical annotations: age, selection, frequency and proxy stay separate
def load_canonical(weights, path, pops, chunk_size, include_proxy, qc):
	annotation = empty_annotation(weights.index)
	lookup = dict(zip(weights["_key"], weights.index))
	seen = set()
	available = set()
	for raw in read_chunks(path, chunk_size):
		frame = normalize_variants(raw, label="canonical annotation")
		available.update(frame.columns)
		frame = frame.loc[frame["_key"].isin(lookup)].copy()
		if frame.empty:
			continue
		for key in frame["_key"]:
			if key in seen:
				raise ValueError(f"Duplicate canonical annotation: {key}")
			seen.add(key)
		frame.index = [lookup[key] for key in frame["_key"]]
		if "AGE_GEN" in frame:
			for column in ["AGE_GEN", "AGE_LO", "AGE_HI", "AGE_QUAL"]:
				frame[column] = numbers(frame[column], f"annotation.{column}") if column in frame else np.nan
			for column in ["AGE_SOURCE", "AGE_METHOD", "AGE_UNCERTAINTY", "ANC"]:
				frame[column] = clean_text(frame[column]) if column in frame else ""
			has_age = frame["AGE_GEN"].notna()
			for column in ["AGE_SOURCE", "AGE_METHOD", "AGE_UNCERTAINTY"]:
				if (has_age & (frame[column] == "")).any():
					raise ValueError(f"Canonical ages require nonempty {column}; unknown uncertainty must be explicit.")
			frame["AGE_REF"] = frame["REF"]
			frame["AGE_ALT"] = frame["ALT"]
			frame["ANC"] = frame["ANC"].str.upper()
			annotation.loc[frame.index, AGE_FIELDS] = frame[AGE_FIELDS]
		selection_columns = [f"SEL_LOG10P_{pop}" for pop in pops if f"SEL_LOG10P_{pop}" in frame]
		for column in selection_columns:
			values = numbers(frame[column], column)
			if (values.dropna() > 0).any():
				raise ValueError(f"{column} is log10(P), not -log10(P), and must be <=0.")
			annotation.loc[frame.index, column] = values
		if selection_columns:
			has_selection = annotation.loc[frame.index, selection_columns].notna().any(axis=1)
			for column in ["SEL_SOURCE", "SEL_METHOD"]:
				values = clean_text(frame[column]) if column in frame else pd.Series("", index=frame.index)
				if (has_selection & (values == "")).any():
					raise ValueError(f"Selection statistics require {column}.")
				annotation.loc[frame.index, column] = values
		frequency_columns = [f"AF_{pop}" for pop in pops if f"AF_{pop}" in frame]
		for column in frequency_columns:
			values = numbers(frame[column], column)
			if ((values.dropna() < 0) | (values.dropna() > 1)).any():
				raise ValueError(f"{column} must be an ALT-allele frequency in [0,1].")
			annotation.loc[frame.index, column] = values
		if frequency_columns:
			has_frequency = annotation.loc[frame.index, frequency_columns].notna().any(axis=1)
			values = clean_text(frame["FREQ_SOURCE"]) if "FREQ_SOURCE" in frame else pd.Series("", index=frame.index)
			if (has_frequency & (values == "")).any():
				raise ValueError("Frequency annotations require FREQ_SOURCE.")
			annotation.loc[frame.index, "FREQ_SOURCE"] = values
			annotation.loc[frame.index, "FREQ_ALT"] = frame["ALT"]
		if include_proxy:
			proxy_columns = [x for x in frame if x.startswith("PROXY_") and x not in ["PROXY_SOURCE", "PROXY_METHOD"]]
			for column in proxy_columns:
				if not re.fullmatch(r"PROXY_[A-Za-z][A-Za-z0-9_]*", column):
					raise ValueError(f"Unsafe proxy feature name: {column}")
				values = numbers(frame[column], column)
				if (values.dropna() < 0).any():
					raise ValueError("Proxy multipliers must be nonnegative to preserve effect-allele signs.")
				annotation.loc[frame.index, column] = values
			if proxy_columns:
				has_proxy = annotation.loc[frame.index, proxy_columns].notna().any(axis=1)
				for column in ["PROXY_SOURCE", "PROXY_METHOD"]:
					values = clean_text(frame[column]) if column in frame else pd.Series("", index=frame.index)
					if (has_proxy & (values == "")).any():
						raise ValueError(f"Proxy annotations require {column}.")
					annotation.loc[frame.index, column] = values
	qc.append(("matching", "canonical_matched_variants", "", len(seen),
		"Canonical age, selection and frequency can have different sources."))
	return annotation


def combine_annotations(geva, canonical, preference):
	answer = geva.copy()
	other_columns = [x for x in canonical if x not in AGE_FIELDS]
	for column in other_columns:
		answer[column] = canonical[column]
	canonical_age = canonical["AGE_GEN"].notna()
	if preference == "geva":
		canonical_age &= geva["AGE_GEN"].isna()
	answer.loc[canonical_age, AGE_FIELDS] = canonical.loc[canonical_age, AGE_FIELDS]
	return answer


# 🚩 Age labels: no benign/pathogenic inference and no beta sign changes
def classify_age(annotation, cutoff1, cutoff2, min_quality):
	answer = annotation.copy()
	for column in ["AGE_GEN", "AGE_LO", "AGE_HI", "AGE_QUAL"]:
		answer[column] = pd.to_numeric(answer[column], errors="raise").astype(float)
	age, lower, upper, quality = [answer[column].to_numpy() for column in ["AGE_GEN", "AGE_LO", "AGE_HI", "AGE_QUAL"]]
	if np.any(np.isfinite(age) & (age <= 0)):
		raise ValueError("Known AGE_GEN estimates must be >0 generations.")
	if np.any(np.isfinite(quality) & ((quality < 0) | (quality > 1))):
		raise ValueError("AGE_QUAL must be in [0,1].")
	if np.any((np.isfinite(lower) & (lower < 0)) | (np.isfinite(upper) & (upper < 0))):
		raise ValueError("Age bounds cannot be negative.")
	if np.any(np.isfinite(lower) & np.isfinite(upper) & (lower > upper)):
		raise ValueError("Age lower bound exceeds upper bound.")
	oriented = (answer["ANC"].fillna("") == answer["AGE_REF"].fillna("")) & (answer["AGE_REF"].fillna("") != "")
	known = np.isfinite(age) & oriented.to_numpy()
	qualified = known & np.isfinite(quality) & (quality >= min_quality)
	bounds = np.isfinite(lower) & np.isfinite(upper)
	compatible_bounds = bounds & (lower <= age) & (age <= upper)
	label = np.full(len(answer), "unknown", dtype=object)
	label[known] = "uncertain"
	label[qualified & compatible_bounds & (upper < cutoff1)] = "young"
	label[qualified & compatible_bounds & (lower >= cutoff1) & (upper < cutoff2)] = "middle"
	label[qualified & compatible_bounds & (lower >= cutoff2)] = "old"
	status = np.full(len(answer), "no_age", dtype=object)
	status[np.isfinite(age) & ~oriented.to_numpy()] = "ancestral_state_unconfirmed"
	status[known] = "missing_quality"
	status[known & np.isfinite(quality) & (quality < min_quality)] = "low_quality"
	status[qualified & ~bounds] = "missing_bounds"
	status[qualified & bounds] = "boundary_uncertain"
	status[qualified & bounds & ~compatible_bounds] = "point_outside_bounds"
	status[np.isin(label, ["young", "middle", "old"])] = "assigned"
	point_label = np.full(len(answer), "unknown", dtype=object)
	point_label[known & (age < cutoff1)] = "young"
	point_label[known & (age >= cutoff1) & (age < cutoff2)] = "middle"
	point_label[known & (age >= cutoff2)] = "old"
	answer["AGE_RAW_GEN"] = age
	answer.loc[~oriented, ["AGE_GEN", "AGE_LO", "AGE_HI"]] = np.nan
	answer["age_bin"] = label
	answer["age_point_bin"] = point_label
	answer["age_status"] = status
	answer["age_usable"] = qualified
	return answer


def add_frequency_features(weights, annotation, pops, cutoffs):
	answer = annotation.copy()
	columns = []
	for pop in pops:
		name = f"AF_{pop}"
		if name in answer:
			columns.append(name)
			is_effect_alt = weights["A1"] == answer["FREQ_ALT"].fillna("")
			answer[f"EAF_{pop}"] = answer[name].where(is_effect_alt, 1 - answer[name])
	if not columns:
		answer["freq_bin"] = "missing"
		return answer, False
	freq = answer[columns].astype(float)
	known_count = freq.notna().sum(axis=1)
	mean_frequency = freq.mean(axis=1)
	answer["reference_maf"] = np.minimum(mean_frequency, 1 - mean_frequency)
	answer["frequency_range"] = (freq.max(axis=1) - freq.min(axis=1)).where(known_count >= 2)
	bins = np.searchsorted(np.array([0.01, 0.05, 0.10, 0.25]), answer["reference_maf"].fillna(0.5), side="right")
	answer["freq_bin"] = pd.Series([f"maf{int(x)}" for x in bins], index=answer.index).where(known_count >= 1, "missing")
	answer["diff_bin"] = "unknown"
	valid = answer["frequency_range"].notna()
	answer.loc[valid & (answer["frequency_range"] < cutoffs[0]), "diff_bin"] = "low"
	answer.loc[valid & (answer["frequency_range"] >= cutoffs[0]) & (answer["frequency_range"] < cutoffs[1]), "diff_bin"] = "middle"
	answer.loc[valid & (answer["frequency_range"] >= cutoffs[1]), "diff_bin"] = "high"
	return answer, True


def permute_age(weights, annotation, seed):
	answer = annotation.copy()
	rng = np.random.default_rng(seed)
	labels = answer["age_bin"].to_numpy(dtype=object)
	permuted = labels.copy()
	# Missing/uncertain annotations have exactly the same SNP support in both
	# arms. Only qualified age-bin membership is randomized; otherwise the
	# control would also destroy annotation coverage and uncertainty patterns.
	eligible = np.isin(labels, AGE_BINS[:3])
	strata = pd.DataFrame({"CHR": weights["CHR"], "freq_bin": answer["freq_bin"]})
	for indices in strata.groupby(["CHR", "freq_bin"], sort=True).groups.values():
		idx = np.asarray(list(indices), dtype=np.int64)
		qualified = idx[eligible[idx]]
		if len(qualified) > 1:
			permuted[qualified] = rng.permutation(labels[qualified])
		if sorted(permuted[idx].tolist()) != sorted(labels[idx].tolist()):
			raise AssertionError("Permutation did not preserve age-module counts within strata.")
	if not np.array_equal(permuted[~eligible], labels[~eligible]):
		raise AssertionError("Permutation changed missing/uncertain age annotation positions.")
	answer["perm_age_bin"] = permuted
	answer["perm_age_eligible"] = eligible
	answer["perm_stratum"] = weights["CHR"].astype(str) + ":" + answer["freq_bin"].astype(str)
	return answer


# 🚩 Module definitions and phenotype-free PLINK weight matrices
def module_specs(annotation, pops, selection_log10_threshold, use_frequency, include_proxy):
	modules = []
	def append(name, group, pop, multiplier, description):
		modules.append({"module": name, "group": group, "population": pop,
			"description": description, "multiplier": np.asarray(multiplier, dtype=float)})
	for pop in pops:
		append(f"total_{pop}", "total", pop, np.ones(len(annotation)), "Original additive population weight.")
	for group, prefix, column in [("age", "age", "age_bin"), ("permuted_age", "perm_age", "perm_age_bin")]:
		for pop in pops:
			for label in AGE_BINS:
				append(f"{prefix}_{label}_{pop}", group, pop, annotation[column] == label,
					"Inferred derived-allele age partition; fixed feature boundaries, not population split dates." if group == "age"
					else "Fixed-seed control: qualified ages permuted within CHR/reference MAF; uncertain/unknown SNPs stay fixed. Not a permutation P value.")
	for pop in pops:
		column = f"SEL_LOG10P_{pop}"
		if column not in annotation:
			continue
		values = annotation[column].astype(float)
		for label, mask in [
			("detected", values.notna() & (values <= selection_log10_threshold)),
			("not_detected", values.notna() & (values > selection_log10_threshold)),
			("unknown", values.isna()),
		]:
			append(f"sel_{label}_{pop}", "selection", pop, mask,
				"Selection evidence as defined by the declared source; no disease-effect direction implied.")
	if use_frequency:
		for pop in pops:
			for label in ["low", "middle", "high", "unknown"]:
				append(f"diff_{label}_{pop}", "differentiation", pop, annotation["diff_bin"] == label,
					"Across-source-population allele-frequency range; not FST, age, or a selection test.")
	if include_proxy:
		columns = sorted(x for x in annotation if x.startswith("PROXY_") and x not in ["PROXY_SOURCE", "PROXY_METHOD"])
		for column in columns:
			name = column[6:].lower()
			for pop in pops:
				append(f"proxy_{name}_{pop}", "proxy", pop, annotation[column].fillna(0),
					"Nonnegative externally supplied proxy multiplier, not validated mutation age.")
				append(f"proxy_{name}_missing_{pop}", "proxy", pop, annotation[column].isna(),
					"Missingness indicator for the named proxy.")
	names = [x["module"] for x in modules]
	if len(names) != len(set(names)):
		raise ValueError("Module names collide after normalization.")
	return modules


def emit_outputs(weights, annotation, modules, out_dir, mode, qc, force):
	out_dir = Path(out_dir)
	out_dir.parent.mkdir(parents=True, exist_ok=True)
	owned_names = ["modules.tsv.gz", "variants.tsv.gz", "qc.tsv"]
	old = [out_dir / name for name in owned_names if (out_dir / name).exists()]
	old += list(out_dir.glob("chr*.weights.tsv.gz")) if out_dir.is_dir() else []
	if old and not force:
		raise ValueError(f"Output files already exist in {out_dir}; use --force to replace this builder's files.")
	with tempfile.TemporaryDirectory(prefix="grid-evolution-", dir=out_dir.parent) as temporary:
		stage = Path(temporary)
		dictionary = []
		for module in modules:
			beta = weights[f"beta_{module['population']}"].fillna(0).to_numpy(dtype=float)
			multiplier = module["multiplier"]
			dictionary.append({k: v for k, v in module.items() if k != "multiplier"} | {
				"n_variants": int(np.count_nonzero(multiplier)),
				"n_weighted_variants": int(np.count_nonzero(beta * multiplier)),
				"mode": mode,
			})
		pd.DataFrame(dictionary).to_csv(stage / "modules.tsv.gz", sep="\t", index=False)
		variants = pd.concat([weights.drop(columns="_key"), annotation], axis=1)
		variants["evolution_mode"] = mode
		variants.to_csv(stage / "variants.tsv.gz", sep="\t", index=False, na_rep="NA")
		for chromosome, sub in weights.groupby("CHR", sort=True):
			idx = sub.index.to_numpy()
			matrix = pd.DataFrame({"SNP": sub["SNP"].to_numpy(), "A1": sub["A1"].to_numpy()})
			for module in modules:
				beta = sub[f"beta_{module['population']}"].fillna(0).to_numpy(dtype=float)
				matrix[module["module"]] = beta * module["multiplier"][idx]
			for pop in sorted(set(x["population"] for x in modules)):
				for prefix in ["age", "perm_age"]:
					partition = matrix[[f"{prefix}_{label}_{pop}" for label in AGE_BINS]].sum(axis=1).to_numpy()
					if not np.array_equal(partition, matrix[f"total_{pop}"].to_numpy()):
						raise AssertionError(f"chr{chromosome}: {prefix}_{pop} does not reconstruct original total.")
			matrix.to_csv(stage / f"chr{chromosome}.weights.tsv.gz", sep="\t", index=False, float_format="%.17g")
		pd.DataFrame(qc, columns=["category", "metric", "population", "value", "note"]).to_csv(
			stage / "qc.tsv", sep="\t", index=False)
		out_dir.mkdir(parents=True, exist_ok=True)
		new_names = {x.name for x in stage.iterdir()}
		for path in stage.iterdir():
			os.replace(path, out_dir / path.name)
		for path in old:
			if path.name not in new_names:
				path.unlink()
	return {"mode": mode, "n_variants": len(weights), "n_modules": len(modules),
		"chromosomes": sorted(int(x) for x in weights["CHR"].unique()),
		"score_col_nums": f"3-{len(modules) + 2}", "out_dir": str(out_dir)}


def build(args):
	if args.build != "GRCh37":
		raise ValueError("This builder requires GRCh37 coordinates; it never performs an implicit liftover.")
	if not (0 <= args.min_age_quality <= 1):
		raise ValueError("--min-age-quality must be in [0,1].")
	if not (0 < args.age_cutoffs[0] < args.age_cutoffs[1]):
		raise ValueError("--age-cutoffs must be two ascending positive generation counts.")
	if not (0 < args.diff_cutoffs[0] < args.diff_cutoffs[1] <= 1):
		raise ValueError("--diff-cutoffs must be ascending values in (0,1].")
	if not (0 < args.selection_p <= 1) or args.chunk_size < 1:
		raise ValueError("--selection-p must be in (0,1], and --chunk-size positive.")
	pops = population_list(args.populations)
	selected_chromosomes = chromosomes(args.chromosomes)
	weights = load_weights(args.weights, pops, selected_chromosomes)
	qc = [
		("input", "variants", "", len(weights), "No target phenotype or split label is read."),
		("definition", "build", "", "GRCh37", "1-based exact coordinate plus alleles; no strand complement."),
		("definition", "age_cutoffs_generations", "", ",".join(map(str, args.age_cutoffs)), "Feature boundaries, not population separation dates."),
		("definition", "permutation_seed", "", args.seed, "One fixed negative control, not a permutation P value; qualified age labels only, within CHR/reference MAF."),
		("definition", "age_interpretation", "", "inferred_derived_allele_origin_generations", "Estimated from genetic data; not directly observed historical dates or AFR/EUR/EAS/SAS divergence dates."),
	]
	geva = empty_annotation(weights.index)
	canonical = empty_annotation(weights.index)
	if args.geva_dir:
		geva = load_geva(weights, args.geva_dir, args.chunk_size, selected_chromosomes, qc)
	if args.annotations:
		canonical = load_canonical(weights, args.annotations, pops, args.chunk_size, args.include_proxy, qc)
	annotation = combine_annotations(geva, canonical, args.age_priority)
	annotation = classify_age(annotation, *args.age_cutoffs, args.min_age_quality)
	annotation["effect_allele_is_dated"] = (weights["A1"] == annotation["AGE_ALT"]).where(
		annotation["AGE_GEN"].notna(), None)
	usable = int(annotation["age_usable"].sum())
	if usable == 0 and not args.allow_proxy_only:
		raise ValueError("No usable genuine age annotations matched. Provide GEVA/canonical ages with confirmed ancestral orientation and quality, or explicitly choose --allow-proxy-only.")
	mode = "human_evolution" if usable else "proxy_only"
	annotation, use_frequency = add_frequency_features(weights, annotation, pops, args.diff_cutoffs)
	annotation = permute_age(weights, annotation, args.seed)
	modules = module_specs(annotation, pops, math.log10(args.selection_p), use_frequency, args.include_proxy)
	qc += [
		("coverage", "usable_age_variants", "", usable, "Finite age, externally confirmed derived orientation, and quality threshold."),
		("coverage", "usable_age_fraction", "", usable / len(weights), "Coverage is reported; it is not prediction accuracy."),
		("definition", "evolution_mode", "", mode, "proxy_only is not an evolutionary reconstruction."),
		("definition", "min_age_quality", "", args.min_age_quality, "GEVA pair-retention quality, not pathogenicity probability."),
		("definition", "selection_threshold_p", "", args.selection_p, "Feature bin threshold, not corrected genome-wide discovery."),
	]
	for label in AGE_BINS:
		qc.append(("age_modules", label, "", int((annotation["age_bin"] == label).sum()), "Mutually exclusive age partitions include all weighted SNPs."))
	for status, count in annotation["age_status"].value_counts().items():
		qc.append(("age_quality", str(status), "", int(count), ""))
	for pop in pops:
		qc.append(("weights", "available_beta_variants", pop, int(weights[f"beta_{pop}"].notna().sum()), "Missing source-population betas score as zero and remain NA in variants table."))
	qc.append(("negative_control", "permutable_qualified_age_variants", "", int(annotation["perm_age_eligible"].sum()), "Only young/middle/old qualified age labels can change; their SNP support is fixed."))
	qc.append(("negative_control", "fixed_uncertain_or_unknown_variants", "", int((~annotation["perm_age_eligible"]).sum()), "Each uncertain/unknown SNP retains its original label and module weights."))
	qc.append(("negative_control", "permuted_labels_changed", "", int((annotation["age_bin"] != annotation["perm_age_bin"]).sum()), "Module counts preserved within strata; qualified singleton/homogeneous strata may not change. One control is not a significance test."))
	return emit_outputs(weights, annotation, modules, args.out_dir, mode, qc, args.force)


# 🚩 Explicit external-data download, integrity checks and atomic completion
def check_geva_file(path, expected_md5=None):
	md5 = hashlib.md5()
	sha256 = hashlib.sha256()
	size = 0
	with open(path, "rb") as stream:
		for block in iter(lambda: stream.read(1024 * 1024), b""):
			md5.update(block)
			sha256.update(block)
			size += len(block)
	if expected_md5 and md5.hexdigest() != expected_md5:
		raise ValueError(f"Official GEVA MD5 mismatch for {path.name}: got {md5.hexdigest()}, expected {expected_md5}")
	with gzip.open(path, "rt") as stream:
		header = None
		for line in stream:
			if not line.startswith("#"):
				header = [x.strip() for x in next(csv.reader([line]))]
				break
		if not header or not set(GEVA_AGE_COLUMNS).issubset(header):
			raise ValueError(f"Not a GEVA Atlas summary gzip: {path}")
		while stream.read(1024 * 1024):
			pass
	return {"bytes": size, "md5": md5.hexdigest(), "sha256": sha256.hexdigest(),
		"official_md5_verified": bool(expected_md5), "gzip_integrity_verified": True}


def download(args):
	directory = Path(args.geva_dir)
	directory.mkdir(parents=True, exist_ok=True)
	records = []
	for chromosome in chromosomes(args.chromosomes):
		name = f"atlas.chr{chromosome}.csv.gz"
		path = directory / name
		expected = GEVA_MD5.get(chromosome)
		url = GEVA_URL.format(chromosome=chromosome)
		if path.is_file() and not args.force:
			info = check_geva_file(path, expected)
			records.append({"chromosome": chromosome, "path": str(path), "status": "verified_existing", **info})
			print(json.dumps(records[-1]), flush=True)
			continue
		fd, temporary_name = tempfile.mkstemp(prefix=name + ".", suffix=".part", dir=directory)
		os.close(fd)
		temporary = Path(temporary_name)
		try:
			request = urllib.request.Request(url, headers={"User-Agent": "GRID-evolution/1.0"})
			with urllib.request.urlopen(request, timeout=args.timeout) as source, temporary.open("wb") as target:
				shutil.copyfileobj(source, target, length=1024 * 1024)
			info = check_geva_file(temporary, expected)
			os.replace(temporary, path)
			records.append({"chromosome": chromosome, "path": str(path), "status": "downloaded", **info})
			print(json.dumps(records[-1]), flush=True)
		finally:
			if temporary.exists():
				temporary.unlink()
	return {"n_chromosomes": len(records), "geva_dir": str(directory), "source": GEVA_DOWNLOAD_PAGE}


def parser():
	main = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	sub = main.add_subparsers(dest="command", required=True)
	build_parser = sub.add_parser("build", help="Create phenotype-free PLINK score modules.")
	build_parser.add_argument("--weights", required=True, help="TSV[.gz]: CHR BP SNP A1 A2 beta_<POP>.")
	build_parser.add_argument("--out-dir", required=True, help="Temporary exchange/output directory.")
	build_parser.add_argument("--geva-dir", help="Existing atlas.chr*.csv.gz directory; no automatic download.")
	build_parser.add_argument("--annotations", help="Canonical annotation TSV[.gz], see module doc/README.")
	build_parser.add_argument("--build", choices=["GRCh37"], default="GRCh37")
	build_parser.add_argument("--populations", default="EUR,AFR,EAS,SAS")
	build_parser.add_argument("--chromosomes", default="1-22")
	build_parser.add_argument("--chunk-size", type=int, default=100000)
	build_parser.add_argument("--age-cutoffs", type=float, nargs=2, default=[1000, 4000], metavar=("YOUNG", "OLD"))
	build_parser.add_argument("--min-age-quality", type=float, default=0.8)
	build_parser.add_argument("--age-priority", choices=["canonical", "geva"], default="canonical")
	build_parser.add_argument("--selection-p", type=float, default=0.001)
	build_parser.add_argument("--diff-cutoffs", type=float, nargs=2, default=[0.05, 0.20])
	build_parser.add_argument("--seed", type=int, default=12345)
	build_parser.add_argument("--include-proxy", action="store_true", help="Include nonnegative PROXY_* fields as explicitly separate modules.")
	build_parser.add_argument("--allow-proxy-only", action="store_true", help="Explicitly permit zero usable age matches; every result is marked proxy_only.")
	build_parser.add_argument("--force", action="store_true", help="Replace only files owned by this module builder.")
	download_parser = sub.add_parser("download", help="Explicitly download GEVA summary data to a reference directory.")
	download_parser.add_argument("--geva-dir", required=True)
	download_parser.add_argument("--chromosomes", default="1-22")
	download_parser.add_argument("--timeout", type=float, default=60)
	download_parser.add_argument("--force", action="store_true")
	return main


def main():
	args = parser().parse_args()
	try:
		result = build(args) if args.command == "build" else download(args)
		print(json.dumps(result, sort_keys=True))
	except (ValueError, OSError, EOFError, pd.errors.ParserError) as error:
		print(f"GRID evolution: {error}", file=sys.stderr)
		return 2
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
