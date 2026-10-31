#!/usr/bin/env python3
"""Build phenotype-free GRID evolutionary score modules from external annotations.

Input weights are additive effect-allele weights; age labels never change their
sign. Genome coordinates must already be normalized to GRCh37. Output is a
temporary exchange format for PLINK and the GRID fitting stage, not a table of
individual phenotypes or predictions.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile

import numpy as np
import pandas as pd


# 🚩 Public data and fixed, explicitly non-demographic feature boundaries
POPULATIONS = ("EUR", "AFR", "EAS", "SAS")
AGE_BINS = ("young", "middle", "old", "uncertain", "unknown")
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
	result = pd.to_numeric(clean.mask(clean == "", np.nan), errors="coerce")
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
	# A high-quality point estimate is not sufficient for an age-stratified
	# score: missing/incompatible bounds or bounds spanning a cut point still
	# leave the SNP in the uncertain module. Keep that weaker quality check
	# separately so it cannot make an uncertainty-only analysis look dated.
	answer["age_quality_pass"] = qualified
	answer["age_usable"] = np.isin(label, AGE_BINS[:3])
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


def permute_age(weights, annotation, seed, mode="auto"):
	answer = annotation.copy()
	if mode not in {"auto", "chr-maf", "chr-only"}:
		raise ValueError("Age permutation mode must be auto, chr-maf, or chr-only.")
	rng = np.random.default_rng(seed)
	labels = answer["age_bin"].to_numpy(dtype=object)
	permuted = labels.copy()
	# Missing/uncertain annotations have exactly the same SNP support in both
	# arms. Only qualified age-bin membership is randomized; otherwise the
	# control would also destroy annotation coverage and uncertainty patterns.
	eligible = np.isin(labels, AGE_BINS[:3])
	maf_known = answer["freq_bin"].fillna("missing").ne("missing").to_numpy()
	n_eligible, n_maf = int(eligible.sum()), int((eligible & maf_known).sum())
	if mode == "chr-maf" and n_maf < n_eligible:
		raise ValueError(
			"Strict CHR+MAF age permutation requires external reference AF for every "
			f"usable age variant; {n_eligible - n_maf} are missing. Supply reference "
			"frequencies, or choose auto/chr-only and report the exploratory limitation."
		)
	if not n_eligible:
		scheme = "not_applicable"
	elif mode == "chr-only" or n_maf == 0:
		scheme = "chr_only_exploratory"
	elif n_maf < n_eligible:
		scheme = "chr_reference_maf_partial"
	else:
		scheme = "chr_reference_maf"
	# Age-only annotations may supply no population AF. In auto mode
	# it therefore gets an explicitly named chromosome-only exploratory
	# control. Never silently claim MAF matching when it did not occur.
	freq_stratum = (pd.Series("all", index=answer.index) if mode == "chr-only"
		else answer["freq_bin"].fillna("missing"))
	strata = pd.DataFrame({"CHR": weights["CHR"], "freq_bin": freq_stratum})
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
	answer["perm_stratum"] = weights["CHR"].astype(str) + ":" + freq_stratum.astype(str)
	answer["permutation_scheme"] = scheme
	answer.attrs["permutation"] = {
		"permutation_requested_mode": mode,
		"permutation_scheme": scheme,
		"permutation_eligible_variants": n_eligible,
		"permutation_maf_covered_variants": n_maf,
		"permutation_maf_covered_fraction": n_maf / n_eligible if n_eligible else None,
		"permutation_ld_controlled": False,
	}
	return answer


# 🚩 Module definitions and phenotype-free PLINK weight matrices
def module_specs(annotation, pops, selection_log10_threshold, use_frequency, include_proxy):
	modules = []
	scheme = annotation.attrs.get("permutation", {}).get("permutation_scheme", "unspecified")
	permutation_descriptions = {
		"chr_reference_maf": "qualified ages permuted within CHR/reference MAF",
		"chr_reference_maf_partial": "qualified ages permuted within CHR/reference MAF where available; missing-AF SNPs use chromosome-only strata",
		"chr_only_exploratory": "exploratory chromosome-only qualified-age permutation; not MAF-matched",
		"not_applicable": "no usable dated SNPs to permute",
	}
	permutation_description = permutation_descriptions.get(scheme, "permutation strata not specified")
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
					else f"Fixed-seed control: {permutation_description}; uncertain/unknown SNPs stay fixed. Not LD-matched or a permutation P value.")
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
	minimum_age_variants = getattr(args, "min_age_variants", 1)
	if minimum_age_variants < 1 or int(minimum_age_variants) != minimum_age_variants:
		raise ValueError("--min-age-variants must be a positive integer.")
	pops = population_list(args.populations)
	selected_chromosomes = chromosomes(args.chromosomes)
	weights = load_weights(args.weights, pops, selected_chromosomes)
	qc = [
		("input", "variants", "", len(weights), "No target phenotype or split label is read."),
		("definition", "build", "", "GRCh37", "1-based exact coordinate plus alleles; no strand complement."),
		("definition", "age_cutoffs_generations", "", ",".join(map(str, args.age_cutoffs)), "Feature boundaries, not population separation dates."),
		("definition", "permutation_seed", "", args.seed, "One fixed negative control, not a permutation P value; qualified age labels only. Actual matching scheme is reported below."),
		("definition", "age_interpretation", "", "inferred_derived_allele_origin_generations", "Estimated from genetic data; not directly observed historical dates or AFR/EUR/EAS/SAS divergence dates."),
	]
	annotation = load_canonical(weights, args.annotations, pops, args.chunk_size, args.include_proxy, qc)
	annotation = classify_age(annotation, *args.age_cutoffs, args.min_age_quality)
	annotation["effect_allele_is_dated"] = (weights["A1"] == annotation["AGE_ALT"]).where(
		annotation["AGE_GEN"].notna(), None)
	usable = int(annotation["age_usable"].sum())
	weighted = weights[[f"beta_{pop}" for pop in pops]].fillna(0).ne(0).any(axis=1)
	usable_weighted = int((annotation["age_usable"] & weighted).sum())
	if usable_weighted < minimum_age_variants and not args.allow_proxy_only:
		raise ValueError(
			f"Only {usable_weighted} usable, nonzero-weight age variants matched; "
			f"--min-age-variants requires {minimum_age_variants}. Provide canonical "
			"ages with confirmed derived orientation, adequate quality, and compatible "
			"bounds wholly inside one age bin, or explicitly choose --allow-proxy-only."
		)
	mode = "human_evolution" if usable_weighted >= minimum_age_variants else "proxy_only"
	annotation, use_frequency = add_frequency_features(weights, annotation, pops, args.diff_cutoffs)
	annotation = permute_age(weights, annotation, args.seed, getattr(args, "age_permutation_mode", "auto"))
	permutation_info = annotation.attrs["permutation"]
	changed = annotation["age_bin"] != annotation["perm_age_bin"]
	bins_present = [label for label in AGE_BINS[:3] if ((annotation["age_bin"] == label) & weighted).any()]
	quality_pass = int(annotation["age_quality_pass"].sum())
	age_info = {
		"age_quality_pass_variants": quality_pass,
		"age_usable_variants": usable,
		"age_usable_fraction": usable / len(weights),
		"age_usable_weighted_variants": usable_weighted,
		"age_bins_present": bins_present,
		"age_bin_count": len(bins_present),
		"minimum_age_variants": int(minimum_age_variants),
		"age_partition_status": "multiple_dated_bins" if len(bins_present) >= 2 else "single_dated_bin" if bins_present else "no_dated_bins",
		"permutation_changed_variants": int(changed.sum()),
		"permutation_changed_weighted_variants": int((changed & weighted).sum()),
		"permutation_has_weighted_contrast": bool((changed & weighted).any()),
		**permutation_info,
	}
	modules = module_specs(annotation, pops, math.log10(args.selection_p), use_frequency, args.include_proxy)
	qc += [
		("coverage", "age_quality_pass_variants", "", quality_pass, "Finite age, externally confirmed derived orientation, and quality threshold; alone not a usable age-bin assignment."),
		("coverage", "usable_age_variants", "", usable, "Adequate quality plus finite compatible bounds entirely inside young/middle/old; uncertain/unknown do not count."),
		("coverage", "usable_age_fraction", "", usable / len(weights), "Coverage is reported; it is not prediction accuracy."),
		("coverage", "usable_weighted_age_variants", "", usable_weighted, "Usable age assignment and nonzero beta in at least one requested source population."),
		("coverage", "weighted_age_bins_present", "", ",".join(bins_present), "One bin alone does not establish a young-versus-old contrast."),
		("definition", "minimum_age_variants", "", minimum_age_variants, "A data-contract floor, not a scientifically sufficient sample-size or coverage requirement."),
		("definition", "evolution_mode", "", mode, "proxy_only is not an evolutionary reconstruction."),
		("definition", "min_age_quality", "", args.min_age_quality, "Declared source quality, not pathogenicity probability."),
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
	qc.append(("negative_control", "permuted_labels_changed", "", age_info["permutation_changed_variants"], "Module counts preserved within strata; qualified singleton/homogeneous strata may not change. One control is not a significance test."))
	qc.append(("negative_control", "permuted_weighted_labels_changed", "", age_info["permutation_changed_weighted_variants"], "A value of zero means no weighted age negative-control contrast was created."))
	qc.append(("negative_control", "permutation_scheme", "", permutation_info["permutation_scheme"], "chr_only_exploratory does not control for MAF; partial reference MAF leaves missing-AF SNPs chromosome-stratified only."))
	qc.append(("negative_control", "permutation_maf_covered_fraction", "", permutation_info["permutation_maf_covered_fraction"], "Among usable age variants; external source frequencies only, never inferred from missing data."))
	qc.append(("negative_control", "permutation_ld_controlled", "", False, "CHR+MAF stratification does not match LD, recombination rate, or other genomic annotations."))
	return emit_outputs(weights, annotation, modules, args.out_dir, mode, qc, args.force) | age_info


def parser():
	main = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	sub = main.add_subparsers(dest="command", required=True)
	build_parser = sub.add_parser("build", help="Create phenotype-free PLINK score modules.")
	build_parser.add_argument("--weights", required=True, help="TSV[.gz]: CHR BP SNP A1 A2 beta_<POP>.")
	build_parser.add_argument("--out-dir", required=True, help="Temporary exchange/output directory.")
	build_parser.add_argument("--annotations", required=True, help="Canonical annotation TSV[.gz], see module doc/README.")
	build_parser.add_argument("--build", choices=["GRCh37"], default="GRCh37")
	build_parser.add_argument("--populations", default="EUR,AFR,EAS,SAS")
	build_parser.add_argument("--chromosomes", default="1-22")
	build_parser.add_argument("--chunk-size", type=int, default=100000)
	build_parser.add_argument("--age-cutoffs", type=float, nargs=2, default=[1000, 4000], metavar=("YOUNG", "OLD"))
	build_parser.add_argument("--min-age-quality", type=float, default=0.8)
	build_parser.add_argument("--min-age-variants", type=int, default=1,
		help="Minimum usable nonzero-weight age variants; a data-contract floor, not a scientific sufficiency threshold.")
	build_parser.add_argument("--age-permutation-mode", choices=["auto", "chr-maf", "chr-only"], default="auto",
		help="auto names any missing-MAF fallback explicitly; chr-maf requires reference AF for every usable age SNP; chr-only is exploratory.")
	build_parser.add_argument("--selection-p", type=float, default=0.001)
	build_parser.add_argument("--diff-cutoffs", type=float, nargs=2, default=[0.05, 0.20])
	build_parser.add_argument("--seed", type=int, default=12345)
	build_parser.add_argument("--include-proxy", action="store_true", help="Include nonnegative PROXY_* fields as explicitly separate modules.")
	build_parser.add_argument("--allow-proxy-only", action="store_true", help="Explicitly permit fewer usable weighted ages than required; the result is marked proxy_only.")
	build_parser.add_argument("--force", action="store_true", help="Replace only files owned by this module builder.")
	return main


def main():
	args = parser().parse_args()
	try:
		result = build(args)
		print(json.dumps(result, sort_keys=True))
	except (ValueError, OSError, EOFError, pd.errors.ParserError) as error:
		print(f"GRID evolution: {error}", file=sys.stderr)
		return 2
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
