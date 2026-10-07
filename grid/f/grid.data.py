"""GRID data alignment, outcome-blind shared splits and training-centred SNP modules.

Public API: prepare_data(args), evolution_data(args, prepared).
The caller owns model fitting, cache lifecycle and formal output publication.
"""
from __future__ import annotations

import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


# 🚩 Paths, tables and identities
POPS = ("EUR", "AFR", "EAS", "SAS")
TRAITS = ("height", "ldl", "t2dm")
CSX = [f"csx.{pop}" for pop in POPS]


def option(args, name, default=None):
	value = getattr(args, name, default)
	return default if value is None else value


def enabled(value):
	return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def optional(value):
	return None if value is None or str(value).strip().lower() in {"", "none", "null"} else Path(value).expanduser()


def required(value):
	path = optional(value)
	if path is None or not path.is_file() or not path.stat().st_size:
		raise ValueError(f"Missing or empty input: {value}")
	return path.resolve()


def names(value):
	if value is None or str(value).lower() in {"", "none"}:
		return []
	return [str(x).strip() for x in value] if isinstance(value, (list, tuple)) else [x.strip() for x in str(value).split(",") if x.strip()]


def selected_traits(args):
	result = names(option(args, "traits", option(args, "trait", "height,ldl,t2dm")))
	if not result or len(set(result)) != len(result) or set(result) - set(TRAITS):
		raise ValueError("Traits must be unique members of height,ldl,t2dm")
	return result


def selected_chromosomes(args):
	result = []
	for token in re.split(r"[,;\s]+", str(option(args, "chrs", "1-22"))):
		if not token:
			continue
		token = re.sub(r"^chr", "", token, flags=re.I)
		if "-" in token:
			a, b = map(int, token.split("-"))
			if b < a:
				raise ValueError("Descending chromosome range")
			result.extend(range(a, b + 1))
		else:
			result.append(int(token))
	result = sorted(set(result))
	if not result or min(result) < 1 or max(result) > 22:
		raise ValueError("GRID supports autosomes 1-22")
	return result


def expand_trait(value, trait):
	return str(value).replace("{trait}", trait)


def log(message):
	print(f"[GRID data] {message}", flush=True)


@contextlib.contextmanager
def text_stream(path):
	path = Path(path)
	if path.suffix == ".gz":
		with gzip.open(path, "rt") as handle:
			yield handle
	elif path.suffix == ".zst":
		try:
			import zstandard
		except ImportError:
			command = shutil.which("zstd")
			if command is None:
				raise ValueError("Reading .zst requires the zstandard Python package or zstd command")
			with subprocess.Popen([command, "-dc", str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
				handle = io.TextIOWrapper(process.stdout)
				partial = True
				try:
					yield handle
					partial = bool(handle.read(1))
				finally:
					# Header-only readers intentionally stop the decompressor early.
					if partial and process.poll() is None:
						process.terminate()
					handle.close()
					stderr = process.stderr.read().decode(errors="replace")
					returncode = process.wait()
				if not partial and returncode != 0:
					raise ValueError(f"zstd failed for {path}: {stderr}")
		else:
			with path.open("rb") as source, zstandard.ZstdDecompressor().stream_reader(source) as reader:
				with io.TextIOWrapper(reader) as handle:
					yield handle
	else:
		with path.open() as handle:
			yield handle


def read_table(value):
	path = required(value)
	if path.suffix.lower() == ".rds":
		import rdata
		data = rdata.read_rds(path)
		if isinstance(data, dict) and len(data) == 1 and isinstance(next(iter(data.values())), pd.DataFrame):
			data = next(iter(data.values()))
		if not isinstance(data, pd.DataFrame):
			raise ValueError(f"Expected one data.frame/data.table in {path}")
		return data.reset_index(drop=True)
	with text_stream(path) as stream:
		lines = 0
		for line in stream:
			if line.startswith("##") or not line.strip():
				lines += 1
				continue
			sep = "\t" if "\t" in line else "," if "," in line else r"\s+"
			break
		else:
			raise ValueError(f"No header in {path}")
	with text_stream(path) as stream:
		return pd.read_csv(stream, sep=sep, skiprows=lines, dtype=str, keep_default_na=True)


def write_rds(data, path):
	import rdata
	path = Path(path)
	path.parent.mkdir(parents=True, exist_ok=True)
	frame = data.copy()
	for column in frame:
		if str(frame[column].dtype) in {"string", "category"}:
			frame[column] = frame[column].astype(object)
		elif frame[column].dtype == np.dtype("float32"):
			frame[column] = frame[column].astype(float)
	with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".rds", delete=False) as stream:
		temporary = Path(stream.name)
	try:
		rdata.write_rds(temporary, frame, compression="gzip")
		os.replace(temporary, path)
	finally:
		temporary.unlink(missing_ok=True)


def normalize_ids(values):
	series = pd.Series(values, copy=False).astype("string").str.strip()
	series = series.str.replace(r"^([+-]?\d+)\.0+$", r"\1", regex=True)
	if series.isna().any() or series.isin(["", "NA", "NaN", "nan", "None"]).any() or series.str.contains(r"\s", regex=True).any():
		raise ValueError("Missing, blank or whitespace-containing individual ID")
	return series.astype(str).to_numpy()


def identified(data, label):
	data = data.copy()
	column = next((x for x in ("eid", "IID", "#IID", "ID_2", "id") if x in data), None)
	if column is None:
		raise ValueError(f"{label}: missing eid/IID column")
	data["eid"] = normalize_ids(data[column])
	if data.eid.duplicated().any():
		raise ValueError(f"{label}: duplicate IDs")
	return data


def numeric(series, label):
	clean = series.replace({"": np.nan, "NA": np.nan, "NaN": np.nan, ".": np.nan})
	values = pd.to_numeric(clean, errors="coerce")
	if (clean.notna() & values.isna()).any() or np.isinf(values.to_numpy(dtype=float, na_value=np.nan)).any():
		raise ValueError(f"{label}: invalid numeric value; encode categorical covariates explicitly")
	return values.astype(float)


def id_list(value):
	path = optional(value)
	if path is None:
		return None
	if not path.is_file():
		raise ValueError(f"Missing ID list: {path}")
	values = []
	with text_stream(path) as stream:
		for line in stream:
			fields = line.split()
			if not fields or fields[0].startswith("#") or fields[0].upper() in {"EID", "IID", "FID", "ID"}:
				continue
			values.append(fields[1] if len(fields) > 1 else fields[0])
	return set(normalize_ids(values)) if values else set()


def filter_ids(data, args):
	keep = id_list(option(args, "keep"))
	remove = id_list(option(args, "remove", "/mnt/d/files/ukb.exclude.id"))
	mask = ~data.eid.str.startswith("-")
	if keep is not None:
		mask &= data.eid.isin(keep)
	if remove is not None:
		mask &= ~data.eid.isin(remove)
	return data.loc[mask].copy()


def source_stamp(path, allow_empty=False):
	path = Path(path).expanduser().resolve() if allow_empty else required(path)
	if not path.is_file():
		raise ValueError(f"Missing input: {path}")
	stat = path.stat()
	return {"path": str(path), "size": int(stat.st_size), "mtime_ns": int(stat.st_mtime_ns)}


# 🚩 Genotype metadata; no allele frequencies are estimated here
def genotype_files(args):
	root = Path(option(args, "dir_gen", "/mnt/f/gen/ukb/37/hap"))
	result = {}
	for chrom in selected_chromosomes(args):
		prefix = root / f"chr{chrom}"
		if Path(str(prefix) + ".pgen").is_file():
			variant = next((Path(str(prefix) + suffix) for suffix in (".pvar", ".pvar.zst") if Path(str(prefix) + suffix).is_file()), None)
			if variant is None:
				raise ValueError(f"Missing PVAR: {prefix}")
			sample = required(str(prefix) + ".psam")
			command = ["--pfile", str(prefix)] + (["vzs"] if variant.suffix == ".zst" else [])
			binary = required(str(prefix) + ".pgen")
			mode = "pfile"
		else:
			variant = required(str(prefix) + ".bim")
			sample = required(str(prefix) + ".fam")
			binary = required(str(prefix) + ".bed")
			command = ["--bfile", str(prefix)]
			mode = "bfile"
		result[chrom] = {"prefix": str(prefix), "variant": variant, "sample": sample,
			"binary": binary, "mode": mode, "command": command}
	return result


def genotype_ids(info):
	if info["mode"] == "bfile":
		data = pd.read_csv(info["sample"], sep=r"\s+", header=None, usecols=[1], dtype=str)
		data.columns = ["eid"]
	else:
		data = read_table(info["sample"])
	return identified(data, str(info["sample"])).eid.tolist()


# 🚩 A shared, outcome-blind family split
def family_assignment(families, fraction, seed):
	counts = pd.Series(families, dtype=str).value_counts()
	if len(counts) < 2:
		raise ValueError("At least two independent family groups are required for a split")
	order = sorted(counts.index, key=lambda value: (-int(counts[value]), hashlib.sha256(f"{seed}:{value}".encode()).digest()))
	membership = {}
	first = second = 0
	for family in order:
		choose_first = first / fraction <= second / (1 - fraction)
		membership[family] = bool(choose_first)
		if choose_first:
			first += int(counts[family])
		else:
			second += int(counts[family])
	if min(first, second) == 0:
		raise ValueError("Family groups cannot form two nonempty splits")
	return membership


def normalize_split(series):
	return series.astype(str).str.strip().str.lower().replace({"training": "train", "testing": "test"})


def attach_families(cohort, args):
	path = optional(option(args, "split_group_file"))
	if path is not None:
		groups = identified(read_table(path), "Family roster")
		column = next((x for x in ("family_id", "group") if x in groups), None)
		if column is None:
			raise ValueError("Family roster requires eid,family_id or eid,group")
		groups["family_id"] = normalize_ids(groups[column])
		cohort = cohort.drop(columns=["family_id"], errors="ignore").merge(groups[["eid", "family_id"]], on="eid", how="left", validate="one_to_one")
		if cohort.family_id.isna().any():
			raise ValueError("Family roster does not cover every retained individual")
		status = "supplied_family_groups"
	elif "family_id" in cohort:
		cohort["family_id"] = normalize_ids(cohort.family_id)
		status = "prepared_family_groups"
	else:
		cohort["family_id"] = cohort.eid
		status = "individual_ID_groups; no claim of unrelatedness"
	return cohort, status


def freeze_split(cohort, args, supplied=None):
	cohort = cohort.sort_values("eid").reset_index(drop=True)
	cohort, family_status = attach_families(cohort, args)
	if optional(option(args, "split_file")) is not None:
		external = identified(read_table(option(args, "split_file")), "Outer split")
		if "split" not in external:
			raise ValueError("Outer split file requires eid,split")
		external["split"] = normalize_split(external.split)
		if supplied is not None:
			comparison = supplied[["eid", "split"]].merge(external[["eid", "split"]], on="eid", suffixes=("_data", "_file"), how="inner")
			if not (comparison.split_data == comparison.split_file).all():
				raise ValueError("Prepared table split conflicts with --split-file")
		supplied = external
	if supplied is not None:
		cohort = cohort.drop(columns=["split"], errors="ignore").merge(supplied[["eid", "split"]], on="eid", how="left", validate="one_to_one")
		cohort["split"] = normalize_split(cohort.split)
		if not cohort.split.isin(["train", "test"]).all():
			raise ValueError("Every retained person needs explicit split=train or test")
	else:
		membership = family_assignment(cohort.family_id, .5, int(option(args, "seed", 20260904)))
		cohort["split"] = np.where(cohort.family_id.map(membership), "train", "test")
	if set(cohort.split) != {"train", "test"} or cohort.groupby("family_id").split.nunique().max() != 1:
		raise ValueError("Outer split must contain train/test with disjoint family groups")
	fraction = float((cohort.split == "train").mean())
	if abs(fraction - .5) > float(option(args, "half_split_tolerance", .02)):
		raise ValueError(f"Expected a 50/50 shared roster; training fraction is {fraction:.4f}. Large families or supplied splits need review.")
	cache = Path(option(args, "cache_dir", "/tmp/grid-cache/grid"))
	cache.mkdir(parents=True, exist_ok=True)
	split_path = cache / "split.rds"
	if split_path.is_file() and not enabled(option(args, "replace", False)):
		previous = identified(read_table(split_path), "Saved outer split")
		if set(previous.eid) != set(cohort.eid):
			raise ValueError("Existing shared split has a different cohort; use a new cache or explicit replacement")
		joined = cohort.merge(previous[["eid", "split", "family_id", "ancestry"]], on="eid", suffixes=("", "_old"), validate="one_to_one")
		if not ((joined.split == joined.split_old) & (joined.family_id == joined.family_id_old) & (joined.ancestry == joined.ancestry_old)).all():
			raise ValueError("Existing shared split/families/ancestry differ; use a new cache or explicit replacement")
	write_rds(cohort[["eid", "ancestry", "family_id", "split"]], split_path)
	cohort[["eid"]].rename(columns={"eid": "#IID"}).to_csv(cache / "common.keep", sep="\t", index=False)
	cohort.loc[cohort.split == "train", ["eid"]].rename(columns={"eid": "#IID"}).to_csv(cache / "train.keep", sep="\t", index=False)
	cohort[["eid", "family_id"]].rename(columns={"family_id": "group"}).to_csv(cache / "families.tsv.gz", sep="\t", index=False)
	internal = cohort[["eid", "family_id", "split"]].copy()
	train = internal.split == "train"
	membership = family_assignment(internal.loc[train, "family_id"], .8, int(option(args, "seed", 20260904)) + 1)
	internal.loc[train, "split"] = np.where(internal.loc[train, "family_id"].map(membership), "train", "validation")
	if internal.groupby("family_id").split.nunique().max() != 1:
		raise AssertionError("PRSformer family groups crossed fit/validation/test")
	internal[["eid", "split"]].to_csv(cache / "prsformer.split.tsv.gz", sep="\t", index=False)
	log(f"Frozen common roster: N={len(cohort)}, train={int(train.sum())}, test={int((~train).sum())}; {family_status}")
	return cohort, {"family_status": family_status, "train_fraction": fraction,
		"split_uses_outcome_values": False, "prsformer_split_counts": {str(k): int(v) for k, v in internal.split.value_counts().items()}}


# 🚩 Phenotypes, original source scores and projected matching coordinates
def outcome(data, trait, args):
	column = option(args, f"{trait}_col", "auto" if trait == "t2dm" else trait)
	if trait == "t2dm" and column == "auto":
		if not {"t2dm.Yr2e", "t2dm.Yt2e"}.issubset(data.columns):
			raise ValueError("T2DM auto needs t2dm.Yr2e and t2dm.Yt2e; supply --t2dm-col for an explicit baseline 0/1 endpoint")
		pre = numeric(data["t2dm.Yr2e"], "t2dm.Yr2e")
		inc = numeric(data["t2dm.Yt2e"], "t2dm.Yt2e")
		values = pd.Series(np.where(pre == 1, 1., np.where((pre == 0) | inc.isin([0., 1.]), 0., np.nan)), index=data.index)
	else:
		if column not in data:
			raise ValueError(f"Missing phenotype column: {column}")
		values = numeric(data[column], column)
	if trait == "t2dm" and not values.dropna().isin([0., 1.]).all():
		raise ValueError("T2DM must be a defined baseline 0/1 outcome, not an incident event indicator or follow-up time")
	return values


def feature_groups(data, covariates, evolution_mode=None):
	shared_evolution = [column for column in data if column.startswith(("evo.sel", "evo.proxy"))]
	groups = {
		"covariates": list(covariates),
		"csx": [column for column in CSX if column in data],
		"ancestry": [column for column in data if re.fullmatch(r"match_PC\d+", column)],
		"evolution": [column for column in data if column.startswith("evo.age")] + shared_evolution,
		"frequency": [column for column in data if column.startswith("evo.diff")],
		"permuted_evolution": [column for column in data if column.startswith("evo.perm_age")] + shared_evolution
			if any(column.startswith("evo.perm_age") for column in data) else [],
		"disco": ["disco"] if "disco" in data and data.disco.notna().all() else [],
	}
	if evolution_mode == "proxy_only":
		# The builder still emits all-unknown age partitions for coverage accounting.
		# These cannot establish an evolution or an age-permutation comparison.
		groups["evolution"] = []
		groups["permuted_evolution"] = []
		groups["frequency"] += [column for column in data if column.startswith("evo.proxy")]
	return groups


def validate_predictors(data, covariates, label):
	missing = set(CSX + covariates) - set(data.columns)
	if missing:
		raise ValueError(f"{label}: missing predictors {sorted(missing)}")
	columns = [x for x in data if x in CSX + covariates + ["disco"] or x.startswith(("match_PC", "evo."))]
	for column in columns:
		data[column] = numeric(data[column], f"{label}.{column}")
	return data


def prepared_tables(args, traits, covariates):
	tables, sources = {}, []
	pattern = str(option(args, "data_file"))
	for trait in traits:
		path = required(expand_trait(pattern, trait))
		data = read_table(path)
		if "trait" in data:
			data = data.loc[data.trait.astype(str).str.lower() == trait].copy()
		elif len(traits) > 1 and "{trait}" not in pattern:
			if trait in data:
				data["y"] = data[trait]
			else:
				raise ValueError("Multi-trait --data-file needs {trait}, a trait column, or separate named phenotype columns")
		data = identified(data, f"Prepared {trait}")
		if not {"y", "split"}.issubset(data.columns):
			raise ValueError(f"Prepared {trait}: requires eid,split,y")
		data = filter_ids(data, args)
		data["y"] = numeric(data.y, f"{trait}.y")
		if trait == "t2dm" and not data.y.dropna().isin([0., 1.]).all():
			raise ValueError("Prepared T2DM y must be 0/1")
		data["split"] = normalize_split(data.split)
		if "ancestry" not in data:
			column = next((x for x in ("target", "genetic_ancestry") if x in data), None)
			data["ancestry"] = data[column] if column else "UNASSIGNED"
		data["ancestry"] = data.ancestry.fillna("UNASSIGNED").astype(str)
		data = validate_predictors(data, covariates, trait)
		core = CSX + covariates + [x for x in data if re.fullmatch(r"match_PC\d+", x)]
		data = data.loc[data[core].notna().all(axis=1)].copy()
		tables[trait] = data
		sources.append(source_stamp(path))
	return tables, sources


def raw_tables(args, traits, covariates):
	path = required(option(args, "pheno_file", "/mnt/d/data/ukb/phe/Rdata/all.rds"))
	data = filter_ids(identified(read_table(path), "Phenotype"), args)
	source_traits = {option(args, f"{t}_col", "auto" if t == "t2dm" else t) for t in traits}
	if set(covariates) & (set(traits) | source_traits | {"y", "outcome", "t2dm.Yr2e", "t2dm.Yt2e"}):
		raise ValueError("An outcome cannot be used as a covariate")
	for trait in traits:
		data[trait] = outcome(data, trait, args)
	for column in covariates:
		if column not in data:
			continue
		data[column] = numeric(data[column], column)
	missing = set(covariates) - set(data)
	if missing:
		raise ValueError(f"Missing phenotype covariates: {sorted(missing)}")
	data = data[["eid", *covariates, *traits]].copy()
	sources = [source_stamp(path)]
	ancestry_path = required(option(args, "ancestry_file", "/mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz"))
	ancestry = identified(read_table(ancestry_path), "Ancestry")
	column = next((x for x in (option(args, "group_col", "genetic_ancestry"), "ancestry", "target", "predicted_ancestry") if x in ancestry), None)
	if column is None:
		raise ValueError("No ancestry label column")
	ancestry["ancestry"] = ancestry[column].fillna("UNASSIGNED").astype(str)
	data = data.merge(ancestry[["eid", "ancestry"]], on="eid", how="inner", validate="one_to_one")
	pca_path = required(option(args, "pca_file", "/mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz"))
	pca = identified(read_table(pca_path), "Projected PCA")
	npc = int(option(args, "distance_pcs", 10))
	if not 1 <= npc <= 20:
		raise ValueError("distance_pcs must be between 1 and 20")
	pcs = [f"PC{j}" for j in range(1, npc + 1)]
	if set(pcs) - set(pca):
		raise ValueError("Projected PCA lacks requested matching PCs")
	for pc in pcs:
		pca[pc] = numeric(pca[pc], f"projected.{pc}")
	pca = pca[["eid", *pcs]].rename(columns={pc: "match_" + pc for pc in pcs})
	data = data.merge(pca, on="eid", how="inner", validate="one_to_one")
	data = data.loc[data[covariates + ["match_" + pc for pc in pcs]].notna().all(axis=1)].copy()
	sources += [source_stamp(ancestry_path), source_stamp(pca_path)]
	genotypes = genotype_files(args)
	common = set(data.eid)
	for info in genotypes.values():
		common.intersection_update(genotype_ids(info))
		sources += [source_stamp(info["sample"]), source_stamp(info["variant"]), source_stamp(info["binary"])]
	data = data.loc[data.eid.isin(common)].copy()
	tables = {}
	for trait in traits:
		score_path = required(Path(option(args, "score_dir", "/mnt/d/data/ukb/pgs")) / trait / "1csx.scores.rds")
		score = identified(read_table(score_path), f"{trait} CSx")
		if set(CSX) - set(score):
			raise ValueError(f"{score_path}: requires all four csx.<POP> scores")
		columns = CSX + [x for x in ("csx.auto", "csx.meta") if x in score]
		for name in columns:
			score[name] = numeric(score[name], f"{trait}.{name}")
		frame = data[["eid", "ancestry", *covariates, *["match_" + pc for pc in pcs], trait]].rename(columns={trait: "y"})
		frame = frame.merge(score[["eid", *columns]], on="eid", how="inner", validate="one_to_one")
		frame = frame.loc[frame[CSX].notna().all(axis=1)].copy()
		disco_path = score_path.with_name("2disco.scores.rds")
		if disco_path.is_file():
			disco = identified(read_table(disco_path), f"{trait} genotype-only Disco")
			if "disco" not in disco:
				raise ValueError(f"{disco_path}: missing disco score")
			disco["disco"] = numeric(disco.disco, "disco")
			frame = frame.merge(disco[["eid", "disco"]], on="eid", how="left", validate="one_to_one")
			sources.append(source_stamp(disco_path))
		frame = validate_predictors(frame, covariates, trait)
		tables[trait] = frame
		sources.append(source_stamp(score_path))
	return tables, sources


def prepare_data(args):
	traits = selected_traits(args)
	covariates = names(option(args, "covariates", "age,sex,PC1,PC2"))
	if len(set(covariates)) != len(covariates) or set(covariates) & {"eid", "split", "family_id", "y", "outcome", *TRAITS}:
		raise ValueError("Covariates must be unique numeric predictors, not outcomes/identifiers")
	external = optional(option(args, "data_file")) is not None
	tables, sources = prepared_tables(args, traits, covariates) if external else raw_tables(args, traits, covariates)
	for name in ("keep", "remove", "split_file", "split_group_file"):
		path = optional(option(args, name))
		if path is not None:
			sources.append(source_stamp(path, allow_empty=name in {"keep", "remove"}))
	common = set.intersection(*(set(frame.eid) for frame in tables.values()))
	any_label = set.union(*(set(frame.loc[frame.y.notna(), "eid"]) for frame in tables.values()))
	common.intersection_update(any_label)
	if len(common) < 10:
		raise ValueError(f"Only {len(common)} shared individuals with at least one observed trait")
	first = tables[traits[0]].set_index("eid").loc[sorted(common)].reset_index()
	cohort = first[["eid", "ancestry"] + (["family_id"] if "family_id" in first else [])].copy()
	supplied = first[["eid", "split"]].copy() if external else None
	for trait, frame in tables.items():
		frame = frame.set_index("eid").loc[cohort.eid].reset_index()
		if not np.array_equal(frame.ancestry.astype(str), cohort.ancestry.astype(str)):
			raise ValueError(f"Ancestry differs between prepared traits: {trait}")
		if external and not np.array_equal(frame.split, supplied.split):
			raise ValueError(f"Outer split differs between prepared traits: {trait}")
		if ("family_id" in cohort) != ("family_id" in frame):
			raise ValueError(f"Family groups must be present consistently in prepared traits: {trait}")
		if "family_id" in cohort and not np.array_equal(frame.family_id.astype(str), cohort.family_id.astype(str)):
			raise ValueError(f"Family groups differ between traits: {trait}")
		tables[trait] = frame
	cohort, split_metadata = freeze_split(cohort, args, supplied)
	groups = {}
	for trait, frame in tables.items():
		frame = frame.drop(columns=["split", "family_id", "target"], errors="ignore").merge(cohort[["eid", "split", "family_id"]], on="eid", validate="one_to_one")
		frame["target"] = frame.ancestry
		tables[trait] = frame
		groups[trait] = feature_groups(frame, covariates)
	t2dm_definition = "not_requested"
	if "t2dm" in traits:
		column = option(args, "t2dm_col", "auto")
		t2dm_definition = "externally_prepared_baseline_0_1; endpoint definition supplied by provider" if external else (
			"baseline_0_1; Yr2e=1 case, Yr2e=0 or valid Yt2e control" if column == "auto" else f"baseline_0_1; explicit phenotype column {column}")
	metadata = {"input_mode": "external_prepared_features" if external else "UKB_CSx_and_reference_projected_PC",
		"sources": sources, "covariates": covariates, "traits": traits, **split_metadata,
		"n_shared": len(cohort), "labels_kept_missing_until_evaluation": True,
		"t2dm_endpoint": t2dm_definition,
		"disco_comparator": {trait: {"score_column_present": "disco" in tables[trait],
			"n_observed_scores": int(tables[trait].disco.notna().sum()) if "disco" in tables[trait] else 0,
			"used_on_complete_common_cohort": bool(groups[trait]["disco"])} for trait in traits},
		"match_PC_source": "external_prepared_features" if external else "reference_projected_PCs_separate_from_covariates"}
	metadata["evolution_verified"] = False
	return {"tables": tables, "feature_groups": groups, "cohort": cohort, "metadata": metadata}


# 🚩 Strict canonical weights; no old transportability re-shrinkage
def load_module(path, name):
	path = required(path)
	spec = importlib.util.spec_from_file_location(name, path)
	module = importlib.util.module_from_spec(spec)
	sys.modules[name] = module
	spec.loader.exec_module(module)
	return module


def normalize_weight_table(frame, label):
	needed = {"CHR", "BP", "SNP", "A1", "A2"}
	if needed - set(frame):
		raise ValueError(f"{label}: missing weight columns {sorted(needed - set(frame))}")
	frame = frame.copy()
	frame["CHR"] = numeric(frame.CHR.astype(str).str.replace(r"^chr", "", case=False, regex=True), f"{label}.CHR")
	frame["BP"] = numeric(frame.BP, f"{label}.BP")
	if frame[["CHR", "BP"]].isna().any().any() or (frame.CHR % 1 != 0).any() or (frame.BP % 1 != 0).any() or not frame.CHR.between(1, 22).all() or (frame.BP < 1).any():
		raise ValueError(f"{label}: invalid autosomal 1-based coordinates")
	frame["CHR"], frame["BP"] = frame.CHR.astype(int), frame.BP.astype(int)
	frame["SNP"] = frame.SNP.astype(str).str.strip()
	for column in ("A1", "A2"):
		frame[column] = frame[column].astype(str).str.upper().str.strip()
		if not frame[column].str.fullmatch("[ACGT]").all():
			raise ValueError(f"{label}: only normalized biallelic SNPs are supported")
	if (frame.A1 == frame.A2).any() or frame.SNP.isin(["", "NA", "nan", "."]).any() or frame.SNP.duplicated().any():
		raise ValueError(f"{label}: identical alleles or missing/duplicate SNP IDs")
	if "BUILD" in frame and not frame.BUILD.astype(str).str.lower().isin(["37", "grch37", "hg19"]).all():
		raise ValueError(f"{label}: BUILD conflicts with GRCh37")
	return frame


def native_weights(args, trait):
	root = Path(option(args, "gwas_dir", "/mnt/f/gwas/4grid/common"))
	chrs = selected_chromosomes(args)
	common = load_module(Path(__file__).with_name("0.common.py"), "grid_data_common")
	snpinfo = required(option(args, "snpinfo", "/mnt/f/refLD/csx/snpinfo_mult_1kg_hm3"))
	canonical, signatures, sources = None, [], []
	for pop in POPS:
		tags = [f"{trait}.{pop}"] + (["t2dm.AFA"] if trait == "t2dm" and pop == "AFR" else [])
		paths = [candidate for tag in tags for candidate in (root / tag / "gwas" / f"{tag}.gz", root / f"{tag}.gz")]
		source = next((path for path in paths if path.is_file()), None)
		if source is None:
			raise ValueError(f"No source GWAS for {trait}.{pop} under {root}")
		suffix = "" if chrs == list(range(1, 23)) else ".chr" + ",".join(map(str, chrs))
		weight_options = [Path(str(source)[:-3] + suffix + ".csx.gz"), Path(str(source)[:-3] + ".csx.gz")]
		weight = next((path for path in weight_options if path.is_file()), None)
		if weight is None:
			raise ValueError(f"Missing completed CSx weights for {trait}.{pop}; run 1.csx.sh first")
		metadata_path = required(str(weight) + ".metadata.json")
		signature_path = required(str(weight) + ".signature")
		metadata = json.loads(metadata_path.read_text())
		if metadata.get("preparation_signature") != common.preparation_signature(source, snpinfo, trait, pop):
			raise ValueError(f"CSx weights/source preparation mismatch: {weight}")
		signatures.append(signature_path.read_text().strip())
		sources.extend(source_stamp(path) for path in (source, weight, metadata_path, signature_path))
		frame = normalize_weight_table(read_table(weight), str(weight))
		if "BETA" not in frame:
			raise ValueError(f"Missing CSx BETA: {weight}")
		frame[f"beta_{pop}"] = numeric(frame.BETA, str(weight) + ".BETA")
		frame = frame.loc[frame.CHR.isin(chrs), ["CHR", "BP", "SNP", "A1", "A2", f"beta_{pop}"]]
		if frame.empty or frame[f"beta_{pop}"].isna().any():
			raise ValueError(f"No finite selected CSx effects: {weight}")
		if canonical is None:
			canonical = frame
			continue
		merged = canonical.merge(frame, on="SNP", how="outer", suffixes=("", "_source"), validate="one_to_one")
		both = merged.CHR.notna() & merged.CHR_source.notna()
		same = (merged.A1 == merged.A1_source) & (merged.A2 == merged.A2_source)
		swap = (merged.A1 == merged.A2_source) & (merged.A2 == merged.A1_source)
		bad = both & ((merged.CHR != merged.CHR_source) | (merged.BP != merged.BP_source) | ~(same | swap))
		if bad.any():
			raise ValueError(f"CSx populations disagree on exact coordinate/alleles: {merged.loc[bad, 'SNP'].head().tolist()}; normalize sources explicitly")
		merged.loc[both & swap, f"beta_{pop}"] *= -1
		for column in ("CHR", "BP", "A1", "A2"):
			merged[column] = merged[column].fillna(merged[column + "_source"])
		canonical = merged.drop(columns=[column + "_source" for column in ("CHR", "BP", "A1", "A2")])
	if len(set(signatures)) != 1 or not signatures[0]:
		raise ValueError("CSx population weights are not from the same completed joint run")
	return normalize_weight_table(canonical, "Aligned CSx weights"), sources


def canonical_weights(args, trait):
	explicit = optional(option(args, "weights_file"))
	if explicit is not None:
		path = required(expand_trait(explicit, trait))
		frame = normalize_weight_table(read_table(path), str(path))
		for pop in POPS:
			if f"beta_{pop}" not in frame:
				raise ValueError(f"Canonical weights need beta_{pop}")
			frame[f"beta_{pop}"] = numeric(frame[f"beta_{pop}"], f"beta_{pop}")
		sources = [source_stamp(path)]
	else:
		frame, sources = native_weights(args, trait)
	frame = frame.loc[frame.CHR.isin(selected_chromosomes(args))].copy()
	if frame.empty or frame[[f"beta_{pop}" for pop in POPS]].isna().all(axis=1).any():
		raise ValueError("Empty canonical weights or SNP with no source-population effect")
	missing_chromosomes = set(selected_chromosomes(args)) - set(frame.CHR)
	if missing_chromosomes:
		raise ValueError(f"No weights on requested chromosomes {sorted(missing_chromosomes)}; choose an explicit --chrs subset or complete the weights")
	if optional(option(args, "snpinfo")) is not None:
		path = required(option(args, "snpinfo"))
		reference = read_table(path)
		lookup = {str(column).upper().lstrip("#"): column for column in reference}
		if not {"CHR", "BP", "SNP"}.issubset(lookup):
			raise ValueError("SNP reference needs CHR,BP,SNP for exact build validation")
		reference = reference[[lookup[x] for x in ("SNP", "CHR", "BP")]].copy()
		reference.columns = ["SNP", "CHR_ref", "BP_ref"]
		if reference.SNP.duplicated().any():
			raise ValueError("Duplicate IDs in SNP build reference")
		x = frame.merge(reference, on="SNP", how="left", validate="one_to_one")
		match = (x.CHR == pd.to_numeric(x.CHR_ref, errors="coerce")) & (x.BP == pd.to_numeric(x.BP_ref, errors="coerce"))
		if not match.all():
			raise ValueError(f"Weights have {int((~match).sum())} SNPs absent from or inconsistent with the GRCh37 SNP reference")
		sources.append(source_stamp(path))
	return frame.sort_values(["CHR", "BP", "SNP"]).reset_index(drop=True), sources


def validate_target_variants(info, weights, chromosome):
	wanted = set(weights.SNP)
	pieces = []
	with text_stream(info["variant"]) as stream:
		if info["mode"] == "bfile":
			reader = pd.read_csv(stream, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "ALT", "REF"], dtype=str, chunksize=250000)
		else:
			skip = 0
			for line in stream:
				if not line.startswith("##"):
					break
				if re.match(r"##reference=", line, re.I) and re.search(r"GRCh38|hg38", line, re.I):
					raise ValueError("PVAR declares GRCh38 but GRID annotations use GRCh37")
				skip += 1
			# PVAR headers are small; reopening also works for compressed streams.
			reader = None
		if reader is not None:
			for chunk in reader:
				pieces.append(chunk.loc[chunk.SNP.isin(wanted), ["CHR", "BP", "SNP", "REF", "ALT"]])
	if info["mode"] == "pfile":
		with text_stream(info["variant"]) as stream:
			for chunk in pd.read_csv(stream, sep=r"\s+", skiprows=skip, dtype=str, chunksize=250000, usecols=["#CHROM", "POS", "ID", "REF", "ALT"]):
				chunk = chunk.rename(columns={"#CHROM": "CHR", "POS": "BP", "ID": "SNP"})
				pieces.append(chunk.loc[chunk.SNP.isin(wanted)])
	variants = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=["CHR", "BP", "SNP", "REF", "ALT"])
	if variants.SNP.duplicated().any():
		raise ValueError(f"chr{chromosome}: a scored SNP ID is duplicated in the target genotype")
	x = weights.merge(variants, on="SNP", how="left", suffixes=("", "_target"), validate="one_to_one")
	chrom = pd.to_numeric(x.CHR_target.astype(str).str.replace(r"^chr", "", regex=True), errors="coerce")
	pos = pd.to_numeric(x.BP_target, errors="coerce")
	alleles = ((x.A1 == x.REF) & (x.A2 == x.ALT)) | ((x.A1 == x.ALT) & (x.A2 == x.REF))
	match = (chrom == x.CHR) & (pos == x.BP) & alleles
	if not match.all():
		raise ValueError(f"chr{chromosome}: {int((~match).sum())} weighted SNPs have missing/mismatched target ID, coordinate or alleles; first {x.loc[~match, 'SNP'].head().tolist()}")
	return variants


# 🚩 PLINK scoring uses only training frequencies, with strict coverage checks
def run_logged(command, output):
	output = Path(output)
	output.parent.mkdir(parents=True, exist_ok=True)
	log("RUN " + shlex.join(map(str, command)))
	with output.open("w") as stream:
		stream.write(shlex.join(map(str, command)) + "\n")
		stream.flush()
		result = subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT, check=False)
	if result.returncode:
		lines = output.read_text(errors="replace").splitlines()
		raise ValueError(f"Command failed ({result.returncode}); {output}\n" + "\n".join(lines[-16:]))


def read_sscore(path, modules, cohort):
	path = required(path)
	columns = [name + "_SUM" for name in modules]
	with text_stream(path) as stream:
		header = stream.readline().split()
	if set(columns) - set(header) or any(name.endswith("_AVG") for name in header):
		raise ValueError("PLINK did not provide the exact requested score SUM columns")
	id_column = next((name for name in ("#IID", "IID", "eid") if name in header), None)
	if id_column is None:
		raise ValueError("PLINK score header has no IID")
	# Numeric inference prevents tens of millions of temporary Python strings
	# when a full UKB cohort has dozens of module score columns.
	with text_stream(path) as stream:
		frame = pd.read_csv(stream, sep=r"\s+", usecols=[id_column, *columns], dtype={id_column: str})
	frame = identified(frame, "PLINK score")
	if set(frame.eid) != set(cohort.eid):
		raise ValueError("PLINK score sample set differs from the fixed common roster")
	frame = frame.set_index("eid").loc[cohort.eid, columns]
	frame.columns = modules
	for name in modules:
		frame[name] = numeric(frame[name], name)
	if not np.isfinite(frame.to_numpy()).all():
		raise ValueError("Nonfinite PLINK module scores")
	return frame


def check_partition(scores, dictionary):
	for pop in POPS:
		total = scores[f"total_{pop}"].to_numpy(float)
		for group in ("age", "permuted_age", "selection", "differentiation"):
			columns = dictionary.loc[(dictionary.population == pop) & (dictionary.group == group), "module"].tolist()
			if not columns:
				continue
			parts = scores[columns].to_numpy(float)
			error = np.abs(parts.sum(axis=1) - total)
			# PLINK text scores are rounded; tolerance scales with summand magnitudes.
			if np.any(error > 2e-5 * (1 + np.abs(parts).sum(axis=1))):
				raise ValueError(f"Scored {group} modules for {pop} do not reconstruct total_{pop}")


def score_modules(args, directory, weights, cohort, genotypes):
	cache = Path(option(args, "cache_dir", "/tmp/grid-cache/grid"))
	dictionary = read_table(directory / "modules.tsv.gz")
	if not {"module", "group", "population"}.issubset(dictionary):
		raise ValueError("Evolution module dictionary is incomplete")
	modules = dictionary.module.tolist()
	if len(set(modules)) != len(modules):
		raise ValueError("Duplicate score module names")
	plink = str(option(args, "plink2", "plink2"))
	if shutil.which(plink) is None:
		raise ValueError(f"PLINK 2 executable not found: {plink}")
	total = pd.DataFrame(0., index=cohort.eid, columns=modules)
	coverage, frequencies = [], []
	threads = int(option(args, "threads", 4))
	memory = int(option(args, "score_memory", 2048))
	if threads < 1 or memory < 128:
		raise ValueError("PLINK threads must be positive and score_memory >=128 MiB")
	for chrom, subset in weights.groupby("CHR", sort=True):
		info = genotypes[int(chrom)]
		variants = validate_target_variants(info, subset, chrom)
		work = directory / f"chr{chrom}"
		work.mkdir(parents=True, exist_ok=True)
		extract = work / "snps.txt"
		subset.SNP.to_csv(extract, index=False, header=False)
		base = [plink, *info["command"], "--extract", str(extract), "--threads", str(threads), "--memory", str(memory)]
		freq_prefix = work / "train"
		frequency = work / "train.afreq"
		run_logged([*base, "--keep", str(cache / "train.keep"), "--nonfounders", "--freq", "--out", str(freq_prefix)], work / "frequency.command.log")
		freq = read_table(frequency).rename(columns={"ID": "SNP", "#CHROM": "CHR"})
		if not {"SNP", "REF", "ALT", "ALT_FREQS", "OBS_CT"}.issubset(freq) or freq.SNP.duplicated().any() or set(freq.SNP) != set(subset.SNP):
			raise ValueError(f"chr{chrom}: training allele frequencies do not cover all scored SNPs")
		freq = freq.set_index("SNP").loc[subset.SNP].reset_index()
		target_alleles = variants.set_index("SNP").loc[subset.SNP, ["REF", "ALT"]].reset_index()
		if not ((freq.REF == target_alleles.REF) & (freq.ALT == target_alleles.ALT)).all():
			raise ValueError(f"chr{chrom}: training frequency alleles differ from target genotypes")
		altfreq = numeric(freq.ALT_FREQS, "training ALT_FREQS")
		nobs = numeric(freq.OBS_CT, "training OBS_CT")
		if not altfreq.between(0, 1).all() or nobs.isna().any() or (nobs <= 0).any():
			raise ValueError(f"chr{chrom}: unavailable training frequency; no target-test frequency fallback is allowed")
		frequencies.append(freq)
		wide = directory / f"chr{chrom}.weights.tsv.gz"
		with text_stream(required(wide)) as stream:
			wide_head = stream.readline().rstrip("\r\n").split("\t")
		if wide_head != ["SNP", "A1", *modules]:
			raise ValueError(f"chr{chrom}: wide score columns differ from module dictionary")
		prefix = work / "scores"
		run_logged([*base, "--keep", str(cache / "common.keep"), "--read-freq", str(frequency), "--error-on-freq-calc",
			"--score", str(wide), "1", "2", "header-read", "center", "no-mean-imputation", "cols=+scoresums,-scoreavgs", "list-variants",
			"--score-col-nums", f"3-{len(modules) + 2}", "--out", str(prefix)], work / "scoring.command.log")
		used = [line.strip() for line in required(str(prefix) + ".sscore.vars").read_text().splitlines() if line.strip()]
		if used and used[0] in {"ID", "#ID"} and used[0] not in set(subset.SNP):
			used = used[1:]
		if len(used) != len(set(used)) or set(used) != set(subset.SNP):
			raise ValueError(f"chr{chrom}: PLINK skipped or duplicated scored variants")
		scores = read_sscore(str(prefix) + ".sscore", modules, cohort)
		check_partition(scores, dictionary)
		total += scores
		coverage.append({"chromosome": int(chrom), "requested_snps": len(subset), "scored_snps": len(used),
			"train_frequency_n": int((cohort.split == "train").sum()), "n_scored_people": len(scores)})
	check_partition(total, dictionary)
	all_frequencies = pd.concat(frequencies, ignore_index=True)
	all_frequencies.to_csv(directory / "training_frequencies.tsv.gz", sep="\t", index=False)
	total = total.reset_index().rename(columns={"index": "eid"})
	if "eid" not in total:
		total = total.rename(columns={total.columns[0]: "eid"})
	total = total.rename(columns={name: "evo." + name for name in modules})
	total.to_csv(directory / "module_scores.tsv.gz", sep="\t", index=False)
	return total, dictionary, pd.DataFrame(coverage)


def evolution_data(args, prepared):
	if prepared["metadata"]["input_mode"] == "external_prepared_features":
		for trait, frame in prepared["tables"].items():
			prepared["feature_groups"][trait] = feature_groups(frame, prepared["metadata"]["covariates"])
		prepared["metadata"]["evolution_mode"] = "external_prepared_features; provenance and no-test-label construction must be established by provider"
		prepared["metadata"]["evolution_verified"] = False
		prepared["scoring_artifacts"] = {}
		return prepared
	if str(option(args, "build", "GRCh37")) != "GRCh37":
		raise ValueError("GRID requires explicitly harmonized GRCh37 inputs")
	helper = load_module(Path(__file__).with_name("grid.evolution.py"), "grid_data_evolution")
	genotypes = genotype_files(args)
	cache = Path(option(args, "cache_dir", "/tmp/grid-cache/grid"))
	metadata, scoring_artifacts = {}, {}
	for trait, frame in prepared["tables"].items():
		weights, sources = canonical_weights(args, trait)
		directory = cache / "evolution" / trait
		directory.mkdir(parents=True, exist_ok=True)
		weight_path = directory / "canonical_weights.tsv.gz"
		weights.to_csv(weight_path, sep="\t", index=False, float_format="%.17g")
		command = ["build", "--weights", str(weight_path), "--out-dir", str(directory), "--build", "GRCh37",
			"--populations", ",".join(POPS), "--chromosomes", ",".join(map(str, selected_chromosomes(args))),
			"--seed", str(option(args, "seed", 20260904)), "--force"]
		annotations = optional(option(args, "annotation_file"))
		geva = optional(option(args, "geva_dir", "/mnt/f/ref/GEVA"))
		if annotations is not None:
			annotation_path = required(expand_trait(annotations, trait))
			command += ["--annotations", str(annotation_path)]
			sources.append(source_stamp(annotation_path))
		# An explicit canonical file can replace a locally absent GEVA atlas.
		if geva is not None and geva.is_dir():
			command += ["--geva-dir", str(geva)]
			sources += [source_stamp(geva / f"atlas.chr{chrom}.csv.gz") for chrom in sorted(weights.CHR.unique())]
		elif annotations is None and not enabled(option(args, "allow_proxy_only", False)):
			raise ValueError("Missing GEVA atlas/canonical age annotations; acquire source data or explicitly allow proxy-only analysis")
		if enabled(option(args, "allow_proxy_only", False)):
			command.append("--allow-proxy-only")
		if enabled(option(args, "include_proxy", False)):
			command.append("--include-proxy")
		for name in ("min_age_quality", "age_priority", "selection_p", "chunk_size"):
			value = getattr(args, name, None)
			if value is not None:
				command += ["--" + name.replace("_", "-"), str(value)]
		for name in ("age_cutoffs", "diff_cutoffs"):
			value = getattr(args, name, None)
			if value is not None:
				values = names(value)
				if len(values) != 2:
					raise ValueError(f"{name} needs exactly two boundaries")
				command += ["--" + name.replace("_", "-"), *map(str, values)]
		build_result = helper.build(helper.parser().parse_args(command))
		scores, dictionary, coverage = score_modules(args, directory, weights, prepared["cohort"], genotypes)
		frame = frame.drop(columns=[column for column in frame if column.startswith("evo.")], errors="ignore").merge(scores, on="eid", validate="one_to_one")
		prepared["tables"][trait] = frame
		prepared["feature_groups"][trait] = feature_groups(frame, prepared["metadata"]["covariates"], build_result["mode"])
		scoring_artifacts[trait] = {
			"canonical_weights": str(weight_path), "modules": str(directory / "modules.tsv.gz"),
			"training_frequencies": str(directory / "training_frequencies.tsv.gz"),
			"chromosomes": {str(chrom): {
				"wide_weights": str(directory / f"chr{chrom}.weights.tsv.gz"),
				"train_afreq": str(directory / f"chr{chrom}" / "train.afreq"),
				"extract": str(directory / f"chr{chrom}" / "snps.txt"),
			} for chrom in build_result["chromosomes"]},
		}
		metadata[trait] = {"builder": build_result, "weights_sources": sources, "modules": dictionary,
			"coverage": coverage, "annotation_qc": read_table(directory / "qc.tsv"),
			"frequency_source": "fixed_training_half_only", "missing_genotypes": "zero centred contribution",
			"old_grid_transport_or_reshrink_used": False}
	prepared["metadata"]["evolution"] = metadata
	prepared["metadata"]["evolution_verified"] = all(item["builder"]["mode"] == "human_evolution" for item in metadata.values())
	prepared["scoring_artifacts"] = scoring_artifacts
	return prepared
