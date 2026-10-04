#!/usr/bin/env python3
"""0.common.py: prepare-sumstats, publish, io, sumstats-cache, split-sumstats, combine-scores, merge-scores. Use a subcommand followed by --help."""

from __future__ import annotations
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

sys.dont_write_bytecode = True


def load_module(filename):
	"""Load a neighbouring workflow module whose filename contains a dot/digit."""
	key = "grid_" + Path(filename).stem.replace(".", "_")
	if key not in sys.modules:
		spec = spec_from_file_location(key, Path(__file__).with_name(filename))
		module = module_from_spec(spec)
		sys.modules[key] = module
		try:
			spec.loader.exec_module(module)
		except BaseException:
			sys.modules.pop(key, None)
			raise
	return sys.modules[key]


# 🚩 Shared result storage and temporary workspaces
_result_spec = spec_from_file_location("shared_results", Path(__file__).resolve().parents[2] / "0f/results.py")
_result_io = module_from_spec(_result_spec)
_result_spec.loader.exec_module(_result_io)
read_result_table = _result_io.read_table
write_rds = _result_io.write_rds
rds_metadata = _result_io.rds_metadata
write_result_workbook = _result_io.write_workbook


def cache_directory(path):
	path = Path(path).resolve()
	key = hashlib.sha256(str(path).encode()).hexdigest()[:16]
	return Path("/tmp/grid-cache") / key / path.name


# 🚩 prepare-sumstats: prepare_sumstats
"""Normalize heterogeneous GWAS summary statistics to PRS-CSx input."""
import argparse, csv, gzip, json, math, re, hashlib, os
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import norm

# Share the preflight threshold, but exclude every conflicting SNP from inference.
MIN_COORDINATE_MATCH = 0.98
ALIASES = {
	"SNP": ["SNP", "RSID", "RS_ID", "MARKERNAME", "VARIANT_ID", "ID", "RS_NUMBER"],
	"CHR": ["CHR", "CHROM", "CHROMOSOME", "CHROMSOME", "CHR_ID"],
	"BP": ["BP", "POS", "POS_B37", "POSITION", "BASE_PAIR_LOCATION"],
	"A1": ["A1", "EA", "EFFECT_ALLELE", "EFFECTALLELE", "ALT", "ALLELE1", "TESTED_ALLELE", "CODED_ALLELE"],
	"A2": ["A2", "NEA", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "NONEFFECTALLELE", "REF", "ALLELE0", "REFERENCE_ALLELE"],
	"BETA": ["BETA", "EFFECT", "EFFECT_SIZE", "ES", "LOG_ODDS", "B"],
	"OR": ["OR", "ODDS_RATIO", "ODDSRATIO"],
	"SE": ["SE", "STDERR", "STANDARD_ERROR", "BETA_SE"],
	"P": ["P", "PVAL", "PVALUE", "P_VALUE", "P-VALUE"],
	"NEFF": ["NEFF", "N_EFF", "N_EFFECTIVE", "EFFECTIVE_N"],
	"N": ["N", "N_TOTAL", "TOTAL_N", "OBS_CT", "N_SAMPLES"],
	"NCASE": ["N_CASE", "NCASE", "NCASES", "CASES", "N_CASES"],
	"NCTRL": ["N_CONTROL", "N_CONTROLS", "NCTRL", "NCONTROLS", "CONTROLS"],
	"EAF": ["EAF", "AF", "POOLED_ALT_AF", "EFFECT_ALLELE_FREQUENCY", "A1FREQ", "FREQ1", "ALT_FREQ"],
}


def canon(s):
	return re.sub(r"[^A-Z0-9]+", "_", str(s).upper()).strip("_")


def choose(cols, key):
	d = {canon(c): c for c in cols}
	return next((d[canon(x)] for x in ALIASES[key] if canon(x) in d), None)


def read_snpinfo(path):
	d = pd.read_csv(path, sep=r"\s+", engine="python", dtype=str)
	sc = choose(d.columns, "SNP") or d.columns[1 if len(d.columns) > 1 else 0]
	cc = choose(d.columns, "CHR")
	bc = choose(d.columns, "BP")
	snps = set(d[sc].dropna().astype(str))
	locus = {}
	if cc and bc:
		for s, c, b in zip(d[sc], d[cc], d[bc]):
			try:
				locus[(str(c).replace("chr", "").replace("CHR", ""), int(float(b)))] = str(s)
			except:
				pass
	return snps, locus


def sep_for(path):
	op = gzip.open if str(path).endswith(".gz") else open
	with op(path, "rt", errors="replace") as h:
		line = next((x for x in h if x.strip() and not x.startswith("##")), "")
	if "\t" in line:
		return "\t"
	if "," in line:
		return ","
	return r"\s+"


def reader_options(path):
	"""Keep #CHROM as the header; skip only leading ## metadata/blank lines."""
	op = gzip.open if str(path).endswith(".gz") else open
	skip = 0
	with op(path, "rt") as h:
		for line in h:
			if line.startswith("##") or not line.strip():
				skip += 1
			else:
				break
	sep = sep_for(path)
	return dict(sep=sep, engine="c", skiprows=skip, compression="infer", low_memory=False)


def preparation_signature(source, snpinfo, trait, pop):
	def stamp(p):
		p = Path(p)
		s = p.stat()
		return [str(p.resolve()), s.st_size, s.st_mtime_ns]

	return hashlib.sha256(
		json.dumps(
			[
				trait.lower(),
				pop.upper(),
				stamp(source),
				stamp(snpinfo),
				hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
			]
		).encode()
	).hexdigest()


def prepare_sumstats_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--input", required=True)
	ap.add_argument("--output", required=True)
	ap.add_argument("--metadata", required=True)
	ap.add_argument("--snpinfo", required=True)
	ap.add_argument("--trait", required=True)
	ap.add_argument("--pop", required=True)
	ap.add_argument("--chunk", type=int, default=500000)
	a = ap.parse_args()
	if a.chunk < 1:
		raise SystemExit("--chunk must be positive")
	snps, locus = read_snpinfo(a.snpinfo)
	positions = {s: (c, b) for (c, b), s in locus.items()}
	Path(a.output).parent.mkdir(parents=True, exist_ok=True)
	seen = set()
	nvals = []
	stats = {
		"input_rows": 0,
		"kept_rows": 0,
		"bad_allele": 0,
		"ambiguous": 0,
		"missing_effect": 0,
		"not_hm3": 0,
		"duplicates": 0,
		"coordinate_checked": 0,
		"coordinate_mismatch": 0,
	}
	coordinate_examples = []
	first = True
	staging = str(a.output) + f".tmp.{os.getpid()}"
	Path(a.metadata).unlink(missing_ok=True)
	it = pd.read_csv(a.input, chunksize=a.chunk, **reader_options(a.input))
	for d in it:
		stats["input_rows"] += len(d)
		cols = {k: choose(d.columns, k) for k in ALIASES}
		if not cols["A1"] or not cols["A2"]:
			raise SystemExit(f"Allele columns not found in {a.input}: {list(d.columns)}")
		chrom = (
			d[cols["CHR"]]
			.astype(str)
			.str.strip()
			.str.replace(r"^chr", "", case=False, regex=True)
			.str.replace(r"\.0$", "", regex=True)
			if cols["CHR"]
			else pd.Series("", index=d.index)
		)
		bp = pd.to_numeric(d[cols["BP"]], errors="coerce") if cols["BP"] else pd.Series(np.nan, index=d.index)
		sid = d[cols["SNP"]].astype(str) if cols["SNP"] else pd.Series("", index=d.index)
		# Recover rsID from SNPINFO by build-37 position when necessary.
		sid = [
			s if s in snps else locus.get((str(c), int(b)), s) if np.isfinite(b) else s
			for s, c, b in zip(sid, chrom, bp)
		]
		chrom = pd.Series(
			[positions.get(s, (c, None))[0] if c in ("", "nan", "NA") else c for s, c in zip(sid, chrom)], index=d.index
		)
		bp = pd.Series(
			[positions.get(s, (None, b))[1] if not np.isfinite(b) else b for s, b in zip(sid, bp)], index=d.index
		)
		disagreement = [s in positions and (str(c), b) != positions[s] for s, c, b in zip(sid, chrom, bp)]
		stats["coordinate_checked"] += sum(s in positions for s in sid)
		stats["coordinate_mismatch"] += sum(disagreement)
		for s, c, b, bad in zip(sid, chrom, bp, disagreement):
			if bad and len(coordinate_examples) < 5:
				coordinate_examples.append(
					{
						"SNP": s,
						"gwas_chr": str(c),
						"gwas_bp": float(b),
						"reference_chr": positions[s][0],
						"reference_bp": positions[s][1],
					}
				)
		a1 = d[cols["A1"]].astype(str).str.upper()
		a2 = d[cols["A2"]].astype(str).str.upper()
		beta = (
			pd.to_numeric(d[cols["BETA"]], errors="coerce")
			if cols["BETA"]
			else (
				np.log(pd.to_numeric(d[cols["OR"]], errors="coerce"))
				if cols["OR"]
				else pd.Series(np.nan, index=d.index)
			)
		)
		se = pd.to_numeric(d[cols["SE"]], errors="coerce") if cols["SE"] else pd.Series(np.nan, index=d.index)
		pv = pd.to_numeric(d[cols["P"]], errors="coerce") if cols["P"] else pd.Series(np.nan, index=d.index)
		need = se.isna() & beta.notna() & pv.notna() & (pv > 0) & (pv <= 1) & (beta != 0)
		se.loc[need] = np.abs(beta.loc[need]) / norm.isf(np.clip(pv.loc[need] / 2, 1e-323, 0.5))
		if cols["NEFF"]:
			n = pd.to_numeric(d[cols["NEFF"]], errors="coerce")
			n_source = "effective_N"
		elif a.trait.lower() in ("t2dm", "t2d") and cols["NCASE"] and cols["NCTRL"]:
			ca = pd.to_numeric(d[cols["NCASE"]], errors="coerce")
			co = pd.to_numeric(d[cols["NCTRL"]], errors="coerce")
			n = 4 / (1 / ca + 1 / co)
			n = n.where((ca > 0) & (co > 0))
			n_source = "4/(1/Ncase+1/Ncontrol)"
		elif cols["N"]:
			n = pd.to_numeric(d[cols["N"]], errors="coerce")
			n_source = "N (verify effective N for case-control GWAS)"
		else:
			n = pd.Series(np.nan, index=d.index)
			n_source = "unavailable; override required"
		n = n.where(np.isfinite(n) & (n > 0))
		eaf = pd.to_numeric(d[cols["EAF"]], errors="coerce") if cols.get("EAF") else pd.Series(np.nan, index=d.index)
		out = pd.DataFrame(
			{
				"SNP": sid,
				"A1": a1,
				"A2": a2,
				"BETA": beta,
				"SE": se,
				"P": pv,
				"N": n,
				"EAF": eaf,
				"CHR": chrom,
				"BP": bp,
			}
		)
		# Do not rewrite conflicting coordinates or let these rows affect N/deduplication.
		out = out.loc[~np.asarray(disagreement, dtype=bool)]
		valid = out.A1.isin(list("ACGT")) & out.A2.isin(list("ACGT")) & (out.A1 != out.A2)
		stats["bad_allele"] += int((~valid).sum())
		out = out[valid]
		amb = (out.A1 + out.A2).isin(["AT", "TA", "CG", "GC"])
		stats["ambiguous"] += int(amb.sum())
		out = out[~amb]
		eff = np.isfinite(out.BETA) & np.isfinite(out.SE) & (out.SE > 0)
		stats["missing_effect"] += int((~eff).sum())
		out = out[eff]
		hm = out.SNP.isin(snps)
		stats["not_hm3"] += int((~hm).sum())
		out = out[hm]
		dup = out.SNP.isin(seen) | out.SNP.duplicated()
		stats["duplicates"] += int(dup.sum())
		out = out[~dup]
		seen.update(out.SNP.tolist())
		if len(out):
			nvals.extend(out.N.dropna().to_numpy().tolist())
			out.to_csv(
				staging,
				sep="\t",
				index=False,
				mode="wt" if first else "at",
				header=first,
				compression="gzip",
				float_format="%.10g",
			)
			first = False
			stats["kept_rows"] += len(out)
	checked = stats["coordinate_checked"]
	matched = checked - stats["coordinate_mismatch"]
	if checked and matched / checked < MIN_COORDINATE_MATCH:
		Path(staging).unlink(missing_ok=True)
		raise SystemExit(
			f"GRCh37 coordinate check failed in {a.input}: {matched}/{checked} matched; require at least {MIN_COORDINATE_MATCH:.0%}. Examples: {coordinate_examples}"
		)
	if first:
		raise SystemExit(f"No usable HapMap3 variants in {a.input}")
	nmed = float(np.nanmedian(nvals)) if nvals else None
	Path(staging).replace(a.output)
	meta = {
		**stats,
		"trait": a.trait,
		"pop": a.pop,
		"input": str(Path(a.input).resolve()),
		"output": str(Path(a.output).resolve()),
		"n_gwas_median": nmed,
		"n_source": n_source,
		"coordinate_match_min": MIN_COORDINATE_MATCH,
		"coordinate_mismatch_examples": coordinate_examples,
		"preparation_signature": preparation_signature(a.input, a.snpinfo, a.trait, a.pop),
	}
	Path(a.metadata).write_text(json.dumps(meta, indent=2) + "\n")
	print(json.dumps(meta))


def prepare_sumstats_cli():
	prepare_sumstats_main()


# 🚩 publish: score_output
"""Publish compact score tables, excluding withdrawn individuals."""
import argparse
import fcntl, os
from pathlib import Path
import pandas as pd
import numpy as np


def excluded_ids(path):
	if not path:
		return set()
	if not Path(path).is_file():
		raise FileNotFoundError(f"Missing withdrawal list: {path}")
	return {
		parts[1] if len(parts) > 1 else parts[0]
		for line in Path(path).read_text().splitlines()
		if (parts := line.split()) and not parts[0].startswith("#")
	}


def filter_samples(d, idcol, remove):
	ids = d[idcol].astype(str)
	return d.loc[~ids.str.startswith("-") & ~ids.isin(excluded_ids(remove))].copy()


def update_csx_table(d, dest, remove=""):
	"""Atomically replace only supplied score columns, preserving other models."""
	d = d.rename(columns={"IID": "eid"}).copy()
	dest = Path(dest)
	dest.parent.mkdir(parents=True, exist_ok=True)
	lock_path = cache_directory(dest) / "write.lock"
	lock_path.parent.mkdir(parents=True, exist_ok=True)
	with lock_path.open("a") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		if dest.exists():
			old = read_result_table(dest, dtype={"eid": str, "IID": str}).rename(columns={"IID": "eid"})
			old = filter_samples(old, "eid", remove)
			if old.eid.isna().any() or old.eid.duplicated().any():
				raise ValueError("Invalid existing CSx IDs")
			if set(old.eid) != set(d.eid):
				raise ValueError(
					"New and existing CSx sample sets differ; use a separate --score-dir for a different cohort"
				)
			extra = [c for c in old if c != "eid" and c not in d]
			d = d.merge(old[["eid"] + extra], on="eid", validate="one_to_one")
		if d.empty or d.eid.isna().any() or d.eid.duplicated().any():
			raise ValueError("Empty/missing/duplicate CSx IDs")
		if not np.isfinite(d.drop(columns="eid").to_numpy(float)).all():
			raise ValueError("Nonfinite CSx scores")
		write_rds(d, dest)
	return d


def publish_table(src, dest, method, remove):
	d = read_result_table(src, dtype={"eid": str, "IID": str})
	idcol = "eid" if "eid" in d else "IID"
	before = len(d)
	d = filter_samples(d, idcol, remove)
	mapping = {f"CSX_{p}": f"csx.{p}" for p in ["AFR", "EAS", "EUR", "SAS"]}
	if method == "disco":
		mapping = {"PRS": "disco"}
	d = d.rename(columns=mapping)
	if d.empty or d[idcol].isna().any() or d[idcol].duplicated().any():
		raise ValueError("Empty/missing/duplicate score IDs")
	if method == "csx":
		update_csx_table(d, dest, remove)
		print(f"{method}: retained={len(d)}; excluded={before - len(d)}; output={dest}")
		return
	if not np.isfinite(d.drop(columns=idcol).to_numpy(float)).all():
		raise ValueError("Nonfinite scores")
	dest = Path(dest)
	dest.parent.mkdir(parents=True, exist_ok=True)
	write_rds(d, dest)
	print(f"{method}: retained={len(d)}; excluded={before - len(d)}; output={dest}")


def score_output_cli():
	p = argparse.ArgumentParser()
	p.add_argument("method", choices=["csx", "disco"])
	p.add_argument("input")
	p.add_argument("output")
	p.add_argument("--remove", default="/mnt/d/files/ukb.exclude.id")
	a = p.parse_args()
	publish_table(a.input, a.output, a.method, a.remove)


# 🚩 io: pipeline_io
"""Validation, provenance and atomic-output staging for the three launchers."""
import argparse, gzip, hashlib, json, os
from pathlib import Path
import numpy as np
import pandas as pd


def stamp(path):
	p = Path(path)
	if not p.is_file():
		raise ValueError(f"Missing file: {p}")
	z = p.stat()
	return [str(p.resolve()), z.st_size, z.st_mtime_ns]


def inspect(snpinfo, files):
	ref = pd.read_csv(snpinfo, sep=r"\s+", usecols=["SNP", "CHR", "BP"], dtype={"SNP": str})
	ref = ref.drop_duplicates("SNP").set_index("SNP")
	for path in files:
		marker = Path(path + ".grch")
		if marker.exists() and marker.read_text().strip() != "37":
			raise ValueError(f"GRCh37 required: {marker}")
		d = pd.read_csv(path, nrows=100000, **reader_options(path))
		columns = {key: choose(d.columns, key) for key in ("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "N")}
		missing = [key for key in ("A1", "A2", "SE") if columns[key] is None]
		if columns["BETA"] is None and choose(d.columns, "OR") is None:
			missing.append("BETA or OR")
		if columns["SNP"] is None and (columns["CHR"] is None or columns["BP"] is None):
			missing.append("rsID or CHR/BP")
		if missing:
			raise ValueError(f"{path}: missing {missing}")
		if columns["SNP"] is None or columns["CHR"] is None or columns["BP"] is None:
			print(f"input GWAS: {path}; coordinate completion uses the explicitly supplied GRCh37 SNPINFO")
			continue
		check = pd.DataFrame(
			{
				"SNP": d[columns["SNP"]].astype(str),
				"CHR": pd.to_numeric(
					d[columns["CHR"]].astype(str).str.strip().str.replace(r"^chr", "", regex=True), errors="coerce"
				),
				"POS": pd.to_numeric(d[columns["BP"]], errors="coerce"),
			}
		)
		z = check.merge(ref, on="SNP", suffixes=("", "_ref"))
		ok = (z.CHR == z.CHR_ref) & (z.POS == z.BP)
		n = len(z)
		rate = float(ok.mean()) if n else 0
		if n < min(100, len(ref)) or rate < MIN_COORDINATE_MATCH:
			raise ValueError(f"{path}: GRCh37 reference-coordinate check failed ({int(ok.sum())}/{n})")
		print(
			f"input GWAS: {path}\n  GRCh37 coordinate sample: {int(ok.sum())}/{n}; N is median of usable HM3 variants (override: --n-gwas)"
		)


def coverage(chrs, files):
	"""Read the BGZF index when available; plain raw gzip is checked after preparation."""
	import struct

	wanted = set(chrs.split())
	for path in files:
		index = Path(path + ".tbi")
		if not index.is_file():
			continue
		with gzip.open(index, "rb") as stream:
			if stream.read(4) != b"TBI\x01":
				raise ValueError(f"Invalid tabix index: {index}")
			header = stream.read(32)
			if len(header) != 32:
				raise ValueError(f"Truncated tabix index: {index}")
			length = struct.unpack("<8i", header)[7]
			names = stream.read(length).decode().strip("\x00").split("\x00")
		present = {s.removeprefix("chr") for s in names}
		missing = sorted(wanted - present, key=int)
		if missing:
			raise ValueError(
				f"{path}: source GWAS is missing chromosomes {','.join(missing)} (tabix index). "
				"Repair the formatted common GWAS from the original source data first; "
				"do not silently omit missing chromosomes."
			)


def signature(values, files):
	obj = {"settings": values, "files": [stamp(p) for p in files]}
	return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def weights(inputs, output):
	d = pd.concat([pd.read_csv(p, sep="\t") for p in inputs], ignore_index=True)
	if d.empty or d.SNP.duplicated().any() or not np.isfinite(d.BETA).all():
		raise ValueError("Invalid/duplicate posterior weights")
	d.to_csv(output, sep="\t", index=False, compression="gzip")
	print(f"posterior variants: {len(d)}")


def disco_inputs(pca, centers, score_dir, outdir, npc, remove="/mnt/d/files/ukb.exclude.id"):
	if not 5 <= npc <= 20:
		raise ValueError("Disco distance PCs must be 5..20")
	pc = [f"PC{i}" for i in range(1, npc + 1)]
	d = pd.read_csv(pca, sep="\t", dtype={"IID": str, "eid": str, "#IID": str})
	idcol = next(x for x in ("IID", "eid", "#IID") if x in d)
	d = d.rename(columns={idcol: "IID"})[["IID"] + pc]
	if d.IID.isna().any() or d.IID.duplicated().any() or not np.isfinite(d[pc].to_numpy()).all():
		raise ValueError("Invalid PCA IDs or PCs")
	med = pd.read_csv(centers, sep="\t")
	med = med[[med.columns[0]] + pc]
	pops = ["AFR", "EAS", "EUR", "SAS"]
	if list(med.iloc[:, 0]) != pops or not np.isfinite(med[pc].to_numpy()).all():
		raise ValueError("Centers must have finite PCs in AFR,EAS,EUR,SAS order")
	scores = []
	merged = Path(score_dir) / "1csx.scores.rds"
	combined = read_result_table(merged, dtype={"eid": str})
	for pop in pops:
		p = merged
		source = combined
		z = source.rename(columns={"eid": "IID", f"CSX_{pop}": "PRS", f"csx.{pop}": "PRS"})[["IID", "PRS"]]
		if z.IID.isna().any() or z.IID.duplicated().any() or not np.isfinite(z.PRS).all():
			raise ValueError(f"Invalid scores: {p}")
		if z.PRS.std() == 0:
			raise ValueError(f"Constant scores: {p}")
		if scores and set(z.IID) != set(scores[0].IID):
			raise ValueError("Population PRS sample sets differ")
		scores.append(z)
	before = len(scores[0])
	scores = [filter_samples(z, "IID", remove) for z in scores]
	print(f"Disco withdrawn filter: excluded={before - len(scores[0])}; retained={len(scores[0])}")
	ids = set(scores[0].IID)
	retained = ids.intersection(d.IID)
	print(f"Disco sample filter: scored={len(ids)}; missing PCA={len(ids - retained)}; retained={len(retained)}")
	if len(retained) < 2:
		raise ValueError("Fewer than two scored samples remain after PCA filtering")
	d = d[d.IID.isin(retained)].sort_values("IID")
	scores = [z[z.IID.isin(retained)].sort_values("IID") for z in scores]
	for pop, z in zip(pops, scores):
		if z.PRS.std() == 0:
			raise ValueError(f"Constant {pop} scores after PCA filtering")
	for _, row in med.iterrows():
		dist = np.linalg.norm(d[pc].to_numpy() - row[pc].to_numpy(dtype=float), axis=1)
		if (dist <= 0).any():
			raise ValueError("Sample exactly at reference center; official interpolation is undefined")
	target = Path(outdir)
	target.mkdir(parents=True, exist_ok=True)
	d.to_csv(target / "pca.tsv", sep="\t", index=False)
	med.to_csv(target / "centers.tsv", sep="\t", index=False)
	for pop, z in zip(pops, scores):
		z.to_csv(target / f"{pop}.tsv", sep="\t", index=False)
	print(f"Disco aligned samples: {len(d)}; distance PCs: {npc}")


def validate_disco(path, inputs):
	d = pd.read_csv(path, sep="\t", dtype={"IID": str})
	ref = pd.read_csv(Path(inputs) / "AFR.tsv", sep="\t", dtype={"IID": str})
	if d.empty or d.IID.duplicated().any() or set(d.IID) != set(ref.IID) or not np.isfinite(pd.to_numeric(d.PRS)).all():
		raise ValueError("Official DiscoDivas returned invalid/incomplete scores")
	co = pd.read_csv(path.replace(".tsv.gz", ".coef.tsv.gz"), sep="\t", dtype={"IID": str})
	vals = co.drop(columns="IID").to_numpy()
	if (
		co.IID.duplicated().any()
		or set(co.IID) != set(ref.IID)
		or not np.isfinite(vals).all()
		or not np.allclose(vals.sum(axis=1), 1)
	):
		raise ValueError("Invalid interpolation coefficients")
	print(f"validated Disco scores: {len(d)}")


def pipeline_io_cli():
	import sys

	action, *args = sys.argv[1:]
	try:
		if action == "inspect":
			inspect(args[0], args[1:])
		elif action == "coverage":
			coverage(args[0], args[1:])
		elif action == "signature":
			split = args.index("--files")
			print(signature(args[:split], args[split + 1 :]))
		elif action == "weights":
			weights(args[1:], args[0])
		elif action == "disco-inputs":
			disco_inputs(*args[:4], int(args[4]), *args[5:])
		elif action == "disco-output":
			validate_disco(*args)
		else:
			raise ValueError(f"Unknown action: {action}")
	except (ValueError, KeyError, OSError) as e:
		raise SystemExit(f"ERROR: {e}")


# 🚩 sumstats-cache: sumstats_cache
"""Store normalized GWAS beside the source; stage independent working copies."""
import argparse
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


FORMAT_VERSION = 1
COLUMNS = ["SNP", "A1", "A2", "BETA", "SE", "P", "N", "EAF", "CHR", "BP"]


def cache_paths(source):
	source = Path(source).resolve()
	stem = source.name[:-3] if source.name.endswith(".gz") else source.name
	prefix = source.parent / (stem + ".csx.sumstats")
	return Path(str(prefix) + ".gz"), Path(str(prefix) + ".json")


def read_metadata(path):
	try:
		data = json.loads(Path(path).read_text())
		return data if isinstance(data, dict) else None
	except (OSError, ValueError):
		return None


def sha256(path):
	result = hashlib.sha256()
	with Path(path).open("rb") as stream:
		for block in iter(lambda: stream.read(1024 * 1024), b""):
			result.update(block)
	return result.hexdigest()


def validate_table(table, metadata):
	"""Check the entire gzip stream and row count before adopting/publishing it."""
	with gzip.open(table, "rt") as stream:
		if stream.readline().rstrip("\r\n").split("\t") != COLUMNS:
			raise ValueError(f"Not a normalized GWAS table: {table}")
		rows = sum(1 for line in stream if line.strip())
	if rows <= 0 or rows != metadata.get("kept_rows"):
		raise ValueError(f"Incomplete normalized GWAS: {table} ({rows} rows)")
	return sha256(table)


def cached_metadata(table, metadata, key):
	data = read_metadata(metadata)
	if not data or data.get("preparation_signature") != key or data.get("cache_format") != FORMAT_VERSION:
		return None
	try:
		stat = table.stat()
		if stat.st_size != data["table_size"] or stat.st_mtime_ns != data["table_mtime_ns"] or stat.st_size == 0:
			return None
		with gzip.open(table, "rt") as stream:
			if stream.readline().rstrip("\r\n").split("\t") != COLUMNS:
				return None
	except (OSError, EOFError, KeyError, ValueError):
		return None
	return data


def legacy_candidates(work, trait, pop, key, source):
	if work is None:
		return []
	found = []
	for path in sorted((Path(work) / trait).glob(f"*/sumstats/{pop}.json")):
		data = read_metadata(path)
		table = path.with_suffix(".tsv.gz")
		if data and data.get("preparation_signature") == key and data.get("input") == str(source) and table.is_file():
			found.append((table, path, data))
	return found


def atomic_json(path, data):
	path = Path(path)
	path.parent.mkdir(parents=True, exist_ok=True)
	fd, temporary = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
	try:
		with os.fdopen(fd, "w") as stream:
			json.dump(data, stream, indent=2)
			stream.write("\n")
		os.replace(temporary, path)
	finally:
		Path(temporary).unlink(missing_ok=True)


def stage_working_copy(table, data, output, metadata):
	"""Keep per-run phi/N annotations and readers separate from the shared cache."""
	output, metadata = Path(output), Path(metadata)
	if output.resolve() == table.resolve():
		raise ValueError("Working output must differ from the permanent cache")
	output.parent.mkdir(parents=True, exist_ok=True)
	old = read_metadata(metadata)
	same = (
		old
		and old.get("preparation_signature") == data["preparation_signature"]
		and old.get("table_sha256") == data["table_sha256"]
		and output.is_file()
		and output.stat().st_size == data["table_size"]
		and output.stat().st_mtime_ns == data["table_mtime_ns"]
	)
	if not same:
		fd, temporary = tempfile.mkstemp(prefix="." + output.name + ".", dir=output.parent)
		os.close(fd)
		try:
			shutil.copy2(table, temporary)
			os.replace(temporary, output)
		finally:
			Path(temporary).unlink(missing_ok=True)
	run_data = dict(data, output=str(output.resolve()), permanent_output=str(table))
	atomic_json(metadata, run_data)


def ensure_prepared(
	source,
	snpinfo,
	trait,
	pop,
	*,
	work=None,
	output=None,
	metadata=None,
	chunk=500000,
	replace=False,
	migrate_only=False,
	remove_legacy=False,
):
	source, snpinfo = Path(source).resolve(), Path(snpinfo).resolve()
	trait, pop = trait.lower(), pop.upper()
	if bool(output) != bool(metadata):
		raise ValueError("Working output and metadata must be supplied together")
	if remove_legacy and not migrate_only:
		raise ValueError("Removing old tables requires --migrate-only")
	table, info = cache_paths(source)
	key = preparation_signature(source, snpinfo, trait, pop)
	lock_path = info.with_suffix(".lock")
	with lock_path.open("a") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		data = None if replace else cached_metadata(table, info, key)
		candidates = legacy_candidates(work, trait, pop, key, source)
		action = "SKIP"
		if data is None:
			chosen = None
			if not replace:
				for old_table, old_info, old_data in candidates:
					try:
						digest = validate_table(old_table, old_data)
					except (OSError, EOFError, ValueError):
						continue
					chosen = (old_table, old_data, digest)
					break
			if chosen is None and migrate_only:
				print(f"MISSING preprocessing {trait}.{pop}: no complete matching table", flush=True)
				return None
			with tempfile.TemporaryDirectory(prefix="." + table.name + ".", dir=table.parent) as temp:
				temporary = Path(temp) / "sumstats.tsv.gz"
				temporary_info = Path(temp) / "sumstats.json"
				if chosen:
					old_table, data, digest = chosen
					shutil.copy2(old_table, temporary)
					if sha256(temporary) != digest:
						raise ValueError(f"Cache copy verification failed: {old_table}")
					action = "MIGRATED"
				else:
					print(f"RUN preprocessing {trait}.{pop}: {source}", flush=True)
					subprocess.run(
						[
							sys.executable,
							str(Path(__file__).with_name("0.common.py")),
							"prepare-sumstats",
							"--input",
							str(source),
							"--output",
							str(temporary),
							"--metadata",
							str(temporary_info),
							"--snpinfo",
							str(snpinfo),
							"--trait",
							trait,
							"--pop",
							pop,
							"--chunk",
							str(chunk),
						],
						check=True,
					)
					data = read_metadata(temporary_info)
					if not data or data.get("preparation_signature") != key:
						raise ValueError("Preprocessing inputs changed during preparation")
					digest = validate_table(temporary, data)
					action = "SAVED"
				if preparation_signature(source, snpinfo, trait, pop) != key:
					raise ValueError("Preprocessing inputs changed during cache publication")
				data = dict(data)
				# These settings belong to a particular inference run, not normalization.
				for name in ("inference_phi", "n_gwas_used", "permanent_output"):
					data.pop(name, None)
				stat = temporary.stat()
				data.update(
					output=str(table),
					cache_format=FORMAT_VERSION,
					table_sha256=digest,
					table_size=stat.st_size,
					table_mtime_ns=stat.st_mtime_ns,
				)
				# The JSON is the commit marker; interrupted publication is never a cache hit.
				info.unlink(missing_ok=True)
				os.replace(temporary, table)
				atomic_json(info, data)
		if remove_legacy:
			for old_table, old_info, old_data in candidates:
				if old_table.resolve() == table.resolve():
					continue
				if old_table.stat().st_size == data["table_size"] and sha256(old_table) == data["table_sha256"]:
					old_table.unlink()
					# Retain per-run metadata (including phi/N) as provenance.
					old_data["output"] = str(table)
					old_data["permanent_output"] = str(table)
					atomic_json(old_info, old_data)
		if output:
			stage_working_copy(table, data, output, metadata)
		print(f"{action} preprocessing {trait}.{pop}: {table}", flush=True)
		return table, info, data


def sumstats_cache_main():
	parser = argparse.ArgumentParser(description=__doc__)
	for key in ("input", "snpinfo", "trait", "pop"):
		parser.add_argument("--" + key, required=True)
	for key in ("work", "output", "metadata"):
		parser.add_argument("--" + key)
	parser.add_argument("--chunk", type=int, default=500000)
	parser.add_argument("--replace", choices=("TRUE", "FALSE"), default="FALSE")
	parser.add_argument("--migrate-only", action="store_true")
	parser.add_argument("--remove-legacy", action="store_true")
	args = parser.parse_args()
	if args.chunk < 1:
		parser.error("--chunk must be positive")
	ensure_prepared(
		args.input,
		args.snpinfo,
		args.trait,
		args.pop,
		work=args.work,
		output=args.output,
		metadata=args.metadata,
		chunk=args.chunk,
		replace=args.replace == "TRUE",
		migrate_only=args.migrate_only,
		remove_legacy=args.remove_legacy,
	)


def sumstats_cache_cli():
	sumstats_cache_main()


# 🚩 split-sumstats: split_sumstats
import argparse, gzip
from pathlib import Path
import pandas as pd


def split_sumstats_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--input", required=True)
	ap.add_argument("--out-dir", required=True)
	ap.add_argument("--prefix", required=True)
	ap.add_argument("--chrs", required=True)
	a = ap.parse_args()
	chrs = {str(x) for x in a.chrs.replace(",", " ").split()}
	Path(a.out_dir).mkdir(parents=True, exist_ok=True)
	first = {c: True for c in chrs}
	for d in pd.read_csv(a.input, sep="\t", compression="infer", dtype={"CHR": str}, chunksize=200000):
		d["CHR"] = d.CHR.str.replace("chr", "", case=False, regex=False)
		for c, x in d[d.CHR.isin(chrs)].groupby("CHR"):
			p = f"{a.out_dir}/{a.prefix}.chr{c}.tsv"
			x.to_csv(p, sep="\t", index=False, mode="wt" if first[c] else "at", header=first[c])
			first[c] = False
	for c in chrs:
		if first[c]:
			raise SystemExit(f"No variants for chromosome {c}")


def split_sumstats_cli():
	split_sumstats_main()


# 🚩 combine-scores: combine_scores
import argparse, gzip, re
from pathlib import Path
import pandas as pd
import numpy as np


def read_score(path):
	d = pd.read_csv(path, sep=r"\s+", dtype=str)
	idc = next((c for c in ["IID", "#IID", "ID_2", "eid"] if c in d.columns), None)
	if idc is None:
		raise SystemExit(f"No IID in {path}")
	cols = [c for c in d.columns if c.endswith("_SUM") and c not in {"NAMED_ALLELE_DOSAGE_SUM"}]
	if len(cols) != 1:
		raise SystemExit(f"Expected exactly one score SUM in {path}: {list(d.columns)}")
	if d[idc].isna().any():
		raise SystemExit(f"Missing score IDs in {path}")
	return pd.DataFrame({"eid": d[idc].astype(str), "score": pd.to_numeric(d[cols[-1]], errors="raise")})


def combine_scores_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--inputs", nargs="+", required=True)
	ap.add_argument("--name", required=True)
	ap.add_argument("--output", required=True)
	a = ap.parse_args()
	z = None
	for p in a.inputs:
		d = read_score(p)
		if d.eid.duplicated().any() or not np.isfinite(d.score).all():
			raise SystemExit(f"Invalid/duplicate score IDs in {p}")
		if z is not None and set(z.eid) != set(d.eid):
			raise SystemExit(f"Chromosome sample sets differ: {p}")
		z = (
			d
			if z is None
			else z.merge(d, on="eid", how="inner", suffixes=("", "_x"), validate="one_to_one")
			.assign(score=lambda x: x["score"] + x["score_x"])
			.drop(columns="score_x")
		)
	z = z.rename(columns={"score": a.name})
	z.to_csv(a.output, sep="\t", index=False, compression="gzip")
	print(len(z))


def combine_scores_cli():
	combine_scores_main()


# 🚩 merge-scores: merge_scores
import argparse
from pathlib import Path
import pandas as pd


def merge_scores_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--inputs", nargs="+", required=True)
	ap.add_argument("--output", required=True)
	a = ap.parse_args()
	z = None
	for p in a.inputs:
		d = pd.read_csv(p, sep="\t", compression="infer", dtype={"eid": str})
		if d.eid.isna().any() or d.eid.duplicated().any():
			raise ValueError(f"Missing/duplicate sample IDs: {p}")
		if z is not None and (set(z.columns) & set(d.columns)) != {"eid"}:
			raise ValueError(f"Duplicate score columns: {p}")
		z = d if z is None else z.merge(d, on="eid", how="outer", validate="one_to_one")
	Path(a.output).parent.mkdir(parents=True, exist_ok=True)
	z.to_csv(a.output, sep="\t", index=False, compression="gzip")
	print(len(z))


def merge_scores_cli():
	merge_scores_main()


def cache_path_cli():
	if len(sys.argv) != 2:
		raise SystemExit("Usage: 0.common.py cache-path PATH")
	print(cache_directory(sys.argv[1]))


COMMANDS = {
	"cache-path": cache_path_cli,
	"prepare-sumstats": prepare_sumstats_cli,
	"publish": score_output_cli,
	"io": pipeline_io_cli,
	"sumstats-cache": sumstats_cache_cli,
	"split-sumstats": split_sumstats_cli,
	"combine-scores": combine_scores_cli,
	"merge-scores": merge_scores_cli,
}




def main():
	import sys

	if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help", "help"):
		print(__doc__)
		print("Commands: " + ", ".join(COMMANDS))
		return
	command = sys.argv.pop(1)
	if command in ("inspect", "coverage", "signature", "weights", "disco-inputs", "disco-output"):
		sys.argv.insert(1, command)
		return pipeline_io_cli()
	if command not in COMMANDS:
		raise SystemExit("Unknown command: " + command)
	sys.argv[0] += " " + command
	COMMANDS[command]()


if __name__ == "__main__":
	main()
