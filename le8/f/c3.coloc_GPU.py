#!/usr/bin/env python3


# 🚩 c3.coloc_GPU
import csv
import hashlib
import importlib.metadata
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

try:
	import numpy as np
	import pandas as pd
	import pyarrow.feather as feather
except Exception as e:
	raise SystemExit(
		"This script needs pandas, numpy and pyarrow in the active Python environment. Install with: python -m pip install pandas numpy pyarrow\n"
		+ str(e)
	)


def opener(path, mode="rt"):
	return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def norm(x):
	return re.sub(r"[^A-Z0-9]+", "", str(x).upper())


def norm_chrom(x):
	chrom = re.sub(r"^chr", "", str(x), flags=re.IGNORECASE).upper()
	chrom = re.sub(r"\.0$", "", chrom)
	chrom = {"X": "23", "Y": "24", "M": "25", "MT": "25"}.get(chrom, chrom)
	return chrom.lstrip("0") or "0"


def pick(cols, names):
	mapping = {norm(c): c for c in cols}
	for n in names:
		if norm(n) in mapping:
			return mapping[norm(n)]
	return None


def safe_name(x):
	return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


def read_header(path):
	with opener(path) as f:
		line = f.readline()
	delim = "\t" if line.count("\t") >= line.count(",") else ","
	cols = line.rstrip("\n\r").split(delim)
	return cols, delim


def sumstat_reader(path):
	cols, delim = read_header(path)
	c_snp = pick(cols, ["SNP", "RSID", "RS_NUMBER", "VARIANT_ID", "ID", "MARKERNAME"])
	c_chr = pick(cols, ["CHR", "CHROM", "CHROMOSOME", "chromosome"])
	c_pos = pick(cols, ["POS", "BP", "POSITION", "BASE_PAIR_LOCATION", "position"])
	c_ea = pick(cols, ["EA", "A1", "ALT", "ALLELE1", "EFFECT_ALLELE", "effect_allele"])
	c_nea = pick(cols, ["NEA", "A2", "REF", "ALLELE0", "OTHER_ALLELE", "reference_allele"])
	# GPU-coloc receives the full regional QTL/GWAS files, so prefer marginal effects.
	c_beta = pick(cols, ["BETA", "beta", "B", "EFFECT", "LOG_ODDS", "BJ", "bJ"])
	c_se = pick(cols, ["SE", "se", "SEBETA", "STDERR", "BJ_SE", "bJ_se"])
	c_p = pick(cols, ["P", "pval", "PVALUE", "P_VALUE", "PJ", "pJ"])
	required = [c_chr, c_pos, c_ea, c_nea, c_beta, c_se]
	if any(x is None for x in required):
		raise RuntimeError(f"Missing required columns in {path}. Need CHR POS EA NEA BETA SE.")
	usecols = [x for x in [c_snp, c_chr, c_pos, c_ea, c_nea, c_beta, c_se, c_p] if x]
	return {
		"delimiter": delim,
		"usecols": set(usecols),
		"snp": c_snp,
		"chr": c_chr,
		"pos": c_pos,
		"ea": c_ea,
		"nea": c_nea,
		"beta": c_beta,
		"se": c_se,
		"p": c_p,
	}


def normalize_sumstat_chunk(ch, spec):
	ch = ch.rename(
		columns={
			spec["chr"]: "CHR",
			spec["pos"]: "POS",
			spec["ea"]: "EA",
			spec["nea"]: "NEA",
			spec["beta"]: "BETA",
			spec["se"]: "SE",
		}
	)
	if spec["snp"]:
		ch = ch.rename(columns={spec["snp"]: "SNP"})
	else:
		ch["SNP"] = ""
	if spec["p"]:
		ch = ch.rename(columns={spec["p"]: "P"})
	else:
		ch["P"] = pd.NA
	ch["CHR"] = ch["CHR"].map(norm_chrom)
	ch["POS"] = pd.to_numeric(ch["POS"], errors="coerce")
	ch["BETA"] = pd.to_numeric(ch["BETA"], errors="coerce")
	ch["SE"] = pd.to_numeric(ch["SE"], errors="coerce")
	ch["P"] = pd.to_numeric(ch["P"], errors="coerce")
	ch = ch.dropna(subset=["CHR", "POS", "BETA", "SE"])
	ch = ch[ch["SE"] > 0]
	ch["EA"] = ch["EA"].astype(str).str.upper()
	ch["NEA"] = ch["NEA"].astype(str).str.upper()
	return ch[(ch["EA"] != "") & (ch["NEA"] != "")]


def iter_sumstat_chunks(path, chunksize=500000):
	spec = sumstat_reader(path)
	for ch in pd.read_csv(
		path,
		sep=spec["delimiter"],
		usecols=lambda c: c in spec["usecols"],
		chunksize=chunksize,
		compression="infer",
	):
		ch = normalize_sumstat_chunk(ch, spec)
		if not ch.empty:
			yield ch


def tabix_index(path):
	for suffix in (".tbi", ".csi"):
		index = str(path) + suffix
		if os.path.isfile(index) and os.path.getsize(index) > 0 and os.path.getmtime(index) >= os.path.getmtime(path):
			return index
	return None


def load_regions_tabix(path, interval_map, label=None, chunksize=500000):
	tabix = shutil.which("tabix")
	if not tabix or not tabix_index(path):
		return None
	contigs = subprocess.check_output([tabix, "-l", str(path)], text=True).splitlines()
	contig_map = {norm_chrom(x): x for x in contigs}
	region_args = []
	for chrom, (starts, ends) in interval_map.items():
		source_chrom = contig_map.get(norm_chrom(chrom))
		if source_chrom is None:
			continue
		region_args.extend(f"{source_chrom}:{int(start)}-{int(end)}" for start, end in zip(starts, ends))
	if label:
		print(
			f"{label}: tabix query {path} for {len(region_args)} merged intervals",
			flush=True,
		)
	if not region_args:
		return empty_sumstats()
	spec = sumstat_reader(path)
	cols, _ = read_header(path)
	proc = subprocess.Popen(
		[tabix, str(path), *region_args],
		stdout=subprocess.PIPE,
		stderr=subprocess.DEVNULL,
		text=True,
	)
	chunks = []
	try:
		for ch in pd.read_csv(
			proc.stdout,
			sep=spec["delimiter"],
			names=cols,
			header=None,
			usecols=lambda c: c in spec["usecols"],
			chunksize=chunksize,
		):
			ch = normalize_sumstat_chunk(ch, spec)
			if not ch.empty:
				chunks.append(ch)
	except pd.errors.EmptyDataError:
		pass
	finally:
		if proc.stdout is not None:
			proc.stdout.close()
	rc = proc.wait()
	if rc:
		raise RuntimeError(f"tabix exited with status {rc}: {path}")
	result = pd.concat(chunks, ignore_index=True) if chunks else empty_sumstats()
	if label:
		print(f"{label}: tabix retained {len(result):,} usable rows", flush=True)
	return result


def empty_sumstats():
	return pd.DataFrame(columns=["SNP", "CHR", "POS", "EA", "NEA", "BETA", "SE", "P"])


def merge_regions(regions):
	"""Return non-overlapping intervals grouped by canonical chromosome."""
	grouped = {}
	for chrom, start, end in regions:
		chrom = norm_chrom(chrom)
		start, end = int(float(start)), int(float(end))
		if end < start:
			start, end = end, start
		grouped.setdefault(chrom, []).append((start, end))
	merged = {}
	for chrom, intervals in grouped.items():
		out = []
		for start, end in sorted(intervals):
			if out and start <= out[-1][1] + 1:
				out[-1] = (out[-1][0], max(out[-1][1], end))
			else:
				out.append((start, end))
		merged[chrom] = (
			np.asarray([x[0] for x in out], dtype=np.int64),
			np.asarray([x[1] for x in out], dtype=np.int64),
		)
	return merged


def load_regions(path, regions, label=None):
	"""Scan a summary-statistics file once and retain the union of regions.

	The input may be ordinary gzip, so this is the bounded-memory fallback for
	files that cannot be queried with tabix. Merged intervals and binary
	searches avoid testing every input row against every requested locus.
	"""
	interval_map = merge_regions(regions)
	if not interval_map:
		return empty_sumstats()
	try:
		indexed = load_regions_tabix(path, interval_map, label=label)
		if indexed is not None:
			return indexed
	except Exception as exc:
		print(
			f"WARNING: tabix query failed; falling back to one full scan: {exc}",
			file=sys.stderr,
			flush=True,
		)
	if label:
		n_intervals = sum(len(x[0]) for x in interval_map.values())
		print(
			f"{label}: scanning {path} once for {n_intervals} merged intervals",
			flush=True,
		)
	chunks = []
	input_rows = 0
	for ch in iter_sumstat_chunks(path):
		input_rows += len(ch)
		chrom_values = ch["CHR"].to_numpy(dtype=str)
		positions = ch["POS"].to_numpy(dtype=np.float64)
		keep = np.zeros(len(ch), dtype=bool)
		for chrom in np.unique(chrom_values):
			interval = interval_map.get(chrom)
			if interval is None:
				continue
			row_index = np.flatnonzero(chrom_values == chrom)
			pos = positions[row_index]
			starts, ends = interval
			interval_index = np.searchsorted(starts, pos, side="right") - 1
			valid = interval_index >= 0
			selected = np.zeros(len(row_index), dtype=bool)
			selected[valid] = pos[valid] <= ends[interval_index[valid]]
			keep[row_index[selected]] = True
		if keep.any():
			chunks.append(ch.loc[keep].copy())
	result = pd.concat(chunks, ignore_index=True) if chunks else empty_sumstats()
	if label:
		print(
			f"{label}: retained {len(result):,} of {input_rows:,} usable rows",
			flush=True,
		)
	return result


def select_region(d, chrom, start, end):
	if d.empty:
		return d.copy()
	chrom = norm_chrom(chrom)
	return d[(d["CHR"].astype(str) == chrom) & (d["POS"] >= start) & (d["POS"] <= end)].copy()


def load_region(path, chrom=None, start=None, end=None):
	if chrom is not None:
		return select_region(load_regions(path, [(chrom, start, end)]), chrom, start, end)
	chunks = list(iter_sumstat_chunks(path))
	return pd.concat(chunks, ignore_index=True) if chunks else empty_sumstats()


def find_lead_and_region(path, region, window_kb):
	if region and re.match(r"^[A-Za-z0-9]+:\d+-\d+$", region):
		chrom, rest = region.split(":", 1)
		start, end = rest.split("-", 1)
		return (
			re.sub(r"^chr", "", chrom, flags=re.IGNORECASE),
			int(float(start)),
			int(float(end)),
		)
	d = load_region(path)
	if d.empty:
		raise RuntimeError("No usable QTL rows in " + path)
	if d["P"].notna().any():
		lead = d.loc[d["P"].idxmin()]
	else:
		z = (d["BETA"] / d["SE"]).abs()
		lead = d.loc[z.idxmax()]
	chrom = re.sub(r"^chr", "", str(lead["CHR"]), flags=re.IGNORECASE)
	pos = int(lead["POS"])
	w = window_kb * 1000
	return chrom, max(1, pos - w), pos + w


def to_signal(d, signal, typ, sdY=None):
	"""Standalone BF conversion; variant identities must already be reconciled."""
	if d.empty:
		raise ValueError("Empty signal")
	if typ == "quant" and (sdY is None or not np.isfinite(sdY) or sdY <= 0):
		raise ValueError("Quantitative BF requires documented sdY")
	if "variant" not in d:
		raise ValueError("A shared prepared variant key is required; allele sorting is not normalization")
	v = d.SE.to_numpy(float)**2
	w = (.15 * sdY)**2 if typ == "quant" else .2**2
	lbf = .5 * (np.log(v) - np.log(v+w) + w/(v+w)*(d.BETA.to_numpy(float)/np.sqrt(v))**2)
	if not np.isfinite(lbf).all():
		raise ValueError("Nonfinite BF")
	z = pd.DataFrame({"variant": d.variant.astype(str), "lbf": lbf}).drop_duplicates()
	if z.variant.duplicated().any():
		raise ValueError("Conflicting duplicate variants")
	wide = pd.DataFrame([z.lbf.to_numpy()], columns=z.variant)
	i = int(np.argmax(z.lbf))
	return wide, str(z.variant.iloc[i]), float(z.lbf.iloc[i])


def reference_posterior(x, y, p1=1e-4, p2=1e-4, p12=1e-5):
	"""Stable H3 without subtracting nearly equal exponentials."""
	from scipy.special import logsumexp
	x, y = np.asarray(x,dtype=float), np.asarray(y,dtype=float)
	if x.shape != y.shape or x.ndim != 1 or not len(x) or not np.isfinite([x,y]).all():
		raise ValueError("Invalid BF vectors")
	pre = np.r_[-np.inf, np.logaddexp.accumulate(y)]
	suf = np.r_[np.logaddexp.accumulate(y[::-1])[::-1], -np.inf]
	off = np.logaddexp(pre[:-1], suf[1:])
	l = np.array([0,np.log(p1)+logsumexp(x),np.log(p2)+logsumexp(y),np.log(p1)+np.log(p2)+logsumexp(x+off),np.log(p12)+logsumexp(x+y)])
	return np.exp(l-logsumexp(l))


def prepare_signals(manifest, cad_gwas, outdir, window_kb, outcome_type="cc"):
	"""Use the exact ordered, harmonized BF vectors prepared once by C3 R."""
	outdir = Path(outdir)
	qdir, ydir = outdir/"qtl_signals", outdir/"cad_signals"
	qdir.mkdir(parents=True, exist_ok=True); ydir.mkdir(parents=True, exist_ok=True)
	for f in [*qdir.glob("*.feather"), *ydir.glob("*.feather")]: f.unlink()
	rows = pd.read_csv(manifest, sep="\t").fillna("")
	qs, ys, pairs, statuses = [], [], [], []
	from scipy.special import logsumexp
	for i, row in rows.iterrows():
		base = {k: str(row.get(k,"")) for k in ["omics","trait","region","snp_hash","BF_model","prior_config"]}
		try:
			path = Path(str(row.get("pair_file","")))
			if not path.is_file(): raise ValueError("Shared CPU pair unavailable: " + str(row.get("input_status","not prepared")))
			d = pd.read_csv(path,sep="\t",dtype={"variant":str})
			if not {"variant","lbf1","lbf2"} <= set(d): raise ValueError("Invalid prepared pair columns")
			if d.empty or d.variant.isna().any() or d.variant.duplicated().any(): raise ValueError("Invalid prepared variant identities")
			if not np.isfinite(d[["lbf1","lbf2"]].to_numpy()).all(): raise ValueError("Nonfinite prepared BF")
			key = hashlib.sha256("\n".join(d.variant).encode()).hexdigest()
			if key != base["snp_hash"]: raise ValueError("Shared SNP hash mismatch")
			chrom, start, end = re.match(r"chr([^:]+):(\d+)-(\d+)$",base["region"]).groups()
			qx = safe_name(f"QTL__{base['omics']}__{base['trait']}__{i}")
			yx = safe_name(f"OUTCOME__{base['trait']}__{i}")
			for col, sig, directory, summaries in [("lbf1",qx,qdir,qs),("lbf2",yx,ydir,ys)]:
				# Explicit SNP intersection per pair: no unmeasured SNPs become zero evidence.
				feather.write_feather(pd.DataFrame([d[col].to_numpy()],columns=d.variant),directory/(sig+".feather"))
				j=int(np.argmax(d[col]))
				summaries.append(dict(signal=sig,chromosome=chrom,location_min=int(start),location_max=int(end),signal_strength=float(d[col].iloc[j]),lead_variant=str(d.variant.iloc[j])))
			pp=reference_posterior(d.lbf1.to_numpy(),d.lbf2.to_numpy(),p12=float(os.getenv("GPU_COLOC_P12","1e-5")))
			pairs.append(dict(**base,qtl_signal=qx,cad_signal=yx,n_snps=len(d),pair_file=str(path),**{f"reference_PP.H{j}":pp[j] for j in range(5)}))
			statuses.append(dict(**base,status="ok",message=f"{len(d)} identical shared CPU/GPU variants"))
		except Exception as exc:
			statuses.append(dict(**base,status="input_invalid",message=str(exc)))
	pd.DataFrame(qs).to_csv(outdir/"qtl_summary.tsv",sep="\t",index=False)
	pd.DataFrame(ys).to_csv(outdir/"cad_summary.tsv",sep="\t",index=False)
	pd.DataFrame(pairs,columns=list(pairs[0]) if pairs else ["qtl_signal","cad_signal"]).to_csv(outdir/"expected_pairs.tsv",sep="\t",index=False)
	pd.DataFrame(statuses).to_csv(outdir/"signal_preparation_status.tsv",sep="\t",index=False)
	if not pairs: raise SystemExit("No eligible shared CPU/GPU pair; see signal_preparation_status.tsv")
	print(f"Prepared {len(pairs)} exact common-variant pairs",flush=True)


def filter_results(raw_results, expected_pairs, output_file, h4_threshold=0.70):
	"""Keep only the QTL/outcome pair requested by each C3 manifest row.

	gpu-coloc intentionally tests every overlapping signal pair. C3 creates a
	disease signal for each requested region, so overlapping regions can also
	create scientifically unintended cross-pairs. Retain the exact manifest
	pairs and explicitly record pairs for which gpu-coloc returned no result.
	"""
	try:
		h4_threshold = float(h4_threshold)
	except (TypeError, ValueError):
		raise SystemExit("GPU-coloc H4 threshold must be numeric")
	if not 0 <= h4_threshold <= 1:
		raise SystemExit("GPU-coloc H4 threshold must be between 0 and 1")
	expected = pd.read_csv(expected_pairs, sep="\t", dtype=str)
	if expected.empty:
		raise SystemExit("No expected GPU-coloc pairs were recorded: " + expected_pairs)
	raw = pd.DataFrame()
	if os.path.exists(raw_results) and os.path.getsize(raw_results) > 0:
		try:
			raw = pd.read_csv(raw_results, sep="\t")
		except pd.errors.EmptyDataError:
			raw = pd.DataFrame()
	if not raw.empty and {"signal1", "signal2"}.issubset(raw.columns):
		pp_col = "PP.H4" if "PP.H4" in raw.columns else None
		if pp_col:
			raw[pp_col] = pd.to_numeric(raw[pp_col], errors="coerce")
			raw = raw.drop_duplicates()
		conflict = raw.duplicated(["signal1","signal2"],keep=False)
		bad = set(map(tuple,raw.loc[conflict,["signal1","signal2"]].to_numpy()))
		raw = raw.drop_duplicates(["signal1","signal2"])
		raw["duplicate_conflict"] = [tuple(v) in bad for v in raw[["signal1","signal2"]].to_numpy()]
		if pp_col: raw.loc[raw.duplicate_conflict,pp_col] = np.nan
		out = expected.merge(
			raw,
			how="left",
			left_on=["qtl_signal", "cad_signal"],
			right_on=["signal1", "signal2"],
			indicator="_gpu_merge",
		)
	else:
		out = expected.copy()
		out["signal1"] = out["qtl_signal"]
		out["signal2"] = out["cad_signal"]
		out["PP.H4"] = np.nan
		out["_gpu_merge"] = "left_only"
	if "signal1" in out:
		out["signal1"] = out["signal1"].fillna(out["qtl_signal"])
	if "signal2" in out:
		out["signal2"] = out["signal2"].fillna(out["cad_signal"])
	if "PP.H4" not in out:
		out["PP.H4"] = np.nan
	pp = pd.to_numeric(out["PP.H4"], errors="coerce")
	returned = out["_gpu_merge"].astype(str) == "both"
	finite = np.isfinite(pp)
	out["GPU_status"] = np.select(
		[~returned, returned & finite],
		["no_gpu_result", "ok"],
		default="gpu_result_na",
	)
	threshold_label = f"{h4_threshold:g}"
	out["GPU_H4_class"] = np.select(
		[~returned, returned & ~finite, pp >= h4_threshold],
		["no_gpu_result", "no_finite_pp_h4", f"H4>={threshold_label}"],
		default=f"H4<{threshold_label}",
	)
	out["GPU_H4_threshold"] = h4_threshold
	if "duplicate_conflict" in out:
		out.loc[out.duplicate_conflict.fillna(False).astype(bool),"GPU_status"] = "validation_failed_duplicate"
	if "reference_PP.H4" in out:
		ref = pd.to_numeric(out["reference_PP.H4"],errors="coerce")
		out["absolute_difference_H4"] = abs(pp-ref)
		out["parity_tolerance"] = float(os.getenv("GPU_COLOC_PARITY_TOL","1e-4"))
		out["parity_status"] = np.where(finite,np.where(out.absolute_difference_H4 <= out.parity_tolerance,"pass","numerical_disagreement"),"native_result_unavailable")
		out["threshold_disagreement"] = np.where(finite,(pp>=h4_threshold)!=(ref>=h4_threshold),None)
		out["reference_role"] = "double precision diagnostic; native result retained"
	try: out["backend_version"] = importlib.metadata.version("gpu-coloc")
	except importlib.metadata.PackageNotFoundError: out["backend_version"] = "unavailable"
	out["native_dtype"] = "float32"
	out = out.drop(columns=["_gpu_merge"])
	out.to_csv(output_file, sep="\t", index=False)


def main():
	if len(sys.argv) in (5, 6) and sys.argv[1] == "filter-results":
		threshold = sys.argv[5] if len(sys.argv) == 6 else 0.70
		filter_results(sys.argv[2], sys.argv[3], sys.argv[4], threshold)
		return
	if len(sys.argv) not in (5, 6):
		raise SystemExit(
			"Usage: c3.coloc_GPU.py <manifest.tsv> <outcome_gwas> <outdir> <window_kb> [cc|quant]\n"
			"   or: c3.coloc_GPU.py filter-results <raw.tsv> <expected_pairs.tsv> <output.tsv> [H4_threshold]"
		)
	outcome_type = sys.argv[5].lower() if len(sys.argv) == 6 else "cc"
	if outcome_type not in {"cc", "quant"}:
		raise SystemExit("outcome_type must be cc or quant")
	prepare_signals(sys.argv[1], sys.argv[2], sys.argv[3], int(float(sys.argv[4])), outcome_type)


if __name__ == "__main__":
	main()
