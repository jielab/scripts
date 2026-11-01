#!/usr/bin/env python3
"""Normalize current method results into SQLite and export introgression density."""

from __future__ import annotations
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
sys.dont_write_bytecode = True

if "gu_0_common" not in sys.modules:
	_spec = spec_from_file_location("gu_0_common", Path(__file__).with_name("0.common.py"))
	_common = module_from_spec(_spec)
	sys.modules[_spec.name] = _common
	try:
		_spec.loader.exec_module(_common)
	except BaseException:
		sys.modules.pop(_spec.name, None)
		raise
from gu_0_common import load_module, reference_role, evidence_sql, x_nonpar_intervals, gu_dataset_from_path, gu_method_runs, gu_method_directory


# 🚩 normalize_results
"""Normalize GU methods into a common, UKB-friendly segment database.

The database is the primary artifact for large cohorts. Exact raw tract boundaries are
kept, while a cohort-level canonical segment_code is assigned by reciprocal-overlap
clustering. This gives each participant 0-N stable-ish tract codes analogous to an
ICD event table, while preserving the original calls for sensitivity analyses.
"""

import argparse, bisect, gzip, hashlib, json, math, os, re, shutil, sqlite3, tempfile
from pathlib import Path
import pandas as pd

load_module("0.common.py")
from gu_0_common import resolve_tsv_path, DEFAULT_ARCHAIC_REFERENCES, SUPPORTED_ARCHAIC_REFERENCES

SEG_COLS = [
	"dataset_id",
	"sample_id",
	"method",
	"source",
	"source_class",
	"reference_role",
	"chr",
	"start",
	"end",
	"length_bp",
	"haplotype",
	"method_haplotype_index",
	"score",
	"posterior",
	"trait",
	"locus_id",
	"genome_build",
	"batch_id",
	"raw_file",
]
POPS_1KG = {
	"ACB",
	"ASW",
	"ESN",
	"GWD",
	"LWK",
	"MSL",
	"YRI",
	"CLM",
	"MXL",
	"PEL",
	"PUR",
	"CDX",
	"CHB",
	"CHS",
	"JPT",
	"KHV",
	"CEU",
	"FIN",
	"GBR",
	"IBS",
	"TSI",
	"BEB",
	"GIH",
	"ITU",
	"PJL",
	"STU",
}


def chrom_norm(x):
	s = str(x)
	s = re.sub(r"^chr", "", s, flags=re.I)
	s = re.sub(r"_part\d+$", "", s, flags=re.I)
	return "X" if s == "23" else s


def source_class(x):
	s = str(x).lower()
	if any(k in s for k in ("altai", "vindija", "chagyr", "neander")):
		return "Neanderthal"
	if "denis" in s:
		return "Denisovan"
	if "mosaic" in s:
		return "Mosaic"
	if "trace" in s or "ghost" in s or "unknown" in s:
		return "Ghost/Unknown"
	return "Archaic/Other"


def empty_segments():
	return pd.DataFrame(columns=SEG_COLS)


def build_norm(x, fallback):
	if pd.isna(x) or str(x).strip() in ("", "nan", "None"):
		x = fallback
	s = str(x).strip()
	if re.fullmatch(r"b?3[78]", s, re.I):
		return f"GRCh{s[-2:]}"
	if re.fullmatch(r"GRCh3[78]", s, re.I):
		return f"GRCh{s[-2:]}"
	return s


def dataset_from_path(path, method, default="1kg"):
	return gu_dataset_from_path(path, method, default)


def finish(df, method, build, raw, batch="", dataset=None):
	if df.empty:
		return empty_segments()
	for c in SEG_COLS:
		if c not in df:
			df[c] = pd.NA
	df["dataset_id"] = dataset or dataset_from_path(raw, method)
	df["method"] = method
	df["genome_build"] = df["genome_build"].map(lambda x: build_norm(x, build))
	df["raw_file"] = str(raw)
	df["batch_id"] = batch
	df["chr"] = df["chr"].map(chrom_norm)
	df["start"] = pd.to_numeric(df["start"], errors="coerce")
	df["end"] = pd.to_numeric(df["end"], errors="coerce")
	df = df[df.start.notna() & df.end.notna() & (df.end > df.start)].copy()
	df["start"] = df.start.astype("int64")
	df["end"] = df.end.astype("int64")
	df["length_bp"] = (df.end - df.start).astype("int64")
	df["source"] = df["source"].fillna("unknown").astype(str)
	df["source_class"] = df["source"].map(source_class)
	df["reference_role"] = df["source"].map(reference_role)
	return df[SEG_COLS]


def is_transient_output(path):
	"""Return True for atomic-work directories such as .chr22.run.<pid>."""
	return any(part.startswith(".") and (".run." in part or ".old." in part) for part in path.parts)


def find_first(root, patterns):
	out = []
	for p in patterns:
		out.extend(root.glob("/".join(gu_method_directory(x) if x in ("ibdmix", "trace", "phyml") else x for x in p.split("/"))))
	out = [x for x in out if not is_transient_output(x) and x.is_file() and x.stat().st_size > 0]
	return sorted(set(out))


# Results are read from the persistent method directories. Shiny's relative
# artifact paths resolve through a link to the same native PhyML results.
def package_core_results(analysis_root, output_dir):
	output_dir.mkdir(parents = True, exist_ok = True)
	link = output_dir / 'phyml'
	target = (analysis_root / gu_method_directory("phyml")).resolve()
	if link.is_symlink():
		if link.resolve() != target:
			link.unlink()
	elif link.exists():
		raise ValueError('Expected a shared PhyML result symlink, but a directory already exists: ' + str(link))
	if not link.is_symlink():
		link.symlink_to(target, target_is_directory = True)
	return {'phyml': 'shared method results'}


def portable_paths(df, package_root, columns):
	"""Store paths inside the share package relative to normalize/, leaving external provenance intact."""
	root = package_root.resolve()
	for column in columns:
		if column not in df:
			continue

		def convert(value):
			if pd.isna(value) or not str(value):
				return value
			try:
				return Path(str(value)).resolve().relative_to(root).as_posix()
			except ValueError:
				return str(value)

		# Millions of rows commonly share one provenance path. Resolve each
		# distinct path once, especially on mounted Windows filesystems.
		mapped = {value: convert(value) for value in df[column].dropna().unique()}
		df[column] = df[column].map(lambda value: mapped.get(value, value))
	return df


def phyml_artifact_in_run(value, run_dir):
	"""Resolve the current run's explicitly recorded artifact path."""
	if value is None or pd.isna(value) or not str(value):
		return ""
	path = Path(str(value))
	return str(path if path.is_absolute() else run_dir / path)


def active_ibdmix_run(run):
	"""MAC2/3 are separate X sensitivity runs, excluded from the primary summary."""
	return not re.search(r"[.]nonpar-v2[.]mac[23]$", run.name)


def require_ibdmix_metadata(meta, path):
	if not meta.get("pipeline_version"):
		raise ValueError(f"IBDmix result lacks pipeline_version: {path}")
	if meta.get("coordinate_system") != "0-based-half-open":
		raise ValueError(f"IBDmix result requires coordinate_system=0-based-half-open: {path}")
	refs = meta.get("refs", "").split()
	if not refs or len(refs) != len(set(refs)):
		raise ValueError(f"IBDmix result requires a unique nonempty refs record: {path}")
	summary_references(refs)


def iter_ibdmix(root, build, completed_runs=None, filter_output=None, daf_mask=None, audits=None):
	load_module("ibdmix.py")
	from gu_ibdmix import prepare

	files = find_first(
		root,
		[
			"ibdmix/**/final/all_archaic_refs*.segments.tsv.gz",
		],
	)
	for f in files:
		if not active_ibdmix_run(f.parent.parent):
			continue
		if completed_runs is not None and f.parent.parent not in completed_runs:
			continue  # jobs finishing after this final snapshot belong to the next one
		meta = read_run_meta(f.parent.parent / "run.meta.tsv")
		if not (f.parent.parent / ".complete").exists():
			continue
		require_ibdmix_metadata(meta, f.parent.parent)
		batch = meta.get("batch_id", "")
		ukb_candidate = dataset_from_path(f, "ibdmix").startswith("ukb-")
		if ukb_candidate:
			# Imputed/site-restricted UKB calls have a different ascertainment
			# and no AFR5 background panel. Keep their native candidate calls;
			# reference/background filtering requires the matching 1KG cohort.
			if not (f.parent.parent / "inputs/ukb.preparation.qc.json").is_file():
				raise ValueError(f"UKB preparation QC is missing: {f.parent.parent}")
			filtered = f
			chrom_match = re.search(r"(?:^|[/_.-])chr([0-9]+|X)(?:$|[/_.-])", str(f), re.I)
			audit = {"status": "ukb_exploratory_reference_affinity", "daf": {"status": "not_applied_ukb_pilot"},
				"chrom": chrom_norm(chrom_match.group(1)) if chrom_match else "",
				"reference_affinity_only": True, "strict_evidence_eligible": False}
		else:
			filtered, audit = prepare(f, filter_output or root / "final/ibdmix_filters", meta, daf_mask)
		if audits is not None:
			audits.append(
				(
					dataset_from_path(f, "ibdmix"),
					build_norm(meta.get("genome_build"), build),
					audit.get("chrom", "X" if "chrX" in str(f) else ""),
					str(f),
					str(filtered),
					audit["status"],
					audit["daf"]["status"],
					json.dumps(audit),
				)
			)
		print(f"[GU FINAL] IBDmix {f.parent.parent.name}: {audit['status']}; DAF={audit['daf']['status']}", flush=True)
		read_options = {"dtype": {"ID": str, "sample_id": str}}
		for d in pd.read_csv(filtered, sep="\t", compression="infer", chunksize=200000, **read_options):
			ren = {}
			for a, b in [
				("ID", "sample_id"),
				("chrom", "chr"),
				("length", "length_bp"),
				("slod", "score"),
				("anc", "source"),
			]:
				if a in d:
					ren[a] = b
			d = d.rename(columns=ren)
			if not {"sample_id", "chr", "start", "end", "source", "score"} <= set(d.columns):
				raise ValueError(f"Invalid IBDmix final segment columns: {filtered}")
			d["haplotype"] = pd.NA
			d["posterior"] = pd.NA
			d["trait"] = pd.NA
			if "locus_id" not in d:
				d["locus_id"] = pd.NA
			if "source" not in d:
				d["source"] = "IBDmix_archaic"
			yield finish(d, "ibdmix", build, filtered, batch=batch, dataset=dataset_from_path(f, "ibdmix"))


def iter_trace(root, build):
	files = find_first(
		root,
		[
			"trace/**/final/trace_haplotype_segments.tsv.gz",
		],
	)
	for f in files:
		batch = read_run_meta(find_parent_file(f, "run.meta.tsv")).get("batch_id", "")
		for d in pd.read_csv(f, sep="\t", compression="infer", chunksize=200000, dtype={"sample": str}):
			d = d.rename(columns={"sample": "sample_id", "chrom": "chr", "mean_posterior": "posterior"})
			d["source"] = "TRACE_ghost_or_unknown"
			d["score"] = pd.NA
			d["trait"] = pd.NA
			if "locus_id" not in d:
				d["locus_id"] = pd.NA
			yield finish(d, "trace", build, f, batch=batch, dataset=dataset_from_path(f, "trace"))


def find_parent_file(f, name):
	for parent in [f.parent, *f.parents]:
		candidate = parent / name
		if candidate.is_file() and candidate.stat().st_size:
			return candidate
	return None




def read_run_meta(path):
	values = {}
	if not path or not Path(path).is_file():
		return values
	try:
		with Path(path).open(errors="replace") as handle:
			for line in handle:
				fields = line.rstrip("\n").split("\t", 1)
				if len(fields) == 2 and fields[0] not in values:
					values[fields[0]] = fields[1]
	except OSError:
		pass
	return values


def discover_method_runs(root, fallback_build):
	"""Record completed chromosome units, including valid zero-call runs.

	Segment presence is not a completion signal: a successful caller can emit
	only a header.  These rows let Shiny distinguish a true negative result
	from a method that was never run (or is unsupported on chrX).
	"""
	rows = {}

	def add(method, path, build, chrom, status="complete"):
		chrom = chrom_norm(chrom)
		build = build_norm(build, fallback_build)
		if not chrom:
			return
		eligible = 1
		note = "completed native analysis"
		if method == "ibdmix":
			meta = read_run_meta(path)
			require_ibdmix_metadata(meta, path)
			note = "audited IBDmix " + meta["pipeline_version"] + "; profile=" + meta.get("profile", "unknown")
		if chrom == "X" and method == "ibdmix":
			eligible = 0
			note = "experimental male-X pseudo-diploid adaptation; excluded from GU evidence score"
		elif method == "trace":
			note = "haplotype-node call; validity depends on the full-chromosome ARG backend and calibration"
		if method == "ibdmix" and dataset_from_path(path, method).startswith("ukb-"):
			eligible = 0
			qualification = "experimental UKB reference-affinity candidates; imputation and site ascertainment require calibration; no AFR5 background subtraction"
			if chrom == "X":
				qualification += "; male non-PAR X pseudo-diploid LOD model is uncalibrated"
			note = "audited IBDmix " + meta["pipeline_version"] + "; " + qualification + "; profile=" + meta.get("profile", "unknown")
		key = (dataset_from_path(path, method), build, chrom, method)
		if key in rows and rows[key][4] == "complete" and status != "complete":
			return
		rows[key] = (*key, status, eligible, note, str(path))

	for marker in sorted(root.glob(gu_method_directory("ibdmix")+"/**/.complete")):
		if not active_ibdmix_run(marker.parent):
			continue
		meta = marker.parent / "run.meta.tsv"
		values = read_run_meta(meta)
		text = meta.read_text(errors="replace") if meta.is_file() else ""
		build = values.get("genome_build", fallback_build)
		chroms = set(re.findall(r"(?:^|[/\t_.-])chr([0-9]+|X)(?=$|[/\t_.-])", text, flags=re.I | re.M))
		if not chroms:
			chroms = set(re.findall(r"chr([0-9]+|X)", str(marker.parent), flags=re.I))
		for chrom in chroms:
			add("ibdmix", meta if meta.is_file() else marker, build, chrom)

	trace_files = find_first(
		root, ["trace/**/final/trace_haplotype_segments.tsv.gz"]
	)
	for output in trace_files:
		meta = find_parent_file(output, "run.meta.tsv")
		values = read_run_meta(meta)
		for chrom in re.split(r"[, ]+", values.get("chromosomes", "").strip()):
			if chrom:
				add("trace", meta or output, values.get("build", fallback_build), chrom)

	return list(rows.values())










def sample_population_table(path, dataset_id="1kg", analysis_root=None):
	cols = ["dataset_id", "sample_id", "population", "super_population", "raw_file"]
	if dataset_id != "1kg":
		# A UKB population label is not a 1KG population. Keep its own dataset
		# key and preserve IDs as strings; missing metadata must not become EUR.
		if not path or not Path(path).is_file():
			panels = []
			if analysis_root and re.fullmatch(r"[A-Za-z0-9_.-]+", dataset_id):
				for method in ("ibdmix", "trace", "phyml"):
					for run in gu_method_runs(analysis_root, method):
						if dataset_from_path(run / 'run.meta.tsv', method) == dataset_id:
							panels.extend(p for name in ('ukb.samples.tsv', 'sample.panel.tsv') if (p := run / 'inputs' / name).is_file())
			if not panels:
				raise ValueError(f"Sample metadata is required for dataset {dataset_id}; restore its native archives or provide --sample-panel")
			out = pd.concat([sample_population_table(p, dataset_id) for p in panels], ignore_index=True)
			labels = out[["sample_id", "population", "super_population"]].drop_duplicates()
			if labels.sample_id.duplicated().any():
				raise ValueError("Archived UKB sample panels disagree on population labels")
			return out.drop_duplicates("sample_id")[cols]
		d = pd.read_csv(path, sep=r"\s+", dtype=str, keep_default_na=False)
		low = {str(c).lower().lstrip("#"): c for c in d.columns}
		sample = next((low[x] for x in ("sample", "sample_id", "iid") if x in low), None)
		pop = next((low[x] for x in ("population", "pop") if x in low), None)
		super_pop = next((low[x] for x in ("super_population", "super_pop", "superpopulation") if x in low), None)
		if sample is None or pop is None or super_pop is None or d.empty:
			raise ValueError("Sample panel requires nonempty sample, pop and super_pop columns")
		out = pd.DataFrame({"dataset_id": dataset_id, "sample_id": d[sample], "population": d[pop].str.upper(),
			"super_population": d[super_pop].str.upper(), "raw_file": str(path)})
		if out.sample_id.duplicated().any() or out[["sample_id", "population", "super_population"]].isin(["", ".", "NA", "NAN"]).any().any():
			raise ValueError("Sample metadata contains duplicate IDs or missing identity/population labels")
		return out[cols]
	if not path or not Path(path).is_file():
		return pd.DataFrame(columns=cols)
	try:
		d = pd.read_csv(path, sep=r"\s+", dtype=str, comment="#")
	except Exception:
		return pd.DataFrame(columns=cols)
	low = {str(c).lower().lstrip("#"): c for c in d.columns}
	sample = next((low[x] for x in ("sample", "sample_id", "iid") if x in low), None)
	pop = next((low[x] for x in ("population", "pop") if x in low), None)
	super_pop = next((low[x] for x in ("super_population", "super_pop", "superpopulation") if x in low), None)
	if sample is None or pop is None:
		return pd.DataFrame(columns=cols)
	out = pd.DataFrame(
		{"dataset_id": "1kg", "sample_id": d[sample].astype(str), "population": d[pop].astype(str).str.upper()}
	)
	out["super_population"] = d[super_pop].astype(str).str.upper() if super_pop else pd.NA
	out["raw_file"] = str(path)
	return out[out.sample_id.ne("") & out.population.isin(POPS_1KG)].drop_duplicates("sample_id")[cols]


def truth_int(value):
	if pd.isna(value):
		return 0
	return int(str(value).strip().upper() in {"1", "TRUE", "T", "YES", "Y", "PASS", "DIRECT_MATCH_PASS"})


def iter_phyml_carriers(root, fallback_build):
	"""Normalize every tested PhyML copy and keep tree and pairwise flags separate."""
	for hp in sorted((root / gu_method_directory("phyml")).glob("**/final/haplotypes.tsv")):
		if not hp.is_file() or not hp.stat().st_size:
			continue
		haps = pd.read_csv(hp, sep="\t", dtype=str)
		if not {"locus_id", "hap_id", "copies"} <= set(haps.columns):
			raise ValueError(f"Invalid current PhyML haplotype columns: {hp}")
		dataset = dataset_from_path(hp, "phyml")
		# The copies column records every tested sample/haplotype and is the
		# prevalence denominator.  A header-only table is a valid result for a
		# locus with no selected haplotypes (for example, low LD information).
		if haps.empty:
			continue
		if not haps["copies"].notna().any():
			raise ValueError(f"PhyML haplotypes table has no populated copies values: {hp}")
		tree_tips = {}
		tree_pass = {}
		# Keep carriers from every passing candidate-focused lineage tree.
		# trees.tsv is only the representative topology for the locus display.
		hap_lineages = {}
		evidence_trees = hp.parent / "evidence_trees.tsv"
		if not evidence_trees.is_file():
			raise ValueError(f"Current PhyML lineage results are missing: {evidence_trees}")
		if evidence_trees.stat().st_size:
			ed = pd.read_csv(evidence_trees, sep="\t", dtype=str, keep_default_na=False)
			for lid in ed["locus_id"].astype(str).unique():
				tree_tips[lid] = set()
				tree_pass[lid] = 0
			for _, row in ed.iterrows():
				if not truth_int(row.get("candidate_clade_pass", 0)):
					continue
				lid = str(row.get("locus_id", ""))
				lin = str(row.get("expected_lineage", ""))
				tips = {x for x in str(row.get("candidate_tips_in_clade", "")).split(",") if x.startswith("H")}
				tree_tips.setdefault(lid, set()).update(tips)
				tree_pass[lid] = 1
				for hid in tips:
					hap_lineages.setdefault((lid, hid), set()).add(lin)
		buffer = []
		for _, hap in haps.iterrows():
			lid = str(hap["locus_id"])
			hap_id = str(hap["hap_id"])
			in_supported_clade = int(bool(tree_pass.get(lid, 0) and hap_id in tree_tips.get(lid, set())))
			copies = "" if pd.isna(hap.get("copies", pd.NA)) else str(hap.get("copies", ""))
			for copy in copies.split(";"):
				if not copy or ":" not in copy:
					continue
				sample, haplotype = copy.rsplit(":", 1)
				if not sample:
					continue
				buffer.append(
					{
						"dataset_id": dataset,
						"genome_build": build_norm(hap.get("genome_build", None), fallback_build),
						"locus_id": lid,
						"sample_id": sample,
						"haplotype": haplotype,
						"hap_id": hap_id,
						"best_archaic": hap.get("best_archaic", pd.NA),
						"best_lineage": hap.get("best_lineage", pd.NA),
						"candidate_lineage": ",".join(sorted(hap_lineages.get((lid, hap_id), set()))),
						"n_compared": pd.to_numeric(hap.get("n_compared", pd.NA), errors="coerce"),
						"n_match": pd.to_numeric(hap.get("n_match", pd.NA), errors="coerce"),
						"prop_match": pd.to_numeric(hap.get("prop_match", pd.NA), errors="coerce"),
						"direct_match_pass": truth_int(hap.get("direct_match_pass", pd.NA)),
						"tree_candidate_pass": in_supported_clade,
						"raw_file": str(hp),
					}
				)
				if len(buffer) >= 100000:
					yield pd.DataFrame(buffer)
					buffer = []
		if buffer:
			yield pd.DataFrame(buffer)


def loci_tables(root, fallback_build):
	def scalar_int(value):
		try:
			if value is None or pd.isna(value):
				return None
			return int(float(value))
		except (TypeError, ValueError):
			return None

	def scalar_float(value):
		try:
			if value is None or pd.isna(value):
				return None
			return float(value)
		except (TypeError, ValueError):
			return None

	def scalar_text(value):
		return None if value is None or pd.isna(value) or not str(value).strip() else str(value)

	seen = set()
	for base in (root / gu_method_directory("phyml"),):
		for pat in ["**/final/loci.tsv"]:
			for f in base.glob(pat):
				if not f.is_file() or f.stat().st_size == 0:
					continue
				d = pd.read_csv(f, sep="\t")
				required = {"locus_id", "chr", "core_start", "core_end", "analysis_start", "analysis_end", "genome_build", "lineage", "selection_method"}
				if not required <= set(d.columns):
					raise ValueError(f"Invalid current PhyML locus columns: {f}")
				if d.empty:
					continue
				trees = {}
				tree_file = f.parent / "trees.tsv"
				if tree_file.is_file() and tree_file.stat().st_size:
					try:
						td = pd.read_csv(tree_file, sep="\t", dtype=str, keep_default_na=False)
						trees = {str(row.get("locus_id", "")): row for _, row in td.iterrows()}
					except Exception:
						pass
				dataset = dataset_from_path(f, "phyml")
				for _, r in d.iterrows():
					chrom = r["chr"]
					input_st, input_en = scalar_int(r["core_start"]), scalar_int(r["core_end"])
					ast, aen = scalar_int(r["analysis_start"]), scalar_int(r["analysis_end"])
					if any(value is None for value in (input_st, input_en, ast, aen)) or input_st < 0 or input_en <= input_st or ast < 0 or aen <= ast:
						raise ValueError(f"Invalid PhyML interval coordinates: {f}")
					st, en = input_st, input_en
					method = "phyml"
					lid = r["locus_id"]
					genome_build = build_norm(r["genome_build"], fallback_build)
					key = (dataset, method, genome_build, str(lid), chrom, st, en)
					if key in seen:
						continue
					seen.add(key)
					tree = trees.get(str(lid), {})
					run_dir = f.parent.parent
					candidate_lineage = scalar_text(tree.get("candidate_lineage", None))
					source = candidate_lineage or scalar_text(r["lineage"]) or "unknown"
					audit = {k: None if pd.isna(v) else v for k, v in r.items()}
					if r.get("selection_method") == "gwas_lead_ld_core":
						audit.update(load_module("phyml.py").saved_locus_audit(f, audit))
					yield {
						'n_panel': scalar_int(audit.get('n_panel')),
						'n_panel_copies': scalar_int(audit.get('n_panel_copies')),
						'n_lead_called': scalar_int(audit.get('n_lead_called')),
						'n_lead_called_copies': scalar_int(audit.get('n_lead_called_copies')),
						'n_recurrent_copies': scalar_int(audit.get('n_recurrent_copies')),
						'n_singleton_copies': scalar_int(audit.get('n_singleton_copies')),
						'input_status': scalar_text(audit.get('input_status')),
						'ld_status': scalar_text(audit.get('ld_status')),
						'alignment_status': scalar_text(audit.get('alignment_status')),
						'tree_run_status': scalar_text(audit.get('tree_run_status')),
						'topology_status': scalar_text(audit.get('topology_status')),
						'sequence_qc': scalar_text(audit.get('sequence_qc')),
						'lead_retained_for_tree': scalar_int(audit.get('lead_retained_for_tree')),
						'high_ld_sites_retained_for_tree': scalar_int(audit.get('high_ld_sites_retained_for_tree')),
						'unique_archaic_sequences': scalar_int(audit.get('unique_archaic_sequences')),
						'high_ld_max_gap_bp': scalar_int(audit.get('high_ld_max_gap_bp')),
						'site_missingness_audit': scalar_text(audit.get('site_missingness_audit')),
						"dataset_id": dataset,
						"method": method,
						"trait": None,
						"locus_id": lid,
						"chr": chrom_norm(chrom),
						"start": st,
						"end": en,
						"input_start": input_st,
						"input_end": input_en,
						"analysis_start": ast,
						"analysis_end": aen,
						"flank_bp": max(0, max(input_st - ast, aen - input_en)),
						"anchor_pos": scalar_int(r.get("lead_pos", None)),
						"selection_method": scalar_text(r.get("selection_method", None)),
						"ld_r2_threshold": scalar_float(r.get("ld_r2_threshold", None)),
						"n_search_sites": scalar_int(r.get("n_search_sites", None)),
						"n_ld_sites": scalar_int(r.get("n_ld_sites", None)),
						"n_ancestral_sites": scalar_int(r.get("n_ancestral_sites", None)),
						"n_sites": scalar_int(r.get("n_sites", None)),
						"best_archaic": scalar_text(r.get("best_archaic", None)),
						"best_lineage": scalar_text(r.get("best_lineage", None)),
						"n_compared": scalar_int(r.get("n_compared", None)),
						"n_match": scalar_int(r.get("n_match", None)),
						"prop_match": scalar_float(r.get("prop_match", None)),
						"source": str(source),
						"source_class": source_class(source),
						"status": r.get("status", None),
						"direct_match_pass": r.get("direct_match_pass", None),
						"n_carriers": None,
						"genome_build": genome_build,
						"tree_status": tree.get("tree_status", "not_run"),
						"tree_file": phyml_artifact_in_run(tree.get("tree_file", ""), run_dir),
						"stats_file": phyml_artifact_in_run(tree.get("stats_file", ""), run_dir),
						"plot_file": phyml_artifact_in_run(tree.get("plot_file", ""), run_dir),
						"n_bootstrap_nodes": tree.get("n_bootstrap_nodes", None),
						"bootstrap_min": tree.get("bootstrap_min", None),
						"bootstrap_median": tree.get("bootstrap_median", None),
						"bootstrap_max": tree.get("bootstrap_max", None),
						"candidate_lineage": candidate_lineage,
						"tree_has_ancestral_outgroup": truth_int(tree.get("tree_has_ancestral_outgroup", None)),
						"candidate_clade_pass": truth_int(tree.get("candidate_clade_pass", None)),
						"candidate_clade_rule": scalar_text(tree.get("candidate_clade_rule", None)),
						"candidate_clade_bootstrap": scalar_float(tree.get("candidate_clade_bootstrap", None)),
						"candidate_clade_node": scalar_text(tree.get("candidate_clade_node", None)),
						"candidate_clade_side": scalar_text(tree.get("candidate_clade_side", None)),
						"candidate_clade_n_tips": scalar_int(tree.get("candidate_clade_n_tips", None)),
						"candidate_clade_modern_tips": scalar_int(tree.get("candidate_clade_modern_tips", None)),
						"candidate_clade_archaic_tips": scalar_int(tree.get("candidate_clade_archaic_tips", None)),
						"candidate_clade_specificity": scalar_float(tree.get("candidate_clade_specificity", None)),
						"candidate_clade_tips": scalar_text(tree.get("candidate_clade_tips", None)),
						"tree_newick": tree.get("tree_newick", ""),
						"raw_file": str(f),
					}


def init_db(con):
	con.executescript("""
    PRAGMA journal_mode=DELETE; PRAGMA synchronous=FULL; PRAGMA temp_store=FILE; PRAGMA user_version=2026;
    DROP VIEW IF EXISTS segments_evidence;
    DROP TABLE IF EXISTS segments_raw; DROP TABLE IF EXISTS loci; DROP TABLE IF EXISTS coord_catalog;
    DROP TABLE IF EXISTS segment_catalog; DROP TABLE IF EXISTS segments; DROP TABLE IF EXISTS carriers;
    DROP TABLE IF EXISTS database_summary;
    DROP TABLE IF EXISTS phyml_haplotype_carriers; DROP TABLE IF EXISTS locus_sample_support;
    DROP TABLE IF EXISTS locus_method_support; DROP TABLE IF EXISTS locus_trajectory; DROP TABLE IF EXISTS locus_evidence;
    DROP TABLE IF EXISTS method_runs;
    DROP TABLE IF EXISTS sample_populations; DROP TABLE IF EXISTS sample_burden; DROP TABLE IF EXISTS burden_union_chr;
    CREATE TABLE segments_raw(dataset_id TEXT,sample_id TEXT,method TEXT,source TEXT,source_class TEXT,reference_role TEXT,chr TEXT,start INTEGER,end INTEGER,length_bp INTEGER,haplotype TEXT,method_haplotype_index TEXT,score REAL,posterior REAL,trait TEXT,locus_id TEXT,genome_build TEXT,batch_id TEXT,raw_file TEXT);
    CREATE TABLE loci(dataset_id TEXT,method TEXT,trait TEXT,locus_id TEXT,chr TEXT,start INTEGER,end INTEGER,input_start INTEGER,input_end INTEGER,analysis_start INTEGER,analysis_end INTEGER,flank_bp INTEGER,anchor_pos INTEGER,selection_method TEXT,ld_r2_threshold REAL,n_search_sites INTEGER,n_ld_sites INTEGER,n_ancestral_sites INTEGER,n_sites INTEGER,best_archaic TEXT,best_lineage TEXT,n_compared INTEGER,n_match INTEGER,prop_match REAL,source TEXT,source_class TEXT,status TEXT,direct_match_pass INTEGER,n_carriers INTEGER,genome_build TEXT,tree_status TEXT,tree_file TEXT,stats_file TEXT,plot_file TEXT,n_bootstrap_nodes INTEGER,bootstrap_min REAL,bootstrap_median REAL,bootstrap_max REAL,candidate_lineage TEXT,tree_has_ancestral_outgroup INTEGER,candidate_clade_pass INTEGER,candidate_clade_rule TEXT,candidate_clade_bootstrap REAL,candidate_clade_node TEXT,candidate_clade_side TEXT,candidate_clade_n_tips INTEGER,candidate_clade_modern_tips INTEGER,candidate_clade_archaic_tips INTEGER,candidate_clade_specificity REAL,candidate_clade_tips TEXT,tree_newick TEXT,n_panel INTEGER,n_panel_copies INTEGER,n_lead_called INTEGER,n_lead_called_copies INTEGER,n_recurrent_copies INTEGER,n_singleton_copies INTEGER,input_status TEXT,ld_status TEXT,alignment_status TEXT,tree_run_status TEXT,topology_status TEXT,sequence_qc TEXT,lead_retained_for_tree INTEGER,high_ld_sites_retained_for_tree INTEGER,unique_archaic_sequences INTEGER,high_ld_max_gap_bp INTEGER,site_missingness_audit TEXT,raw_file TEXT);
    CREATE TABLE coord_catalog(dataset_id TEXT,genome_build TEXT,source_class TEXT,chr TEXT,start INTEGER,end INTEGER,segment_code TEXT);
    CREATE TABLE segment_catalog(segment_code TEXT PRIMARY KEY,dataset_id TEXT,source_class TEXT,chr TEXT,start INTEGER,end INTEGER,length_bp INTEGER,n_distinct_boundaries INTEGER,genome_build TEXT);
    CREATE TABLE sample_populations(dataset_id TEXT,sample_id TEXT,population TEXT,super_population TEXT,raw_file TEXT,PRIMARY KEY(dataset_id,sample_id));
    CREATE TABLE phyml_haplotype_carriers(dataset_id TEXT,genome_build TEXT,locus_id TEXT,sample_id TEXT,haplotype TEXT,hap_id TEXT,best_archaic TEXT,best_lineage TEXT,candidate_lineage TEXT,n_compared INTEGER,n_match INTEGER,prop_match REAL,direct_match_pass INTEGER,tree_candidate_pass INTEGER,raw_file TEXT);
    CREATE TABLE method_runs(dataset_id TEXT,genome_build TEXT,chr TEXT,method TEXT,status TEXT,evidence_eligible INTEGER,availability_note TEXT,raw_file TEXT);
    CREATE TABLE database_summary(dataset_id TEXT,genome_build TEXT,metric TEXT,value INTEGER,PRIMARY KEY(dataset_id,genome_build,metric));
    """)


def insert_df(con, table, df):
	if df is None or df.empty:
		return 0
	df.to_sql(table, con, if_exists="append", index=False, method="multi", chunksize=1000)
	return len(df)




def reciprocal(a, b, c, d):
	ov = max(0, min(b, d) - max(a, c))
	if ov <= 0:
		return 0.0
	return min(ov / max(1, b - a), ov / max(1, d - c))


def code_for(dataset, build, src, chrom, start, end):
	tag = {
		"Neanderthal": "NEA",
		"Denisovan": "DEN",
		"Mosaic": "MOS",
		"Ghost/Unknown": "GST",
		"Archaic/Other": "ARC",
	}.get(src, "ARC")
	dig = (
		hashlib.sha1(f"{dataset}|{build}|{src}|{chrom}|{round(start, -3)}|{round(end, -3)}".encode())
		.hexdigest()[:10]
		.upper()
	)
	return f"GU{build.replace('GRCh', '')}-{tag}-{chrom}-{dig}"


def build_catalog(con, build, threshold):
	groups = con.execute(
		'SELECT DISTINCT dataset_id,genome_build,source_class,chr FROM segments_raw WHERE ' + evidence_sql() + ' ORDER BY dataset_id,genome_build,source_class, CASE WHEN chr GLOB "[0-9]*" THEN CAST(chr AS INT) ELSE 99 END, chr'
	).fetchall()
	ins_coord = []
	ins_cat = []
	for dataset, group_build, src, chrom in groups:
		cur = con.execute(
			"SELECT DISTINCT start,end FROM segments_raw WHERE dataset_id=? AND genome_build=? AND source_class=? AND chr=? AND " + evidence_sql() + " ORDER BY start,end",
			(dataset, group_build, src, chrom),
		)
		coords = cur.fetchall()
		cluster = []
		mean_s = mean_e = None

		def flush():
			nonlocal cluster, mean_s, mean_e, ins_coord, ins_cat
			if not cluster:
				return
			rs = int(round(mean_s))
			re = int(round(mean_e))
			code = code_for(dataset, group_build, src, chrom, rs, re)
			ins_cat.append((code, dataset, src, chrom, rs, re, re - rs, len(cluster), group_build))
			ins_coord.extend((dataset, group_build, src, chrom, s, e, code) for s, e in cluster)
			if len(ins_coord) >= 100000:
				con.executemany("INSERT INTO coord_catalog VALUES(?,?,?,?,?,?,?)", ins_coord)
				ins_coord = []
			if len(ins_cat) >= 10000:
				con.executemany("INSERT OR REPLACE INTO segment_catalog VALUES(?,?,?,?,?,?,?,?,?)", ins_cat)
				ins_cat = []
			cluster = []
			mean_s = mean_e = None

		for s, e in coords:
			if not cluster:
				cluster = [(s, e)]
				mean_s = float(s)
				mean_e = float(e)
				continue
			if reciprocal(s, e, mean_s, mean_e) >= threshold:
				n = len(cluster)
				cluster.append((s, e))
				mean_s = (mean_s * n + s) / (n + 1)
				mean_e = (mean_e * n + e) / (n + 1)
			else:
				flush()
				cluster = [(s, e)]
				mean_s = float(s)
				mean_e = float(e)
		flush()
	if ins_coord:
		con.executemany("INSERT INTO coord_catalog VALUES(?,?,?,?,?,?,?)", ins_coord)
	if ins_cat:
		con.executemany("INSERT OR REPLACE INTO segment_catalog VALUES(?,?,?,?,?,?,?,?,?)", ins_cat)
	con.commit()
	con.executescript("""
      CREATE INDEX idx_coord ON coord_catalog(dataset_id,genome_build,source_class,chr,start,end);
      CREATE TABLE segments AS SELECT s.*, c.segment_code FROM segments_raw s LEFT JOIN coord_catalog c USING(dataset_id,genome_build,source_class,chr,start,end);
      CREATE VIEW segments_evidence AS SELECT * FROM segments WHERE __EVIDENCE_SQL__;
      CREATE INDEX idx_segments_sample ON segments(dataset_id,sample_id); CREATE INDEX idx_segments_region ON segments(dataset_id,genome_build,chr,start,end); CREATE INDEX idx_segments_code ON segments(segment_code);
      CREATE INDEX idx_segments_method_source ON segments(dataset_id,genome_build,method,source_class,chr,start);
      CREATE INDEX idx_segments_build_region ON segments(genome_build,chr,start,end);
      CREATE INDEX idx_segments_build_method ON segments(genome_build,dataset_id,method,source_class,chr,start);
      CREATE TABLE carriers AS
        SELECT dataset_id,sample_id,segment_code,source_class,
               MIN(2, CASE WHEN COUNT(DISTINCT CASE WHEN haplotype NOT IN ('','nan','None') THEN haplotype END)>0 THEN COUNT(DISTINCT CASE WHEN haplotype NOT IN ('','nan','None') THEN haplotype END) ELSE 1 END) AS dosage,
               COUNT(DISTINCT method) AS n_methods,
               GROUP_CONCAT(DISTINCT method) AS methods_support,
               MAX(score) AS max_score, MAX(posterior) AS max_posterior
        FROM segments_evidence WHERE sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None') GROUP BY dataset_id,sample_id,segment_code,source_class;
      CREATE INDEX idx_carriers_sample ON carriers(dataset_id,sample_id); CREATE INDEX idx_carriers_code ON carriers(segment_code);
    """.replace("__EVIDENCE_SQL__", evidence_sql()))
	con.commit()






def interval_summary(intervals, weighted=False):
	if not intervals:
		return 0, 0, 0, 0
	x = sorted(intervals, key=lambda z: (z[0], z[1]))
	merged = 0
	total = 0
	s = e = None
	for start, end, *_ in x:
		if s is None:
			s, e = start, end
		elif start <= e:
			e = max(e, end)
		else:
			merged += 1
			total += e - s
			s, e = start, end
	merged += 1
	total += e - s
	if not weighted:
		return len(x), merged, total, total
	events = {}
	for start, end, dose in x:
		dose = max(1, min(2, int(dose or 1)))
		events.setdefault(start, []).append((1, dose))
		events.setdefault(end, []).append((-1, dose))
	active = {1: 0, 2: 0}
	prev = None
	dosage_bp = 0
	for pos in sorted(events):
		if prev is not None:
			dose = 2 if active[2] else (1 if active[1] else 0)
			dosage_bp += (pos - prev) * dose
		for delta, dose in events[pos]:
			active[dose] += delta
		prev = pos
	return len(x), merged, total, dosage_bp


def insert_interval_burden(con, query, burden_type, weighted=False):
	cur = con.execute(query)
	current = None
	intervals = []
	rows = []

	def flush():
		nonlocal intervals, rows
		if current is None:
			return
		n, merged, total, dosage = interval_summary(intervals, weighted)
		rows.append((*current, burden_type, n, merged, total, dosage))
		intervals = []
		if len(rows) >= 10000:
			con.executemany("INSERT INTO burden_union_chr VALUES(?,?,?,?,?,?,?,?,?,?,?)", rows)
			rows = []

	for row in cur:
		key = tuple(row[:6])
		interval = tuple(row[6:])
		if current is not None and key != current:
			flush()
		current = key
		intervals.append(interval)
	flush()
	if rows:
		con.executemany("INSERT INTO burden_union_chr VALUES(?,?,?,?,?,?,?,?,?,?,?)", rows)


def build_burden(con):
	con.executescript("""
      DROP TABLE IF EXISTS sample_burden;
      CREATE TABLE sample_burden(
        dataset_id TEXT,sample_id TEXT,method TEXT,source_class TEXT,genome_build TEXT,burden_type TEXT,
        n_input_segments INTEGER,n_merged_intervals INTEGER,n_chromosomes INTEGER,total_bp INTEGER,dosage_bp INTEGER
      );
      INSERT INTO sample_burden
        SELECT dataset_id,sample_id,method,source_class,genome_build,'raw_call',COUNT(*),NULL,COUNT(DISTINCT chr),SUM(length_bp),SUM(length_bp)
        FROM segments_evidence WHERE sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None')
        GROUP BY dataset_id,sample_id,method,source_class,genome_build;
      CREATE TEMP TABLE burden_union_chr(
        dataset_id TEXT,sample_id TEXT,method TEXT,source_class TEXT,genome_build TEXT,chr TEXT,burden_type TEXT,
        n_input_segments INTEGER,n_merged_intervals INTEGER,total_bp INTEGER,dosage_bp INTEGER
      );
    """)
	# The biological source_class merges support across Neanderthal references,
	# preventing the same Neanderthal tract from being counted repeatedly.
	insert_interval_burden(
		con,
		"""
      SELECT dataset_id,sample_id,method,source_class,genome_build,chr,start,end
      FROM segments_evidence WHERE sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None')
      ORDER BY dataset_id,sample_id,method,source_class,genome_build,chr,start,end
    """,
		"nonredundant_union",
		False,
	)
	insert_interval_burden(
		con,
		"""
      SELECT c.dataset_id,c.sample_id,'consensus',c.source_class,s.genome_build,s.chr,s.start,s.end,c.dosage
      FROM carriers c JOIN segment_catalog s
        ON s.dataset_id=c.dataset_id AND s.segment_code=c.segment_code AND s.source_class=c.source_class
      WHERE c.n_methods>=2
      ORDER BY c.dataset_id,c.sample_id,c.source_class,s.genome_build,s.chr,s.start,s.end
    """,
		"consensus_catalog",
		True,
	)
	con.execute("""
      INSERT INTO sample_burden
      SELECT dataset_id,sample_id,method,source_class,genome_build,burden_type,
             SUM(n_input_segments),SUM(n_merged_intervals),COUNT(DISTINCT chr),SUM(total_bp),SUM(dosage_bp)
      FROM burden_union_chr GROUP BY dataset_id,sample_id,method,source_class,genome_build,burden_type
    """)
	con.executescript("""
      CREATE INDEX idx_sample_burden_sample ON sample_burden(dataset_id,sample_id);
      CREATE INDEX idx_sample_burden_filter ON sample_burden(genome_build,dataset_id,burden_type,method,source_class);
    """)
	con.commit()


def wilson_interval(successes, total, z=1.959963984540054):
	if not total:
		return (None, None, None)
	p = successes / total
	z2 = z * z
	denominator = 1 + z2 / total
	centre = (p + z2 / (2 * total)) / denominator
	half = z * math.sqrt((p * (1 - p) + z2 / (4 * total)) / total) / denominator
	return p, max(0.0, centre - half), min(1.0, centre + half)


def build_locus_support(con, trajectory_bins=120):
	"""Align PhyML carriers with individual tract calls and build scalable curves."""
	con.executescript("""
      CREATE INDEX idx_phyml_carrier_locus ON phyml_haplotype_carriers(dataset_id,genome_build,locus_id,sample_id);
      CREATE INDEX idx_phyml_carrier_match ON phyml_haplotype_carriers(dataset_id,genome_build,locus_id,tree_candidate_pass);
      DROP TABLE IF EXISTS locus_sample_support;
      CREATE TABLE locus_sample_support(
        dataset_id TEXT,genome_build TEXT,locus_id TEXT,chr TEXT,source_class TEXT,
        sample_id TEXT,population TEXT,super_population TEXT,
        phyml INTEGER,ibdmix INTEGER,trace INTEGER,n_methods INTEGER,
        methods_support TEXT,matched_haplotypes INTEGER,best_lineages TEXT,matched_source_classes TEXT,max_prop_match REAL,
        max_segment_score REAL,max_posterior REAL
      );
      WITH tested AS (
        SELECT pc.dataset_id,pc.genome_build,pc.locus_id,l.chr,l.source_class,pc.sample_id,
               p.population,p.super_population,
               MAX(pc.tree_candidate_pass) AS phyml,
               COUNT(DISTINCT CASE WHEN pc.tree_candidate_pass=1 THEN pc.hap_id END) AS matched_haplotypes,
               GROUP_CONCAT(DISTINCT CASE WHEN pc.tree_candidate_pass=1 THEN pc.candidate_lineage END) AS best_lineages,
               GROUP_CONCAT(DISTINCT CASE WHEN pc.tree_candidate_pass=1 THEN
                 CASE
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%neander%' AND LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%denis%' THEN 'Neanderthal,Denisovan'
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%neander%' OR LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%altai%' OR LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%vindija%' OR LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%chagyr%' THEN 'Neanderthal'
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%denis%' THEN 'Denisovan'
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%mosaic%' THEN 'Mosaic'
                   WHEN LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%ghost%' OR LOWER(COALESCE(pc.candidate_lineage,'')) LIKE '%unknown%' THEN 'Ghost/Unknown'
                   ELSE 'Archaic/Other'
                 END END) AS matched_source_classes,
               MAX(CASE WHEN pc.tree_candidate_pass=1 THEN pc.prop_match END) AS max_prop_match
        FROM phyml_haplotype_carriers pc
        JOIN loci l ON l.dataset_id=pc.dataset_id AND l.genome_build=pc.genome_build
                   AND l.locus_id=pc.locus_id AND l.method='phyml'
        LEFT JOIN sample_populations p ON p.dataset_id=pc.dataset_id AND p.sample_id=pc.sample_id
        GROUP BY pc.dataset_id,pc.genome_build,pc.locus_id,l.chr,l.source_class,pc.sample_id,p.population,p.super_population
      ), flags AS (
        SELECT t.dataset_id,t.genome_build,t.locus_id,t.sample_id,
               MAX(s.method='ibdmix') AS ibdmix,MAX(s.method='trace') AS trace,
               MAX(s.score) AS max_segment_score,MAX(s.posterior) AS max_posterior
        FROM tested t
        JOIN loci l ON l.dataset_id=t.dataset_id AND l.genome_build=t.genome_build AND l.locus_id=t.locus_id
        JOIN segments_evidence s ON s.dataset_id=t.dataset_id AND s.genome_build=t.genome_build AND s.sample_id=t.sample_id
                       AND s.chr=l.chr AND (
                         s.method='trace'
                         OR (t.phyml=1 AND INSTR(','||COALESCE(t.matched_source_classes,'')||',',','||s.source_class||',')>0)
                         OR (t.phyml=0 AND s.source_class=l.source_class)
                       )
                       AND s.end>l.start AND s.start<l.end
                       AND s.method IN ('ibdmix','trace')
        GROUP BY t.dataset_id,t.genome_build,t.locus_id,t.sample_id
      )
      INSERT INTO locus_sample_support
      SELECT t.dataset_id,t.genome_build,t.locus_id,t.chr,t.source_class,t.sample_id,t.population,t.super_population,
             t.phyml,COALESCE(f.ibdmix,0),COALESCE(f.trace,0),
             t.phyml+COALESCE(f.ibdmix,0)+COALESCE(f.trace,0),
             TRIM((CASE WHEN t.phyml=1 THEN 'phyml,' ELSE '' END)||
                  (CASE WHEN COALESCE(f.ibdmix,0)=1 THEN 'ibdmix,' ELSE '' END)||
                  (CASE WHEN COALESCE(f.trace,0)=1 THEN 'trace,' ELSE '' END),','),
             t.matched_haplotypes,t.best_lineages,t.matched_source_classes,t.max_prop_match,f.max_segment_score,f.max_posterior
      FROM tested t LEFT JOIN flags f USING(dataset_id,genome_build,locus_id,sample_id);
      CREATE UNIQUE INDEX idx_locus_sample_support_pk ON locus_sample_support(dataset_id,genome_build,locus_id,sample_id);
      CREATE INDEX idx_locus_sample_support_filter ON locus_sample_support(dataset_id,genome_build,locus_id,population,n_methods);
      UPDATE loci SET n_carriers=(
        SELECT COUNT(*) FROM locus_sample_support s
        WHERE s.dataset_id=loci.dataset_id AND s.genome_build=loci.genome_build
          AND s.locus_id=loci.locus_id AND s.phyml=1
      ) WHERE method='phyml';
      DROP TABLE IF EXISTS locus_method_support;
      CREATE TABLE locus_method_support(
        dataset_id TEXT,genome_build TEXT,locus_id TEXT,chr TEXT,population TEXT,method TEXT,source_class TEXT,
        method_available INTEGER,evidence_eligible INTEGER,availability_note TEXT,
        n_tested INTEGER,n_carriers INTEGER,carrier_fraction REAL,ci_low REAL,ci_high REAL,
        n_phyml_carriers INTEGER,n_phyml_supported INTEGER,phyml_support_fraction REAL
      );
      DROP TABLE IF EXISTS locus_trajectory;
      CREATE TABLE locus_trajectory(
        dataset_id TEXT,genome_build TEXT,locus_id TEXT,chr TEXT,population TEXT,denominator_mode TEXT,method TEXT,source_class TEXT,
        bin_start INTEGER,bin_end INTEGER,bin_mid INTEGER,n_denominator INTEGER,n_carriers INTEGER,
        prevalence REAL,ci_low REAL,ci_high REAL
      );
    """)

	available = {
		(d, b, c, m)
		for d, b, c, m in con.execute(
			"SELECT DISTINCT dataset_id,genome_build,chr,method FROM segments_evidence WHERE method IN ('ibdmix','trace') UNION SELECT dataset_id,genome_build,chr,method FROM method_runs"
		)
	}
	run_info = {
		(d, b, c, m): (int(eligible or 0), note or "")
		for d, b, c, m, eligible, note in con.execute(
			"SELECT dataset_id,genome_build,chr,method,MAX(evidence_eligible),MAX(availability_note) FROM method_runs GROUP BY dataset_id,genome_build,chr,method"
		)
	}
	aggregate_sql = """
      SELECT dataset_id,genome_build,locus_id,chr,source_class,'ALL' AS population,COUNT(*) AS n_tested,
             SUM(phyml),SUM(ibdmix),SUM(trace),
             SUM(phyml*ibdmix),SUM(phyml*trace)
      FROM locus_sample_support GROUP BY dataset_id,genome_build,locus_id,chr,source_class
      UNION ALL
      SELECT dataset_id,genome_build,locus_id,chr,source_class,population,COUNT(*),
             SUM(phyml),SUM(ibdmix),SUM(trace),
             SUM(phyml*ibdmix),SUM(phyml*trace)
      FROM locus_sample_support WHERE population IS NOT NULL AND population!=''
      GROUP BY dataset_id,genome_build,locus_id,chr,source_class,population
    """
	support_rows = []
	for d, b, lid, chrom, src, pop, n, phy, ibd, trc, phy_ibd, phy_trc in con.execute(aggregate_sql):
		values = {"phyml": (phy, phy), "ibdmix": (ibd, phy_ibd), "trace": (trc, phy_trc)}
		for method, (carriers, intersection) in values.items():
			carriers = int(carriers or 0)
			intersection = int(intersection or 0)
			n = int(n or 0)
			phy = int(phy or 0)
			fraction, lo, hi = wilson_interval(carriers, n)
			method_source = "Ghost/Unknown" if method == "trace" else src
			key = (d, b, chrom, method)
			method_available = int(method == "phyml" or key in available)
			if method == "phyml":
				eligible, note = 1, "bootstrap-supported tree clade; pairwise match is descriptive only"
			elif key in run_info:
				eligible, note = run_info[key]
			elif chrom == "X" and method == "ibdmix":
				eligible, note = 0, "experimental male-X pseudo-diploid adaptation; excluded from GU evidence score"
			elif method == "ibdmix":
				eligible, note = 0, "IBDmix completion/provenance unavailable; excluded from GU evidence score"
			elif method == "trace":
				eligible, note = (
					1,
					"haplotype-node call; validity depends on the full-chromosome ARG backend and calibration",
				)
			else:
				eligible, note = 1, "segment result observed; completion metadata unavailable"
			support_rows.append(
				(
					d,
					b,
					lid,
					chrom,
					pop,
					method,
					method_source,
					method_available,
					int(eligible),
					note,
					n,
					carriers,
					fraction,
					lo,
					hi,
					phy,
					intersection,
					(intersection / phy if phy else None),
				)
			)
	con.executemany("INSERT INTO locus_method_support VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", support_rows)
	con.executescript("""
      CREATE UNIQUE INDEX idx_locus_method_support_pk ON locus_method_support(dataset_id,genome_build,locus_id,population,method);
      CREATE INDEX idx_locus_method_support_filter ON locus_method_support(dataset_id,genome_build,locus_id,population);
    """)

	trajectory_rows = []
	loci = con.execute("""
      SELECT dataset_id,genome_build,locus_id,chr,source_class,start,end,analysis_start,analysis_end
      FROM loci WHERE method='phyml' ORDER BY dataset_id,genome_build,locus_id
    """).fetchall()
	for d, b, lid, chrom, src, core_start, core_end, analysis_start, analysis_end in loci:
		ast = int(analysis_start if analysis_start is not None else core_start)
		aen = int(analysis_end if analysis_end is not None else core_end)
		if aen <= ast:
			continue
		bins = max(1, min(int(trajectory_bins), aen - ast))
		edges = [ast + (aen - ast) * i // bins for i in range(bins + 1)]
		denominators = {
			"ALL": con.execute(
				"""
          SELECT COUNT(*) FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=?
        """,
				(d, b, lid),
			).fetchone()[0]
		}
		phy_counts = {
			"ALL": con.execute(
				"""
          SELECT COALESCE(SUM(phyml),0) FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=?
        """,
				(d, b, lid),
			).fetchone()[0]
		}
		for pop, n, carriers in con.execute(
			"""
          SELECT population,COUNT(*),SUM(phyml) FROM locus_sample_support
          WHERE dataset_id=? AND genome_build=? AND locus_id=? AND population IS NOT NULL AND population!=''
          GROUP BY population
        """,
			(d, b, lid),
		):
			denominators[pop] = n
			phy_counts[pop] = carriers or 0
		diffs = {}

		def diff_for(method, pop, mode):
			return diffs.setdefault((method, pop, mode), [0] * (bins + 1))

		current = None
		ranges = []
		current_pop = ""
		current_phyml = 0

		def flush_ranges():
			nonlocal ranges
			if current is None or not ranges:
				return
			method, _sample = current
			merged = []
			for lo, hi in sorted(ranges):
				if merged and lo <= merged[-1][1] + 1:
					merged[-1][1] = max(merged[-1][1], hi)
				else:
					merged.append([lo, hi])
			for pop in ("ALL", current_pop):
				if pop != "ALL" and (not pop or pop not in denominators):
					continue
				modes = ["tested"] + (["phyml_positive"] if current_phyml else [])
				for mode in modes:
					diff = diff_for(method, pop, mode)
					for lo, hi in merged:
						diff[lo] += 1
						diff[hi + 1] -= 1
			ranges = []

		query = """
          SELECT s.method,s.sample_id,COALESCE(ls.population,''),ls.phyml,s.start,s.end
          FROM segments_evidence s JOIN locus_sample_support ls
            ON ls.dataset_id=s.dataset_id AND ls.genome_build=s.genome_build AND ls.sample_id=s.sample_id
           AND ls.locus_id=?
          WHERE s.dataset_id=? AND s.genome_build=? AND s.chr=? AND (
              s.method='trace'
              OR (ls.phyml=1 AND INSTR(','||COALESCE(ls.matched_source_classes,'')||',',','||s.source_class||',')>0)
              OR (ls.phyml=0 AND s.source_class=?)
            )
            AND s.method IN ('ibdmix','trace') AND s.end>? AND s.start<?
          ORDER BY s.method,s.sample_id,s.start,s.end
        """
		for method, sample, pop, is_phyml, start, end in con.execute(query, (lid, d, b, chrom, src, ast, aen)):
			key = (method, sample)
			if current is not None and key != current:
				flush_ranges()
			current = key
			current_pop = pop or ""
			current_phyml = int(is_phyml or 0)
			start = max(ast, int(start))
			end = min(aen, int(end))
			if end <= start:
				continue
			lo = max(0, min(bins - 1, bisect.bisect_right(edges, start) - 1))
			hi = max(0, min(bins - 1, bisect.bisect_left(edges, end) - 1))
			if hi >= lo:
				ranges.append((lo, hi))
		flush_ranges()
		locus_methods = {m for dd, bb, cc, m in available if (dd, bb, cc) == (d, b, chrom)}
		for pop, tested_total in denominators.items():
			phy_n = int(phy_counts.get(pop, 0))
			for mode, total in (("tested", tested_total), ("phyml_positive", phy_n)):
				if not total:
					continue
				for method in sorted(locus_methods):
					diff_for(method, pop, mode)
				for method in sorted(locus_methods):
					running = 0
					diff = diff_for(method, pop, mode)
					for i in range(bins):
						running += diff[i]
						fraction, lo, hi = wilson_interval(running, total)
						method_source = "Ghost/Unknown" if method == "trace" else src
						trajectory_rows.append(
							(
								d,
								b,
								lid,
								chrom,
								pop,
								mode,
								method,
								method_source,
								edges[i],
								edges[i + 1],
								(edges[i] + edges[i + 1]) // 2,
								total,
								running,
								fraction,
								lo,
								hi,
							)
						)
				phy_carriers = phy_n if mode == "tested" else total
				fraction, lo, hi = wilson_interval(phy_carriers, total)
				for i in range(bins):
					if edges[i + 1] <= core_start or edges[i] >= core_end:
						continue
					trajectory_rows.append(
						(
							d,
							b,
							lid,
							chrom,
							pop,
							mode,
							"phyml",
							src,
							edges[i],
							edges[i + 1],
							(edges[i] + edges[i + 1]) // 2,
							total,
							phy_carriers,
							fraction,
							lo,
							hi,
						)
					)
		if len(trajectory_rows) >= 100000:
			con.executemany("INSERT INTO locus_trajectory VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", trajectory_rows)
			trajectory_rows = []
	if trajectory_rows:
		con.executemany("INSERT INTO locus_trajectory VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", trajectory_rows)
	con.executescript("""
      CREATE INDEX idx_locus_trajectory_filter ON locus_trajectory(dataset_id,genome_build,locus_id,population,denominator_mode,method,bin_start);
    """)
	con.commit()


def build_locus_evidence(con):
	"""Build a QC-first evidence row without combining incomparable scores."""
	con.executescript("""
      DROP TABLE IF EXISTS locus_evidence;
      CREATE TABLE locus_evidence(
        dataset_id TEXT,genome_build TEXT,locus_id TEXT,chr TEXT,source_class TEXT,
        n_tested INTEGER,n_phyml_carriers INTEGER,carrier_fraction REAL,
        best_archaic TEXT,best_lineage TEXT,n_compared INTEGER,n_match INTEGER,prop_match REAL,
        sequence_information TEXT,candidate_clade_bootstrap REAL,candidate_clade_specificity REAL,
        independent_methods_available INTEGER,independent_methods_supporting INTEGER,independent_method_names TEXT,
        tract_support_mean REAL,ibdmix_max_lod REAL,trace_max_posterior REAL,
        evidence_score REAL,evidence_completeness REAL,evidence_tier TEXT,
        probability_calibrated INTEGER,score_interpretation TEXT,
        candidate_clade_pass INTEGER,region_coverage_fraction REAL,phyml_qc TEXT,evidence_state TEXT,evidence_summary TEXT
      );
    """)
	rows = []
	loci = con.execute("""
      SELECT dataset_id,genome_build,locus_id,chr,source_class,start,end,input_start,input_end,
             best_archaic,best_lineage,n_compared,n_match,prop_match,tree_status,candidate_lineage,tree_has_ancestral_outgroup,
             candidate_clade_pass,candidate_clade_bootstrap,candidate_clade_specificity,candidate_clade_rule
      FROM loci WHERE method='phyml' ORDER BY dataset_id,genome_build,locus_id
    """).fetchall()
	for (
		d,
		b,
		lid,
		chrom,
		src,
		start,
		end,
		input_start,
		input_end,
		best_archaic,
		best_lineage,
		n_compared,
		n_match,
		prop_match,
		tree_status,
		candidate_lineage,
		has_outgroup,
		clade_pass,
		clade_boot,
		clade_specificity,
		clade_rule,
	) in loci:
		if candidate_lineage:
			src = source_class(candidate_lineage)
		n_tested, n_carriers = con.execute(
			"""
          SELECT COUNT(*),COALESCE(SUM(phyml),0) FROM locus_sample_support
          WHERE dataset_id=? AND genome_build=? AND locus_id=?
        """,
			(d, b, lid),
		).fetchone()
		carrier_fraction = (n_carriers / n_tested) if n_tested else None
		method_rows = con.execute(
			"""
          SELECT method,phyml_support_fraction FROM locus_method_support
          WHERE dataset_id=? AND genome_build=? AND locus_id=? AND population='ALL'
            AND method!='phyml' AND method_available=1 AND evidence_eligible=1
          ORDER BY method
        """,
			(d, b, lid),
		).fetchall()
		fractions = [float(x[1]) for x in method_rows if x[1] is not None]
		tract_support = sum(fractions) / len(fractions) if fractions else None
		supporting = sum(float(x[1] or 0) > 0 for x in method_rows)
		names = ",".join(x[0] for x in method_rows)
		native = con.execute(
			"""
          SELECT MAX(CASE WHEN method='ibdmix' THEN score END),
                 MAX(CASE WHEN method='trace' THEN posterior END)
          FROM segments_evidence WHERE dataset_id=? AND genome_build=? AND chr=? AND end>? AND start<?
            AND (method='trace' OR source_class=?)
        """,
			(d, b, chrom, start, end, src),
		).fetchone()
		n_compared = int(n_compared or 0)
		n_match = int(n_match or 0)
		sequence_information = "adequate" if n_compared >= 10 else ("low" if n_compared > 0 else "missing")
		input_bp = max(0, int(input_end or 0) - int(input_start or 0))
		selected_bp = max(0, int(end or 0) - int(start or 0))
		region_coverage = (selected_bp / input_bp) if input_bp else None
		if tree_status != "complete":
			phyml_qc = "not_evaluable"
		elif sequence_information != "adequate":
			phyml_qc = "low_sequence_information"
		elif region_coverage is not None and region_coverage < 0.25:
			phyml_qc = "narrow_selected_region"
		elif int(clade_pass or 0) != 1 and clade_rule == "supported_mixed_archaic_edge_lineage_unresolved":
			phyml_qc = "archaic_lineage_unresolved"
		elif int(clade_pass or 0) != 1:
			phyml_qc = "no_supported_archaic_clade"
		elif int(has_outgroup or 0) != 1:
			phyml_qc = "unrooted_tree"
		else:
			phyml_qc = "tree_supported"
		if phyml_qc in {"not_evaluable", "low_sequence_information"}:
			state = "not_evaluable"
		elif phyml_qc == "narrow_selected_region":
			state = "requires_core_region_rerun"
		elif phyml_qc == "archaic_lineage_unresolved":
			state = "mixed_archaic_tree_support_lineage_unresolved"
		elif phyml_qc == "no_supported_archaic_clade":
			state = "phyml_not_supported"
		elif phyml_qc == "unrooted_tree":
			state = "provisional_tree_support_with_independent_overlap" if supporting else "provisional_tree_support"
		elif not method_rows:
			state = "tree_supported_no_independent_method"
		elif supporting:
			state = "tree_supported_with_independent_overlap"
		else:
			state = "tree_supported_without_independent_overlap"
		summary = (
			f"{phyml_qc}; tree carriers={n_carriers}/{n_tested}; "
			f"eligible callers with overlap={supporting}/{len(method_rows)}"
		)
		rows.append(
			(
				d,
				b,
				lid,
				chrom,
				src,
				n_tested,
				n_carriers,
				carrier_fraction,
				best_archaic,
				best_lineage,
				n_compared,
				n_match,
				prop_match,
				sequence_information,
				clade_boot,
				clade_specificity,
				len(method_rows),
				supporting,
				names,
				tract_support,
				*native,
				None,
				None,
				state,
				0,
				"No composite score: inspect tree QC and each caller in its native scale",
				int(clade_pass or 0),
				region_coverage,
				phyml_qc,
				state,
				summary,
			)
		)
	con.executemany(
		"INSERT INTO locus_evidence VALUES(" + ",".join("?" for _ in range(32)) + ")", rows
	)
	con.execute("ALTER TABLE locus_evidence ADD COLUMN n_panel INTEGER")
	con.execute("UPDATE locus_evidence SET n_panel=(SELECT n_panel FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN n_panel_copies INTEGER")
	con.execute("UPDATE locus_evidence SET n_panel_copies=(SELECT n_panel_copies FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN n_lead_called INTEGER")
	con.execute("UPDATE locus_evidence SET n_lead_called=(SELECT n_lead_called FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN n_lead_called_copies INTEGER")
	con.execute("UPDATE locus_evidence SET n_lead_called_copies=(SELECT n_lead_called_copies FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN n_recurrent_copies INTEGER")
	con.execute("UPDATE locus_evidence SET n_recurrent_copies=(SELECT n_recurrent_copies FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN n_singleton_copies INTEGER")
	con.execute("UPDATE locus_evidence SET n_singleton_copies=(SELECT n_singleton_copies FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN input_status TEXT")
	con.execute("UPDATE locus_evidence SET input_status=(SELECT input_status FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN ld_status TEXT")
	con.execute("UPDATE locus_evidence SET ld_status=(SELECT ld_status FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN alignment_status TEXT")
	con.execute("UPDATE locus_evidence SET alignment_status=(SELECT alignment_status FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN tree_run_status TEXT")
	con.execute("UPDATE locus_evidence SET tree_run_status=(SELECT tree_run_status FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN topology_status TEXT")
	con.execute("UPDATE locus_evidence SET topology_status=(SELECT topology_status FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN sequence_qc TEXT")
	con.execute("UPDATE locus_evidence SET sequence_qc=(SELECT sequence_qc FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN lead_retained_for_tree INTEGER")
	con.execute("UPDATE locus_evidence SET lead_retained_for_tree=(SELECT lead_retained_for_tree FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN high_ld_sites_retained_for_tree INTEGER")
	con.execute("UPDATE locus_evidence SET high_ld_sites_retained_for_tree=(SELECT high_ld_sites_retained_for_tree FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN unique_archaic_sequences INTEGER")
	con.execute("UPDATE locus_evidence SET unique_archaic_sequences=(SELECT unique_archaic_sequences FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN high_ld_max_gap_bp INTEGER")
	con.execute("UPDATE locus_evidence SET high_ld_max_gap_bp=(SELECT high_ld_max_gap_bp FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.execute("ALTER TABLE locus_evidence ADD COLUMN site_missingness_audit TEXT")
	con.execute("UPDATE locus_evidence SET site_missingness_audit=(SELECT site_missingness_audit FROM loci l WHERE l.dataset_id=locus_evidence.dataset_id AND l.genome_build=locus_evidence.genome_build AND l.locus_id=locus_evidence.locus_id AND l.method='phyml')")
	con.executescript("""
      CREATE UNIQUE INDEX idx_locus_evidence_pk ON locus_evidence(dataset_id,genome_build,locus_id);
      CREATE INDEX idx_locus_evidence_state ON locus_evidence(dataset_id,genome_build,evidence_state,phyml_qc);
    """)
	con.commit()


def build_database_summary(con):
	con.executescript("""
      DELETE FROM database_summary;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'segment_calls',COUNT(*) FROM segments GROUP BY dataset_id,genome_build;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'samples',COUNT(DISTINCT sample_id) FROM segments
        WHERE sample_id IS NOT NULL AND sample_id!='' GROUP BY dataset_id,genome_build;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'methods',COUNT(DISTINCT method) FROM (
          SELECT dataset_id,genome_build,method FROM segments UNION ALL SELECT dataset_id,genome_build,method FROM loci
          UNION ALL SELECT dataset_id,genome_build,method FROM method_runs WHERE status='complete'
        ) GROUP BY dataset_id,genome_build;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'completed_method_units',COUNT(*) FROM method_runs WHERE status='complete' GROUP BY dataset_id,genome_build;
      INSERT INTO database_summary
        SELECT dataset_id,genome_build,'phyml_loci',COUNT(*) FROM loci WHERE method='phyml' GROUP BY dataset_id,genome_build;
    """)
	con.commit()


def export_query(con, query, path):
	path.parent.mkdir(parents=True, exist_ok=True)
	tmp = path.with_name(f".{path.name}.part.{os.getpid()}")
	cur = con.execute(query)
	cols = [x[0] for x in cur.description]
	with gzip.open(tmp, "wt") as h:
		h.write("\t".join(cols) + "\n")
		while True:
			rows = cur.fetchmany(50000)
			if not rows:
				break
			for r in rows:
				h.write("\t".join("" if x is None else str(x) for x in r) + "\n")
	os.replace(tmp, path)


def normalize_results_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--analysis-root", required=True, type=Path)
	ap.add_argument("--output-dir", required=True, type=Path)
	ap.add_argument("--build", default="GRCh37")
	ap.add_argument("--reciprocal-overlap", type=float, default=0.5)
	ap.add_argument("--trajectory-bins", type=int, default=120)
	ap.add_argument("--database", type=Path, help="SQLite destination; defaults to <analysis-root>/final/gu.sqlite")
	ap.add_argument("--sample-panel", type=Path, help="Sample metadata with sample/pop/super_pop columns")
	ap.add_argument("--sample-panel-dataset", default="1kg", help="Dataset owning this sample panel; default preserves the 1KG interface")
	ap.add_argument(
		"--ibdmix-daf-mask",
		type=Path,
		default=os.environ.get("GU_IBDMIX_DAF_MASK"),
		help="Altai or Altai.2013 GRCh37 top-0.1%% DAF exclusion BED with exact-reference .json provenance sidecar",
	)
	args = ap.parse_args()
	if args.trajectory_bins < 20 or args.trajectory_bins > 1000:
		ap.error("--trajectory-bins must be between 20 and 1000")
	args.output_dir.mkdir(parents=True, exist_ok=True)
	if args.sample_panel_dataset.startswith("ukb-"):
		# The shared PhyML summary link needs a directory even when this UKB
		# pilot has only IBDmix calls. An empty directory certifies no method.
		(args.analysis_root / gu_method_directory("phyml")).mkdir(parents=True, exist_ok=True)
	packaged = package_core_results(args.analysis_root, args.output_dir)
	parse_root = args.output_dir
	db = (args.database or (args.analysis_root / "final" / "gu.sqlite")).resolve()
	db.parent.mkdir(parents=True, exist_ok=True)
	fd, tmp_name = tempfile.mkstemp(prefix="gu-normalize-", suffix=".sqlite")
	os.close(fd)
	tmp_db = Path(tmp_name)
	db_part = db.with_name(f".{db.name}.part.{os.getpid()}")
	con = None
	loci = []
	try:
		con = sqlite3.connect(tmp_db)
		init_db(con)
		counts = {}
		method_runs = discover_method_runs(args.analysis_root, args.build)
		ukb_datasets = {row[0] for row in method_runs if row[0].startswith("ukb-")}
		if ukb_datasets:
			ukb_prep = load_module("ukb_prepare.py")
			for dataset in sorted(ukb_datasets):
				if ukb_prep.validate_result_cohort(args.analysis_root, dataset) is None:
					raise ValueError(f"UKB method results lack archived preparation QC: {dataset}")
		if method_runs:
			con.executemany("INSERT INTO method_runs VALUES(?,?,?,?,?,?,?,?)", method_runs)
		con.execute("CREATE INDEX idx_method_runs_unit ON method_runs(dataset_id,genome_build,chr,method)")
		counts["completed_method_units"] = len(method_runs)
		populations = sample_population_table(args.sample_panel, args.sample_panel_dataset, args.analysis_root)
		if not populations.empty:
			insert_df(con, "sample_populations", populations)
		completed_ibdmix = {Path(row[7]).parent for row in method_runs if row[3] == "ibdmix" and row[4] == "complete"}
		ibdmix_audits = []
		if any(not dataset_from_path(path, "ibdmix").startswith("ukb-") for path in completed_ibdmix) and args.ibdmix_daf_mask is None:
			print(
				"[GU FINAL] WARNING: Altai DAF filter NOT applied: author exclusion BED unavailable; set GU_IBDMIX_DAF_MASK. Denisovan filters remain enabled.",
				flush=True,
			)
		for name, it in [
			(
				"ibdmix",
				iter_ibdmix(
					args.analysis_root,
					args.build,
					completed_ibdmix,
					args.output_dir / "ibdmix_filters",
					args.ibdmix_daf_mask,
					ibdmix_audits,
				),
			),
			("trace", iter_trace(args.analysis_root, args.build)),
		]:
			n = 0
			if it is not None:
				for d in it:
					n += insert_df(con, "segments_raw", portable_paths(d, db.parent, ["raw_file"]))
			counts[name] = n
			con.commit()
		con.execute(
			"CREATE TABLE ibdmix_filter_runs(dataset_id TEXT,genome_build TEXT,chr TEXT,native_file TEXT,filtered_file TEXT,status TEXT,daf_status TEXT,audit_json TEXT)"
		)
		con.executemany("INSERT INTO ibdmix_filter_runs VALUES(?,?,?,?,?,?,?,?)", ibdmix_audits)
		con.commit()
		loci = list(loci_tables(parse_root, args.build))
		if loci:
			insert_df(
				con,
				"loci",
				portable_paths(pd.DataFrame(loci), db.parent, ["raw_file", "tree_file", "stats_file", "plot_file"]),
			)
		n_phyml_carriers = 0
		for d in iter_phyml_carriers(parse_root, args.build):
			n_phyml_carriers += insert_df(con, "phyml_haplotype_carriers", portable_paths(d, db.parent, ["raw_file"]))
		counts["phyml_haplotype_copies"] = n_phyml_carriers
		con.execute("CREATE INDEX idx_loci_region ON loci(dataset_id,genome_build,chr,start,end)")
		con.commit()
		build_catalog(con, args.build, args.reciprocal_overlap)
		build_burden(con)
		build_locus_support(con, args.trajectory_bins)
		build_locus_evidence(con)
		build_database_summary(con)
		summary = args.output_dir / "summary"
		export_query(
			con,
			"SELECT * FROM segments ORDER BY dataset_id,genome_build,chr,start,end,sample_id",
			summary / "segments.tsv.gz",
		)
		export_query(
			con,
			"SELECT * FROM segment_catalog ORDER BY dataset_id,genome_build,chr,start,end",
			summary / "segment_catalog.tsv.gz",
		)
		export_query(
			con, "SELECT * FROM carriers ORDER BY dataset_id,sample_id,segment_code", summary / "carriers.tsv.gz"
		)
		export_query(con, "SELECT * FROM loci ORDER BY dataset_id,genome_build,chr,start,end", summary / "loci.tsv.gz")
		export_query(
			con,
			"SELECT * FROM phyml_haplotype_carriers ORDER BY dataset_id,genome_build,locus_id,sample_id,haplotype",
			summary / "phyml_haplotype_carriers.tsv.gz",
		)
		export_query(
			con,
			"SELECT * FROM locus_sample_support ORDER BY dataset_id,genome_build,locus_id,n_methods DESC,sample_id",
			summary / "locus_sample_support.tsv.gz",
		)
		export_query(
			con,
			"SELECT * FROM locus_method_support ORDER BY dataset_id,genome_build,locus_id,population,method",
			summary / "locus_method_support.tsv.gz",
		)
		export_query(
			con,
			"SELECT * FROM locus_evidence ORDER BY dataset_id,genome_build,locus_id",
			summary / "locus_evidence.tsv.gz",
		)
		export_query(
			con,
			"SELECT * FROM locus_trajectory ORDER BY dataset_id,genome_build,locus_id,population,denominator_mode,method,bin_start",
			summary / "locus_trajectory.tsv.gz",
		)
		export_query(
			con, "SELECT * FROM method_runs ORDER BY dataset_id,genome_build,chr,method", summary / "method_runs.tsv.gz"
		)
		export_query(
			con,
			"SELECT * FROM sample_populations ORDER BY dataset_id,population,sample_id",
			summary / "sample_populations.tsv.gz",
		)
		export_query(
			con,
			"SELECT * FROM sample_burden ORDER BY dataset_id,genome_build,burden_type,sample_id,method,source_class",
			summary / "sample_burden.tsv.gz",
		)
		export_query(
			con,
			"SELECT dataset_id,genome_build,chr,native_file,filtered_file,status,daf_status FROM ibdmix_filter_runs",
			summary / "ibdmix_filters.tsv",
		)
		for method in ("trace", "ibdmix"):
			export_query(
				con,
				f"SELECT * FROM segments WHERE method='{method}' ORDER BY dataset_id,genome_build,chr,start,end,sample_id",
				args.output_dir / method / "segments.tsv.gz",
			)
		export_query(
			con,
			"SELECT * FROM loci WHERE method='phyml' ORDER BY dataset_id,genome_build,chr,start,end",
			args.output_dir / "phyml/loci.tsv.gz",
		)
		check = con.execute("PRAGMA quick_check").fetchone()
		if not check or check[0] != "ok":
			raise RuntimeError(f"SQLite quick_check failed: {check}")
		con.close()
		con = None
		if db_part.exists():
			db_part.unlink()
		shutil.copy2(tmp_db, db_part)
		os.replace(db_part, db)
	finally:
		if con is not None:
			con.close()
		if tmp_db.exists():
			tmp_db.unlink()
		if db_part.exists():
			db_part.unlink()
	print("Normalized:", counts, "loci=", len(loci), "packaged_files=", packaged, "database=", db)


def normalize_results_cli():
	normalize_results_main()


# 🚩 introgression_density
"""Cache per-individual genomic coverage, merging overlapping IBDmix references.

Fractions are lineage-specific union bp / diploid physical span. Unphased calls do
not resolve homozygous tract dosage. Male non-PAR X has a haploid denominator.
Missing chromosomes and untested female X rows are NA, never a negative call.
"""

import argparse, csv, gzip, hashlib, json, os, re, sqlite3, tempfile, shutil, atexit, errno
from pathlib import Path
import numpy as np

load_module("0.common.py")
from gu_0_common import CHROM_LENGTHS, write_tsv_rows

SCHEMA = 2026
SUMMARY_REFERENCES = tuple(DEFAULT_ARCHAIC_REFERENCES)


def summary_references(references):
	"""Order only references explicitly recorded by completed runs."""
	actual = set(references)
	unknown = actual - set(SUPPORTED_ARCHAIC_REFERENCES)
	if unknown:
		raise ValueError("Unsupported archaic references in completed run: " + ", ".join(sorted(unknown)))
	return tuple(ref for ref in SUMMARY_REFERENCES if ref in actual) + tuple(sorted(actual - set(SUMMARY_REFERENCES)))


def density_input_stat(path):
	try:
		return path.stat() if path.is_file() else None
	except OSError as error:
		# Archived manifests can point to a drive that is no longer mounted.
		if error.errno not in (errno.ENOENT, errno.ENOTDIR, errno.ENODEV, errno.ENXIO):
			raise
		return None


def density_run_path(raw, database):
	"""Resolve current provenance paths relative to the result database."""
	path = Path(raw)
	return path if path.is_absolute() else database.parent / path


def density_run_signature(paths):
	"""Include exact saved calling rosters in the cache identity."""
	inputs = set(paths)
	for path in paths:
		inputs.update(path.parent / "samples" / name / "ALL.txt" for name in ("X_MALE", "X_PAR"))
		inputs.update(path.parent / "samples" / ("C" + str(chrom)) / "ALL.txt" for chrom in range(1, 23))
	records = []
	for path in sorted(inputs):
		stat = density_input_stat(path)
		records.append((str(path), stat.st_mtime_ns if stat else None, stat.st_size if stat else None))
	return hashlib.sha256(json.dumps(records).encode()).hexdigest()


def whole_chromosome_run(text, chrom):
	"""Recognize both native-VCF and cached modern_source run records."""
	records = [line.split("\t") for line in text.splitlines()]
	if any(row[0] in {"loci_file", "loci"} for row in records):
		return False
	if any(row[0] == "modern_vcf" for row in records):
		return True
	source = [row[2:] for row in records if len(row) >= 4 and row[:2] == ["modern_source", "chr" + str(chrom)]]
	scopes = [row[1] for row in source if row[0] == "loci"]
	if scopes:
		return all(scope == "whole_chromosome" for scope in scopes)
	# Native VCF source manifests omit a loci row for whole chromosomes.
	return any(row[:2] == ["format", "vcf"] for row in source)




def shared_reference_targets(reference_targets, references, chromosomes=tuple(map(str, range(1, 23)))):
	"""Only certify sample/chromosome pairs tested against every group member."""
	common = {}
	for chrom in chromosomes:
		samples = set.intersection(*(reference_targets.get(ref, {}).get(chrom, set()) for ref in references))
		if samples:
			common[chrom] = samples
	return common


def reference_summary(con, dataset, build, samples, lengths, targets, reference, x_male=False):
	"""Exact per-reference interval union on autosomes or male non-PAR X.

	Targets certify whole-chromosome runs.

	Untested individuals remain missing, including when provenance is unavailable.
	Do not infer a denominator from positive calls or extrapolate partial genomes.
	"""
	burden = dict.fromkeys(samples, 0)
	chromosomes = ("X",) if x_male else tuple(map(str, range(1, 23)))
	spans = {
		ch: x_nonpar_intervals(build, lengths[ch]) if x_male else [(0, lengths[ch])]
		for ch in chromosomes
		if ch in lengths
	}
	sources = (reference,)
	placeholders = ",".join("?" for _ in sources)
	scope_sql = "chr='X'" if x_male else "chr IN (" + ",".join("'" + ch + "'" for ch in chromosomes) + ")"
	rows = con.execute(
		f"SELECT sample_id,chr,start,end FROM segments WHERE dataset_id=? AND genome_build=? AND method='ibdmix' AND source IN ({placeholders}) AND {scope_sql} ORDER BY sample_id,chr,start,end",
		(dataset, build, *sources),
	)
	for sample, chrom, left, right in merged_intervals(rows):
		if chrom in spans and sample in burden and sample in targets.get(chrom, set()):
			burden[sample] += sum(max(0, min(hi, right) - max(lo, left)) for lo, hi in spans[chrom])
	result = []
	for sample in samples:
		chroms = [ch for ch in chromosomes if ch in spans and sample in targets.get(ch, set())]
		physical = sum(hi - lo for ch in chroms for lo, hi in spans[ch])
		denominator = physical if x_male else 2 * physical
		result.append(
			dict(
				sample_id=sample,
				chromosomes=",".join(chroms),
				n_chromosomes=len(chroms),
				tested_bp=denominator,
				haploid_bp=physical,
				archaic_bp=burden[sample] if denominator else "",
				coverage_pct=100 * burden[sample] / denominator if denominator else "",
				coordinate_coverage_pct=100 * burden[sample] / physical if physical else "",
			)
		)
	return result


def merged_intervals(rows):
	"""Input sorted by sample, chromosome, start, end; yield union intervals."""
	key = None
	start = end = 0
	for sample, chrom, left, right in rows:
		if right <= left:
			continue
		current = (sample, chrom)
		if current == key and left <= end:
			end = max(end, right)
		else:
			if key is not None:
				yield (*key, start, end)
			key, start, end = current, left, right
	if key is not None:
		yield (*key, start, end)


def add_interval(covered, row, offset, length, bin_bp, left, right):
	left = max(0, left)
	right = min(length, right)
	if right <= left:
		return
	for b in range(left // bin_bp, (right - 1) // bin_bp + 1):
		covered[row, offset + b] += min(right, (b + 1) * bin_bp) - max(left, b * bin_bp)



def population_burdens(rows, metadata, dataset, build, scope, population_samples, filters):
	"""Saved descriptive data: zero calls require a certified tested denominator."""
	by_ref = {}
	for row in rows: by_ref.setdefault(row["reference"], {})[row["sample_id"]] = row
	groups = {"ALL": set(population_samples)}
	for sample in population_samples:
		pop, superpop = metadata.get(sample, ("UNKNOWN", "UNKNOWN"))
		groups.setdefault(pop, set()).add(sample)
		groups.setdefault(superpop, set()).add(sample)
		if pop in {"ESN", "GWD", "LWK", "MSL", "YRI"}: groups.setdefault("AFR (5)", set()).add(sample)
	result = []
	for ref, values in by_ref.items():
		for group, members in groups.items():
			measured = [values[s] for s in members if s in values and values[s]["tested_bp"] > 0 and values[s]["archaic_bp"] != ""]
			n = len(measured)
			chroms = sorted({r["chromosomes"] for r in measured})
			denoms = [r["tested_bp"] for r in measured]
			total = sum(r["archaic_bp"] for r in measured)
			result.append(dict(dataset_id=dataset, genome_build=build, scope=scope, group=group, reference=ref,
				reference_role=reference_role(ref),
				n_panel=len(members), n_tested=n, n_positive=sum(r["archaic_bp"] > 0 for r in measured), n_zero=sum(r["archaic_bp"] == 0 for r in measured), n_untested=len(members)-n,
				total_union_bp=total if n else None, mean_mb=total/n/1e6 if n else None,
				coverage_pct=sum(r["coverage_pct"] for r in measured)/n if n else None,
				coordinate_coverage_pct=sum(100*r["archaic_bp"]/r["haploid_bp"] for r in measured)/n if n else None,
				coordinate_coverage_definition="interval_union_bp / physical_haploid_span; not ancestry dosage or callable fraction",
				physical_denominator_bp=denoms[0] if denoms and len(set(denoms)) == 1 else None,
				denominator_min_bp=min(denoms) if denoms else None, denominator_max_bp=max(denoms) if denoms else None,
				chromosomes=";".join(chroms), comparable_scope=int(len(chroms) == 1),
				complete_chromosomes=int(bool(measured) and all(r["n_chromosomes"] == (1 if scope == "male_nonpar_X" else 22) for r in measured)),
				filter_status=json.dumps(filters, sort_keys=True), denominator_definition="physical_nonPAR_haploid" if scope == "male_nonpar_X" else "physical_autosomes_diploid_not_dosage_or_callable",
				mask_history="see_native_qc;unknown_if_not_saved", model_calibration="experimental_pseudodiploid_LOD_uncalibrated" if scope == "male_nonpar_X" else "existing_autosomal_profile"))
	return result


def paired_male_burdens(x_rows, auto_rows, metadata, dataset, build):
	"""Same people, same reference, certified full X and 22-autosome coverage."""
	autos = {(r["reference"], r["sample_id"]): r for r in auto_rows if r["n_chromosomes"] == 22 and r["tested_bp"] > 0}
	pairs = []
	for x in x_rows:
		a = autos.get((x["reference"], x["sample_id"]))
		if a is None or x["n_chromosomes"] != 1 or x["tested_bp"] <= 0: continue
		pairs.append((x, a))
	groups = {}
	for x, a in pairs:
		pop, sp = metadata.get(x["sample_id"], ("UNKNOWN", "UNKNOWN"))
		for group in set(("ALL", pop, sp)) | ({"AFR (5)"} if pop in {"ESN", "GWD", "LWK", "MSL", "YRI"} else set()):
			groups.setdefault((x["reference"], group), []).append((x, a))
	result = []
	for (ref, group), rr in groups.items():
		x_pct = sum(x["coverage_pct"] for x,a in rr)/len(rr)
		a_pct = sum(a["coverage_pct"] for x,a in rr)/len(rr)
		x_coord = sum(100*x["archaic_bp"]/x["haploid_bp"] for x,a in rr)/len(rr)
		a_coord = sum(100*a["archaic_bp"]/a["haploid_bp"] for x,a in rr)/len(rr)
		result.append(dict(dataset_id=dataset, genome_build=build, reference=ref, group=group, n_same_males=len(rr),
			x_mean_mb=sum(x["archaic_bp"] for x,a in rr)/len(rr)/1e6, autosome_mean_mb=sum(a["archaic_bp"] for x,a in rr)/len(rr)/1e6,
			x_coverage_pct=x_pct, autosome_coverage_pct=a_pct,
			x_coordinate_coverage_pct=x_coord, autosome_coordinate_coverage_pct=a_coord,
			x_to_autosome_reported_coverage_ratio=x_pct/a_pct if a_pct else None,
			x_to_autosome_coordinate_coverage_ratio=x_coord/a_coord if a_coord else None,
			x_physical_bp=rr[0][0]["tested_bp"], autosome_diploid_bp=rr[0][1]["tested_bp"],
			definition="same_men;22_autosomes_and_nonPAR_X;coordinate_coverage=union/one_physical_span;reported_autosome_coverage=union/two_spans;not_ancestry_dosage;X_LOD_uncalibrated;filters_differ"))
	return result


def prepare(database, output, sample_panel=None, bin_bp=5_000_000):
	stamp = database.stat()
	con = sqlite3.connect(f"file:{database.resolve()}?mode=ro", uri=True)
	run_paths = {
		raw: density_run_path(raw, database)
		for (raw,) in con.execute("SELECT DISTINCT raw_file FROM method_runs WHERE method='ibdmix' AND status='complete'")
	}
	signature = {
		"schema": SCHEMA,
		"database": str(database.resolve()),
		"mtime_ns": stamp.st_mtime_ns,
		"size": stamp.st_size,
		"bin_bp": bin_bp,
		"sample_panel": str(sample_panel) if sample_panel else "",
		"panel_mtime_ns": sample_panel.stat().st_mtime_ns if sample_panel and sample_panel.exists() else 0,
		"run_inputs_sha256": density_run_signature(run_paths.values()),
	}
	version = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:16]
	output.mkdir(parents=True, exist_ok=True)
	units = con.execute(
		"SELECT DISTINCT dataset_id,genome_build FROM method_runs WHERE method='ibdmix' AND status='complete'"
	).fetchall()
	missing = any(
		not (output / re.sub(r"[^A-Za-z0-9_.-]", "_", dataset + "." + build) / version / "manifest.json").exists()
		for dataset, build in units
	)
	local_database = None
	if missing:
		# Random indexed reads on /mnt/* are much slower than one sequential
		# copy followed by local SQLite reads. Never mutate the source database.
		fd, name = tempfile.mkstemp(prefix="gu-density-", suffix=".sqlite")
		os.close(fd)
		local_database = Path(name)
		atexit.register(lambda: local_database.unlink(missing_ok=True))
		print("DENSITY staging SQLite for local interval aggregation", flush=True)
		shutil.copyfile(database, local_database)
		current = database.stat()
		if (current.st_mtime_ns, current.st_size) != (stamp.st_mtime_ns, stamp.st_size):
			raise RuntimeError("database changed while preparing density; retry")
		con.close()
		con = sqlite3.connect(local_database)
	for dataset, build in units:
		lengths = CHROM_LENGTHS.get(build.replace("GRCh", ""))
		if not lengths:
			continue
		slug = re.sub(r"[^A-Za-z0-9_.-]", "_", dataset + "." + build)
		root = output / slug
		dest = root / version
		pointer = root / "current.tsv"
		if (dest / "manifest.json").exists():
			write_tsv_rows(pointer, ["directory"], [{"directory": version}])
			print(f"DENSITY cached {dataset} {build}", flush=True)
			continue
		root.mkdir(parents=True, exist_ok=True)
		stage = Path(tempfile.mkdtemp(prefix=".building-", dir=root))
		metadata = {
			r[0]: (r[1] or "UNKNOWN", r[2] or "UNKNOWN")
			for r in con.execute(
				"SELECT sample_id,population,super_population FROM sample_populations WHERE dataset_id=?", (dataset,)
			)
		}
		called_pairs = list(
			con.execute(
				"SELECT DISTINCT sample_id,chr FROM segments WHERE dataset_id=? AND genome_build=? AND method='ibdmix'",
				(dataset, build),
			)
		)
		called = {r[0] for r in called_pairs}
		runs = con.execute(
			"SELECT chr,evidence_eligible,raw_file,availability_note FROM method_runs WHERE dataset_id=? AND genome_build=? AND method='ibdmix' AND status='complete'",
			(dataset, build),
		).fetchall()
		tested = {r[0] for r in runs}
		x_male = False
		targets = {}
		neanderthal_tested = set()
		lineage_targets = {"Neanderthal": {}, "Denisovan": {}}
		reference_targets = {}
		# Use the database's provenance snapshot. Raw metadata may already have
		# changed during a rerun while this database still contains old calls.
		profiles = sorted({str(r[3]).split("profile=", 1)[1] for r in runs if "profile=" in str(r[3])})
		for chrom, _, raw, note in runs:
			p = run_paths[raw]
			if not density_input_stat(p):
				print(f"DENSITY WARNING {dataset} {build} chr{chrom}: run metadata missing: {p}; coverage remains N/A", flush=True)
				continue
			text = p.read_text()
			run_targets = set()
			if chrom == "X" and "male_haploid_nonpar" in text:
				x_male = True
			# Saved caller rosters survive relocation and include zero-call people.
			# Prefer them over external input files, which may no longer exist.
			unit = ("X_MALE" if "male_haploid_nonpar" in text else "X_PAR") if chrom == "X" else "C" + chrom
			roster = p.parent / "samples" / unit / "ALL.txt"
			if roster.is_file():
				run_targets = set(roster.read_text().split())
			if run_targets:
				targets.setdefault(chrom, set()).update(run_targets)
			if not run_targets:
				print(f"DENSITY WARNING {dataset} {build} chr{chrom}: no target sample roster for {p}; coverage remains N/A", flush=True)
			refs = next((line.split("\t", 1)[1].split() for line in text.splitlines() if line.startswith("refs\t")), [])
			# Locus-only calls cannot certify an entire chromosome denominator.
			whole = whole_chromosome_run(text, chrom)
			neand_refs = {ref for ref in summary_references(refs) if source_class(ref) == "Neanderthal"}
			if neand_refs:
				neanderthal_tested.add(chrom)
			# Background-only Denisova runs export no Denisovan ancestry calls.
			# Their absence in segments is unmeasured, not zero coverage.
			background_only = any(
				line.split("\t")[:2] == ["background_filter", "1"] for line in text.splitlines()
			) and not any(line.split("\t")[:2] == ["export_denisovan", "1"] for line in text.splitlines())
			if (
				whole
				and run_targets
				and (chrom in map(str, range(1, 23)) or (chrom == "X" and "male_haploid_nonpar" in text))
			):
				for ref in summary_references(refs):
					if background_only and ref.startswith("Denisova"):
						continue
					reference_targets.setdefault(ref, {}).setdefault(chrom, set()).update(run_targets)
			for lineage, has_refs in [
				("Neanderthal", bool(neand_refs)),
				("Denisovan", not background_only and any("denis" in ref.lower() for ref in refs)),
			]:
				if has_refs and whole and run_targets:
					lineage_targets[lineage].setdefault(chrom, set()).update(run_targets)
		tested &= neanderthal_tested  # a Denisova-only run did not test Neanderthal ancestry
		# Use the actual caller target, not the larger annotation panel.
		samples = sorted(
			called | set().union(*targets.values()),
			key=lambda s: (metadata.get(s, ("UNKNOWN", "UNKNOWN"))[1], metadata.get(s, ("UNKNOWN", "UNKNOWN"))[0], s),
		)
		sample_index = {s: i for i, s in enumerate(samples)}
		for chrom in tested:
			if chrom not in targets:
				targets[chrom] = {sample for sample, c in called_pairs if c == chrom}
		bins = []
		offsets = {}
		for chrom, length in lengths.items():
			offsets[chrom] = len(bins)
			for left in range(0, length, bin_bp):
				bins.append(
					dict(
						bin_index=len(bins),
						chr=chrom,
						start=left,
						end=min(length, left + bin_bp),
						tested=int(chrom in tested and (chrom != "X" or not x_male or any(max(left, lo) < min(left + bin_bp, hi) for lo, hi in x_nonpar_intervals(build, length)))),
						physical_bp=sum(max(0, min(length, left + bin_bp, hi) - max(left, lo)) for lo, hi in x_nonpar_intervals(build, length)) if chrom == "X" and x_male else min(length, left + bin_bp) - left,
					)
				)
		print(f"DENSITY building {dataset} {build}: {len(samples)} individuals x {len(bins)} bins", flush=True)
		widths = np.array([b["physical_bp"] for b in bins])
		ploidy = np.array([1 if b["chr"] == "X" and x_male else 2 for b in bins])
		unions = {}
		for lineage, filename in [("Neanderthal", "matrix.tsv.gz"), ("Denisovan", "denisovan_matrix.tsv.gz")]:
			covered = np.zeros((len(samples), len(bins)), dtype=np.int64)
			rows = con.execute(
				"SELECT sample_id,chr,start,end FROM segments WHERE dataset_id=? AND genome_build=? AND method='ibdmix' AND source_class=? ORDER BY sample_id,chr,start,end",
				(dataset, build, lineage),
			)
			unions[lineage] = 0
			for sample, chrom, left, right in merged_intervals(rows):
				if chrom not in lengths:
					continue
				spans = x_nonpar_intervals(build, lengths[chrom]) if chrom == "X" and x_male else [(0, lengths[chrom])]
				for lo, hi in spans:
					add_interval(covered, sample_index[sample], offsets[chrom], lengths[chrom], bin_bp, max(left, lo), min(right, hi))
				unions[lineage] += 1
			if np.any(covered > widths):
				raise ValueError("union coverage exceeds bin length")
			density = np.round(np.divide(100 * covered, widths * ploidy, out=np.full(covered.shape, np.nan), where=widths > 0), 3)
			for j, b in enumerate(bins):
				eligible = lineage_targets[lineage].get(b["chr"], set())
				for i, sample in enumerate(samples):
					if sample not in eligible:
						density[i, j] = np.nan
			with gzip.open(stage / filename, "wt") as f:
				np.savetxt(f, density, delimiter="\t", fmt="%.3f")
		write_tsv_rows(
			stage / "samples.tsv",
			["row_index", "sample_id", "population", "super_population"],
			[
				dict(
					row_index=i,
					sample_id=s,
					population=metadata.get(s, ("UNKNOWN", "UNKNOWN"))[0],
					super_population=metadata.get(s, ("UNKNOWN", "UNKNOWN"))[1],
				)
				for i, s in enumerate(samples)
			],
		)
		write_tsv_rows(stage / "bins.tsv", ["bin_index", "chr", "start", "end", "tested", "physical_bp"], bins)
		reference_rows = []
		for ref in summary_references(reference_targets):
			group_targets = shared_reference_targets(reference_targets, (ref,))
			if not group_targets:
				continue
			print(f"DENSITY summarizing {ref} autosomes: {len(group_targets)}/22 chromosomes", flush=True)
			reference_rows.extend(
				dict(row, reference=ref)
				for row in reference_summary(con, dataset, build, samples, lengths, group_targets, ref)
			)
		x_rows = []
		for ref in summary_references(reference_targets):
			group_targets = shared_reference_targets(reference_targets, (ref,), ("X",))
			if not group_targets:
				continue
			x_rows.extend(
				dict(row, reference=ref)
				for row in reference_summary(con, dataset, build, samples, lengths, group_targets, ref, x_male=True)
			)
		write_tsv_rows(
			stage / "archaic_summary_x.tsv",
			[
				"reference",
				"sample_id",
				"chromosomes",
				"n_chromosomes",
				"tested_bp",
				"haploid_bp",
				"archaic_bp",
				"coverage_pct",
				"coordinate_coverage_pct",
			],
			x_rows,
		)
		write_tsv_rows(
			stage / "archaic_summary.tsv",
			[
				"reference",
				"sample_id",
				"chromosomes",
				"n_chromosomes",
				"tested_bp",
				"haploid_bp",
				"archaic_bp",
				"coverage_pct",
				"coordinate_coverage_pct",
			],
			reference_rows,
		)
		filters = []
		if con.execute("SELECT 1 FROM sqlite_master WHERE name='ibdmix_filter_runs'").fetchone():
			filters = [
				dict(chrom=c, status=s, daf_status=d)
				for c, s, d in con.execute(
					"SELECT chr,status,daf_status FROM ibdmix_filter_runs WHERE dataset_id=? AND genome_build=?",
					(dataset, build),
				)
			]
		x_population = set().union(*(v.get("X", set()) for v in reference_targets.values()))
		aggregates = population_burdens(reference_rows, metadata, dataset, build, "autosomes", set(samples), filters)
		aggregates += population_burdens(x_rows, metadata, dataset, build, "male_nonpar_X", x_population, filters)
		paired = paired_male_burdens(x_rows, reference_rows, metadata, dataset, build)
		write_tsv_rows(stage / "population_summary.tsv", list(aggregates[0]) if aggregates else ["dataset_id", "scope"], aggregates)
		write_tsv_rows(stage / "same_males_X_autosomes.tsv", list(paired[0]) if paired else ["dataset_id", "reference", "n_same_males"], paired)
		(stage / "manifest.json").write_text(
			json.dumps(
				dict(
					signature,
					dataset=dataset,
					build=build,
					n_samples=len(samples),
					n_bins=len(bins),
					union_intervals=unions,
					x_male_only=x_male,
					profiles=profiles,
					diploid_autosome_bp=2 * sum(lengths.get(str(ch), 0) for ch in range(1, 23)),
					lineages=["Neanderthal", "Denisovan"],
					definition="Per-lineage union bp / diploid bin span (male non-PAR X: haploid); percent; untested entries are missing",
					coordinate_coverage_definition="Exact per-individual interval union / haploid physical span; available for both autosomes and X; not ancestry dosage",
					archaic_summary_references=list(summary_references(reference_targets)),
					x_summary_definition="Male non-PAR X only; exact interval union / haploid physical non-PAR span; untested samples missing; experimental calls; chromosome-specific reference and background filters",
					archaic_summary_definition="Per-individual interval union within each explicitly tested reference on certified autosomes; overlapping bp counted once; untested entries are missing",
					ibdmix_filters=filters,
				),
				indent=2,
			)
		)
		os.replace(stage, dest)
		write_tsv_rows(pointer, ["directory"], [{"directory": version}])
		print(f"DENSITY ready {dest}", flush=True)
	con.close()
	if local_database is not None:
		local_database.unlink(missing_ok=True)


def introgression_density_cli():
	ap = argparse.ArgumentParser()
	ap.add_argument("--database", required=True, type=Path)
	ap.add_argument("--output", type=Path)
	ap.add_argument("--sample-panel", type=Path)
	ap.add_argument("--bin-bp", type=int, default=5000000)
	args = ap.parse_args()
	if args.bin_bp <= 0:
		ap.error("--bin-bp must be positive")
	prepare(args.database, args.output or args.database.parent / "normalize/density", args.sample_panel, args.bin_bp)


SUBCOMMANDS = {
	"density": introgression_density_cli,
}


def main():
	import sys

	if len(sys.argv) > 1 and sys.argv[1] in SUBCOMMANDS:
		command = sys.argv.pop(1)
		return SUBCOMMANDS[command]()
	return normalize_results_cli()


if __name__ == "__main__":
	main()
