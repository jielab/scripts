#!/usr/bin/env python3
"""Inspect UKB PGENs and prepare a bounded, explicitly QC'd GU pilot.

Preparation exchange files stay under /tmp. Verified scientific QC and sample
metadata are attached to the native method archive. No genotype imputation or
phasing is performed here. IBDmix calling code is not modified.
"""

from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager, nullcontext
import csv
from dataclasses import dataclass
import fcntl
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile

sys.dont_write_bytecode = True


# 🚩 Input and sample contracts

def fail(message):
	raise ValueError(message)


def boolean(value):
	if isinstance(value, bool):
		return value
	if str(value).lower() in {"1", "true", "yes"}:
		return True
	if str(value).lower() in {"0", "false", "no"}:
		return False
	raise argparse.ArgumentTypeError("expected true or false")


def chromosome(value):
	value = re.sub(r"^chr", "", str(value), flags = re.I).upper()
	return "X" if value == "23" else value


def chromosomes(value):
	result = []
	for item in str(value).replace(",", " ").split():
		if re.fullmatch(r"\d+-\d+", item):
			start, end = map(int, item.split("-"))
			if not 1 <= start <= end <= 22:
				fail("chromosome ranges must be within 1-22")
			items = map(str, range(start, end + 1))
		else:
			items = [chromosome(item)]
		for chrom in items:
			if chrom not in {str(x) for x in range(1, 23)} | {"X"}:
				fail(f"unsupported chromosome: {chrom}")
			if chrom not in result:
				result.append(chrom)
	if not result:
		fail("at least one chromosome is required")
	return result


def nonpar(chrom, pos, build):
	if chromosome(chrom) != "X":
		return False
	par1, par2, length = ((60001, 2699520), (154931044, 155260560), 155270560) if build == 37 else ((10001, 2781479), (155701383, 156030895), 156040895)
	return 1 <= pos <= length and not (par1[0] <= pos <= par1[1] or par2[0] <= pos <= par2[1])


def temporary_path(value):
	path = Path(value).expanduser().resolve()
	if Path("/tmp") not in path.parents:
		fail(f"UKB exchange files must be under /tmp (repository convention): {path}")
	path.mkdir(parents = True, exist_ok = True)
	return path


def permanent_path(value):
	path = Path(value).expanduser().resolve()
	for root in map(Path, ["/tmp", "/var/tmp", "/dev/shm", "/run"]):
		if path == root or root in path.parents:
			fail(f"scientific results need a permanent directory: {path}")
	legacy = Path("/mnt/d/analysis/gu")
	if path == legacy:
		fail("UKB results require a dataset subdirectory, e.g. /mnt/d/analysis/gu/ukb003")
	return path


@contextmanager
def open_text(path):
	path = Path(path)
	if path.suffix == ".zst":
		command = shutil.which("zstd")
		if not command:
			fail("zstd is required to read .pvar.zst")
		process = subprocess.Popen([command, "-dc", str(path)], stdout = subprocess.PIPE, text = True)
		try:
			yield process.stdout
		finally:
			process.stdout.close()
			code = process.wait()
			if code:
				fail(f"could not decompress {path}: exit {code}")
	else:
		with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as handle:
			yield handle


def pfile(root, chrom):
	base = Path(root) / f"chr{chrom}"
	if chrom == "X" and not base.with_suffix(".pgen").is_file():
		base = Path(root) / "chrX.male"
	paths = [Path(str(base) + suffix) for suffix in [".pgen", ".psam", ".pvar"]]
	if not paths[2].is_file():
		paths[2] = Path(str(base) + ".pvar.zst")
	for path in paths:
		if not path.is_file() or not path.stat().st_size:
			fail(f"missing or empty PGEN input: {path}")
	return base, paths


def sex_value(value):
	value = str(value).lower()
	known = {"1": "male", "m": "male", "male": "male", "2": "female", "f": "female", "female": "female"}
	if value in known:
		return known[value]
	if value in {"", "0", "na", "nan", ".", "-9", "n", "u", "unknown", "none"}:
		return "unknown"
	fail("unrecognized sex code in sample metadata")


def read_psam(path):
	rows = []
	seen = set()
	with Path(path).open() as handle:
		header = None
		for number, line in enumerate(handle, 1):
			if not line.strip() or line.startswith("##"):
				continue
			fields = line.split()
			if header is None:
				header = [x.lstrip("#").upper() for x in fields]
				if "IID" not in header:
					fail(f"PSAM has no IID header: {path}")
				continue
			if len(fields) != len(header):
				fail(f"malformed PSAM row {number}: {path}")
			data = dict(zip(header, fields))
			iid = data["IID"]
			if iid in {"", ".", "0", "NA"} or iid in seen:
				fail(f"PSAM IID must be nonmissing and unique for VCF export (row {number}): {path}")
			seen.add(iid)
			rows.append({"sample": iid, "pop": "UKB", "super_pop": "ALL", "sex": sex_value(data.get("SEX")), "fid": data.get("FID", "0")})
	if not rows:
		fail(f"no samples in {path}")
	return rows


def read_ids(path):
	result = []
	seen = set()
	column = None
	with open_text(path) as handle:
		for number, line in enumerate(handle, 1):
			if not line.strip() or line.startswith("##"):
				continue
			fields = line.split()
			if column is None:
				header = [x.lstrip("#").upper() for x in fields]
				id_headers = {"IID", "SAMPLE", "ID", "SNP"}
				if any(x in id_headers for x in header):
					column = next(i for i, x in enumerate(header) if x in id_headers)
					continue
				if line.startswith("#"):
					continue
				if len(fields) not in {1, 2}:
					fail(f"use one ID column, or a header with IID, in {path}")
				column = len(fields) - 1
			if line.startswith("#"):
				continue
			if column >= len(fields):
				fail(f"missing ID column at {path}:{number}")
			iid = fields[column]
			if iid in seen:
				fail(f"duplicate ID at {path}:{number}")
			seen.add(iid)
			result.append(iid)
	if not result:
		fail(f"no IDs in {path}")
	return result


def attach_panel(samples, path):
	if path is None:
		return samples
	with Path(path).open() as handle:
		header = handle.readline().strip().split()
		header = [x.lstrip("#").lower() for x in header]
		key = next((x for x in ["sample", "iid", "id"] if x in header), None)
		if key is None:
			fail("sample panel needs a sample or IID column")
		panel = {}
		for line in handle:
			if not line.strip():
				continue
			values = line.split()
			if len(values) != len(header):
				fail("sample panel must be a rectangular whitespace-delimited table")
			row = dict(zip(header, values))
			if row[key] in panel:
				fail("sample panel contains duplicate IDs")
			panel[row[key]] = row
	for sample in samples:
		sample["panel_present"] = sample["sample"] in panel
		if sample["sample"] not in panel:
			continue
		row = panel[sample["sample"]]
		for field in ["pop", "super_pop"]:
			if row.get(field) not in {None, "", ".", "NA"}:
				sample[field] = row[field]
		sex = sex_value(row.get("sex"))
		if sex != "unknown":
			if sample["sex"] not in {"unknown", sex}:
				fail("PSAM and supplied sample panel disagree on sex; reconcile the source metadata")
			sample["sex"] = sex
	return samples


def write_panel(path, samples, protected = None):
	handle = io.StringIO()
	writer = csv.DictWriter(handle, fieldnames = ["sample", "pop", "super_pop", "sex"], delimiter = "\t", lineterminator = "\n", extrasaction = "ignore")
	writer.writeheader()
	writer.writerows(sorted(samples, key = lambda row: row["sample"]))
	content = handle.getvalue()
	path = Path(path)
	if path.is_file() and path.read_text() == content:
		return
	if protected is not None and path.resolve() == Path(protected).resolve():
		fail("the prepared sample-panel path would overwrite the supplied panel; use a separate --ukb-work directory")
	with tempfile.NamedTemporaryFile("w", dir = path.parent, prefix = ".samples-", delete = False) as handle:
		handle.write(content)
		stage = Path(handle.name)
	stage.replace(path)


def sample_hash(samples):
	return hashlib.sha256("".join(x["sample"] + "\n" for x in samples).encode()).hexdigest()


def cohort_hash(samples):
	fields = ["sample", "pop", "super_pop", "sex"]
	text = "".join("\t".join(x[key] for key in fields) + "\n" for x in sorted(samples, key = lambda x: x["sample"]))
	return hashlib.sha256(text.encode()).hexdigest()


def validate_cohort(samples, args):
	if not samples:
		fail("the selected cohort contains no samples")
	if args.sample_panel and any(not x.get("panel_present") for x in samples):
		fail("the supplied sample panel is missing selected individuals; provide complete metadata for the keep set")
	labels = {(x["pop"], x["super_pop"]) for x in samples}
	if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", label) for pair in labels for label in pair):
		fail("population labels may contain only letters, numbers, '.', '_' and '-'")
	if len(labels) != 1:
		fail("the pilot must use one fixed population cohort; select --ukb-pop or a population-specific --ukb-keep instead of splitting IBDmix into sample batches")
	pop, super_pop = next(iter(labels))
	return {"population": pop, "super_population": super_pop, "samples": len(samples), "metadata_status": "unlabelled_technical_pilot" if (pop, super_pop) == ("UKB", "ALL") else "supplied_population_labels", "sample_metadata_sha256": cohort_hash(samples)}


def compare_sample_metadata(expected, observed):
	by_id = {row["sample"]: row for row in expected}
	for row in observed:
		other = by_id.get(row["sample"])
		if other is None or any(row[field] != other[field] for field in ["sex", "fid", "pop", "super_pop"]):
			fail("selected IID/FID/sex/population metadata differ between chromosome PSAMs; reconcile the source files before export")


def metadata_chromosome(args):
	return next((x for x in args.chroms if x != "X"), "1" if (args.root / "chr1.psam").is_file() else "X")


def metadata_psam(args):
	chrom = metadata_chromosome(args)
	path = args.root / f"chr{chrom}.psam"
	return path if path.is_file() else pfile(args.root, chrom)[1][1]


# 🚩 Variant quality: dense data, exact allele matching, no LD pruning

@dataclass
class Variant:
	chrom: str
	pos: int
	id: str
	ref: str
	alt: str
	info: str
	ordinal: int = 0


def read_pvar(path, indexed_copy=None):
	# A PGEN indexes variants by row. Rewrite only IDs in a disposable PVAR,
	# preserving every row, allele and header so duplicate rsIDs cannot make
	# --extract include the wrong variant. The source files remain untouched.
	with open_text(path) as handle, (Path(indexed_copy).open('w') if indexed_copy else nullcontext()) as rewritten:
		header = None
		ordinal = 0
		for number, line in enumerate(handle, 1):
			if not line.strip() or line.startswith("##"):
				if rewritten: rewritten.write(line)
				continue
			if line.startswith("#"):
				header = [x.lstrip("#").upper() for x in line.split()]
				if header[:5] != ["CHROM", "POS", "ID", "REF", "ALT"]:
					fail(f"invalid PVAR header: {path}")
				if rewritten: rewritten.write(line)
				continue
			fields = line.split()
			if header is None or len(fields) != len(header):
				fail(f"malformed PVAR row {number}: {path}")
			try:
				pos = int(fields[1])
			except ValueError:
				fail(f"invalid PVAR position at row {number}: {path}")
			ordinal += 1
			row = Variant(chromosome(fields[0]), pos, fields[2], fields[3].upper(), fields[4].upper(), fields[header.index("INFO")] if "INFO" in header else ".", ordinal)
			if rewritten:
				fields[2] = f'gu_row_{ordinal}'
				rewritten.write('\t'.join(fields) + '\n')
			yield row


def biallelic_snp(row):
	return row.ref in {"A", "C", "G", "T"} and row.alt in {"A", "C", "G", "T"} and row.ref != row.alt


def quality_number(value):
	try:
		value = float(value)
	except (TypeError, ValueError):
		return None
	return value if math.isfinite(value) and 0 <= value <= 1 else None


def info_quality(row, key):
	for token in row.info.split(";"):
		name, _, value = token.partition("=")
		if name == key:
			return quality_number(value)
	return None


def quality_rows(path, chrom):
	"""Read sorted UKB 8-column MFI or CHROM POS REF ALT INFO TSV.

	MFI A1/A2 order is irrelevant to a site-level quality score; allele sets must
	still match exactly. Never match solely on position or silently complement.
	"""
	with open_text(path) as handle:
		header = None
		last_pos = -1
		for number, line in enumerate(handle, 1):
			if not line.strip() or line.startswith("##"):
				continue
			fields = line.split()
			if header is None:
				candidate = [x.lstrip("#").upper() for x in fields]
				if {"CHROM", "POS", "REF", "ALT", "INFO"} <= set(candidate):
					header = {x: candidate.index(x) for x in ["CHROM", "POS", "REF", "ALT", "INFO"]}
					continue
				if line.startswith("#"):
					continue
				if len(fields) != 8:
					fail(f"quality file needs UKB 8-column MFI or CHROM POS REF ALT INFO header: {path}")
				header = {}
			if header:
				if len(fields) <= max(header.values()):
					fail(f"malformed quality row {number}: {path}")
				if chromosome(fields[header["CHROM"]]) != chrom:
					continue
				pos, ref, alt, value = (fields[header[x]] for x in ["POS", "REF", "ALT", "INFO"])
			else:
				if len(fields) != 8:
					fail(f"malformed MFI row {number}: {path}")
				pos, ref, alt, value = fields[2], fields[3], fields[4], fields[7]
			try:
				pos = int(pos)
			except ValueError:
				fail(f"invalid quality-file position at {path}:{number}")
			if pos < 1:
				fail(f"quality-file positions must be positive: {path}:{number}")
			if pos < last_pos:
				fail(f"quality file must be sorted by position for streaming lookup: {path}")
			last_pos = pos
			yield pos, tuple(sorted([ref.upper(), alt.upper()])), quality_number(value)


class QualityStream:
	def __init__(self, path, chrom):
		self.rows = iter(quality_rows(path, chrom))
		self.next = next(self.rows, None)
		self.pos = None
		self.values = {}

	def get(self, row):
		if self.pos != row.pos:
			while self.next is not None and self.next[0] < row.pos:
				self.next = next(self.rows, None)
			self.pos, self.values = row.pos, {}
			while self.next is not None and self.next[0] == row.pos:
				_, key, value = self.next
				if key in self.values and self.values[key] != value:
					fail(f"conflicting quality scores at {row.chrom}:{row.pos}")
				self.values[key] = value
				self.next = next(self.rows, None)
		return self.values.get(tuple(sorted([row.ref, row.alt])))


def formatted_path(value, chrom):
	return Path(value.replace("{chr}", chrom)).expanduser().resolve() if value else None


def select_variants(path, chrom, args, output, indexed_pvar=None):
	quality_path = formatted_path(args.ukb_info_file, chrom)
	quality = QualityStream(quality_path, chrom) if quality_path else None
	pass_path = formatted_path(args.ukb_qc_pass_variants, chrom)
	pass_ids = set(read_ids(pass_path)) if pass_path else None
	regions = []
	if args.regions:
		regions = [(int(f[1]), int(f[2])) for f in (line.split() for line in args.regions.read_text().splitlines() if line.strip() and not line.startswith("#")) if chromosome(f[0]) == chrom]
	counts = Counter()
	seen = set()
	last_pos = -1
	with Path(output).open("w") as handle:
		pending = []
		pending_biallelic = 0
		def flush_position():
			nonlocal pending_biallelic
			if pending_biallelic > 1:
				# IBDmix joins on position. Choosing one of multiple records at
				# the same position can silently change the observed allele;
				# a low-INFO second ALT does not prove that allele is absent.
				counts["ambiguous_position_variants_excluded"] += len(pending)
				counts["ambiguous_positions_excluded"] += 1
			elif len(pending) == 1:
				handle.write(pending[0] + "\n")
				counts["selected_variants"] += 1
			pending.clear()
			pending_biallelic = 0
		for row in read_pvar(path, indexed_pvar):
			counts["input_variants"] += 1
			if args.regions and not any(lo < row.pos <= hi for lo,hi in regions):
				continue
			if not args.existing_pgen and row.id not in {".", "", "NA"}:
				if row.id in seen:
					fail("duplicate PVAR IDs prevent unambiguous --extract selection; make unique source IDs first")
				seen.add(row.id)
			if row.chrom != chrom or (chrom == "X" and not nonpar(row.chrom, row.pos, args.grch)):
				counts["outside_chromosome_or_nonpar"] += 1
				continue
			if row.pos < last_pos:
				fail("PVAR must be sorted by coordinate before quality filtering")
			if row.pos < 1:
				fail("PVAR positions must be positive")
			if row.pos != last_pos:
				flush_position()
			last_pos = row.pos
			if not biallelic_snp(row):
				counts["not_biallelic_snp"] += 1
				continue
			pending_biallelic += 1
			if row.id in {".", "", "NA"}:
				counts["missing_variant_id"] += 1
				continue
			if pass_ids is not None and row.id not in pass_ids:
				counts["not_in_qc_pass_list"] += 1
				continue
			if not args.existing_pgen and args.ukb_source == "imp" and (pass_ids is None or quality is not None):
				value = quality.get(row) if quality is not None else info_quality(row, args.ukb_info_field)
				if value is None:
					counts["quality_missing_or_invalid"] += 1
					continue
				if value < args.ukb_info_min:
					counts["quality_below_threshold"] += 1
					continue
			pending.append(f'gu_row_{row.ordinal}' if args.existing_pgen else row.id)
		flush_position()
	if not counts["selected_variants"]:
		reason = "no unique biallelic SNPs in the requested region" if args.existing_pgen else "imp requires INFO in PVAR, --ukb-info-file (MFI/TSV), or an explicit --ukb-qc-pass-variants list"
		fail("no eligible SNPs: " + reason + "; QC=" + json.dumps(counts))
	return {**dict(counts), 'extraction_id_scheme': 'original_PVAR_row_1based' if args.existing_pgen else 'source_variant_ID'}


# 🚩 Tool execution and genotype/phase audit

def run(command, log = None):
	command = list(map(str, command))
	print("Running: " + shlex.join(command), flush = True)
	result = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, text = True)
	if log:
		Path(log).write_text(result.stdout)
	if result.returncode:
		fail(f"command failed (exit {result.returncode}): {shlex.join(command)}\n{result.stdout[-6000:]}")
	return result.stdout


def plink_input(base, paths):
	return ["--pfile", str(base)] + (["vzs"] if str(paths[2]).endswith(".zst") else [])


def pgen_information(base, paths, args, out):
	raw = run(["plink2", *plink_input(base, paths), "--pgen-info", "--threads", str(args.ukb_threads), "--memory", str(args.ukb_memory_mb), "--out", out], str(out) + ".txt")
	lower = raw.lower()
	phase = False if "no hardcalls are explicitly phased" in lower or "no phased hardcalls" in lower else (True if re.search(r"phased hardcalls[^\n]*(?:present|yes)", lower) else None)
	dosage = False if re.search(r"no dosages?|dosages?[^\n]*(?:absent|not present)", lower) else (True if re.search(r"dosages?[^\n]*(?:present|yes)", lower) else None)
	phase_dosage = True if "explicitly phased dosages present" in lower else (False if dosage is False or "none explicitly phased" in lower or "unphased dosages present" in lower else None)
	return {"phased_hardcalls_present": phase, "dosages_present": dosage, "phased_dosages_present": phase_dosage}


def audit_vcf(path, expected_samples, chrom, require_phase = False, max_records = None, build = None):
	if not expected_samples:
		fail("a genotype audit requires at least one explicitly selected sample")
	counts = Counter()
	samples = None
	last_pos = 0
	site_hash = hashlib.sha256()
	max_site_missing = 0
	with open_text(path) as handle:
		for line in handle:
			if line.startswith("##"):
				continue
			if line.startswith("#CHROM"):
				if samples is not None:
					fail("VCF contains more than one sample header")
				samples = line.rstrip("\n").split("\t")[9:]
				if samples != [x["sample"] for x in expected_samples]:
					fail("VCF sample IDs/order differ from the requested PSAM subset")
				continue
			if samples is None:
				fail("VCF lacks a sample header")
			fields = line.rstrip("\n").split("\t")
			if len(fields) != 9 + len(samples) or chromosome(fields[0]) != chrom:
				fail("VCF chromosome or sample width is inconsistent")
			try:
				pos = int(fields[1])
			except ValueError:
				fail("VCF has an invalid position")
			if pos <= last_pos:
				fail("VCF positions must be positive, unique and increasing for IBDmix")
			last_pos = pos
			ref, alt = fields[3:5]
			if ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"} or ref == alt:
				fail("VCF contains a non-biallelic SNP record")
			if chrom == "X" and build is not None and not nonpar(chrom, pos, build):
				fail("chrX export contains a PAR or out-of-range coordinate")
			site_hash.update(f"{chrom}\t{pos}\t{ref}\t{alt}\n".encode())
			formats = fields[8].split(":")
			if formats.count("GT") != 1:
				fail("VCF has no GT; dosage-only output cannot be passed to GU callers")
			index = formats.index("GT")
			if fields[8] == "GT":
				calls = Counter(fields[9:])
			else:
				# Trailing missing FORMAT values are legal VCF. A missing GT
				# remains missing; it must never be interpreted as reference.
				values = (x.split(":") for x in fields[9:])
				calls = Counter(x[index] if len(x) > index else "." for x in values)
			counts["records_checked"] += 1
			missing_here = 0
			for gt, number in calls.items():
				alleles = re.split(r"[/|]", gt)
				if any(x not in {"0", "1", "."} for x in alleles) or len(alleles) not in {1, 2}:
					fail("VCF contains a non-biallelic or invalid genotype")
				counts["genotypes_checked"] += number
				if any(x != "." for x in alleles):
					if chrom == "X" and len(alleles) != 1:
						fail("male non-PAR chrX export contains diploid GT; do not duplicate the male X allele")
					if chrom != "X" and len(alleles) != 2:
						fail("autosomal export contains a non-diploid GT")
				if "." in alleles:
					counts["missing_genotypes"] += number
					missing_here += number
					if any(x != "." for x in alleles):
						counts["partially_missing_genotypes"] += number
					continue
				counts["nonmissing_genotypes"] += number
				if len(alleles) == 2 and alleles[0] != alleles[1]:
					counts["heterozygous_genotypes"] += number
					if "/" in gt:
						counts["unphased_heterozygous_genotypes"] += number
						if require_phase:
							fail("unphased heterozygous GT found: AS3/TRACE require genuinely phased dense data; never replace / with |")
			max_site_missing = max(max_site_missing, missing_here)
			if max_records is not None and counts["records_checked"] >= max_records:
				break
	if not counts["records_checked"]:
		fail("VCF contains no records")
	if require_phase and chrom != "X" and not counts["heterozygous_genotypes"]:
		fail("phase audit is inconclusive: no heterozygous GTs in the export")
	phase_status = "haploid_male_X" if chrom == "X" else ("unphased_heterozygotes_observed" if counts["unphased_heterozygous_genotypes"] else ("all_observed_heterozygotes_phased" if counts["heterozygous_genotypes"] else "inconclusive_no_heterozygotes"))
	return {**dict(counts), "phase_audit_scope": "full_export" if max_records is None else f"first_{max_records}_records", "phase_status": phase_status, "phase_required": require_phase, "postqc_sites_sha256": site_hash.hexdigest(), "missing_fraction_in_audit": counts["missing_genotypes"] / counts["genotypes_checked"], "max_site_missing_fraction": max_site_missing / len(expected_samples)}


def sample_missingness(path, samples, variants, threshold):
	"""Read PLINK statistics produced after variant QC, with explicit ID/count checks."""
	with Path(path).open() as handle:
		header = [x.lstrip("#") for x in handle.readline().split()]
		missing_key = next((x for x in ["MISSING_CT", "N_MISS"] if x in header), None)
		observed_key = next((x for x in ["OBS_CT", "N_GENO"] if x in header), None)
		if "IID" not in header or "F_MISS" not in header or not missing_key or not observed_key:
			fail("PLINK sample missingness report lacks IID, count or F_MISS columns")
		records = []
		for line in handle:
			if not line.strip():
				continue
			fields = line.split()
			if len(fields) != len(header):
				fail("malformed PLINK sample missingness report")
			records.append(dict(zip(header, fields)))
	if len(records) != len(samples) or [row["IID"] for row in records] != [row["sample"] for row in samples]:
		fail("sample missingness report does not describe the exact exported sample order")
	fractions, missing_total = [], 0
	for row, sample in zip(records, samples):
		if "FID" in row and row["FID"] != sample["fid"]:
			fail("sample missingness FID differs from the selected PSAM")
		missing, observed = int(row[missing_key]), int(row[observed_key])
		value = float(row["F_MISS"])
		if observed != variants or not 0 <= missing <= observed or not math.isfinite(value) or not 0 <= value <= 1 or abs(value - missing / observed) > 1e-5:
			fail("sample missingness statistics do not describe the final exported variants")
		fraction = missing / observed
		if fraction > threshold:
			fail("export has excessive sample hardcall missingness; inspect the imputation/hardcall threshold and eligible sample QC")
		fractions.append(fraction)
		missing_total += missing
	return {"mean_sample_missingness": sum(fractions) / len(fractions), "max_sample_missingness": max(fractions), "sample_missing_genotypes": missing_total}


def phase_certified(audit, chrom):
	if audit.get("phase_audit_scope") != "full_export" or not audit.get("nonmissing_genotypes"):
		return False
	if chrom == "X":
		return audit.get("phase_status") == "haploid_male_X"
	return audit.get("phase_status") == "all_observed_heterozygotes_phased" and audit.get("heterozygous_genotypes", 0) > 0 and audit.get("unphased_heterozygous_genotypes", 0) == 0


def selected_samples(samples, keep, chrom):
	by_id = {x["sample"]: x for x in samples}
	if keep is not None:
		missing = set(keep) - by_id.keys()
		if missing:
			fail(f"{len(missing)} requested IDs are absent from the chromosome PSAM; harmonize keep sets by IID")
		wanted = set(keep)
		samples = [x for x in samples if x["sample"] in wanted]
	if chrom == "X":
		if keep is not None and any(x["sex"] != "male" for x in samples):
			fail("the fixed male chrX keep set conflicts with chromosome PSAM/panel sex metadata")
		samples = [x for x in samples if x["sex"] == "male"]
	if not samples:
		fail("no eligible samples; chrX requires explicitly male PSAM/panel metadata")
	return samples


def write_keep_and_sex(prefix, samples):
	keep, sex = Path(str(prefix) + ".keep.txt"), Path(str(prefix) + ".sex.txt")
	keep.write_text("#IID\n" + "".join(x["sample"] + "\n" for x in samples))
	sex.write_text("#IID\tSEX\n" + "".join(x["sample"] + "\t" + {"male": "1", "female": "2", "unknown": "0"}[x["sex"]] + "\n" for x in samples))
	return keep, sex


# 🚩 Inspect and define a reproducible pilot

def inspect(args):
	probe_dir = args.work / "probe"
	probe_dir.mkdir(exist_ok = True)
	reports = []
	first_hash = None
	for chrom in args.chroms:
		base, paths = pfile(args.root, chrom)
		samples = attach_panel(read_psam(paths[1]), args.sample_panel)
		stats, buckets, probes = Counter(), set(), []
		for row in read_pvar(paths[2]):
			stats["variants"] += 1
			if row.chrom != chrom or (chrom == "X" and not nonpar(row.chrom, row.pos, args.grch)) or not biallelic_snp(row):
				continue
			stats["biallelic_snps_in_scope"] += 1
			if info_quality(row, args.ukb_info_field) is not None:
				stats["snps_with_quality"] += 1
			bucket = row.pos // 250000
			if bucket not in buckets and row.id not in {".", "NA"}:
				buckets.add(bucket)
				probes.append(row.id)
		hash_value = sample_hash(samples)
		if first_hash is None:
			first_hash = hash_value
		report = {"chromosome": chrom, "source": args.ukb_source, "samples": len(samples), **dict(stats), "sample_order_matches_first_chr": hash_value == first_hash, "sample_order_sha256": hash_value, **{f"sex_{key}": value for key, value in Counter(x["sex"] for x in samples).items()}}
		if not args.ukb_no_probe:
			report.update(pgen_information(base, paths, args, probe_dir / f"chr{chrom}.pgen_info"))
			eligible = [x for x in samples if chrom != "X" or x["sex"] == "male"]
			if eligible and probes:
				eligible = sorted(eligible, key = lambda x: hashlib.sha256((str(args.ukb_seed) + x["sample"]).encode()).digest())[:64]
				wanted = {x["sample"] for x in eligible}
				eligible = [x for x in samples if x["sample"] in wanted]
				prefix = probe_dir / f"chr{chrom}"
				keep, sex = write_keep_and_sex(prefix, eligible)
				extract = Path(str(prefix) + ".sites.txt")
				extract.write_text("\n".join(probes) + "\n")
				run(["plink2", *plink_input(base, paths), "--keep", keep, "--update-sex", sex, "--extract", extract, "--snps-only", "just-acgt", "--max-alleles", "2", "--threads", args.ukb_threads, "--memory", args.ukb_memory_mb, "--export", "vcf", "bgz", "id-paste=iid", "--out", prefix], str(prefix) + ".command.txt")
				report.update(audit_vcf(Path(str(prefix) + ".vcf.gz"), eligible, chrom, build = args.grch))
				report["phase_audit_scope"] = "spatial_probe_up_to_64_samples_one_SNP_per_250kb"
			else:
				report["phase_audit_scope"] = "not_run_no_eligible_samples_or_sites"
		reports.append(report)
		print(json.dumps(report, ensure_ascii = False), flush = True)
	(args.work / "inspection.json").write_text(json.dumps(reports, indent = 2) + "\n")
	print(f"Inspection saved under {args.work}; probe results do not certify full phasing.")


def panel_or_pilot(args):
	samples = attach_panel(read_psam(metadata_psam(args)), args.sample_panel)
	if args.ukb_keep:
		samples = selected_samples(samples, read_ids(args.ukb_keep), "1")
	if args.ukb_group or args.ukb_pop:
		if args.sample_panel is None:
			fail("--ukb-group/--ukb-pop requires an ancestry-labelled --sample-panel; UKB is not automatically EUR")
		samples = [x for x in samples if (not args.ukb_group or x["super_pop"] == args.ukb_group) and (not args.ukb_pop or x["pop"] == args.ukb_pop)]
	if not samples:
		fail("no eligible samples after keep/group selection")
	if args.action == "pilot":
		samples = sorted(samples, key = lambda x: hashlib.sha256((str(args.ukb_seed) + ":" + x["sample"]).encode()).digest())[:args.ukb_pilot_n]
		cohort = validate_cohort(samples, args)
		for chrom in args.chroms:
			_, other = pfile(args.root, chrom)
			chrom_keep = [x["sample"] for x in samples if chrom != "X" or x["sex"] == "male"]
			observed = selected_samples(attach_panel(read_psam(other[1]), args.sample_panel), chrom_keep, chrom)
			validate_cohort(observed, args)
			compare_sample_metadata(samples, observed)
		prefix = args.work / "pilot"
		keep, _ = write_keep_and_sex(prefix, samples)
		path = args.work / "pilot.samples.tsv"
		write_panel(path, samples, protected = args.sample_panel)
		print(f"Pilot: {json.dumps(cohort)}; keep={keep}; sample_panel={path}")
		print("Use the same keep set on every chromosome. An unlabelled pilot is UKB/ALL, not an ancestry-homogeneous cohort.")
	else:
		path = args.work / "samples.txt"
		write_panel(path, samples, protected = args.sample_panel)
		print(f"Sample panel: n={len(samples)}; {path}")


# 🚩 Bounded export, preserving real genotype phase

def fasta_output_code(path, build, chrom):
	index = Path(str(path) + ".fai")
	if not index.is_file():
		run(["samtools", "faidx", path])
	lengths = {}
	with index.open() as handle:
		for line in handle:
			fields = line.split()
			lengths[fields[0]] = int(fields[1])
	name = next((x for x in [chrom, "chr" + chrom] if x in lengths), None)
	if name is None:
		fail(f"FASTA lacks chromosome {chrom}")
	sentinels = {37: {"1": 249250621, "22": 51304566, "X": 155270560}, 38: {"1": 248956422, "22": 50818468, "X": 156040895}}
	for target, expected in sentinels[build].items():
		for key in [target, "chr" + target]:
			if key in lengths and lengths[key] != expected:
				fail(f"FASTA length disagrees with GRCh{build}: {key}={lengths[key]}, expected {expected}")
	return "chrMT" if name.startswith("chr") else "MT"


def signature(paths):
	paths = dict.fromkeys(Path(x).resolve() for x in paths if x is not None)
	return [{"path": str(path), "bytes": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns} for path in paths]


def export_vcf(args):
	if not args.ukb_keep:
		fail("pgen-vcf requires --ukb-keep; export of the entire UKB cohort is not a default operation")
	if not args.existing_pgen and args.ukb_source in {"hap", "typ"} and not args.ukb_sparse_control:
		fail("hap/typ are sparse array inputs; use --ukb-sparse-control true only for the density-control comparison")
	if not args.ukb_ref_fasta or not args.ukb_ref_fasta.is_file():
		fail("pgen-vcf requires a build-matched --ukb-ref-fasta")
	keep_ids = read_ids(args.ukb_keep)
	if len(keep_ids) > args.ukb_max_export_samples:
		fail(f"requested {len(keep_ids)} samples exceeds --ukb-max-export-samples {args.ukb_max_export_samples}; validate the pilot and storage plan before raising the bound")
	out = temporary_path(args.ukb_vcf_out or args.work / "vcf")
	if out != args.work and args.work not in out.parents:
		fail("--ukb-vcf-out must be inside --ukb-work so its cohort identity and preparation lock are shared")
	results_root = permanent_path(args.ukb_results_root or f"/mnt/d/analysis/gu/ukb-{args.ukb_source}-pilot")
	env = args.work / "gu-target.env"
	env.unlink(missing_ok = True)
	metadata_path = metadata_psam(args)
	all_samples = selected_samples(attach_panel(read_psam(metadata_path), args.sample_panel), keep_ids, "1")
	if args.male_only:
		all_samples = [x for x in all_samples if x["sex"] == "male"]
	if args.existing_pgen:
		psam_sexes = {x["sample"]: x["sex"] for x in read_psam(metadata_path)}
		if any(x["sex"] != psam_sexes[x["sample"]] for x in all_samples):
			fail("direct --keep-males uses PSAM sex; a panel cannot fill or override source sex")
	cohort = validate_cohort(all_samples, args)
	if (args.ukb_pop and cohort["population"] != args.ukb_pop) or (args.ukb_group and cohort["super_population"] != args.ukb_group):
		fail("the fixed --ukb-keep contains samples outside the requested population; rerun pilot with the population filter")
	# Write before taking signatures. If the user's input is a previously
	# prepared panel, an identical canonical panel keeps its original mtime.
	write_panel(args.work / "samples.txt", all_samples, protected = args.sample_panel)
	write_panel(out / "samples.txt", all_samples, protected = args.sample_panel)
	if out.parent != args.work and args.work in out.parent.parents:
		write_panel(out.parent / "samples.txt", all_samples, protected = args.sample_panel)
	reports = []
	for chrom in args.chroms:
		base, paths = pfile(args.root, chrom)
		full_samples = attach_panel(read_psam(paths[1]), args.sample_panel)
		chrom_keep = [x["sample"] for x in all_samples if x["sex"] == "male"] if chrom == "X" else [x["sample"] for x in all_samples]
		missing_x = sorted(set(chrom_keep) - {x['sample'] for x in full_samples}) if args.existing_pgen and chrom == 'X' else []
		if missing_x:
			print(f'chrX: {len(missing_x)} requested males absent from source PSAM; excluded on X only and recorded in preparation QC', flush=True)
		samples = selected_samples(full_samples, [x for x in chrom_keep if x not in set(missing_x)], chrom)
		validate_cohort(samples, args)
		compare_sample_metadata(all_samples, samples)
		write_panel(args.work / f'chr{chrom}.samples.txt', samples)
		name = f"chr{chrom}" + (".male" if chrom == "X" else "")
		final = out / (name + ".vcf.gz")
		receipt = args.work / ("." + name + ".complete.json")
		code = fasta_output_code(args.ukb_ref_fasta, args.grch, chrom)
		inputs = signature([*paths, metadata_path, args.ukb_keep, args.sample_panel, args.ukb_ref_fasta, Path(str(args.ukb_ref_fasta) + ".fai"), formatted_path(args.ukb_info_file, chrom), formatted_path(args.ukb_qc_pass_variants, chrom), args.regions, Path(__file__), shutil.which("plink2"), shutil.which("bcftools")])
		options = {key: getattr(args, key) for key in ["grch", "ukb_source", "ukb_info_min", "ukb_info_field", "ukb_hardcall_threshold", "ukb_trust_hardcalls", "ukb_require_phase", "ukb_sparse_control", "ukb_max_missing", "ukb_max_sample_missing"]}
		options["ukb_quality_source"] = ("external_INFO" if args.ukb_info_file else "PVAR_INFO") if args.ukb_source == "imp" else "sparse_array_hardcalls"
		if args.ukb_qc_pass_variants:
			options["ukb_quality_source"] = options["ukb_quality_source"] + "_and_explicit_pre_qc_IDs" if args.ukb_info_file and args.ukb_source == "imp" else "explicit_pre_qc_IDs"
		if args.existing_pgen:
			options.update(ukb_info_min=None, ukb_info_field=None, ukb_hardcall_threshold=None, ukb_trust_hardcalls=True, ukb_sparse_control=args.ukb_source != "imp", ukb_max_missing=None, ukb_max_sample_missing=None, ukb_quality_source="existing_PGEN_hardcalls;INFO_not_reassessed", existing_pgen=True)
			options.update(variant_extraction='PVAR_row_ordinal', reference_alleles='FASTA_reference_with_GT_reindexed', missing_x_policy='exclude_absent_X_only')
			if chrom == 'X': options['invalid_haploid_hardcalls'] = 'PLINK_set_invalid_haploid_missing'
		options["male_only"] = args.male_only
		if os.environ.get("GU_KEEP_PSAM_INFO"):
			options["keep_psam"] = json.loads(Path(os.environ["GU_KEEP_PSAM_INFO"]).read_text())
		contract = {"inputs": inputs, "options": options, "pgen_root": str(args.root), "sample_order_sha256": sample_hash(samples), "cohort_metadata_sha256": cohort["sample_metadata_sha256"], "scored_sample_metadata_sha256": cohort_hash(samples)}
		if missing_x:
			contract.update(missing_x_sample_ids=missing_x, chromosome_psam=str(paths[1].resolve()))
		if receipt.is_file() and final.is_file() and Path(str(final) + ".tbi").is_file():
			try:
				old = json.loads(receipt.read_text())
			except (json.JSONDecodeError, OSError):
				old = {}
			if old.get("contract") == contract and old.get("outputs") == signature([final, Path(str(final) + ".tbi")]):
				reports.append(old["qc"])
				print(f"SKIP verified UKB preparation: {final}", flush = True)
				continue
		receipt.unlink(missing_ok = True)
		with tempfile.TemporaryDirectory(prefix = ".prepare-", dir = out) as temporary:
			stage = Path(temporary)
			info = pgen_information(base, paths, args, stage / "pgen_info")
			if args.ukb_source == "imp" and info["dosages_present"] is not True and not args.ukb_trust_hardcalls and not args.existing_pgen:
				fail("imputed PGEN has no verified dosages; the original GT threshold cannot be reconstructed. Inspect the conversion, then explicitly use --ukb-trust-hardcalls true if its hardcalls are suitable")
			if args.ukb_require_phase and chrom != "X" and info["phased_hardcalls_present"] is False and info["phased_dosages_present"] is not True:
				fail("PGEN has no phased hardcalls: IBDmix can use unphased data; PhyML/AS3/TRACE require real phase; recover dense phased data before requesting phase")
			ids = stage / "selected.snps.txt"
			indexed_pvar = stage / 'source.indexed.pvar' if args.existing_pgen else None
			selection = select_variants(paths[2], chrom, args, ids, indexed_pvar)
			print(f"chr{chrom} SNP quality selection: " + json.dumps(selection), flush = True)
			keep, sex = write_keep_and_sex(stage / "selected", samples)
			filtered = stage / "filtered"
			input_args = ["--pgen", paths[0], "--pvar", indexed_pvar, "--psam", paths[1]] if indexed_pvar else plink_input(base, paths)
			command = ["plink2", *input_args, "--keep", args.ukb_keep if args.existing_pgen else keep, "--extract", ids, "--snps-only", "just-acgt", "--max-alleles", "2", "--threads", args.ukb_threads, "--memory", args.ukb_memory_mb, "--fa", args.ukb_ref_fasta, "--ref-from-fa"]
			if args.existing_pgen:
				# Imported UKB allele order may be marked known although ALT is
				# the FASTA reference. PLINK swaps allele indexes and GTs together;
				# the subsequent strict FASTA audit still rejects true mismatches.
				command += ["force"]
			if args.male_only or chrom == "X":
				command += ["--keep-males"]
			if not args.existing_pgen and args.ukb_source == "imp" and info["dosages_present"] is True:
				command += ["--hard-call-threshold", args.ukb_hardcall_threshold]
			if not args.existing_pgen:
				command += ["--update-sex", sex]
			command += ["--make-pgen", "--out", filtered]
			run(command, stage / "filter.command.txt")
			# Regenerate hardcalls first, then evaluate their missingness in a new
			# invocation. REF reassignment cannot share a PLINK stats command.
			quality_filtered = stage / "quality_filtered"
			quality_args = ([] if args.existing_pgen else ["--geno", args.ukb_max_missing])
			if args.existing_pgen and chrom == 'X':
				quality_args += ['--missing', 'sample-only', 'scols=maybefid,nmiss,nobs,hethap', '--set-invalid-haploid-missing']
			run(["plink2", "--pfile", filtered, *quality_args, "--threads", args.ukb_threads, "--memory", args.ukb_memory_mb, "--make-pgen", "--out", quality_filtered], stage / "missing_filter.command.txt")
			ploidy_qc = {}
			if args.existing_pgen and chrom == 'X':
				with Path(str(quality_filtered) + '.smiss').open() as handle:
					before = list(csv.DictReader(handle, delimiter='\t'))
				ploidy_qc = {'input_missing_genotypes': sum(int(row['MISSING_CT']) for row in before),
					'invalid_haploid_genotypes_set_missing': sum(int(row['HETHAP_CT']) for row in before)}
			# PLINK's --missing reports statistics from before --geno in the same
			# invocation. Reload the filtered file so sample QC describes exactly
			# the variants which are exported, not rejected imputation sites.
			export = stage / "export"
			run(["plink2", "--pfile", quality_filtered, "--threads", args.ukb_threads, "--memory", args.ukb_memory_mb, "--output-chr", code, "--missing", "sample-only", "--export", "vcf", "bgz", "id-paste=iid", "--out", export], stage / "export.command.txt")
			checked = stage / "checked.vcf.gz"
			run(["bcftools", "norm", "-f", args.ukb_ref_fasta, "-c", "e", "-Oz", "-o", checked, Path(str(export) + ".vcf.gz")], stage / "reference_check.txt")
			run(["bcftools", "index", "-t", checked])
			n_variants = int(run(["bcftools", "index", "-n", checked]).strip())
			if n_variants <= 0:
				fail("no variants after hardcall missingness filtering")
			missingness = sample_missingness(Path(str(export) + ".smiss"), samples, n_variants, (1.0 if args.existing_pgen else args.ukb_max_sample_missing))
			if ploidy_qc and sum(ploidy_qc.values()) != missingness['sample_missing_genotypes']:
				fail('male X missingness disagrees with input missing calls plus invalid haploid hardcalls')
			audit = audit_vcf(checked, samples, chrom, args.ukb_require_phase, build = args.grch)
			if audit["records_checked"] != n_variants or audit.get("missing_genotypes", 0) != missingness["sample_missing_genotypes"]:
				fail("full VCF genotype counts disagree with PLINK's post-QC sample missingness report")
			if not args.existing_pgen and audit["max_site_missing_fraction"] > args.ukb_max_missing:
				fail("exported genotypes violate the requested site missingness threshold")
			qc = {
				"chromosome": chrom, "source": args.ukb_source, "samples": len(samples), "variants": n_variants,
				"purpose": "existing_PGEN_subset" if args.existing_pgen else "sparse_array_control" if args.ukb_source != "imp" else "imputed_dense_pilot",
				"cohort": cohort, "quality_selection": selection, "pgen": info, "genotype_audit": audit, **missingness,
				"modern_indel_mask": "unavailable_after_SNP_only_export",
				"ibdmix_interpretation": {"status": "ukb_exploratory_not_cell2020", "reference_affinity_only": True, "strict_evidence_eligible": False},
				"hardcall_uncertainty": "DS distance is not an original genotype posterior probability" if info["dosages_present"] else "existing hardcalls explicitly accepted",
				"sample_selection": {"requested": len(chrom_keep), "retained": len(samples), "missing_from_chromosome_psam": missing_x},
				"ploidy_preparation": ploidy_qc,
			}
			checked.replace(final)
			Path(str(checked) + ".tbi").replace(Path(str(final) + ".tbi"))
			receipt.write_text(json.dumps({"contract": contract, "outputs": signature([final, Path(str(final) + ".tbi")]), "qc": qc}, indent = 2) + "\n")
			reports.append(qc)
			print(json.dumps(qc), flush = True)
	(args.work / "preparation.qc.json").write_text(json.dumps(reports, indent = 2) + "\n")
	phase_checked = all(phase_certified(report["genotype_audit"], report["chromosome"]) for report in reports)
	environment = {
		"GU_TARGET_NATIVE_VCF_PREFIX": out / "chr", "GU_TARGET_GEN_PREFIX": out / "chr",
		"GU_SAMPLE_PANEL": args.work / "samples.txt", "GU_NORMALIZE_SAMPLE_PANEL": args.work / "samples.txt",
		"GU_UKB_SCORED_SAMPLE_PANEL": args.work / 'chrX.samples.txt' if 'X' in args.chroms else '',
		"GU_TARGET_ROOT": args.work, "GU_BUILD_CHECK_INPUT": os.environ.get("GU_BUILD_CHECK_INPUT") or str(metadata_path).removesuffix(".psam") + ".pvar",
		"GU_UKB_PGEN_PREPARED": 1, "GU_UKB_PHASE_CHECKED": int(phase_checked),
		"GU_TARGET": args.dataset or f"ukb-{args.ukb_source}-pilot", "GU_BUILD": args.grch,
		"GU_ANALYSIS_ROOT": results_root, "GU_PUBLISHED_ROOT": results_root, "GU_UKB_RESULTS_ROOT": results_root,
	}
	env.write_text("\n".join("export " + key + "=" + shlex.quote(str(value)) for key, value in environment.items()) + "\n")
	print(f"Prepared {len(reports)} chromosomes. source {env}; target prefix {out / 'chr'}")
	print("IBDmix frequencies are not frozen by exporting separate batches. Use one defined cohort for the pilot; validate the production frequency strategy separately.")


# 🚩 Verify the fixed cohort and carry scientific QC into the method archive

def prepared_cohort(args):
	path = args.work / "samples.txt"
	if not path.is_file():
		fail("prepared sample metadata are missing; rerun pgen-vcf")
	samples = [{"sample": iid, "fid": "0", "sex": "unknown", "pop": "UKB", "super_pop": "ALL"} for iid in read_ids(path)]
	return attach_panel(samples, path)


def scored_cohort(cohort, contract, chrom):
	eligible = [x for x in cohort if chrom != "X" or x["sex"] == "male"]
	missing = contract.get('missing_x_sample_ids', [])
	if missing:
		options = contract.get('options', {})
		if chrom != 'X' or not options.get('existing_pgen') or options.get('missing_x_policy') != 'exclude_absent_X_only':
			fail('missing sample exclusions are only valid for direct male X preparation')
		psam = contract.get('chromosome_psam')
		if not psam or psam not in {x['path'] for x in contract.get('inputs', [])}:
			fail('chrX sample exclusions lack a fingerprinted source PSAM')
		present = {x['sample'] for x in read_psam(psam)}
		if missing != sorted({x['sample'] for x in eligible} - present):
			fail('chrX sample exclusions disagree with the source PSAM')
		missing = set(missing)
		eligible = [x for x in eligible if x['sample'] not in missing]
	return eligible


def validated_preparation(args, chrom):
	"""Validate a completed export before either calling or preparing references."""
	out = Path(args.ukb_vcf_out or args.work / "vcf").expanduser().resolve()
	name = f"chr{chrom}" + (".male" if chrom == "X" else "")
	receipt = args.work / ("." + name + ".complete.json")
	if not receipt.is_file():
		fail(f"missing preparation receipt for chr{chrom}; rerun pgen-vcf")
	saved = json.loads(receipt.read_text())
	contract = saved.get("contract", {})
	inputs = contract.get("inputs", [])
	if not inputs or signature([item["path"] for item in inputs]) != inputs:
		fail(f"chr{chrom} preparation inputs changed; rerun pgen-vcf before calling")
	outputs = [out / (name + ".vcf.gz"), out / (name + ".vcf.gz.tbi")]
	if saved.get("outputs") != signature(outputs):
		fail(f"chr{chrom} prepared VCF/index changed; rerun pgen-vcf before calling")
	options = contract.get("options", {})
	qc = saved.get("qc", {})
	if options.get("grch") != args.grch or qc.get("chromosome") != chrom or qc.get("source") != options.get("ukb_source") or not qc.get("samples") or not qc.get("variants"):
		fail(f"chr{chrom} preparation receipt is inconsistent with this run")
	cohort = prepared_cohort(args)
	cohort_digest = cohort_hash(cohort)
	if contract.get("cohort_metadata_sha256") != cohort_digest or qc.get("cohort", {}).get("sample_metadata_sha256") != cohort_digest or qc["cohort"].get("samples") != len(cohort):
		fail(f"chr{chrom} receipt no longer matches the prepared cohort metadata; use a separate --ukb-work for each fixed cohort")
	eligible = scored_cohort(cohort, contract, chrom)
	if len(eligible) != qc["samples"] or cohort_hash(eligible) != contract.get("scored_sample_metadata_sha256"):
		fail(f"chr{chrom} sample count no longer matches the prepared cohort")
	audit = qc.get("genotype_audit", {})
	if audit.get("phase_audit_scope") != "full_export" or audit.get("records_checked") != qc["variants"] or audit.get("genotypes_checked") != qc["variants"] * qc["samples"] or not re.fullmatch("[0-9a-f]{64}", audit.get("postqc_sites_sha256", "")):
		fail(f"chr{chrom} lacks a complete post-QC genotype/site audit; rerun pgen-vcf")
	if (options.get("ukb_require_phase") or args.ukb_require_phase) and not phase_certified(audit, chrom):
		fail(f"chr{chrom} receipt does not certify the requested genotype phase")
	if qc["source"] == "imp" and qc.get("pgen", {}).get("dosages_present") is not True and options.get("ukb_trust_hardcalls") is not True:
		fail(f"chr{chrom} imputation hardcall provenance is not verified")
	return saved


def _preparation_identity(payload):
	"""The scientific cohort/quality contract shared by every chromosome/method."""
	rows = payload.get("chromosomes", [])
	if not isinstance(rows, list) or not rows:
		fail("UKB archived preparation QC has no chromosome contracts")
	identity = None
	keys = ["ukb_info_min", "ukb_info_field", "ukb_hardcall_threshold", "ukb_trust_hardcalls", "ukb_sparse_control", "ukb_max_missing", "ukb_max_sample_missing", "ukb_quality_source"]
	for row in rows:
		options = row.get("quality_thresholds", {})
		cohort = row.get("cohort_metadata_sha256", "")
		root = row.get("pgen_root")
		if not re.fullmatch("[0-9a-f]{64}", cohort) or not root or not Path(root).is_absolute() or any(key not in options for key in keys):
			fail("UKB archived QC lacks a complete fixed-cohort/source/quality contract; restage verified preparation before reuse")
		if options.get("ukb_source") not in {"imp", "hap", "typ"} or options.get("grch") not in {37, 38} or row.get("genome_build") != options["grch"] or row.get("preparation", {}).get("source") != options["ukb_source"]:
			fail("UKB archived QC source/build contract is inconsistent")
		current = {"cohort_metadata_sha256": cohort, "pgen_root": str(root), "source": options["ukb_source"], "genome_build": options["grch"], "quality_thresholds": {key: options[key] for key in keys}}
		# ukb_require_phase verifies existing GTs; changing it does not change
		# the cohort, sites or hardcalls, so it is deliberately excluded.
		if identity is not None and current != identity:
			fail("UKB QC mixes cohorts, PGEN sources, builds or genotype/site quality thresholds within one archive")
		identity = current
	return identity


def archived_preparation_identity(path):
	"""Public read-only helper for finalization after /tmp has been cleared."""
	try:
		payload = json.loads(Path(path).read_text())
		if not isinstance(payload, dict):
			fail("UKB archived preparation QC must be a JSON object")
		return _preparation_identity(payload)
	except (ValueError, OSError, TypeError, AttributeError) as error:
		fail(f"cannot validate archived UKB preparation QC {path}: {error}")


def _native_archive_preparation_identity(path):
	"""Read only the QC member of a compact archive, without extraction."""
	try:
		with tarfile.open(path, "r|gz") as bundle:
			for member in bundle:
				if Path(member.name) != Path("inputs/ukb.preparation.qc.json"):
					continue
				if not member.isfile() or member.size > 32 * 1024 * 1024:
					fail("native UKB preparation QC must be a bounded regular JSON member")
				with bundle.extractfile(member) as handle:
					payload = json.load(handle)
				if not isinstance(payload, dict):
					fail("native UKB preparation QC must be a JSON object")
				return _preparation_identity(payload)
		fail("native UKB archive has no inputs/ukb.preparation.qc.json")
	except (ValueError, OSError, TypeError, AttributeError, tarfile.TarError) as error:
		fail(f"cannot validate compact UKB archive {path}: {error}")


def validate_result_cohort(root, dataset, expected = None, expected_method = None):
	"""Compare all archived/restored method QC belonging to one UKB dataset.

	Returns its one scientific identity, or None when no results exist. Compact
	*.raw.tar.gz files are streamed only as far as their QC member. PGEN/VCF
	files and temporary input directories are never needed for this check.
	"""
	if not re.fullmatch(r"ukb-[A-Za-z0-9_.-]+", str(dataset)):
		fail("fixed UKB cohort validation requires a safe ukb-* dataset name")
	root = Path(root).expanduser().resolve()
	identity = expected
	method_identities = {expected_method: expected} if expected is not None and expected_method else {}
	def compare(observed, path, method):
		nonlocal identity
		if identity is not None and observed["cohort_metadata_sha256"] != identity["cohort_metadata_sha256"]:
			fail(f"UKB dataset {dataset} already contains a different cohort, PGEN source, build or quality contract ({path}); use a separate results root/dataset for a different analysis")
		if method in method_identities and observed != method_identities[method]:
			fail(f"UKB {method} results mix PGEN sources/builds/quality contracts; use a separate results root")
		method_identities[method] = observed
		identity = observed
	for method in ["ibdmix", "trace", "as3", "phyml"]:
		for base in [root / method / dataset, root / "ukb" / method / dataset]:
			native_runs = set()
			for path in sorted(base.glob("**/inputs/ukb.preparation.qc.json")):
				try:
					observed = archived_preparation_identity(path)
				except ValueError:
					if path.exists():
						raise
					# Compaction publishes its archive before removing native
					# files. The archive pass below covers that transition.
					continue
				compare(observed, path, method)
				native_runs.add(path.parent.parent)
			for archive in sorted(base.glob("**/*.raw.tar.gz")):
				if archive.parent in native_runs:
					continue
				compare(_native_archive_preparation_identity(archive), archive, method)
	return identity


@contextmanager
def _result_cohort_lock(root):
	# Lock the results directory itself: no permanent manifest or ID list is
	# needed just to serialize competing first-chromosome preparations.
	root.mkdir(parents = True, exist_ok = True)
	fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
	try:
		fcntl.flock(fd, fcntl.LOCK_EX)
		yield
	finally:
		os.close(fd)


def stage_qc(args):
	"""Attach verified preparation QC to the native GU run for scientific archiving."""
	if not args.ukb_run_dir:
		fail("stage-qc requires --ukb-run-dir")
	rows, scored, excluded = [], [], []
	cohort = prepared_cohort(args)
	for chrom in args.chroms:
		saved = validated_preparation(args, chrom)
		contract = saved["contract"]
		scored.extend(dict(chromosome=chrom, **row) for row in scored_cohort(cohort, contract, chrom))
		excluded.extend(dict(chromosome=chrom, sample=iid, reason='absent_from_source_PSAM') for iid in contract.get('missing_x_sample_ids', []))
		rows.append({"chromosome": chrom, "genome_build": args.grch, "preparation": saved["qc"], "quality_thresholds": contract["options"], "pgen_root": contract["pgen_root"], "sample_order_sha256": contract["sample_order_sha256"], "cohort_metadata_sha256": contract["cohort_metadata_sha256"], "scored_sample_metadata_sha256": contract["scored_sample_metadata_sha256"], "input_fingerprints": contract["inputs"]})
	run_dir = Path(args.ukb_run_dir).expanduser().resolve()
	destination = run_dir / "inputs" / "ukb.preparation.qc.json"
	samples = prepared_cohort(args)
	if all(chrom == "X" for chrom in args.chroms):
		samples = [row for row in samples if row["sex"] == "male"]
	payload = {"cohort_scope": "fixed UKB pilot; IBDmix frequencies are not frozen across different sample batches", "chromosomes": rows}
	identity = _preparation_identity(payload)
	dataset = os.environ.get("GU_TARGET") or f"ukb-{identity['source']}-pilot"
	results_root = permanent_path(args.ukb_results_root or f"/mnt/d/analysis/gu/ukb-{identity['source']}-pilot")
	content = json.dumps(payload, indent = 2, sort_keys = True) + "\n"
	with _result_cohort_lock(results_root):
		validate_result_cohort(results_root, dataset, identity, run_dir.relative_to(results_root).parts[0])
		if destination.is_file() and archived_preparation_identity(destination) != identity:
			fail("the current method run already contains a different UKB cohort/quality contract; do not overwrite it with a new pilot")
		destination.parent.mkdir(parents = True, exist_ok = True)
		# Keep the requested cohort, including unavailable X samples, so final
		# can distinguish untested people from measured zero-call individuals.
		write_panel(destination.parent / "ukb.samples.tsv", samples)
		for name, records, fields in [('ukb.scored.samples.tsv', scored, ['chromosome', 'sample', 'pop', 'super_pop', 'sex']),
			('ukb.excluded.samples.tsv', excluded, ['chromosome', 'sample', 'reason'])]:
			with (destination.parent / name).open('w') as handle:
				writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
				writer.writeheader()
				writer.writerows(records)
		if not destination.is_file() or destination.read_text() != content:
			with tempfile.NamedTemporaryFile("w", dir = destination.parent, prefix = ".ukb-qc-", delete = False) as handle:
				handle.write(content)
				staged = Path(handle.name)
			staged.replace(destination)
	print(f"Verified UKB preparation QC staged for method-result archiving: {destination}")


# 🚩 UKB archaic scoring at observed modern sites

def _ukb_scoring_chrom(value):
	value = re.sub(r'^chr', '', str(value), flags=re.I).upper()
	return 'X' if value == '23' else value


def _ukb_scoring_blocks(chrom, build):
	if chrom != 'X':
		return []
	if build == 37:
		return [(0, 60000), (2699520, 154931043), (155260560, 155270560)]
	if build == 38:
		return [(0, 10000), (2781479, 155701382), (156030895, 156040895)]
	raise ValueError('UKB scoring requires GRCh37 or GRCh38')


def _ukb_scoring_source(root, ref, chrom):
	patterns = [f'*chr{chrom}_*.vcf.gz', f'*chr{chrom}.*.vcf.gz', f'*chr{chrom}.vcf.gz', f'*.{chrom}.mod.vcf.gz']
	for directory in [root / ref, root / 'avcf' / ref]:
		if not directory.is_dir():
			continue
		found = sorted({p for pattern in patterns for p in directory.glob(pattern) if p.is_file() and not Path(str(p) + '.aria2').exists()})
		if len(found) > 1:
			raise ValueError(f'Ambiguous archaic input for {ref} chr{chrom}; provide an unambiguous source directory')
		if found:
			return found[0].resolve()
	raise ValueError(f'Missing archaic VCF for {ref} chr{chrom} under {root}')


def _ukb_scoring_sites(modern, chrom, build, destination, expected_sha256=None):
	"""Stream just four fields; never materialize all modern GTs in Python."""
	blocks = [{'start': lo, 'end': hi, 'name': f'UKB_X_NP{i + 1}', 'modern_sites': 0} for i, (lo, hi) in enumerate(_ukb_scoring_blocks(chrom, build))]
	process = subprocess.Popen(['bcftools', 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\n', str(modern)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	count, previous = 0, 0
	digest = hashlib.sha256()
	try:
		with destination.open('w') as out:
			for line in process.stdout:
				fields = line.rstrip('\n').split('\t')
				if len(fields) != 4:
					raise ValueError('Malformed post-QC modern site table')
				c, position, ref, alt = fields
				c, position, ref, alt = _ukb_scoring_chrom(c), int(position), ref.upper(), alt.upper()
				if c != chrom or position <= previous:
					raise ValueError('Modern scoring sites must be one chromosome, sorted and unique by position')
				if len(ref) != 1 or len(alt) != 1 or ref not in 'ACGT' or alt not in 'ACGT' or ref == alt:
					raise ValueError('Modern scoring sites must be unambiguous biallelic SNPs')
				if blocks:
					selected = [b for b in blocks if b['start'] < position <= b['end']]
					if len(selected) != 1:
						raise ValueError('Modern X scoring site is outside male non-PAR scope')
					selected[0]['modern_sites'] += 1
				canonical = f'{c}\t{position}\t{ref}\t{alt}\n'
				out.write(canonical)
				digest.update(canonical.encode())
				previous = position
				count += 1
		stderr = process.stderr.read()
		if process.wait():
			raise ValueError('bcftools failed while reading post-QC sites: ' + stderr[-2000:])
	except BaseException:
		process.kill()
		process.wait()
		raise
	finally:
		process.stdout.close()
		process.stderr.close()
	if not count:
		raise ValueError('No post-QC modern scoring sites')
	value = digest.hexdigest()
	if expected_sha256 is not None and value != expected_sha256:
		raise ValueError('Post-QC modern site hash disagrees with the validated preparation receipt')
	for block in blocks:
		block['status'] = 'scored_in_separate_call_unit' if block['modern_sites'] else 'no_information_no_modern_sites'
	return {'modern_sites': count, 'postqc_sites_sha256': value, 'blocks': blocks}


def _ukb_scoring_site_rows(path):
	with path.open() as handle:
		for line in handle:
			c, pos, ref, alt = line.rstrip('\n').split('\t')
			yield int(pos), ref, alt


def _ukb_scoring_filter_reference(source, destination, sites, chrom, modern_site_count, blocks=()):
	"""Restrict archaic records to observed modern sites, keeping original GTs.

	Do not prefilter on archaic quality masks: upstream intentionally retains
	discordant-homozygote evidence at masked sites. Its existing masks still apply.
	Archaic ALT='.' / 0/0 is retained. A missing archaic record is never invented.
	"""
	header = subprocess.run(['bcftools', 'view', '-h', str(source)], check=True, capture_output=True, text=True).stdout
	contigs = re.findall(r'^##contig=<ID=([^,>]+)', header, flags=re.M)
	matches = [c for c in contigs if _ukb_scoring_chrom(c) == chrom]
	if len(matches) != 1:
		raise ValueError(f'Archaic input needs exactly one declared alias for chr{chrom}: {source}')
	source_contig = matches[0]
	sample_line = next((line for line in header.splitlines() if line.startswith('#CHROM')), '')
	if len(sample_line.split('\t')) != 10:
		raise ValueError(f'Archaic input must have exactly one sample: {source}')
	targets = destination.parent / ('targets.' + chrom + '.bed')
	with targets.open('w') as out:
		for pos, _, _ in _ukb_scoring_site_rows(sites):
			out.write(f'{source_contig}\t{pos - 1}\t{pos}\n')
	log = destination.parent / 'prepare.log'
	stats = {'retained_archaic_records': 0, 'allele_compatible_records': 0, 'nonmissing_allele_compatible_records': 0, 'missing_archaic_GT_records': 0, 'non_SNP_archaic_records_ignored': 0, 'incompatible_ALT_records': 0}
	stats['blocks'] = [{**b, 'allele_compatible_records': 0, 'nonmissing_allele_compatible_records': 0} for b in blocks]
	modern_rows = iter(_ukb_scoring_site_rows(sites))
	current = next(modern_rows, None)
	previous = 0
	reader = writer = None
	try:
		with log.open('w') as errors:
			reader = subprocess.Popen(['bcftools', 'view', '-T', str(targets), '--targets-overlap', '0', '-Ov', str(source)], stdout=subprocess.PIPE, stderr=errors, text=True)
			writer = subprocess.Popen(['bcftools', 'view', '-Oz', '-o', str(destination)], stdin=subprocess.PIPE, stderr=errors, text=True)
			for line in reader.stdout:
				if line.startswith('##contig='):
					match = re.match(r'##contig=<ID=([^,>]+)', line)
					if match and match.group(1) == source_contig:
						writer.stdin.write(line.replace('ID=' + source_contig, 'ID=' + chrom, 1))
					continue
				if line.startswith('#'):
					writer.stdin.write(line)
					continue
				f = line.rstrip('\n').split('\t')
				if len(f) != 10 or f[0] != source_contig:
					raise ValueError('Unexpected archaic chromosome or sample width')
				pos = int(f[1])
				if pos <= previous:
					raise ValueError('Archaic scoring records must be sorted and unique by position')
				previous = pos
				while current is not None and current[0] < pos:
					current = next(modern_rows, None)
				if current is None or current[0] != pos:
					raise ValueError('Archaic scoring record escaped the exact modern position set')
				ref, alt = f[3].upper(), f[4].upper()
				if len(ref) != 1 or ref not in 'ACGT' or (alt != '.' and (len(alt) != 1 or alt not in 'ACGT')):
					stats['non_SNP_archaic_records_ignored'] += 1
					continue
				if ref != current[1]:
					raise ValueError(f'Archaic REF disagrees with the FASTA-checked modern REF at chr{chrom}:{pos}')
				fmt = f[8].split(':')
				if 'GT' not in fmt or len(f[9].split(':')) <= fmt.index('GT'):
					raise ValueError(f'Missing archaic GT field at chr{chrom}:{pos}')
				gt = f[9].split(':')[fmt.index('GT')]
				alleles = re.split(r'[/|]', gt)
				if len(alleles) not in {1, 2} or any(a not in {'.', '0', '1'} for a in alleles) or (alt == '.' and '1' in alleles):
					raise ValueError(f'Invalid biallelic archaic GT at chr{chrom}:{pos}')
				stats['missing_archaic_GT_records'] += int('.' in alleles)
				compatible = alt in {'.', current[2]}
				nonmissing_compatible = compatible and '.' not in alleles
				stats['allele_compatible_records'] += int(compatible)
				stats['nonmissing_allele_compatible_records'] += int(nonmissing_compatible)
				for block in stats['blocks']:
					if block['start'] < pos <= block['end']:
						block['allele_compatible_records'] += int(compatible)
						block['nonmissing_allele_compatible_records'] += int(nonmissing_compatible)
				stats['incompatible_ALT_records'] += int(not compatible)
				stats['retained_archaic_records'] += 1
				f[0], f[3], f[4] = chrom, ref, alt
				writer.stdin.write('\t'.join(f) + '\n')
			reader.stdout.close()
			if reader.wait():
				raise ValueError('bcftools failed when subsetting archaic input: ' + log.read_text()[-2000:])
			writer.stdin.close()
			if writer.wait():
				raise ValueError('bcftools failed when writing scoped archaic VCF: ' + log.read_text()[-2000:])
		subprocess.run(['bcftools', 'index', '-t', str(destination)], check=True)
	except BaseException:
		for process in (reader, writer):
			if process is not None and process.poll() is None:
				process.kill()
				process.wait()
		destination.unlink(missing_ok=True)
		Path(str(destination) + '.tbi').unlink(missing_ok=True)
		raise
	finally:
		modern_rows.close()
		for process, attr in ((reader, 'stdout'), (writer, 'stdin')):
			if process is not None:
				stream = getattr(process, attr)
				if stream is not None and not stream.closed:
					try:
						stream.close()
					except BrokenPipeError:
						pass
	stats['modern_sites_without_compatible_archaic_record'] = modern_site_count - stats['allele_compatible_records']
	stats['unknown_sites_are_not_homozygous_reference'] = True
	for block in stats['blocks']:
		block['status'] = 'eligible_separate_call_unit' if block['nonmissing_allele_compatible_records'] else 'no_information_no_nonmissing_compatible_archaic_records'
	return stats


def prepare_archaic(args):
	"""CLI action: prepare-archaic; emit a receipt path per requested chromosome."""
	if not args.ukb_archaic_root:
		raise ValueError('prepare-archaic requires --ukb-archaic-root with original VCFs')
	source_root = Path(args.ukb_archaic_root).expanduser().resolve()
	mask_root = Path(args.ukb_mask_root).expanduser().resolve() if args.ukb_mask_root else source_root.parent / 'mask'
	fallback_mask_root = source_root.parent / 'mask'
	if not fallback_mask_root.is_dir():
		fallback_mask_root = None
	if not source_root.is_dir() or not mask_root.is_dir():
		raise ValueError('Original archaic root and original mask root must exist')
	canonical = {x.lower(): x for x in ('Altai', 'Chagyr', 'Vindija', 'Denisova', 'Denisova25')}
	names = str(args.ukb_archaic_refs or 'Altai Vindija').replace(',', ' ').split()
	if not names or any(x.lower() not in canonical for x in names):
		raise ValueError('UKB scoring accepts only the five current references; 2013 references are reserved for IBDmix replication')
	refs = list(dict.fromkeys(canonical[x.lower()] for x in names))
	vcf_root = Path(args.ukb_vcf_out or args.work / 'vcf').resolve()
	for chrom in args.chroms:
		prepared = validated_preparation(args, chrom)
		modern = vcf_root / (f'chr{chrom}' + ('.male' if chrom == 'X' else '') + '.vcf.gz')
		base = Path(args.work) / 'scoring' / ('chr' + chrom)
		base.mkdir(parents=True, exist_ok=True)
		if Path('/tmp') not in base.resolve().parents:
			raise ValueError('UKB scoring exchange must be under /tmp')
		with (base / '.prepare.lock').open('a') as lock:
			fcntl.flock(lock, fcntl.LOCK_EX)
			sources = {ref: _ukb_scoring_source(source_root, ref, chrom) for ref in refs}
			contract = {'schema': 'ukb-observed-sites-v1', 'chromosome': chrom, 'build': args.grch, 'references': refs,
				'original_mask_root': str(mask_root), 'original_fallback_mask_root': str(fallback_mask_root) if fallback_mask_root else None,
				'preparation': prepared['contract'], 'postqc_sites_sha256': prepared['qc']['genotype_audit']['postqc_sites_sha256'], 'inputs': signature([*sources.values(), Path(__file__)])}
			receipt = base / 'target.json'
			if receipt.is_file():
				old = json.loads(receipt.read_text())
				outputs = old.get('outputs', [])
				try:
					valid_outputs = bool(outputs) and signature([x['path'] for x in outputs]) == outputs
					if fallback_mask_root:
						link = Path(old.get('archaic_root', '')).parent / 'mask'
						valid_outputs = valid_outputs and link.is_symlink() and link.resolve() == fallback_mask_root
				except OSError:
					valid_outputs = False
				if old.get('contract') == contract and valid_outputs:
					print(receipt)
					continue
			stage = base / ('inputs-' + hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest()[:24])
			if stage.exists(): shutil.rmtree(stage)
			stage.mkdir()
			try:
				if fallback_mask_root:
					# The frozen caller also searches archaic_root.parent/mask.
					# Preserve that original fallback when an explicit mask root
					# contains only some of the requested masks.
					(stage / 'mask').symlink_to(fallback_mask_root, target_is_directory=True)
				sites = stage / 'modern.sites.tsv'
				expected_hash = prepared.get('qc', {}).get('genotype_audit', {}).get('postqc_sites_sha256')
				qc = _ukb_scoring_sites(modern, chrom, args.grch, sites, expected_hash)
				outputs = [sites]
				ref_qc = {}
				for ref, source in sources.items():
					destination = stage / 'archaic' / ref / ('chr' + chrom + '.vcf.gz')
					destination.parent.mkdir(parents=True, exist_ok=True)
					ref_qc[ref] = _ukb_scoring_filter_reference(source, destination, sites, chrom, qc['modern_sites'], qc['blocks'])
					for file in (destination, Path(str(destination) + '.tbi')):
						os.utime(file, ns=(source.stat().st_mtime_ns, source.stat().st_mtime_ns))
					if chrom != 'X' and not ref_qc[ref]['nonmissing_allele_compatible_records']:
						raise ValueError(f'No nonmissing compatible archaic records for {ref} chr{chrom}; this chromosome cannot be scored')
					outputs.extend([destination, Path(str(destination) + '.tbi')])
				bed = None
				if chrom == 'X':
					for index, block in enumerate(qc['blocks']):
						block['reference_compatible_records'] = {ref: ref_qc[ref]['blocks'][index]['nonmissing_allele_compatible_records'] for ref in refs}
						block['missing_reference_callable_sites'] = [ref for ref, count in block['reference_compatible_records'].items() if not count]
						block['status'] = ('eligible_separate_call_unit' if block['modern_sites'] and not block['missing_reference_callable_sites'] else
							'no_information_no_modern_sites' if not block['modern_sites'] else 'omitted_missing_reference_callable_sites')
					if not any(b['status'] == 'eligible_separate_call_unit' for b in qc['blocks']):
						raise ValueError('No X non-PAR block has compatible observed sites in all requested references; choose references with a common informative block')
					bed = stage / 'nonpar.bed'
					bed.write_text(''.join(f"X\t{b['start']}\t{b['end']}\t{b['name']}\n" for b in qc['blocks'] if b['status'] == 'eligible_separate_call_unit'))
					outputs.append(bed)
				result = {'contract': contract, 'archaic_root': str(stage / 'archaic'), 'mask_root': str(mask_root),
					'original_fallback_mask_root': str(fallback_mask_root) if fallback_mask_root else None,
					'nonpar_bed': str(bed) if bed else None, 'modern_vcf': str(modern), 'scope': 'UKB_postqc_observed_sites',
					'modern_indel_mask': 'unavailable_after_SNP_only_export',
					'status': 'ukb_exploratory_not_cell2020', 'reference_affinity_only': True, 'strict_evidence_eligible': False,
					'supported_nonpar_block_bp': sum(b['end'] - b['start'] for b in qc['blocks'] if b['status'] == 'eligible_separate_call_unit') if chrom == 'X' else None,
					'unknown_sites': 'Not scored; absent modern or archaic records are never filled with reference genotype',
					'x_call_state': 'Each informative non-PAR block is a separate caller unit; no-information blocks are declared, not called',
					'full_X_physical_coverage_certified': False if chrom == 'X' else None, **qc, 'reference_qc': ref_qc, 'outputs': signature(outputs)}
				fd, name = tempfile.mkstemp(prefix='.target-', suffix='.json', dir=base)
				with os.fdopen(fd, 'w') as out:
					json.dump(result, out, indent=2)
					out.write('\n')
				os.replace(name, receipt)
			except BaseException:
				shutil.rmtree(stage)
				raise
			print(receipt)


# 🚩 CLI

def parser():
	ap = argparse.ArgumentParser(description = __doc__)
	ap.add_argument("action", choices = ["inspect-pgen", "make-panel-pgen", "pilot", "pgen-vcf", "stage-qc", "prepare-archaic"])
	ap.add_argument("--chr", default = os.environ.get("GU_CHRS", "22"))
	ap.add_argument("--grch", type = int, choices = [37, 38], default = int(os.environ.get("GU_BUILD", "37")))
	ap.add_argument("--sample-panel", type = Path, default = os.environ.get("GU_SAMPLE_PANEL") or None)
	ap.add_argument("--ukb-source", choices = ["imp", "hap", "typ"], default = os.environ.get("UKB_SOURCE", "imp"))
	ap.add_argument("--ukb-pgen-root", type = Path, default = os.environ.get("UKB_PGEN_ROOT") or None)
	ap.add_argument("--ukb-work", type = Path, default = os.environ.get("UKB_WORK") or None)
	ap.add_argument("--ukb-results-root", type = Path, default = os.environ.get("UKB_RESULTS_ROOT") or os.environ.get("GU_UKB_RESULTS_ROOT") or None)
	ap.add_argument("--ukb-run-dir", type = Path, help = "internal stage-qc native method-run destination")
	ap.add_argument("--ukb-archaic-root", type = Path, help = "internal prepare-archaic original reference VCF directory")
	ap.add_argument("--ukb-archaic-refs", default = "Altai Vindija", help = "internal prepare-archaic reference names, separated by commas or spaces")
	ap.add_argument("--ukb-mask-root", type = Path, help = "original archaic quality-mask directory; filtering VCFs does not redefine quality masks")
	ap.add_argument("--ukb-vcf-out", type = Path, default = os.environ.get("UKB_VCF_OUT") or None)
	ap.add_argument("--ukb-ref-fasta", type = Path, default = os.environ.get("UKB_REF_FASTA") or None)
	ap.add_argument("--ukb-keep", "--keep", type = Path, default = os.environ.get("UKB_KEEP") or None)
	ap.add_argument("--ukb-info-file", default = os.environ.get("UKB_INFO_FILE") or None, help = "UKB .mfi or CHROM POS REF ALT INFO TSV; {chr} expands per chromosome")
	ap.add_argument("--ukb-info-field", default = os.environ.get("UKB_INFO_FIELD", "INFO"))
	ap.add_argument("--ukb-qc-pass-variants", default = os.environ.get("UKB_QC_PASS_VARIANTS") or None, help = "explicit already-QC-passed variant IDs, as an alternative to missing INFO; {chr} supported")
	ap.add_argument("--ukb-info-min", type = float, default = float(os.environ.get("UKB_INFO_MIN", "0.8")))
	ap.add_argument("--ukb-hardcall-threshold", type = float, default = float(os.environ.get("UKB_HARDCALL_THRESHOLD", "0.1")))
	ap.add_argument("--ukb-max-missing", type = float, default = float(os.environ.get("UKB_MAX_MISSING", "0.02")))
	ap.add_argument("--ukb-max-sample-missing", type = float, default = float(os.environ.get("UKB_MAX_SAMPLE_MISSING", "0.1")))
	ap.add_argument("--ukb-threads", type = int, default = int(os.environ.get("UKB_THREADS", "8")))
	ap.add_argument("--ukb-memory-mb", type = int, default = int(os.environ.get("GU_PREP_MEMORY_MB", "8192")))
	ap.add_argument("--ukb-pilot-n", type = int, default = int(os.environ.get("UKB_PILOT_N", "1000")))
	ap.add_argument("--ukb-max-export-samples", type = int, default = int(os.environ.get("UKB_MAX_EXPORT_SAMPLES", "10000")))
	ap.add_argument("--ukb-seed", type = int, default = int(os.environ.get("UKB_SEED", "20261006")))
	ap.add_argument("--ukb-group", default = os.environ.get("UKB_GROUP") or None)
	ap.add_argument("--ukb-pop", default = os.environ.get("UKB_POP") or None, help = "select one population from a supplied sample panel when defining the fixed pilot")
	for option in ["require-phase", "trust-hardcalls", "sparse-control"]:
		ap.add_argument("--ukb-" + option, type = boolean, default = boolean(os.environ.get("UKB_" + option.replace("-", "_").upper(), "false")))
	ap.add_argument("--keep-males", dest="male_only", action="store_true", help="PLINK --keep-males; sex comes from source PSAM")
	ap.add_argument("--existing-pgen", action="store_true", help="Use existing PGEN hardcalls; report INFO and original imputation QC as unverified")
	ap.add_argument("--regions", type=Path, help="Limit direct PGEN export to BED0 windows")
	ap.add_argument("--dataset", help="Result dataset identity")
	ap.add_argument("--ukb-no-probe", action = "store_true", help = "inspect text indexes without calling PLINK; does not assess phase/dosage")
	return ap


def main(argv = None):
	args = parser().parse_args(argv)
	args.chroms = chromosomes(args.chr)
	args.root = (args.ukb_pgen_root or Path(f"/mnt/f/gen/ukb/{args.grch}/{args.ukb_source}")).expanduser().resolve()
	args.work = temporary_path(args.ukb_work or Path("/tmp/gu-ukb") / args.ukb_source)
	for name in ["ukb_threads", "ukb_memory_mb", "ukb_pilot_n", "ukb_max_export_samples"]:
		if getattr(args, name) <= 0:
			fail(f"{name} must be positive")
	if not 0 <= args.ukb_info_min <= 1:
		fail("ukb_info_min must be in [0, 1]")
	for name in ["ukb_max_missing", "ukb_max_sample_missing"]:
		if not 0 <= getattr(args, name) < 1:
			fail(f"{name} must be in [0, 1)")
	if not 0 <= args.ukb_hardcall_threshold < 0.5:
		fail("hardcall threshold must be in [0, 0.5)")
	with (args.work / ".prepare.lock").open("w") as lock:
		try:
			fcntl.flock(lock, (fcntl.LOCK_SH if args.action in {"stage-qc", "prepare-archaic"} else fcntl.LOCK_EX) | fcntl.LOCK_NB)
		except BlockingIOError:
			fail(f"another UKB preparation is using {args.work}; use a separate work directory")
		if args.action == "inspect-pgen":
			inspect(args)
		elif args.action in {"make-panel-pgen", "pilot"}:
			panel_or_pilot(args)
		elif args.action == "stage-qc":
			stage_qc(args)
		elif args.action == "prepare-archaic":
			prepare_archaic(args)
		else:
			export_vcf(args)


if __name__ == "__main__":
	try:
		main()
	except (ValueError, OSError, subprocess.SubprocessError) as error:
		print(f"ERROR: {error}", file = sys.stderr)
		raise SystemExit(2)
