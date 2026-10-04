#!/usr/bin/env python3
"""GU phyml utilities. See --help for commands."""

from __future__ import annotations
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

if "gu_0_common" not in sys.modules:
	_spec = spec_from_file_location("gu_0_common", Path(__file__).with_name("0.common.py"))
	_common = module_from_spec(_spec)
	sys.modules[_spec.name] = _common
	try:
		_spec.loader.exec_module(_common)
	except BaseException:
		sys.modules.pop(_spec.name, None)
		raise
from gu_0_common import load_module


# 🚩 phyml_contract
"""Scientific workflow identity shared by generation and resume checks."""
WORKFLOW = "gwas_lead_ld_core_archaic5_v2"
LINEAGE_REFS = {"Neanderthal": ("Altai", "Chagyr", "Vindija"), "Denisovan": ("Denisova", "Denisova25")}
REFS = tuple(ref for refs in LINEAGE_REFS.values() for ref in refs)


# 🚩 phyml_core
"""Shared VCF and haplotype readers for the GU GWAS risk-core workflow."""


import csv
import re
import subprocess
import tempfile
from array import array
from collections import Counter
from typing import Iterator
from functools import lru_cache
from pathlib import Path

load_module("0.common.py")
from gu_0_common import enable_wide_csv_fields


enable_wide_csv_fields()

BASES = {"A", "C", "G", "T"}


# 🚩 Locus errors and command execution
class SkipLocus(RuntimeError):
	"""A chromosome-window has no usable data and may be skipped."""

	def __init__(self, message: str, code: str = "low_sequence_information", details: dict | None = None):
		super().__init__(message)
		self.code = code
		self.details = details or {}


def fail(message: str) -> None:
	raise SystemExit(f"ERROR: {message}")


def run(args: list[str]) -> str:
	p = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if p.returncode:
		fail(f"command failed ({p.returncode}): {' '.join(args)}\n{p.stderr.strip()}")
	return p.stdout


# 🚩 Sample metadata and VCF lookup
def read_sexes(path: Path) -> dict[str, str]:
	with path.open() as handle:
		reader = csv.DictReader(handle, delimiter="\t")
		if not reader.fieldnames or "sample" not in reader.fieldnames or "sex" not in reader.fieldnames:
			fail(f"{path} must contain tab-separated sample and sex columns")
		out = {}
		for row in reader:
			sex = str(row["sex"]).strip().lower()
			if sex in {"1", "m", "male"}:
				sex = "male"
			elif sex in {"2", "f", "female"}:
				sex = "female"
			elif sex in {"", "0", "na", "n/a", ".", "unknown", "u"}:
				sex = "unknown"
			else:
				fail(f"invalid sex for sample {row['sample']}: {row['sex']}")
			out[row["sample"]] = sex
	return out


@lru_cache(maxsize=None)
def vcf_path(vcf_dir: Path, chrom: str) -> Path:
	candidates = [
		vcf_dir / f"chr{chrom}.vcf.gz",
		vcf_dir / f"{chrom}.vcf.gz",
		vcf_dir / f"chr{chrom}.bcf",
		vcf_dir / f"{chrom}.bcf",
	]
	for path in candidates:
		if path.is_file() and path.stat().st_size:
			return path
	fail(f"modern VCF for chr{chrom} not found under {vcf_dir}")


@lru_cache(maxsize=None)
def contigs(vcf: Path) -> dict[str, int | None]:
	header = run(["bcftools", "view", "-h", str(vcf)])
	out: dict[str, int | None] = {}
	for line in header.splitlines():
		m = re.match(r"##contig=<ID=([^,>]+)(?:,length=([0-9]+))?", line, re.I)
		if m:
			out[m.group(1)] = int(m.group(2)) if m.group(2) else None
	return out


@lru_cache(maxsize=None)
def vcf_contig(vcf: Path, chrom: str) -> str:
	names = contigs(vcf)
	for candidate in (chrom, f"chr{chrom}", "23" if chrom == "X" else chrom):
		if candidate in names:
			return candidate
	fail(f"VCF {vcf} has no contig for chr{chrom}")


@lru_cache(maxsize=None)
def archaic_vcf(root: Path, ref: str, chrom: str) -> Path:
	aliases = [ref]
	if ref.lower() == "denisova":
		aliases += ["Denisovan"]
	dirs = []
	for alias in aliases:
		dirs += [root / alias, root / "avcf" / alias]
	patterns = [f"*chr{chrom}_*.vcf.gz", f"*chr{chrom}.*.vcf.gz", f"*chr{chrom}.vcf.gz", f"*{chrom}.vcf.gz"]
	for directory in dirs:
		if not directory.is_dir():
			continue
		for pattern in patterns:
			found = sorted(x for x in directory.glob(pattern) if x.is_file() and x.stat().st_size)
			if found:
				return found[0]
	fail(f"archaic VCF ref={ref} chr={chrom} not found under {root}")


@lru_cache(maxsize=None)
def vcf_samples(vcf: Path) -> tuple[str, ...]:
	return tuple(x for x in run(["bcftools", "query", "-l", str(vcf)]).splitlines() if x)


# 🚩 Streaming VCF queries
def query_rows(vcf: Path, region: str, include_aa: bool = False) -> tuple[list[str], Iterator[list[str]]]:
	samples = list(vcf_samples(vcf))
	fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT" + ("\t%INFO/AA" if include_aa else "") + "[\t%GT]\n"
	command = ["bcftools", "query", "-r", region, "-f", fmt, str(vcf)]

	def stream():
		# Do not retain the full sample-by-variant text or millions of GT strings.
		# A file for stderr avoids pipe deadlocks when bcftools emits warnings.
		with tempfile.TemporaryFile(mode="w+t") as errors:
			proc = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=errors)
			try:
				for line in proc.stdout:
					if line.strip():
						yield line.rstrip("\r\n").split("\t")
				if proc.wait():
					errors.seek(0)
					fail(f"command failed ({proc.returncode}): {' '.join(command)}\n{errors.read().strip()}")
			finally:
				proc.stdout.close()
				if proc.poll() is None:
					proc.terminate()
				proc.wait()

	return samples, stream()


# 🚩 Genotype and ancestral allele decoding
def called_base(ref: str, alt: str, gt: str, phased_required: bool) -> tuple[str, str]:
	alleles = [ref] + alt.split(",")
	if gt in {".", "./.", ".|."}:
		return "N", "N"
	if phased_required and "/" in gt and gt.split("/")[0] != gt.split("/")[-1]:
		return "N", "N"
	parts = re.split(r"[/|]", gt)
	if len(parts) == 1:
		parts *= 2
	if len(parts) != 2 or any(not x.isdigit() or int(x) >= len(alleles) for x in parts):
		return "N", "N"
	return alleles[int(parts[0])].upper(), alleles[int(parts[1])].upper()


def called_haploid_base(ref: str, alt: str, gt: str) -> str:
	"""Return one male-X allele and reject every diploid representation."""
	alleles = [ref] + alt.split(",")
	if gt in {".", "", "./.", ".|."}:
		return "N"
	parts = re.split(r"[/|]", gt)
	if len(parts) != 1 or not parts[0].isdigit() or int(parts[0]) >= len(alleles):
		fail(f"chrX.male contains a non-haploid/heterozygous GT: {gt}")
	return alleles[int(parts[0])].upper()


def ancestral_allele(value: str, ref: str, alt: str, allow_third_allele: bool = False) -> str:
	"""Normalize 1KG INFO/AA values such as A, A|, or A|||."""
	token = re.split(r"[|,]", str(value).strip().upper(), maxsplit=1)[0]
	return token if len(token) == 1 and token in BASES and (allow_third_allele or token in {ref, alt}) else "N"


def alt_dosage(gt: str, ploidy: int) -> int | None:
	"""Return ALT-allele dosage without requiring phased genotypes."""
	if gt in {".", "", "./.", ".|."}:
		return None
	parts = re.split(r"[/|]", gt)
	if len(parts) != ploidy or any(not x.isdigit() or int(x) > 1 for x in parts):
		return None
	return sum(int(x) for x in parts)


# 🚩 Modern haplotypes
def modern_data(
	vcf: Path,
	contig: str,
	locus: dict,
	min_mac: int,
	sexes: dict[str, str],
	x_male_only: bool,
	x_par_diploid: bool,
	ancestral_any_base: bool = False,
):
	samples, rows = query_rows(vcf, f"{contig}:{locus['start'] + 1}-{locus['end']}", include_aa=True)
	if not samples:
		fail(f"modern VCF has no samples: {vcf}")
	missing = [sample for sample in samples if sample not in sexes] if sexes else []
	if sexes and missing:
		fail(f"samples absent from samples.txt (first examples): {','.join(missing[:10])}")
	sites = []
	anchors = []
	coordinate_anchor = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)", locus["name"], re.I)
	if locus["chrom"] == "X" and x_male_only and sexes:
		nonmale = [sample for sample in samples if sexes.get(sample) != "male"]
		if nonmale:
			fail(f"chrX.male VCF contains samples not marked male (first examples): {','.join(nonmale[:10])}")
	for fields in rows:
		if len(fields) < 6 + len(samples):
			continue
		_, pos, vid, ref, alt, aa, *gts = fields
		pos_i = int(pos)
		coordinate_match = coordinate_anchor and (
			coordinate_anchor[1],
			int(coordinate_anchor[2]),
			coordinate_anchor[3].upper(),
			coordinate_anchor[4].upper(),
		) == (locus["chrom"], pos_i, ref.upper(), alt.upper())
		if locus["name"] in vid.split(";") or coordinate_match:
			anchors.append((pos_i, ref, alt))
		if len(ref) != 1 or len(alt) != 1 or ref.upper() not in BASES or alt.upper() not in BASES:
			continue
		haps = []
		ploidy = 1 if locus["chrom"] == "X" and x_male_only else 2
		# GT vocabulary is tiny; decode once per site, not once per sample.
		decoded = {}
		dosage_values = array("b")
		for sample, gt in zip(samples, gts):
			if gt not in decoded:
				if locus["chrom"] == "X" and x_male_only:
					pair = (called_haploid_base(ref.upper(), alt.upper(), gt), "N")
				elif locus["chrom"] != "X" or x_par_diploid:
					pair = called_base(ref.upper(), alt.upper(), gt, phased_required=True)
				else:
					fail("chrX input requires an explicit male non-PAR or diploid PAR mode")
				dose = alt_dosage(gt, ploidy)
				decoded[gt] = (pair, -1 if dose is None else dose)
			pair, dose = decoded[gt]
			haps.extend(pair)
			dosage_values.append(dose)
		haps = "".join(haps)
		counts = Counter(x for x in haps if x in BASES)
		if len(counts) < 2 or min(counts.values()) < min_mac:
			continue
		sites.append(
			dict(
				pos=pos_i,
				vid=vid,
				ref=ref.upper(),
				alt=alt.upper(),
				ancestral=ancestral_allele(aa, ref.upper(), alt.upper(), ancestral_any_base),
				haps=haps,
				dosages=dosage_values,
			)
		)
	if not sites:
		raise SkipLocus(f"no polymorphic biallelic SNPs with MAC >= {min_mac}")
	target_pos = anchors[0][0] if anchors else (locus["core_start"] + locus["core_end"]) // 2
	# LD selection needs an anchor that survived the biallelic/MAC filters.  An
	# exact named anchor is preferred; otherwise use the closest eligible SNP
	# and report its actual position downstream.
	anchor_pos = min(sites, key=lambda x: (abs(x["pos"] - target_pos), x["pos"]))["pos"]
	return samples, sites, anchor_pos


# 🚩 Archaic haplotypes
def archaic_calls(
	vcf: Path, contig: str, locus: dict, modern_sites: list[dict], allow_third_allele: bool = False
) -> dict[int, str]:
	_, rows = query_rows(vcf, f"{contig}:{locus['start'] + 1}-{locus['end']}")
	wanted = {x["pos"]: x for x in modern_sites}
	calls: dict[int, str] = {}
	for fields in rows:
		if len(fields) < 6:
			continue
		_, pos, _, ref, alt, *gts = fields
		pos_i = int(pos)
		if pos_i not in wanted:
			continue
		observed = []
		for gt in gts:
			observed.extend(called_base(ref.upper(), alt.upper(), gt, phased_required=False))
		# A heterozygote containing one allele outside the modern pair is
		# still a heterozygote; filtering that allele first creates a false
		# homozygous match. Missing/heterozygous archaic calls remain unknown.
		calls[pos_i] = (
			observed[0]
			if (
				observed
				and len(set(observed)) == 1
				and observed[0] in BASES
				and (allow_third_allele or observed[0] in {wanted[pos_i]["ref"], wanted[pos_i]["alt"]})
			)
			else "N"
		)
	return calls


# 🚩 phyml_thresholds
"""ILS probability used by the GWAS risk-core workflow.
"""

import math


# 🚩 Incomplete lineage sorting probability
def ils_probability(length_bp, recomb_cm_mb=0.53, split_years=550000, archaic_age_years=50000, generation_years=29):
	"""Gamma(shape=2) survival; same model as old ils_p, without 1-CDF cancellation.

	This is a model-based screen, not a calibrated introgression probability.
	The caller supplies the observed diagnostic span, never the flanked region.
	"""
	x = max(0, length_bp) * recomb_cm_mb * 1e-8 * ((2 * split_years - archaic_age_years) / generation_years)
	return math.exp(-x) * (1 + x)


# 🚩 phyml_tree_summary
"""Shared Newick parser and bootstrap support for GWAS risk trees."""


import re
from dataclasses import dataclass, field

load_module("0.common.py")
from gu_0_common import enable_wide_csv_fields


enable_wide_csv_fields()


# 🚩 Bootstrap support
def bootstrap_values(newick: str) -> list[float]:
	return [float(x) for x in re.findall(r"\)([0-9]+(?:\.[0-9]+)?)(?=[:),;])", newick)]


# 🚩 Newick tree structure and parsing
@dataclass
class NewickNode:
	label: str = ""
	children: list["NewickNode"] = field(default_factory=list)
	support: float | None = None
	node_id: str = ""


def parse_newick(text: str) -> NewickNode:
	"""Parse the simple, unquoted Newick emitted by PhyML.

	Keeping this tiny parser local avoids adding a Biopython dependency to the
	normalization path.  Tip labels have already been sanitized by safe_label.
	"""
	source = "".join(text.split())
	index = 0
	serial = 0

	def token() -> str:
		nonlocal index
		start = index
		while index < len(source) and source[index] not in ":,();":
			index += 1
		return source[start:index]

	def branch_length() -> None:
		nonlocal index
		if index < len(source) and source[index] == ":":
			index += 1
			while index < len(source) and source[index] not in ",();":
				index += 1

	def subtree() -> NewickNode:
		nonlocal index, serial
		if index >= len(source):
			raise ValueError("unexpected end of Newick")
		if source[index] != "(":
			label = token()
			if not label:
				raise ValueError("empty Newick tip")
			node = NewickNode(label=label)
			branch_length()
			return node
		index += 1
		children = [subtree()]
		while index < len(source) and source[index] == ",":
			index += 1
			children.append(subtree())
		if index >= len(source) or source[index] != ")":
			raise ValueError("unterminated Newick clade")
		index += 1
		label = token()
		support = None
		try:
			support = float(label) if label else None
		except ValueError:
			pass
		serial += 1
		node = NewickNode(label=label, children=children, support=support, node_id=f"N{serial}")
		branch_length()
		return node

	root = subtree()
	if index < len(source) and source[index] == ";":
		index += 1
	if index != len(source):
		raise ValueError(f"unexpected Newick content at offset {index}")
	return root


# 🚩 phyml_gwas_input
"""Prepare original COJO leads, never distance-clump independent GWAS signals."""
import argparse, csv, hashlib, json, math, os, re, fcntl
from pathlib import Path

load_module("0.common.py")
from gu_0_common import CHROM_LENGTHS

load_module("cojo.py")
from gu_cojo import lift


def write(path, data, fields=None):
	path.parent.mkdir(parents=True, exist_ok=True)
	fields = fields or list(dict.fromkeys(k for r in data for k in r)) or ["status"]
	tmp = path.with_name(path.name + ".tmp")
	with tmp.open("w") as f:
		w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
		w.writeheader()
		w.writerows(data)
	tmp.replace(path)


def read(path):
	with path.open() as f:
		return list(csv.DictReader(f, delimiter="\t"))


def complement(a):
	return a.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def parse_cojo(path, build):
	lines = [x.split() for x in path.read_text().splitlines() if x.strip()]
	required = {"SNP", "Chr", "bp", "refA", "bJ", "pJ"}
	if not lines or not required <= set(lines[0]):
		raise ValueError("COJO requires columns: " + ", ".join(sorted(required)))
	accepted, audit, seen = [], [], set()
	for n, vals in enumerate(lines[1:], 1):
		r = dict(zip(lines[0], vals))
		name = r.get("SNP", "")
		try:
			if len(vals) != len(lines[0]):
				raise ValueError("wrong_column_count")
			ch = re.sub("^chr", "", r["Chr"], flags=re.I).upper()
			ch = "X" if ch == "23" else ch
			pos = int(r["bp"])
			beta = float(r["bJ"])
			p = float(r["pJ"])
			effect = r["refA"].upper()
			if ch not in CHROM_LENGTHS[build] or not 1 <= pos <= CHROM_LENGTHS[build][ch]:
				raise ValueError("invalid_position")
			if not math.isfinite(beta) or beta == 0 or not math.isfinite(p) or not 0 <= p <= 1:
				raise ValueError("invalid_effect_or_P")
			if not re.fullmatch("[ACGT]+", effect):
				raise ValueError("invalid_effect_allele")
			m = re.fullmatch(r"(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)", name, re.I)
			ref, alt = ("", "")
			if m:
				mc = "X" if m[1] == "23" else m[1].upper()
				if mc != ch or int(m[2]) != pos:
					raise ValueError("SNP_ID_disagrees_with_CHR_POS")
				ref, alt = m[3].upper(), m[4].upper()
				if ref == alt or effect not in (ref, alt):
					raise ValueError("effect_allele_disagrees_with_SNP_ID")
			elif not re.fullmatch("rs[0-9]+", name):
				raise ValueError("unsupported_SNP_ID")
			if name in seen:
				raise ValueError("duplicate_SNP_ID")
			seen.add(name)
			accepted.append(
				dict(
					input_row=n,
					locus_id=name,
					index_snp=name,
					source_build="GRCh" + build,
					source_chr=ch,
					source_pos=pos,
					chr=ch,
					lead_pos=pos,
					ref=ref,
					alt=alt,
					effect_allele=effect,
					beta_j=beta,
					p_j=p,
					strand="+",
				)
			)
		except (ValueError, KeyError) as e:
			audit.append(dict(input_row=n, index_snp=name, status="skipped", reason=str(e)))
	return accepted, audit


def prepare(a):
	out = a.output
	out.mkdir(parents=True, exist_ok=True)
	leads, audit = parse_cojo(a.input, a.build)
	provenance = dict(
		input=str(a.input.resolve()),
		input_sha256=hashlib.sha256(a.input.read_bytes()).hexdigest(),
		source_build=a.build,
		analysis_build="37",
		method="original_COJO_leads_no_distance_collapse",
		search_flank_bp=a.window,
	)
	if a.build == "38" and leads:
		binary = Path(os.environ.get("GU_LIFTOVER", "/mnt/d/software/bin/liftOver"))
		print(
			f"[GU GWAS] liftOver={binary.resolve()} sha256={hashlib.sha256(binary.read_bytes()).hexdigest()}",
			flush=True,
		)
		source = out / "leads.GRCh38.bed"
		source.write_text(
			"".join(
				f"chr{r['chr']}\t{r['lead_pos'] - 1}\t{r['lead_pos'] - 1 + max(1, len(r['ref']))}\t{r['input_row']}\t0\t+\n"
				for r in leads
			)
		)
		hits = lift(binary, a.chain, source, out / "leads.lifted.GRCh37.bed", out / "leads.unmapped.bed")
		provenance.update(
			liftover_binary=str(binary.resolve()),
			liftover_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
			chain=str(a.chain),
			chain_sha256=hashlib.sha256(a.chain.read_bytes()).hexdigest(),
		)
		mapped = []
		for r in leads:
			h = hits[str(r["input_row"])]
			reason = ""
			if len(h) != 1:
				reason = "lead_unmapped_or_multiple_mapping"
			else:
				x = h[0]
				ch = x[0].removeprefix("chr")
				pos = int(x[1]) + 1
				if ch != r["chr"] or ch not in CHROM_LENGTHS["37"] or int(x[2]) - int(x[1]) != max(1, len(r["ref"])):
					reason = "lead_mapping_changes_chromosome_or_allele_span"
				else:
					r.update(chr=ch, lead_pos=pos, strand=x[5])
					if x[5] == "-":
						for key in ("ref", "alt", "effect_allele"):
							r[key] = complement(r[key])
			if reason:
				audit.append(dict(input_row=r["input_row"], index_snp=r["index_snp"], status="skipped", reason=reason))
			else:
				mapped.append(r)
		leads = mapped
	for r in leads:
		r.update(
			search_start=max(0, r["lead_pos"] - 1 - a.window),
			search_end=min(CHROM_LENGTHS["37"][r["chr"]], r["lead_pos"] + a.window),
		)
	write(out / "gwas_leads.GRCh37.tsv", leads)
	write(out / "input.audit.tsv", audit, ["input_row", "index_snp", "status", "reason"])
	for r in audit:
		print(f"[GU GWAS] SKIP row={r['input_row']} SNP={r['index_snp']} reason={r['reason']}", flush=True)
	if not leads:
		raise ValueError("No valid original COJO leads remain")
	# These are extraction windows, not inferred core haplotypes or standalone inputs.
	(out / "search_windows.GRCh37.bed").write_text(
		"".join(f"{r['chr']}\t{r['search_start']}\t{r['search_end']}\t{r['locus_id']}\n" for r in leads)
	)
	provenance.update(n_leads=len(leads), n_skipped=len(audit))
	(out / "manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")
	print(
		f"[GU GWAS] {len(leads)} original leads retained; {len(audit)} skipped; LD=1KG EUR; cores will be defined by r² > 0.98",
		flush=True,
	)
	print(
		f"[GU GWAS] lead audit: {out / 'gwas_leads.GRCh37.tsv'}; extraction BED is NOT a standalone PhyML input",
		flush=True,
	)


def phyml_gwas_input_cli():
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("--input", type=Path, required=True)
	p.add_argument("--output", type=Path, required=True)
	p.add_argument("--build", choices=["37", "38"], required=True)
	p.add_argument("--window", type=int, default=500000)
	p.add_argument("--chain", type=Path, default=Path("/mnt/d/files/liftOver/hg38ToHg19.over.chain.gz"))
	args = p.parse_args()
	args.output.mkdir(parents=True, exist_ok=True)
	with (args.output / ".prepare.lock").open("w") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		prepare(args)


# 🚩 phyml_run
"""Run one PhyML tree; a pre-bootstrap main tree is never a completion marker."""
import argparse, fcntl, hashlib, json, os, re, shutil, signal, subprocess, sys, time
from pathlib import Path

SUFFIXES = [
	"_phyml_tree.txt",
	"_phyml_stats.txt",
	"_phyml_boot_trees.txt",
	"_phyml_boot_stats.txt",
	"_phyml_tree.png",
	"_phyml_tree.pdf",
	".phyml.log",
	".phyml.complete.json",
]


def phyml_run_digest(p):
	h = hashlib.sha256()
	with Path(p).open("rb") as f:
		for block in iter(lambda: f.read(1024 * 1024), b""):
			h.update(block)
	return h.hexdigest()


def phyml_run_outputs(phy, boot):
	names = ["_phyml_tree.txt", "_phyml_stats.txt", ".phyml.log"]
	if boot > 0:
		names += ["_phyml_boot_trees.txt", "_phyml_boot_stats.txt"]
	return [Path(str(phy) + s) for s in names]


def completion_error(phy, boot):
	files = phyml_run_outputs(phy, boot)
	if any(not p.is_file() or not p.stat().st_size for p in files):
		return "missing/empty tree, stats, log or bootstrap outputs"
	log = Path(str(phy) + ".phyml.log").read_text(errors="replace")
	finish = list(re.finditer(r"\. Time used\s+\d+h\d+m\d+s", log))
	if not finish:
		return "no final runtime footer (interrupted run)"
	if not re.search(r"Printing the most likely tree", log):
		return "missing final tree report"
	if boot > 0:
		progress = list(re.finditer(r"\b" + str(boot) + r"\s*/\s*" + str(boot) + r"\b", log))
		if not progress or progress[-1].end() > finish[-1].start():
			return f"bootstrap did not finish {boot}/{boot}"
		with Path(str(phy) + "_phyml_boot_trees.txt").open() as f:
			count = 0
			for line in f:
				line = line.strip()
				if not line:
					continue
				if not line.endswith(";") or line.count("(") != line.count(")"):
					return "invalid bootstrap Newick"
				count += 1
		if count != boot:
			return f"bootstrap trees={count}, expected={boot}"
	tree = Path(str(phy) + "_phyml_tree.txt").read_text().strip()
	if not tree.endswith(";") or tree.count("(") != tree.count(")"):
		return "invalid final Newick"
	return None


def phyml_run_request(phy, boot, binary, seed=None):
	data = dict(
		schema=1,
		input_sha256=phyml_run_digest(phy),
		binary=str(binary),
		binary_sha256=phyml_run_digest(binary),
		model="HKY85",
		categories=4,
		alpha="e",
		invariant="e",
		bootstrap=boot,
	)
	if seed is not None:
		data["seed"] = seed
	return data


def seal(phy, req, execution=None):
	p = Path(str(phy) + ".phyml.complete.json")
	q = p.with_suffix(".next")
	data = dict(request=req, outputs={str(x): phyml_run_digest(x) for x in phyml_run_outputs(phy, req["bootstrap"])})
	if execution is not None:
		data["execution"] = execution
	q.write_text(json.dumps(data, indent=2) + "\n")
	q.replace(p)


def reusable(phy, req, adopt=True):
	why = completion_error(phy, req["bootstrap"])
	if why:
		return False, why
	receipt = Path(str(phy) + ".phyml.complete.json")
	if receipt.exists():
		try:
			data = json.loads(receipt.read_text())
			if data["request"] != req:
				return False, "input/software/options changed"
			if set(data["outputs"]) != {str(x) for x in phyml_run_outputs(phy, req["bootstrap"])}:
				return False, "incomplete output receipt"
			if any(phyml_run_digest(p) != d for p, d in data["outputs"].items()):
				return False, "completed output changed"
		except (OSError, ValueError, KeyError):
			return False, "invalid completion receipt"
		return True, "verified completion receipt"
	# Adopt historical results only when their own full log proves completion
	# and the input predates the main tree. Subsequent reuse uses content hashes.
	log = Path(str(phy) + ".phyml.log").read_text(errors="replace")
	import shlex

	command = re.search(r"\. Command line:\s*(.*)", log)
	try:
		args = shlex.split(command.group(1)) if command else []
	except ValueError:
		args = []
	if not args or Path(args[0]).name != Path(req["binary"]).name:
		return False, "historical executable does not match"
	options = [("-i", str(phy)), ("-m", "HKY85"), ("-c", "4"), ("-a", "e"), ("-v", "e"), ("-b", str(req["bootstrap"]))]
	if "seed" in req:
		options.append(("--r_seed", str(req["seed"])))
	for key, val in options:
		if key not in args or args.index(key) + 1 >= len(args) or args[args.index(key) + 1] != val:
			return False, "historical command does not match requested input/options"
	if phy.stat().st_mtime_ns > Path(str(phy) + "_phyml_tree.txt").stat().st_mtime_ns:
		return False, "input newer than historical tree"
	if adopt:
		seal(phy, req)
	return True, "verified complete historical bootstrap"


def protected_result(phy, boot=100):
	# A receipt is also evidence of prior completion when a product was lost.
	return Path(str(phy) + ".phyml.complete.json").is_file() or completion_error(phy, boot) is None


def require_replace(phy, reason):
	raise RuntimeError(
		f"completed PhyML result preserved: {reason}; input={phy}; "
		"use --replace-phyml TRUE only to explicitly replace it"
	)


def prepare_input(phy, text, replace=False):
	old = phy.read_text() if phy.exists() else None
	# Whitespace formatting alone is not a different alignment. Keep the exact
	# original bytes and mtime so its completion hashes remain valid.
	if old is not None and old.split() == text.split():
		return
	if not replace and protected_result(phy):
		require_replace(phy, "alignment changed or original input is missing")
	clean(phy)
	phy.write_text(text)


def clean(phy):
	for suffix in SUFFIXES:
		Path(str(phy) + suffix).unlink(missing_ok=True)
	# Only products for this exact input, never other loci/trees.
	for p in list(phy.parent.glob(phy.name + "_phyml_tree.panelB*")) + list(
		phy.parent.glob(phy.name + "_phyml_tree.*.panelB*")
	):
		if p.is_file():
			p.unlink()


def reuse_identical_sibling(phy, req):
	"""Only byte-identical alignments and identical settings may share a tree."""
	for source in sorted(phy.parent.glob("haplotypes.evidence.*.phy")):
		if source == phy or phyml_run_digest(source) != req["input_sha256"]:
			continue
		ok, _ = reusable(source, req)
		if not ok:
			continue
		clean(phy)
		for src, dst in zip(phyml_run_outputs(source, req["bootstrap"]), phyml_run_outputs(phy, req["bootstrap"])):
			shutil.copy2(src, dst)
		with Path(str(phy) + ".phyml.log").open("a") as f:
			f.write(f"\n[GU PHYML] Reused identical alignment and model from {source}\n")
		seal(phy, req)
		return source
	return None


def state(phy, scope, status, rc, elapsed=0):
	p = Path(str(phy) + ".phyml.run.status.tsv")
	q = p.with_suffix(".next")
	q.write_text(
		"scope\tstatus\trc\telapsed_seconds\tphy\ttree\tlog\n"
		+ f"{scope}\t{status}\t{rc}\t{elapsed}\t{phy}\t{phy}_phyml_tree.txt\t{phy}.phyml.log\n"
	)
	q.replace(p)


def stop_process_group(proc):
	# Distribution launchers may spawn MPI/worker children.
	try:
		os.killpg(proc.pid, signal.SIGTERM)
	except ProcessLookupError:
		return
	try:
		proc.wait(timeout=10)
	except subprocess.TimeoutExpired:
		pass
	try:
		os.killpg(proc.pid, signal.SIGKILL)
	except ProcessLookupError:
		pass
	proc.wait()


def interrupted(signum, frame):
	raise RuntimeError(f"PhyML interrupted by signal {signum}")


def discard_partial(phy):
	for suffix in SUFFIXES:
		if suffix != ".phyml.log":
			Path(str(phy) + suffix).unlink(missing_ok=True)


def run_attempt(phy, cmd, deadline):
	remaining = None if deadline is None else deadline - time.monotonic()
	proc = None
	try:
		with Path(str(phy) + ".phyml.log").open("w") as f:
			if remaining is not None and remaining <= 0:
				f.write("[GU PHYML] Timeout budget exhausted before starting this attempt\n")
				return 124
			# MPI otherwise consumes the caller's remaining task-list lines.
			proc = subprocess.Popen(
				cmd, stdin=subprocess.DEVNULL, stdout=f, stderr=subprocess.STDOUT, start_new_session=True
			)
			while True:
				remaining = None if deadline is None else deadline - time.monotonic()
				if remaining is not None and remaining <= 0:
					stop_process_group(proc)
					return 124
				try:
					return proc.wait(timeout=60 if remaining is None else min(60, remaining))
				except subprocess.TimeoutExpired:
					if deadline is not None and time.monotonic() >= deadline:
						stop_process_group(proc)
						return 124
					# Continue enforcing the deadline silently.
	finally:
		if proc is not None and proc.poll() is None:
			stop_process_group(proc)


def archive_failure(phy, req, attempt, why):
	prefix = str(phy) + f".phyml.failed.{time.time_ns()}"
	log = Path(prefix + ".log")
	shutil.copy2(str(phy) + ".phyml.log", log)
	Path(prefix + ".json").write_text(
		json.dumps(dict(request=req, attempt=attempt, completion_error=why, log=str(log)), indent=2) + "\n"
	)
	return str(log)


def phyml_run_main():
	p = argparse.ArgumentParser()
	p.add_argument("--phy", type=Path, required=True)
	p.add_argument("--scope", required=True)
	p.add_argument("--bootstrap", type=int, default=100)
	p.add_argument("--timeout", type=float, default=0)
	p.add_argument("--cpus", type=int, default=1)
	p.add_argument("--seed", type=int)
	p.add_argument("--mpi-fallback", choices=["serial", "error"], default="serial")
	p.add_argument("--verify-only", action="store_true")
	a = p.parse_args()
	phy = a.phy.resolve()
	if a.bootstrap < 0 or a.timeout < 0 or a.cpus < 1:
		p.error("bootstrap/timeout must be nonnegative and cpus positive")
	serial = shutil.which("phyml")
	mpi = shutil.which("phyml-mpi")
	mpirun = shutil.which("mpirun")
	if not serial:
		raise RuntimeError("phyml executable unavailable")
	with Path(str(phy) + ".phyml.lock").open("a") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
		return run_locked(a, phy, serial, mpi, mpirun)


def run_locked(a, phy, serial, mpi, mpirun):
	force = os.environ.get("PHYML_REPLACE") == "TRUE" and not a.verify_only
	# Changing the execution CPU count alone must not replace a complete tree.
	for binary in dict.fromkeys(x for x in (serial, mpi) if x):
		req = phyml_run_request(phy, a.bootstrap, Path(binary).resolve(), a.seed)
		ok, why = reusable(phy, req, adopt=not a.verify_only)
		if ok and not force:
			if not a.verify_only:
				state(phy, a.scope, "COMPLETE_VERIFIED", 0)
			print(f"[GU PHYML] tree={'VERIFIED' if a.verify_only else 'SKIP'} reason={why} input={phy}", flush=True)
			return 0
	if a.verify_only:
		print(f"ERROR: planned tree is not complete: {why}; input={phy}", file=sys.stderr)
		return 1
	if not force and protected_result(phy, a.bootstrap):
		require_replace(phy, why)
	binary = serial
	launcher = []
	if a.cpus > 1 and a.bootstrap > 0:
		if mpi and mpirun:
			# PhyML rounds replicates UP to a multiple of MPI ranks. Choose a
			# divisor so requesting 100 always produces exactly 100, not 112.
			a.cpus = max(n for n in range(1, min(a.cpus, a.bootstrap) + 1) if a.bootstrap % n == 0)
			binary = mpi
			launcher = [mpirun, "--bind-to", "none", "-np", str(a.cpus)]
		else:
			print("[GU PHYML] MPI unavailable; using one CPU", flush=True)
	req = phyml_run_request(phy, a.bootstrap, Path(binary).resolve(), a.seed)
	reused = None if force else reuse_identical_sibling(phy, req)
	if reused:
		state(phy, a.scope, "COMPLETE_REUSED", 0)
		print(f"[GU PHYML] tree=REUSE identical_input={reused} input={phy}", flush=True)
		return 0
	if force:
		why = "explicit --replace-phyml TRUE"
	print(f"[GU PHYML] tree=REPLACE reason={why} input={phy}", flush=True)
	clean(phy)
	state(phy, a.scope, "RUNNING", "")
	start = time.monotonic()
	deadline = start + a.timeout if a.timeout else None
	cmd = launcher + [binary, "-i", str(phy), "-m", "HKY85", "-c", "4", "-a", "e", "-v", "e", "-b", str(a.bootstrap)]
	if a.seed is not None:
		cmd += ["--r_seed", str(a.seed)]
	print(f"[GU PHYML] bootstrap_cpus={a.cpus if launcher else 1} input={phy}", flush=True)
	attempts = []
	rc = 1
	signal.signal(signal.SIGTERM, interrupted)
	signal.signal(signal.SIGINT, interrupted)
	try:
		rc = run_attempt(phy, cmd, deadline)
		why = completion_error(phy, a.bootstrap)
		log = Path(str(phy) + ".phyml.log").read_text(errors="replace")
		match = re.search(r"\. Random seed:\s*(\d+)", log)
		seed = int(match.group(1)) if match else a.seed
		attempts.append(dict(backend="mpi" if launcher else "serial", command=cmd, rc=rc, seed=seed))
		# One computational-error recovery, never a retry for weak support,
		# incomplete zero-exit output, timeout, or interruption. MPI and serial
		# have different random streams even with the same initial seed.
		mpi_error = (
			rc > 0 and rc not in (124, 130, 137, 143) and bool(re.search(r"MPI_ERR_[A-Z_]+|MPI_ERRORS_ARE_FATAL", log))
		)
		if launcher and (rc == 1 or mpi_error) and a.mpi_fallback == "serial" and seed is not None:
			failed_log = archive_failure(phy, req, attempts[-1], why)
			attempts[-1]["failed_log"] = failed_log
			if deadline is not None and time.monotonic() >= deadline:
				rc = 124
			else:
				print(f"[GU PHYML] mpi=FAILED rc={rc} fallback=serial seed={seed} failed_log={failed_log}", flush=True)
				clean(phy)
				req = phyml_run_request(phy, a.bootstrap, Path(serial).resolve(), a.seed)
				cmd = [
					serial,
					"-i",
					str(phy),
					"-m",
					"HKY85",
					"-c",
					"4",
					"-a",
					"e",
					"-v",
					"e",
					"-b",
					str(a.bootstrap),
					"--r_seed",
					str(seed),
				]
				rc = run_attempt(phy, cmd, deadline)
				why = completion_error(phy, a.bootstrap)
				attempts.append(dict(backend="serial", command=cmd, rc=rc, seed=seed))
		if rc == 0 and why is None:
			seal(phy, req, dict(attempts=attempts))
			state(phy, a.scope, "COMPLETE", 0, int(time.monotonic() - start))
			return 0
		print(f"ERROR: incomplete PhyML rc={rc}: {why}; input={phy}", file=sys.stderr)
		# Prevent this invocation's summary stages from reading a partial tree.
		discard_partial(phy)
		state(phy, a.scope, "FAILED", rc or 1, int(time.monotonic() - start))
		return rc or 1
	except Exception:
		discard_partial(phy)
		state(phy, a.scope, "FAILED", 1, int(time.monotonic() - start))
		raise


def phyml_run_cli():
	try:
		sys.exit(phyml_run_main())
	except (OSError, RuntimeError, ValueError) as e:
		print(f"ERROR: {e}", file=sys.stderr)
		sys.exit(1)


# 🚩 phyml_locus_cache
"""Validate whole-locus completion before launching a GU worker.

Large reference files use path/size/mtime fingerprints; small results use SHA256.
Historical runs require an explicit successful worker footer, never just TSVs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex


def phyml_locus_cache_digest(path):
	h = hashlib.sha256()
	with path.open("rb") as f:
		for block in iter(lambda: f.read(1024 * 1024), b""):
			h.update(block)
	return h.hexdigest()


def rows(path):
	with path.open() as f:
		return list(csv.DictReader(f, delimiter="\t"))


def command(path):
	env, argv = {}, []
	for line in path.read_text().splitlines():
		words = shlex.split(line)
		if words[:1] == ["export"]:
			for word in words[1:]:
				key, value = word.split("=", 1)
				env[key] = value
		elif words[:1] == ["exec"]:
			argv = words[1:]
	if len(argv) < 2 or argv[1] != "phyml":
		raise ValueError("not a phyml worker")
	args = argv[2:]
	if args[:1] == ["run"]:
		args = args[1:]
	if len(args) % 2 or any(not k.startswith("--") for k in args[::2]):
		raise ValueError("unsupported worker arguments")
	opts = dict(zip(args[::2], args[1::2]))
	return env, opts


def phyml_locus_cache_request(cmd, archaic_root):
	env, opts = command(cmd)
	bed = Path(opts["--loci"]).read_text().split()
	if len(bed) != 4:
		raise ValueError("expected one locus")
	lead = [r for r in rows(Path(env["GU_PHYML_LEAD_TABLE"])) if r["locus_id"] == bed[3]]
	if len(lead) != 1:
		raise ValueError("missing original lead")
	prefix = Path(opts["--target-dir"])
	ch = bed[0].removeprefix("chr")
	sources = [
		p
		for p in prefix.parent.glob(prefix.name + ch + ".*")
		if p.name.endswith(
			(".pgen", ".pvar", ".pvar.zst", ".psam", ".bed", ".bim", ".fam", ".vcf.gz", ".bcf", ".tbi", ".csi")
		)
	]
	if not sources:
		raise ValueError("missing target input")
	panel = opts.get("--sample-panel") or os.environ.get("GU_SAMPLE_PANEL")
	if not panel:
		panel = next(
			(str(p) for p in (prefix.parent / "samples.txt", prefix.parent.parent / "samples.txt") if p.is_file()), None
		)
	if not panel:
		raise ValueError("sample panel unavailable")
	sources.append(Path(panel))
	root = Path(archaic_root)
	if not root.is_dir():
		raise ValueError("archaic root unavailable")
	sources += [
		p for p in root.rglob("*") if re.search(r"(?<![0-9])(?:chr)?" + re.escape(ch) + r"[_.]", p.name) and p.is_file()
	]
	stamps = {}
	for p in sorted(set(sources)):
		if p.is_file():
			s = p.stat()
			stamps[str(p.resolve())] = [s.st_size, s.st_mtime_ns]
	ignored = {
		"--memory-cap",
		"--foreground",
		"--auto-final",
		"--replace-phyml",
		"--phyml-jobs",
		"--loci",
		"--sample-panel",
	}
	code = {name: phyml_locus_cache_digest(Path(__file__).with_name(name)) for name in ("phyml.py", "phyml.R")}
	return dict(
		schema=1,
		workflow=WORKFLOW,
		lead=lead[0],
		bed=bed,
		options={k: v for k, v in opts.items() if k not in ignored},
		sources=stamps,
		sample_panel=str(Path(panel).resolve()),
		archaic_root=str(root.resolve()),
		code=code,
		environment={k: os.environ.get(k, "") for k in ("GU_CHRX_MALE_ONLY", "GU_CHRX_PAR_DIPLOID")},
	)


def phyml_locus_cache_outputs(out):
	required = (
		"gwas_loci.tsv",
		"gwas_haplotypes.tsv",
		"gwas_copies.tsv",
		"gwas_lead.tsv",
		"gwas_parameters.json",
		"loci.tsv",
		"trees.tsv",
		"evidence_trees.tsv",
		"haplotypes.tsv",
		"haplotype_samples.tsv",
		"skipped_loci.tsv",
	)
	if any(not (out / "final" / name).is_file() or not (out / "final" / name).stat().st_size for name in required):
		raise ValueError("missing final outputs")
	return {
		str(p.relative_to(out)): phyml_locus_cache_digest(p)
		for folder in ("final", "loci")
		for p in sorted((out / folder).glob("*"))
		if p.is_file() and not p.name.endswith(".lock") and ".failed." not in p.name
	}


def successful(cmd, req):
	out = cmd.parent
	log = cmd.with_suffix(".log").read_text()
	if cmd.with_suffix(".err").exists() or f"analysis_unit={cmd.stem} status=complete" not in log:
		raise ValueError("no successful worker completion")
	if "status=failed" in log or "plot export failed" in log:
		raise ValueError("failed worker")
	if f"archaic reference root={req['archaic_root']}\n" not in log:
		raise ValueError("reference root changed")
	if f"sample metadata={req['sample_panel']}\n" not in log:
		raise ValueError("sample panel changed")
	if rows(out / "final/gwas_lead.tsv") != [req["lead"]]:
		raise ValueError("lead changed")
	params = json.loads((out / "final/gwas_parameters.json").read_text())
	if params["workflow"] != WORKFLOW:
		raise ValueError("historical workflow changed")
	summaries = rows(out / "final/gwas_loci.tsv")
	if len(summaries) != 2 or any(r["status"] in ("tree_failed", "not_evaluable") for r in summaries):
		raise ValueError("incomplete locus")
	trees = rows(out / "final/trees.tsv")
	if len(trees) != 1 or trees[0]["tree_status"] not in ("complete", "not_run"):
		raise ValueError("incomplete tree summary")
	if (
		any(r["status"] in ("tree_supported", "tree_not_supported") for r in summaries)
		and trees[0]["tree_status"] != "complete"
	):
		raise ValueError("tree result without completed tree")
	for tree in trees:
		if tree["tree_status"] == "complete":
			if completion_error(out / "loci/haplotypes.phy", 100):
				raise ValueError("incomplete tree")
			for lineage in ("Neanderthal", "Denisovan"):
				for suffix in (".png", ".pdf", ".full.png", ".full.pdf"):
					plot = out / "loci" / f"haplotypes.phy_phyml_tree.{lineage}.panelB{suffix}"
					if not plot.is_file() or not plot.stat().st_size:
						raise ValueError("missing plot")


def process(mode, cmd, archaic_root):
	receipt = cmd.parent / ".phyml.locus.complete.json"
	if mode == "adopt" and receipt.exists():
		return False
	if mode == "adopt":
		# Reject interrupted historical workers before scanning references.
		log = cmd.with_suffix(".log").read_text()
		if cmd.with_suffix(".err").exists() or f"analysis_unit={cmd.stem} status=complete" not in log:
			raise ValueError("no successful worker completion")
	req = phyml_locus_cache_request(cmd, archaic_root)
	if mode == "check":
		if command(cmd)[1].get("--replace-phyml", "FALSE") == "TRUE":
			return False
		if not receipt.is_file():
			return process("adopt", cmd, archaic_root)
		data = json.loads(receipt.read_text())
		# Code hashes record provenance, not permission to replace completed
		# analyses. Actual inputs/options and output integrity still must match.
		saved = {k: v for k, v in data["request"].items() if k != "code"}
		current = {k: v for k, v in req.items() if k != "code"}
		return saved == current and data["outputs"] == phyml_locus_cache_outputs(cmd.parent)
	successful(cmd, req)
	if mode == "adopt":
		# Old runs have no source fingerprints. Only adopt sources older than
		# their successful log, with the original command and exact saved lead.
		finished = cmd.with_suffix(".log").stat().st_mtime_ns
		if any(stamp[1] > finished for stamp in req["sources"].values()):
			return False
		if cmd.stat().st_mtime_ns > finished:
			return False
	data = dict(request=req, outputs=phyml_locus_cache_outputs(cmd.parent))
	tmp = receipt.with_name(receipt.name + f".{os.getpid()}.tmp")
	tmp.write_text(json.dumps(data, indent=2) + "\n")
	tmp.replace(receipt)
	return True


def phyml_locus_cache_main():
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("mode", choices=("check", "seal", "adopt", "partition"))
	p.add_argument("cmd", type=Path)
	p.add_argument("--pending", type=Path)
	p.add_argument("--archaic-root", required=True)
	a = p.parse_args()
	if a.mode == "partition":
		if a.pending is None:
			p.error("partition requires --pending")
		commands = [Path(line) for line in a.cmd.read_text().splitlines() if line]
		pending = []
		skipped = 0
		skipped_results = 0
		for index, cmd in enumerate(commands, 1):
			try:
				ok = process("check", cmd, a.archaic_root)
			except (OSError, ValueError, KeyError, TypeError):
				ok = False
			if ok:
				skipped += 1
				skipped_results += bool(rows(cmd.parent / "final/skipped_loci.tsv"))
			else:
				pending.append(str(cmd))
		a.pending.write_text("".join(cmd + "\n" for cmd in pending))
		print(
			f"[GU CMD] RESUME total={len(commands)} reused={skipped} pending={len(pending)} skipped={skipped_results}",
			flush=True,
		)
		return
	try:
		ok = process(a.mode, a.cmd, a.archaic_root)
	except (OSError, ValueError, KeyError, TypeError):
		ok = False
	raise SystemExit(0 if ok else 1)


def phyml_locus_cache_cli():
	phyml_locus_cache_main()


# 🚩 phyml_gwas
"""GWAS-anchored LD-core phylogeny following Zeberg & Pääbo (2020).

EUR defines LD. All sampled populations supply recurrent modern haplotypes.
No archaic-similarity search, replacement index SNP, or candidate tip sampling.
"""
import argparse, hashlib, json, math, os, subprocess, sys, fcntl
from collections import Counter, defaultdict
from pathlib import Path


def reference_lineage(ref):
	return next(lineage for lineage, refs in LINEAGE_REFS.items() if ref in refs)


def risk_from_effect(ref, alt, effect, beta):
	if effect not in (ref, alt) or "," in alt or ref == alt:
		raise SkipLocus("effect allele does not match exact target variant", "lead_allele_mismatch")
	return effect if beta > 0 else (alt if effect == ref else ref)


def exact_lead(row, records, n_samples, haploid=False):
	found = []
	for f in records:
		if int(f[1]) != int(row["lead_pos"]):
			continue
		ref, alt = f[3].upper(), f[4].upper()
		if row["ref"]:
			if (ref, alt) != (row["ref"], row["alt"]):
				continue
		elif row["index_snp"] not in f[2].split(";"):
			continue
		if "," in alt or len(f) != 5 + n_samples:
			continue
		found.append(f)
	if len(found) != 1:
		raise SkipLocus("exact original lead absent or ambiguous; no replacement SNP", "lead_absent_or_ambiguous")
	f = found[0]
	risk = risk_from_effect(f[3], f[4], row["effect_allele"], float(row["beta_j"]))
	copies = []
	for gt in f[5:]:
		alleles = (called_haploid_base(f[3], f[4], gt), "N") if haploid else called_base(f[3], f[4], gt, True)
		copies.extend(1 if x == risk else 0 if x in (f[3], f[4]) else -1 for x in alleles)
	return dict(pos=int(f[1]), ref=f[3], alt=f[4], vcf_id=f[2], risk_allele=risk, copies=copies)


def phased_ld(risk, site, indexes):
	n = sx = sy = sxy = 0
	for i in indexes:
		a, b = risk[i], site["haps"][i]
		if a < 0 or b not in (site["ref"], site["alt"]):
			continue
		y = int(b == site["alt"])
		n += 1
		sx += a
		sy += y
		sxy += a * y
	denom = sx * (n - sx) * sy * (n - sy)
	if n < 4 or denom == 0:
		return None, n, ""
	cov = n * sxy - sx * sy
	return min(1.0, cov * cov / denom), n, site["alt"] if cov > 0 else site["ref"]


def define_core(sites, lead, indexes):
	ld = []
	high = []
	for s in sites:
		r2, n, allele = phased_ld(lead["copies"], s, indexes)
		passed = r2 is not None and r2 > 0.98
		ld.append(
			dict(
				pos=s["pos"],
				id=s["vid"],
				ref=s["ref"],
				alt=s["alt"],
				ld_r2=r2,
				n_EUR_copies=n,
				risk_linked_allele=allele,
				core_marker=int(passed),
			)
		)
		if passed:
			high.append(s)
	if len(high) < 2:
		raise SkipLocus("fewer than two SNPs with EUR r² > 0.98", "insufficient_high_LD_markers", dict(ld=ld))
	start = min([lead["pos"]] + [s["pos"] for s in high]) - 1
	end = max([lead["pos"] + len(lead["ref"]) - 1] + [s["pos"] for s in high])
	return [s for s in sites if start < s["pos"] <= end], ld, start, end


def risk_clade(newick, risk_tips, modern_tips, refs=LINEAGE_REFS["Neanderthal"]):
	"""Test all risk tips + the specified lineage, excluding other archaics."""
	root = parse_newick(newick)
	edges = []

	def visit(node):
		tips = {node.label} if not node.children else set().union(*(visit(c) for c in node.children))
		edges.append((node, tips))
		return tips

	all_tips = visit(root)
	required = set(risk_tips) | set(refs)
	if (
		not risk_tips
		or not required <= all_tips
		or "Ancestral" not in all_tips
		or not (set(modern_tips) - set(risk_tips))
	):
		return None
	matches = []
	for node, tips in edges:
		if node is root:
			continue
		for side in (tips, all_tips - tips):
			if side == required:
				matches.append(dict(bootstrap=node.support, tips=",".join(sorted(side)), node=node.node_id))
	return max(matches, key=lambda x: -1 if x["bootstrap"] is None else x["bootstrap"], default=None)


def prepare_sequences(sites, calls, samples, lead, haploid):
	# Use the same sites and modern tips for both lineage tests. All five
	# archaic references must be callable. Missing ancestry
	# remains N (never substitute REF); constant columns are retained after
	# singleton removal, as required for HKY+I parameter estimation.
	sites = [s for s in sites if all(calls[r].get(s["pos"], "N") in BASES for r in REFS)]
	indexes = [2 * i + h for i in range(len(samples)) for h in range(1 if haploid else 2)]
	grouped = defaultdict(list)
	for i in indexes:
		if lead["copies"][i] < 0:
			continue
		seq = "".join(s["haps"][i] for s in sites)
		grouped[seq].append(i)
	recurrent = []
	for seq, ix in sorted(grouped.items(), key=lambda x: (-len(x[1]), x[0])):
		if len(ix) < 2:
			continue
		flags = {lead["copies"][i] for i in ix}
		role = "risk" if flags == {1} else "nonrisk" if flags == {0} else "mixed"
		recurrent.append(dict(hap_id=f"H{len(recurrent) + 1:05d}", seq=seq, indices=ix, role=role, n=len(ix)))
	return sites, recurrent, grouped


def run_verified_tree(command):
	result = subprocess.run(command)
	if result.returncode == 0:
		result = subprocess.run(command + ["--verify-only"])
	return result


def run_locus(a, row):
	out = a.out
	final = out / "final"
	loc = out / "loci"
	final.mkdir(parents=True, exist_ok=True)
	loc.mkdir(parents=True, exist_ok=True)
	row = dict(row)
	lid = row["locus_id"]
	ch = row["chr"]
	pos = int(row["lead_pos"])
	# Defer intermediate writes until the existing completed tree is safe.

	preserve = os.environ.get("PHYML_REPLACE") != "TRUE" and protected_result(loc / "haplotypes.phy")
	pending_tables = {}
	haploid = ch == "X" and a.x_male_only
	ident = dict(
		locus_id=lid,
		index_snp=row["index_snp"],
		source_build=row["source_build"],
		source_pos=row["source_pos"],
		dataset_id=a.dataset,
		genome_build="GRCh37",
		chr=ch,
		lead_pos=pos,
		p_j=row["p_j"],
		beta_j=row["beta_j"],
		effect_allele=row["effect_allele"],
		ld_population="1KG EUR",
		ld_r2_threshold=0.98,
		selection_method="gwas_lead_ld_core",
		analysis_start=int(row["search_start"]),
		analysis_end=int(row["search_end"]),
		core_start=pos - 1,
		core_end=pos,
		lineage="Neanderthal",
		tree_pass=0,
		tree_bootstrap=None,
		status="not_evaluable",
		reason="",
		risk_allele="",
		n_ld_sites=None,
		n_sites=None,
		n_risk_copies=None,
		n_candidate_haplotypes=None,
		n_candidate_copies=None,
	)
	locus = dict(
		chrom=ch, name=lid, start=int(row["search_start"]), end=int(row["search_end"]), core_start=pos - 1, core_end=pos
	)
	haps = []
	detail = []
	copyrows = []
	allcopies = []
	tree = dict(
		locus_id=lid,
		tree_status="not_run",
		candidate_lineage="Neanderthal",
		expected_lineage="Neanderthal",
		candidate_clade_pass=0,
		candidate_clade_bootstrap=None,
		tree_call_reason="not_run",
		candidate_clade_tips="",
		candidate_tips_in_clade="",
		expected_archaic_tips_in_clade="",
		control_tips_in_clade="",
		candidate_context_tips_in_clade="",
		n_candidate_tips_in_clade=0,
		n_expected_archaic_tips_in_clade=0,
		candidate_clade_n_tips=0,
		candidate_clade_modern_tips=0,
		candidate_clade_archaic_tips=0,
		candidate_clade_specificity=None,
		tree_scope="gwas_risk_core",
		tree_newick="",
		tree_file="",
		stats_file="",
		plot_file="",
		tree_has_ancestral_outgroup=0,
	)
	failed = False
	try:
		vcf = vcf_path(a.vcf_dir, ch)
		contig = vcf_contig(vcf, ch)
		panel = read(a.sample_file)
		meta = {r["sample"]: r for r in panel}
		samples, records = query_rows(vcf, f"{contig}:{pos}-{pos}")
		if any(s not in meta for s in samples):
			raise ValueError("sample panel must include every VCF sample and population")
		eur = [
			i
			for i, s in enumerate(samples)
			if meta[s].get("super_pop", "") == "EUR" or meta[s].get("pop", "") in {"CEU", "GBR", "FIN", "IBS", "TSI"}
		]
		if not eur:
			raise ValueError("No EUR samples in --sample-panel; LD cannot use all populations as a fallback")
		lead = exact_lead(row, records, len(samples), haploid)
		ident.update(
			risk_allele=lead["risk_allele"],
			target_variant=f"{ch}:{pos}:{lead['ref']}:{lead['alt']}",
			n_EUR_individuals=len(eur),
			n_risk_copies=sum(x == 1 for x in lead["copies"]),
		)
		indexes = [2 * i + h for i in eur for h in range(1 if haploid else 2)]
		eur_called = [lead["copies"][i] for i in indexes if lead["copies"][i] >= 0]
		ident["risk_frequency_EUR"] = sum(eur_called) / len(eur_called) if eur_called else None
		_, sites, _ = modern_data(
			vcf, contig, locus, 2, read_sexes(a.sample_file), a.x_male_only, a.x_par_diploid, ancestral_any_base=True
		)
		# Multiple records at one coordinate cannot be matched to an archaic
		# base unambiguously; exclude them rather than merge unrelated alleles.
		counts = Counter(s["pos"] for s in sites)
		sites = [s for s in sites if counts[s["pos"]] == 1]
		core, ld, start, end = define_core(sites, lead, indexes)
		pending_tables["ld.tsv"] = ld
		ident.update(
			core_start=start,
			core_end=end,
			core_kb=(end - start) / 1000,
			n_ld_sites=sum(r["core_marker"] for r in ld),
			n_search_sites=len(sites),
			anchor_pos=pos,
			selected_start=start,
			selected_end=end,
		)
		high = [r for r in ld if r["core_marker"]]
		# Search limits are explicit: the inferred span is conditional on this
		# finite window, and a marker near an edge requests a wider search.
		ident["search_edge_warning"] = int(start - locus["start"] < 10000 or locus["end"] - end < 10000)
		calls = {}
		for ref in REFS:
			av = archaic_vcf(a.archaic_root, ref, ch)
			calls[ref] = archaic_calls(
				av, vcf_contig(av, ch), dict(locus, start=start, end=end), core, allow_third_allele=True
			)
		markers = {r["pos"]: r["risk_linked_allele"] for r in high}
		for ref in REFS:
			called = [(p, b) for p, b in markers.items() if calls[ref].get(p, "N") in BASES]
			ident[ref + "_LD_matches"] = sum(calls[ref][p] == b for p, b in called)
			ident[ref + "_LD_called"] = len(called)
		sites, haps, grouped = prepare_sequences(core, calls, samples, lead, haploid)
		ident.update(
			n_sites=len(sites),
			n_compared=len(sites),
			n_ancestral_sites=sum(s["ancestral"] in BASES for s in sites),
			n_tree_haplotypes=len(haps),
			n_singleton_copies=sum(len(ix) for ix in grouped.values() if len(ix) == 1),
			n_candidate_haplotypes=sum(h["role"] == "risk" for h in haps),
			n_candidate_copies=sum(h["n"] for h in haps if h["role"] == "risk"),
			n_nonrisk_haplotypes=sum(h["role"] == "nonrisk" for h in haps),
		)
		if len(sites) < 2:
			raise SkipLocus("insufficient five-reference callable core sites", "insufficient_tree_sites")
		arch = {r: "".join(calls[r][s["pos"]] for s in sites) for r in REFS}
		ancestor = "".join(s["ancestral"] for s in sites)
		pending_tables["sites.tsv"] = [
			dict(chr=ch, pos=s["pos"], id=s["vid"], ref=s["ref"], alt=s["alt"]) for s in sites
		]
		pending_tables["archaic.tsv"] = [dict(archaic=r, lineage=reference_lineage(r), seq=arch[r]) for r in REFS]
		pending_tables["ancestral.tsv"] = [
			dict(reference="Ancestral", n_callable=ident["n_ancestral_sites"], seq=ancestor)
		]
		for h in haps:
			h["copies"] = ";".join(f"{samples[i // 2]}:{i % 2 + 1}" for i in h["indices"])
			comp = []
			for ref, seq in arch.items():
				pairs = [(x, y) for x, y in zip(h["seq"], seq) if x in BASES and y in BASES]
				nm = sum(x == y for x, y in pairs)
				nc = len(pairs)
				comp.append((nm / nc if nc else 0, nc, nm, ref))
			prop, nc, nm, ref = max(comp)
			d = dict(
				ident,
				hap_id=h["hap_id"],
				role=h["role"],
				n_copies=h["n"],
				n_individuals=len({i // 2 for i in h["indices"]}),
				n_compared=nc,
				n_match=nm,
				archaic=ref,
				prop_match=prop,
				candidate_start=start,
				candidate_end=end,
				call="risk_haplotype"
				if h["role"] == "risk"
				else "nonrisk_control"
				if h["role"] == "nonrisk"
				else "risk_nonrisk_sequence_unresolved",
				tree_pass=0,
			)
			d["superpopulation_copy_counts"] = ",".join(
				f"{p}:{n}"
				for p, n in sorted(
					Counter(meta[samples[i // 2]].get("super_pop", "UNKNOWN") for i in h["indices"]).items()
				)
			)
			detail.append(d)
			for i in h["indices"]:
				cp = dict(
					locus_id=lid,
					hap_id=h["hap_id"],
					sample=samples[i // 2],
					sample_id=samples[i // 2],
					haplotype=i % 2 + 1,
					candidate_start=start,
					candidate_end=end,
					role=h["role"],
				)
				allcopies.append(cp)
				if h["role"] == "risk":
					copyrows.append(cp)
			h.update(
				locus_id=lid,
				genome_build="GRCh37",
				best_archaic=ref,
				best_lineage=reference_lineage(ref),
				n_compared=nc,
				n_match=nm,
				prop_match=prop,
				direct_match_pass=0,
			)
		pending_tables["haplotypes.tsv"] = [{k: v for k, v in h.items() if k != "indices"} for h in haps]
		if any(h["role"] == "mixed" for h in haps):
			ident.update(n_candidate_haplotypes=None, n_candidate_copies=None)
			raise SkipLocus("lead alleles share identical callable core sequences", "risk_nonrisk_sequence_unresolved")
		if not ident["n_candidate_haplotypes"] or not ident["n_nonrisk_haplotypes"]:
			raise SkipLocus("both recurrent risk and nonrisk haplotypes required", "insufficient_recurrent_haplotypes")
		if not ident["n_ancestral_sites"]:
			raise SkipLocus("ancestral human bases absent; REF is not an outgroup", "ancestral_sequence_unavailable")
		phy = loc / "haplotypes.phy"
		tree["phy_file"] = str(phy)
		seqs = [(h["hap_id"], h["seq"]) for h in haps] + list(arch.items()) + [("Ancestral", ancestor)]
		phytext = f"{len(seqs)} {len(sites)}\n" + "".join(f"{label:<10} {seq}\n" for label, seq in seqs)
		prepare_input(phy, phytext, replace=os.environ.get("PHYML_REPLACE") == "TRUE")
		pending_tables["haplotypes.phy.meta.tsv"] = [
			dict(
				phy_label=label,
				label=label,
				role=next(
					(h["role"] for h in haps if h["hap_id"] == label),
					"ancestral" if label == "Ancestral" else "archaic",
				),
			)
			for label, seq in seqs
		]
		ident.update(status="tree_not_requested", reason="sequence_prepared")
		if a.action == "run" and a.plot_phy == "TRUE":
			print(
				f"[GU PHYML] {lid}: EUR={len(eur)}; core={ch}:{start + 1}-{end}; LD markers={ident['n_ld_sites']}; tree sites={len(sites)}; recurrent haplotypes={len(haps)}; bootstrap=100",
				flush=True,
			)
			cmd = [
				sys.executable,
				str(Path(__file__).with_name("phyml.py")),
				"run",
				"--phy",
				str(phy),
				"--scope",
				"gwas_risk_core",
				"--bootstrap",
				"100",
				"--timeout",
				str(a.timeout),
				"--cpus",
				str(a.cpus),
				"--mpi-fallback",
				"serial",
			]
			proc = run_verified_tree(cmd)
			treepath = Path(str(phy) + "_phyml_tree.txt")
			statpath = Path(str(phy) + "_phyml_stats.txt")
			if proc.returncode:
				if preserve:
					require_replace(phy, "existing tree could not be verified; see runner diagnostic")
				failed = True
				logpath = Path(str(phy) + ".phyml.log")
				failure_log = logpath.read_text(errors="replace") if logpath.exists() else ""
				reason = (
					"tree_timeout"
					if proc.returncode == 124
					else "numerical_model_fit_failed"
					if "Cannot work out eigen vectors" in failure_log
					else "see_PhyML_run_log"
				)
				ident.update(status="tree_failed", reason=reason)
				tree.update(tree_status="failed", tree_call_reason=reason)
			else:
				newick = treepath.read_text().strip()
				risk = {h["hap_id"] for h in haps if h["role"] == "risk"}
				match = risk_clade(newick, risk, {h["hap_id"] for h in haps})
				bs = match["bootstrap"] if match else None
				passed = bool(bs is not None and bs >= 70)
				reason = "risk_Neanderthal_split" if match else "no_exclusive_risk_Neanderthal_split"
				ident.update(
					status="tree_supported" if passed else "tree_not_supported",
					reason=reason,
					tree_pass=int(passed),
					tree_bootstrap=bs,
				)
				tree.update(
					tree_status="complete",
					tree_newick=newick,
					tree_file=str(treepath),
					stats_file=str(statpath),
					tree_has_ancestral_outgroup=1,
					candidate_clade_pass=int(passed),
					candidate_clade_bootstrap=bs,
					candidate_clade_tips=match["tips"] if match else "",
					candidate_tips_in_clade=",".join(sorted(risk)) if match else "",
					candidate_clade_rule="all_recurrent_risk_plus_three_Neanderthals_no_nonrisk_or_ancestor",
					tree_call_reason=reason,
					n_bootstrap_nodes=len(bootstrap_values(newick)),
					expected_archaic_tips_in_clade=",".join(REFS) if match else "",
					n_candidate_tips_in_clade=len(risk) if match else 0,
					n_expected_archaic_tips_in_clade=3 if match else 0,
					candidate_clade_n_tips=len(risk) + 3 if match else 0,
					candidate_clade_modern_tips=len(risk) if match else 0,
					candidate_clade_archaic_tips=3 if match else 0,
					candidate_clade_specificity=1 if match else None,
				)
		tree["tree_call_reason"] = ident["reason"]
	except SkipLocus as e:
		if preserve:
			require_replace(loc / "haplotypes.phy", f"new analysis would skip completed locus: {e.code}")
		ident.update(status=e.code, reason=str(e))
		tree["tree_call_reason"] = str(e)
		if e.code == "risk_nonrisk_sequence_unresolved":
			copyrows = []
		if e.details and "ld" in e.details:
			pending_tables["ld.tsv"] = e.details["ld"]
			ident.update(n_ld_sites=sum(r["core_marker"] for r in e.details["ld"]), n_search_sites=len(e.details["ld"]))
		print(f"[GU PHYML] SKIP {lid}: {e.code}: {e}", flush=True)
	for name in ("sites.tsv", "archaic.tsv", "ancestral.tsv", "haplotypes.tsv", "ld.tsv"):
		(loc / name).unlink(missing_ok=True)
	for name, table in pending_tables.items():
		write(loc / name, table)
	ident["call"] = ident["status"]
	# A length-model sensitivity statistic, not a calibrated locus-specific P.
	if ident.get("core_kb"):
		ident["ils_probability"] = ils_probability(ident["core_end"] - ident["core_start"])
		ident["ils_model"] = "assumed_0.53cM/Mb_29yr_550k_split_50k_archaic_age;not_local_map;uncorrected"
	summaries, lineage_trees, lineage_details = lineage_results(ident, tree, detail, haps, arch if detail else {})
	write(final / "gwas_loci.tsv", summaries)
	write(final / "gwas_haplotypes.tsv", lineage_details)
	write(final / "gwas_copies.tsv", [dict(c, lineage=lineage) for lineage in LINEAGE_REFS for c in copyrows])
	# Keep one locus row for database/browser identity; retain both tests in
	# the evidence/report tables. Prefer a supported lineage as representative.
	representative = max(range(len(summaries)), key=lambda i: summaries[i]["tree_pass"])
	write(final / "loci.tsv", [summaries[representative]])
	write(final / "trees.tsv", [lineage_trees[representative]])
	write(final / "evidence_trees.tsv", lineage_trees)
	completed_haps = [{k: v for k, v in h.items() if k != "indices"} for h in haps if "copies" in h]
	write(
		final / "haplotypes.tsv",
		completed_haps,
		list(completed_haps[0]) if completed_haps else ["locus_id", "hap_id", "copies"],
	)
	write(
		final / "haplotype_samples.tsv",
		allcopies,
		["locus_id", "hap_id", "sample", "sample_id", "haplotype", "candidate_start", "candidate_end", "role"],
	)
	write(final / "skipped_loci.tsv", [] if ident["status"].startswith("tree_") else [ident])
	write(final / "gwas_lead.tsv", [row])
	parameters = dict(
		workflow=WORKFLOW,
		ld_population="1KG EUR",
		ld_rule="phased_r2 > 0.98",
		lead=row,
		tree_populations="all target 1KG samples",
		minimum_haplotype_copies=2,
		minimum_minor_allele_copies=2,
		ancestral_source="target VCF INFO/AA; unknown remains N; not verified as Ensembl release 100",
		references=list(REFS),
		heterozygous_archaic_policy="mask_as_N",
		bootstrap=100,
		model="HKY85+G4+I;estimated",
		lineage_references=LINEAGE_REFS,
		tree_rule="one five-reference tree; separately test all recurrent risk + each lineage, excluding other archaics, nonrisk and ancestor; BS >=70",
		method_source="https://www.nature.com/articles/s41586-020-2818-3",
	)
	(final / "gwas_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
	if tree["tree_status"] == "complete":
		plot = subprocess.run(["Rscript", "--vanilla", str(Path(__file__).with_name("phyml.R")), "--out", str(out)])
		if plot.returncode:
			print("[GU PHYML] WARNING: tree completed but plot export failed", flush=True)
			failed = True
	print(
		f"[GU PHYML] {lid}: " + "; ".join(f"{s['lineage']}={s['status']} ({s['reason']})" for s in summaries),
		flush=True,
	)
	return 1 if failed else 0


def lineage_results(ident, tree, details, haps, arch):
	"""Report both prespecified lineage tests from the same inferred tree."""
	summaries, trees, rows = [], [], []
	risk = {h["hap_id"] for h in haps if h["role"] == "risk"}
	modern = {h["hap_id"] for h in haps}
	sequences = {h["hap_id"]: h["seq"] for h in haps}
	for lineage, refs in LINEAGE_REFS.items():
		s = dict(ident, lineage=lineage)
		t = dict(tree, candidate_lineage=lineage, expected_lineage=lineage)
		if tree["tree_status"] == "complete":
			match = risk_clade(tree["tree_newick"], risk, modern, refs)
			bs = match["bootstrap"] if match else None
			passed = bool(bs is not None and bs >= 70)
			reason = f"risk_{lineage}_split" if match else f"no_exclusive_risk_{lineage}_split"
			s.update(
				status="tree_supported" if passed else "tree_not_supported",
				reason=reason,
				tree_pass=int(passed),
				tree_bootstrap=bs,
			)
			t.update(
				candidate_clade_pass=int(passed),
				candidate_clade_bootstrap=bs,
				candidate_clade_rule="all_recurrent_risk_plus_lineage_no_other_archaics_nonrisk_or_ancestor",
				tree_call_reason=reason,
				candidate_clade_tips=match["tips"] if match else "",
				candidate_tips_in_clade=",".join(sorted(risk)) if match else "",
				expected_archaic_tips_in_clade=",".join(refs) if match else "",
				n_candidate_tips_in_clade=len(risk) if match else 0,
				n_expected_archaic_tips_in_clade=len(refs) if match else 0,
				candidate_clade_n_tips=len(risk) + len(refs) if match else 0,
				candidate_clade_modern_tips=len(risk) if match else 0,
				candidate_clade_archaic_tips=len(refs) if match else 0,
				candidate_clade_specificity=1 if match else None,
			)
		s["call"] = s["status"]
		if lineage == "Denisovan":
			s.update(ils_probability=None, ils_model="not_parameterized_for_Denisovan")
		summaries.append(s)
		trees.append(t)
		for detail in details:
			d = dict(
				detail,
				lineage=lineage,
				tree_bootstrap=s["tree_bootstrap"],
				tree_pass=int(s["tree_pass"] and detail["role"] == "risk"),
				call=s["status"] if detail["role"] == "risk" else detail["call"],
			)
			comparisons = []
			for ref in refs:
				pairs = [(x, y) for x, y in zip(sequences[d["hap_id"]], arch[ref]) if x in BASES and y in BASES]
				nc = len(pairs)
				nm = sum(x == y for x, y in pairs)
				comparisons.append((nm / nc if nc else 0, nc, nm, ref))
			prop, nc, nm, ref = max(comparisons)
			d.update(archaic=ref, n_compared=nc, n_match=nm, prop_match=prop)
			rows.append(d)
	return summaries, trees, rows


def phyml_gwas_main():
	p = argparse.ArgumentParser(description=__doc__)
	for name in ("lead-table", "loci", "vcf-dir", "archaic-root", "sample-file", "out"):
		p.add_argument("--" + name, type=Path, required=True)
	p.add_argument("--dataset", default="1kg")
	p.add_argument("--action", choices=["run", "match", "check"], default="run")
	p.add_argument("--plot-phy", choices=["TRUE", "FALSE"], default="TRUE")
	p.add_argument("--timeout", type=int, default=86400)
	p.add_argument("--cpus", type=int, default=4)
	p.add_argument("--x-male-only", action="store_true")
	p.add_argument("--x-par-diploid", action="store_true")
	a = p.parse_args()
	if a.dataset != "1kg":
		raise ValueError(
			"This EUR LD workflow currently requires target 1kg; another target needs a separate fixed 1KG EUR LD reference"
		)
	names = [line.split()[3] for line in a.loci.read_text().splitlines() if line.strip() and not line.startswith("#")]
	rows = [r for r in read(a.lead_table) if r["locus_id"] in names]
	if len(rows) != 1 or len(names) != 1:
		raise ValueError("Each worker requires exactly one original GWAS lead")
	if a.action == "check":
		print("[GU PHYML] original lead, EUR LD and recurrent-haplotype workflow configured")
		return
	a.out.mkdir(parents=True, exist_ok=True)
	with (a.out / ".gwas.lock").open("w") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		raise SystemExit(run_locus(a, rows[0]))


def phyml_gwas_cli():
	try:
		phyml_gwas_main()
	except (OSError, RuntimeError, ValueError) as e:
		print(f"ERROR: {e}", file=sys.stderr)
		sys.exit(1)


# 🚩 phyml_report
"""Report predefined GWAS risk cores and cross-check their carriers with IBDmix."""
import argparse, hashlib, json, os, sqlite3
from pathlib import Path
from collections import defaultdict
from datetime import datetime, timezone

dl = load_module("review.py")
load_module("review.py")
from gu_review import write_tsv, atomic_text, coverage, union


def local_path(value):
	value = value or ""
	if os.name == "nt" and str(value).startswith("/mnt/"):
		return Path(str(value)[5] + ":" + str(value)[6:])
	return Path(value)


def validate(copies, con, threshold=0.8):
	"""Same sample and actual tract; TRACE supports unknown ancestry, not lineage."""
	con.row_factory = sqlite3.Row
	runs = {
		(r["dataset_id"], r["genome_build"], dl.dual_lead_chrom(r["chr"]), r["method"]): dict(r)
		for r in con.execute("SELECT * FROM method_runs")
	}
	for unit, run in runs.items():
		path = local_path(run.get("raw_file", ""))
		if not path.is_file():
			continue
		text = path.read_text(errors="replace").replace("\\t", "\t")
		# A chromosome-level completion row must not turn an uncovered locus
		# into a negative when only a subset of that chromosome was analyzed.
		if any(line.startswith(("loci_file\t", "loci\t")) for line in text.splitlines()):
			bed = path.parent / "request.loci.analysis.bed"
			run["_scope"] = []
			if bed.is_file():
				for line in bed.read_text().splitlines():
					fields = line.split()
					if len(fields) >= 3 and dl.dual_lead_chrom(fields[0]) == unit[2]:
						run["_scope"].append((int(fields[1]), int(fields[2])))
		if unit[3] == "ibdmix":
			meta = dict(line.split("\t", 1) for line in text.splitlines() if "\t" in line)
			refs = meta.get("refs", "").split()
			run["_tested_lineages"] = set()
			if set(refs) & {"Altai", "Vindija", "Chagyr", "Chagyrskaya"}:
				run["_tested_lineages"].add("Neanderthal")
			if any("denis" in ref.lower() for ref in refs) and (
				meta.get("background_filter") == "0" or meta.get("export_denisovan") == "1"
			):
				run["_tested_lineages"].add("Denisovan")
			chrom = unit[2]
			scope = ("X_MALE" if "male_haploid_nonpar" in text else "X_PAR") if chrom == "X" else "C" + chrom
			for roster in (path.parent / "samples" / scope / "ALL.txt", path.parent / "samples" / "ALL.txt"):
				if roster.is_file():
					run["_tested_samples"] = set(roster.read_text().split())
					break
		if unit[3] == "trace":
			sm = path.parent / "samples/trace_sample_map.tsv"
			if sm.is_file():
				nodes = path.parent / "samples/tree_nodes.txt"
				chosen = set(nodes.read_text().split()) if nodes.is_file() else None
				run["_tested_copies"] = {
					str(r["sample"]) + ":" + str(r["haplotype"])
					for r in dl.rows(sm)
					if chosen is None or str(r["tree_node_id"]) in chosen
				}
	windows = defaultdict(list)
	for c in copies:
		windows[(c["dataset_id"], c["genome_build"], c["chr"])].append((c["candidate_start"], c["candidate_end"]))
	segments = defaultdict(list)
	for unit, ww in windows.items():
		seen = set()
		for st, en in union(ww):
			for r in con.execute(
				"SELECT * FROM segments WHERE dataset_id=? AND genome_build=? AND chr=? AND method IN ('ibdmix','trace') AND start<? AND end>?",
				(*unit, en, st),
			):
				r = dict(r)
				key = tuple(r.values())
				if key in seen:
					continue
				seen.add(key)
				segments[(*unit, r["sample_id"], r["method"])].append(r)
	result = []
	for c in copies:
		unit = (c["dataset_id"], c["genome_build"], c["chr"])
		interval = (c["candidate_start"], c["candidate_end"])
		for method in ["ibdmix", "trace"]:
			run = runs.get((*unit, method), {})
			complete = run.get("status") == "complete"
			in_scope = "_scope" not in run or coverage(interval, run["_scope"])["union_fraction"] >= 1 - 1e-9
			in_panel = (
				("_tested_copies" not in run or c["sample_id"] + ":" + str(c["haplotype"]) in run["_tested_copies"])
				and ("_tested_samples" not in run or c["sample_id"] in run["_tested_samples"])
				and (method != "ibdmix" or c["lineage"] in run.get("_tested_lineages", set()))
			)
			compared = complete and in_scope and in_panel
			eligible = compared and dl.dual_lead_truth(run.get("evidence_eligible"))
			pool = segments.get((*unit, c["sample_id"], method), [])
			pool = [s for s in pool if method == "trace" or s["source_class"] == c["lineage"]]
			samephase = [s for s in pool if str(s.get("haplotype")) == str(c["haplotype"])] if method == "trace" else []
			m = coverage(interval, [(s["start"], s["end"]) for s in pool])
			p = coverage(interval, [(s["start"], s["end"]) for s in samephase])
			result.append(
				dict(
					c,
					method=method,
					method_complete=int(complete),
					comparison_available=int(compared),
					evidence_eligible=int(eligible),
					scope_status="in_scope" if in_scope else "outside_or_unknown_run_scope",
					panel_status="included_or_panel_not_recorded" if in_panel else "not_in_method_sample_map",
					availability_note=run.get("availability_note", "no_completed_run"),
					coverage_threshold=threshold,
					overlap_fraction=m["best_single_fraction"] if complete or pool else None,
					union_fraction=m["union_fraction"] if complete or pool else None,
					overlap_pass=int(m["best_single_fraction"] >= threshold) if compared else None,
					phase_overlap_pass=int(p["best_single_fraction"] >= threshold)
					if compared and method == "trace"
					else None,
					phase_note="same_sample_only_unphased"
					if method == "ibdmix"
					else "same_sample_and_stored_haplotype_index",
					lineage_note="same_lineage" if method == "ibdmix" else "ghost_unknown_not_lineage_validation",
					callability_note="individual_callable_bases_not_available",
				)
			)
	return result, runs


def attach_validation(rows, validation, runs, by_hap=False):
	grouped = defaultdict(list)
	for v in validation:
		grouped[v["candidate_id"] if by_hap else (v["locus_key"], v["lineage"])].append(v)
	for row in rows:
		rr = grouped[row.get("candidate_id") if by_hap else (row["locus_key"], row["lineage"])]
		for method in ["ibdmix", "trace"]:
			run = runs.get((row["dataset_id"], row["genome_build"], row["chr"], method), {})
			candidates = [v for v in rr if v["method"] == method]
			vv = [v for v in candidates if v.get("comparison_available", 1)]
			complete = run.get("status") == "complete"
			seen = {v["sample_id"] for v in vv}
			passing = {v["sample_id"] for v in vv if v["overlap_pass"] == 1}
			any_overlap = {v["sample_id"] for v in vv if (v.get("overlap_fraction") or 0) > 0}
			status = (
				"not_run"
				if not complete
				else "no_candidate_to_test"
				if not candidates
				else "not_evaluable_scope_or_panel"
				if not vv
				else "exploratory"
				if not dl.dual_lead_truth(run.get("evidence_eligible"))
				else "overlap_detected"
				if passing
				else "partial_overlap"
				if any_overlap
				else "not_detected"
			)
			row.update(
				{
					method + "_status": status,
					method + "_individuals": len(seen) if complete and vv else None,
					method + "_unassessed_copies": len(candidates) - len(vv) if complete else None,
					method + "_any_overlap_individuals": len(any_overlap) if complete and vv else None,
					method + "_supported_individuals": len(passing) if complete and vv else None,
					method + "_support_fraction": len(passing) / len(seen) if complete and seen else None,
					method + "_supported_copies": sum(v["phase_overlap_pass"] == 1 for v in vv)
					if method == "trace" and complete and vv
					else None,
					method + "_candidate_copies": len(vv) if complete and vv else None,
				}
			)


def gwas_rows(path):
	summaries = read(path)
	details = read(path.with_name("gwas_haplotypes.tsv"))
	copies = read(path.with_name("gwas_copies.tsv"))
	summaries = [s for s in summaries if s.get("locus_id")]
	result = []
	all_details = []
	all_copies = []
	for s in summaries:
		s = dict(s)
		lineage = s.get("lineage") or "Neanderthal"
		if lineage == "Denisova":
			lineage = "Denisovan"
		for key in (
			"core_start",
			"core_end",
			"lead_pos",
			"tree_pass",
			"n_candidate_haplotypes",
			"n_candidate_copies",
			"n_sites",
			"n_ld_sites",
		):
			if s.get(key) not in (None, ""):
				s[key] = int(s[key])
		key = hashlib.sha256(
			"|".join(
				str(s.get(k, "")) for k in ("dataset_id", "genome_build", "locus_id", "core_start", "core_end")
			).encode()
		).hexdigest()[:20]
		s.update(
			locus_key=key,
			lineage=lineage,
			core_interval=f"{s['core_start'] + 1}–{s['core_end']}" if s.get("core_kb") else "未定义",
			ld_rule="> 0.98 (EUR)",
			risk_haplotypes=f"{s['n_candidate_haplotypes']} / {s['n_candidate_copies']}"
			if s.get("n_candidate_haplotypes") not in (None, "")
			else None,
		)
		result.append(s)
		for d in details:
			if d.get("locus_id") != s["locus_id"] or d.get("lineage", "Neanderthal") != lineage:
				continue
			d = dict(d, lineage=lineage, locus_key=key, candidate_id=f"{key}|{lineage}|{d['hap_id']}")
			all_details.append(d)
			if d["role"] == "risk":
				for c in copies:
					if (
						c.get("locus_id") != s["locus_id"]
						or c.get("hap_id") != d["hap_id"]
						or c.get("lineage", "Neanderthal") != lineage
					):
						continue
					all_copies.append(
						dict(
							c,
							dataset_id=s["dataset_id"],
							genome_build=s["genome_build"],
							chr=s["chr"],
							locus_key=key,
							lineage=lineage,
							candidate_id=d["candidate_id"],
							candidate_start=s["core_start"],
							candidate_end=s["core_end"],
						)
					)
	return result, all_details, all_copies


def phyml_report_main():
	ap = argparse.ArgumentParser(description=__doc__)
	ap.add_argument("--normalize", required=True, type=Path)
	ap.add_argument("--database", required=True, type=Path)
	ap.add_argument("--output", required=True, type=Path)
	ap.add_argument("--loci", type=Path)
	args = ap.parse_args()
	summary = []
	details = []
	copies = []
	files = []
	trees = []
	for path in sorted((args.normalize / "phyml").glob("**/final/gwas_loci.tsv")):
		ss, dd, cc = gwas_rows(path)
		summary += ss
		details += dd
		copies += cc
		files.extend(path.parent.glob("gwas_*.tsv"))
		trees += read(path.with_name("evidence_trees.tsv"))
	with sqlite3.connect(args.database.resolve().as_uri() + "?mode=ro", uri=True) as con:
		validation, runs = validate(copies, con)
	attach_validation(summary, validation, runs)
	attach_validation(details, validation, runs, True)
	out = args.output
	out.mkdir(parents=True, exist_ok=True)
	for name, rr in [
		("phyml_locus_report", summary),
		("phyml_haplotype_report", details),
		("phyml_copy_validation", validation),
		("phyml_lineage_trees", trees),
	]:
		write_tsv(out / (name + ".tsv"), rr)
	# This generated report no longer has two competing lead roles.
	(out / "phyml_lead_report.tsv").unlink(missing_ok=True)
	atomic_text(
		out / "phyml_report_manifest.json",
		json.dumps(
			dict(
				created_utc=datetime.now(timezone.utc).isoformat(),
				workflow="gwas_lead_ld_core",
				n_loci=len({r["locus_key"] for r in summary}),
				n_lineage_tests=len(summary),
				n_haplotype_rows=len(details),
				ld_population="1KG EUR",
				ld_rule="r2 > 0.98",
				risk_definition="COJO refA and bJ; original lead retained",
				coverage_threshold=0.8,
				validation_denominator="all recurrent risk-haplotype carriers, regardless of tree support; individuals counted once",
				tree_support_note="per-lineage predefined risk+archaic split, bootstrap >=70; not a high-confidence introgression classification",
				input_files=[
					dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(set(files))
				],
			),
			indent=2,
		),
	)
	print(
		f"PhyML GWAS report ready: {len({r['locus_key'] for r in summary})} original leads, {len(summary)} lineage tests, {len(details)} recurrent haplotype rows",
		flush=True,
	)


def phyml_report_cli():
	phyml_report_main()


SUBCOMMANDS = {
	"input": phyml_gwas_input_cli,
	"run": phyml_run_cli,
	"cache": phyml_locus_cache_cli,
	"gwas": phyml_gwas_cli,
	"report": phyml_report_cli,
}


def main():
	import sys

	if len(sys.argv) > 1 and sys.argv[1] in SUBCOMMANDS:
		command = sys.argv.pop(1)
		return SUBCOMMANDS[command]()
	if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help", "help"):
		print(__doc__)
		print("Commands: " + ", ".join(SUBCOMMANDS))
		return
	raise SystemExit("Unknown command: " + sys.argv[1])


if __name__ == "__main__":
	main()
