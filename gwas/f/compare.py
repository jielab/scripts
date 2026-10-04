#!/usr/bin/env python3
"""compare workflow utilities; use --help for commands."""

from __future__ import annotations


# 🚩 gwas_ldsc_summary
"""Summarize only the requested run's LDSC logs, retaining unavailable traits/pairs."""
import argparse
import csv
import itertools
import math
from pathlib import Path
import re
import subprocess

RG_COLUMNS = "p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se".split()


def read_rg(path):
	text = Path(path).read_text()
	if "Analysis finished at" not in text:
		raise ValueError("Incomplete LDSC log: " + str(path))
	marker = "Summary of Genetic Correlation Results"
	if marker not in text:
		raise ValueError("Missing rg summary: " + str(path))
	lines = text.rsplit(marker, 1)[1].strip().splitlines()
	header = lines[0].split()
	if header != RG_COLUMNS:
		raise ValueError("Unexpected rg columns: " + str(path))
	rows = []
	for line in lines[1:]:
		if not line.strip():
			break
		values = line.split()
		if len(values) != len(header):
			raise ValueError("Malformed rg summary: " + line)
		row = dict(zip(header, values))
		for key in ("p1", "p2"):
			row[key] = Path(row[key]).name.removesuffix(".sumstats.gz")
		rows.append(row)
	return rows


def read_h2(path):
	text = Path(path).read_text()
	if "Analysis finished at" not in text:
		raise ValueError("Incomplete LDSC log: " + str(path))
	fields = {}
	for name, pattern in [("h2", "Total Observed scale h2"), ("intercept", "Intercept"), ("ratio", "Ratio")]:
		match = re.search(re.escape(pattern) + r":\s*([-+\deE.]+)\s*\(([-+\deE.]+)\)", text)
		fields[name], fields[name + "_se"] = match.groups() if match else ("NA", "NA")
	for name, pattern in [("lambda_gc", "Lambda GC"), ("mean_chi2", r"Mean Chi\^2")]:
		match = re.search(pattern + r":\s*([-+\deE.]+)", text)
		fields[name] = match.group(1) if match else "NA"
	if not finite(fields["h2"]) or not finite(fields["h2_se"]):
		raise ValueError("Missing/non-finite h2 estimate: " + str(path))
	return fields


def finite(value):
	try:
		return math.isfinite(float(value))
	except (ValueError, TypeError):
		return False


def write_tsv(path, columns, rows):
	tmp = path.with_suffix(".tmp")
	with tmp.open("w") as handle:
		writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", restval="NA")
		writer.writeheader()
		writer.writerows(rows)
	tmp.replace(path)


def summarize(out, run_h2=True, run_rg=True):
	out = Path(out)
	with (out / "inputs.status.tsv").open() as handle:
		statuses = list(csv.DictReader(handle, delimiter="\t"))
	traits = [Path(row["FILE"]).name.removesuffix(".gz") for row in statuses]
	available = [trait for trait, row in zip(traits, statuses) if row["STATUS"] == "SCHEDULED"]
	reasons = {trait: row["REASON"] for trait, row in zip(traits, statuses) if row["STATUS"] != "SCHEDULED"}
	h2 = []
	for trait in traits:
		row = dict(
			trait=trait, scope="1-22", scale="observed", status="UNAVAILABLE", reason=reasons.get(trait, "h2 disabled")
		)
		if run_h2 and trait in available:
			row.update(read_h2(out / "h2.log" / (trait + ".h2.log")))
			row.update(status="ESTIMATED", reason="Standard LDSC; X/Y/MT not estimated, see trait qc chromosome audit")
		h2.append(row)
	write_tsv(
		out / "h2.tsv",
		[
			"trait",
			"h2",
			"h2_se",
			"intercept",
			"intercept_se",
			"ratio",
			"ratio_se",
			"lambda_gc",
			"mean_chi2",
			"scope",
			"scale",
			"status",
			"reason",
		],
		h2,
	)
	pairs = {}
	if run_rg:
		for trait in available[:-1]:
			for row in read_rg(out / "rg.log" / (trait + ".rg.log")):
				key = tuple(sorted((row["p1"], row["p2"])))
				if key in pairs:
					raise ValueError("Duplicate rg pair: " + str(key))
				pairs[key] = row
	rg = []
	for a, b in itertools.combinations(traits, 2):
		key = tuple(sorted((a, b)))
		row = dict(
			p1=a,
			p2=b,
			status="UNAVAILABLE",
			reason="; ".join(reasons[t] for t in (a, b) if t in reasons) or "rg disabled",
		)
		if run_rg and a in available and b in available:
			if key not in pairs:
				raise ValueError("Missing rg pair: " + str(key))
			row.update(pairs[key])
			valid = finite(row["rg"]) and finite(row["p"]) and finite(row["se"])
			row.update(
				status="ESTIMATED" if valid else "UNAVAILABLE",
				reason="Standard LDSC, chromosomes 1-22"
				if valid
				else "LDSC returned a non-finite estimate; see rg log",
			)
		rg.append(row)
	write_tsv(out / "rg.tsv", RG_COLUMNS + ["status", "reason"], rg)
	if run_rg:
		subprocess.run(["Rscript", str(Path(__file__).with_name("compare.R")), "ldsc", str(out)], check=True)
	else:
		(out / "rg.png").unlink(missing_ok=True)
	print(
		f"Summarized {sum(r['status'] == 'ESTIMATED' for r in h2)} h2 estimates and "
		f"{sum(r['status'] == 'ESTIMATED' for r in rg)} rg pairs: {out}",
		flush=True,
	)


def gwas_ldsc_summary_cli():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--output-dir", required=True)
	parser.add_argument("--h2", choices=("True", "False"), default="True")
	parser.add_argument("--rg", choices=("True", "False"), default="True")
	args = parser.parse_args()
	summarize(args.output_dir, args.h2 == "True", args.rg == "True")


# 🚩 gwas_ldsc
"""Run the user's bulik/ldsc checkout in its Python 2 conda environment."""
import argparse
import gzip
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import json
import csv
import hashlib
import tempfile
import fcntl
from contextlib import contextmanager
import sys
from collections import Counter

DEFAULT_REF = "/mnt/f/refLD/ldsc/1000G/1000G_Phase3_ldscores/LDscore."
DEFAULT_WEIGHTS = "/mnt/f/refLD/ldsc/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
DEFAULT_ALLELES = "/mnt/f/refLD/ldsc/hm3/w_hm3.snplist"
AUTOSOMES = {str(c) for c in range(1, 23)}


def validate_references(ref, weights):
	paths = []
	for c in range(1, 23):
		for prefix, suffix in ((ref, ".l2.ldscore.gz"), (ref, ".l2.M_5_50"), (weights, ".l2.ldscore.gz")):
			stem = prefix.replace("@", str(c)) if "@" in prefix else prefix + str(c)
			f = Path(stem + suffix)
			if not f.is_file() or not f.stat().st_size:
				raise ValueError("Missing LDSC reference: " + str(f))
			paths.append(f)
	return paths


def boolean(x):
	if x.upper() not in ("TRUE", "FALSE"):
		raise argparse.ArgumentTypeError("Use TRUE/FALSE")
	return x.upper() == "TRUE"


@contextmanager
def directory_lock(path):
	fd = os.open(path, os.O_RDONLY)
	try:
		fcntl.flock(fd, fcntl.LOCK_EX)
		yield
	finally:
		os.close(fd)


def cache_paths(source):
	source = Path(source)
	trait = source.name.removesuffix(".gz")
	qc = source.parent.parent / "qc" if source.parent.name == "gwas" else source.parent / "qc"
	return (
		source.parent / (trait + ".sumstats.gz"),
		qc / (trait + ".sumstats.cache.tsv"),
		qc / (trait + ".sumstats.chromosomes.tsv"),
	)


def file_identity(path):
	path = Path(path).resolve()
	st = path.stat()
	return [str(path), st.st_size, st.st_mtime_ns]


def validate_sumstats(cache):
	if not cache.is_file():
		raise ValueError("Missing pre-generated sumstats: " + str(cache))
	with gzip.open(cache, "rt") as handle:
		if not {"SNP", "A1", "A2", "Z", "N"}.issubset(handle.readline().split()) or not handle.readline().strip():
			raise ValueError("Empty/invalid pre-generated sumstats: " + str(cache))


def prepare(source, merge_alleles, software, interpreter, fallback_n=None, allow_build="True"):
	source = Path(source).resolve()
	cache, meta, auditfile = cache_paths(source)
	if allow_build != "True":
		validate_sumstats(cache)
		print("Reuse pre-generated sumstats: " + str(cache), flush=True)
		return
	meta.parent.mkdir(parents=True, exist_ok=True)
	# Lock the directory inode, avoiding persistent per-trait lock files.
	with directory_lock(source.parent):
		signature = hashlib.sha256(
			json.dumps(
				dict(
					version=2,
					source=file_identity(source),
					alleles=file_identity(merge_alleles),
					munge=file_identity(Path(software) / "munge_sumstats.py"),
					interpreter=interpreter,
					N=fallback_n,
				),
				sort_keys=True,
			).encode()
		).hexdigest()
		prior = {}
		if meta.exists():
			with meta.open() as handle:
				prior = dict(csv.reader(handle, delimiter="\t"))
		if (
			cache.is_file()
			and auditfile.is_file()
			and prior.get("signature") == signature
			and prior.get("output") == str(file_identity(cache))
		):
			print("Reuse munged GWAS: " + str(cache), flush=True)
			return
		print("Prepare munged GWAS: " + str(cache), flush=True)
		with tempfile.TemporaryDirectory(prefix=".ldsc-", dir=source.parent) as work:
			prefix = Path(work) / source.name.removesuffix(".gz")
			fixed = str(prefix) + ".input.gz"
			fix_p(str(source), fixed, merge_alleles)
			audit = json.loads(Path(fixed + ".chromosomes.json").read_text())
			with gzip.open(fixed, "rt") as handle:
				header = handle.readline().split()
			if "N" not in header and fallback_n is None:
				raise ValueError("Missing N: " + str(source))
			n = ["--N-col", "N"] if "N" in header else ["--N", str(fallback_n)]
			cmd = (
				interpreter
				+ [
					str(Path(software) / "munge_sumstats.py"),
					"--sumstats",
					fixed,
					"--out",
					str(prefix),
					"--merge-alleles",
					merge_alleles,
					"--snp",
					"SNP",
					"--a1",
					"EA",
					"--a2",
					"NEA",
					"--signed-sumstats",
					"BETA,0",
					"--p",
					"P",
					"--chunksize",
					"10000",
				]
				+ n
			)
			log = meta.parent / (source.name.removesuffix(".gz") + ".sumstats.log")
			with log.open("w") as handle:
				handle.write(audit["reason"] + "\n")
				handle.flush()
				subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, check=True)
			munged = Path(str(prefix) + ".sumstats.gz")
			with gzip.open(munged, "rt") as handle:
				if not {"SNP", "A1", "A2", "Z", "N"}.issubset(handle.readline().split()) or not handle.readline():
					raise ValueError("Empty/invalid munged GWAS: " + str(munged))
			with auditfile.open("w") as handle:
				writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
				writer.writerow(["CHR", "INPUT", "KEPT", "REASON"])
				for chrom, count in audit["input"].items():
					writer.writerow(
						[
							chrom,
							count,
							audit["kept"].get(chrom, 0),
							audit["reason"] if chrom not in AUTOSOMES else "HM3 and valid P filter before munge QC",
						]
					)
			os.replace(munged, cache)
			tmp = meta.with_suffix(".tmp")
			with tmp.open("w") as handle:
				csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(
					[("signature", signature), ("output", str(file_identity(cache)))]
				)
			os.replace(tmp, meta)


def gwas_ldsc_main():
	p = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog="""Examples (replace /path/... with your files):
  compare.sh ldsc --gwas-files /path/A.gz,/path/B.gz --output-dir /path/ldsc-results \\
    --merge-alleles /path/w_hm3.snplist \\
    --ref-ld-chr /path/eur_w_ld_chr/ --w-ld-chr /path/eur_w_ld_chr/
  Add --run FALSE to write commands only, --run-rg FALSE for h2 only,
  or --N 10000 if the actual sample size is 10000 and no N column is present.
""",
	)
	p.add_argument("--gwas-files", required=True, help="Comma-separated standardized GWAS files")
	p.add_argument("--output-dir", default="/mnt/d/analysis/gwas/ldsc")
	p.add_argument("--ldsc-software-dir", default="/mnt/d/software/ldsc")
	p.add_argument("--conda", default=shutil.which("conda") or str(Path.home() / "anaconda3/bin/conda"))
	p.add_argument("--conda-env", default="ldsc")
	p.add_argument("--python", help="Explicit LDSC-compatible interpreter; bypass conda")
	p.add_argument("--merge-alleles", default=DEFAULT_ALLELES, help="w_hm3.snplist")
	p.add_argument("--ref-ld-chr", default=DEFAULT_REF, help="Reference LD-score prefix, including trailing / or dot")
	p.add_argument("--w-ld-chr", default=DEFAULT_WEIGHTS, help="Regression-weight prefix")
	p.add_argument("--N", type=float, help="Explicit fallback only when N is absent")
	p.add_argument(
		"--missing-n",
		choices=("error", "skip"),
		default="error",
		help="skip writes an explicit unavailable status; never invents N",
	)
	p.add_argument("--run-munge", type=boolean, default=True)
	p.add_argument("--run-h2", type=boolean, default=True)
	p.add_argument("--run-rg", type=boolean, default=True)
	p.add_argument("--run", type=boolean, default=True, help="FALSE writes commands only")
	a = p.parse_args()
	files = [Path(x.strip()).resolve() for x in a.gwas_files.split(",")]
	if any(not f.is_file() for f in files):
		p.error("Input file missing")
	if a.N is not None and a.N <= 0:
		p.error("--N must be positive")
	if a.run_munge and not Path(a.merge_alleles).is_file():
		p.error("Missing --merge-alleles file")
	validate_references(a.ref_ld_chr, a.w_ld_chr)
	out = Path(a.output_dir).resolve()
	out.mkdir(parents=True, exist_ok=True)
	python = [a.python] if a.python else [a.conda, "run", "--no-capture-output", "-n", a.conda_env, "python"]
	script = [
		"#!/usr/bin/env bash",
		"set -euo pipefail",
		shlex.join(["rm", "-f", "--"] + [str(out / name) for name in ("h2.tsv", "rg.tsv", "rg.png")]),
	]
	names = []
	statuses = []
	labels = [f.name.removesuffix(".gz") for f in files]
	if len(set(labels)) != len(labels):
		p.error("Trait names must be unique")
	helper = str(Path(__file__).resolve())
	for folder in ("h2.log", "rg.log"):
		(out / folder).mkdir(exist_ok=True)
	for f, trait in zip(files, labels):
		cache, _, _ = cache_paths(f)
		if not a.run_munge:
			validate_sumstats(cache)
			names.append((trait, str(cache)))
			statuses.append(
				(str(f), "SCHEDULED", "Pre-generated adjacent sumstats; standard LDSC/reference covers 1-22 only")
			)
			continue
		with gzip.open(f, "rt") if f.suffix == ".gz" else f.open() as handle:
			header = handle.readline().strip().split()
		if "N" not in header and a.N is None:
			if a.missing_n != "skip":
				p.error(str(f) + ": missing N; provide --N (no assumed sample size)")
			statuses.append(
				(str(f), "UNAVAILABLE", "Missing N: supply a verified sample size; excluded from h2 and rg")
			)
			print("UNAVAILABLE LDSC (missing N): " + str(f), flush=True)
			continue
		if not set(("SNP", "CHR", "EA", "NEA", "BETA", "P")).issubset(header):
			p.error(str(f) + ": need SNP CHR EA NEA BETA P")
		statuses.append(
			(
				str(f),
				"SCHEDULED",
				"Standard LDSC/reference covers 1-22 only; X/Y/MT exclusions recorded in trait qc chromosome audit",
			)
		)
		cache, _, _ = cache_paths(f)
		names.append((trait, str(cache)))
		cmd = [
			sys.executable,
			helper,
			"prepare",
			"--source",
			str(f),
			"--merge-alleles",
			a.merge_alleles,
			"--software",
			a.ldsc_software_dir,
			"--allow-build",
			str(a.run_munge),
			"--interpreter",
		] + python
		if a.N is not None:
			cmd[3:3] = ["--fallback-n", str(a.N)]
		script.append(shlex.join(cmd))
	ref = ["--ref-ld-chr", a.ref_ld_chr, "--w-ld-chr", a.w_ld_chr]
	ldsc = python + [str(Path(a.ldsc_software_dir) / "ldsc.py")]
	if a.run_h2:
		for trait, f in names:
			script.append(shlex.join(ldsc + ["--h2", f, "--out", str(out / "h2.log" / (trait + ".h2"))] + ref))
	if a.run_rg:
		for i in range(len(names) - 1):
			script.append(
				shlex.join(
					ldsc
					+ ["--rg", ",".join(f for _, f in names[i:]), "--out", str(out / "rg.log" / (names[i][0] + ".rg"))]
					+ ref
				)
			)
	script.append(
		shlex.join(
			[
				sys.executable,
				str(Path(__file__).with_name("compare.py")),
				"ldsc-summary",
				"--output-dir",
				str(out),
				"--h2",
				str(a.run_h2),
				"--rg",
				str(a.run_rg),
			]
		)
	)
	cmdfile = out / "ldsc.cmd.sh"
	cmdfile.write_text("\n".join(script) + "\n")
	(out / "inputs.status.tsv").write_text(
		"FILE\tSTATUS\tREASON\n" + "\n".join("\t".join(row) for row in statuses) + "\n"
	)
	if not names:
		p.error("No GWAS with usable N; see inputs.status.tsv")
	print(cmdfile, flush=True)
	if a.run:
		subprocess.run(["bash", str(cmdfile)], check=True)


def fix_p(src, dst, merge_alleles=None):
	import math

	tmp = dst + ".tmp"
	chromosomes = Counter()
	kept = Counter()
	invalid = 0
	allowed = None
	if merge_alleles:
		with open(merge_alleles) as inp:
			allowed = {line.split()[0] for line in inp if line.strip()}
	with gzip.open(src, "rt") if src.endswith(".gz") else open(src) as inp, gzip.open(tmp, "wt") as out:
		header = inp.readline().split()
		pi = header.index("P")
		ci = header.index("CHR")
		si = header.index("SNP")
		out.write("\t".join(header) + "\n")
		for line in inp:
			row = line.split()
			if len(row) != len(header):
				raise ValueError("Malformed GWAS row")
			chr = row[ci].upper().removeprefix("CHR")
			chr = {"X": "23", "Y": "24", "MT": "25", "M": "25"}.get(chr, chr)
			chromosomes[chr] += 1
			# Standard LDSC is autosomal. Audit omitted X/Y explicitly.
			if chr not in AUTOSOMES:
				continue
			if allowed is not None and row[si] not in allowed:
				continue
			try:
				pv = float(row[pi])
			except ValueError:
				invalid += 1
				continue
			if not math.isfinite(pv) or pv < 0 or pv > 1:
				invalid += 1
				continue
			if pv < 1e-300:
				row[pi] = "1e-300"
			out.write("\t".join(row) + "\n")
			kept[chr] += 1
	os.replace(tmp, dst)
	Path(dst + ".chromosomes.json").write_text(
		json.dumps(
			dict(
				input=dict(chromosomes),
				kept=dict(kept),
				invalid_p=invalid,
				scope="autosomes",
				reason="Standard LDSC and installed LD scores cover chromosomes 1-22; X/Y/MT are not estimated",
			),
			indent=2,
		)
	)


def gwas_ldsc_cli():
	import sys

	if len(sys.argv) > 1 and sys.argv[1] == "fix-p":
		fix_p(*sys.argv[2:])
	elif len(sys.argv) > 1 and sys.argv[1] == "prepare":
		parser = argparse.ArgumentParser()
		parser.add_argument("--source", required=True)
		parser.add_argument("--merge-alleles", required=True)
		parser.add_argument("--software", required=True)
		parser.add_argument("--fallback-n", type=float)
		parser.add_argument("--allow-build", choices=("True", "False"), default="True")
		parser.add_argument("--interpreter", nargs=argparse.REMAINDER, required=True)
		args = parser.parse_args(sys.argv[2:])
		prepare(**vars(args))
	else:
		gwas_ldsc_main()


# 🚩 ld_blocks
"""Download the published GRCh37 LDetect blocks and record their provenance."""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil
from urllib.request import urlopen


def read_bed(text, omitted=None):
	rows = []
	for line in text.splitlines():
		fields = line.split()
		if not fields or fields[0].lower() in ("chr", "chrom", "#chrom"):
			continue
		chromosome = fields[0].removeprefix("chr")
		if chromosome not in {str(i) for i in range(1, 23)}:
			raise ValueError("Unexpected LDetect chromosome: " + chromosome)
		if "None" in fields[1:3] and omitted is not None:
			omitted.append(line)
			continue
		start, end = map(int, fields[1:3])
		if not 0 <= start < end:
			raise ValueError("Invalid BED interval: " + line)
		rows.append((int(chromosome), start, end))
	rows.sort()
	if {row[0] for row in rows} != set(range(1, 23)):
		raise ValueError("Expected all 22 autosomes")
	for left, right in zip(rows, rows[1:]):
		if left[0] == right[0] and left[2] > right[1]:
			raise ValueError("Overlapping BED intervals")
	return rows


def ld_blocks_main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--output-dir", default="/mnt/f/refLD/block")
	args = parser.parse_args()
	root = Path(args.output_dir)
	root.mkdir(parents=True, exist_ok=True)
	records = []
	for race in ("AFR", "ASN", "EUR"):
		url = (
			f"https://api.bitbucket.org/2.0/repositories/nygcresearch/ldetect-data/src/master/{race}/fourier_ls-all.bed"
		)
		with urlopen(url, timeout=60) as response:
			raw = response.read()
		omitted = []
		rows = read_bed(raw.decode(), omitted)
		# Upstream AFR chr11 contains two records with an unknown shared edge.
		# Omit these invalid intervals; do not invent a replacement breakpoint.
		if omitted:
			(root / f"{race}.37.omitted.tsv").write_text("\n".join(omitted) + "\n")
		data = "".join(f"chr{chrom}\t{start}\t{end}\n" for chrom, start, end in rows).encode()
		targets = ("EAS", "SAS") if race == "ASN" else (race,)
		# Do not overwrite a different locally curated block set on reruns.
		for target in targets:
			path = root / f"{target}.37.bed"
			if path.exists() and path.read_bytes() != data:
				raise FileExistsError("Existing BED differs: " + str(path))
		stage = root / f"{race}.37.bed"
		if stage.exists() and stage.read_bytes() != data:
			raise FileExistsError("Existing BED differs: " + str(stage))
		stage.write_bytes(data)
		if race == "ASN":
			for target in targets:
				shutil.copyfile(stage, root / f"{target}.37.bed")
			stage.unlink()
		for target in targets:
			records.append(
				[
					f"{target}.37.bed",
					"37",
					target,
					race,
					len(rows),
					"0-based half-open; chromosomes 1-22",
					url,
					hashlib.sha256(raw).hexdigest(),
					hashlib.sha256(data).hexdigest(),
				]
			)
			print(f"{target}.37.bed: {len(rows)} blocks; source={race}", flush=True)
	for race in ("AFR", "EAS", "EUR", "SAS"):
		path = root / f"{race}.38.bed"
		if path.exists():
			rows = read_bed(path.read_text())
			records.append(
				[
					path.name,
					"38",
					race,
					race,
					len(rows),
					"0-based half-open; chromosomes 1-22",
					"https://github.com/jmacdon/LDblocks_GRCh38",
					"",
					hashlib.sha256(path.read_bytes()).hexdigest(),
				]
			)
	with (root / "block_sources.tsv").open("w", newline="") as handle:
		writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
		writer.writerow(
			[
				"file",
				"grch",
				"race",
				"source_population",
				"blocks",
				"coordinates",
				"source_url",
				"download_sha256",
				"file_sha256",
			]
		)
		writer.writerows(records)
	(root / "README.bplot.md").write_text(
		"# LD block references for bplot\n\n"
		"GRCh37: published LDetect fourier_ls-all.bed files, with BED coordinates retained.\n"
		"AFR: two upstream chr11 intervals contain a None endpoint and are omitted, "
		"recorded in AFR.37.omitted.tsv. No replacement boundary was inferred.\n"
		"EAS.37.bed and SAS.37.bed are identical copies of ASN, as requested; "
		"they are not independently inferred EAS/SAS block sets.\n\n"
		"GRCh38: existing pyrho files renamed to [race].38.bed without changing contents.\n"
		"Both sets contain autosomes 1-22 only. No HIS, ALL or chrX boundary set is supplied.\n"
		"See block_sources.tsv for source URLs, population mappings and checksums.\n"
	)


def ld_blocks_cli():
	ld_blocks_main()


# 🚩 gwas_compare_project
"""Resolve project GWAS files and create non-destructive build-aligned comparison inputs."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def build_of(file):
	trait = file.name.removesuffix(".gz").removesuffix(".thin")
	candidates = [Path(str(file) + ".grch"), file.parent.parent / "qc" / (trait + ".grch")]
	builds = {p.read_text().strip() for p in candidates if p.is_file()}
	if len(builds) != 1 or not builds.issubset({"37", "38"}):
		raise ValueError("Missing/conflicting genome-build metadata for " + str(file))
	return next(iter(builds))


def resolve(project, category, anchor):
	base = Path(project) / category
	files = [p / "gwas" / (p.name + ".gz") for p in base.iterdir() if p.is_dir()]
	files = sorted(
		(p for p in files if p.is_file()),
		key=lambda p: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", p.parent.parent.name)],
	)
	if anchor:
		first = [p for p in files if p.parent.parent.name == anchor]
		if not first:
			raise ValueError(f"Anchor {anchor} not found; finish processing that GWAS first")
		files = first + [p for p in files if p not in first]
	return files


def require_build(files, required):
	builds = {str(f): build_of(f) for f in files}
	wrong = [f"{Path(f).stem}=GRCh{b}" for f, b in builds.items() if b != required]
	if wrong:
		raise ValueError(f"All inputs must be GRCh{required}; finish format.sh liftover first: " + ", ".join(wrong))
	return builds


def gwas_compare_project_main():
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("mode", choices=("compare", "ldsc", "shiny"))
	p.add_argument("--dir-gwas", "--project-dir", dest="project_dir")
	p.add_argument("--category", default="common")
	p.add_argument("--anchor")
	p.add_argument("--gwas-files")
	p.add_argument("--grch", choices=("37", "38"))
	p.add_argument(
		"--require-grch",
		choices=("37", "38"),
		help="Require every finalized input to have this build; do not lift comparison copies",
	)
	p.add_argument("--dir-out", "--output-dir", dest="output_dir", help="Default: /mnt/d/analysis/gwas/<mode>")
	p.add_argument("--check-only", choices=("TRUE", "FALSE"), default="FALSE")
	a, other = p.parse_known_args()
	if bool(a.project_dir) == bool(a.gwas_files):
		p.error("Specify one of --dir-gwas or --gwas-files")
	files = (
		resolve(a.project_dir, a.category, a.anchor)
		if a.project_dir
		else [Path(x).resolve() for x in a.gwas_files.split(",")]
	)
	minimum = 1 if a.mode == "shiny" else 2
	if len(files) < minimum or any(not f.is_file() for f in files):
		p.error(f"{a.mode} requires at least {minimum} existing GWAS file(s)")
	if a.mode == "shiny":
		files = [
			f.with_name(f.name[:-3] + ".thin.gz")
			if not f.name.endswith(".thin.gz") and f.with_name(f.name[:-3] + ".thin.gz").is_file()
			else f
			for f in files
		]
	labels = [
		f.name.removesuffix(".gz").removesuffix(".thin") if a.mode == "shiny" else f.name.removesuffix(".gz")
		for f in files
	]
	checked_builds = {}
	if a.require_grch:
		if a.grch and a.grch != a.require_grch:
			p.error("--grch and --require-grch must agree")
		checked_builds = require_build(files, a.require_grch)
		a.grch = a.require_grch
	out = Path(a.output_dir or ("/mnt/d/analysis/gwas/" + a.mode)).resolve()
	if a.mode == "shiny":
		a.grch = a.grch or "37"
		builds = {str(f): build_of(f) for f in files}
		if len(set(builds.values())) != 1:
			p.error("shiny inputs must share one source GRCh build; harmonize the inputs first")
		if len(set(labels)) != len(labels):
			p.error("shiny requires unique GWAS filenames")
		source_build = next(iter(builds.values()))
		out.mkdir(parents=True, exist_ok=True)
		config = out / "shiny.inputs.json"
		config.write_text(
			json.dumps(
				dict(
					source_build=source_build,
					tracks=[dict(file=str(f.resolve()), trait=label) for f, label in zip(files, labels)],
				),
				indent=2,
			)
		)
		print(f"shiny: {len(files)} GWAS; source GRCh{source_build}; initial GRCh{a.grch or source_build}", flush=True)
		if a.check_only == "TRUE":
			print(config)
			return
		os.execvp(
			"Rscript",
			[
				"Rscript",
				"--vanilla",
				str(Path(__file__).resolve().parent.parent / "shiny" / "app.R"),
				"--manifest",
				str(config),
				"--output-dir",
				str(out),
				"--grch",
				a.grch or source_build,
			]
			+ other,
		)
		return
	manifest = []
	aligned = []
	for f in files:
		build = checked_builds.get(str(f)) or (build_of(f) if a.mode == "compare" else "not_required_for_rsID_LDSC")
		target = f
		if a.mode == "compare" and a.grch and build != a.grch:
			target = out / "inputs" / ("grch" + a.grch) / f.name
			if a.check_only != "TRUE":
				chain = {"37": "hg19ToHg38.over.chain.gz", "38": "hg38ToHg19.over.chain.gz"}[build]
				subprocess.run(
					[
						sys.executable,
						str(Path(__file__).with_name("format.py")),
						"liftover",
						"--input",
						str(f),
						"--output",
						str(target),
						"--qc-prefix",
						str(target.parent / "qc" / f.stem),
						"--source-build",
						build,
						"--target-build",
						a.grch,
						"--chain",
						"/mnt/d/files/liftOver/" + chain,
						"--liftOver",
						"/mnt/d/software/bin/liftOver",
					],
					check=True,
				)
		aligned.append(target)
		manifest.append(dict(trait=f.stem, input=str(f), source_build=build, comparison_input=str(target)))
	if a.mode == "compare" and not a.grch and len({x["source_build"] for x in manifest}) > 1:
		p.error("Mixed builds; supply --grch 38 for cached harmonization")
	out.mkdir(parents=True, exist_ok=True)
	import csv

	with (out / "inputs.tsv").open("w") as handle:
		writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter="\t", lineterminator="\n")
		writer.writeheader()
		writer.writerows(manifest)
	print(f"Discovered {len(files)} GWAS: " + ",".join(labels), flush=True)
	if a.check_only == "TRUE":
		print(out / "inputs.tsv")
		return
	args = ["--gwas-files", ",".join(map(str, aligned)), "--output-dir", str(out)] + other
	if a.mode == "compare":
		if "--labels" not in other:
			args += ["--labels", ",".join(labels)]
		if a.grch:
			args += ["--grch", a.grch]
		cmd = ["Rscript", str(Path(__file__).with_name("compare.R")), "compare"] + args
	else:
		cmd = [sys.executable, str(Path(__file__).with_name("compare.py")), "ldsc-run"] + args
	subprocess.run(cmd, check=True)


def gwas_compare_project_cli():
	gwas_compare_project_main()


SUBCOMMANDS = {
	"ldsc-summary": gwas_ldsc_summary_cli,
	"ldsc-run": gwas_ldsc_cli,
	"ld-blocks": ld_blocks_cli,
}


def main():
	import sys

	if len(sys.argv) > 1 and sys.argv[1] in SUBCOMMANDS:
		return SUBCOMMANDS[sys.argv.pop(1)]()
	if sys.argv[1:2] in (["--help"], ["-h"]):
		print("Helper commands: " + ", ".join(SUBCOMMANDS))
	return gwas_compare_project_cli()


if __name__ == "__main__":
	main()
