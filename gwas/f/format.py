#!/usr/bin/env python3
"""format workflow utilities; use --help for commands."""

from __future__ import annotations


# 🚩 gwas_alleles
"""Indexed FASTA access and sequence allele normalization (1-based positions)."""
from bisect import bisect_right
from functools import lru_cache
import gzip
from pathlib import Path
import struct


class IndexedFasta:
	def __init__(self, path):
		self.index = {}
		for line in Path(str(path) + ".fai").read_text().splitlines():
			name, length, offset, bases, width = line.split()[:5]
			self.index[name] = tuple(map(int, (length, offset, bases, width)))
		self._compressed = str(path).endswith((".gz", ".bgz", ".bgzf"))
		if self._compressed:
			# GZI entries map uncompressed byte offsets to BGZF block offsets;
			# the first (0, 0) block is implicit in the on-disk index.
			data = Path(str(path) + ".gzi").read_bytes()
			if len(data) < 8 or len(data) != 8 + 16 * struct.unpack_from("<Q", data)[0]:
				raise ValueError(f"Invalid BGZF index: {path}.gzi")
			pairs = [(0, 0), *struct.iter_unpack("<QQ", data[8:])]
			self._compressed_offsets, self._plain_offsets = map(list, zip(*pairs))
			self._compressed_size = Path(path).stat().st_size
			if (
				any(b[0] <= a[0] or b[1] <= a[1] for a, b in zip(pairs, pairs[1:]))
				or pairs[-1][0] >= self._compressed_size
			):
				raise ValueError(f"Invalid BGZF block offsets: {path}.gzi")
		self.handle = open(path, "rb")
		if self._compressed:
			header = self.handle.read(18)
			if header[:4] != b"\x1f\x8b\x08\x04" or header[12:16] != b"BC\x02\x00":
				self.handle.close()
				raise ValueError(f"Expected BGZF-compressed FASTA (use bgzip): {path}")
			self._block = lru_cache(maxsize=32)(self._read_block)

	def _read_block(self, block):
		start = self._compressed_offsets[block]
		end = (
			self._compressed_offsets[block + 1] if block + 1 < len(self._compressed_offsets) else self._compressed_size
		)
		# The final indexed block may be followed by the 28-byte BGZF EOF marker.
		if not 0 < end - start <= 65536 + 28:
			raise ValueError("Invalid BGZF compressed block size")
		self.handle.seek(start)
		compressed = self.handle.read(end - start)
		if len(compressed) != end - start:
			raise ValueError("Truncated BGZF FASTA")
		data = gzip.decompress(compressed)
		if len(data) > 65536 or (
			block + 1 < len(self._plain_offsets)
			and len(data) != self._plain_offsets[block + 1] - self._plain_offsets[block]
		):
			raise ValueError("BGZF index does not match FASTA blocks")
		return data

	def _read_bytes(self, offset, length):
		if not self._compressed:
			self.handle.seek(offset)
			return self.handle.read(length)
		block = bisect_right(self._plain_offsets, offset) - 1
		parts = []
		while length:
			data = self._block(block)
			within = offset - self._plain_offsets[block]
			if not 0 <= within < len(data):
				raise ValueError("FASTA interval exceeds indexed BGZF data")
			part = data[within : within + length]
			parts.append(part)
			offset += len(part)
			length -= len(part)
			block += 1
			if length and block >= len(self._plain_offsets):
				raise ValueError("Truncated indexed BGZF FASTA")
		return b"".join(parts)

	def fetch(self, chrom, pos, length):
		name = {"23": "X", "24": "Y", "25": "MT"}.get(str(chrom), str(chrom))
		if name not in self.index:
			name = "chr" + ("M" if name == "MT" else name)
		size, offset, bases, width = self.index[name]
		if pos < 1 or length < 1 or pos + length - 1 > size:
			raise ValueError(f"FASTA interval out of bounds: {chrom}:{pos}+{length}")
		start = pos - 1
		byte_start = offset + start // bases * width + start % bases
		end = start + length - 1
		byte_end = offset + end // bases * width + end % bases
		sequence = self._read_bytes(byte_start, byte_end - byte_start + 1)
		sequence = sequence.replace(b"\n", b"").replace(b"\r", b"").decode().upper()
		if len(sequence) != length:
			raise ValueError("FASTA index does not match sequence length")
		return sequence

	def close(self):
		if self._compressed:
			self._block.cache_clear()
		self.handle.close()


def normalize(pos, ref, alt, fasta=None, chrom=None):
	"""Trim identical sequence and, with FASTA, left-align repeat indels."""
	if ref == alt or not ref or not alt or not set(ref + alt) <= set("ACGT"):
		raise ValueError("Expected two distinct DNA alleles")
	while ref[-1] == alt[-1]:
		if min(len(ref), len(alt)) == 1:
			if fasta is None or pos <= 1 or len(ref) == len(alt):
				break
			previous = fasta.fetch(chrom, pos - 1, 1)
			if previous not in "ACGT":
				break
			ref, alt, pos = previous + ref, previous + alt, pos - 1
		ref, alt = ref[:-1], alt[:-1]
	while min(len(ref), len(alt)) > 1 and ref[0] == alt[0]:
		ref, alt, pos = ref[1:], alt[1:], pos + 1
	return pos, ref, alt


def id_alleles(snp, pos):
	fields = snp.split(":")
	if len(fields) != 4 or fields[1] != str(pos):
		return None
	ref, alt = fields[2:]
	if not ref or not alt or ref == alt or not set(ref + alt) <= set("ACGT"):
		return None
	return ref, alt


# 🚩 gwas_liftover
"""Stream a complete standardized GWAS through UCSC liftOver, retaining X/Y."""
import argparse
from collections import Counter
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


PRIMARY_CHROMOSOMES = frozenset(str(c) for c in range(1, 26))

# GRC PAR intervals, inclusive and 1-based. Keep source-X PAR associations on
# X when a chain picks the homologous Y representation. True source-Y and
# non-PAR mappings are not changed.
PAR = {
	37: ((60001, 2699520, 10001, 2649520), (154931044, 155260560, 59034050, 59363566)),
	38: ((10001, 2781479, 10001, 2781479), (155701383, 156030895, 56887903, 57217415)),
}


def source_x_par_target(source_chr, source_pos, source_build, target_build, mapping):
	if chrom(source_chr) != "23" or chrom(mapping[0]) != "24":
		return mapping
	start, end = int(mapping[1]) + 1, int(mapping[2])
	for old, new in zip(PAR[source_build], PAR[target_build]):
		if old[0] <= int(source_pos) <= old[1] and new[2] <= start <= end <= new[3]:
			result = list(mapping)
			offset = new[0] - new[2]
			result[0], result[1], result[2] = "chrX", str(start - 1 + offset), str(end + offset)
			return result
	return mapping


def identity(path):
	p = Path(path).resolve()
	s = p.stat()
	return [str(p), s.st_size, s.st_mtime_ns]


def chrom(value):
	value = value.upper().removeprefix("CHR")
	return {"X": "23", "Y": "24", "M": "25", "MT": "25"}.get(value, value)


def ucsc(value):
	return "chr" + {"23": "X", "24": "Y", "25": "M"}.get(chrom(value), chrom(value))


def reverse_complement(value):
	if not value or any(c not in "ACGT" for c in value.upper()):
		raise ValueError("non-DNA allele on reverse strand")
	return value.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def gwas_liftover_run(a):
	src, dst = Path(a.input), Path(a.output)
	if src.resolve() == dst.resolve():
		raise ValueError("Source and target must be distinct; retain the source build")
	qc = Path(a.qc_prefix)
	qc.parent.mkdir(parents=True, exist_ok=True)
	dst.parent.mkdir(parents=True, exist_ok=True)
	state = Path(str(qc) + ".liftover.state.json")
	# Reference files are optional for legacy SNV-only use, but mandatory for
	# preserving sequence-resolved indels. Never interpret a point lift as an
	# indel lift without checking its full reference interval.
	fasta_paths = []
	for option, build in (("source_fasta", a.source_build), ("target_fasta", a.target_build)):
		path = getattr(a, option, None)
		if not path:
			candidate = Path(f"/mnt/f/gen/1hgp/{build}/GRCH{build}.fasta.gz")
			path = str(candidate) if all(Path(str(candidate) + ext).is_file() for ext in ("", ".fai", ".gzi")) else None
		fasta_paths.append(path)
	signature = dict(
		version=5,
		source=identity(src),
		chain=identity(a.chain),
		source_build=a.source_build,
		target_build=a.target_build,
		fastas=[identity(p) if p else None for p in fasta_paths],
		fasta_indexes=[
			[
				identity(str(p) + ext)
				for ext in ((".fai", ".gzi") if str(p).endswith((".gz", ".bgz", ".bgzf")) else (".fai",))
			]
			if p
			else None
			for p in fasta_paths
		],
	)
	if a.replace == "FALSE" and dst.exists() and state.exists():
		saved = json.loads(state.read_text())
		if saved.get("signature") == signature and saved.get("output") == identity(dst):
			print("Reuse verified liftover: " + str(dst), flush=True)
			return
	# Invalidate before work; a failed run must not leave a success state.
	state.unlink(missing_ok=True)
	env = dict(os.environ, LC_ALL="C")
	counts = Counter()
	by_chr = Counter()
	source_fasta, target_fasta = [IndexedFasta(p) if p else None for p in fasta_paths]
	sequence_rows = {}
	with tempfile.TemporaryDirectory(prefix="liftover.", dir=qc.parent) as temp:
		work = Path(temp)
		bed, old = work / "input.bed", work / "input.tsv"
		opener = gzip.open if str(src).endswith((".gz", ".bgz")) else open
		with opener(src, "rt") as inp, bed.open("w") as bp, old.open("w") as op:
			header = inp.readline().rstrip("\r\n").split("\t")
			ci, pi, ai, bi = [header.index(c) for c in ("CHR", "POS", "EA", "NEA")]
			for n, line in enumerate(inp, 1):
				row = line.rstrip("\r\n").split("\t")
				if len(row) != len(header):
					raise ValueError(f"Malformed source row {n}")
				pos = int(row[pi])
				if pos < 1:
					raise ValueError(f"Invalid position at row {n}")
				if chrom(row[ci]) == "25" and pos > 16569:
					raise ValueError(
						"CHR=25 exceeds mitochondrial length. If the source uses PLINK "
						"25=XY/PAR, standardize it to X (23) before liftover."
					)
				key = f"{n:012d}"
				length = 1
				if max(len(row[ai]), len(row[bi])) > 1 and source_fasta and target_fasta:
					pair = id_alleles(row[header.index("SNP")], pos)
					if pair and set(pair) == {row[ai], row[bi]}:
						ref, alt = pair
						if source_fasta.fetch(row[ci], pos, len(ref)) == ref:
							sequence_rows[key] = (ref, alt, row[ai] == ref)
							length = len(ref)
						else:
							sequence_rows[key] = "source_reference_mismatch"
				bp.write(f"{ucsc(row[ci])}\t{pos - 1}\t{pos + length - 1}\t{key}\t0\t+\n")
				op.write(key + "\t" + "\t".join(row) + "\n")
				counts["input"] += 1
				by_chr[(chrom(row[ci]), "input")] += 1
		lifted, unmapped = work / "lifted.bed", work / "unmapped.bed"
		with Path(str(qc) + ".liftover.log").open("w") as log:
			subprocess.run(
				[a.liftOver, str(bed), a.chain, str(lifted), str(unmapped)],
				stdout=log,
				stderr=subprocess.STDOUT,
				check=True,
			)
		ordered = work / "mapped.sorted.bed"
		with ordered.open("w") as out:
			subprocess.run(["sort", "-T", temp, "-S", "256M", "-k4,4", str(lifted)], stdout=out, env=env, check=True)
		mapped = work / "mapped.tsv"
		with (
			ordered.open() as lp,
			old.open() as op,
			mapped.open("w") as out,
			gzip.open(str(qc) + ".liftover.par.tsv.gz", "wt") as par_audit,
			gzip.open(str(qc) + ".liftover.nonprimary.tsv.gz", "wt") as nonprimary,
			gzip.open(str(qc) + ".liftover.alleles.tsv.gz", "wt") as allele_audit,
		):
			allele_audit.write("SNP\tCHR_SOURCE\tPOS_SOURCE\tEA_SOURCE\tNEA_SOURCE\tREASON\n")
			nonprimary.write(
				"\t".join(header + ["TARGET_CONTIG", "TARGET_START_0", "TARGET_END_0", "TARGET_STRAND", "REASON"])
				+ "\n"
			)
			par_audit.write("SNP\tCHR_SOURCE\tPOS_SOURCE\tCHAIN_TARGET\tCHAIN_POS\tFINAL_TARGET\tFINAL_POS\tREASON\n")
			current = next(lp, "").split()
			for line in op:
				key, rest = line.rstrip("\n").split("\t", 1)
				row = rest.split("\t")
				source_chr = chrom(row[ci])
				if not current or current[3] != key:
					counts["unmapped"] += 1
					by_chr[(source_chr, "unmapped")] += 1
					continue
				if len(current) != 6:
					raise ValueError("Invalid mapped BED6 interval")
				canonical = source_x_par_target(source_chr, row[pi], a.source_build, a.target_build, current)
				if canonical != current:
					par_audit.write(
						"\t".join(
							[
								row[header.index("SNP")],
								source_chr,
								row[pi],
								current[0],
								str(int(current[1]) + 1),
								canonical[0],
								str(int(canonical[1]) + 1),
								"source_X_PAR_retained_on_homologous_target_X_PAR",
							]
						)
						+ "\n"
					)
					current = canonical
					counts["source_x_par_retained_on_x"] += 1
				if chrom(current[0]) not in PRIMARY_CHROMOSOMES:
					# The project schema/LD references support primary chromosomes.
					# Numeric sorting would interleave chr1 with chr1_* contigs.
					# Preserve the entire source association and mapping separately;
					# do not mislabel these mappings as an autosomal position.
					nonprimary.write(
						"\t".join(
							row
							+ [
								current[0],
								current[1],
								current[2],
								current[5],
								"nonprimary_target_unsupported_by_project_schema_and_LD_references",
							]
						)
						+ "\n"
					)
					counts["nonprimary_excluded"] += 1
					by_chr[(source_chr, "nonprimary_excluded")] += 1
					previous = current[3]
					current = next(lp, "").split()
					if current and current[3] == previous:
						raise ValueError("Multiple mappings for the same source variant")
					continue
				seq = sequence_rows.get(key)
				reject = None
				new_pos = int(current[1]) + 1
				expected_length = len(seq[0]) if isinstance(seq, tuple) else 1
				if int(current[2]) - int(current[1]) != expected_length:
					reject = "reference_interval_changed_length"
				if isinstance(seq, str):
					reject = seq
				if isinstance(seq, tuple) and reject is None:
					ref, alt, effect_is_ref = seq
					if current[5] == "-":
						ref, alt = reverse_complement(ref), reverse_complement(alt)
					target_chr = chrom(current[0])
					if target_fasta.fetch(target_chr, new_pos, len(ref)) != ref:
						reject = "target_reference_mismatch"
					else:
						new_pos, ref, alt = normalize(new_pos, ref, alt, target_fasta, target_chr)
						row[ai], row[bi] = (ref, alt) if effect_is_ref else (alt, ref)
						counts["sequence_alleles_lifted"] += 1
				if reject:
					allele_audit.write(
						"\t".join(
							[
								row[header.index("SNP")],
								row[ci],
								row[pi],
								row[ai],
								row[bi],
								reject + "; association_excluded",
							]
						)
						+ "\n"
					)
					counts["alleles_excluded"] += 1
					by_chr[(source_chr, "alleles_excluded")] += 1
					previous = current[3]
					current = next(lp, "").split()
					if current and current[3] == previous:
						raise ValueError("Multiple mappings for the same source variant")
					continue
				# Unknown allele codes or indels lacking explicit REF/ALT and
				# FASTA validation retain association statistics with NA alleles.
				unresolved = any(len(row[i]) != 1 or row[i].upper() not in "ACGT" for i in (ai, bi))
				if isinstance(seq, tuple):
					unresolved = False
				if unresolved:
					allele_audit.write(
						"\t".join(
							[
								row[header.index("SNP")],
								row[ci],
								row[pi],
								row[ai],
								row[bi],
								"alleles_not_resolved_by_point_liftover; association_retained",
							]
						)
						+ "\n"
					)
					row[ai] = row[bi] = "NA"
					counts["alleles_unresolved"] += 1
					by_chr[(source_chr, "alleles_unresolved")] += 1
				try:
					if current[5] == "-" and not unresolved and not isinstance(seq, tuple):
						row[ai], row[bi] = reverse_complement(row[ai]), reverse_complement(row[bi])
						counts["reverse_complemented"] += 1
					row[ci], row[pi] = chrom(current[0]), str(new_pos)
					out.write("\t".join(row) + "\n")
					counts["lifted"] += 1
					by_chr[(source_chr, "lifted")] += 1
				except ValueError:
					counts["unsupported_reverse_allele"] += 1
					by_chr[(source_chr, "unsupported_reverse_allele")] += 1
				previous = current[3]
				current = next(lp, "").split()
				if current and current[3] == previous:
					raise ValueError("Multiple mappings for the same source variant")
			if current:
				raise ValueError("Unconsumed mapping rows")
		if not counts["lifted"]:
			raise ValueError("No variants lifted; inspect chain/build and liftover.log")
		tmp = Path(str(dst) + ".tmp")
		with tmp.open("wb") as out:
			bg = subprocess.Popen(["bgzip", "-@", "4", "-c"], stdin=subprocess.PIPE, stdout=out)
			try:
				bg.stdin.write(("\t".join(header) + "\n").encode())
				with subprocess.Popen(
					[
						"sort",
						"-T",
						temp,
						"-S",
						"512M",
						"-t",
						"\t",
						f"-k{ci + 1},{ci + 1}n",
						f"-k{pi + 1},{pi + 1}n",
						"-k1,1",
						str(mapped),
					],
					stdout=subprocess.PIPE,
					env=env,
				) as sorter:
					shutil.copyfileobj(sorter.stdout, bg.stdin, 1024 * 1024)
					if sorter.wait():
						raise RuntimeError("Coordinate sort failed")
				bg.stdin.close()
				if bg.wait():
					raise RuntimeError("BGZF compression failed")
			except BaseException:
				bg.kill()
				bg.wait()
				raise
		subprocess.run(
			["tabix", "-s", str(ci + 1), "-b", str(pi + 1), "-e", str(pi + 1), "-S", "1", str(tmp)], check=True
		)
		os.replace(tmp, dst)
		os.replace(str(tmp) + ".tbi", str(dst) + ".tbi")
		Path(str(dst) + ".csi").unlink(missing_ok=True)
		with gzip.open(str(qc) + ".liftover.unmapped.bed.gz", "wb") as out, unmapped.open("rb") as inp:
			shutil.copyfileobj(inp, out)
	with Path(str(qc) + ".liftover.n.tsv").open("w") as out:
		out.write(
			"CHR\tN_INPUT\tN_LIFTED\tN_UNMAPPED\tN_ALLELES_UNRESOLVED\tN_ALLELES_EXCLUDED\tN_NONPRIMARY_EXCLUDED\n"
		)
		for c in sorted({c for c, _ in by_chr}, key=lambda c: int(c)):
			out.write(
				c
				+ "\t"
				+ "\t".join(
					str(by_chr[c, k])
					for k in (
						"input",
						"lifted",
						"unmapped",
						"alleles_unresolved",
						"alleles_excluded",
						"nonprimary_excluded",
					)
				)
				+ "\n"
			)
	for fasta in (source_fasta, target_fasta):
		if fasta:
			fasta.close()
	Path(str(src) + ".grch").write_text(str(a.source_build) + "\n")
	Path(str(dst) + ".grch").write_text(str(a.target_build) + "\n")
	Path(str(qc) + ".grch").write_text(str(a.target_build) + "\n")
	state.write_text(json.dumps(dict(signature=signature, output=identity(dst), counts=counts), indent=2))
	print("Liftover completed: " + json.dumps(counts), flush=True)


def gwas_liftover_cli():
	p = argparse.ArgumentParser(description=__doc__)
	for name in ("input", "output", "chain", "qc-prefix"):
		p.add_argument("--" + name, required=True)
	p.add_argument("--liftOver", default="liftOver")
	p.add_argument("--source-fasta", help="Indexed source FASTA for sequence-resolved indels")
	p.add_argument("--target-fasta", help="Indexed target FASTA for validation and left alignment")
	p.add_argument("--source-build", type=int, choices=(37, 38), required=True)
	p.add_argument("--target-build", type=int, choices=(37, 38), required=True)
	p.add_argument("--replace", choices=("TRUE", "FALSE"), default="FALSE")
	gwas_liftover_run(p.parse_args())


# 🚩 gwas_h2
"""Summary-statistic SNP heritability, with explicit estimands and unavailable components."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from compare import DEFAULT_REF, DEFAULT_WEIGHTS, DEFAULT_ALLELES, validate_references, cache_paths


def parse_h2(text):
	pattern = r"Total Observed scale h2:\s*([-+\deE.]+)\s*\(([-+\deE.]+)\)"
	match = re.search(pattern, text)
	if not match:
		raise ValueError("No observed-scale h2 estimate found in LDSC log")
	import math

	estimate, se = map(float, match.groups())
	if not all(map(math.isfinite, (estimate, se))) or se < 0:
		raise ValueError("LDSC returned non-finite h2/SE")
	return estimate, se


def report_rows(trait, sex, estimate, se, chromosomes):
	auto = set(str(c) for c in range(1, 23))
	extra = set(chromosomes) - auto
	rows = []

	def add(metric, value="NA", error="NA", scope="all", status="UNAVAILABLE", reason=""):
		rows.append(
			[
				trait,
				metric,
				value,
				error,
				"LDSC" if value != "NA" else "not_estimated",
				"observed",
				scope,
				sex,
				status,
				reason,
			]
		)

	add(
		"h2_autosome",
		estimate,
		se,
		"1-22",
		"ESTIMATED",
		"SNP heritability tagged by the specified LD reference; not pedigree heritability",
	)
	if extra:
		add(
			"h2_total",
			reason="Cannot call autosomal LDSC whole-genome h2; unestimated chromosomes: " + ",".join(sorted(extra)),
		)
	else:
		add(
			"h2_total", estimate, se, "1-22", "ESTIMATED", "Input contains autosomes only; same estimate as h2_autosome"
		)
	for name, code in [("X", "23"), ("Y", "24")]:
		add(
			"h2_chr" + name,
			scope=name,
			status="UNAVAILABLE" if code in chromosomes else "NOT_PRESENT",
			reason="Standard LDSC/reference does not support this chromosome; requires a validated sex/ploidy-specific LD method"
			if code in chromosomes
			else "No input variants on this chromosome",
		)
	add(
		"h2_sig",
		reason="Do not run LDSC on GWAS-selected significant SNPs: selection bias. Requires an appropriate LD-aware variance/partition model and phenotype-scale information or individual data",
	)
	add(
		"h2_sig_chrX",
		scope="X",
		reason="Same limitation as h2_sig, plus male-X dosage/LD scaling; not the X proportion of significant-SNP heritability",
	)
	add(
		"h2_related",
		reason="Pedigree/close-relative variance requires individual phenotypes and family/genotype data; a reference kinship matrix alone is insufficient",
	)
	for target in ("male", "female"):
		if sex == target:
			add(
				"h2_" + target + "s",
				estimate,
				se,
				"1-22",
				"ESTIMATED",
				"Autosomal estimate for the declared sex-specific GWAS; not an independent estimate or all-chromosome GREML",
			)
		else:
			add(
				"h2_" + target + "s",
				reason="Requires a corresponding sex-specific GWAS or individual-level analysis; cannot split pooled summary statistics",
			)
	return rows


def gwas_h2_run(a):
	out = Path(a.output_dir).resolve()
	out.mkdir(parents=True, exist_ok=True)
	done = out / "h2.done.json"
	table = out / (a.trait + ".h2.tsv")
	refs = validate_references(a.ref_ld_chr, a.w_ld_chr)
	signature = dict(
		version=2,
		source=identity(a.gwas_file),
		sex=a.sex,
		references=[identity(f) for f in refs],
		alleles=identity(a.merge_alleles),
		python=a.python,
		conda_env=a.conda_env,
	)
	if a.replace == "FALSE" and done.exists() and table.exists():
		previous = json.loads(done.read_text())
		if previous.get("signature") == signature and previous.get("table") == identity(table):
			print("Reuse h2: " + str(table))
			return
	done.unlink(missing_ok=True)
	cmd = [
		sys.executable,
		str(Path(__file__).with_name("compare.py")),
		"ldsc-run",
		"--gwas-files",
		a.gwas_file,
		"--output-dir",
		str(out / "ldsc"),
		"--merge-alleles",
		a.merge_alleles,
		"--ref-ld-chr",
		a.ref_ld_chr,
		"--w-ld-chr",
		a.w_ld_chr,
		"--run-rg",
		"FALSE",
		"--conda-env",
		a.conda_env,
	]
	if a.python:
		cmd += ["--python", a.python]
	with (out / "h2.run.log").open("w") as log:
		subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
	stem = Path(a.gwas_file).name.removesuffix(".gz")
	estimate, se = parse_h2((out / "ldsc" / "h2.log" / (stem + ".h2.log")).read_text())
	with cache_paths(Path(a.gwas_file))[2].open() as handle:
		chromosomes = {row["CHR"]: int(row["INPUT"]) for row in csv.DictReader(handle, delimiter="\t")}
	rows = report_rows(a.trait, a.sex, estimate, se, chromosomes)
	tmp = table.with_suffix(".tmp")
	with tmp.open("w") as handle:
		writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
		writer.writerow(
			["GWAS", "METRIC", "ESTIMATE", "SE", "METHOD", "SCALE", "CHROMOSOMES", "SEX", "STATUS", "REASON"]
		)
		writer.writerows(rows)
	os.replace(tmp, table)
	(out / "h2.components.log").write_text("\n".join("\t".join(map(str, row)) for row in rows) + "\n")
	(out / "h2.err").unlink(missing_ok=True)
	done.write_text(
		json.dumps(
			dict(
				signature=signature,
				table=identity(table),
				status="completed_supported_estimates; see component statuses",
			),
			indent=2,
		)
	)
	print(str(table))


def gwas_h2_cli():
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("--gwas-file", required=True)
	p.add_argument("--trait", required=True)
	p.add_argument("--output-dir", required=True)
	p.add_argument("--sex", choices=("unknown", "mixed", "male", "female"), default="unknown")
	p.add_argument("--ref-ld-chr", default=DEFAULT_REF)
	p.add_argument("--w-ld-chr", default=DEFAULT_WEIGHTS)
	p.add_argument("--merge-alleles", default=DEFAULT_ALLELES)
	p.add_argument("--conda-env", default="ldsc")
	p.add_argument("--python")
	p.add_argument("--replace", choices=("TRUE", "FALSE"), default="FALSE")
	args = p.parse_args()
	try:
		gwas_h2_run(args)
	except Exception as exc:
		out = Path(args.output_dir)
		out.mkdir(parents=True, exist_ok=True)
		(out / "h2.done.json").unlink(missing_ok=True)
		(out / "h2.err").write_text(str(exc) + "\nSee h2.run.log for tool output.\n")
		raise


# 🚩 gwas_magma_ids
"""Map coordinate GWAS IDs to MAGMA reference rsIDs in the declared build.

The reusable SQLite cache contains only rsIDs present in the MAGMA BIM.
Coordinate matching requires both alleles; ambiguous matches are discarded.
"""
import argparse
from contextlib import closing
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import sys
import shutil
import tempfile
import time


def chromosome(value):
	value = value.removeprefix("chr").upper()
	return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(value) or int(value)


def signature(*paths):
	return hashlib.sha256(
		json.dumps(
			[(str(Path(p).resolve()), Path(p).stat().st_size, Path(p).stat().st_mtime_ns) for p in paths]
		).encode()
	).hexdigest()[:24]


def reference_cache(bim, dbsnp, cache_dir):
	target = Path(cache_dir) / (signature(bim, dbsnp) + ".sqlite3")
	target.parent.mkdir(parents=True, exist_ok=True)
	with open(str(target) + ".lock", "w") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		if target.exists():
			return target
		print("Build MAGMA rsID cache (one dbSNP scan): " + str(target), file=sys.stderr, flush=True)
		ids = set()
		with open(bim, buffering=8 * 1024 * 1024) as source:
			for line in source:
				snp = line.split()[1]
				if snp.startswith("rs") and snp[2:].isdigit():
					ids.add(int(snp[2:]))
		if not ids:
			raise ValueError("Coordinate mapping requires rsIDs in the MAGMA BIM")
		temp = Path(str(target) + f".tmp.{os.getpid()}")
		workspace = tempfile.TemporaryDirectory(prefix="magma-rsid-build-", dir="/tmp")
		local_db = Path(workspace.name) / "reference.sqlite3"
		try:
			# A Connection context manager commits but does not close the file.
			# Close before publishing: renaming an open database on WSL/DrvFS
			# can leave the destination unavailable to subsequent readers.
			with closing(sqlite3.connect(local_db)) as conn:
				conn.execute("PRAGMA journal_mode=OFF")
				conn.execute("PRAGMA synchronous=OFF")
				conn.execute("CREATE TABLE variants (chr INTEGER, pos INTEGER, a TEXT, b TEXT, rs INTEGER)")
				batch = []
				started = time.monotonic()
				with open(dbsnp, "rb", buffering=8 * 1024 * 1024) as compressed, gzip.open(compressed, "rt") as source:
					for row_number, line in enumerate(source, 1):
						if row_number % 5000000 == 0:
							print(
								f"MAGMA cache: scanned {row_number:,} dbSNP rows in {time.monotonic() - started:.0f}s",
								file=sys.stderr,
								flush=True,
							)
						ch, pos, snp, ref, alt = line.rstrip().split("\t")[:5]
						if not snp.startswith("rs") or not snp[2:].isdigit() or int(snp[2:]) not in ids:
							continue
						try:
							ch = chromosome(ch)
						except ValueError:
							continue
						if not 1 <= ch <= 23:
							continue
						for allele in alt.split(","):
							a, b = sorted((ref.upper(), allele.upper()))
							batch.append((ch, int(pos), a, b, int(snp[2:])))
						if len(batch) >= 10000:
							conn.executemany("INSERT INTO variants VALUES (?,?,?,?,?)", batch)
							batch.clear()
					conn.executemany("INSERT INTO variants VALUES (?,?,?,?,?)", batch)
				conn.execute("CREATE INDEX coordinates ON variants(chr,pos)")
				conn.commit()
			shutil.copyfile(local_db, temp)
			temp.replace(target)
			print("MAGMA rsID cache ready: " + str(target), file=sys.stderr, flush=True)
		finally:
			temp.unlink(missing_ok=True)
			workspace.cleanup()
	return target


def map_rows(rows, output, snploc, audit, cache):
	counts = dict(rsid=0, mapped=0, unmatched=0, ambiguous=0)
	conn = sqlite3.connect(Path(cache).resolve().as_uri() + "?mode=ro&immutable=1", uri=True) if cache else None
	if conn:
		# reference_cache publishes a closed database atomically and never edits
		# it. Immutable readers avoid per-variant locks/stat calls on DrvFS.
		# mmap shares cached pages between workers without a private DB copy.
		conn.execute("PRAGMA mmap_size=2147483648")
	previous = None
	matches = {}
	try:
		with open(rows) as source, open(output, "w") as out, open(snploc, "w") as loc:
			for line in source:
				snp, p, n, ch, pos, ea, nea = line.rstrip("\n").split("\t")
				ch, pos = chromosome(ch), int(pos)
				if snp.startswith("rs") and snp[2:].isdigit():
					counts["rsid"] += 1
				else:
					if conn is None:
						raise ValueError("Non-rsID input requires a dbSNP mapping cache")
					key = (ch, pos)
					if key != previous:
						matches = {}
						for a, b, rs in conn.execute("SELECT a,b,rs FROM variants WHERE chr=? AND pos=?", key):
							matches.setdefault((a, b), set()).add(rs)
						previous = key
					candidates = matches.get(tuple(sorted((ea.upper(), nea.upper()))), set())
					if len(candidates) != 1:
						counts["ambiguous" if candidates else "unmatched"] += 1
						continue
					snp = "rs" + str(next(iter(candidates)))
					counts["mapped"] += 1
				out.write(f"{snp}\t{p}\t{n}\n")
				loc.write(f"{snp}\t{ch}\t{pos}\n")
	finally:
		if conn:
			conn.close()
	with open(audit, "w") as out:
		out.write("status\trows\n")
		for key, value in counts.items():
			out.write(f"{key}\t{value}\n")
	print("MAGMA IDs: " + json.dumps(counts), file=sys.stderr)


def gwas_magma_ids_main():
	parser = argparse.ArgumentParser(description=__doc__)
	for name in ("rows", "output", "snploc", "audit", "bim", "dbsnp", "cache"):
		parser.add_argument("--" + name, required=True)
	args = parser.parse_args()
	with open(args.rows) as source:
		need_mapping = any(not (s := line.split("\t", 1)[0]).startswith("rs") or not s[2:].isdigit() for line in source)
	cache = reference_cache(args.bim, args.dbsnp, args.cache) if need_mapping else None
	map_rows(args.rows, args.output, args.snploc, args.audit, cache)


def gwas_magma_ids_cli():
	gwas_magma_ids_main()


# 🚩 prepare_yap2018_mpb
"""Prepare the audited GCST007020 source; this is not a general Y/Z decoder.

Y=first ID allele and Z=second ID allele was checked against 1000G EUR
frequencies in both orientations. Every recovered REF is checked against hg19.
The downloaded original is retained. Downstream use is through format.sh.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


SOURCE_SHA256 = "14d54ed2a090a86ea2e65b83d83489e4359e019c6e3b9cb5c2e0ab9526798190"
HEADER = "SNP CHR POS EA NEA EAF N BETA SE P LOG10P".split()


def decode(row, fasta):
	snp, c, pos, _, ea, nea, eaf, missing, beta, se, _, p = row
	reasons = []
	if c == "25":
		# PLINK 25 is XY; this study's source uses X coordinates for PAR.
		if not (60001 <= int(pos) <= 2699520 or 154931044 <= int(pos) <= 155260560):
			raise ValueError("Source CHR=25 outside GRCh37 X PAR: " + snp)
		c = "23"
		reasons.append("PLINK_XY_25_to_X_23")
	if {ea, nea} == {"Y", "Z"}:
		pair = id_alleles(snp, pos)
		if pair:
			if fasta.fetch(c, int(pos), len(pair[0])) != pair[0]:
				raise ValueError("ID reference allele disagrees with GRCh37: " + snp)
			alleles = dict(zip(("Y", "Z"), pair))
			ea, nea = alleles[ea], alleles[nea]
			reasons.append("Y_first_ID_allele_Z_second_ID_allele")
		else:
			ea = nea = "NA"
			reasons.append("symbolic_SV_without_sequence_or_END;association_retained_alleles_NA")
	elif not set(ea + nea) <= set("ACGT"):
		raise ValueError("Unexpected allele coding: " + snp)
	miss = float(missing)
	if not 0 <= miss < 1:
		raise ValueError("Invalid F_MISS: " + snp)
	n = str(int(205327 * (1 - miss) + 0.5))
	lp = format(-math.log10(float(p)), ".12g") if float(p) > 0 else "NA"
	return [snp, c, pos, ea, nea, eaf, n, beta, se, p, lp], reasons


def prepare_yap2018_mpb_run(a):
	src, dst, qc = Path(a.input), Path(a.output), Path(a.qc_prefix)
	if src.resolve() == dst.resolve():
		raise ValueError("Preserve the original source")
	if dst.exists():
		raise ValueError("Prepared output already exists; inspect its audit before replacing it")
	dst.parent.mkdir(parents=True, exist_ok=True)
	qc.parent.mkdir(parents=True, exist_ok=True)
	fasta = IndexedFasta(a.fasta)
	counts, chromosomes = Counter(), Counter()
	digest = hashlib.sha256()
	with tempfile.TemporaryDirectory(prefix="mpb.prepare.", dir=dst.parent) as tmpdir:
		tsv = Path(tmpdir) / "rows.tsv"
		audit = Path(tmpdir) / "audit.tsv.gz"
		with src.open("rb") as inp, tsv.open("w") as out, gzip.open(audit, "wt") as log:
			first = next(inp)
			digest.update(first)
			expected = "SNP CHR BP GENPOS ALLELE1 ALLELE0 A1FREQ F_MISS BETA SE P_BOLT_LMM_INF P_BOLT_LMM".split()
			if first.decode().split() != expected:
				raise ValueError("Unexpected source header")
			log.write("SNP\tCHR_SOURCE\tEA_SOURCE\tNEA_SOURCE\tCHR\tEA\tNEA\tREASON\n")
			for n, line in enumerate(inp, 1):
				digest.update(line)
				row = line.decode().rstrip("\r\n").split("\t")
				result, reasons = decode(row, fasta)
				out.write("\t".join(result) + "\n")
				counts["input"] += 1
				chromosomes[result[1]] += 1
				for reason in reasons:
					counts[reason] += 1
				if reasons:
					log.write(
						"\t".join([row[0], row[1], row[4], row[5], result[1], result[3], result[4], ";".join(reasons)])
						+ "\n"
					)
				if n % 2000000 == 0:
					print(f"Prepared {n:,} rows", flush=True)
		if digest.hexdigest() != SOURCE_SHA256:
			raise ValueError("Source differs from the audited GCST007020 download; no final output published")
		packed = Path(tmpdir) / "prepared.gz"
		with packed.open("wb") as out:
			bg = subprocess.Popen(["bgzip", "-@", "4", "-c"], stdin=subprocess.PIPE, stdout=out)
			try:
				bg.stdin.write(("\t".join(HEADER) + "\n").encode())
				with subprocess.Popen(
					["sort", "-T", tmpdir, "-S", "1G", "-t", "\t", "-k2,2n", "-k3,3n", "-k1,1", str(tsv)],
					stdout=subprocess.PIPE,
					env=dict(os.environ, LC_ALL="C"),
				) as sorter:
					shutil.copyfileobj(sorter.stdout, bg.stdin)
					if sorter.wait():
						raise RuntimeError("Sort failed")
				bg.stdin.close()
				if bg.wait():
					raise RuntimeError("BGZF compression failed")
			except BaseException:
				bg.kill()
				bg.wait()
				raise
		subprocess.run(["tabix", "-s", "2", "-b", "3", "-e", "3", "-S", "1", str(packed)], check=True)
		os.replace(packed, dst)
		os.replace(str(packed) + ".tbi", str(dst) + ".tbi")
		os.replace(audit, str(qc) + ".source.alleles.tsv.gz")
	fasta.close()
	Path(str(dst) + ".grch").write_text("37\n")
	report = dict(
		source=str(src),
		source_sha256=digest.hexdigest(),
		output=str(dst),
		counts=counts,
		chromosomes=chromosomes,
		build=37,
		n_total=205327,
		note="BETA, SE, EAF and P_BOLT_LMM preserved. Symbolic SV association rows retained with NA alleles.",
	)
	Path(str(qc) + ".source.prepare.json").write_text(json.dumps(report, indent=2) + "\n")
	print(json.dumps(report, indent=2), flush=True)


def prepare_yap2018_mpb_cli():
	p = argparse.ArgumentParser(description=__doc__)
	for arg in ("input", "output", "qc-prefix"):
		p.add_argument("--" + arg, required=True)
	p.add_argument("--fasta", default="/mnt/f/gen/1hgp/37/GRCH37.fasta.gz")
	prepare_yap2018_mpb_run(p.parse_args())


SUBCOMMANDS = {
	"liftover": gwas_liftover_cli,
	"h2": gwas_h2_cli,
	"magma-ids": gwas_magma_ids_cli,
	"yap2018-mpb": prepare_yap2018_mpb_cli,
}


def main():
	import sys

	if len(sys.argv) > 1 and sys.argv[1] in SUBCOMMANDS:
		return SUBCOMMANDS[sys.argv.pop(1)]()
	if len(sys.argv) < 2 or sys.argv[1] in ("--help", "-h"):
		print("Commands: " + ", ".join(SUBCOMMANDS))
		return
	raise SystemExit("Unknown helper command: " + sys.argv[1])


if __name__ == "__main__":
	main()
