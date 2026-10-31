#!/usr/bin/env python3
"""ibdmix utilities and command entry points."""
from __future__ import annotations
import sys

sys.dont_write_bytecode = True


# 🚩 ibdmix
# GU ibdmix utilities. See --help for commands.
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

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
from gu_0_common import load_module, x_nonpar_intervals


# 🚩 ibdmix_final_filter
"""Non-destructive Cell 2020 postprocessing for GU final.

The published text does not specify the DAF window construction. Accept a
provenanced exclusion BED; never substitute segment LOD or carrier frequency
for archaic derived-allele proportion. All intervals here are BED0 half-open.
"""

import csv, gzip, hashlib, json, math, os, tempfile
from collections import defaultdict, Counter
from bisect import bisect_right

ibdmix_final_filter_VERSION = "2026-09-15.1"
AFR = {"ESN", "GWD", "LWK", "MSL", "YRI"}
NEA = {"Altai", "Altai.2013", "Vindija", "Chagyr", "Chagyrskaya"}
DEN = {"Denisova", "Denisova.2013", "Denisova25"}


def rows(path):
	with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as f:
		yield from csv.DictReader(f, delimiter="\t")


def ibdmix_final_filter_union(intervals):
	result = []
	for lo, hi in sorted(intervals):
		if lo < 0 or hi <= lo:
			raise ValueError("Invalid interval")
		if result and lo <= result[-1][1]:
			result[-1] = (result[-1][0], max(hi, result[-1][1]))
		else:
			result.append((lo, hi))
	return result


def ibdmix_final_filter_subtract(lo, hi, mask):
	i = max(0, bisect_right(mask, (lo, math.inf)) - 1)
	for j in range(i, len(mask)):
		left, right = mask[j]
		if left >= hi:
			break
		if right <= lo:
			continue
		if left > lo:
			yield lo, left
		lo = max(lo, right)
		if lo >= hi:
			return
	if lo < hi:
		yield lo, hi


def carrier_mask(calls, members, threshold=0.30):
	"""Pointwise pooled AFR5 carrier frequency, one vote per person.

    Include zero-call participants in the denominator. This is a carrier
    fraction, not an inferred phased allele frequency (IBDmix is unphased).
    """
	if not members:
		raise ValueError("No African denominator")
	changes = Counter()
	for sample, intervals in calls.items():
		if sample not in members:
			raise ValueError("Control sample outside denominator")
		for lo, hi in ibdmix_final_filter_union(intervals):
			changes[lo] += 1
			changes[hi] -= 1
	count = 0
	previous = None
	result = []
	minimum = math.ceil(threshold * len(members))
	for pos, change in sorted(changes.items()):
		if previous is not None and count >= minimum and pos > previous:
			result.append((previous, pos))
		count += change
		previous = pos
	return ibdmix_final_filter_union(result)


def checksum(path):
	h = hashlib.sha256()
	with open(path, "rb") as f:
		while data := f.read(1024 * 1024):
			h.update(data)
	return h.hexdigest()


def load_daf(path):
	"""BED plus .json declaring build, reference, percentile and provenance."""
	if path is None:
		return {}, {"status": "not_applied_missing_author_mask"}
	path = Path(path)
	meta = json.loads(Path(str(path) + ".json").read_text())
	if meta.get("genome_build") != "GRCh37" or meta.get("reference") != "Altai":
		raise ValueError("DAF mask requires genome_build=GRCh37 and reference=Altai")
	if meta.get("percentile") != 99.9 or not meta.get("source"):
		raise ValueError("DAF mask requires percentile=99.9 and a source citation")
	masks = defaultdict(list)
	with gzip.open(path, "rt") if str(path).endswith(".gz") else path.open() as f:
		for line in f:
			if not line.strip() or line.startswith(("#", "track ", "browser ")):
				continue
			ch, lo, hi, *_ = line.split()
			ch = ch.removeprefix("chr")
			if ch not in set(map(str, range(1, 23))):
				raise ValueError("DAF mask must be autosomal")
			masks[ch].append((int(lo), int(hi)))
	if not masks:
		raise ValueError("DAF mask is empty")
	return {ch: ibdmix_final_filter_union(iv) for ch, iv in masks.items()}, dict(
		meta, status="applied", path=str(path.resolve()), sha256=checksum(path)
	)


def raw_masks(run, chrom, refs, altai_ref="Altai"):
	"""Use original Neanderthal calls, before the AFR-Denisova subtraction."""
	nea = []
	controls = {ref: defaultdict(list) for ref in refs if ref in DEN}
	members = set()
	seen_pops = set()
	input_paths = []
	units = sorted((run / "samples").glob("*/populations.tsv"))
	if not units:
		raise ValueError(f"IBDmix sample denominators missing: {run}")
	for table in units:
		populations = list(rows(table))
		unit = table.parent.name
		for pop in populations:
			name = pop["population"]
			sample_file = table.parent / (name + ".txt")
			ids = set(sample_file.read_text().split())
			input_paths.append(sample_file)
			if len(ids) != int(pop["n"]):
				raise ValueError(f"Population denominator mismatch: {sample_file}")
			if name in AFR:
				members.update(ids)
				seen_pops.add(name)
			# Altai defines the paper comparison. Denisova25 uses the same
			# Altai overlap mask as an explicitly labelled extension.
			for ref in [altai_ref] + sorted(controls):
				if ref != altai_ref and name not in AFR:
					continue
				path = run / "raw" / unit / f"{ref}.{name}.raw.txt.gz"
				if not path.is_file():
					raise ValueError(f"Original IBDmix calls required by final filter: {path}")
				input_paths.append(path)
				for r in rows(path):
					ch = r["chrom"].removeprefix("chr")
					ch = "X" if ch == "23" else ch
					if ch != chrom:
						continue
					if r["ID"] not in ids:
						raise ValueError(f"Wrong population sample: {path}")
					lo, hi = int(r["start"]) - 1, int(r["end"]) - 1
					score = float(r.get("LOD", r.get("slod", "nan")))
					if not math.isfinite(score):
						raise ValueError(f"Invalid LOD: {path}")
					if score < 4 or hi - lo < 50000:
						continue
					if ref == altai_ref:
						nea.append((lo, hi))
					else:
						controls[ref][r["ID"]].append((lo, hi))
		input_paths.append(table)
	if seen_pops != AFR:
		raise ValueError(f"All five African controls required: {run}")
	return (
		ibdmix_final_filter_union(nea),
		{ref: carrier_mask(data, members) for ref, data in controls.items()},
		len(members),
	)


def filter_rows(source, masks, daf, chrom, counts):
	for row in rows(source):
		ref = row["anc"]
		lo, hi = int(row["start"]), int(row["end"])
		if row["chrom"].removeprefix("chr") != chrom:
			raise ValueError("Mixed chromosome final output")
		counts[ref]["input_segments"] += 1
		counts[ref]["input_bp"] += hi - lo
		mask = daf.get(chrom, []) if ref == "Altai" else masks.get(ref, [])
		for left, right in ibdmix_final_filter_subtract(lo, hi, mask):
			if right - left < 50000:
				counts[ref]["short_fragment_bp"] += right - left
				continue
			result = dict(row, start=left, end=right, length=right - left)
			if (left, right) != (lo, hi):
				result["score_scope"] = "parent_native_call_after_final_mask"
				# Parent LOD is retained; site counts do not describe fragments.
				for key in ("sites", "positive_lods", "negative_lods"):
					if key in result:
						result[key] = ""
			counts[ref]["output_segments"] += 1
			counts[ref]["output_bp"] += right - left
			yield result


def prepare(source, output_root, meta, daf_path=None):
	"""Return a filtered snapshot and audit; native chromosome files stay intact."""
	source = Path(source)
	run = source.parent.parent
	daf, daf_info = load_daf(daf_path)
	with gzip.open(source, "rt") as f:
		reader = csv.DictReader(f, delimiter="\t")
		fields = reader.fieldnames
		first = next(reader, None)
	if first is None:
		return source, {"status": "empty_native_result", "daf": daf_info}
	chrom = first["chrom"].removeprefix("chr")
	if chrom not in set(map(str, range(1, 23))):
		return source, {"status": "extension_chrX_not_filtered", "daf": {"status": "not_applicable"}}
	if (
		meta.get("genome_build") not in ("b37", "37", "GRCh37")
		or meta.get("coordinate_system") != "0-based-half-open"
	):
		raise ValueError(f"Cell final filters require GRCh37 BED0 calls: {source}")
	refs = set(meta.get("refs", "").split())
	if refs & {"Denisova", "Denisova25"} and "Altai" not in refs:
		raise ValueError("Denisovan filtering requires original Altai calls")
	if "Denisova.2013" in refs and "Altai.2013" not in refs:
		raise ValueError("Denisova.2013 final filtering requires original Altai.2013 calls")
	dependencies = [
		source,
		run / "run.meta.tsv",
		*sorted((run / "raw").glob("*/*.raw.txt.gz")),
		*sorted((run / "samples").glob("*/*.txt")),
		*sorted((run / "samples").glob("*/populations.tsv")),
	]
	signature = dict(
		version=ibdmix_final_filter_VERSION,
		code_sha256=checksum(__file__),
		daf=daf_info,
		run_meta=meta,
		inputs=[[str(p.resolve()), p.stat().st_size, p.stat().st_mtime_ns] for p in dependencies],
	)
	key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
	dest = Path(output_root) / key
	result = dest / "segments.tsv.gz"
	manifest = dest / "audit.json"
	if result.is_file() and manifest.is_file():
		saved = json.loads(manifest.read_text())
		if saved["signature"] == signature and saved["output_sha256"] == checksum(result):
			return result, saved
	Path(output_root).mkdir(parents=True, exist_ok=True)
	with tempfile.TemporaryDirectory(prefix=".filter-", dir=output_root) as tmp:
		stage = Path(tmp)
		overlap, high, denominator = (
			raw_masks(run, chrom, refs - {"Denisova.2013"})
			if refs & {"Denisova", "Denisova25"}
			else ([], {}, 0)
		)
		masks = {ref: ibdmix_final_filter_union(overlap + iv) for ref, iv in high.items()}
		overlap2013 = []
		if "Denisova.2013" in refs:
			overlap2013, high2013, denominator = raw_masks(
				run, chrom, {"Denisova.2013"}, altai_ref="Altai.2013"
			)
			masks.update({ref: ibdmix_final_filter_union(overlap2013 + iv) for ref, iv in high2013.items()})
		for name, iv in [
			("Altai_overlap", overlap),
			("Altai.2013_overlap", overlap2013),
			*[(r + "_AFR30", iv) for r, iv in high.items()],
			("Altai_DAF", daf.get(chrom, [])),
		]:
			(stage / (name + ".bed")).write_text("".join(f"{chrom}\t{lo}\t{hi}\n" for lo, hi in iv))
		counts = defaultdict(Counter)
		with gzip.open(stage / "segments.tsv.gz", "wt") as f:
			w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
			w.writeheader()
			w.writerows(filter_rows(source, masks, daf, chrom, counts))
		audit = dict(
			signature=signature,
			status="filtered",
			daf=daf_info,
			chrom=chrom,
			denisovan_overlap="union_of_original_Altai_calls_across_all_populations; 2013 pair uses Altai.2013",
			african_frequency="pointwise_distinct_carriers_in_pooled_ESN_GWD_LWK_MSL_YRI / all_tested_members >= 0.30",
			n_african_controls=denominator,
			reference_scope="Altai/Denisova and Altai.2013/Denisova.2013 independent pairs; Denisova25 uses Altai extension",
			counts=dict(counts),
			output_sha256=checksum(stage / "segments.tsv.gz"),
		)
		(stage / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
		# A process-specific cache generation avoids replacing active readers.
		if dest.exists():
			dest = Path(str(dest) + f".{os.getpid()}")
		os.replace(stage, dest)
	return dest / "segments.tsv.gz", audit


# 🚩 ibdmix_provenance
"""Compare IBDmix run contracts using input paths exactly as recorded.

pipeline_version in ibdmix.sh is the semantic contract for the workflow scripts:
changes to genotype preparation, calling or filtering MUST bump that version.
Script file timestamps and mask manifests remain audit information. Actual mask
BED checksums and all input stat records remain part of reuse validation.
"""

import argparse
import difflib

WORKFLOW_SCRIPTS = {
	str(Path(__file__).with_name(name).resolve()) for name in ("ibdmix.sh", "ibdmix.py", "0.common.py")
}


def contract(text):
	# GNU stat's format in older run records emitted literal backslash-t.
	rows = [line.replace("\\t", "\t").split("\t") for line in text.splitlines() if line]
	refs = next((row[1].split() for row in rows if row[0] == "refs"), [])
	versions = [row for row in rows if row[0] == "pipeline_version"]
	masks = [row for row in rows if row[0] == "excluded_mask_sha256"]
	units = {row[1] for row in masks if len(row) == 4}
	units.update(row[1] for row in rows if row[0] == "mask_manifest_sha256" and len(row) >= 2)
	if len(versions) != 1 or len(versions[0]) != 2 or not versions[0][1] or not refs or not units:
		raise ValueError(
			"Incomplete provenance: pipeline version, references and mask checksums are required"
		)
	expected = {(unit, ref) for unit in units for ref in refs}
	actual = {
		(row[1], row[2])
		for row in masks
		if len(row) == 4 and len(row[3]) == 64 and all(c in "0123456789abcdef" for c in row[3])
	}
	if actual != expected or len(masks) != len(expected):
		raise ValueError("Incomplete or duplicate excluded-mask checksums")
	result = []
	for row in rows:
		key = row[0]
		if key in {"mask_root", "mask_manifest_sha256"}:
			continue  # Actual BED hashes, rather than their storage metadata.
		if key == "software" and len(row) == 2:
			path = row[1].rsplit(":", 2)[0]
			if path in WORKFLOW_SCRIPTS:
				continue  # Explicit pipeline_version controls semantic compatibility.
		result.append("\t".join(row))
	return sorted(result)


def differences(previous, current):
	before, after = contract(previous), contract(current)

	def references(lines):
		return set(next(line.split("\t")[1].split() for line in lines if line.startswith("refs\t")))

	common = references(before) & references(after)

	def project(lines, refs):
		result = []
		for line in lines:
			row = line.split("\t")
			if row[0] == "refs":
				continue
			if row[0] == "excluded_mask_sha256" and row[2] not in common:
				continue
			if row[0] in {"archaic_vcf", "excluded_mask"}:
				path = row[1].rsplit(":", 2)[0]
				owners = set(Path(path).parts) & refs
				if len(owners) != 1:
					raise ValueError(f"Cannot identify reference in provenance input: {path}")
				if not owners & common:
					continue
			result.append(line)
		return result

	before, after = project(before, references(before)), project(after, references(after))
	return "".join(
		difflib.unified_diff(
			[s + "\n" for s in before],
			[s + "\n" for s in after],
			fromfile="previous run contract",
			tofile="current run contract",
		)
	)


def merge_cache(previous, current):
	"""Retain validated per-reference contracts when references are deselected."""
	before, after = contract(previous), contract(current)

	def refs(lines):
		return set(next(line.split("\t")[1].split() for line in lines if line.startswith("refs\t")))

	old_refs, new_refs = refs(before), refs(after)
	retained = old_refs - new_refs
	result = [line for line in after if not line.startswith("refs\t")]
	for line in before:
		row = line.split("\t")
		if row[0] == "excluded_mask_sha256" and row[2] in retained:
			result.append(line)
		elif row[0] in {"archaic_vcf", "excluded_mask"}:
			path = row[1].rsplit(":", 2)[0]
			if set(Path(path).parts) & retained:
				result.append(line)
	result.append("refs\t" + " ".join(sorted(old_refs | new_refs)))
	merged = "\n".join(sorted(result)) + "\n"
	contract(merged)  # Validate before committing the retained contracts.
	return merged


def ibdmix_provenance_main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--write-cache", action="store_true")
	parser.add_argument("previous", type=Path)
	parser.add_argument("current", type=Path)
	args = parser.parse_args()
	try:
		if args.write_cache:
			current = args.current.read_text()
			merged = merge_cache(args.previous.read_text(), current) if args.previous.exists() else current
			staged = args.previous.with_name(args.previous.name + ".part")
			staged.write_text(merged)
			staged.replace(args.previous)
			return 0
		diff = differences(args.previous.read_text(), args.current.read_text())
	except ValueError as error:
		print(str(error), file=sys.stderr)
		return 1
	if diff:
		print(diff, end="", file=sys.stderr)
		return 1
	return 0


def ibdmix_provenance_cli():
	sys.exit(ibdmix_provenance_main())


# 🚩 ibdmix_workflow
"""IBDmix preparation and final filtering following the authors' Cell workflow.

The caller remains the unmodified upstream executable. BED files here are
0-based half-open; native IBDmix coordinates are converted exactly once.
"""
from contextlib import contextmanager
import fcntl
import re
import subprocess
import numpy as np

load_module("0.common.py")
from gu_0_common import CHROM_LENGTHS, write_tsv_rows

ibdmix_workflow_VERSION = "2026-09-14.1"
AFRICAN_CONTROLS = {"ESN", "GWD", "LWK", "MSL", "YRI"}
NEANDERTHALS = {"Altai", "Altai.2013", "Chagyr", "Chagyrskaya", "Vindija"}
DENISOVANS = {"Denisova", "Denisova.2013", "Denisova25"}


REFERENCE_ALIASES = {ref.lower(): ref for ref in sorted(NEANDERTHALS | DENISOVANS)}
REFERENCE_ALIASES.update(vindijia="Vindija", denisovan="Denisova", den25="Denisova25", chagyrskaya="Chagyr")
CALL_REFERENCES = ("Altai", "Chagyr", "Denisova", "Denisova25", "Vindija")


def normalize_refs(value):
	tokens = value.replace(",", " ").split()
	if not tokens or ",," in value or value.startswith(",") or value.endswith(","):
		raise ValueError("Archaic reference list must be nonempty without empty entries")
	try:
		return list(dict.fromkeys(REFERENCE_ALIASES[token.lower()] for token in tokens))
	except KeyError as e:
		raise ValueError(f"Unknown archaic reference: {e.args[0]}") from None


def opener(path):
	return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def fingerprint(path):
	p = Path(path).resolve()
	s = p.stat()
	return [str(p), s.st_size, s.st_mtime_ns]


def file_checksum(path, algorithm="sha256"):
	h = hashlib.new(algorithm)
	with open(path, "rb") as handle:
		while block := handle.read(1024 * 1024):
			h.update(block)
	return h.hexdigest()


def sha256_file(path):
	return file_checksum(path)


def samples(panel, actual, output, background_required=True):
	metadata = {}
	with open(panel) as handle:
		header = handle.readline().lower().split()

		def index(names):
			return next((header.index(n) for n in names if n in header), None)

		sample_col = index(["sample", "sample_id", "iid", "#iid"])
		pop_col = index(["pop", "population"])
		super_col = index(["super_pop", "super_population", "superpop"])
		if sample_col is None or pop_col is None:
			raise ValueError(
				"IBDmix requires sample and pop columns; provide --sample-panel with population labels."
			)
		for line in handle:
			fields = line.split()
			if not fields:
				continue
			sample = fields[sample_col]
			pop = fields[pop_col]
			if sample in metadata:
				raise ValueError(f"Duplicate metadata sample: {sample}")
			if not re.fullmatch(r"[A-Za-z0-9_-]+", pop) or pop.upper() in {"NA", "NONE", "UNKNOWN", "ALL"}:
				raise ValueError(f"Missing or unsafe population label for {sample}: {pop}")
			metadata[sample] = (pop, fields[super_col] if super_col is not None else "UNKNOWN")
	ids = Path(actual).read_text().split()
	if len(set(ids)) != len(ids) or not ids:
		raise ValueError("Target sample IDs must be nonempty and unique")
	groups = {}
	for sample in ids:
		if sample not in metadata:
			raise ValueError(f"Target sample missing from metadata: {sample}")
		pop, superpop = metadata[sample]
		groups.setdefault(pop, []).append(sample)
	if background_required and not AFRICAN_CONTROLS.issubset(groups):
		raise ValueError(
			"Cell 2020 filtering needs all five African control populations: ESN GWD LWK MSL YRI. Use a complete 1KG cohort, or an explicitly configured custom profile."
		)
	output = Path(output)
	output.mkdir(parents=True, exist_ok=True)
	rows = []
	for pop, members in sorted(groups.items()):
		if len(members) < 10:
			raise ValueError(
				f"Population {pop} has {len(members)} targets; IBDmix requires at least 10 individuals."
			)
		path = output / f"{pop}.txt"
		path.write_text("".join(s + "\n" for s in members))
		supers = {metadata[s][1] for s in members}
		if len(supers) != 1:
			raise ValueError(f"Conflicting super-population labels for {pop}")
		rows.append(
			dict(
				population=pop,
				super_population=next(iter(supers)),
				n=len(members),
				sample_file=str(path.resolve()),
			)
		)
	write_tsv_rows(output / "populations.tsv", ["population", "super_population", "n", "sample_file"], rows)
	Path(output / "ALL.txt").write_text("".join(s + "\n" for s in ids))
	return rows


def bed_intervals(path, chrom, length, columns=(0, 1, 2)):
	with opener(path) as handle:
		for line in handle:
			if line.startswith(("#", "track", "browser")) or not line.strip():
				continue
			f = line.split()
			if f[columns[0]].removeprefix("chr") != chrom:
				continue
			lo, hi = int(f[columns[1]]), int(f[columns[2]])
			if not 0 <= lo < hi <= length:
				raise ValueError(f"Invalid BED interval in {path}: {line.strip()}")
			yield lo, hi


def write_boolean_bed(path, chrom, mask):
	# Emit runs without allocating one Python object per masked base.
	edges = np.flatnonzero(np.diff(mask.astype(np.int8), prepend=0, append=0))
	with open(path, "w") as out:
		for lo, hi in zip(edges[::2], edges[1::2]):
			out.write(f"{chrom}\t{lo}\t{hi}\n")


def cpg_sites(reference, variants, alignments, chrom):
	"""Potential CpGs in hg19, 1KG alleles and the three ancestral alignments.

    Same union-of-C/G rule as upstream generate_cpg_mask.py, vectorized per
    chromosome. AXT target starts are 1-based inclusive.
    """
	sequence = np.frombuffer(reference.upper(), dtype="S1")
	can_c = sequence == b"C"
	can_g = sequence == b"G"
	with opener(variants) as handle:
		for line in handle:
			ch, pos, ref, alt = line.split()
			if ch.removeprefix("chr") != chrom:
				continue
			i = int(pos) - 1
			if not 0 <= i < len(sequence):
				raise ValueError("CpG variant outside chromosome")
			if ref not in {"N", "M"} and sequence[i] not in (ref.encode(), b"N", b"M"):
				raise ValueError(f"CpG variant/reference mismatch at chr{chrom}:{pos}")
			if alt == "C":
				can_c[i] = True
			if alt == "G":
				can_g[i] = True
	for path in alignments:
		with opener(path) as handle:
			for line in handle:
				if not line.strip() or line.startswith("#"):
					continue
				header = line.split()
				human = handle.readline().strip().upper()
				other = handle.readline().strip().upper()
				if len(header) != 9 or len(human) != len(other):
					raise ValueError(f"Malformed AXT: {path}")
				if header[1].removeprefix("chr") != chrom:
					raise ValueError(f"Wrong AXT chromosome: {path}")
				h = np.frombuffer(human.encode(), dtype="S1")
				o = np.frombuffer(other.encode(), dtype="S1")
				keep = h != b"-"
				mapped = o[keep]
				start = int(header[2]) - 1
				end = start + len(mapped)
				if end != int(header[3]) or not 0 <= start < end <= len(sequence):
					raise ValueError(f"Invalid AXT coordinates: {path}")
				can_c[start:end] |= mapped == b"C"
				can_g[start:end] |= mapped == b"G"
	pairs = can_c[:-1] & can_g[1:]
	result = np.zeros(len(sequence), dtype=bool)
	result[:-1] |= pairs
	result[1:] |= pairs
	return result


def axt_directory():
	annot = Path(os.environ.get("GU_ANNOT_ROOT", "/mnt/f/annot"))
	return Path(os.environ.get("IBDMIX_AXT_DIR", str(annot / "axt/37")))


@contextmanager
def cached_mask_artifacts(root, signature, names):
	"""Serialize builders and publish complete, checksum-verified cache entries."""
	root = Path(root)
	root.mkdir(parents=True, exist_ok=True)
	key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
	destination = root / key
	with (root / f"{key}.lock").open("w") as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		try:
			meta = json.loads((destination / "cache.json").read_text())
			valid = meta["signature"] == signature and all(
				(destination / name).is_file() and sha256_file(destination / name) == meta["sha256"][name]
				for name in names
			)
		except (OSError, ValueError, KeyError, TypeError):
			valid = False
		if valid:
			yield destination, None
			return
		with tempfile.TemporaryDirectory(prefix=f".{key}.", dir=root) as temporary:
			staging = Path(temporary)
			yield destination, staging
			checksums = {name: sha256_file(staging / name) for name in names}
			(staging / "cache.json").write_text(
				json.dumps(dict(signature=signature, sha256=checksums), indent=2) + "\n"
			)
			destination.mkdir(exist_ok=True)
			for name in [*names, "cache.json"]:
				os.replace(staging / name, destination / name)


def cached_cpg_sites(fasta, variants, alignments, chrom, length, root):
	"""Cache the CpG bitmap under the analysis output, in packed-bit form."""
	reference_files = [fasta, str(fasta) + ".fai"]
	if str(fasta).endswith((".gz", ".bgz", ".bgzf")):
		reference_files.append(str(fasta) + ".gzi")
	signature = dict(
		schema=1,
		chrom=chrom,
		length=length,
		code=sha256_file(__file__),
		inputs=[fingerprint(p) for p in [*reference_files, variants, *alignments]],
	)
	with cached_mask_artifacts(Path(root) / "cpg" / f"chr{chrom}", signature, ["cpg.exclude.packed.npy"]) as (
		cache,
		staging,
	):
		if staging is None:
			packed = np.load(cache / "cpg.exclude.packed.npy", allow_pickle=False)
			if packed.dtype != np.uint8 or packed.shape != ((length + 7) // 8,):
				raise ValueError(f"Invalid cached CpG bitmap: {cache}")
			print(f"IBDMIX chr{chrom}: reusing analysis CpG cache {cache}", flush=True)
			return np.unpackbits(packed, count=length).astype(bool)
		with open(str(fasta) + ".fai") as handle:
			contigs = {line.split()[0]: int(line.split()[1]) for line in handle}
		contig = chrom if chrom in contigs else "chr" + chrom
		if contigs.get(contig) != length:
			raise ValueError("Reference FASTA chromosome length does not match GRCh37")
		fa = subprocess.check_output(["samtools", "faidx", str(fasta), contig])
		sequence = b"".join(fa.splitlines()[1:])
		del fa
		if len(sequence) != length:
			raise ValueError("Incomplete FASTA sequence")
		mask = cpg_sites(sequence, variants, alignments, chrom)
		np.save(staging / "cpg.exclude.packed.npy", np.packbits(mask), allow_pickle=False)
		print(f"IBDMIX chr{chrom}: built analysis CpG cache {cache}", flush=True)
	return mask


def link_cached_masks(cache, output, names):
	"""Link each analysis unit to generated masks within the analysis tree."""
	output = Path(output)
	output.mkdir(parents=True, exist_ok=True)
	with tempfile.TemporaryDirectory(prefix=".mask-links.", dir=output) as temporary:
		for name in names:
			link = Path(temporary) / name
			link.symlink_to((cache / name).resolve())
			os.replace(link, output / name)


def axt_reference(chrom, species, root=None):
	filename = f"chr{chrom}.hg19.{species}.synNet.axt.gz"
	url = f"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/vs{species[0].upper() + species[1:]}/syntenicNet/{filename}"
	# Explicit roots retain the flat-directory interface. The default cache
	# separates species, matching the user's 37/vsPanTro2 layout.
	directory = Path(root or axt_directory())
	if root is None and not os.environ.get("IBDMIX_AXT_DIR"):
		directory /= f"vs{species[0].upper() + species[1:]}"
	return url, directory / filename


def local_axt(chrom, species, root=None):
	"""Use only a completed local AXT matching its publisher MD5."""
	_, path = axt_reference(chrom, species, root)
	checksums = path.parent / "md5sum.txt"
	if not local_mask_ready(checksums):
		checksums = path.parent / "checksums" / f"{species}.md5sum.txt"
	for required in [path, checksums]:
		if not local_mask_ready(required):
			raise FileNotFoundError(
				f"Missing, empty or unfinished AXT resource: {required}; "
				"automatic reference downloads are disabled. Download it manually."
			)
	expected = {}
	for line in checksums.read_text().splitlines():
		if not line.strip():
			continue
		digest, name = line.split()
		expected[Path(name.lstrip("*")).name] = digest
	if path.name not in expected:
		raise ValueError(f"AXT absent from publisher checksum list: {path.name}")
	if file_checksum(path, "md5") != expected[path.name]:
		raise ValueError(f"Publisher checksum mismatch: {path}")
	return path


def published_mask_name(ref, chrom):
	"""The publisher's included-base mask for this individual and chromosome."""
	if ref == "Chagyr":
		return f"Chagyr/chr{chrom}_mask.bed.gz"
	elif ref == "Vindija":
		return f"Vindija/chr{chrom}_mask.bed.gz"
	elif ref == "Denisova25":
		suffix = "" if chrom == "X" else ".min10x"
		filename = f"Denisova25.chr{chrom}.hg19.L35MQ25.map35_100.GCcov.noSimpleRepeat.noIndel{suffix}.bed.gz"
		return "Denisova25/" + filename
	else:
		raise ValueError(f"No published mask configured for {ref}")


MINIMAL_MASK_NAMES = {
	"Altai": "Altai/AltaiNea.map35_50.MQ30.Cov.indels.TRF.bed.bgz",
	"Denisova": "Denisova/DenisovaPinky.map35_50.MQ30.Cov.indels.TRF.bed.bgz",
	"Altai.2013": "Altai/AltaiNea.map35_50.MQ30.Cov.indels.TRF.bed.bgz",
	"Denisova.2013": "Denisova/DenisovaPinky.map35_50.MQ30.Cov.indels.TRF.bed.bgz",
}


def local_mask_candidates(local_name, assets, archaic_root=None):
	candidates = [Path(assets) / local_name]
	if archaic_root:
		candidates.append(Path(archaic_root).parent / "mask" / local_name)
	return list(dict.fromkeys(candidates))


def local_mask_ready(path, allow_empty=False):
	path = Path(path)
	# aria2 preallocates its final filename before the transfer is complete.
	return (
		path.is_file()
		and (allow_empty or path.stat().st_size > 0)
		and not Path(str(path) + ".aria2").exists()
	)


def require_local_mask(candidates, allow_empty=False):
	for path in candidates:
		if local_mask_ready(path, allow_empty):
			return Path(path)
	raise FileNotFoundError(
		"Missing, empty or unfinished mask; automatic mask downloads are disabled. "
		"Expected one of: " + ", ".join(map(str, candidates))
	)


def local_reference(path):
	"""Read an existing mask resource without any network fallback."""
	path = require_local_mask([Path(path)])
	provenance = Path(str(path) + ".source.json")
	if provenance.exists() and sha256_file(path) != json.loads(provenance.read_text())["sha256"]:
		raise ValueError(f"Reference checksum mismatch: {path}")
	return path


def published_reference_mask(ref, chrom, assets, archaic_root=None):
	return reference_mask_resource(published_mask_name(ref, chrom), assets, archaic_root)


def reference_mask_resource(relative_name, root, archaic_root=None):
	path = require_local_mask(local_mask_candidates(relative_name, root, archaic_root))
	return local_reference(path)


def check_mask_inputs(root, chroms, refs, archaic_root=None, custom_masks=None):
	"""List every missing mask before target export or reference downloads."""
	refs = normalize_refs(" ".join(refs))
	chroms = list(
		dict.fromkeys(
			"X" if str(c).removeprefix("chr") == "23" else str(c).removeprefix("chr") for c in chroms
		)
	)
	if not chroms or any(c not in CHROM_LENGTHS["37"] for c in chroms):
		raise ValueError("Mask checks require supported chromosomes")
	assets = Path(root)
	required = []
	if custom_masks:
		required = [[Path(custom_masks) / ref / f"chr{chrom}.bed"] for chrom in chroms for ref in refs]
	else:
		if any(c != "X" for c in chroms):
			required.append(
				local_mask_candidates("common/1kg.strict_mask.autosomes.bed", assets, archaic_root)
			)
		required.append(local_mask_candidates("common/genomicSuperDups.txt.gz", assets, archaic_root))
		required.extend(
			local_mask_candidates(MINIMAL_MASK_NAMES[ref], assets, archaic_root)
			for ref in refs
			if ref in MINIMAL_MASK_NAMES
		)
		required.extend(
			local_mask_candidates(published_mask_name(ref, chrom), assets, archaic_root)
			for chrom in chroms
			for ref in refs
			if ref not in MINIMAL_MASK_NAMES
		)
	missing = []
	for candidates in required:
		try:
			require_local_mask(candidates, allow_empty=bool(custom_masks))
		except FileNotFoundError:
			missing.append("  " + " OR ".join(map(str, candidates)))
	if missing:
		raise FileNotFoundError(
			"Required IBDmix mask files are missing, empty or unfinished (.aria2):\n"
			+ "\n".join(missing)
			+ "\nAutomatic mask downloads are disabled; "
			"finish manual downloads before rerunning."
		)
	return len(required)


def prepare_masks(root, chrom, modern, fasta, upstream, output, axt_root=None, refs=None, archaic_root=None):
	"""Build hg19/GRCh37 masks, with an explicit non-blocking X extension.

    map35_50 mappability is already present in the published minimal masks;
    recomputing and intersecting that identical mask is unnecessary.
    """
	length = CHROM_LENGTHS["37"][chrom]
	root = Path(root)
	output = Path(output)
	output.mkdir(parents=True, exist_ok=True)
	assets = root
	refs = refs or ["Altai", "Denisova"]
	check_mask_inputs(root, [chrom], refs, archaic_root)
	# This published accessibility BED contains autosomes only. Applying its
	# empty X subset would exclude all of X, independently of the AXT issue.
	strict = (
		None
		if chrom == "X"
		else reference_mask_resource("common/1kg.strict_mask.autosomes.bed", assets, archaic_root)
	)
	dups = reference_mask_resource("common/genomicSuperDups.txt.gz", assets, archaic_root)
	minimal = {}
	for ref, filename in MINIMAL_MASK_NAMES.items():
		if ref not in refs:
			continue
		minimal[ref] = reference_mask_resource(filename, assets, archaic_root)
	extra = {
		ref: published_reference_mask(ref, chrom, assets, archaic_root) for ref in refs if ref not in minimal
	}
	alignments = []
	species_list = ["panTro2", "ponAbe2", "rheMac2"]
	missing_axt = [
		axt_reference(chrom, species, axt_root)[1]
		for species in species_list
		if not local_mask_ready(axt_reference(chrom, species, axt_root)[1])
	]
	cpg_enabled = chrom != "X" or not missing_axt
	if cpg_enabled:
		for species in species_list:
			alignments.append(local_axt(chrom, species, axt_root))
		print(
			f"IBDMIX chr{chrom}: CpG filter enabled; verified AXT: " + ", ".join(map(str, alignments)),
			flush=True,
		)
	else:
		print(
			"IBDMIX chrX: CpG filter SKIPPED (X AXT unavailable); continuing X analysis with archaic minimal, segmental-duplication and modern-indel masks.",
			flush=True,
		)
	variants = Path(upstream) / "j.cell.2020.01.012-workflow/abridged_variants.gz"
	if cpg_enabled and not variants.exists():
		raise ValueError(f"Upstream CpG variant resource missing: {variants}")
	source_meta = Path(os.environ.get("GU_TARGET_TMP_DIR", "/nonexistent")) / f"chr{chrom}" / "source.tsv"
	modern_signature = sha256_file(source_meta) if source_meta.exists() else fingerprint(modern)
	inputs = [dups, *minimal.values(), *extra.values()]
	if strict is not None:
		inputs.append(strict)
	if cpg_enabled:
		inputs.extend([*alignments, variants, fasta])
	signature = dict(
		version=ibdmix_workflow_VERSION,
		chrom=chrom,
		modern=modern_signature,
		refs=refs,
		cpg_enabled=cpg_enabled,
		inputs=[fingerprint(p) for p in inputs],
	)
	# Final masks also depend on modern indels; do not share them by chromosome alone.
	signature["cache_schema"] = 1
	signature["code"] = sha256_file(__file__)
	# Reference resources are read-only inputs. Generated masks belong beside
	# the analysis units, including when the caller overrides --root.
	cache_root = output.parent / "derived"
	analysis_output = output
	names = [f"{ref}.exclude.bed" for ref in refs] + ["manifest.json"]
	with cached_mask_artifacts(cache_root / "combined" / f"chr{chrom}", signature, names) as (cache, staging):
		if staging is not None:
			output = staging
			marker = output / "manifest.json"
			if cpg_enabled:
				exclude = cached_cpg_sites(fasta, variants, alignments, chrom, length, cache_root)
			else:
				exclude = np.zeros(length, dtype=bool)
			for lo, hi in bed_intervals(dups, chrom, length, (1, 2, 3)):
				exclude[lo:hi] = True
			with subprocess.Popen(
				["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(modern)],
				stdout=subprocess.PIPE,
				text=True,
			) as proc:
				for line in proc.stdout:
					ch, pos, ref, alt = line.split()
					pos = int(pos)
					if ch.removeprefix("chr") != chrom:
						raise ValueError("Modern VCF must contain only the requested chromosome")
					if len(ref) > 1:
						exclude[max(0, pos - 6) : min(length, pos + 5 + len(ref))] = True
					if len(alt) > 1:
						exclude[max(0, pos - 6) : min(length, pos + 5)] = True
				if proc.wait():
					raise RuntimeError("bcftools failed while building indel mask")
			accessible = np.ones(length, dtype=bool) if strict is None else np.zeros(length, dtype=bool)
			if strict is not None:
				for lo, hi in bed_intervals(strict, chrom, length):
					accessible[lo:hi] = True
			accessible &= ~exclude
			totals = {}
			for ref in refs:
				allowed = np.zeros(length, dtype=bool)
				for lo, hi in bed_intervals(minimal.get(ref, extra.get(ref)), chrom, length):
					allowed[lo:hi] = True
				allowed &= accessible
				totals[ref] = int(allowed.sum())
				if not totals[ref]:
					raise ValueError(f"No callable bases after masking: {ref} chr{chrom}")
				write_boolean_bed(output / f"{ref}.exclude.bed", "23" if chrom == "X" else chrom, ~allowed)
			components = [
				"reference-specific published quality/callability masks",
				"segmental duplications",
				"modern indels +/-5bp",
			]
			if strict is not None:
				components.append("1KG strict accessibility")
			if cpg_enabled:
				components.append("CpG: hg19 + upstream modern variants + panTro2/ponAbe2/rheMac2")
			marker.write_text(
				json.dumps(
					dict(
						signature=signature,
						profile=(
							"multi_reference"
							if extra
							else "cell2020_chrX_extension" if chrom == "X" else "cell2020"
						),
						reference_callability={
							ref: (
								"published_minimal_mask"
								if ref in minimal
								else "published_individual_FilterBed"
							)
							for ref in refs
						},
						callable_bp=totals,
						cpg_filter="applied" if cpg_enabled else "skipped_missing_chrX_axt",
						strict_accessibility="applied" if strict is not None else "unavailable_for_chrX",
						mask_semantics="excluded; BED0; native X contig=23",
						components=components,
					),
					indent=2,
				)
			)
	link_cached_masks(cache, analysis_output, names)
	print(f"IBDMIX chr{chrom}: analysis masks {cache}", flush=True)


def ibdmix_workflow_union(intervals):
	result = []
	for lo, hi in sorted(intervals):
		if result and lo <= result[-1][1]:
			result[-1] = (result[-1][0], max(hi, result[-1][1]))
		else:
			result.append((lo, hi))
	return result


def ibdmix_workflow_subtract(lo, hi, mask):
	import bisect

	i = max(0, bisect.bisect_right(mask, (lo, float("inf"))) - 1)
	for j in range(i, len(mask)):
		left, right = mask[j]
		if left >= hi:
			break
		if right <= lo:
			continue
		if left > lo:
			yield lo, left
		lo = max(lo, right)
		if lo >= hi:
			return
	if lo < hi:
		yield lo, hi


FIELDS = [
	"ID",
	"chrom",
	"start",
	"end",
	"length",
	"slod",
	"sites",
	"positive_lods",
	"negative_lods",
	"sample_set",
	"super_pop",
	"anc",
	"locus_id",
	"genome_build",
	"parent_start",
	"parent_end",
	"score_scope",
]


def finalize(
	raw_dir,
	populations,
	refs,
	output,
	chrom,
	build,
	locus,
	core_start=0,
	core_end=0,
	min_bp=50000,
	lod=4,
	background=True,
	export_denisovan=False,
):
	with open(populations) as handle:
		groups = list(csv.DictReader(handle, delimiter="\t"))
	calls = []
	control = []
	control2013 = []
	for ref in refs:
		for pop in groups:
			if (
				background
				and not export_denisovan
				and ref in DENISOVANS
				and pop["population"] not in AFRICAN_CONTROLS
			):
				continue
			path = Path(raw_dir) / f"{ref}.{pop['population']}.raw.txt.gz"
			members = set(Path(pop["sample_file"]).read_text().split())
			with opener(path) as handle:
				for row in csv.DictReader(handle, delimiter="\t"):
					if row["ID"] not in members or row["chrom"].removeprefix("chr") not in {
						chrom,
						"23" if chrom == "X" else chrom,
					}:
						raise ValueError(f"Wrong population sample or chromosome in {path}")
					start = int(row["start"]) - 1
					end = int(row["end"]) - 1
					score = float(row.get("LOD", row.get("slod", "nan")))
					if start < 0 or end <= start:
						raise ValueError(f"Invalid native IBDmix interval: {path}")
					if not np.isfinite(score):
						raise ValueError(f"Invalid LOD score in {path}")
					if score < lod or end - start < min_bp:
						continue
					if ref == "Denisova" and pop["population"] in AFRICAN_CONTROLS:
						control.append((start, end))
					if ref == "Denisova.2013" and pop["population"] in AFRICAN_CONTROLS:
						control2013.append((start, end))
					if ref in NEANDERTHALS or (ref in DENISOVANS and (not background or export_denisovan)):
						calls.append((ref, pop, row, start, end, score))
	mask = ibdmix_workflow_union(control) if background else []
	mask2013 = ibdmix_workflow_union(control2013) if background else []
	Path(output).parent.mkdir(parents=True, exist_ok=True)
	with gzip.open(output, "wt") as handle:
		writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
		writer.writeheader()
		for ref, pop, row, start, end, score in calls:
			# The African Denisova control mask belongs to the Neanderthal
			# protocol. Applying it to Denisovan calls would erase the signal
			# from the controls that generated the mask.
			for lo, hi in ibdmix_workflow_subtract(
				start, end, (mask2013 if ref == "Altai.2013" else mask) if ref in NEANDERTHALS else []
			):
				if hi - lo < min_bp:
					continue  # authors reapply length filter after subtraction
				if core_end:
					lo = max(lo, core_start)
					hi = min(hi, core_end)
				if hi <= lo:
					continue
				writer.writerow(
					dict(
						ID=row["ID"],
						chrom=chrom,
						start=lo,
						end=hi,
						length=hi - lo,
						slod=score,
						sites=row.get("sites", ""),
						positive_lods=row.get("positive_lods", ""),
						negative_lods=row.get("negative_lods", ""),
						sample_set=pop["population"],
						super_pop=pop["super_population"],
						anc=ref,
						locus_id=locus,
						genome_build=build,
						parent_start=start,
						parent_end=end,
						score_scope="parent_native_call",
					)
				)
	Path(str(output) + ".afr_denisovan.bed").write_text("".join(f"{chrom}\t{lo}\t{hi}\n" for lo, hi in mask))


def validate_genotypes(path, output, x_nonpar_build=None):
	rows = ref_rows = 0
	last = None
	stream = gzip.open if str(path).endswith((".gz", ".bgz")) else open
	with stream(path, "rb") as handle:
		header = handle.readline().strip().split(b"\t")
		expected = len(header) - 4
		if expected < 2 or len(set(header[4:])) != expected:
			raise ValueError("Invalid genotype sample header")
		for line in handle:
			fields = line.split(b"\t", 4)
			if len(fields) != 5:
				raise ValueError("Malformed genotype row")
			body = fields[4].strip()
			key = (fields[0], int(fields[1]))
			if x_nonpar_build and (
				fields[0].lower().removeprefix(b"chr") not in (b"x", b"23")
				or not any(lo <= key[1] - 1 < hi for lo, hi in x_nonpar_intervals(x_nonpar_build))
			):
				raise ValueError(f"Scoring position outside X non-PAR: {key}")
			if body.count(b"\t") != expected - 1 or body.translate(None, b"0129\t"):
				raise ValueError(f"Unsupported genotype encoding at {key}; expected diploid 0/1/2/9 values")
			if last is not None and (key[0] != last[0] or key[1] <= last[1]):
				raise ValueError("Genotype positions must be strictly increasing on one chromosome")
			last = key
			rows += 1
			ref_rows += body.startswith(b"0\t")
	if not rows:
		raise ValueError("No informative genotype records")
	Path(output).write_text(
		json.dumps(
			dict(
				rows=rows,
				archaic_hom_ref_rows=ref_rows,
				modern_samples=expected - 1,
				scoring_build=x_nonpar_build,
				allowed_intervals_bed0=x_nonpar_intervals(x_nonpar_build) if x_nonpar_build else None,
			),
			indent=2,
		)
	)


def split_x_genotypes(path, output, build):
	"""Each connected non-PAR component gets an independent caller state."""
	from contextlib import ExitStack

	spans = x_nonpar_intervals(build)
	with gzip.open(path, "rt") as source, ExitStack() as stack:
		header = source.readline()
		outputs = [
			stack.enter_context(gzip.open(str(output) + f".{i+1}.gz", "wt")) for i in range(len(spans))
		]
		for handle in outputs:
			handle.write(header)
		for line in source:
			fields = line.split("\t", 2)
			pos = int(fields[1]) - 1
			part = next((i for i, (lo, hi) in enumerate(spans) if lo <= pos < hi), None)
			if fields[0].lower().removeprefix("chr") not in ("x", "23") or part is None:
				raise ValueError("PAR or non-X position in X genotype table")
			outputs[part].write(line)


def ibdmix_workflow_main():
	parser = argparse.ArgumentParser()
	subs = parser.add_subparsers(dest="action", required=True)
	p = subs.add_parser("samples")
	for flag in ["panel", "actual", "output"]:
		p.add_argument("--" + flag, required=True)
	p.add_argument("--no-background", action="store_true")
	p = subs.add_parser("mask")
	for flag in ["root", "chrom", "modern", "fasta", "upstream", "output"]:
		p.add_argument("--" + flag, required=True)
	p.add_argument("--axt-root")
	p.add_argument("--refs", nargs="+")
	p.add_argument("--archaic-root")
	p = subs.add_parser("check-masks")
	p.add_argument("--root", required=True)
	p.add_argument("--chroms", nargs="+", required=True)
	p.add_argument("--refs", nargs="+", required=True)
	p.add_argument("--archaic-root")
	p.add_argument("--custom-masks")
	p = subs.add_parser("finalize")
	for flag in ["raw-dir", "populations", "output", "chrom", "build", "locus"]:
		p.add_argument("--" + flag, required=True)
	p.add_argument("--refs", nargs="+", required=True)
	p.add_argument("--core-start", type=int, default=0)
	p.add_argument("--core-end", type=int, default=0)
	p.add_argument("--min-bp", type=int, default=50000)
	p.add_argument("--lod", type=float, default=4)
	p.add_argument("--no-background", action="store_true")
	p.add_argument("--export-denisovan", action="store_true")
	p = subs.add_parser("validate-mask")
	for flag in ["path", "chrom", "build"]:
		p.add_argument("--" + flag, required=True)
	p = subs.add_parser("split-x-genotypes")
	p.add_argument("--path", required=True)
	p.add_argument("--output", required=True)
	p.add_argument("--build", required=True)
	p = subs.add_parser("validate-genotypes")
	p.add_argument("--path", required=True)
	p.add_argument("--output", required=True)
	p.add_argument("--x-nonpar-build", choices=["37", "38"])
	args = vars(parser.parse_args())
	action = args.pop("action")
	if action == "samples":
		args["background_required"] = not args.pop("no_background")
		samples(**args)
	elif action == "mask":
		prepare_masks(**args)
	elif action == "check-masks":
		try:
			count = check_mask_inputs(**args)
		except (FileNotFoundError, ValueError, KeyError) as exc:
			parser.exit(1, f"ERROR: {exc}\n")
		print(f"IBDMIX mask preflight passed: {count} local resources; automatic mask downloads disabled")
	elif action == "validate-mask":
		length = CHROM_LENGTHS[args["build"].removeprefix("GRCh").removeprefix("b")][args["chrom"]]
		previous = 0
		for lo, hi in bed_intervals(args["path"], args["chrom"], length):
			if lo < previous:
				raise ValueError("Excluded mask must be sorted and merged")
			previous = hi
		# Reject wrong-contig rows rather than silently interpreting an include
		# mask or a mismatched chromosome as an empty excluded mask.
		with opener(args["path"]) as handle:
			for line in handle:
				if (
					line.strip()
					and not line.startswith(("#", "track", "browser"))
					and line.split()[0] != args["chrom"]
				):
					raise ValueError(
						"Use a chromosome-specific BED with numeric contig names (X for custom X)"
					)
	elif action == "split-x-genotypes":
		split_x_genotypes(**args)
	elif action == "validate-genotypes":
		validate_genotypes(**args)
	else:
		args["background"] = not args.pop("no_background")
		finalize(**args)


def ibdmix_workflow_cli():
	import sys

	if len(sys.argv) == 3 and sys.argv[1] == "normalize-refs":
		try:
			refs = normalize_refs(sys.argv[2])
			if any(ref not in CALL_REFERENCES for ref in refs):
				raise ValueError(
					"IBDmix calling references are limited to "
					+ ",".join(CALL_REFERENCES)
					+ "; .2013 replication references are retired (historical results remain readable)"
				)
			print(" ".join(refs))
		except ValueError as error:
			print(str(error), file=sys.stderr)
			sys.exit(2)
	else:
		ibdmix_workflow_main()


WORKFLOW_SUBCOMMANDS = {
	"provenance": ibdmix_provenance_cli,
}


def workflow_main():
	import sys

	if len(sys.argv) > 1 and sys.argv[1] in WORKFLOW_SUBCOMMANDS:
		command = sys.argv.pop(1)
		return WORKFLOW_SUBCOMMANDS[command]()
	return ibdmix_workflow_cli()


import io
import shlex
import tarfile
import zlib
from functools import lru_cache

sys.dont_write_bytecode = True
import ukb_prepare as ukb

resume_F = Path(__file__).resolve().parent


def resume_require(condition, message):
	if not condition:
		raise ValueError(message)


def resume_digest(data):
	return hashlib.sha256(data).hexdigest()


@lru_cache(maxsize=128)
def resume_file_digest(path, size, mtime_ns):
	with open(path, "rb") as handle:
		return hashlib.file_digest(handle, "sha256").hexdigest()


def resume_sha(path):
	path = Path(path)
	info = path.stat()
	return resume_file_digest(str(path), info.st_size, info.st_mtime_ns)


def resume_atomic_json(path, value):
	path.parent.mkdir(parents=True, exist_ok=True)
	fd, name = tempfile.mkstemp(prefix=".resume-", dir=path.parent)
	try:
		with os.fdopen(fd, "w") as handle:
			json.dump(value, handle)
		os.replace(name, path)
	finally:
		Path(name).unlink(missing_ok=True)


def resume_wanted(name):
	return (
		name in {".complete", "run.meta.tsv", "cache.meta.tsv", "inputs/ukb.samples.tsv"}
		or (name.startswith("inputs/ukb-scoring/") and name.endswith("/target.json"))
		or (name.startswith("qc/masks/") and name.endswith(".preparation.json"))
	)


def resume_evidence(run):
	"""Freeze small receipts before stage-qc can overwrite historical inputs."""
	record = (run / "run.meta.tsv").read_text()
	key = resume_digest(record.encode())
	receipt = run / ".ibdmix.resume.json"
	archive = run / "ibdmix.raw.tar.gz"
	archive_sig = ukb.signature([archive]) if archive.is_file() else []
	completion = (run / ".complete").exists()
	if receipt.is_file():
		saved = json.loads(receipt.read_text())
		if (
			saved.get("schema") == 1
			and saved.get("record_sha256") == key
			and saved.get("archive") == archive_sig
			and (archive_sig or saved.get("completion") == completion)
		):
			return saved
	texts, outputs = {}, {}
	if archive.is_file():
		observed = {}
		manifest = None
		with tarfile.open(archive, "r|gz") as bundle:
			for member in bundle:
				name = member.name
				if not (resume_wanted(name) or name.startswith("final/") or name == "GU-MANIFEST.json"):
					continue
				resume_require(member.isfile(), f"nonregular archived receipt: {name}")
				with bundle.extractfile(member) as handle:
					data = handle.read()
				observed[name] = resume_digest(data)
				if resume_wanted(name):
					texts[name] = data.decode()
				elif name == "GU-MANIFEST.json":
					manifest = json.loads(data)
				elif name.endswith(".segments.tsv.gz"):
					with gzip.GzipFile(fileobj=io.BytesIO(data)) as handle:
						while handle.read(1024 * 1024):
							pass
					outputs[name] = observed[name]
		resume_require(manifest is not None, "archive lacks GU-MANIFEST.json")
		declared = {row["path"]: row for row in manifest["files"]}
		for name, checksum in observed.items():
			if name != "GU-MANIFEST.json":
				resume_require(
					declared.get(name, {}).get("sha256") == checksum, f"archive checksum mismatch: {name}"
				)
		resume_require(texts.get("run.meta.tsv") == record, "published and working run contracts differ")
	else:
		for path in run.rglob("*"):
			name = path.relative_to(run).as_posix()
			if path.is_file() and resume_wanted(name):
				texts[name] = path.read_text()
		if completion:
			for path in (run / "final").glob("*.segments.tsv.gz"):
				with gzip.open(path, "rb") as handle:
					while handle.read(1024 * 1024):
						pass
				outputs[path.relative_to(run).as_posix()] = resume_sha(path)
	saved = dict(
		schema=1, record_sha256=key, archive=archive_sig, completion=completion, texts=texts, outputs=outputs
	)
	resume_atomic_json(receipt, saved)
	return saved


def resume_command(path):
	env, argv = dict(os.environ), []
	for line in path.read_text().splitlines():
		words = shlex.split(line)
		if words[:1] == ["export"]:
			env.update(word.split("=", 1) for word in words[1:])
		elif words[:1] == ["exec"]:
			argv = words[1:]
	resume_require(len(argv) >= 2 and argv[1] in {"ibdmix", "ibdmix-GPU"}, "not an IBDmix worker command")
	env["IBDMIX_BACKEND"] = "cuda" if argv[1] == "ibdmix-GPU" else "cpu"
	args = argv[2:]
	if args[:1] == ["run"]:
		args = args[1:]
	if "--keep-males" in args:
		args.insert(args.index("--keep-males") + 1, "TRUE")
	resume_require(
		len(args) % 2 == 0 and all(x.startswith("--") for x in args[::2]), "unsupported worker arguments"
	)
	return env, dict(zip(args[::2], args[1::2]))


def resume_records(text):
	return [line.replace("\\t", "\t").split("\t") for line in text.splitlines() if line]


def resume_scoring(saved, chrom):
	return json.loads(saved["texts"][f"inputs/ukb-scoring/chr{chrom}/target.json"])


def resume_disposable(path):
	p = Path(path)
	return (p.parent.parent == Path("/tmp/gu-psam-chunks") and p.name == "chunk.keep") or (
		p.parent.parent == Path("/tmp/gu-cohorts") and p.name in {"samples.txt", "cohort.keep"}
	)


def resume_stable_scoring(value):
	"""Only generated selection files have a content-based replacement below."""
	value = json.loads(json.dumps(value))
	value["preparation"]["inputs"] = [
		x for x in value["preparation"]["inputs"] if not resume_disposable(x["path"])
	]
	return value


def resume_verify_stats(items):
	for item in items:
		resume_require(ukb.signature([item["path"]])[0] == item, f'input changed: {item["path"]}')


@lru_cache(maxsize=32)
def resume_psam_rows(path, size, mtime_ns):
	return ukb.read_psam(path)


def resume_samples_at(path):
	path = Path(path).resolve()
	info = path.stat()
	# attach_panel mutates rows; cached PSAM data must stay pristine.
	return [dict(row) for row in resume_psam_rows(str(path), info.st_size, info.st_mtime_ns)]


def resume_verify_cohort(prep, env, opts):
	chrom = opts["--chr"]
	root = Path(opts["--target-dir"]).parent.resolve()
	resume_require(prep["pgen_root"] == str(root), "target PGEN directory changed")
	resume_require(Path(opts["--target-dir"]).name == "chr", "unsupported UKB PGEN prefix")
	options = prep["options"]
	male = opts.get("--keep-males") == "TRUE"
	resume_require(
		options.get("existing_pgen") is True and options["male_only"] == male,
		"PGEN mode or sex selection changed",
	)
	resume_require(
		options["grch"] == int(opts["--grch"]) and options["ukb_source"] == env["GU_UKB_SOURCE"],
		"build or UKB source changed",
	)
	resume_require(
		options["ukb_require_phase"] == ukb.boolean(env.get("UKB_REQUIRE_PHASE", "false")),
		"phase requirement changed",
	)
	resume_require(
		not env.get("UKB_INFO_FILE") and not env.get("UKB_QC_PASS_VARIANTS"),
		"external UKB variant filters require full validation",
	)
	metadata = root / ("chr1.psam" if chrom == "X" else f"chr{chrom}.psam")
	samples = resume_samples_at(metadata)
	if "--keep-psam" in opts:
		name, size, index = opts["--keep-psam"].rsplit(",", 2)
		path, size, index = Path(name).resolve(), int(size), int(index)
		resume_require(size > 0 and index > 0, "invalid PSAM chunk")
		info = options.get("keep_psam", {})
		resume_require(
			(str(path), size, index, resume_sha(path))
			== (
				info.get("source"),
				info.get("chunk_size"),
				info.get("chunk_index"),
				info.get("source_sha256"),
			),
			"PSAM chunk selection changed",
		)
		keep = [s["sample"] for s in resume_samples_at(path)[(index - 1) * size : index * size]]
	elif opts.get("--keep") or env.get("GU_KEEP"):
		keep = ukb.read_ids(opts.get("--keep") or env["GU_KEEP"])
		resume_require(not options.get("keep_psam"), "PSAM chunk selection changed")
	else:
		raise ValueError("missing fixed-cohort selection")
	panel = opts.get("--sample-panel")
	samples = ukb.selected_samples(ukb.attach_panel(samples, panel), keep, "1")
	if male:
		samples = [s for s in samples if s["sex"] == "male"]
	ukb.validate_cohort(samples, argparse.Namespace(sample_panel=panel))
	resume_require(
		ukb.cohort_hash(samples) == prep["cohort_metadata_sha256"],
		"sample cohort/sex/population metadata changed",
	)
	_, paths = ukb.pfile(root, chrom)
	full = ukb.attach_panel(resume_samples_at(paths[1]), panel)
	keep = [s["sample"] for s in samples if chrom != "X" or s["sex"] == "male"]
	if chrom == "X":
		present = {s["sample"] for s in full}
		keep = [s for s in keep if s in present]
	selected = ukb.selected_samples(full, keep, chrom)
	resume_require(
		ukb.cohort_hash(selected) == prep["scored_sample_metadata_sha256"]
		and ukb.sample_hash(selected) == prep["sample_order_sha256"],
		"chromosome sample set/order changed",
	)
	recorded = {x["path"] for x in prep["inputs"]}
	resume_require(all(str(p.resolve()) in recorded for p in paths), "PGEN input paths changed")
	fasta = Path(
		env.get(
			"UKB_REF_FASTA",
			str(
				Path(env.get("GU_REF_ROOT", "/mnt/f/gen"))
				/ f'1hgp/{opts["--grch"]}/GRCH{opts["--grch"]}.fasta.gz'
			),
		)
	).resolve()
	resume_require(
		str(fasta) in recorded and str(fasta) + ".fai" in recorded, "preparation reference FASTA changed"
	)
	resume_verify_stats(x for x in prep["inputs"] if not resume_disposable(x["path"]))


# The original file SHA remains immutable audit evidence. This exact migration
# accepts its masks only if the AST of every mask routine and constant is unchanged.
# New caches still record the current file SHA; unknown or changed code is rejected.
MASK_IMPLEMENTATION_FUNCTIONS = (
	"opener",
	"fingerprint",
	"file_checksum",
	"sha256_file",
	"bed_intervals",
	"write_boolean_bed",
	"cpg_sites",
	"axt_directory",
	"cached_mask_artifacts",
	"cached_cpg_sites",
	"link_cached_masks",
	"axt_reference",
	"local_axt",
	"published_mask_name",
	"local_mask_candidates",
	"local_mask_ready",
	"require_local_mask",
	"local_reference",
	"published_reference_mask",
	"reference_mask_resource",
	"check_mask_inputs",
	"prepare_masks",
)
MASK_IMPLEMENTATION_CONSTANTS = ("ibdmix_workflow_VERSION", "MINIMAL_MASK_NAMES")
MASK_PRECONSOLIDATION_FILE_SHA256 = "800b1b29ac695c197724bf12b83148eaa55131ef75185d02f904386babfb2ab4"
MASK_PRECONSOLIDATION_IMPLEMENTATION = "7cf342f562afd3a7bc5e2f0d92a4281e8b16e1f520b8a579831a44630e6d0005"


def mask_implementation_digest(source):
	"""Fingerprint mask calculations independently of file layout and formatting."""
	import ast

	class StripNonComputational(ast.NodeTransformer):
		def visit_Import(self, node):
			return None

		def visit_ImportFrom(self, node):
			return None

		def visit_Expr(self, node):
			if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
				return None
			return self.generic_visit(node)

	def structure(value):
		if isinstance(value, ast.AST):
			return [
				type(value).__name__,
				{
					name: structure(item)
					for name, item in ast.iter_fields(value)
					if name != "type_params" or item
				},
			]
		if isinstance(value, list):
			return [structure(item) for item in value]
		if isinstance(value, bytes):
			return {"bytes": value.hex()}
		if value is Ellipsis:
			return {"ellipsis": True}
		return value

	nodes = {}
	for node in ast.parse(source).body:
		if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in MASK_IMPLEMENTATION_FUNCTIONS:
			nodes[node.name] = structure(StripNonComputational().visit(node))
		elif isinstance(node, ast.Assign):
			for target in node.targets:
				if isinstance(target, ast.Name) and target.id in MASK_IMPLEMENTATION_CONSTANTS:
					nodes[target.id] = structure(node.value)
	if set(nodes) != set(MASK_IMPLEMENTATION_FUNCTIONS) | set(MASK_IMPLEMENTATION_CONSTANTS):
		return None
	return hashlib.sha256(json.dumps(nodes, sort_keys=True).encode()).hexdigest()


def mask_code_matches(recorded):
	"""Accept the audited pre-consolidation masks only while their code is identical."""
	if recorded == sha256_file(__file__):
		return True
	return (
		recorded == MASK_PRECONSOLIDATION_FILE_SHA256
		and mask_implementation_digest(Path(__file__).read_text()) == MASK_PRECONSOLIDATION_IMPLEMENTATION
	)


def resume_verify_request(saved, env, opts, complete):
	chrom = opts["--chr"]
	old = resume_records(saved["texts"]["run.meta.tsv"])
	contract(saved["texts"]["run.meta.tsv"])
	values = {row[0]: row[1] for row in old if len(row) == 2}
	cuda = env["IBDMIX_BACKEND"] == "cuda"
	version = re.search(r"^pipeline_version=([^\n]+)", (resume_F / "ibdmix.sh").read_text(), re.M)[1]
	if cuda:
		version += ".cuda-segment-v1"
	expected = dict(
		pipeline_version=version,
		profile=env.get("IBDMIX_PROFILE", "multi_reference"),
		refs=env["IBDMIX_REFS"],
		genome_build="b" + opts["--grch"],
		background_filter=env.get("IBDMIX_AFR_DENISOVAN_FILTER", "0"),
		custom_masks=env.get("IBDMIX_MASK_DIR", ""),
		lod_cut=env.get("IBDMIX_LOD", "4"),
		len_cut=env.get("IBDMIX_MIN_BP", "50000"),
		emit_lod_cut=env.get("emit_lod_cut", env.get("IBDMIX_LOD", "4")),
		minor_allele_count=env.get("minor_allele_count", "1"),
		archaic_error=env.get("archaic_error", "0.01"),
		modern_error_max=env.get("modern_error_max", "0.002"),
		modern_error_proportion=env.get("modern_error_proportion", "2"),
	)
	if chrom == "X":
		expected.update(
			x_scoring_profile=env.get("IBDMIX_X_PROFILE", "nonpar-v2"),
			x_minor_allele_count=env.get("IBDMIX_X_MINOR_ALLELE_COUNT", "1"),
		)
	for key, value in expected.items():
		resume_require(
			values.get(key) == value, f"parameter changed: {key} ({values.get(key)!r} -> {value!r})"
		)
	soft = Path(env.get("GU_SOFT", "/mnt/d/software/gu"))
	prepare = Path(env.get("IBDMIX_PREPARE_BIN", soft / "ibdmix-prepare/prepare"))
	resume_require(
		resume_sha(prepare) == values["genotype_binary_sha256"], "genotype preparation binary changed"
	)
	if cuda:
		resume_require(
			resume_sha(env.get("IBDMIX_GPU_BIN", soft / "ibdmix-gpu/ibdmix-gpu"))
			== values["gpu_binary_sha256"],
			"GPU caller binary changed",
		)
	else:
		runtime = Path(env.get("IBDMIX_RUNTIME", soft / "ibdmix"))
		if not runtime.is_dir():
			runtime = soft / "IBDmix"
		caller = runtime / "build/src/ibdmix"
		info = caller.stat()
		resume_require(
			["software", f"{caller}:{info.st_size}:{int(info.st_mtime)}"] in old, "CPU caller binary changed"
		)
	value = resume_scoring(saved, chrom)
	recorded_contract = value["contract"]
	resume_verify_cohort(recorded_contract["preparation"], env, opts)
	resume_require(
		recorded_contract["references"] == env["IBDMIX_REFS"].split(), "archaic reference selection changed"
	)
	source_root = Path(opts["--archaic-path"]).resolve()
	sources = [
		ukb._ukb_scoring_source(source_root, ref, chrom).resolve() for ref in recorded_contract["references"]
	]
	resume_require(
		[x["path"] for x in recorded_contract["inputs"] if Path(x["path"]) != resume_F / "ukb_prepare.py"]
		== list(map(str, sources)),
		"archaic source paths changed",
	)
	resume_verify_stats(recorded_contract["inputs"])
	mask_root = str(Path(env.get("IBDMIX_MASK_ROOT", source_root.parent / "mask")).resolve())
	resume_require(recorded_contract["original_mask_root"] == mask_root, "mask root changed")
	fallback = source_root.parent / "mask"
	resume_require(
		recorded_contract.get("original_fallback_mask_root") in {None, str(fallback)},
		"fallback mask root changed",
	)
	fasta = env.get(
		"IBDMIX_REFERENCE_FASTA", str(Path(env.get("GU_REF_ROOT", "/mnt/f/gen")) / "1hgp/37/GRCH37.fasta.gz")
	)
	resume_require(values["reference_fasta"].rsplit(":", 2)[0] == fasta, "mask reference FASTA changed")
	# The scoring receipt must belong to these calls, including after a failed
	# rerun has overwritten the working inputs directory.
	prep_hash = resume_digest(json.dumps(recorded_contract["preparation"], sort_keys=True).encode())
	resume_require(
		["modern_source", "chr" + chrom, "preparation_contract_sha256", prep_hash] in old,
		"scoring receipt does not match the original call contract",
	)
	if complete:
		masks = [json.loads(text) for name, text in saved["texts"].items() if name.startswith("qc/masks/")]
		resume_require(masks, "completed run lacks mask provenance")
		for mask in masks:
			sig = mask["signature"]
			resume_require(mask_code_matches(sig["code"]), "mask implementation changed")
			resume_require(sig["refs"] == recorded_contract["references"], "mask references changed")
			resume_verify_stats(dict(path=p, bytes=n, mtime_ns=t) for p, n, t in sig["inputs"])
			if chrom == "X":
				enabled = all(
					local_mask_ready(axt_reference(chrom, species, None)[1])
					for species in ["panTro2", "ponAbe2", "rheMac2"]
				)
				resume_require(sig["cpg_enabled"] == enabled, "chrX CpG resources changed")


def resume_check(cmd):
	env, opts = resume_command(cmd)
	if (
		env.get("GU_UKB_DIRECT") != "1"
		or "--chr" not in opts
		or opts.get("--replace-ibdmix", "FALSE") != "FALSE"
	):
		return 1
	run = Path(env["IBDMIX_OUT"])
	if not (run / "run.meta.tsv").is_file():
		return 1
	saved = resume_evidence(run)
	complete = ".complete" in saved["texts"]
	resume_verify_request(saved, env, opts, complete)
	if not complete:
		print(f'[GU CMD] RESUME unit=chr{opts["--chr"]} reason=partial_calls_verified', file=sys.stderr)
		return 1
	final = f'final/all_archaic_refs.lod{env.get("IBDMIX_LOD", "4")}.len{env.get("IBDMIX_MIN_BP", "50000")}.segments.tsv.gz'
	resume_require(final in saved["outputs"], "completed run lacks verified final segments")
	if (run / final).exists():
		resume_require(
			resume_sha(run / final) == saved["outputs"][final], "final segments changed or corrupted"
		)
	else:
		resume_require(bool(saved["archive"]), "final segments are missing")
	print(run / final if (run / final).exists() else run / "ibdmix.raw.tar.gz")
	return 0


def resume_provenance(previous, current, run):
	before, after = previous.read_text(), current.read_text()
	if not differences(before, after):
		return 0
	receipt = run / ".ibdmix.resume.json"
	if receipt.is_file():
		saved = json.loads(receipt.read_text())
		resume_require(
			saved["texts"]["run.meta.tsv"] == (run / "run.meta.tsv").read_text(),
			"resume receipt no longer matches original run",
		)
		old_paths = [
			n for n in saved["texts"] if n.startswith("inputs/ukb-scoring/") and n.endswith("/target.json")
		]
		resume_require(len(old_paths) == 1, "expected one chromosome scoring receipt")
		old = json.loads(saved["texts"][old_paths[0]])["contract"]
		new = json.loads((run / old_paths[0]).read_text())["contract"]
		resume_require(
			resume_stable_scoring(old) == resume_stable_scoring(new),
			"original sources, cohort, preparation options or scored sites changed",
		)
		for text, value in [(before, old), (after, new)]:
			prep_hash = resume_digest(json.dumps(value["preparation"], sort_keys=True).encode())
			resume_require(
				["modern_source", "chr" + value["chromosome"], "preparation_contract_sha256", prep_hash]
				in resume_records(text),
				"scoring receipt does not match the compared call contract",
			)

		# Only exchange-file identities are disposable. Masks, caller binaries,
		# thresholds and all other provenance still use the existing comparator.
		def normalize(text):
			return (
				"\n".join(
					"\t".join(row)
					for row in resume_records(text)
					if row[0] not in {"sample_file", "archaic_vcf"}
					and not (row[0] == "modern_source" and row[2] == "preparation_contract_sha256")
				)
				+ "\n"
			)

		diff = differences(normalize(before), normalize(after))
	else:
		diff = differences(before, after)
	if diff:
		print(diff, end="", file=sys.stderr)
		return 1
	return 0


def resume_main():
	parser = argparse.ArgumentParser(description=__doc__)
	sub = parser.add_subparsers(dest="action", required=True)
	p = sub.add_parser("check")
	p.add_argument("command", type=Path)
	p = sub.add_parser("provenance")
	p.add_argument("previous", type=Path)
	p.add_argument("current", type=Path)
	p.add_argument("--run", type=Path, required=True)
	args = parser.parse_args()
	try:
		return (
			resume_check(args.command)
			if args.action == "check"
			else resume_provenance(args.previous, args.current, args.run)
		)
	except (
		ValueError,
		OSError,
		KeyError,
		IndexError,
		TypeError,
		EOFError,
		tarfile.TarError,
		zlib.error,
	) as error:
		print(f"ERROR: IBDmix resume check: {error}", file=sys.stderr)
		return 2


import bisect
import itertools
import shutil
import sqlite3


loci_VERSION = 1
loci_LINEAGES = ("Neanderthal", "Denisovan")
loci_SQL = """SELECT source_class,chr,sample_id,start,end,source,score
    FROM segments_evidence
    WHERE dataset_id=? AND genome_build=? AND method='ibdmix'
      AND reference_role='scientific_reference'
      AND source_class IN ('Neanderthal','Denisovan')
      AND sample_id IS NOT NULL AND sample_id NOT IN ('','nan','None')
    ORDER BY source_class,chr,sample_id,start,end"""
loci_FIELDS = (
	"rank",
	"segment_code",
	"chr",
	"start",
	"end",
	"length_bp",
	"source_class",
	"genome_build",
	"n_tested",
	"peak_start",
	"peak_end",
	"peak_fraction",
	"peak_raw_fraction",
	"peak_effective_copies",
	"mean_fraction",
	"background_fraction",
	"extension_threshold",
	"core_threshold",
	"core_windows",
	"n_windows",
	"excess_support_kb",
	"n_carriers",
	"n_references",
	"reference_names",
	"max_lod",
	"n_calls",
)


def loci_canonical_ref(ref):
	return "Chagyr" if ref == "Chagyrskaya" else ref


def loci_tested_rosters(con, database, dataset, build):
	from results_normalize import density_run_path, reference_role, source_class, whole_chromosome_run

	"""Certify full-autosome denominators from saved caller rosters only.

    Missing provenance is not inferred from positive calls or annotation panels.
    Historical controls and locus-only runs cannot certify these windows.
    """
	runs = con.execute(
		"""SELECT chr,raw_file FROM method_runs
        WHERE dataset_id=? AND genome_build=? AND method='ibdmix'
        AND status='complete' AND evidence_eligible=1""",
		(dataset, build),
	).fetchall()
	paths = {raw: density_run_path(raw, database) for _, raw in runs}
	targets, warnings = {}, []
	for chrom, raw in runs:
		if chrom not in {str(i) for i in range(1, 23)}:
			continue
		path = paths[raw]
		if not path.is_file():
			warnings.append(f"chr{chrom}: missing run metadata; excluded")
			continue
		text = path.read_text().replace("\\t", "\t")
		if not whole_chromosome_run(text, chrom):
			warnings.append(f"chr{chrom}: locus-only/unverified genome scope; excluded")
			continue
		records = dict(line.split("\t", 1) for line in text.splitlines() if "\t" in line)
		roster = path.parent / "samples" / ("C" + chrom) / "ALL.txt"
		if not roster.is_file():
			roster = path.parent / "samples" / "ALL.txt"
		members = set(roster.read_text().split()) if roster.is_file() else set()
		if not members:
			warnings.append(f"chr{chrom}: missing/empty actual caller roster; excluded")
			continue
		background_only = records.get("background_filter") == "1" and records.get("export_denisovan") != "1"
		for ref in records.get("refs", "").split():
			ref = loci_canonical_ref(ref)
			if reference_role(ref) != "scientific_reference":
				continue
			lineage = source_class(ref)
			if lineage not in loci_LINEAGES or (lineage == "Denisovan" and background_only):
				continue
			targets.setdefault((lineage, chrom), {}).setdefault(ref, set()).update(members)
	return targets, paths, sorted(set(warnings))


def loci_union_intervals(intervals):
	"""Merge sorted per-person intervals across all scientific references."""
	left = right = None
	for lo, hi in intervals:
		if hi <= lo:
			continue
		if left is not None and lo <= right:
			right = max(right, hi)
		else:
			if left is not None:
				yield left, right
			left, right = lo, hi
	if left is not None:
		yield left, right


def loci_add_interval(boundary_bp, full_windows, start, end, window_bp):
	"""Accumulate exact covered bp using O(1) updates per union interval."""
	first, last = start // window_bp, (end - 1) // window_bp
	if first == last:
		boundary_bp[first] += end - start
	else:
		boundary_bp[first] += (first + 1) * window_bp - start
		boundary_bp[last] += end - last * window_bp
		full_windows[first + 1] += 1
		full_windows[last] -= 1


def loci_window_profiles(rows, targets, lengths, window_bp):
	profiles = []
	# Include completely tested chromosomes with no positive calls in background.
	for (lineage, chrom), refs in sorted(targets.items()):
		if chrom not in lengths:
			continue
		length = lengths[chrom]
		n = len(set().union(*refs.values()))
		starts = np.arange(0, length, window_bp, dtype=np.int64)
		profiles.append(
			dict(
				source_class=lineage,
				chr=chrom,
				n_tested=n,
				starts=starts,
				ends=np.minimum(starts + window_bp, length),
				boundary=np.zeros(len(starts)),
				full=np.zeros(len(starts) + 1),
			)
		)
	lookup = {(p["source_class"], p["chr"]): p for p in profiles}
	for (lineage, chrom, sample), group in itertools.groupby(rows, key=lambda r: r[:3]):
		p = lookup.get((lineage, chrom))
		if p is None:
			# groupby advances past this group's rows on the next iteration.
			continue
		refs = targets[(lineage, chrom)]
		intervals = (
			(max(0, r[3]), min(lengths[chrom], r[4]))
			for r in group
			if sample in refs.get(loci_canonical_ref(r[5]), set())
		)
		for lo, hi in loci_union_intervals(intervals):
			loci_add_interval(p["boundary"], p["full"], lo, hi, window_bp)
	for p in profiles:
		widths = p["ends"] - p["starts"]
		p["covered_bp"] = p.pop("boundary") + np.cumsum(p.pop("full")[:-1]) * widths
		p["raw"] = p["covered_bp"] / (2 * p["n_tested"] * widths)
	return profiles


def loci_median_smooth(values, size=3):
	"""No padding/reflection across chromosome edges; use available neighbors."""
	radius = size // 2
	return np.array([np.median(values[max(0, i - radius) : i + radius + 1]) for i in range(len(values))])


def loci_true_runs(mask):
	changes = np.diff(np.r_[False, mask, False].astype(np.int8))
	return zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1))


def loci_call_loci(profiles, dataset, build, window_bp=50_000, floor=0.01, smooth_windows=3):
	backgrounds = {}
	for lineage in loci_LINEAGES:
		group = [p for p in profiles if p["source_class"] == lineage]
		numerator = sum(float(p["covered_bp"].sum()) for p in group)
		denominator = sum(2 * p["n_tested"] * int(p["ends"][-1]) for p in group)
		backgrounds[lineage] = numerator / denominator if denominator else None
	loci = []
	for p in profiles:
		lineage, chrom = p["source_class"], p["chr"]
		baseline = max(floor, backgrounds[lineage])
		core_threshold = 2 * baseline
		smooth = loci_median_smooth(p["raw"], smooth_windows)
		# Raw gate prevents smoothing from manufacturing support in a valley.
		extent = (smooth >= baseline) & (p["raw"] >= baseline)
		core = (smooth >= core_threshold) & (p["raw"] >= core_threshold)
		# A short terminal bin cannot certify a full 50 kb core window.
		core &= (p["ends"] - p["starts"]) == window_bp
		for a, b in loci_true_runs(extent):
			if not any(hi - lo >= 2 for lo, hi in loci_true_runs(core[a:b])):
				continue
			start, end = int(p["starts"][a]), int(p["ends"][b - 1])
			if end - start < 2 * window_bp:
				continue
			peak = a + int(np.argmax(smooth[a:b]))
			widths = p["ends"][a:b] - p["starts"][a:b]
			# Area is in fraction*kb (not overlapping-window counts or a p-value).
			excess = float(np.sum(np.maximum(0, smooth[a:b] - baseline) * widths) / 1000)
			identity = f"{dataset}|{build}|{lineage}|{chrom}|{start}|{end}"
			code = "IBD-" + hashlib.sha256(identity.encode()).hexdigest()[:12].upper()
			loci.append(
				dict(
					segment_code=code,
					chr=chrom,
					start=start,
					end=end,
					length_bp=end - start,
					source_class=lineage,
					genome_build=build,
					n_tested=p["n_tested"],
					peak_start=int(p["starts"][peak]),
					peak_end=int(p["ends"][peak]),
					peak_fraction=float(smooth[peak]),
					peak_raw_fraction=float(p["raw"][peak]),
					peak_effective_copies=float(p["raw"][peak] * 2 * p["n_tested"]),
					mean_fraction=float(np.average(p["raw"][a:b], weights=widths)),
					background_fraction=backgrounds[lineage],
					extension_threshold=baseline,
					core_threshold=core_threshold,
					core_windows=int(core[a:b].sum()),
					n_windows=int(b - a),
					excess_support_kb=excess,
				)
			)
	loci.sort(
		key=lambda r: (
			-r["excess_support_kb"],
			-r["peak_fraction"],
			int(r["chr"]),
			r["start"],
			r["source_class"],
		)
	)
	for rank, row in enumerate(loci, 1):
		row["rank"] = rank
	return loci, backgrounds


def loci_annotate_loci(rows, loci, targets):
	"""One streaming pass: no hundred repeated chromosome range scans."""
	index = {}
	for row in loci:
		row.update(_carriers=set(), _references=set(), max_lod=None, n_calls=0)
		index.setdefault((row["source_class"], row["chr"]), []).append(row)
	index = {key: sorted(group, key=lambda r: r["start"]) for key, group in index.items()}
	bounds = {key: ([r["start"] for r in group], [r["end"] for r in group]) for key, group in index.items()}
	for lineage, chrom, sample, start, end, source, score in rows:
		key = (lineage, chrom)
		ref = loci_canonical_ref(source)
		if key not in index or sample not in targets[key].get(ref, set()):
			continue
		starts, ends = bounds[key]
		for i in range(bisect.bisect_right(ends, start), bisect.bisect_left(starts, end)):
			row = index[key][i]
			row["_carriers"].add(sample)
			row["_references"].add(ref)
			row["n_calls"] += 1
			if score is not None and math.isfinite(score):
				row["max_lod"] = score if row["max_lod"] is None else max(score, row["max_lod"])
	for row in loci:
		row["n_carriers"] = len(row.pop("_carriers"))
		refs = sorted(row.pop("_references"))
		row["n_references"] = len(refs)
		row["reference_names"] = ",".join(refs)


def loci_database_stamp(database):
	stat = database.stat()
	return (str(database.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def loci_prepare(database, cache_dir, dataset, build, window_bp=50_000, floor=0.01, smooth_windows=3):
	from results_normalize import density_run_signature

	if window_bp <= 0 or not 0.01 <= floor <= 0.25 or smooth_windows < 3 or smooth_windows % 2 == 0:
		raise ValueError("Require positive window, floor in [0.01,0.25], and odd smoothing width >=3")
	lengths = CHROM_LENGTHS.get(build.replace("GRCh", ""))
	if lengths is None:
		raise ValueError(f"Unsupported genome build: {build}")
	before = loci_database_stamp(database)
	with sqlite3.connect(f"file:{database.resolve()}?mode=ro", uri=True) as con:
		targets, paths, warnings = loci_tested_rosters(con, database, dataset, build)
	provenance = density_run_signature(paths.values())
	settings = dict(
		window_bp=window_bp,
		smooth_windows=smooth_windows,
		floor_fraction=floor,
		core_multiplier=2,
		min_core_windows=2,
	)
	signature = dict(
		version=loci_VERSION,
		database=before,
		provenance=provenance,
		dataset=dataset,
		build=build,
		settings=settings,
		code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
	)
	key = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()[:24]
	cache_dir.mkdir(parents=True, exist_ok=True)
	dest = cache_dir / key
	if (dest / "loci.json").is_file():
		return dest / "loci.json"
	with tempfile.TemporaryDirectory(prefix="gu-ibdmix-loci-") as work:
		# Indexed SQLite access on /mnt is expensive; stage one immutable copy.
		local = Path(work) / "input.sqlite"
		shutil.copyfile(database, local)
		if loci_database_stamp(database) != before:
			raise RuntimeError("Database changed during snapshot; retry")
		with sqlite3.connect(f"file:{local}?mode=ro", uri=True) as con:
			profiles = loci_window_profiles(
				con.execute(loci_SQL, (dataset, build)), targets, lengths, window_bp
			)
			loci, backgrounds = loci_call_loci(profiles, dataset, build, window_bp, floor, smooth_windows)
			if loci:
				loci_annotate_loci(con.execute(loci_SQL, (dataset, build)), loci, targets)
		if loci_database_stamp(database) != before or density_run_signature(paths.values()) != provenance:
			raise RuntimeError("Database or caller roster changed during aggregation; retry")
		meta = dict(
			signature=signature,
			settings=settings,
			dataset=dataset,
			build=build,
			backgrounds=backgrounds,
			candidate_count=len(loci),
			warnings=warnings,
			tested=[
				dict(source_class=p["source_class"], chr=p["chr"], n_tested=p["n_tested"]) for p in profiles
			],
			definition="Per-person union bp / (window bp * 2 * tested individuals); unphased conservative copy support, not dosage or a calibrated probability",
		)
		stage = Path(tempfile.mkdtemp(prefix=".building-", dir=cache_dir))
		try:
			(stage / "loci.json").write_text(json.dumps(dict(metadata=meta, loci=loci), allow_nan=False))
			with (stage / "candidates.tsv").open("w") as handle:
				writer = csv.DictWriter(handle, fieldnames=loci_FIELDS, delimiter="\t")
				writer.writeheader()
				writer.writerows(loci)
			with (stage / "windows.tsv").open("w") as handle:
				writer = csv.writer(handle, delimiter="\t")
				writer.writerow(
					("source_class", "chr", "start", "end", "n_tested", "raw_fraction", "smooth_fraction")
				)
				for p in profiles:
					smooth = loci_median_smooth(p["raw"], smooth_windows)
					writer.writerows(
						(p["source_class"], p["chr"], int(lo), int(hi), p["n_tested"], float(raw), float(sm))
						for lo, hi, raw, sm in zip(p["starts"], p["ends"], p["raw"], smooth)
					)
			try:
				os.rename(stage, dest)
			except OSError:
				if not (dest / "loci.json").is_file():
					raise
		finally:
			if stage.exists():
				shutil.rmtree(stage)
	return dest / "loci.json"


def loci_main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--database", type=Path, required=True)
	parser.add_argument("--cache-dir", type=Path, default=Path("/tmp/gu-shiny-evidence"))
	parser.add_argument("--dataset", required=True)
	parser.add_argument("--build", required=True)
	parser.add_argument("--window-bp", type=int, default=50_000)
	parser.add_argument("--floor", type=float, default=0.01)
	parser.add_argument("--smooth-windows", type=int, default=3)
	args = parser.parse_args()
	print(
		loci_prepare(
			args.database,
			args.cache_dir,
			args.dataset,
			args.build,
			args.window_bp,
			args.floor,
			args.smooth_windows,
		)
	)


def main():
	if sys.argv[1:2] == ["resume"]:
		sys.argv.pop(1)
		return resume_main()
	if sys.argv[1:2] == ["loci"]:
		sys.argv.pop(1)
		return loci_main()
	return workflow_main()


if __name__ == "__main__":
	sys.exit(main())
