#!/usr/bin/env python3
"""arg utilities and command entry points."""
from __future__ import annotations
import sys

sys.dont_write_bytecode = True


# 🚩 arg
# GU arg utilities. See --help for commands.
from typing import TYPE_CHECKING, NoReturn

if TYPE_CHECKING:
	import tskit
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
from gu_0_common import load_module


# 🚩 arg_input
"""Validate the shared male-only, one-haplotype-per-person ARG X contract."""
import argparse
import csv


def pfile_prefix(directory, chrom):
	return Path(directory) / ("chrX.male" if chrom == "X" else f"chr{chrom}")


def male_sample_ids(psam):
	with Path(psam).open() as handle:
		header = handle.readline().lstrip("#").split()
		if not {"IID", "SEX"} <= set(header):
			raise ValueError(f"chrX.male PSAM requires IID and SEX: {psam}")
		rows = [dict(zip(header, line.split())) for line in handle if line.strip()]
	names = [row["IID"] for row in rows]
	if not names or len(set(names)) != len(names) or any(row.get("SEX") != "1" for row in rows):
		raise ValueError(f"chrX.male requires unique male samples (SEX=1): {psam}")
	return set(names)


def validate_x_map(sample_map, psam):
	males = male_sample_ids(psam)
	with Path(sample_map).open() as handle:
		rows = list(csv.DictReader(handle, delimiter="\t"))
	names = [row["sample"] for row in rows]
	if (
		not names
		or not set(names) <= males
		or len(set(names)) != len(names)
		or any(row["haplotype"] != "1" for row in rows)
	):
		raise ValueError(
			"chrX ARG must contain only chrX.male samples, each with haplotype=1; "
			"rebuild chrX with arg.sh build --chr X --format trace"
		)
	return len(names)


def arg_input_cli():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--psam", type=Path, required=True)
	parser.add_argument("--sample-map", type=Path)
	args = parser.parse_args()
	count = validate_x_map(args.sample_map, args.psam) if args.sample_map else len(male_sample_ids(args.psam))
	print(f"ARG chrX male-only sample check PASS: samples={count} source={args.psam}")


# 🚩 arg
"""ARG preparation, inference, validation and TRACE export utilities."""
import gzip

# ---- x ----
"""Non-PAR X input handling and explicit biological haplotype identities."""


def x_nonpar(position, build):
	# One-based inclusive non-PAR interval between PAR1 and PAR2.
	start, end = {"37": (2699521, 154931043), "38": (2781480, 155701382)}[str(build)]
	return start <= position <= end


def x_read_sexes(path):
	with open(path) as handle:
		header = handle.readline().lstrip("#").split()
		if "IID" not in header or "SEX" not in header:
			raise ValueError(f"X requires IID and SEX in PSAM: {path}")
		result = {}
		for line in handle:
			if not line.strip():
				continue
			row = dict(zip(header, line.split()))
			if row["IID"] in result or row["SEX"] not in ("1", "2"):
				raise ValueError(f"Duplicate IID or unknown X sex: {row['IID']}")
			result[row["IID"]] = int(row["SEX"])
	return result


def x_haploid_zarr(source, destination):
	"""Flatten only real GT copies, never infer the -2 padding as a sample.

    Uniform haploid storage also works with tsinfer's partitioned matcher.
    Scan every GT chunk to enforce constant per-individual ploidy.
    """
	import numpy as np
	import numcodecs
	import zarr

	src = zarr.open(str(source), mode="r")
	gt = src["call_genotype"]
	names = [str(x) for x in src["sample_id"][:]]
	present = np.asarray(gt[0]) != -2
	ploidies = present.sum(axis=1)
	if np.any((ploidies < 1) | (ploidies > 2)) or np.any(~present[:, 0]):
		raise ValueError("Invalid per-sample X ploidy in VCF-Zarr")
	identities = [(name, h + 1) for name, p in zip(names, ploidies) for h in range(int(p))]
	dest = zarr.open_group(str(destination), mode="w")
	for key in (
		"variant_position",
		"variant_allele",
		"variant_AA",
		"variant_contig",
		"contig_id",
		"contig_length",
	):
		if key in src:
			zarr.copy(src[key], dest, name=key)
	ids = np.asarray([f"arg_hap_{i}" for i in range(len(identities))], dtype=object)
	dest.create_dataset("sample_id", data=ids, object_codec=numcodecs.VLenUTF8())
	width = len(identities)
	chunk = min(gt.chunks[0], max(1, 32_000_000 // (gt.shape[1] * gt.shape[2])))
	out = dest.create_dataset(
		"call_genotype", shape=(gt.shape[0], width, 1), chunks=(chunk, min(width, 256), 1), dtype=gt.dtype
	)
	for start in range(0, gt.shape[0], chunk):
		values = np.asarray(gt[start : start + chunk])
		if np.any((values != -2) != present[None, :, :]):
			raise ValueError(f"X ploidy changes within chromosome at variant chunk {start}")
		if np.any(values < -2) or np.any(values > 1):
			raise ValueError("Non-biallelic GT in prepared X VCF")
		out[start : start + len(values), :, 0] = values[:, present]
	zarr.consolidate_metadata(dest.store)
	return identities


def x_restore_individuals(ts, identities):
	"""Replace storage-only haploid identities with original biological samples."""
	import numpy as np
	import tskit

	if ts.num_samples != len(identities):
		raise ValueError("X ARG sample count differs from real input haplotypes")
	tables = ts.dump_tables()
	individuals = {}
	tables.individuals.clear()
	tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
	node_individual = np.full(ts.num_nodes, tskit.NULL, dtype=np.int32)
	for node, (name, hap) in zip(ts.samples(), identities):
		if name not in individuals:
			individuals[name] = tables.individuals.add_row(metadata={"sample": name})
		node_individual[node] = individuals[name]
	tables.nodes.individual = node_individual
	return tables.tree_sequence()


# ---- needle_x ----
"""Pack real non-PAR X copies into ASMC's two-column storage containers."""
import json


def needle_x_panel(keep, psam, policy, qc):
	sexes = x_read_sexes(psam)
	path = Path(keep)
	rows = path.read_text().splitlines()
	names = [r.split()[-1] for r in rows[1:] if r.strip()]
	if not set(names) <= sexes.keys():
		raise ValueError("Selected X panel has samples absent from PSAM")
	males = [name for name in names if sexes[name] == 1]
	dropped = []
	if len(males) % 2:
		if policy != "drop-last":
			raise ValueError(
				"ASMC requires an even number of real haplotypes. "
				"Use --x-odd-male drop-last to exclude the last selected male, "
				"or provide an even-male sample panel with --keep."
			)
		dropped = males[-1:]
		rows = [rows[0]] + [r for r in rows[1:] if r.split()[-1] not in dropped]
		path.write_text("\n".join(rows) + "\n")
	result = dict(
		policy=policy,
		excluded_males=dropped,
		individuals=len(names) - len(dropped),
		haplotypes=sum(1 if sexes[n] == 1 else 2 for n in names) - len(dropped),
	)
	qc = Path(qc)
	text = json.dumps(result, indent=2) + "\n"
	if not qc.is_file() or qc.read_text() != text:
		temporary = qc.with_name(qc.name + ".next")
		temporary.write_text(text)
		temporary.replace(qc)
	print("Needle X panel: " + json.dumps(result))


def needle_x_pack(prefix, psam, build, keep):
	import numpy as np

	prefix = str(prefix)
	sexes = x_read_sexes(psam)
	sample = Path(prefix + ".sample")
	names = [line.split()[1] for line in sample.read_text().splitlines()[2:] if line.strip()]
	wanted = [
		line.split()[-1]
		for line in Path(keep).read_text().splitlines()
		if line.strip() and not line.startswith("#")
	]
	if names != wanted:
		raise ValueError("X export sample order differs from selected panel")
	identities = [(name, h) for name in names for h in range(1 if sexes[name] == 1 else 2)]
	if len(identities) % 2:
		raise ValueError("Odd real X haplotype count: select an even-male panel before export")
	male = np.array([sexes[name] == 1 for name in names])
	columns = np.array(
		[2 * i + h for i, name in enumerate(names) for h in range(1 if sexes[name] == 1 else 2)]
	)
	haps = Path(prefix + ".haps.gz")
	if not haps.exists():
		haps = Path(prefix + ".haps")
	tmp = haps.with_name(haps.name + ".x.next")
	compressed = haps.suffix == ".gz"
	retained = removed = 0
	try:
		with (
			gzip.open(haps, "rt") if compressed else haps.open() as src,
			gzip.open(tmp, "wt", compresslevel=1) if compressed else tmp.open("w") as dest,
		):
			for line in src:
				z = line.split()
				if len(z) != 5 + 2 * len(names) or z[0].removeprefix("chr") not in ("X", "23"):
					raise ValueError("Malformed X HAPS row")
				if not x_nonpar(int(z[2]), build):
					removed += 1
					continue
				# Vectorize across individuals: full-density X contains billions
				# of calls. ASCII 0/1 are 48/49; PLINK's absent male copy is '-'.
				values = np.frombuffer("".join(z[5:]).encode("ascii"), dtype=np.uint8)
				if len(values) != 2 * len(names):
					raise ValueError(f"Missing/unphased X HAPS at {z[2]}")
				pairs = values.reshape(-1, 2)
				a, b = pairs[:, 0], pairs[:, 1]
				valid_a = (a == 48) | (a == 49)
				valid_b = np.where(male, (b == 45) | (b == a), (b == 48) | (b == 49))
				if not np.all(valid_a & valid_b):
					bad = int(np.flatnonzero(~(valid_a & valid_b))[0])
					raise ValueError(f"Invalid/missing X HAPS: {names[bad]} at {z[2]}")
				gt = values[columns].tobytes().decode("ascii")
				dest.write(" ".join(["X"] + z[1:5]) + " " + " ".join(gt) + "\n")
				retained += 1
		if retained < 2:
			raise ValueError("Too few non-PAR X HAPS variants")
		tmp.replace(haps)
	finally:
		tmp.unlink(missing_ok=True)
	# Pairs are file-format containers only; biological identities live in the sidecar.
	storage_ids = [f"X_pair_{i}" for i in range(len(identities) // 2)]
	sample.write_text("ID_1 ID_2 missing\n0 0 0\n" + "".join(f"{n} {n} 0\n" for n in storage_ids))
	Path(prefix + ".packed.keep").write_text("#IID\n" + "\n".join(storage_ids) + "\n")
	with open(prefix + ".haplotypes.tsv", "w", newline="") as handle:
		writer = csv.writer(handle, delimiter="\t")
		writer.writerow(("eid", "hap"))
		writer.writerows(identities)
	print(
		f"Needle X packed: individuals={len(names)} real_haplotypes={len(identities)} "
		f"nonpar_variants={retained} excluded_outside_nonpar={removed}"
	)


def needle_x_sample_map(tree, identities, output):
	import tskit

	with open(identities) as handle:
		rows = [(r["eid"], int(r["hap"])) for r in csv.DictReader(handle, delimiter="\t")]
	ts = tskit.load(tree)
	if len(rows) != ts.num_samples:
		raise ValueError("X tree sample count differs from biological haplotype map")
	restored = x_restore_individuals(ts, rows)
	temp = Path(str(tree) + ".x.next")
	restored.dump(str(temp))
	temp.replace(tree)
	with open(output, "w", newline="") as handle:
		writer = csv.writer(handle, delimiter="\t")
		writer.writerow(("eid", "hap", "sample_node"))
		writer.writerows((name, hap, int(node)) for (name, hap), node in zip(rows, restored.samples()))


def needle_x_annotate_qc(qc, identities):
	from collections import Counter

	with open(identities) as handle:
		counts = Counter(r["eid"] for r in csv.DictReader(handle, delimiter="\t"))
	path = Path(qc)
	result = json.loads(path.read_text())
	if result["haplotypes"] != sum(counts.values()):
		raise ValueError("X HAPS QC count differs from real haplotypes")
	result["storage_pairs"] = result["individuals"]
	result["individuals"] = len(counts)
	result["haploid_individuals"] = sum(n == 1 for n in counts.values())
	result["diploid_individuals"] = sum(n == 2 for n in counts.values())
	result["storage"] = "ASMC-pairs; biological identities in chrX.haplotypes.tsv"
	path.write_text(json.dumps(result, indent=2) + "\n")


def needle_x_check(tree, identities, output, psam, keep):
	import tskit

	sexes = x_read_sexes(psam)
	names = [r.split()[-1] for r in Path(keep).read_text().splitlines()[1:] if r.strip()]
	expected = [(name, h) for name in names for h in range(1 if sexes[name] == 1 else 2)]
	with open(identities) as handle:
		actual = [(r["eid"], int(r["hap"])) for r in csv.DictReader(handle, delimiter="\t")]
	if actual != expected or len(set(actual)) != len(actual):
		raise ValueError("X biological haplotype identities disagree with panel/sex")
	ts = tskit.load(tree)
	with open(output) as handle:
		mapped = [
			(r["eid"], int(r["hap"]), int(r["sample_node"])) for r in csv.DictReader(handle, delimiter="\t")
		]
	if len(expected) != ts.num_samples or mapped != [
		(n, h, int(node)) for (n, h), node in zip(expected, ts.samples())
	]:
		raise ValueError("X sample map does not match tree samples and real ploidy")
	for name, hap, node in mapped:
		individual = ts.node(node).individual
		if individual < 0 or ts.individual(individual).metadata.get("sample") != name:
			raise ValueError("X tree individual identity disagrees with sample map")
	print(f"Needle X identity CHECK PASS: individuals={len(names)} real_haplotypes={len(expected)}")


def needle_x_main():
	parser = argparse.ArgumentParser(description=__doc__)
	sub = parser.add_subparsers(dest="command", required=True)
	p = sub.add_parser("panel")
	for key in ("keep", "psam", "qc"):
		p.add_argument("--" + key, required=True)
	p.add_argument("--policy", choices=("error", "drop-last"), default="drop-last")
	p = sub.add_parser("pack")
	for key in ("prefix", "psam", "keep", "build"):
		p.add_argument("--" + key, required=True)
	p = sub.add_parser("map")
	for key in ("tree", "identities", "output"):
		p.add_argument("--" + key, required=True)
	p = sub.add_parser("check")
	for key in ("tree", "identities", "output", "psam", "keep"):
		p.add_argument("--" + key, required=True)
	p = sub.add_parser("qc")
	for key in ("qc", "identities"):
		p.add_argument("--" + key, required=True)
	args = vars(parser.parse_args())
	command = args.pop("command")
	{
		"panel": needle_x_panel,
		"pack": needle_x_pack,
		"map": needle_x_sample_map,
		"check": needle_x_check,
		"qc": needle_x_annotate_qc,
	}[command](**args)


# ---- infer ----
"""Build a dated tsinfer/tsdate tree sequence for TRACE from a filtered phased VCF.

The workflow uses explicit generate_ancestors -> match_ancestors -> match_samples
stages rather than tsinfer.infer(). This makes memory use and failure points easier to
inspect, and it allows mmap-backed ancestor generation on local Linux storage.

The input should be produced by ``f/arg.sh prep_vcf``. It must contain only
biallelic variants, phased GTs, and a normalized INFO/AA field. INFO/AA, rather
than VCF REF, is used as the ancestral state.

Large VCF-Zarr intermediates are created below a Linux temporary directory
(default: /tmp), not beside the final output on /mnt/d. This is important in
WSL because Zarr contains many small files and is extremely slow to create or
remove on a Windows-mounted filesystem.
"""
from collections import Counter
import inspect
import os
import pickle
import shutil
import subprocess
import time
import warnings
from typing import Any


def infer_die(message: str) -> "NoReturn":
	raise SystemExit(message)


def infer_command_path(name: str) -> str:
	path = shutil.which(name)
	if path is None:
		infer_die(f"Required command not found: {name}")
	return path


def infer_version_of(module: Any) -> str:
	return str(getattr(module, "__version__", "unknown"))


def infer_extract_sample_ids(vcf_path: Path) -> list[str]:
	"""Read sample IDs with bcftools, avoiding Python gzip/header edge cases."""
	result = subprocess.run(
		[infer_command_path("bcftools"), "query", "-l", str(vcf_path)],
		check=True,
		text=True,
		capture_output=True,
	)
	samples = [line.strip() for line in result.stdout.splitlines() if line.strip()]
	if not samples:
		infer_die(f"No samples found in VCF: {vcf_path}")
	if len(samples) != len(set(samples)):
		infer_die(f"Duplicate sample IDs found in VCF: {vcf_path}")
	return samples


def infer_check_phasing(vcf_path: Path, max_records: int = 2000) -> dict[str, int]:
	"""Check diploid heterozygotes; haploid calls are phase-complete by definition."""
	import pysam

	records = 0
	called = 0
	heterozygous = 0
	unphased_heterozygous = 0
	with pysam.VariantFile(str(vcf_path)) as vcf:
		for record in vcf:
			records += 1
			for call in record.samples.values():
				gt = call.get("GT")
				if gt is None or any(a is None for a in gt):
					continue
				called += 1
				if len(gt) < 2:
					continue
				if len(set(gt)) > 1:
					heterozygous += 1
					if not call.phased:
						unphased_heterozygous += 1
						if unphased_heterozygous <= 5:
							print(
								f"UNPHASED: {record.contig}:{record.pos} GT={gt}",
								flush=True,
							)
			if records >= max_records:
				break
	if records == 0:
		infer_die(f"VCF has no records: {vcf_path}")
	if unphased_heterozygous:
		infer_die(
			f"VCF is not fully phased: found {unphased_heterozygous} unphased "
			f"heterozygous calls in the first {records} records"
		)
	return {
		"records_checked": records,
		"called_genotypes_checked": called,
		"heterozygous_genotypes_checked": heterozygous,
		"unphased_heterozygous_checked": unphased_heterozygous,
	}


def infer_scalar_string(value: Any) -> str:
	"""Convert Zarr string/byte/scalar/length-one arrays to a plain string."""
	if hasattr(value, "tolist"):
		value = value.tolist()
	while isinstance(value, (list, tuple)) and len(value) == 1:
		value = value[0]
	if isinstance(value, bytes):
		value = value.decode("utf-8", errors="replace")
	return str(value)


def infer_normalize_aa(value: Any) -> str:
	aa = infer_scalar_string(value).strip().upper()
	aa = aa.split("|", 1)[0].split(",", 1)[0]
	return aa


def infer_build_ancestral_array(zarr_path: Path, chunk_size: int = 250_000) -> tuple[str, int]:
	"""Create a validated one-base ancestral-state array from variant_AA."""
	import numpy as np
	import zarr

	root = zarr.open(str(zarr_path), mode="r+")
	required = ("variant_AA", "variant_allele")
	missing = [name for name in required if name not in root]
	if missing:
		infer_die(
			"VCF-Zarr is missing "
			+ ", ".join(missing)
			+ ". The VCF must retain INFO/AA; REF is not an acceptable fallback."
		)
	aa_source = root["variant_AA"]
	alleles_source = root["variant_allele"]
	n_variants = int(aa_source.shape[0])
	if n_variants == 0:
		infer_die("VCF-Zarr contains no variants")
	out_name = "variant_ancestral_state"
	if out_name in root:
		del root[out_name]
	chunks = getattr(aa_source, "chunks", None)
	out_chunk = int(chunks[0]) if chunks and chunks[0] else min(chunk_size, n_variants)
	out = root.create_dataset(
		out_name,
		shape=(n_variants,),
		chunks=(max(1, out_chunk),),
		dtype="<U1",
		overwrite=True,
	)

	invalid: list[str] = []
	for start in range(0, n_variants, chunk_size):
		end = min(start + chunk_size, n_variants)
		aa_block = aa_source[start:end]
		allele_block = alleles_source[start:end]
		normalized = np.empty(end - start, dtype="<U1")
		for offset, raw_aa in enumerate(aa_block):
			aa = infer_normalize_aa(raw_aa)
			alleles = {
				infer_scalar_string(x).strip().upper()
				for x in allele_block[offset]
				if infer_scalar_string(x).strip() not in ("", ".")
			}
			if aa not in {"A", "C", "G", "T"} or aa not in alleles:
				if len(invalid) < 10:
					invalid.append(f"variant_index={start + offset} AA={aa!r} alleles={sorted(alleles)}")
				normalized[offset] = "N"
			else:
				normalized[offset] = aa
		out[start:end] = normalized
	if invalid:
		infer_die("Invalid ancestral states remain after VCF preparation. Examples: " + "; ".join(invalid))
	return out_name, n_variants


def infer_vcf_zarr_ploidy(zarr_path: Path) -> int:
	"""Read the fixed GT ploidy dimension written by vcf2zarr."""
	import zarr

	root = zarr.open(str(zarr_path), mode="r")
	if "call_genotype" not in root:
		infer_die("VCF-Zarr is missing call_genotype")
	shape = root["call_genotype"].shape
	if len(shape) != 3 or int(shape[2]) < 1:
		infer_die(f"Unexpected call_genotype shape in VCF-Zarr: {shape}")
	return int(shape[2])


def infer_metadata_dict(obj: Any) -> dict[str, Any]:
	metadata = getattr(obj, "metadata", None)
	if isinstance(metadata, dict):
		return metadata
	if isinstance(metadata, bytes):
		try:
			return json.loads(metadata.decode())
		except Exception:
			return {}
	return {}


def infer_write_sample_map(tree_sequence: Any, sample_ids: list[str], output_path: Path) -> None:
	"""Write actual tree node IDs and map haplotypes to VCF samples."""
	seen: dict[str, int] = {}
	sample_id_set = set(sample_ids)
	with output_path.open("w", encoding="utf-8", newline="") as handle:
		handle.write("tree_node_id\tsample\thaplotype\n")
		for node_id in tree_sequence.samples():
			node = tree_sequence.node(int(node_id))
			candidates: list[str] = []
			objects = [node]
			if node.individual >= 0:
				objects.append(tree_sequence.individual(node.individual))
			for obj in objects:
				metadata = infer_metadata_dict(obj)
				for key in (
					"variant_data_sample_id",
					"sample",
					"sample_id",
					"individual",
					"individual_id",
					"name",
					"id",
				):
					if key in metadata:
						candidates.append(str(metadata[key]))
			sample = next((x for x in candidates if x in sample_id_set), None)
			if sample is None and 0 <= node.individual < len(sample_ids):
				sample = sample_ids[node.individual]
			if sample is None:
				sample = candidates[0] if candidates else f"node_{node_id}"
			seen[sample] = seen.get(sample, 0) + 1
			handle.write(f"{int(node_id)}\t{sample}\t{seen[sample]}\n")


def infer_dump_atomic(tree_sequence: Any, output: Path) -> None:
	temp_output = output.with_name(f"{output.name}.part.{os.getpid()}")
	try:
		tree_sequence.dump(str(temp_output))
		os.replace(temp_output, output)
	finally:
		temp_output.unlink(missing_ok=True)


def infer_parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Infer and date an ARG for TRACE from a phased INFO/AA VCF")
	parser.add_argument("--vcf", required=True, type=Path)
	parser.add_argument("--out", required=True, type=Path)
	parser.add_argument("--sample-map", type=Path)
	parser.add_argument("--qc-json", type=Path)
	parser.add_argument("--chr", default="")
	parser.add_argument("--threads", type=int, default=8)
	parser.add_argument("--tmp-root", type=Path, default=Path("/tmp"))
	parser.add_argument("--recombination-rate", type=float, default=1.0e-8)
	parser.add_argument("--mutation-rate", type=float, default=1.25e-8)
	parser.add_argument("--population-size", type=float, default=30_000)
	parser.add_argument(
		"--tsdate-method",
		default="variational_gamma",
		choices=("variational_gamma", "inside_outside", "maximization"),
	)
	parser.add_argument("--mismatch-ratio", type=float)
	parser.add_argument("--phase-check-records", type=int, default=2000)
	parser.add_argument("--no-progress", action="store_true")
	parser.add_argument(
		"--one-bit",
		action="store_true",
		help="use 1-bit biallelic genotype encoding; requires no missing genotypes",
	)
	parser.add_argument(
		"--sample-match-min-work",
		type=int,
		default=50_000_000,
		help=("genotypes per resumable match-samples partition; 0 disables batching (default: 50000000)"),
	)
	parser.add_argument("--keep-work", action="store_true")
	return parser.parse_args()


def infer_main() -> None:
	args = infer_parse_args()
	if args.threads < 1:
		infer_die("--threads must be positive")
	if args.sample_match_min_work < 0:
		infer_die("--sample-match-min-work must be nonnegative")
	for name, value in (
		("recombination rate", args.recombination_rate),
		("mutation rate", args.mutation_rate),
		("population size", args.population_size),
	):
		if value <= 0:
			infer_die(f"{name} must be positive")
	args.vcf = args.vcf.resolve()
	if not args.vcf.is_file() or args.vcf.stat().st_size == 0:
		infer_die(f"Input VCF is missing or empty: {args.vcf}")
	args.out.parent.mkdir(parents=True, exist_ok=True)
	args.tmp_root.mkdir(parents=True, exist_ok=True)
	tmp_filesystem = subprocess.check_output(
		[infer_command_path("stat"), "-f", "-c", "%T", str(args.tmp_root)], text=True
	).strip()
	if tmp_filesystem in {"v9fs", "9p", "drvfs"}:
		infer_die(
			"ARG --tmp-root must be on the local WSL Linux filesystem, not "
			f"{tmp_filesystem}: {args.tmp_root}. LMDB requires a 1 TiB sparse map; "
			"use /tmp or set GU_ARG_TMP_ROOT to an ext4 path."
		)
	sample_map = args.sample_map or args.out.with_suffix(".sample_map.tsv")
	qc_json = args.qc_json or args.out.with_suffix(".arg_qc.json")
	sample_map.parent.mkdir(parents=True, exist_ok=True)
	qc_json.parent.mkdir(parents=True, exist_ok=True)

	warnings.filterwarnings(
		"ignore",
		message=r"The LMDBStore is deprecated and will be removed in a Zarr-Python version 3.*",
		category=FutureWarning,
	)
	import tsdate
	import tsinfer
	import tskit
	import zarr

	sample_ids = infer_extract_sample_ids(args.vcf)
	phase_qc = infer_check_phasing(args.vcf, args.phase_check_records)
	# Checkpoint identity follows the requested output, not input metadata.
	work = args.tmp_root / f"gu_arg_{args.out.stem}"
	work.mkdir(parents=True, exist_ok=True)
	started = time.time()
	completed = False
	try:
		vcz_path = work / f"{args.out.stem}.vcz"
		inferred_path = work / f"{args.out.stem}.inferred.trees"
		processed_path = work / f"{args.out.stem}.preprocessed.trees"

		print(f"[ARG] Work directory: {work}", flush=True)
		if (vcz_path / ".zmetadata").is_file():
			print(f"[ARG] Reuse VCF-Zarr checkpoint: {vcz_path}", flush=True)
		else:
			if vcz_path.exists():
				shutil.rmtree(vcz_path)
			print(f"[ARG] Convert VCF to VCF-Zarr: {vcz_path}", flush=True)
			subprocess.run(
				[infer_command_path("vcf2zarr"), "convert", str(args.vcf), str(vcz_path)],
				check=True,
			)
		input_ploidy = infer_vcf_zarr_ploidy(vcz_path)
		x_identities = None
		if args.chr.upper().removeprefix("CHR") in ("X", "23"):
			haploid_path = work / f"{args.out.stem}.haploid.vcz"
			x_identities = x_haploid_zarr(vcz_path, haploid_path)
			vcz_path = haploid_path
			print(
				f"[ARG] Non-PAR X: {len(sample_ids)} individuals, {len(x_identities)} real haplotypes",
				flush=True,
			)
		ancestral_array, n_input_variants = infer_build_ancestral_array(vcz_path)
		print(f"[ARG] VCF ploidy: {input_ploidy}", flush=True)

		data = tsinfer.VariantData(str(vcz_path), ancestral_state=ancestral_array)
		# Explicit staged inference is easier to diagnose than the one-call
		# tsinfer.infer() path and permits mmap-backed ancestor generation.
		ancestor_path = work / f"{args.out.stem}.ancestors"
		ancestor_done = work / f"{args.out.stem}.ancestors.done"
		if ancestor_path.is_file() and ancestor_done.is_file():
			print(f"[ARG] Reuse ancestor-data checkpoint: {ancestor_path}", flush=True)
			ancestor_data = tsinfer.load(str(ancestor_path))
		else:
			ancestor_path.unlink(missing_ok=True)
			ancestor_done.unlink(missing_ok=True)
			print("[ARG] tsinfer.generate_ancestors", flush=True)
			gen_kwargs: dict[str, Any] = {
				"path": str(ancestor_path),
				"num_threads": args.threads,
				"mmap_temp_dir": str(work),
			}
			if args.one_bit:
				gen_kwargs["genotype_encoding"] = tsinfer.GenotypeEncoding.ONE_BIT
			ancestor_data = tsinfer.generate_ancestors(data, **gen_kwargs)
			ancestor_done.write_text("complete\n")

		match_kwargs: dict[str, Any] = {
			"recombination_rate": args.recombination_rate,
			"num_threads": args.threads,
		}
		if args.mismatch_ratio is not None:
			match_kwargs["mismatch_ratio"] = args.mismatch_ratio
		ancestors_ts_path = work / f"{args.out.stem}.ancestors.trees"
		if ancestors_ts_path.is_file():
			print(f"[ARG] Reuse ancestor-tree checkpoint: {ancestors_ts_path}", flush=True)
			ancestors_ts = tskit.load(str(ancestors_ts_path))
		else:
			print("[ARG] tsinfer.match_ancestors", flush=True)
			ancestors_ts = tsinfer.match_ancestors(data, ancestor_data, **match_kwargs)
			infer_dump_atomic(ancestors_ts, ancestors_ts_path)

		if inferred_path.is_file():
			print(f"[ARG] Reuse matched-samples checkpoint: {inferred_path}", flush=True)
			inferred_ts = tskit.load(str(inferred_path))
		elif args.sample_match_min_work == 0:
			print("[ARG] tsinfer.match_samples", flush=True)
			inferred_ts = tsinfer.match_samples(data, ancestors_ts, **match_kwargs)
			infer_dump_atomic(inferred_ts, inferred_path)
		else:
			batch_dir = work / "match_samples_batch"
			metadata_path = batch_dir / "metadata.json"
			batch_kwargs: dict[str, Any] = {
				"recombination_rate": args.recombination_rate,
			}
			if args.mismatch_ratio is not None:
				batch_kwargs["mismatch_ratio"] = args.mismatch_ratio
			if not metadata_path.is_file():
				print("[ARG] tsinfer.match_samples batch init", flush=True)
				descriptor = tsinfer.match_samples_batch_init(
					str(batch_dir),
					str(vcz_path),
					ancestral_array,
					str(ancestors_ts_path),
					args.sample_match_min_work,
					**batch_kwargs,
				)
				num_partitions = int(descriptor.num_partitions)
			else:
				metadata = json.loads(metadata_path.read_text())
				num_partitions = int(metadata["num_partitions"])
				print(
					f"[ARG] Reuse match-samples batch checkpoint: {batch_dir}",
					flush=True,
				)
			for partition_index in range(num_partitions):
				partition_path = batch_dir / f"partition_{partition_index}.pkl"
				if partition_path.is_file() and partition_path.stat().st_size > 0:
					try:
						with partition_path.open("rb") as partition_file:
							pickle.load(partition_file)
						continue
					except (EOFError, OSError, pickle.UnpicklingError):
						partition_path.unlink(missing_ok=True)
				print(
					f"[ARG] tsinfer.match_samples partition {partition_index + 1}/{num_partitions}",
					flush=True,
				)
				tsinfer.match_samples_batch_partition(str(batch_dir), partition_index)
			print("[ARG] tsinfer.match_samples batch finalise", flush=True)
			inferred_ts = tsinfer.match_samples_batch_finalise(str(batch_dir))
			infer_dump_atomic(inferred_ts, inferred_path)

		if x_identities is not None:
			inferred_ts = x_restore_individuals(inferred_ts, x_identities)
		if hasattr(tsdate, "preprocess_ts"):
			print("[ARG] tsdate.preprocess_ts", flush=True)
			processed_ts = tsdate.preprocess_ts(inferred_ts)
		else:
			processed_ts = inferred_ts
		infer_dump_atomic(processed_ts, processed_path)

		date_kwargs: dict[str, Any] = {
			"mutation_rate": args.mutation_rate,
			"method": args.tsdate_method,
		}
		# Variational gamma does not accept population_size.
		if args.tsdate_method != "variational_gamma":
			date_kwargs["population_size"] = args.population_size
		date_signature = inspect.signature(tsdate.date)
		if "num_threads" in date_signature.parameters:
			date_kwargs["num_threads"] = args.threads
		if "progress" in date_signature.parameters:
			date_kwargs["progress"] = not args.no_progress
		date_parameters = f"method={args.tsdate_method}, mutation_rate={args.mutation_rate:g}"
		if "population_size" in date_kwargs:
			date_parameters += f", Ne={args.population_size:g}"
		print(f"[ARG] tsdate: {date_parameters}", flush=True)
		dated_ts = tsdate.date(processed_ts, **date_kwargs)
		expected_nodes = len(x_identities) if x_identities is not None else input_ploidy * len(sample_ids)
		if dated_ts.num_samples != expected_nodes:
			infer_die(
				"ARG sample-node count does not match VCF samples × ploidy "
				f"({dated_ts.num_samples} vs {len(sample_ids)} × {input_ploidy} = {expected_nodes})"
			)
		infer_dump_atomic(dated_ts, args.out)

		map_part = sample_map.with_name(f"{sample_map.name}.part.{os.getpid()}")
		infer_write_sample_map(dated_ts, sample_ids, map_part)
		os.replace(map_part, sample_map)

		qc = {
			"input_vcf": str(args.vcf),
			"output_trees": str(args.out.resolve()),
			"sample_map": str(sample_map.resolve()),
			"ancestral_state_source": "INFO/AA (VCF-Zarr variant_AA)",
			"input_variants": n_input_variants,
			"vcf_samples": len(sample_ids),
			"vcf_ploidy": input_ploidy,
			"sample_ploidies": (
				dict(Counter(name for name, _ in x_identities))
				if x_identities is not None
				else dict.fromkeys(sample_ids, input_ploidy)
			),
			"sample_nodes": int(dated_ts.num_samples),
			"sites": int(dated_ts.num_sites),
			"mutations": int(dated_ts.num_mutations),
			"nodes": int(dated_ts.num_nodes),
			"edges": int(dated_ts.num_edges),
			"trees": int(dated_ts.num_trees),
			"sequence_length": float(dated_ts.sequence_length),
			"parameters": {
				"inference_pipeline": "generate_ancestors/match_ancestors/match_samples",
				"one_bit": bool(args.one_bit),
				"threads": args.threads,
				"sample_match_min_work": args.sample_match_min_work,
				"recombination_rate": args.recombination_rate,
				"mutation_rate": args.mutation_rate,
				"population_size": args.population_size,
				"tsdate_method": args.tsdate_method,
				"mismatch_ratio": args.mismatch_ratio,
			},
			"phase_qc": phase_qc,
			"versions": {
				"tsinfer": infer_version_of(tsinfer),
				"tsdate": infer_version_of(tsdate),
				"tskit": infer_version_of(tskit),
				"zarr": infer_version_of(zarr),
			},
			"elapsed_seconds": round(time.time() - started, 3),
		}
		qc_part = qc_json.with_name(f"{qc_json.name}.part.{os.getpid()}")
		qc_part.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
		os.replace(qc_part, qc_json)
		print(f"[ARG] Dated ARG: {args.out}", flush=True)
		print(f"[ARG] TRACE sample map: {sample_map}", flush=True)
		print(f"[ARG] QC: {qc_json}", flush=True)
		completed = True
	finally:
		if args.keep_work:
			print(
				f"[ARG] Keeping work directory by request: {work}",
				flush=True,
			)
		else:
			shutil.rmtree(work, ignore_errors=True)
			if not completed:
				print(f"[ARG] Removed incomplete work directory: {work}", flush=True)


# ---- check ----
"""Fast structural validation for a completed tsinfer ARG and its sidecars."""


load_module("0.common.py")
from gu_0_common import enable_wide_csv_fields


enable_wide_csv_fields()


def check_fail(message: str) -> None:
	import tskit

	raise SystemExit(f"ERROR: ARG post-build check failed: {message}")


def check_resolved(value: str | Path) -> str:
	import tskit

	return str(Path(value).resolve())


def check_main() -> None:
	import tskit

	parser = argparse.ArgumentParser()
	parser.add_argument("--trees", required=True, type=Path)
	parser.add_argument("--sample-map", required=True, type=Path)
	parser.add_argument("--qc-json", required=True, type=Path)
	parser.add_argument("--input-vcf", required=True, type=Path)
	parser.add_argument("--chr", required=True)
	args = parser.parse_args()

	try:
		ts = tskit.load(str(args.trees))
	except Exception as exc:
		check_fail(f"unreadable tree sequence {args.trees}: {exc}")
	if ts.num_trees < 1 or ts.num_sites < 1 or ts.num_samples < 1:
		check_fail(f"empty tree sequence: trees={ts.num_trees} sites={ts.num_sites} samples={ts.num_samples}")
	if ts.sequence_length <= 0:
		check_fail(f"invalid sequence length: {ts.sequence_length}")

	try:
		qc = json.loads(args.qc_json.read_text(encoding="utf-8"))
	except Exception as exc:
		check_fail(f"unreadable QC JSON {args.qc_json}: {exc}")
	expected_paths = {
		"input_vcf": check_resolved(args.input_vcf),
		"output_trees": check_resolved(args.trees),
		"sample_map": check_resolved(args.sample_map),
	}
	for key, expected in expected_paths.items():
		if check_resolved(qc.get(key, "")) != expected:
			check_fail(f"QC {key} does not match {expected}")
	expected_counts = {
		"sample_nodes": ts.num_samples,
		"sites": ts.num_sites,
		"mutations": ts.num_mutations,
		"nodes": ts.num_nodes,
		"edges": ts.num_edges,
		"trees": ts.num_trees,
	}
	for key, expected in expected_counts.items():
		if qc.get(key) != expected:
			check_fail(f"QC {key}={qc.get(key)!r}, tree sequence has {expected}")
	if float(qc.get("sequence_length", -1)) != float(ts.sequence_length):
		check_fail("QC sequence_length does not match tree sequence")
	if qc.get("ancestral_state_source") != "INFO/AA (VCF-Zarr variant_AA)":
		check_fail("QC does not record INFO/AA as the ancestral-state source")

	try:
		with args.sample_map.open(encoding="utf-8", newline="") as handle:
			reader = csv.DictReader(handle, delimiter="\t")
			if reader.fieldnames != ["tree_node_id", "sample", "haplotype"]:
				check_fail(f"unexpected sample-map header: {reader.fieldnames}")
			rows = list(reader)
	except Exception as exc:
		check_fail(f"unreadable sample map {args.sample_map}: {exc}")
	if len(rows) != ts.num_samples:
		check_fail(f"sample-map rows={len(rows)}, expected {ts.num_samples}")
	try:
		mapped_nodes = [int(row["tree_node_id"]) for row in rows]
		if any(not row["sample"] or int(row["haplotype"]) < 1 for row in rows):
			check_fail("sample map contains an empty sample or invalid haplotype")
	except (TypeError, ValueError) as exc:
		check_fail(f"invalid sample-map value: {exc}")
	if len(mapped_nodes) != len(set(mapped_nodes)):
		check_fail("sample map contains duplicate tree node IDs")
	if set(mapped_nodes) != set(map(int, ts.samples())):
		check_fail("sample-map node IDs do not equal tree-sequence sample nodes")
	if args.chr == "X" and "sample_ploidies" not in qc:
		check_fail("X QC lacks biological sample ploidies; rebuild with the updated X workflow")
	if "sample_ploidies" in qc:
		observed = {}
		for row in rows:
			observed.setdefault(row["sample"], []).append(int(row["haplotype"]))
		expected = {name: list(range(1, int(ploidy) + 1)) for name, ploidy in qc["sample_ploidies"].items()}
		if {name: sorted(haps) for name, haps in observed.items()} != expected:
			check_fail("Sample map differs from biological sample ploidies")

	print(
		f"[ARG CHECK] PASS chr={args.chr} trees={ts.num_trees} "
		f"sites={ts.num_sites} samples={ts.num_samples} sequence_length={ts.sequence_length:g}"
	)


# ---- trace ----
"""Export a native ARG tree and sample map to GU TRACE's input contract."""


def trace_atomic_text(path: Path, text: str) -> None:
	import tskit

	path.parent.mkdir(parents=True, exist_ok=True)
	tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
	tmp.write_text(text)
	os.replace(tmp, path)


def trace_read_map(path: Path, method: str) -> list[tuple[int, str, int]]:
	import tskit

	with path.open(newline="") as handle:
		reader = csv.DictReader(handle, delimiter="\t")
		fields = set(reader.fieldnames or [])
		rows = []
		if method == "needle":
			required = {"eid", "hap", "sample_node"}
			if not required <= fields:
				raise SystemExit(f"Needle sample map requires {sorted(required)}: {path}")
			for row in reader:
				rows.append((int(row["sample_node"]), str(row["eid"]), int(row["hap"]) + 1))
		else:
			required = {"tree_node_id", "sample", "haplotype"}
			if not required <= fields:
				raise SystemExit(f"tsinfer sample map requires {sorted(required)}: {path}")
			for row in reader:
				rows.append(
					(
						int(row["tree_node_id"]),
						str(row["sample"]),
						int(row["haplotype"]),
					)
				)
	return rows


def trace_validate(ts: tskit.TreeSequence, rows: list[tuple[int, str, int]], chrom: str) -> dict:
	import tskit

	if min(ts.num_trees, ts.num_samples, ts.num_sites) < 1:
		raise SystemExit(f"Empty tree sequence for chr{chrom}")
	nodes = list(map(int, ts.samples()))
	mapped = [row[0] for row in rows]
	if len(rows) != len(nodes) or set(mapped) != set(nodes) or len(set(mapped)) != len(mapped):
		raise SystemExit(
			f"TRACE sample map does not bijectively cover chr{chrom} sample nodes: "
			f"rows={len(rows)} samples={len(nodes)}"
		)
	positions = [site.position for site in ts.sites()]
	span = (max(positions) - min(positions)) / max(1.0, float(ts.sequence_length))
	if any(hap < 1 for _, _, hap in rows):
		raise SystemExit(f"TRACE haplotype numbers must be one-based for chr{chrom}")
	return {
		"chromosome": str(chrom),
		"trees": ts.num_trees,
		"nodes": ts.num_nodes,
		"samples": ts.num_samples,
		"sites": ts.num_sites,
		"mutations": ts.num_mutations,
		"sequence_length": ts.sequence_length,
		"site_span_fraction": span,
	}


def trace_write_tree(source: Path, dest: Path, ts: tskit.TreeSequence, method: str) -> str:
	import tskit

	dest.parent.mkdir(parents=True, exist_ok=True)
	tmp = dest.with_name(f".{dest.name}.tmp.{os.getpid()}")
	tmp.unlink(missing_ok=True)
	offset = (
		float(ts.metadata.get("offset", 0)) if method == "needle" and isinstance(ts.metadata, dict) else 0
	)
	if offset:
		# Native ARG-Needle coordinates are relative to the first input SNP.
		# GU loci and TRACE maps are chromosome coordinates, so translate once.
		tables = ts.dump_tables()
		tables.sequence_length += offset
		tables.edges.left += offset
		tables.edges.right += offset
		tables.sites.position += offset
		if len(tables.migrations):
			tables.migrations.left += offset
			tables.migrations.right += offset
		tables.metadata_schema = tskit.MetadataSchema.permissive_json()
		tables.metadata = {
			**ts.metadata,
			"offset": 0,
			"source_arg_offset": offset,
			"coordinate_system": "chromosome_bp",
		}
		tables.time_units = "generations"
		tables.tree_sequence().dump(tmp)
		os.replace(tmp, dest)
		return "chromosome_coordinate_copy"
	if ts.time_units == "generations":
		os.symlink(os.path.relpath(source.resolve(), dest.parent.resolve()), tmp)
		os.replace(tmp, dest)
		return "relative_symlink"
	if method != "needle" or ts.time_units not in ("unknown", ""):
		raise SystemExit(f"TRACE requires node times in generations; got {ts.time_units!r}: {source}")
	tables = ts.dump_tables()
	tables.time_units = "generations"
	tables.provenances.add_row(
		record=json.dumps(
			{
				"schema_version": "1.0.0",
				"software": {"name": "arg.trace", "version": "1"},
				"parameters": {
					"source_method": "needle",
					"operation": "declare_ARG_Needle_normalized_node_times_as_generations",
				},
			},
			separators=(",", ":"),
		)
	)
	tables.tree_sequence().dump(tmp)
	os.replace(tmp, dest)
	return "converted_copy"


def trace_main() -> None:
	import tskit

	parser = argparse.ArgumentParser()
	parser.add_argument("--method", required=True, choices=("needle", "tsinfer", "threads"))
	parser.add_argument("--tree", required=True, type=Path)
	parser.add_argument("--sample-map", required=True, type=Path)
	parser.add_argument("--out-dir", required=True, type=Path)
	parser.add_argument("--chr", required=True)
	parser.add_argument("--replace", action="store_true")
	parser.add_argument("--check", action="store_true")
	args = parser.parse_args()

	for path in (args.tree, args.sample_map):
		if not path.is_file() or path.stat().st_size == 0:
			raise SystemExit(f"Missing ARG input: {path}")
	out_tree = args.out_dir / f"chr{args.chr}.trees"
	out_map = args.out_dir / f"chr{args.chr}.sample_map.tsv"
	out_qc = args.out_dir / f"chr{args.chr}.arg_qc.json"
	out_meta = args.out_dir / f"chr{args.chr}.trace.meta.json"
	source_signature = {
		"tree": [str(args.tree.resolve()), args.tree.stat().st_size, args.tree.stat().st_mtime_ns],
		"sample_map": [
			str(args.sample_map.resolve()),
			args.sample_map.stat().st_size,
			args.sample_map.stat().st_mtime_ns,
		],
		"method": args.method,
		"format": "trace",
	}
	outputs = (out_tree, out_map, out_qc, out_meta)
	if args.check:
		for path in outputs:
			if not path.exists() or path.stat().st_size == 0:
				raise SystemExit(f"Missing TRACE-format ARG output: {path}")
		old = json.loads(out_meta.read_text())
		if old.get("source") != source_signature:
			raise SystemExit(
				f"TRACE-format ARG source provenance differs for chr{args.chr}; rerun build with --format trace"
			)
		out_ts = tskit.load(out_tree)
		rows = trace_read_map(out_map, "tsinfer")
		stats = trace_validate(out_ts, rows, args.chr)
		if out_ts.time_units != "generations":
			raise SystemExit(f"TRACE tree time units are not generations: {out_tree}")
		print(
			f"TRACE format CHECK PASS method={args.method} chr{args.chr} "
			f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
		)
		return

	if not args.replace and all(path.exists() and path.stat().st_size > 0 for path in outputs):
		try:
			old = json.loads(out_meta.read_text())
			if old.get("source") == source_signature:
				out_ts = tskit.load(out_tree)
				rows = trace_read_map(out_map, "tsinfer")
				stats = trace_validate(out_ts, rows, args.chr)
				if out_ts.time_units == "generations":
					print(
						f"TRACE format SKIP method={args.method} chr{args.chr} "
						f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
					)
					return
		except (OSError, ValueError, json.JSONDecodeError, tskit.FileFormatError):
			pass

	ts = tskit.load(args.tree)
	rows = trace_read_map(args.sample_map, args.method)
	stats = trace_validate(ts, rows, args.chr)
	storage = trace_write_tree(args.tree, out_tree, ts, args.method)
	map_text = "tree_node_id\tsample\thaplotype\n" + "".join(
		f"{node}\t{sample}\t{hap}\n" for node, sample, hap in rows
	)
	trace_atomic_text(out_map, map_text)
	out_ts = tskit.load(out_tree)
	stats = trace_validate(out_ts, trace_read_map(out_map, "tsinfer"), args.chr)
	if out_ts.time_units != "generations":
		raise SystemExit(f"TRACE tree time units are not generations: {out_tree}")
	qc = {
		**stats,
		"method": args.method,
		"format": "trace",
		"time_units": out_ts.time_units,
		"storage": storage,
	}
	trace_atomic_text(out_qc, json.dumps(qc, indent=2, sort_keys=True) + "\n")
	trace_atomic_text(
		out_meta,
		json.dumps({"source": source_signature, "output": qc}, indent=2, sort_keys=True) + "\n",
	)
	print(
		f"TRACE format PASS method={args.method} chr{args.chr} storage={storage} "
		f"trees={stats['trees']} samples={stats['samples']} sites={stats['sites']}"
	)


# ---- unique ----
"""Drop all records at repeated positions from sorted PLINK Oxford HAPS.

Keep only one pending row in memory. Do not select an arbitrary allele from
split multiallelic sites. Replace the input only after a successful full scan.
"""


def unique_clean(path, chromosome):
	path = Path(path)
	temporary = path.with_name(path.name + ".unique.next")
	total = kept = duplicate_sites = duplicate_records = 0
	previous = None
	pending = None
	count = 0
	try:
		compressed = path.suffix == ".gz"
		with (
			gzip.open(path, "rt") if compressed else path.open() as source,
			(
				gzip.open(temporary, "wt", compresslevel=1)
				if compressed
				else temporary.open("w", newline="\n")
			) as output,
		):
			for number, line in enumerate(source, 1):
				if not line.strip():
					continue
				fields = line.split(maxsplit=5)
				if len(fields) != 6 or fields[0].removeprefix("chr") != str(chromosome):
					raise ValueError(f"HAPS line {number}: malformed row or wrong chromosome")
				position = int(fields[2])
				if position < 1 or (previous is not None and position < previous):
					raise ValueError(
						f"HAPS line {number}: invalid/unsorted position {position} after {previous}"
					)
				total += 1
				if position == previous:
					if count == 1:
						duplicate_sites += 1
						duplicate_records += 1
					duplicate_records += 1
					count += 1
					pending = None
				else:
					if pending is not None:
						output.write(pending)
						kept += 1
					previous, pending, count = position, line, 1
			if pending is not None:
				output.write(pending)
				kept += 1
		if kept < 2:
			raise ValueError(f"Too few unique-position HAPS variants: {kept}")
		os.replace(temporary, path)
	finally:
		temporary.unlink(missing_ok=True)
	return dict(
		input_variants=total,
		retained_variants=kept,
		duplicate_positions=duplicate_sites,
		removed_variants=duplicate_records,
		policy="exclude-all-repeated-positions",
	)


def unique_main():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--haps", required=True)
	parser.add_argument("--chr", required=True)
	parser.add_argument("--out", required=True)
	args = parser.parse_args()
	result = unique_clean(args.haps, args.chr)
	Path(args.out).write_text(json.dumps(result, indent=2) + "\n")
	print("HAPS position QC: " + json.dumps(result))


def cache_main():
	#!/usr/bin/env python3
	"""Stable request fingerprints. Files use canonical names and nanosecond stat data."""
	import argparse
	import hashlib
	import json
	from pathlib import Path

	p = argparse.ArgumentParser()
	p.add_argument("--file", action="append", default=[])
	p.add_argument("--text-file", action="append", default=[])
	p.add_argument("--value", action="append", default=[])
	a = p.parse_args()
	files = []
	for name in a.file + a.text_file:
		path = Path(name).resolve(strict=True)
		s = path.stat()
		row = [str(path), s.st_size, s.st_mtime_ns]
		if name in a.text_file:
			row.append(hashlib.sha256(path.read_bytes()).hexdigest())
		files.append(row)
	print(json.dumps({"schema": 2, "values": a.value, "files": files}, sort_keys=True))


def validate_vcf_main():
	#!/usr/bin/env python3
	"""Full, streaming phase/AA/position validation before committing prepared VCF."""
	import argparse
	import pysam

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("vcf")
	parser.add_argument("--chr", default="")
	parser.add_argument("--build", choices=("37", "38"), default="37")
	parser.add_argument("--x-psam")
	args = parser.parse_args()
	with pysam.VariantFile(args.vcf) as vcf:
		ploidies = None
		if args.chr == "X":
			if args.x_psam:
				sexes = x_read_sexes(args.x_psam)
				if not set(vcf.header.samples) <= sexes.keys():
					raise SystemExit("X VCF contains samples absent from PSAM")
				ploidies = {n: 1 if sexes[n] == 1 else 2 for n in vcf.header.samples}
		previous = None
		count = 0
		for r in vcf:
			count += 1
			if previous and (r.contig != previous[0] or r.pos <= previous[1]):
				raise SystemExit(f"Non-unique/unsorted positions: {r.contig}:{r.pos}")
			previous = (r.contig, r.pos)
			if args.chr == "X":
				if r.contig.removeprefix("chr") not in ("X", "23") or not x_nonpar(r.pos, args.build):
					raise SystemExit(f"Expected non-PAR X: {r.contig}:{r.pos}")
				if ploidies is None:
					ploidies = {n: len(c["GT"]) for n, c in r.samples.items()}
			aa = r.info.get("AA")
			if (
				len(r.alleles) != 2
				or aa not in r.alleles
				or any(x not in ("A", "C", "G", "T") for x in r.alleles)
			):
				raise SystemExit(f"Invalid alleles/AA: {r.contig}:{r.pos}")
			for name, call in r.samples.items():
				gt = call.get("GT")
				if gt is None or len(gt) not in (1, 2):
					raise SystemExit(f"Invalid ploidy at {r.pos}: {name}")
				if ploidies is not None and len(gt) != ploidies[name]:
					raise SystemExit(
						f"Incorrect/changing non-PAR X ploidy at {r.pos}: {name}; expected {ploidies[name]}, got {len(gt)}. Supply a VCF with haploid males and diploid females."
					)
				if any(x is not None and x not in (0, 1) for x in gt):
					raise SystemExit("Invalid GT allele")
				if len(gt) == 2 and None not in gt and gt[0] != gt[1] and not call.phased:
					raise SystemExit(f"Unphased heterozygote at {r.contig}:{r.pos}: {name}")
		if count == 0:
			raise SystemExit("No variants")
	print(f"Validated all {count} VCF records: phased GT, valid AA, unique positions")


arg_COMMANDS = {
	"cache": cache_main,
	"validate-vcf": validate_vcf_main,
	"unique-haps": unique_main,
	"needle-x": needle_x_main,
	"infer-tsinfer": infer_main,
	"check-tsinfer": check_main,
	"trace": trace_main,
}


def arg_main():
	if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
		print("Usage: arg.py {" + ",".join(arg_COMMANDS) + "} [options]")
		return
	command = sys.argv.pop(1)
	if command not in arg_COMMANDS:
		raise SystemExit("Unknown ARG command: " + command)
	arg_COMMANDS[command]()


def arg_cli():
	arg_main()


arg_COMMANDS["input"] = arg_input_cli


BASE_SUBCOMMANDS = {}


def base_main():
	import sys

	return arg_cli()


sys.dont_write_bytecode = True


# 🚩 arg_singer
"""Phased WGS -> SINGER posterior windows -> chromosome-coordinate ARG samples.

SINGER is an ARG backend, not an independent introgression caller. Keep each
retained posterior draw separate; never treat windows as posterior replicates.
"""

import fcntl, hashlib, math, random, re, runpy
from itertools import groupby
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading, signal
import io, tempfile
from contextlib import contextmanager

singer_HERE = Path(__file__).resolve().parent


def singer_stamp(p):
	p = Path(p).resolve(strict=True)
	s = p.stat()
	return [str(p), s.st_size, s.st_mtime_ns]


def singer_save(p, data):
	p = Path(p)
	p.parent.mkdir(parents=True, exist_ok=True)
	q = p.with_name(p.name + ".next")
	q.write_text(json.dumps(data, indent=2) + "\n")
	q.replace(p)


def singer_valid(p, request=None):
	try:
		d = json.loads(Path(p).read_text())
		return (request is None or d["request"] == request) and all(
			singer_stamp(x[0]) == x for x in d["outputs"]
		)
	except (OSError, ValueError, KeyError):
		return False


def singer_record(p, request, files):
	singer_save(p, dict(request=request, outputs=[singer_stamp(x) for x in files]))


def singer_run(cmd, log, stop=None):
	import shlex

	print("[SINGER] " + shlex.join(map(str, cmd)), flush=True)
	with Path(log).open("a") as f:
		f.write("\n" + shlex.join(map(str, cmd)) + "\n")
		f.flush()
		if stop is not None and stop.is_set():
			raise RuntimeError("SINGER cancelled")
		proc = subprocess.Popen(list(map(str, cmd)), stdout=f, stderr=subprocess.STDOUT)
		try:
			while True:
				if stop is not None and stop.is_set():
					raise RuntimeError("SINGER cancelled")
				try:
					rc = proc.wait(timeout=1)
					break
				except subprocess.TimeoutExpired:
					pass  # Keep cancellation responsive without periodic console output.
			if rc:
				raise RuntimeError(f"command failed ({rc}); see {log}")
		finally:
			if proc.poll() is None:
				proc.terminate()
				try:
					proc.wait(timeout=10)
				except subprocess.TimeoutExpired:
					proc.kill()
					proc.wait()


def singer_parse():
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("--method", choices=["singer"], default="singer")
	p.add_argument("--action", choices=["all", "build", "prepare", "check"], default="build")
	p.add_argument("--format", choices=["native", "trace"], default="native")
	p.add_argument("--dir-gen", default="/mnt/f/gen/1kg/37")
	p.add_argument("--dir-vcf", "--target-vcf-dir", dest="dir_vcf")
	p.add_argument("--arg-dir")
	p.add_argument("--gen4arg-dir")
	p.add_argument("--sample-file", "--sample-panel", dest="sample_file")
	p.add_argument("--keep")
	p.add_argument("--chr", "--chrs", dest="chr", default="1-22,X")
	p.add_argument("--grch", default="37")
	p.add_argument("--replace", type=str.upper, choices=["TRUE", "FALSE"], default="FALSE")
	p.add_argument("--threads", type=int, default=1, help="total CPU budget; each SINGER window uses one CPU")
	p.add_argument(
		"--jobs",
		type=int,
		default=1,
		help="maximum simultaneous SINGER windows, also limited by CPU and memory budgets",
	)
	p.add_argument(
		"--singer-worker-memory-gb",
		type=float,
		default=24,
		help="reserved GiB per active window including conversion (default: 24 for 5 Mb / ~1000 haplotypes); raise for larger inputs",
	)
	p.add_argument("--features", type=str.upper, choices=["FALSE"], default="FALSE")
	p.add_argument("--sample-match-min-work", type=int, default=50000000)
	p.add_argument("--singer-populations", default="GBR,CHB,ITU,YRI,LWK")
	p.add_argument("--singer-per-pop", type=int, default=100)
	p.add_argument("--singer-samples", type=int, default=50)
	p.add_argument(
		"--singer-burnin",
		type=int,
		default=50,
		help="discard this many saved draws, each after --singer-thin iterations",
	)
	p.add_argument("--singer-thin", type=int, default=20)
	p.add_argument("--singer-window-bp", type=int, default=5000000)
	p.add_argument("--singer-min-sites", type=int, default=20)
	p.add_argument("--singer-ne", type=float, default=10000)
	p.add_argument("--mutation-rate", type=float, default=1.25e-8)
	p.add_argument("--singer-recombination-rate", type=float, default=1.25e-8)
	p.add_argument("--singer-start", type=int, default=0)
	p.add_argument("--singer-end", type=int)
	p.add_argument("--seed", type=int, default=20260909)
	p.add_argument("--mac-min", type=int, default=1)
	p.add_argument("--missing-max", type=float, choices=[0.0], default=0)
	a = p.parse_args()
	if os.environ.get("REFGEN_GEN4ARG_ONLY") == "1":
		a.action = "prepare"
	for k in [
		"threads",
		"jobs",
		"singer_samples",
		"singer_thin",
		"singer_window_bp",
		"singer_min_sites",
		"mac_min",
	]:
		if getattr(a, k) < 1:
			p.error(k + " must be positive")
	for k in ["singer_ne", "mutation_rate", "singer_recombination_rate", "singer_worker_memory_gb"]:
		if not math.isfinite(getattr(a, k)) or getattr(a, k) <= 0:
			p.error(k + " must be finite and positive")
	if min(a.singer_burnin, a.singer_per_pop, a.singer_start, a.seed) < 0:
		p.error("counts, coordinates and seed must be nonnegative")
	if a.singer_end is not None and a.singer_end <= a.singer_start:
		p.error("end must exceed start")
	a.grch = {"b37": "37", "grch37": "37", "hg19": "37", "b38": "38", "grch38": "38", "hg38": "38"}.get(
		a.grch.lower(), a.grch
	)
	if a.grch not in ("37", "38"):
		p.error("invalid genome build")
	chrs = []
	for s in re.split("[,; ]+", a.chr):
		s = re.sub("^chr", "", s, flags=re.I).upper()
		if not s:
			continue
		(
			chrs.extend(map(str, range(int(s.split("-")[0]), int(s.split("-")[1]) + 1)))
			if "-" in s
			else chrs.append("X" if s == "23" else s)
		)
	a.chrs = list(dict.fromkeys(chrs))
	if not a.chrs or not set(a.chrs) <= set(map(str, range(1, 23))) | {"X"}:
		p.error("invalid chromosome selection")
	root = Path(a.dir_gen).resolve()
	if root.name in ("vcf", "pfile", "hap"):
		root = root.parent
	a.root = root
	a.vcf = Path(a.dir_vcf or root / "vcf").resolve()
	a.prep = Path(a.gen4arg_dir or str(a.vcf) + ".4arg.singer").resolve()
	a.output = Path(a.arg_dir or root / "arg.singer").resolve()
	a.panel = Path(a.sample_file or root / "samples.txt").resolve()
	if len({a.vcf, a.prep, a.output}) != 3:
		p.error("source, prepared and ARG directories must differ")
	return a


def singer_sample_panel(a, vcf, chrom):
	available = list(vcf.header.samples)
	meta = {}
	if a.panel.is_file():
		with a.panel.open() as f:
			for r in csv.DictReader(f, delimiter="\t"):
				name = r.get("sample", r.get("IID", r.get("#IID", "")))
				if name:
					meta[name] = r
	if a.keep:
		wanted = [s.strip().split()[0] for s in Path(a.keep).read_text().splitlines() if s.strip()]
		if len(set(wanted)) != len(wanted):
			raise ValueError("duplicate --keep sample")
		if set(wanted) - set(available):
			raise ValueError("--keep contains samples absent from VCF")
		keep = set(wanted)
	else:
		pops = a.singer_populations.split(",")
		groups = {}
		for name in available:
			pop = meta.get(name, {}).get("pop", "UNKNOWN")
			if pops == ["ALL"] or pop in pops:
				groups.setdefault(pop, []).append(name)
		keep = set()
		for pop, names in sorted(groups.items()):
			# Hash ranking is stable across chromosomes and independent of VCF order.
			names.sort(key=lambda n: hashlib.sha256(f"{a.seed}:{n}".encode()).hexdigest())
			keep.update(names[: a.singer_per_pop] if a.singer_per_pop else names)
	selected = [s for s in available if s in keep]
	if not selected:
		raise ValueError("no selected samples; provide --sample-file / --keep / --singer-populations ALL")
	if chrom == "X":
		bad = [s for s in selected if str(meta.get(s, {}).get("sex", "")).lower() not in ("1", "m", "male")]
		if bad:
			raise ValueError("X requires male-only VCF/panel; nonmale or unknown sex: " + ",".join(bad[:5]))
	if len(selected) * (1 if chrom == "X" else 2) < 4:
		raise ValueError("SINGER requires at least four selected haplotypes in this workflow")
	return selected, meta


def singer_source_vcf(a, chrom):
	return a.vcf / ("chrX.male.vcf.gz" if chrom == "X" else f"chr{chrom}.vcf.gz")


def singer_unique_snps(records, start, end, dropped):
	"""Exclude every PASS biallelic SNP at a repeated position, before GT filters."""

	def candidates():
		for r in records:
			# Indexed fetch can include an earlier structural variant overlapping this window.
			if not start < r.pos <= end:
				continue
			reason = None
			if set(r.filter) not in (set(), {"PASS"}):
				reason = "filter"
			elif len(r.alleles) != 2 or any(x not in ("A", "C", "G", "T") for x in r.alleles):
				reason = "not_biallelic_snp"
			if reason:
				dropped[reason] = dropped.get(reason, 0) + 1
			else:
				yield r

	for _, group in groupby(candidates(), key=lambda r: r.pos):
		first = next(group)
		count = 1 + sum(1 for _ in group)
		if count > 1:
			dropped["duplicate_position"] = dropped.get("duplicate_position", 0) + count
		else:
			yield first


def singer_prep_request(a, chrom):
	src = singer_source_vcf(a, chrom)
	return dict(
		source=singer_stamp(src),
		index=singer_stamp(str(src) + (".tbi" if Path(str(src) + ".tbi").is_file() else ".csi")),
		code=singer_stamp(__file__),
		panel=singer_stamp(a.panel) if a.panel.is_file() else None,
		keep=singer_stamp(a.keep) if a.keep else None,
		populations=a.singer_populations,
		per_pop=a.singer_per_pop,
		seed=a.seed,
		window=a.singer_window_bp,
		start=a.singer_start,
		end=a.singer_end,
		min_sites=a.singer_min_sites,
		mac=a.mac_min,
		build=a.grch,
	)


@contextmanager
def singer_compressed_vcf_output(path):
	temporary = path.with_name(path.name + ".next")
	try:
		with io.TextIOWrapper(pysam.BGZFile(str(temporary), "w")) as out:
			yield out
		temporary.replace(path)
	finally:
		temporary.unlink(missing_ok=True)


def singer_prepared_vcf(prefix):
	compressed = Path(str(prefix) + ".vcf.gz")
	return compressed if compressed.exists() else Path(str(prefix) + ".vcf")


def singer_remove_native_text(wd):
	for path in wd.glob("arg_*.txt*"):
		path.unlink(missing_ok=True)


def singer_prepare(a, chrom):
	dest = a.prep / f"chr{chrom}"
	done = dest / "complete.json"
	req = singer_prep_request(a, chrom)
	if a.replace == "FALSE" and singer_valid(done, req):
		return json.loads((dest / "windows.json").read_text())
	done.unlink(missing_ok=True)
	dest.mkdir(parents=True, exist_ok=True)
	v = pysam.VariantFile(str(singer_source_vcf(a, chrom)))
	matches = [
		c for c in v.header.contigs if c.removeprefix("chr") in (chrom, "23" if chrom == "X" else chrom)
	]
	if len(matches) != 1:
		raise ValueError(f"expected one chr{chrom} contig")
	contig = matches[0]
	length = v.header.contigs[contig].length
	if not length:
		raise ValueError("VCF contig must declare chromosome length")
	samples, meta = singer_sample_panel(a, v, chrom)
	v.subset_samples(samples)
	lo = a.singer_start
	hi = min(a.singer_end or length, length)
	if chrom == "X":
		lo = max(lo, 2699520 if a.grch == "37" else 2781479)
		hi = min(hi, 154931043 if a.grch == "37" else 155701382)
	if hi <= lo:
		raise ValueError("empty selected interval")
	ploidy = 1 if chrom == "X" else 2
	windows = []
	files = []
	qc = []
	for start in range(lo, hi, a.singer_window_bp):
		end = min(start + a.singer_window_bp, hi)
		prefix = dest / f"w{start}_{end}"
		path = Path(str(prefix) + ".vcf.gz")
		count = 0
		dropped = {}
		with singer_compressed_vcf_output(path) as out:
			out.write("##fileformat=VCFv4.2\n" + f"##contig=<ID={contig},length={end - start}>\n")
			out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n')
			out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
			for r in singer_unique_snps(v.fetch(contig, start, end), start, end, dropped):
				reason = None
				gts = []
				flat = []
				for name in samples:
					call = r.samples[name]
					gt = call.get("GT")
					if gt is None or any(x is None for x in gt):
						reason = "missing"
						break
					if len(gt) != ploidy:
						raise ValueError(f"wrong ploidy {name} at {r.pos}")
					if ploidy == 2 and len(set(gt)) > 1 and not call.phased:
						raise ValueError(f"unphased heterozygote {name} at {r.pos}")
					gts.append("|".join(map(str, gt)))
					flat.extend(gt)
				if not reason and min(sum(flat), len(flat) - sum(flat)) < a.mac_min:
					reason = "mac"
				if reason:
					dropped[reason] = dropped.get(reason, 0) + 1
					continue
				out.write(
					f"{contig}\t{r.pos - start}\t.\t{r.ref}\t{r.alts[0]}\t.\tPASS\t.\tGT\t"
					+ "\t".join(gts)
					+ "\n"
				)
				count += 1
		qc.append(dict(start=start, end=end, sites=count, dropped=dropped))
		print(f"[SINGER] chr{chrom} window {start}-{end}: sites={count} dropped={dropped}", flush=True)
		if count < a.singer_min_sites:
			path.unlink()
			continue
		pysam.tabix_index(str(path), preset="vcf", force=True)
		files.extend([path, Path(str(path) + ".tbi")])
		windows.append(dict(start=start, end=end, sites=count, prefix=str(prefix)))
	v.close()
	if not windows:
		raise ValueError("no windows with sufficient polymorphic sites")
	info = dict(
		chrom=chrom,
		length=length,
		ploidy=ploidy,
		samples=samples,
		populations={s: meta.get(s, {}).get("pop", "UNKNOWN") for s in samples},
		windows=windows,
		qc=qc,
		regional=a.singer_start != 0 or a.singer_end is not None,
	)
	singer_save(dest / "windows.json", info)
	singer_record(done, req, files + [dest / "windows.json"])
	for old in dest.glob("w*.vcf*"):
		if old not in files:
			old.unlink()
	print(f"[SINGER] prepared chr{chrom}: samples={len(samples)} windows={len(windows)}", flush=True)
	return info


def singer_load_expected(prefix, ploidy):
	"""Keep the full expected genotype matrix compact and load it once per window."""
	positions = []
	rows = []
	with pysam.VariantFile(str(singer_prepared_vcf(prefix))) as f:
		n = len(f.header.samples) * ploidy
		for r in f:
			positions.append(r.pos - 1)
			rows.append(np.fromiter((x for c in r.samples.values() for x in c["GT"]), dtype=np.int8, count=n))
	if not rows:
		raise ValueError("empty SINGER input")
	return np.asarray(positions, dtype=np.float64), np.stack(rows)


def singer_validate_window(ts, prefix, ploidy, expected=None):
	"""Check every leaf and site; explicit alleles make integer comparison exact."""
	positions, genotypes = singer_load_expected(prefix, ploidy) if expected is None else expected
	if ts.num_samples != genotypes.shape[1]:
		raise ValueError("SINGER sample count changed")
	if ts.num_sites != len(positions):
		raise ValueError("SINGER lost input sites")
	for i, v in enumerate(ts.variants(alleles=("0", "1"))):
		if v.site.position != positions[i]:
			raise ValueError(f"SINGER coordinate mismatch {v.site.position}")
		if not np.array_equal(v.genotypes, genotypes[i]):
			raise ValueError(f"SINGER leaf genotype mismatch at {positions[i]}")


def singer_merge_windows(chrom, info, sources, dest, validated=False):
	tables = tskit.TableCollection(info["length"])
	tables.time_units = "generations"
	tables.metadata_schema = tskit.MetadataSchema.permissive_json()
	tables.metadata = dict(
		method="singer",
		chromosome=chrom,
		coordinate_system="0-based chromosome bp",
		regional=info["regional"],
	)
	tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()
	tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
	n = len(info["samples"]) * info["ploidy"]
	for name in info["samples"]:
		ind = tables.individuals.add_row(metadata=dict(sample=name, population=info["populations"][name]))
		for h in range(info["ploidy"]):
			tables.nodes.add_row(
				flags=tskit.NODE_IS_SAMPLE,
				time=0,
				individual=ind,
				metadata=dict(sample=name, haplotype=h + 1),
			)
	for w, source in zip(info["windows"], sources):
		ts = tskit.load(source)
		if not validated:
			singer_validate_window(ts, w["prefix"], info["ploidy"])
		if ts.sequence_length != w["end"] - w["start"]:
			raise ValueError("SINGER window length mismatch")
		mapping = {int(node): j for j, node in enumerate(ts.samples())}
		for node in ts.nodes():
			if not node.is_sample():
				mapping[node.id] = tables.nodes.add_row(flags=node.flags, time=node.time)
		for e in ts.edges():
			tables.edges.add_row(
				e.left + w["start"], e.right + w["start"], mapping[e.parent], mapping[e.child]
			)
		for site in ts.sites():
			s = tables.sites.add_row(site.position + w["start"], site.ancestral_state)
			for mut in site.mutations:
				tables.mutations.add_row(s, mapping[mut.node], mut.derived_state)
	tables.sort()
	tables.build_index()
	tables.compute_mutation_parents()
	tables.compute_mutation_times()
	ts = tables.tree_sequence()
	if ts.num_samples != n or not ts.num_mutations:
		raise ValueError("invalid merged ARG")
	tmp = dest.with_suffix(".next")
	try:
		ts.dump(tmp)
		tmp.replace(dest)
	finally:
		tmp.unlink(missing_ok=True)
	return dict(
		sites=ts.num_sites, trees=ts.num_trees, samples=ts.num_samples, sequence_length=ts.sequence_length
	)


def singer_memory_budget_gb():
	"""Use available host RAM and all visible cgroup ancestor limits."""
	gib = 1024**3
	limits = []
	for line in Path("/proc/meminfo").read_text().splitlines():
		if line.startswith("MemAvailable:"):
			limits.append(int(line.split()[1]) * 1024 * 0.85 / gib)
	if os.environ.get("ARG_MEMORY_LIMIT_GB"):
		limits.append(float(os.environ["ARG_MEMORY_LIMIT_GB"]) * 0.85)
	for line in Path("/proc/self/cgroup").read_text().splitlines():
		if line.startswith("0::"):
			root = Path("/sys/fs/cgroup")
			path = root / line[3:].lstrip("/")
			while path == root or root in path.parents:
				try:
					cap = (path / "memory.max").read_text().strip()
					used = int((path / "memory.current").read_text())
					if cap != "max":
						limits.append(max(0, int(cap) - used) * 0.85 / gib)
				except (OSError, ValueError):
					pass
				if path == root:
					break
				path = path.parent
	return min(limits) if limits else 24.0


def singer_window_workers(a, count):
	capacity = singer_memory_budget_gb()
	cpu = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1
	workers = min(count, a.jobs, a.threads, cpu, max(1, int(capacity / a.singer_worker_memory_gb)))
	if capacity < a.singer_worker_memory_gb:
		print(
			f"[SINGER] WARNING: available budget {capacity:.1f} GiB is below the per-window reservation; one window may exceed the memory cap",
			flush=True,
		)
	return workers


def singer_infer_window(a, info, home, root, req, i, w, stop, read_ts):
	if stop.is_set():
		raise RuntimeError("SINGER cancelled")
	n = a.singer_burnin + a.singer_samples
	wd = root / f"w{w['start']}_{w['end']}"
	wd.mkdir(exist_ok=True)
	prefix = wd / "arg"
	mark = wd / "complete.json"
	request = dict(req, window=w)
	if a.replace == "FALSE" and singer_valid(mark, request):
		return
	mark.unlink(missing_ok=True)
	for partial in wd.glob("posterior_*.trees"):
		partial.unlink()
	raw = wd / "raw.complete.json"
	if a.replace == "TRUE" or not singer_valid(raw, request):
		raw.unlink(missing_ok=True)
		singer_remove_native_text(wd)
		try:
			# SINGER accepts a filename prefix, so expand only this active window.
			with tempfile.TemporaryDirectory(prefix="input.", dir=wd) as temporary:
				input_prefix = Path(temporary) / "window"
				source = singer_prepared_vcf(w["prefix"])
				with (
					gzip.open(source, "rb") if source.suffix == ".gz" else source.open("rb") as inp,
					Path(str(input_prefix) + ".vcf").open("wb") as out,
				):
					shutil.copyfileobj(inp, out)
				singer_run(
					[
						home / "singer",
						"-Ne",
						a.singer_ne,
						"-m",
						a.mutation_rate,
						"-r",
						a.singer_recombination_rate,
						"-input",
						input_prefix,
						"-output",
						prefix,
						"-start",
						0,
						"-end",
						w["end"] - w["start"] + 1,
						"-polar",
						0.5,
						"-n",
						n,
						"-thin",
						a.singer_thin,
						"-ploidy",
						info["ploidy"],
						"-seed",
						a.seed + i,
					],
					wd / "run.log",
					stop,
				)
			raw_files = []
			for j in range(n):
				for kind in ("nodes", "branches", "muts"):
					path = Path(str(prefix) + f"_{kind}_{j}.txt")
					if j < a.singer_burnin:
						path.unlink(missing_ok=True)
						continue
					packed = Path(str(path) + ".gz")
					with path.open("rb") as inp, gzip.open(packed, "wb", compresslevel=1) as out:
						shutil.copyfileobj(inp, out)
					path.unlink()
					raw_files.append(packed)
			for unused in wd.glob("arg_*.txt"):
				unused.unlink()
			singer_record(raw, request, raw_files)
		except BaseException:
			singer_remove_native_text(wd)
			raise
	# Compressed native draws allow conversion to resume after a failure.
	try:
		files = []
		expected = singer_load_expected(w["prefix"], info["ploidy"])
		for j in range(a.singer_burnin, n):
			if stop.is_set():
				raise RuntimeError("SINGER cancelled")
			with tempfile.TemporaryDirectory(prefix="convert.", dir=wd) as temporary:
				paths = {}
				for kind in ("nodes", "branches", "muts"):
					paths[kind] = Path(temporary) / kind
					source = Path(str(prefix) + f"_{kind}_{j}.txt.gz")
					with gzip.open(source, "rb") as inp, paths[kind].open("wb") as out:
						shutil.copyfileobj(inp, out)
				tables = read_ts(str(paths["nodes"]), str(paths["branches"]))
				mutations = np.atleast_2d(np.loadtxt(paths["muts"]))
			previous = None
			for row in mutations:
				if row[0] != previous:
					site = tables.sites.add_row(position=row[0], ancestral_state="0")
					previous = row[0]
				tables.mutations.add_row(site=site, node=int(row[1]), derived_state=str(int(row[3])))
			tables.sort()
			tables.build_index()
			tables.compute_mutation_parents()
			tables.compute_mutation_times()
			ts = (
				tables.tree_sequence().keep_intervals([[1, w["end"] - w["start"] + 1]], simplify=False).trim()
			)
			singer_validate_window(ts, w["prefix"], info["ploidy"], expected)
			file = wd / f"posterior_{j}.trees"
			ts.dump(file)
			files.append(file)
		(wd / "convert.log").write_text(
			f"Validated {len(files)} retained posterior draws against every input genotype.\n"
		)
		singer_record(mark, request, files + [wd / "run.log", wd / "convert.log"])
		singer_remove_native_text(wd)
		raw.unlink(missing_ok=True)
	except BaseException:
		for partial in wd.glob("posterior_*.trees"):
			partial.unlink(missing_ok=True)
		raise


def singer_infer(a, chrom, info, home):
	root = a.output / "native" / f"chr{chrom}"
	root.mkdir(parents=True, exist_ok=True)
	done = root / "complete.json"
	req = dict(
		prepared=singer_stamp(a.prep / f"chr{chrom}" / "complete.json"),
		code=singer_stamp(__file__),
		binary=singer_stamp(home / "singer"),
		converter=singer_stamp(home / "convert_to_tskit"),
		samples=a.singer_samples,
		burnin=a.singer_burnin,
		thin=a.singer_thin,
		Ne=a.singer_ne,
		mutation=a.mutation_rate,
		recombination=a.singer_recombination_rate,
		seed=a.seed,
	)
	if a.replace == "FALSE" and singer_valid(done, req):
		print(f"[SINGER] reuse chr{chrom}", flush=True)
		return
	done.unlink(missing_ok=True)
	n = a.singer_burnin + a.singer_samples
	workers = singer_window_workers(a, len(info["windows"]))
	print(
		f"[SINGER] chr{chrom}: {workers} concurrent windows; CPU budget={a.threads}, requested jobs={a.jobs}, reserve={a.singer_worker_memory_gb:g} GiB/window",
		flush=True,
	)
	stop = threading.Event()
	read_ts = runpy.run_path(str(home / "convert_to_tskit"))["read_ts"]
	pool = ThreadPoolExecutor(max_workers=workers)
	futures = []
	try:
		for i, w in enumerate(info["windows"]):
			futures.append(pool.submit(singer_infer_window, a, info, home, root, req, i, w, stop, read_ts))
		for completed, future in enumerate(as_completed(futures), 1):
			future.result()
			print(f"[SINGER] chr{chrom}: windows complete {completed}/{len(futures)}", flush=True)
	except BaseException:
		stop.set()
		for future in futures:
			future.cancel()
		raise
	finally:
		pool.shutdown(wait=True, cancel_futures=True)
	# Only merge outputs whose complete genotype validation and file stamps pass.
	for w in info["windows"]:
		mark = root / f"w{w['start']}_{w['end']}" / "complete.json"
		if not singer_valid(mark, dict(req, window=w)):
			raise ValueError(f"changed window outputs: {mark}")
	merged = a.output / "posterior"
	merged.mkdir(exist_ok=True)
	files = []
	stats = []
	# Remove only stale posterior exports for this chromosome on explicit rebuild.
	for p in merged.glob(f"chr{chrom}.sample*.trees"):
		p.unlink()
	for j in range(a.singer_burnin, n):
		sources = [root / f"w{w['start']}_{w['end']}" / f"posterior_{j}.trees" for w in info["windows"]]
		dest = merged / f"chr{chrom}.sample{j - a.singer_burnin:03d}.trees"
		stats.append(singer_merge_windows(chrom, info, sources, dest, validated=True))
		files.append(dest)
	smap = merged / f"chr{chrom}.sample_map.tsv"
	with smap.open("w") as f:
		f.write("tree_node_id\tsample\thaplotype\n")
		for i, s in enumerate(info["samples"]):
			for h in range(info["ploidy"]):
				f.write(f"{i * info['ploidy'] + h}\t{s}\t{h + 1}\n")
	mask = merged / f"chr{chrom}.inferred.bed"
	mask.write_text("".join(f"chr{chrom}\t{w['start']}\t{w['end']}\n" for w in info["windows"]))
	singer_save(
		root / "qc.json",
		dict(
			stats=stats,
			regional=info["regional"],
			retained=a.singer_samples,
			burnin=a.singer_burnin,
			thin=a.singer_thin,
			convergence="not assessed: inspect MCMC logs before scientific use",
		),
	)
	singer_record(done, req, files + [smap, mask, root / "qc.json"])
	# The validated chromosome posterior is now the durable scientific output.
	for w in info["windows"]:
		wd = root / f"w{w['start']}_{w['end']}"
		for path in wd.glob("posterior_*.trees"):
			path.unlink()
		singer_remove_native_text(wd)
		(wd / "raw.complete.json").unlink(missing_ok=True)


def singer_export(a, chrom):
	dest = a.output / "trace/singer"
	dest.mkdir(parents=True, exist_ok=True)
	for p in dest.glob(f"chr{chrom}.sample*.trees"):
		p.unlink()
	files = list((a.output / "posterior").glob(f"chr{chrom}.*"))
	for p in files:
		q = dest / p.name
		if q.is_symlink() or q.exists():
			q.unlink()
		q.symlink_to(p.resolve())
	marker = dest / "ARG_TRACE_BUILD.tsv"
	text = (
		"method\tsinger\nformat\ttrace\ncoordinate_system\t0-based chromosome bp\ntime_units\tgenerations\n"
	)
	if not marker.exists() or marker.read_text() != text:
		marker.write_text(text)
	singer_record(
		dest / f"chr{chrom}.complete.json",
		dict(native=singer_stamp(a.output / "native" / f"chr{chrom}" / "complete.json")),
		[dest / p.name for p in files] + [dest / "ARG_TRACE_BUILD.tsv"],
	)


def singer_check(a, chrom):
	p = a.output / "native" / f"chr{chrom}" / "complete.json"
	if not singer_valid(p):
		raise ValueError(f"incomplete or changed SINGER outputs: {p}")
	d = json.loads(p.read_text())
	req = d["request"]
	for key in ["prepared", "code", "binary", "converter"]:
		if singer_stamp(req[key][0]) != req[key]:
			raise ValueError(f"SINGER {key} changed; rebuild required")
	prep = Path(req["prepared"][0])
	data = json.loads(prep.read_text())
	if not singer_valid(prep):
		raise ValueError("prepared input changed")
	for key in ["source", "index", "panel", "keep"]:
		row = data["request"].get(key)
		if row and singer_stamp(row[0]) != row:
			raise ValueError(f"source {key} changed")
	tsfiles = [Path(x[0]) for x in d["outputs"] if x[0].endswith(".trees")]
	if len(tsfiles) != req["samples"]:
		raise ValueError("posterior count mismatch")
	for f in tsfiles:
		ts = tskit.load(f)
		if ts.time_units != "generations" or not ts.num_sites or not ts.num_mutations:
			raise ValueError("invalid tree")
	if a.format == "trace" and not singer_valid(a.output / "trace/singer" / f"chr{chrom}.complete.json"):
		raise ValueError("TRACE export missing or changed")
	print(f"[SINGER] CHECK PASS chr{chrom}: {len(tsfiles)} posterior samples", flush=True)


def singer_arg_singer_main():
	def interrupted(signum, frame):
		raise KeyboardInterrupt(f"signal {signum}")

	for sig in (signal.SIGTERM, signal.SIGHUP):
		signal.signal(sig, interrupted)
	a = singer_parse()
	home = Path(os.environ.get("GU_SINGER_HOME", str(Path.home() / "software/singer-0.1.9")))
	if a.action == "check":
		for chrom in a.chrs:
			singer_check(a, chrom)
		return
	if not (home / "singer").is_file():
		raise ValueError("SINGER missing; run bash install.sh --singer")
	locks = []
	for path in (a.prep, a.output):
		path.mkdir(parents=True, exist_ok=True)
		lock = (path / ".refgen.lock").open("a")
		fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
		locks.append(lock)
	for chrom in a.chrs:
		info = singer_prepare(a, chrom)
		if a.action == "prepare":
			continue
		singer_infer(a, chrom, info, home)
		if a.format == "trace":
			singer_export(a, chrom)
		singer_check(a, chrom)
	print("[SINGER] COMPLETE", flush=True)


def singer_arg_singer_cli():
	try:
		singer_arg_singer_main()
	except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as e:
		print("ERROR: " + str(e), file=sys.stderr)
		sys.exit(2)


sys.dont_write_bytecode = True


# 🚩 arg_threads_worker
"""Launch the installed Threads CLI with a bounded, private Ray runtime."""
import importlib.metadata


def threads_threads_worker_main():
	if sys.argv[1:2] == ["infer"]:
		# Ray otherwise rewrites even an explicit loopback IP to an external IP.
		# Despite its name, this switch controls that behavior on Linux too.
		os.environ["RAY_ENABLE_WINDOWS_OR_OSX_CLUSTER"] = "0"
		import ray

		original_init = ray.init

		def bounded_init(*args, **kwargs):
			# This pipeline is single-host. Avoid WSL/VPN interface changes breaking
			# Ray heartbeats, and never attach to a cluster from RAY_ADDRESS.
			kwargs["address"] = "local"
			kwargs["_node_ip_address"] = "127.0.0.1"
			kwargs.setdefault("num_cpus", int(os.environ["REFGEN_THREADS"]))
			kwargs.setdefault("object_store_memory", 512 * 1024**2)
			kwargs.setdefault("include_dashboard", False)
			return original_init(*args, **kwargs)

		ray.init = bounded_init

	(entry,) = [
		e
		for e in importlib.metadata.distribution("threads_arg").entry_points
		if e.group == "console_scripts" and e.name == "threads"
	]
	entry.load()()


def threads_arg_threads_worker_cli():
	threads_threads_worker_main()


# 🚩 arg_threads
"""PGEN -> Threads -> chromosome-coordinate, mutation-bearing TRACE trees.

Runs in an isolated environment. Stage manifests bind checkpoints to inputs,
options, software versions and output stamps; incomplete stages are never reused.
"""
import contextlib
from functools import wraps


threads_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(threads_HERE))


class threads_StageTimings:
	def __init__(self, **metadata):
		self.started = time.perf_counter()
		self.data = dict(
			metadata,
			**{s + "_seconds": 0.0 for s in ("prepare", "infer", "convert", "check")},
			cache_hits={},
			peak_rss_mb=None,
			peak_rss_scope="sampled sum of process-tree RSS; shared pages may be counted more than once",
		)
		self.stack = []
		self.stop = threading.Event()
		self.sampler = threading.Thread(target=self._sample_memory, daemon=True)
		self.sampler.start()

	def _sample_memory(self):
		try:
			page_size = os.sysconf("SC_PAGE_SIZE")
			while not self.stop.is_set():
				total = 0
				pending, seen = [os.getpid()], set()
				while pending:
					pid = pending.pop()
					if pid in seen:
						continue
					seen.add(pid)
					proc = Path("/proc") / str(pid)
					try:
						total += int((proc / "statm").read_text().split()[1]) * page_size
						for task in (proc / "task").iterdir():
							try:
								pending.extend(map(int, (task / "children").read_text().split()))
							except (OSError, ValueError):
								pass
					except (OSError, ValueError, IndexError):
						pass
				self.data["peak_rss_mb"] = max(self.data["peak_rss_mb"] or 0, total / 2**20)
				self.stop.wait(0.25)
		except (OSError, ValueError):
			pass  # Profiling must not change inference success/failure.

	@contextlib.contextmanager
	def stage(self, name):
		frame = [time.perf_counter(), 0.0]
		self.stack.append(frame)
		try:
			yield
		finally:
			elapsed = time.perf_counter() - frame[0]
			self.stack.pop()
			self.data[name + "_seconds"] += elapsed - frame[1]
			if self.stack:
				self.stack[-1][1] += elapsed

	def finish(self, path, status, **metadata):
		self.data.update(metadata, status=status, total_seconds=time.perf_counter() - self.started)
		self.stop.set()
		self.sampler.join(timeout=2)
		path = Path(path)
		path.parent.mkdir(parents=True, exist_ok=True)
		tmp = path.with_name(path.name + f".{os.getpid()}.next")
		try:
			tmp.write_text(json.dumps(self.data, indent=2, sort_keys=True) + "\n")
			tmp.replace(path)
		finally:
			tmp.unlink(missing_ok=True)


threads_TIMING = None


def threads_stage(name):
	return threads_TIMING.stage(name) if threads_TIMING is not None else contextlib.nullcontext()


def threads_timed_stage(name):
	def decorate(function):
		@wraps(function)
		def wrapped(*args, **kwargs):
			with threads_stage(name):
				return function(*args, **kwargs)

		return wrapped

	return decorate


def threads_cache_hit(name):
	if threads_TIMING is not None:
		threads_TIMING.data["cache_hits"][name] = True


def threads_env(key, default=None):
	value = os.environ.get("REFGEN_" + key, default)
	if value is None:
		raise ValueError(f"Missing REFGEN_{key}")
	return value


def threads_stamp(path):
	p = Path(path).resolve(strict=True)
	s = p.stat()
	return [str(p), s.st_size, s.st_mtime_ns]


def threads_atomic_json(path, data):
	path = Path(path)
	tmp = path.with_suffix(path.suffix + ".next")
	tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
	tmp.replace(path)


def threads_checkpoint(path, request, files):
	threads_atomic_json(path, dict(request=request, outputs=[threads_stamp(p) for p in files]))


def threads_verified(path, request=None):
	try:
		data = json.loads(Path(path).read_text())
		return (request is None or data["request"] == request) and all(
			threads_stamp(row[0]) == row for row in data["outputs"]
		)
	except (OSError, ValueError, KeyError, TypeError):
		return False


def threads_run(args, log):
	import shlex

	args = list(map(str, args))
	print("[THREADS] " + shlex.join(args), flush=True)
	with Path(log).open("a") as handle:
		handle.write("\nCOMMAND " + shlex.join(args) + "\n")
		handle.flush()
		proc = subprocess.Popen(args, stdout=handle, stderr=subprocess.STDOUT)
		try:
			code = proc.wait()
			if code:
				with Path(log).open("rb") as tail:
					tail.seek(max(0, tail.seek(0, 2) - 5000))
					print(tail.read().decode(errors="replace"), file=sys.stderr)
				raise subprocess.CalledProcessError(code, args)
		finally:
			if proc.poll() is None:
				proc.terminate()
				try:
					proc.wait(timeout=10)
				except subprocess.TimeoutExpired:
					proc.kill()
					proc.wait()


def threads_map_for(chrom):
	pattern = threads_env("MAP_PATTERN", "")
	if pattern:
		path = Path(pattern.replace("{chr}", chrom).replace("%CHR%", chrom))
		if path.is_file():
			return path
	root = Path(threads_env("MAP_DIR"))
	build = threads_env("GRCH")
	for name in (
		f"chr{chrom}.map",
		f"chr{chrom}.txt",
		f"chr{chrom}.map.gz",
		f"chr{chrom}.b{build}.gmap.gz",
		f"genetic_map_GRCh{build}_chr{chrom}.txt",
		f"genetic_map_chr{chrom}_combined_b{build}.txt",
		f"genetic_map_chr{chrom}_combined_b{build}.txt.map",
	):
		for suffix in ("", ".gz"):
			if (root / (name + suffix)).is_file():
				return root / (name + suffix)
	raise ValueError(f"No genetic map for chr{chrom} in {root}; use --map-pattern")


def threads_read_psam(path):
	with Path(path).open() as h:
		fields = h.readline().lstrip("#").split()
		rows = [dict(zip(fields, line.split())) for line in h if line.strip()]
	names = [r["IID"] for r in rows]
	if len(names) != len(set(names)):
		raise ValueError("Duplicate sample IID")
	return rows


def threads_repack(source, output, chrom, source_map):
	"""Stream PGEN; enforce phase and unique positions, pack only real X copies."""
	from importlib.util import module_from_spec, spec_from_file_location

	grid = Path(os.environ.get("REFGEN_GRID_ROOT", threads_HERE.parents[1] / "grid")) / "f"
	_spec = spec_from_file_location("grid_arg_helpers", grid / "0.arg.py")
	_grid_arg = module_from_spec(_spec)
	_spec.loader.exec_module(_grid_arg)
	read_source = _grid_arg.read_source
	samples = threads_read_psam(str(source) + ".psam")
	ids, columns = [], []
	for i, row in enumerate(samples):
		if chrom == "X" and row.get("SEX") not in ("1", "2"):
			raise ValueError("Non-PAR X requires known sex for every selected sample")
		copies = 1 if chrom == "X" and row["SEX"] == "1" else 2
		for hap in range(copies):
			ids.append((row["IID"], hap + 1))
			columns.append(2 * i + hap)
	if len(ids) < 4 or len(ids) % 2:
		raise ValueError("Threads needs an even number of >=4 real haplotypes; check X panel policy")
	variants = []
	with Path(str(source) + ".pvar").open() as h:
		for line in h:
			if line.startswith("##"):
				continue
			if line.startswith("#"):
				fields = line.lstrip("#").split()
				continue
			variants.append(dict(zip(fields, line.split())))
	bp, cm = read_source(source_map)
	selected, positions, seen = [], [], set()
	previous = -1
	for i, row in enumerate(variants):
		pos = int(row["POS"])
		if pos < previous:
			raise ValueError("Unsorted PVAR")
		previous = pos
		if pos not in seen:
			selected.append(i)
			positions.append(pos)
			seen.add(pos)
	if len(selected) < 2:
		raise ValueError("Need at least two distinct polymorphic sites")
	cols = np.asarray(columns, dtype=np.int32)
	alleles = np.empty(2 * len(samples), dtype=np.int32)
	phase = np.empty(len(samples), dtype=np.uint8)
	with pgenlib.PgenReader(os.fsencode(str(source) + ".pgen")) as reader:
		with pgenlib.PgenWriter(
			os.fsencode(str(output) + ".pgen"),
			len(ids) // 2,
			variant_ct=len(selected),
			hardcall_phase_present=True,
		) as writer:
			for index in selected:
				reader.read_alleles_and_phasepresent(index, alleles, phase)
				if np.any((alleles < 0) | (alleles > 1)):
					raise ValueError(
						f"Missing/non-biallelic genotype at variant {index}; use --missing-max 0"
					)
				het = alleles[::2] != alleles[1::2]
				if np.any(het & (phase == 0)):
					raise ValueError(f"Unphased heterozygote at variant {index}")
				if chrom == "X" and any(het[i] for i, row in enumerate(samples) if row["SEX"] == "1"):
					raise ValueError("Heterozygous male non-PAR X genotype")
				writer.append_alleles(np.ascontiguousarray(alleles[cols]), all_phased=True)
	with Path(str(output) + ".pvar").open("w") as h:
		h.write("#CHROM\tPOS\tID\tREF\tALT\n")
		for i in selected:
			r = variants[i]
			h.write(f"{chrom}\t{r['POS']}\t{r['ID']}\t{r['REF']}\t{r['ALT']}\n")
	with Path(str(output) + ".psam").open("w") as h:
		h.write("#IID\tSEX\n")
		for i in range(len(ids) // 2):
			h.write(f"pair_{i}\t2\n")
	with Path(str(output) + ".haplotypes.tsv").open("w") as h:
		h.write("sample\thaplotype\n")
		for sample, hap in ids:
			h.write(f"{sample}\t{hap}\n")
	with Path(str(output) + ".map").open("w") as h:
		h.write("pos\tchr\tcM\n")
		for pos, genetic in zip(positions, np.interp(positions, bp, cm)):
			h.write(f"{pos}\t{chrom}\t{genetic:.12g}\n")
	threads_atomic_json(
		str(output) + ".qc.json",
		dict(
			chromosome=chrom,
			haplotypes=len(ids),
			individuals=len(samples),
			variants=len(selected),
			duplicate_positions_removed=len(variants) - len(selected),
			errors=[],
		),
	)


@threads_timed_stage("prepare")
def threads_prepare(chrom, prep, work, grid, replace, require=False):
	p = prep / ("chr" + chrom)
	source = pfile_prefix(threads_env("PFILE_DIR"), chrom)
	if chrom == "X":
		male_sample_ids(str(source) + ".psam")
	var = Path(str(source) + ".pvar")
	if not var.is_file():
		var = Path(str(source) + ".pvar.zst")
	source_map = threads_map_for(chrom)
	inputs = [
		str(source) + ".pgen",
		var,
		str(source) + ".psam",
		source_map,
		threads_env("ANCESTRY_FILE"),
		__file__,
		threads_HERE / "arg.py",
		grid / "0.arg.py",
		grid / "0.arg.py",
		shutil.which("plink2"),
	]
	for option in ("ARG_KEEP", "EXTRACT"):
		if threads_env(option, ""):
			inputs.append(threads_env(option))
	if chrom == "X":
		inputs.append(threads_HERE / "arg.py")
	request = dict(
		inputs=[threads_stamp(x) for x in inputs],
		options={
			k: threads_env(k)
			for k in (
				"GRCH",
				"ARG_MAF_MIN",
				"ARG_GENO_MAX",
				"SEED",
				"MAX_INDIVIDUALS",
				"ANCHORS_PER_POP",
				"ARG_FULL",
				"X_ODD_MALE",
			)
		},
	)
	marker = Path(str(p) + ".prepared.json")
	if not replace and threads_verified(marker, request):
		threads_cache_hit("prepare")
		print(f"SKIP verified Threads PGEN chr{chrom}", flush=True)
		return p
	if require:
		raise ValueError(f"Prepared chr{chrom} missing/stale; run prep_gen with the same options")
	marker.unlink(missing_ok=True)
	log = prep / "log" / f"arg.prepare.chr{chrom}.log"
	with tempfile.TemporaryDirectory(prefix=f"prepare.chr{chrom}.", dir=work) as tmp:
		tmp = Path(tmp)
		keep = tmp / "keep.tsv"
		order = tmp / "order.tsv"
		cmd = [
			sys.executable,
			grid / "0.arg.py",
			"make-arg-keep",
			"--sample",
			str(source) + ".psam",
			"--ancestry",
			threads_env("ANCESTRY_FILE"),
			"--out-keep",
			keep,
			"--out-panel",
			tmp / "panel.tsv",
			"--seed",
			threads_env("SEED"),
			"--anchors-per-pop",
			threads_env("ANCHORS_PER_POP"),
			"--max-individuals",
			threads_env("MAX_INDIVIDUALS"),
		]
		if threads_env("ARG_FULL") == "TRUE":
			cmd += ["--full"]
		if threads_env("ARG_KEEP", ""):
			cmd += ["--custom-keep", threads_env("ARG_KEEP")]
		threads_run(cmd, log)
		if chrom == "X":
			needle_x_panel(keep, str(source) + ".psam", threads_env("X_ODD_MALE"), tmp / "selection.qc.json")
			shutil.copy2(tmp / "selection.qc.json", str(p) + ".selection.qc.json")
		order.write_text("\n".join(keep.read_text().splitlines()[1:]) + "\n")
		subset = tmp / "subset"
		cmd = ["plink2", "--pfile", source]
		if str(var).endswith(".zst"):
			cmd += ["vzs"]
		cmd += [
			"--keep",
			keep,
			"--indiv-sort",
			"f",
			order,
			"--snps-only",
			"just-acgt",
			"--max-alleles",
			"2",
			"--maf",
			threads_env("ARG_MAF_MIN"),
			"--mac",
			"1",
			"--geno",
			threads_env("ARG_GENO_MAX"),
			"--sort-vars",
			"--make-pgen",
			"--threads",
			threads_env("THREADS"),
			"--out",
			subset,
		]
		if chrom == "X":
			lo, hi = {"37": (2699521, 154931043), "38": (2781480, 155701382)}[threads_env("GRCH")]
			cmd += ["--chr", "X", "--from-bp", str(lo), "--to-bp", str(hi)]
		if threads_env("EXTRACT", ""):
			cmd += ["--extract", threads_env("EXTRACT")]
		threads_run(cmd, log)
		threads_repack(subset, tmp / "ready", chrom, source_map)
		extensions = ("pgen", "pvar", "psam", "haplotypes.tsv", "map", "qc.json")
		try:
			for ext in extensions:
				shutil.copy2(tmp / f"ready.{ext}", str(p) + "." + ext + ".next")
				Path(str(p) + "." + ext + ".next").replace(str(p) + "." + ext)
		finally:
			for ext in extensions:
				Path(str(p) + "." + ext + ".next").unlink(missing_ok=True)
		shutil.copy2(keep, str(p) + ".keep.tsv")
	files = [str(p) + "." + ext for ext in (*extensions, "keep.tsv")]
	if chrom == "X":
		files.append(str(p) + ".selection.qc.json")
	threads_checkpoint(marker, request, files)
	print(f"PREPARED Threads chr{chrom}: {p}", flush=True)
	return p


def threads_export_tree(raw, dest, identities, sample_map, pgen_prefix):
	ts = tszip.decompress(raw)
	with Path(identities).open() as h:
		ids = [(r["sample"], int(r["haplotype"])) for r in csv.DictReader(h, delimiter="\t")]
	if ts.num_samples != len(ids):
		raise ValueError("Threads conversion changed sample count")
	tables = ts.dump_tables()
	offset = ts.metadata.get("offset", 0) if isinstance(ts.metadata, dict) else 0
	if offset:
		tables.sequence_length += offset
		tables.edges.left += offset
		tables.edges.right += offset
		tables.sites.position += offset
	tables.time_units = "generations"
	tables.metadata_schema = tskit.MetadataSchema.permissive_json()
	tables.metadata = dict(
		source_method="threads", source_arg_offset=offset, offset=0, coordinate_system="chromosome_bp"
	)
	tables.provenances.add_row(
		record=json.dumps(
			dict(
				software={"name": "gu-threads"},
				parameters={
					"source_offset": offset,
					"time_units": "generations",
					"ancestral_state_source": "tree_parsimony_not_external_AA",
				},
			)
		)
	)
	# arg-needle-lib 1.1.3's mutation converter assumes infinite sites; Threads
	# parsimony can have several mutations at one site. Build one site per PVAR
	# row and preserve parent links for recurrent/back mutations instead.
	topology = tables.tree_sequence()
	marginal = topology.first()
	genotypes = np.empty(len(ids), dtype=np.int32)
	with pgenlib.PgenReader(os.fsencode(str(pgen_prefix) + ".pgen")) as reader:
		with Path(str(pgen_prefix) + ".pvar").open() as variants:
			index = 0
			for line in variants:
				if line.startswith("#"):
					continue
				chrom, position, vid, ref, alt = line.split()[:5]
				position = int(position)
				if not 0 <= position < topology.sequence_length:
					raise ValueError(f"Variant outside inferred ARG: {position}")
				reader.read_alleles(index, genotypes)
				marginal.seek(position)
				ancestor, mutations = marginal.map_mutations(genotypes, alleles=[ref, alt])
				site = tables.sites.add_row(position=position, ancestral_state=ancestor)
				base = len(tables.mutations)
				for mutation in mutations:
					parent = tskit.NULL if mutation.parent == tskit.NULL else base + mutation.parent
					tables.mutations.add_row(
						site=site, node=mutation.node, derived_state=mutation.derived_state, parent=parent
					)
				index += 1
				if index % 100000 == 0:
					print(f"CONVERT mapped {index} sites: {dest}", flush=True)
	tables.sort()
	result = tables.tree_sequence()
	result.dump(str(dest) + ".next")
	Path(str(dest) + ".next").replace(dest)
	Path(sample_map).write_text(
		"tree_node_id\tsample\thaplotype\n"
		+ "".join(f"{int(node)}\t{name}\t{hap}\n" for node, (name, hap) in zip(result.samples(), ids))
	)
	return result, [(int(node), name, hap) for node, (name, hap) in zip(result.samples(), ids)]


def threads_build(chrom, p, out, work, model, replace, action):
	try:
		with tempfile.TemporaryDirectory(prefix=f"build.chr{chrom}.", dir=work) as scratch:
			return threads__build(chrom, p, out, Path(scratch), model, replace, action)
	finally:
		if out is not None:
			for pending in (
				out / "threads" / f"chr{chrom}.threads.next",
				out / "trees" / f"chr{chrom}.trees.next",
			):
				pending.unlink(missing_ok=True)


def threads__build(chrom, p, out, work, model, replace, action):
	with threads_stage("infer"):
		instruction = out / "threads" / f"chr{chrom}.threads"
		marker = out / "threads" / f"chr{chrom}.infer.json"
		request = dict(
			prepared=threads_stamp(str(p) + ".prepared.json"),
			model=model,
			mode=threads_env("THREADS_MODE"),
			fit=threads_env("THREADS_FIT"),
			mutation_rate=threads_env("THREADS_MUTATION_RATE"),
			threads=threads_env("THREADS"),
			code=[threads_stamp(__file__), threads_stamp(threads_HERE / "arg.py")],
			versions={
				k: importlib.metadata.version(k)
				for k in ("threads_arg", "arg-needle-lib", "numpy", "pgenlib", "tskit")
			},
		)
		if replace or not threads_verified(marker, request):
			if action == "convert":
				raise ValueError(f"No compatible inference checkpoint for chr{chrom}")
			marker.unlink(missing_ok=True)
			# Native Linux copies for repeated PGEN scans by all worker processes.
			local = work / f"chr{chrom}"
			for ext in ("pgen", "pvar", "psam", "map"):
				shutil.copyfile(str(p) + "." + ext, str(local) + "." + ext)
			cmd = [
				sys.executable,
				threads_HERE / "arg.py",
				"threads",
				"worker",
				"infer",
				"--pgen",
				str(local) + ".pgen",
				"--map",
				str(local) + ".map",
				"--demography",
				model["path"],
				"--mutation_rate",
				threads_env("THREADS_MUTATION_RATE"),
				"--mode",
				threads_env("THREADS_MODE"),
				"--num_threads",
				threads_env("THREADS"),
				"--out",
				str(instruction) + ".next",
			]
			if threads_env("THREADS_FIT") == "TRUE":
				cmd += ["--fit_to_data"]
			threads_run(cmd, out / "log" / f"arg.infer.chr{chrom}.log")
			Path(str(instruction) + ".next").replace(instruction)
			threads_checkpoint(marker, request, [instruction])
		else:
			threads_cache_hit("infer")
			print(f"SKIP verified Threads inference chr{chrom}", flush=True)
	if action == "infer":
		return
	with threads_stage("convert"):
		tree = out / "trees" / f"chr{chrom}.trees"
		smap = tree.with_suffix(".sample_map.tsv")
		converted = out / "trees" / f"chr{chrom}.converted.json"
		conversion = dict(
			inference=threads_stamp(marker),
			identities=threads_stamp(str(p) + ".haplotypes.tsv"),
			code=threads_stamp(__file__),
		)
		if replace or not threads_verified(converted, conversion):
			converted.unlink(missing_ok=True)
			raw = work / f"chr{chrom}.tsz"
			threads_run(
				["threads", "convert", "--threads", instruction, "--tsz", raw],
				out / "log" / f"arg.convert.chr{chrom}.log",
			)
			ts, rows = threads_export_tree(raw, tree, str(p) + ".haplotypes.tsv", smap, p)
			with threads_stage("check"):
				stats = trace_validate(ts, rows, chrom)
				if ts.num_mutations < 1:
					raise ValueError("Threads tree has no mutations")
			threads_atomic_json(tree.with_suffix(".arg_qc.json"), stats)
			threads_checkpoint(converted, conversion, [tree, smap, tree.with_suffix(".arg_qc.json")])
			raw.unlink()
		else:
			threads_cache_hit("convert")
		if threads_env("ARG_FORMAT") == "trace":
			threads_run(
				[
					sys.executable,
					threads_HERE / "arg.py",
					"trace",
					"--method",
					"threads",
					"--tree",
					tree,
					"--sample-map",
					smap,
					"--out-dir",
					out / "trace/threads",
					"--chr",
					chrom,
				],
				out / "log" / f"arg.trace.chr{chrom}.log",
			)


@threads_timed_stage("check")
def threads_check(chrom, out):
	tree = out / "trees" / f"chr{chrom}.trees"
	smap = tree.with_suffix(".sample_map.tsv")
	if chrom == "X":
		validate_x_map(smap, str(pfile_prefix(threads_env("PFILE_DIR"), chrom)) + ".psam")
	if not threads_verified(out / "threads" / f"chr{chrom}.infer.json") or not threads_verified(
		out / "trees" / f"chr{chrom}.converted.json"
	):
		raise ValueError(f"Missing/modified Threads checkpoint chr{chrom}")
	converted = json.loads((out / "trees" / f"chr{chrom}.converted.json").read_text())
	if converted["request"]["inference"] != threads_stamp(out / "threads" / f"chr{chrom}.infer.json"):
		raise ValueError(f"Tree does not match current Threads inference chr{chrom}")

	ts = tskit.load(tree)
	stats = trace_validate(ts, trace_read_map(smap, "threads"), chrom)
	if ts.time_units != "generations" or ts.num_mutations < 1:
		raise ValueError("Invalid Threads time units or missing mutations")
	if threads_env("ARG_FORMAT") == "trace":
		# Read-only validation: no logs or manifests are written here.
		subprocess.run(
			[
				sys.executable,
				str(threads_HERE / "arg.py"),
				"trace",
				"--method",
				"threads",
				"--tree",
				str(tree),
				"--sample-map",
				str(smap),
				"--out-dir",
				str(out / "trace/threads"),
				"--chr",
				chrom,
				"--check",
			],
			check=True,
		)
	print(f"Threads CHECK PASS chr{chrom}: {stats}", flush=True)


def threads_write_build_metadata(out, chrs, model):
	"""A chrX-only rebuild must retain the completed autosomes in the manifest."""
	path = out / "ARG_THREADS_BUILD.json"
	previous = json.loads(path.read_text()).get("chromosomes", []) if path.is_file() else []
	completed = []
	for chrom in dict.fromkeys(previous + chrs):
		inference = out / "threads" / f"chr{chrom}.infer.json"
		converted = out / "trees" / f"chr{chrom}.converted.json"
		if (
			threads_verified(inference)
			and threads_verified(converted)
			and json.loads(converted.read_text())["request"]["inference"] == threads_stamp(inference)
		):
			completed.append(chrom)
	threads_atomic_json(path, dict(method="threads", chromosomes=completed, model=model))
	if threads_env("ARG_FORMAT") == "trace":
		(out / "trace/threads/ARG_TRACE_BUILD.tsv").write_text(
			f"method\tthreads\nformat\ttrace\nbuild\tGRCh{threads_env('GRCH')}\nchromosomes\t{' '.join(completed)}\n"
		)


def threads_arg_threads_main():
	global threads_TIMING

	def stopped(signum, frame):
		raise SystemExit(128 + signum)

	for sig in (signal.SIGTERM, signal.SIGHUP):
		signal.signal(sig, stopped)
	out = Path(threads_env("ARG_DIR"))
	prep = Path(threads_env("GEN4ARG_DIR"))
	action = threads_env("ARG_ACTION")
	chrs = threads_env("ARG_CHRS").split()
	if action == "check":
		for chrom in chrs:
			threads_check(chrom, out)
		return
	grid = Path(os.environ.get("REFGEN_GRID_ROOT", threads_HERE.parents[1] / "grid")) / "f"
	sys.path.insert(0, str(grid))
	if not shutil.which("plink2"):
		raise ValueError("plink2 is missing from PATH")
	for chrom in chrs:
		source = pfile_prefix(threads_env("PFILE_DIR"), chrom)
		for ext in ("pgen", "psam"):
			threads_stamp(str(source) + "." + ext)
		if not any(Path(str(source) + "." + ext).is_file() for ext in ("pvar", "pvar.zst")):
			raise ValueError(f"Missing PVAR for chr{chrom}")
		if chrom == "X":
			male_sample_ids(str(source) + ".psam")
		threads_map_for(chrom)
	threads_stamp(threads_env("ANCESTRY_FILE"))
	work = Path(threads_env("ARG_SCRATCH")) / "refgen_threads"
	if str(work).startswith("/mnt/"):
		raise ValueError("--scratch must use native Linux storage")
	for path in (work, prep / "log", out / "threads", out / "trees", out / "log"):
		path.mkdir(parents=True, exist_ok=True)
	# Ray uses Unix sockets, whose paths must fit the ~108-byte OS limit.
	ray_tmp = tempfile.TemporaryDirectory(prefix="gu-ray-")
	os.environ["RAY_TMPDIR"] = ray_tmp.name
	import fcntl

	lock = (work / ".lock").open("a")
	fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
	replace = threads_env("REPLACE") == "TRUE"
	mu = float(threads_env("THREADS_MUTATION_RATE"))
	if not math.isfinite(mu) or mu <= 0:
		raise ValueError("--mutation-rate must be finite and positive")
	custom = threads_env("THREADS_DEMOGRAPHY", "")
	if custom:
		demo = Path(custom).resolve()
	else:
		demo = work / "constant.haploid20000.demo"
		if not demo.exists():
			demo.write_text("0\t20000\n")
	values = np.loadtxt(demo, ndmin=2)
	if (
		values.shape[1] != 2
		or not np.isfinite(values).all()
		or values[0, 0] != 0
		or np.any(np.diff(values[:, 0]) <= 0)
		or np.any(values[:, 1] <= 0)
	):
		raise ValueError(
			"Demography must start at generation 0 with increasing times and positive haploid Ne"
		)
	model = dict(
		path=str(demo),
		sha256=hashlib.sha256(demo.read_bytes()).hexdigest(),
		description="custom haploid Ne" if custom else "constant haploid Ne=20000 (diploid Ne=10000)",
	)
	print(f"Threads model: {model}; workers={threads_env('THREADS')}; chromosome jobs=1", flush=True)
	for chrom in chrs:
		threads_TIMING = threads_StageTimings(
			method="threads",
			chromosome=chrom,
			action=action,
			threads=int(threads_env("THREADS")),
			mode=threads_env("THREADS_MODE"),
			fit_to_data=threads_env("THREADS_FIT") == "TRUE",
			samples=None,
			sites=None,
			chromosome_jobs=1,
			started_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
		)
		status, error = "failed", None
		try:
			p = threads_prepare(chrom, prep, work, grid, replace, require=action in ("infer", "convert"))
			qc = json.loads(Path(str(p) + ".qc.json").read_text())
			threads_TIMING.data.update(
				samples=qc["individuals"], sites=qc["variants"], haplotypes=qc["haplotypes"]
			)
			if action != "prepare":
				threads_build(chrom, p, out, work, model, replace, action)
			if action not in ("prepare", "infer"):
				threads_check(chrom, out)
			status = "passed"
		except BaseException as exc:
			error = str(exc)
			raise
		finally:
			try:
				threads_TIMING.finish(out / "log" / f"arg.timing.chr{chrom}.json", status, error=error)
			except OSError as exc:
				print(f"WARNING: could not write Threads timing: {exc}", file=sys.stderr)
			finally:
				threads_TIMING = None
	if action not in ("prepare", "infer"):
		threads_write_build_metadata(out, chrs, model)


def threads_arg_threads_cli():
	try:
		threads_arg_threads_main()
	except Exception as error:
		print(f"ERROR: {error}", file=sys.stderr, flush=True)
		raise


def singer_entry():
	global np, pysam, tskit
	import numpy as np
	import pysam
	import tskit

	return singer_arg_singer_cli()


def threads_entry():
	if sys.argv[1:2] == ["worker"]:
		sys.argv.pop(1)
		return threads_threads_worker_main()
	if sys.argv[1:2] in (["--help"], ["-h"]):
		print(
			"Threads driver: use ./arg.sh build|check --method threads; worker: arg.py threads worker --help"
		)
		return
	global np, pgenlib, tskit, tszip
	import numpy as np
	import pgenlib
	import tskit
	import tszip

	return threads_arg_threads_cli()


arg_COMMANDS.update(singer=singer_entry, threads=threads_entry)

if __name__ == "__main__":
	arg_main()
