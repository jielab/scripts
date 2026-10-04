#!/usr/bin/env python3
"""GU 0.common utilities. See --help for commands."""

from __future__ import annotations
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


def load_module(filename):
	"""Load a neighbouring GU module, including filenames with dots."""
	key = "gu_" + Path(filename).stem.replace(".", "_")
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


# 🚩 comm
"""Compatibility helpers for GU TSV files with very wide fields."""


import csv
import gzip
import os
import sys
import tempfile
from pathlib import Path


def resolve_tsv_path(path: Path) -> Path:
	"""Prefer compressed unfiltered inventories, with legacy plain TSV fallback."""
	path = Path(path)
	if path.name.endswith("unfiltered.tsv.gz"):
		plain, compressed = path.with_suffix(""), path
	elif path.name.endswith("unfiltered.tsv"):
		plain, compressed = path, path.with_name(path.name + ".gz")
	else:
		return path
	if compressed.is_file():
		return compressed
	return plain if plain.is_file() else path


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
	path = resolve_tsv_path(path)
	if not path.is_file() or path.stat().st_size == 0:
		return []
	opener = gzip.open if path.suffix == ".gz" else open
	with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
		return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv_rows(path: Path, fieldnames, rows) -> None:
	"""Stream TSV output atomically; unfiltered inventories are always gzip."""
	path = Path(path)
	if path.name.endswith("unfiltered.tsv"):
		path = path.with_name(path.name + ".gz")
	path.parent.mkdir(parents=True, exist_ok=True)
	fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
	os.close(fd)
	tmp = Path(name)
	try:
		opener = gzip.open if path.suffix == ".gz" else open
		with opener(tmp, "wt", encoding="utf-8", newline="") as handle:
			fields = list(fieldnames)
			writer = csv.DictWriter(
				handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
			)
			writer.writeheader()
			for row in rows:
				writer.writerow({key: row.get(key, "") for key in fields})
		tmp.replace(path)
		if path.name.endswith("unfiltered.tsv.gz"):
			# Retire a previous plain output only after the new gzip is closed
			# and atomically published. Failed writes preserve both old forms.
			path.with_suffix("").unlink(missing_ok=True)
	finally:
		tmp.unlink(missing_ok=True)


def enable_wide_csv_fields() -> int:
	"""Raise the process-wide CSV field limit to the platform maximum."""
	limit = sys.maxsize
	while limit > 0:
		try:
			csv.field_size_limit(limit)
			return limit
		except OverflowError:
			# Some Python builds expose a C long smaller than sys.maxsize.
			limit //= 10
	raise RuntimeError("unable to configure a usable CSV field-size limit")


CHROM_LENGTHS = {
	"37": {
		"1": 249250621,
		"2": 243199373,
		"3": 198022430,
		"4": 191154276,
		"5": 180915260,
		"6": 171115067,
		"7": 159138663,
		"8": 146364022,
		"9": 141213431,
		"10": 135534747,
		"11": 135006516,
		"12": 133851895,
		"13": 115169878,
		"14": 107349540,
		"15": 102531392,
		"16": 90354753,
		"17": 81195210,
		"18": 78077248,
		"19": 59128983,
		"20": 63025520,
		"21": 48129895,
		"22": 51304566,
		"X": 155270560,
	},
	"38": {
		"1": 248956422,
		"2": 242193529,
		"3": 198295559,
		"4": 190214555,
		"5": 181538259,
		"6": 170805979,
		"7": 159345973,
		"8": 145138636,
		"9": 138394717,
		"10": 133797422,
		"11": 135086622,
		"12": 133275309,
		"13": 114364328,
		"14": 107043718,
		"15": 101991189,
		"16": 90338345,
		"17": 83257441,
		"18": 80373285,
		"19": 58617616,
		"20": 64444167,
		"21": 46709983,
		"22": 50818468,
		"X": 156040895,
	},
}

"""Expand normalized BED loci while retaining an auditable core/analysis map."""
import argparse
import re
from pathlib import Path


def parse_bp(value: str) -> int:
	m = re.fullmatch(r"\s*([0-9]+)\s*(bp|k|kb|m|mb)?\s*", value, re.I)
	if not m:
		raise argparse.ArgumentTypeError("use a non-negative size such as 0, 100kb, or 1mb")
	scale = {None: 1, "bp": 1, "k": 1000, "kb": 1000, "m": 1000000, "mb": 1000000}[
		m.group(2).lower() if m.group(2) else None
	]
	return int(m.group(1)) * scale


def expand_loci_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--input", required=True, type=Path)
	ap.add_argument("--output", required=True, type=Path)
	ap.add_argument("--map", required=True, type=Path)
	ap.add_argument("--flank", default="100kb", type=parse_bp)
	ap.add_argument("--build", required=True, choices=["37", "38"])
	args = ap.parse_args()
	rows = []
	for line_no, line in enumerate(args.input.read_text().splitlines(), 1):
		if not line.strip() or line.lstrip().startswith("#"):
			continue
		fields = line.split("\t")
		if len(fields) < 4:
			raise SystemExit(f"ERROR: normalized BED line {line_no} has fewer than four columns")
		chrom, start, end, locus = fields[:4]
		start = int(start)
		end = int(end)
		chrom_len = CHROM_LENGTHS[args.build][chrom]
		analysis_start = max(0, start - args.flank)
		analysis_end = min(chrom_len, end + args.flank)
		rows.append((chrom, start, end, analysis_start, analysis_end, locus, args.flank))
	if not rows:
		raise SystemExit("ERROR: no loci to expand")
	args.output.parent.mkdir(parents=True, exist_ok=True)
	args.map.parent.mkdir(parents=True, exist_ok=True)
	args.output.write_text("".join(f"{c}\t{a0}\t{a1}\t{locus}\n" for c, _, _, a0, a1, locus, _ in rows))
	args.map.write_text(
		"chr\tcore_start\tcore_end\tanalysis_start\tanalysis_end\tlocus_id\tflank_bp\n"
		+ "".join(f"{c}\t{s}\t{e}\t{a0}\t{a1}\t{locus}\t{flank}\n" for c, s, e, a0, a1, locus, flank in rows)
	)
	print(args.flank)


def comm_cli():
	if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
		print("Usage: 0.common.py expand-loci [options]")
	elif sys.argv.pop(1) == "expand-loci":
		expand_loci_main()
	else:
		raise SystemExit("Unknown shared utility command")


# 🚩 genetic_cache
"""Bind reusable genetic outputs to explicit source files and conversion options."""
import argparse
import json
from pathlib import Path


def stamp(name):
	p = Path(name).resolve(strict=True)
	s = p.stat()
	return [str(p), s.st_size, s.st_mtime_ns]


def genetic_cache_main():
	p = argparse.ArgumentParser()
	p.add_argument("action", choices=["check", "record"])
	p.add_argument("--receipt", type=Path, required=True)
	p.add_argument("--input", action="append", default=[])
	p.add_argument("--output", action="append", default=[])
	p.add_argument("--value", action="append", default=[])
	a = p.parse_args()
	try:
		data = dict(schema=1, inputs=[stamp(x) for x in a.input], outputs=[stamp(x) for x in a.output], values=a.value)
		if not all(x[1] > 0 for x in data["outputs"]):
			raise ValueError("empty output")
		if a.action == "check":
			return 0 if json.loads(a.receipt.read_text()) == data else 1
		temp = a.receipt.with_name(a.receipt.name + ".next")
		try:
			temp.write_text(json.dumps(data, indent=2) + "\n")
			temp.replace(a.receipt)
		finally:
			temp.unlink(missing_ok=True)
		return 0
	except (OSError, ValueError):
		if a.action == "check":
			return 1
		raise


def genetic_cache_cli():
	raise SystemExit(genetic_cache_main())


# 🚩 make_batches
import argparse, hashlib
from pathlib import Path


def stable(xs):
	return sorted(map(str, xs), key=lambda x: hashlib.sha256(x.encode()).hexdigest())


def make_batches_main():
	import pandas as pd

	ap = argparse.ArgumentParser(description="Deterministic ancestry-stratified target batches with fixed anchors")
	ap.add_argument("--panel", required=True, type=Path)
	ap.add_argument("--outdir", required=True, type=Path)
	ap.add_argument("--batch-size", type=int, default=1000)
	ap.add_argument("--sample-col", default="sample")
	ap.add_argument("--group-col", default="super_pop")
	ap.add_argument("--anchor-list", type=Path)
	ap.add_argument("--anchors-per-group", type=int, default=0)
	a = ap.parse_args()
	d = pd.read_csv(a.panel, sep=None, engine="python", dtype=str)
	if a.sample_col not in d:
		raise SystemExit(f"missing sample column {a.sample_col}; columns={list(d)}")
	if a.group_col not in d:
		d[a.group_col] = "ALL"
	d = d[[a.sample_col, a.group_col]].dropna(subset=[a.sample_col]).drop_duplicates(a.sample_col)
	a.outdir.mkdir(parents=True, exist_ok=True)
	anchors = []
	if a.anchor_list:
		anchors = [x.split()[0] for x in a.anchor_list.read_text().splitlines() if x.strip() and not x.startswith("#")]
	elif a.anchors_per_group > 0:
		for _, g in d.groupby(a.group_col, dropna=False):
			anchors += stable(g[a.sample_col])[: a.anchors_per_group]
		anchors = sorted(set(anchors))
		(a.outdir / "anchors.samples.txt").write_text("\n".join(anchors) + "\n")
	aset = set(anchors)
	rows = []
	for group, g in d.groupby(a.group_col, dropna=False):
		group = "NA" if pd.isna(group) else str(group)
		ids = [x for x in stable(g[a.sample_col]) if x not in aset]
		for k in range(0, len(ids), a.batch_size):
			target = ids[k : k + a.batch_size]
			bid = f"{group}.b{k // a.batch_size + 1:04d}"
			tf = a.outdir / f"{bid}.targets.txt"
			tf.write_text("\n".join(target) + "\n")
			jf = a.outdir / f"{bid}.joint.txt"
			jf.write_text("\n".join(sorted(set(anchors + target))) + "\n")
			rows.append([bid, group, len(target), len(anchors), str(tf), str(jf)])
	pd.DataFrame(rows, columns=["batch_id", "group", "n_target", "n_anchor", "target_list", "joint_list"]).to_csv(
		a.outdir / "batch_manifest.tsv", sep="\t", index=False
	)
	print(
		f"batches={len(rows)} target_samples={sum(r[2] for r in rows)} anchors={len(anchors)} manifest={a.outdir / 'batch_manifest.tsv'}"
	)


def make_batches_cli():
	make_batches_main()


# 🚩 make_sample_map
"""Create a TRACE tree-node -> sample/haplotype map from tree metadata."""

import argparse, json
from pathlib import Path


def load_ts(path: Path):
	import tskit

	if str(path).endswith(".tsz"):
		import tszip

		return tszip.decompress(str(path))
	return tskit.load(str(path))


def meta(obj):
	x = getattr(obj, "metadata", None)
	if isinstance(x, dict):
		return x
	if isinstance(x, (bytes, bytearray)):
		try:
			return json.loads(x.decode())
		except Exception:
			return {}
	return {}


def make_sample_map_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--tree", required=True, type=Path)
	ap.add_argument("--out", required=True, type=Path)
	args = ap.parse_args()
	ts = load_ts(args.tree)
	nodes = list(map(int, ts.samples()))
	rows = []
	unresolved = []
	seen = {}
	for j, node_id in enumerate(nodes):
		node = ts.node(node_id)
		cand = []
		for obj in [node] + ([ts.individual(node.individual)] if node.individual >= 0 else []):
			m = meta(obj)
			for k in ("variant_data_sample_id", "sample", "sample_id", "individual", "individual_id", "name", "id"):
				if k in m:
					cand.append(str(m[k]))
		sample = next((x for x in cand if x), None)
		if sample is None:
			unresolved.append((j, node_id))
		else:
			seen[sample] = seen.get(sample, 0) + 1
			rows.append([node_id, sample, seen[sample], "metadata"])
	if unresolved:
		examples = ",".join(str(node_id) for _, node_id in unresolved[:10])
		raise SystemExit(f"Cannot map {len(unresolved)} sample nodes from ARG metadata; example node IDs: {examples}")
	args.out.parent.mkdir(parents=True, exist_ok=True)
	with args.out.open("w") as h:
		h.write("tree_node_id\tsample\thaplotype\tmapping_mode\n")
		for r in rows:
			h.write("\t".join(map(str, r)) + "\n")
	print(f"Wrote {len(rows)} nodes to {args.out}; mode={rows[0][3] if rows else 'none'}")


def make_sample_map_cli():
	make_sample_map_main()


# 🚩 vcf_gt_fix
"""Stream a VCF while validating/fixing GT ploidy for X-chromosome units.

This does not invent phase. With --require-phased, unphased heterozygotes fail.
Haploid calls can be duplicated as homozygous diploid calls when a downstream
program requires two alleles per sample.
"""


import argparse
import sys


def vcf_gt_fix_main() -> None:
	p = argparse.ArgumentParser()
	p.add_argument("--chrom", help="replace CHROM with this value (e.g. 23)")
	p.add_argument("--duplicate-haploid", choices=["slash", "pipe"])
	p.add_argument("--require-diploid", action="store_true")
	p.add_argument("--require-phased", action="store_true")
	args = p.parse_args()

	n_records = n_calls = n_haploid = n_unphased_het = 0
	for line_no, raw in enumerate(sys.stdin, start=1):
		if raw.startswith("#"):
			sys.stdout.write(raw)
			continue
		fields = raw.rstrip("\n").split("\t")
		if len(fields) < 10:
			raise SystemExit(f"Malformed VCF line {line_no}: expected >=10 columns")
		n_records += 1
		if args.chrom:
			fields[0] = args.chrom
		fmt = fields[8].split(":")
		try:
			gt_idx = fmt.index("GT")
		except ValueError:
			raise SystemExit(f"VCF line {line_no} has no GT field")
		for i in range(9, len(fields)):
			vals = fields[i].split(":")
			if gt_idx >= len(vals):
				raise SystemExit(f"VCF line {line_no}, sample column {i + 1}: missing GT value")
			gt = vals[gt_idx]
			n_calls += 1
			if gt in {".", ""}:
				if args.duplicate_haploid:
					sep = "|" if args.duplicate_haploid == "pipe" else "/"
					vals[gt_idx] = f".{sep}."
				elif args.require_diploid:
					raise SystemExit(f"Haploid/missing-width GT {gt!r} at VCF line {line_no}")
			elif "/" not in gt and "|" not in gt:
				n_haploid += 1
				if args.duplicate_haploid:
					sep = "|" if args.duplicate_haploid == "pipe" else "/"
					vals[gt_idx] = f"{gt}{sep}{gt}"
				elif args.require_diploid:
					raise SystemExit(f"Haploid GT {gt!r} at VCF line {line_no}; duplication was not enabled")
			else:
				sep = "|" if "|" in gt else "/"
				alleles = gt.split(sep)
				if len(alleles) != 2:
					raise SystemExit(f"Non-diploid GT {gt!r} at VCF line {line_no}")
				if args.require_phased and sep == "/" and len(set(alleles) - {"."}) > 1:
					n_unphased_het += 1
					if n_unphased_het <= 5:
						print(f"UNPHASED_HET line={line_no} gt={gt}", file=sys.stderr)
			fields[i] = ":".join(vals)
		sys.stdout.write("\t".join(fields) + "\n")
	if args.require_phased and n_unphased_het:
		raise SystemExit(
			f"Found {n_unphased_het} unphased heterozygous GT calls; AS3 requires phased target haplotypes"
		)
	print(f"records={n_records} calls={n_calls} haploid_duplicated_or_seen={n_haploid}", file=sys.stderr)


def vcf_gt_fix_cli():
	vcf_gt_fix_main()


SUBCOMMANDS = {
	"genetic-cache": genetic_cache_cli,
	"batches": make_batches_cli,
	"sample-map": make_sample_map_cli,
	"fix-vcf-gt": vcf_gt_fix_cli,
	"expand-loci": expand_loci_main,
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
