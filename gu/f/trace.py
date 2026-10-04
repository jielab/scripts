#!/usr/bin/env python3
"""GU trace utilities. See --help for commands."""

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


# 🚩 trace_output
"""Publish TRACE NPZ outputs only after successful execution and ZIP validation."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import zipfile
import zlib


def valid_npz(path):
	try:
		with zipfile.ZipFile(path) as archive:
			names = archive.namelist()
			return bool(names) and all(n.endswith(".npy") for n in names) and archive.testzip() is None
	except (OSError, ValueError, EOFError, RuntimeError, zipfile.BadZipFile, zlib.error):
		return False


def signature(paths):
	return [[str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in paths]


def reusable(paths, receipt):
	try:
		if json.loads(receipt.read_text()) == signature(paths):
			return True
	except (OSError, ValueError):
		pass
	# Adopt complete pre-existing outputs after validating every ZIP member.
	if not all(valid_npz(p) for p in paths):
		return False
	seal(paths, receipt)
	return True


def seal(paths, receipt):
	part = receipt.with_name(receipt.name + ".part")
	part.write_text(json.dumps(signature(paths)) + "\n")
	part.replace(receipt)


def run(prefix, suffixes, command):
	prefix = Path(prefix)
	paths = [Path(str(prefix) + s) for s in suffixes]
	receipt = Path(str(prefix) + ".complete.json")
	if reusable(paths, receipt):
		print(f"[GU TRACE] SKIP verified={prefix}", flush=True)
		return 0
	receipt.unlink(missing_ok=True)
	# Delete invalid legacy files so a failed replacement cannot be adopted.
	for p in paths:
		p.unlink(missing_ok=True)
	with tempfile.TemporaryDirectory(prefix=".trace-part-", dir=prefix.parent) as tmp:
		staged = Path(tmp) / prefix.name
		print(f"[GU TRACE] START output={prefix}", flush=True)
		proc = subprocess.Popen(command + ["-o", str(staged)], stdin=subprocess.DEVNULL, start_new_session=True)
		try:
			rc = proc.wait()
		finally:
			if proc.poll() is None:
				os.killpg(proc.pid, signal.SIGTERM)
				try:
					proc.wait(timeout=5)
				except subprocess.TimeoutExpired:
					os.killpg(proc.pid, signal.SIGKILL)
					proc.wait()
		if rc:
			return rc if rc > 0 else 128 - rc
		outputs = [Path(str(staged) + s) for s in suffixes]
		if not all(valid_npz(p) for p in outputs):
			print(f"ERROR: incomplete TRACE NPZ output: {prefix}", flush=True)
			return 1
		for source, target in zip(outputs, paths):
			source.replace(target)
		seal(paths, receipt)
		print(f"[GU TRACE] DONE output={prefix}", flush=True)
	return 0


def interrupted(signum, frame):
	raise KeyboardInterrupt


def trace_output_cli():
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--prefix", required=True)
	parser.add_argument("--suffix", action="append", required=True)
	parser.add_argument("command", nargs=argparse.REMAINDER)
	args = parser.parse_args()
	command = args.command[1:] if args.command[:1] == ["--"] else args.command
	if not command:
		parser.error("a TRACE command is required after --")
	signal.signal(signal.SIGTERM, interrupted)
	try:
		raise SystemExit(run(args.prefix, args.suffix, command))
	except KeyboardInterrupt:
		raise SystemExit(130)


# 🚩 trace_combine
"""Normalize TRACE calls and optionally apply requested loci post hoc.

TRACE inference should normally use chromosome-scale/multi-chromosome context.  The
loci BED is therefore a final segment-selection map, not an HMM training region.
"""

import argparse, re
from pathlib import Path
import pandas as pd


def first(cols, names):
	low = {str(c).lower(): c for c in cols}
	return next((low[n.lower()] for n in names if n.lower() in low), None)


def norm_one(path: Path, node: str) -> pd.DataFrame:
	try:
		d = pd.read_csv(path, sep="\t")
	except pd.errors.EmptyDataError:
		return pd.DataFrame()
	if d.empty:
		return d
	c = first(d.columns, ["chrom", "chromosome", "chr"])
	s = first(d.columns, ["start", "start_bp", "left", "begin"])
	e = first(d.columns, ["end", "end_bp", "right", "stop"])
	p = first(d.columns, ["mean_posterior", "posterior", "posterior_mean", "prob", "probability"])
	cm = first(d.columns, ["length(cM)", "length_cM", "length_cm", "genetic_length_cm"])
	if c is None or s is None or e is None:
		raise SystemExit(f"Unrecognized TRACE columns in {path}: {list(d.columns)}")
	out = pd.DataFrame(
		{
			"tree_node_id": str(node),
			"chr": d[c].astype(str).str.replace("chr", "", regex=False),
			"start": pd.to_numeric(d[s], errors="coerce"),
			"end": pd.to_numeric(d[e], errors="coerce"),
		}
	)
	out["posterior"] = pd.to_numeric(d[p], errors="coerce") if p else pd.NA
	out["length_cM"] = pd.to_numeric(d[cm], errors="coerce") if cm else pd.NA
	return out[(out.start.notna()) & (out.end > out.start)].copy()


def clip_loci(data: pd.DataFrame, path: Path) -> pd.DataFrame:
	loci = pd.read_csv(
		path, sep="\t", header=None, names=["chr", "locus_start", "locus_end", "locus_id"], dtype={"chr": str}
	)
	loci["chr"] = loci.chr.str.replace("chr", "", regex=False)
	joined = data.merge(loci, on="chr", how="inner")
	joined["start"] = joined[["start", "locus_start"]].max(axis=1)
	joined["end"] = joined[["end", "locus_end"]].min(axis=1)
	return joined[joined.end > joined.start].drop(columns=["locus_start", "locus_end"])


def trace_combine_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--root", required=True, type=Path)
	ap.add_argument("--sample-map", required=True, type=Path)
	ap.add_argument("--build", default="GRCh37")
	ap.add_argument("--loci", type=Path)
	args = ap.parse_args()
	sample_map = pd.read_csv(args.sample_map, sep="\t", dtype={"tree_node_id": str})
	pieces = []
	for path in sorted((args.root / "calls").glob("hap*.summary.txt")):
		match = re.match(r"hap(.+)\.summary\.txt$", path.name)
		if match:
			part = norm_one(path, match.group(1))
			pieces.extend([part] if not part.empty else [])
	cols = ["tree_node_id", "chr", "start", "end", "posterior", "length_cM"]
	hap = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=cols)
	if args.loci and not hap.empty:
		hap = clip_loci(hap, args.loci)
	hap["length_bp"] = hap.end - hap.start
	hap = hap.merge(sample_map[["tree_node_id", "sample", "haplotype"]], on="tree_node_id", how="left")
	if not hap.empty and hap["sample"].isna().any():
		missing = sorted(hap.loc[hap["sample"].isna(), "tree_node_id"].astype(str).unique())
		raise SystemExit(f"TRACE sample map is missing tree nodes: {','.join(missing[:10])}")
	hap["method"] = "trace"
	hap["source"] = "ghost_or_unknown_archaic"
	hap["genome_build"] = args.build
	order = [
		"sample",
		"haplotype",
		"tree_node_id",
		"method",
		"source",
		"chr",
		"start",
		"end",
		"length_bp",
		"posterior",
		"length_cM",
		"genome_build",
	]
	if "locus_id" in hap:
		order.insert(5, "locus_id")
	final = args.root / "final"
	final.mkdir(parents=True, exist_ok=True)
	hap[order].to_csv(final / "trace_haplotype_segments.tsv.gz", sep="\t", index=False, compression="gzip")
	print(f"TRACE segments: {len(hap)} -> {final / 'trace_haplotype_segments.tsv.gz'}")


def trace_combine_cli():
	trace_combine_main()


SUBCOMMANDS = {
	"output": trace_output_cli,
	"combine": trace_combine_cli,
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
