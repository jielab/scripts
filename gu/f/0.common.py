#!/usr/bin/env python3
"""GU 0.common utilities. See --help for commands."""

from __future__ import annotations
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
sys.dont_write_bytecode = True


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
"""Readers and writers for temporary GU tables with wide fields."""


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


# 🚩 Persistent native results and RDS exports
import hashlib
import shutil
import subprocess
import fcntl
from contextlib import ExitStack


RESULT_NAMES = {"phyml": "phyml.haplotypes.xlsx", "ibdmix": "ibdmix.tracts.xlsx", "trace": "trace.segments.xlsx", "as3": "as3.tracts.xlsx"}
LEGACY_RESULT_NAMES = {method: str(Path(name).with_suffix(".rds")) for method, name in RESULT_NAMES.items()}


def gu_result_digest(path):
	digest = hashlib.sha256()
	with Path(path).open('rb') as handle:
		for block in iter(lambda: handle.read(1024 * 1024), b''):
			digest.update(block)
	return digest.digest()


def gu_work_root(published):
	root = Path(published).resolve()
	for temporary in ('/tmp', '/var/tmp', '/dev/shm', '/run'):
		if root == Path(temporary) or Path(temporary) in root.parents:
			raise ValueError(f'Main GU results require a persistent directory, not {root}')
	return root


def gu_result_lock(run):
 directory = Path('/tmp/gu-locks')
 directory.mkdir(parents=True, exist_ok=True)
 return directory / (hashlib.sha256(str(Path(run).resolve()).encode()).hexdigest()[:16] + '.lock')


def gu_atomic_copy(source, destination):
	"""Verify a same-directory staged copy before replacing a durable result."""
	destination = Path(destination)
	destination.parent.mkdir(parents = True, exist_ok = True)
	fd, name = tempfile.mkstemp(prefix = '.' + destination.name + '.part.', dir = destination.parent)
	os.close(fd)
	stage = Path(name)
	try:
		shutil.copy2(source, stage)
		if gu_result_digest(source) != gu_result_digest(stage):
			raise IOError(f'Result publication verification failed: {destination}')
		with stage.open('rb') as handle:
			os.fsync(handle.fileno())
		os.replace(stage, destination)
	finally:
		stage.unlink(missing_ok = True)


def gu_result_r(*arguments):
	environment = os.environ.copy()
	for name in ('R_HOME', 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE', 'R_ENVIRON_USER'):
		environment.pop(name, None)
	with tempfile.TemporaryDirectory(prefix = 'gu-results-code-', dir = '/tmp') as directory:
		stage = Path(directory)
		for source in (Path(__file__).with_suffix('.R'), Path(__file__).with_name('phyml.R'), Path(__file__).parents[2] / '0f/results.R'):
			shutil.copy2(source, stage / source.name)
		environment['GU_RESULTS_R'] = str(stage / 'results.R')
		environment['GU_PHYML_R'] = str(stage / 'phyml.R')
		paired_r = Path(sys.executable).with_name('Rscript')
		rscript = environment.get('GU_RESULTS_RSCRIPT') or (str(paired_r) if paired_r.is_file() else shutil.which('Rscript'))
		if not rscript:
			raise RuntimeError('Rscript is required to recover or export GU RDS results')
		subprocess.run([rscript, '--vanilla', str(stage / '0.common.R'), *map(str, arguments)], env = environment, check = True)


def gu_run_specs(work, output, methods, run_only = None):
	specs = []
	for method in methods:
		if method not in RESULT_NAMES:
			continue
		root = work / method
		if not root.is_dir():
			continue
		for target in ([run_only.parent] if run_only is not None else sorted(root.iterdir())):
			if not target.is_dir() or target.name in ('tmp', 'log', 'inputs'):
				continue
			for run in ([run_only] if run_only is not None else sorted(target.iterdir())):
				if not run.is_dir() or run.name in ('tmp', 'log', 'inputs'):
					continue
				files, figures = [], []
				for directory, folders, names in os.walk(run, followlinks = False):
					folders[:] = [name for name in folders if name not in ('log', 'logs', 'tmp', 'mask', 'genotype', '__pycache__')]
					for name in folders:
						path = Path(directory) / name
						if path.is_symlink():
							files.append({'name': path.relative_to(run).as_posix(), 'link': os.readlink(path)})
					for name in names:
						path = Path(directory) / name
						relative = path.relative_to(run).as_posix()
						# Exports live beside native results. Never embed an archive
						# in itself or reinterpret already-published figure names.
						if path.parent == run and (name.startswith(method + '.')):
							continue
						if name.startswith(('.published-', '.result-')) or (path.suffix in ('.log', '.err', '.cmd', '.lock', '.pdf', '.xlsx') and not name.endswith('.phyml.log')) or '.tmp' in name or '.part.' in name:
							continue
						if path.suffix == '.png':
							if method != 'phyml' or '.panelB' not in name:
								raise ValueError(f'Add the numerical table mapping for figure: {path}')
							stem = name.split('_phyml_tree.', 1)[1].replace('.panelB', '')
							prefix = '' if Path(relative).parent.as_posix() == 'loci' else Path(relative).parent.name + '.'
							info = path.stat()
							figures.append({'source': relative, 'name': f'phyml.{prefix}{stem}', 'modified': info.st_mtime, 'size': info.st_size})
							continue
						if path.is_symlink():
							files.append({'name': relative, 'link': os.readlink(path)})
						else:
							info = path.stat()
							files.append({'name': relative, 'modified': info.st_mtime, 'size': info.st_size})
				if not files:
					continue
				destination = output / run.relative_to(work) / RESULT_NAMES[method]
				specs.append({'method': method, 'source': str(run), 'source_root': str(work), 'destination': str(destination), 'files': files, 'figures': figures})
	return specs


def gu_review_specs(work, output):
	directory = work / 'final/review'
	files = []
	for name in (
		'phyml_locus_report.tsv', 'phyml_haplotype_report.tsv', 'phyml_copy_validation.tsv', 'phyml_lineage_trees.tsv',
		'loci_review.tsv', 'population_concordance.tsv', 'candidate_copy_support.tsv', 'strict_support_summary.tsv',
		'selected_ibdmix_calls.tsv', 'locus_tag_summary.tsv', 'lead_comparison.tsv', 'haplotype_tag_candidates.tsv',
		'tag_population_metrics.tsv', 'best_tag_by_population.tsv', 'best_haplotypes.tsv', 'gwas_lookup.tsv',
	):
		path = directory / name
		if path.is_file():
			info = path.stat()
			files.append({'name': name, 'modified': info.st_mtime, 'size': info.st_size})
	return [{'method': 'validation', 'source': str(directory), 'source_root': str(work),
		'destination': str(output / 'final/gu.validation.xlsx'), 'files': files, 'figures': []}] if files else []


def gu_publish_results(work, output, methods, run_only = None):
	work, output = Path(work).resolve(), gu_work_root(output)
	with ExitStack() as locks, tempfile.TemporaryDirectory(prefix = 'gu-publish-', dir = '/tmp') as temporary:
		stage = Path(temporary)
		specs = gu_run_specs(work, stage, methods, run_only)
		if 'final' in methods:
			specs.extend(gu_review_specs(work, stage))
		pending = []
		for spec in specs:
			lock_path = gu_result_lock(spec['source'])
			if lock_path.is_file():
				lock = locks.enter_context(lock_path.open('r'))
				try:
					fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
				except BlockingIOError:
					print('Result publication deferred for active run:', spec['source'], file = sys.stderr)
					continue
			marker = Path(spec['source']) / '.published-content'
			signature = hashlib.sha256(json.dumps([str(output), spec['files'], spec['figures']], sort_keys = True).encode()).hexdigest()
			existing = output / Path(spec['destination']).relative_to(stage)
			if existing.is_file() and marker.is_file() and marker.read_text() == signature:
				continue
			spec['signature'] = signature
			pending.append(spec)
		specs = pending
		specification = stage / 'jobs.json'
		specification.write_text(json.dumps(specs))
		gu_result_r('pack', specification)
		for spec in specs:
			for figure in spec['figures']:
				destination = Path(spec['destination']).parent / figure['name']
				shutil.copy2(Path(spec['source']) / figure['source'], destination)
				if not destination.with_suffix('.xlsx').is_file():
					raise ValueError(f'Missing numerical results: {destination}')
		if 'final' in methods and (work / 'final/gu.sqlite').is_file():
			gu_result_r('database-pack', work / 'final/gu.sqlite', stage / 'final/gu.results.rds')
			gu_publish_summary(work, stage)
		for source in stage.rglob('*'):
			if not source.is_file() or source == specification:
				continue
			destination = output / source.relative_to(stage)
			gu_atomic_copy(source, destination)
		for spec in specs:
			(Path(spec['source']) / '.published-content').write_text(spec['signature'])


def gu_publish_summary(work, stage):
	# Final summary tables contain substantial overlap with the participant RDS.
	# Publish only independent aggregate evidence needed for inspection.
	module = spec_from_file_location('gu_shared_results', Path(__file__).parents[2] / '0f/results.py')
	results = module_from_spec(module)
	module.loader.exec_module(results)
	import pandas as pd
	root = work / 'final/normalize/summary'
	validation = {}
	for name in ('phyml_locus_report', 'phyml_haplotype_report'):
		path = work / 'final/review' / (name + '.tsv')
		if not path.is_file():
			continue
		frame = pd.read_csv(path, sep = '\t', low_memory = False, float_precision = 'round_trip')
		if frame.empty:
			continue
		if name == 'phyml_haplotype_report':
			columns = ['locus_id', 'hap_id', 'lineage', 'role', 'n_copies', 'n_individuals', 'archaic', 'prop_match',
				'candidate_start', 'candidate_end', 'call', 'superpopulation_copy_counts', 'ibdmix_status',
				'ibdmix_supported_individuals', 'ibdmix_support_fraction', 'trace_status', 'trace_supported_individuals', 'trace_support_fraction']
			frame = frame[[column for column in columns if column in frame]]
		validation[name.removeprefix('phyml_').removesuffix('_report')] = frame
	if validation:
		results.write_workbook(validation, stage / 'final/gu.validation.summary.xlsx')
	groups = {'gu.loci': ['locus_evidence', 'locus_method_support'], 'gu.segments': ['segment_catalog'], 'gu.trajectory': ['locus_trajectory']}
	for name, names in groups.items():
		tables = {}
		for table in names:
			path = root / (table + '.tsv.gz')
			if path.is_file() and path.stat().st_size:
				frame = pd.read_csv(path, sep = '\t', low_memory = False, float_precision = 'round_trip')
				frame = frame.loc[:, ~frame.columns.str.contains(r'(?:file|path)$')]
				if not frame.empty:
					tables[table.removeprefix('locus_')] = frame
		if tables:
			if name == 'gu.trajectory':
				for chromosome, frame in tables['trajectory'].groupby('chr', sort = False):
					results.write_workbook({'trajectory': frame}, stage / 'final' / f'gu.trajectory.chr{chromosome}.xlsx')
			else:
				results.write_workbook(tables, stage / 'final' / (name + '.xlsx'))


def gu_restore_results(published, method="all", run_only=None):
	published = work = gu_work_root(published)
	work.mkdir(parents = True, exist_ok = True, mode = 0o700)
	with (work / '.restore.lock').open('a') as lock:
		fcntl.flock(lock, fcntl.LOCK_EX)
		jobs = []
		markers = []
		for archive_method in RESULT_NAMES:
			if method not in ("all", "final", archive_method):continue
			archives = [run_only / (archive_method + '.raw.tar.gz')] if run_only else sorted((work / archive_method).glob('*/*/*.raw.tar.gz'))
			for archive in archives:
				if not archive.is_file():continue
				marker = archive.parent / '.native-restored'
				if not marker.exists():
					gu_extract_native(archive,archive.parent,archive.parent,work)
					marker.touch()
		for legacy_method, name in LEGACY_RESULT_NAMES.items():
			for source in sorted((published / legacy_method).glob('*/*/' + name)):
				if run_only and source.parent != run_only:continue
				if (source.parent / (legacy_method + '.raw.tar.gz')).is_file():continue
				stat = source.stat()
				key = f'{stat.st_size}:{stat.st_mtime_ns}'
				target = work / source.parent.relative_to(published)
				marker = target / '.published-result'
				# Native files are authoritative, including partially completed
				# analyses. RDS is a recovery/export format, never a rollback.
				if not (target / '.result-restoring').exists() and any((target / name).exists() for name in ('run.meta.tsv', 'loci', 'final', 'results', 'extract')):
					continue
				jobs.append({'source': str(source), 'target': str(target)})
				markers.append((marker, key))
		review = published / 'final/gu.validation.rds'
		if review.is_file():
			stat = review.stat(); key = f'{stat.st_size}:{stat.st_mtime_ns}'
			target = work / 'final/review'; marker = target / '.published-result'
			if (target / '.result-restoring').exists() or not (target / 'phyml_locus_report.tsv').is_file():
				jobs.append({'source': str(review), 'target': str(target)})
				markers.append((marker, key))
		with tempfile.NamedTemporaryFile(mode = 'w', suffix = '.json', dir = '/tmp') as handle:
			json.dump(jobs, handle); handle.flush()
			if jobs:
				for job in jobs:
					target = Path(job['target'])
					target.mkdir(parents = True, exist_ok = True)
					(target / '.result-restoring').write_text(job['source'])
				gu_result_r('restore', handle.name)
		for job in jobs:
			target = Path(job['target'])
			origin = target / '.result-source-root'
			gu_rebase_paths(target, Path(origin.read_text().strip()) if origin.is_file() else published, work)
			(target / '.result-restoring').unlink(missing_ok = True)

		for marker, key in markers:
			marker.write_text(key)
		source = published / 'final/gu.results.rds'
		if source.is_file():
			stat = source.stat(); key = f'{stat.st_size}:{stat.st_mtime_ns}'
			marker = work / 'final/.published-database'
			if not (work / 'final/gu.sqlite').is_file():
				gu_result_r('database-restore', source, work / 'final/gu.sqlite')
				marker.write_text(key)
		if method == "shiny":gu_prepare_read_view(work)
		else:gu_link_phyml(work)
	return work


def gu_rebase_paths(target, published, work):
	if str(published) == str(work):
		return
	old, new = str(published).encode(), str(work).encode()
	for directory, folders, names in os.walk(target, followlinks = False):
		for name in [*names, *[name for name in folders if (Path(directory) / name).is_symlink()]]:
			path = Path(directory) / name
			if name in ('.result-source-root', 'run.meta.tsv', 'cache.meta.tsv'):
				continue
			if path.is_symlink():
				link = os.readlink(path)
				if link.startswith(str(published) + '/'):
					path.unlink(); path.symlink_to(str(work) + link[len(str(published)):])
			elif path.suffix in ('.json', '.tsv', '.txt', '.cmd', '.list'):
				data = path.read_bytes()
				if old in data:
					stat = path.stat(); path.write_bytes(data.replace(old, new)); os.utime(path, ns = (stat.st_atime_ns, stat.st_mtime_ns))


def gu_link_phyml(work, target=None):
	directory = work / 'final/normalize'
	directory.mkdir(parents = True, exist_ok = True)
	link = directory / 'phyml'
	target = target or work / 'phyml'
	if link.is_symlink() and os.readlink(link) != str(target):link.unlink()
	if not link.exists() and not link.is_symlink():
		link.symlink_to(target, target_is_directory = True)


# 🚩 Compact native results: durable archives, disposable read copies
import tarfile
import io


def gu_native_files(run, method):
	for directory, folders, names in os.walk(run, followlinks=False):
		folders[:] = [n for n in folders if n not in ('tmp', 'mask', 'log', 'logs', '__pycache__')]
		for name in names:
			p = Path(directory) / name
			if p.parent == run and (name.startswith(method + '.') or name.startswith(('.published-', '.result-', '.native-'))):
				continue
			if p.suffix in ('.lock', '.cmd', '.err', '.pdf') or (p.suffix == '.log' and not name.endswith('.phyml.log')) or '.part.' in name:
				continue
			yield p


def gu_native_view_root(work):
	return Path('/tmp/gu-native-view') / hashlib.sha256(str(work).encode()).hexdigest()[:16]


def gu_extract_native(archive, destination, run, work, replace=False):
	with tarfile.open(archive, 'r:gz') as bundle:
		manifest = json.load(bundle.extractfile('GU-MANIFEST.json'))
		for row in manifest['files']:
			rel = Path(row['path'])
			if rel.is_absolute() or '..' in rel.parts:
				raise ValueError('Unsafe native result path')
			p = destination / rel
			if p.exists() or p.is_symlink():
				if not replace:continue
				p.unlink()
			p.parent.mkdir(parents=True, exist_ok=True)
			if 'link' in row:
				link = row['link']
				if destination != run:
					absolute = os.path.abspath(run / rel.parent / link)
					if absolute.startswith(str(run) + '/'):
						link = str(destination) + absolute[len(str(run)):]
					elif absolute.startswith(str(work) + '/'):
						link = str(gu_native_view_root(work)) + absolute[len(str(work)):]
				p.symlink_to(link)
			else:
				h = hashlib.sha256()
				fd,temporary = tempfile.mkstemp(prefix='.'+p.name+'.part.',dir=p.parent)
				try:
					with bundle.extractfile(row['path']) as source, os.fdopen(fd,'wb') as out:
						for block in iter(lambda:source.read(4*1024*1024), b''):
							out.write(block);h.update(block)
					if h.hexdigest() != row['sha256']:raise ValueError('Native archive checksum mismatch: '+str(p))
					os.replace(temporary,p)
					os.utime(p, ns=(row['mtime_ns'],row['mtime_ns']))
				finally:Path(temporary).unlink(missing_ok=True)
	return manifest


def gu_run_read_view(work, run, method):
	view = gu_native_view_root(work)
	archive = run / (method + '.raw.tar.gz')
	link = view / run.relative_to(work)
	link.parent.mkdir(parents=True, exist_ok=True)
	target = run
	if archive.is_file():
		info = archive.stat()
		key = hashlib.sha256(f'{archive}:{info.st_size}:{info.st_mtime_ns}'.encode()).hexdigest()[:24]
		target = view / '.objects' / key
		if not (target / '.gu-view-ready').is_file():
			stage = Path(tempfile.mkdtemp(prefix='.extract-', dir=view))
			try:
				gu_extract_native(archive,stage,run,work)
				(stage / '.gu-view-ready').touch()
				target.parent.mkdir(parents=True,exist_ok=True)
				if target.exists():shutil.rmtree(target)
				os.replace(stage,target)
				# Internal links were rebased to the staging path before rename.
				gu_rebase_paths(target,stage,target)
			finally:
				if stage.exists():shutil.rmtree(stage)
	if not link.is_symlink() or os.readlink(link) != str(target):
		staged = link.with_name('.'+link.name+'.link')
		staged.unlink(missing_ok=True);staged.symlink_to(target, target_is_directory=True)
		os.replace(staged,link)
	return target


def gu_prepare_read_view(work):
	for method in RESULT_NAMES:
		for run in sorted((work / method).glob('*/*')):
			if run.is_dir() and run.name not in ('inputs','log','tmp'):
				gu_run_read_view(work,run,method)
	gu_link_phyml(work,gu_native_view_root(work) / 'phyml')


def gu_archive_run(work, run, method):
	archive = run / (method + '.raw.tar.gz')
	files = list(gu_native_files(run,method))
	# A compact run has no materialized algorithm outputs; keep its archive.
	if archive.is_file() and not any((run / n).is_dir() for n in ('final','loci','extract','results','raw')):
		return archive
	if not files:return None
	manifest = {'format':'gu-native-v1','method':method,'files':[]}
	with tempfile.TemporaryDirectory(prefix='gu-archive-',dir='/tmp') as directory:
		stage = Path(directory) / archive.name
		with tarfile.open(stage,'w:gz',compresslevel=1) as bundle:
			for p in files:
				name = p.relative_to(run).as_posix()
				if p.is_symlink():row = {'path':name,'link':os.readlink(p)}
				else:
					info = p.stat()
					row = {'path':name,'sha256':gu_result_digest(p).hex(),'mtime_ns':info.st_mtime_ns,'bytes':info.st_size}
				bundle.add(p,arcname=name,recursive=False)
				manifest['files'].append(row)
			data = json.dumps(manifest).encode()
			entry = tarfile.TarInfo('GU-MANIFEST.json');entry.size=len(data)
			bundle.addfile(entry,io.BytesIO(data))
		# Every native file must pass independent read-back before publication.
		check = Path(directory) / 'verified'
		gu_extract_native(stage,check,run,work)
		gu_atomic_copy(stage,archive)
	return archive


def _gu_compact_results(work, methods, run_only=None):
	work = gu_work_root(work)
	# Initial view links allow the running Shiny to read each untouched run.
	gu_prepare_read_view(work)
	for method in methods:
		if method not in RESULT_NAMES:continue
		for run in ([run_only] if run_only is not None else sorted((work / method).glob('*/*'))):
			if not run.is_dir() or run.name in ('inputs','log','tmp'):continue
			lock_path = gu_result_lock(run)
			with lock_path.open('a') as lock:
				try:fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
				except BlockingIOError:continue
				archive = gu_archive_run(work,run,method)
				if archive is None:continue
				gu_run_read_view(work,run,method)
				# Archive and read copy both verified. Keep only user exports and
				# the small caller roster/provenance required by density summaries.
				for p in list(run.iterdir()):
					if p.name.startswith(method + '.') or p.name in ('run.meta.tsv','cache.meta.tsv','.published-content') or '.part.' in p.name:continue
					if p.name == 'samples' and method == 'ibdmix':
						for q in list(p.rglob('*')):
							if q.is_file() and q.name not in ('ALL.txt','male.txt'):q.unlink()
						continue
					if p.is_dir() and not p.is_symlink():
						if p.name in ('mask','log','logs','tmp'):
							cache = Path('/tmp/gu-intermediate') / hashlib.sha256(str(run).encode()).hexdigest()[:16] / p.name
							cache.parent.mkdir(parents=True,exist_ok=True)
							if cache.exists():shutil.rmtree(cache)
							shutil.move(str(p),cache)
						else:shutil.rmtree(p)
					else:p.unlink()
				print('COMPACT',method,run.name,flush=True)
	gu_link_phyml(work,gu_native_view_root(work) / 'phyml')


def gu_compact_results(work, methods, run_only=None):
 work = gu_work_root(work)
 with (work / '.restore.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  _gu_compact_results(work,methods,run_only)


def results_cli():
	parser = argparse.ArgumentParser()
	parser.add_argument('action', choices = ['work-root', 'restore', 'publish', 'compact'])
	parser.add_argument('--published', type = Path, required = True)
	parser.add_argument('--work', type = Path)
	parser.add_argument('--run', type = Path)
	parser.add_argument('--method', choices = [*RESULT_NAMES, 'final', 'all', 'shiny', 'ukb'], default = 'all')
	args = parser.parse_args()
	if args.action == 'work-root':
		print(gu_work_root(args.published))
	elif args.action == 'restore':
		print(gu_restore_results(args.published, args.method, args.run))
	else:
		methods = [*RESULT_NAMES, 'final'] if args.method == 'all' else [args.method]
		if args.action == 'compact':gu_compact_results(args.published,methods,args.run)
		else:gu_publish_results(args.work or gu_work_root(args.published), args.published, methods, args.run)


SUBCOMMANDS = {
	"results": results_cli,
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
