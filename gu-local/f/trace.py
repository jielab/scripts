#!/usr/bin/env python3
"""GU trace utilities. See --help for commands."""

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
from gu_0_common import load_module


# 🚩 trace_output
"""Publish TRACE NPZ outputs only after successful execution and ZIP validation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import shutil
import subprocess
import tempfile
import zipfile
import zlib


def valid_npz(path):
	try:
		with zipfile.ZipFile(path) as archive:
			names = archive.namelist()
			if not (names and all(n.endswith(".npy") for n in names) and archive.testzip() is None):
				return False
		if str(path).endswith('.xss.npz'):
			import numpy as np
			with np.load(path, allow_pickle=False) as data:
				for key in ('gammas', 'params'):
					if key not in data or not np.isfinite(data[key]).all():
						return False
				g = data['gammas']
				if g.ndim != 2 or g.shape[0] != 2 or not np.all((g >= 0) & (g <= 1)) or not np.allclose(g.sum(axis=0), 1, atol=1e-7):
					return False
		return True
	except (OSError, ValueError, EOFError, RuntimeError, zipfile.BadZipFile, zlib.error):
		return False


def valid_output(path):
    if str(path).endswith('.npz'):
        return valid_npz(path)
    if str(path).endswith('.summary.txt'):
        import csv
        import math
        try:
            with Path(path).open() as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                required = {'chromosome', 'start', 'end', 'mean_posterior', 'length(bp)', 'length(cM)'}
                if not required <= set(reader.fieldnames or []):
                    return False
                for row in reader:
                    start,end,pp,length,cm = map(float, (row[k] for k in ('start','end','mean_posterior','length(bp)','length(cM)')))
                    if not all(map(math.isfinite,(start,end,pp,length,cm))) or not (0<=start<end and 0<=pp<=1 and length>0 and cm>=0):
                        return False
            return True
        except (OSError, ValueError, KeyError, TypeError):
            return False
    return False


def signature(paths):
	return [[str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in paths]


def receipt_data(paths):
	outputs = signature(paths)
	if any(str(p).endswith('.xss.npz') for p in paths):
		return dict(validation='finite-normalized-posterior-v1', outputs=outputs)
	return outputs


def reusable(paths, receipt):
	try:
		if json.loads(receipt.read_text()) == receipt_data(paths):
			return True
	except (OSError, ValueError):
		pass
	return False


def seal(paths, receipt):
	part = receipt.with_name(receipt.name + ".part")
	part.write_text(json.dumps(receipt_data(paths)) + "\n")
	part.replace(receipt)


def run(prefix, suffixes, command):
	prefix = Path(prefix)
	paths = [Path(str(prefix) + s) for s in suffixes]
	receipt = Path(str(prefix) + ".complete.json")
	if reusable(paths, receipt):
		print(f"[GU TRACE] SKIP verified={prefix}", flush=True)
		return 0
	receipt.unlink(missing_ok=True)
	# Remove unverified files before publishing a new completed result.
	for p in paths:
		p.unlink(missing_ok=True)
	with tempfile.TemporaryDirectory(prefix=".trace-part-", dir=prefix.parent) as tmp:
		staged = Path(tmp) / prefix.name
		command = list(command)
		if command and Path(command[0]).name == 'trace-summarize':
			# Upstream trace-summarize rewrites its input NPZs to append states.
			# Work on copies so interrupted summaries cannot corrupt verified
			# inference checkpoints (nor invalidate another concurrent reader).
			for option in ('-f', '--files'):
				if option in command:
					index = command.index(option) + 1
					copies = []
					for i, name in enumerate(command[index].split(',')):
						copy = Path(tmp) / f'input{i}.xss.npz'
						shutil.copy2(name, copy)
						copies.append(str(copy))
					command[index] = ','.join(copies)
		print(f"[GU TRACE] START output={prefix}", flush=True)
		guard = Path(__file__).resolve().parent / "process_lifetime.py"
		proc = subprocess.Popen([sys.executable, str(guard), "--parent-pid", str(os.getpid()), "--",
		                         *command, "-o", str(staged)], stdin=subprocess.DEVNULL, start_new_session=False)
		try:
			rc = proc.wait()
		finally:
			if proc.poll() is None:
				proc.terminate()
				try:
					proc.wait(timeout=10)
				except subprocess.TimeoutExpired:
					proc.kill()
					proc.wait()
		if rc:
			return rc if rc > 0 else 128 - rc
		outputs = [Path(str(staged) + s) for s in suffixes]
		if not all(valid_output(p) for p in outputs):
			print(f"ERROR: invalid TRACE output (missing/corrupt/non-finite or unnormalized posterior): {prefix}", flush=True)
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


def norm_one(path, node):
	import pandas as pd

	if not valid_output(path):
		raise SystemExit(f"Invalid current TRACE summary: {path}")
	d = pd.read_csv(path, sep="\t")
	out = pd.DataFrame({
		"tree_node_id": str(node),
		"chr": d["chromosome"].astype(str).str.replace("chr", "", regex=False),
		"start": pd.to_numeric(d["start"], errors="raise"),
		"end": pd.to_numeric(d["end"], errors="raise"),
		"posterior": pd.to_numeric(d["mean_posterior"], errors="raise"),
		"length_cM": pd.to_numeric(d["length(cM)"], errors="raise"),
	})
	return out


def clip_loci(data: pd.DataFrame, path: Path) -> pd.DataFrame:
	import pandas as pd

	loci = pd.read_csv(
		path, sep="\t", header=None, names=["chr", "locus_start", "locus_end", "locus_id"], dtype={"chr": str}
	)
	loci["chr"] = loci.chr.str.replace("chr", "", regex=False)
	joined = data.merge(loci, on="chr", how="inner")
	joined["start"] = joined[["start", "locus_start"]].max(axis=1)
	joined["end"] = joined[["end", "locus_end"]].min(axis=1)
	return joined[joined.end > joined.start].drop(columns=["locus_start", "locus_end"])


def trace_combine_main():
	import pandas as pd

	ap = argparse.ArgumentParser()
	ap.add_argument("--root", required=True, type=Path)
	ap.add_argument("--sample-map", required=True, type=Path)
	ap.add_argument("--build", default="GRCh37")
	ap.add_argument("--loci", type=Path)
	ap.add_argument("--filename", default="trace_haplotype_segments.tsv.gz")
	args = ap.parse_args()
	if Path(args.filename).name != args.filename:
		ap.error("--filename must be a filename without directories")
	sample_map = pd.read_csv(args.sample_map, sep="\t", dtype={"tree_node_id": str, "sample": str})
	if sample_map.tree_node_id.duplicated().any() or sample_map[["sample", "haplotype"]].duplicated().any():
		raise SystemExit("TRACE sample map contains duplicate nodes or biological haplotypes")
	wanted_file = args.root / "samples" / "tree_nodes.txt"
	wanted = set(wanted_file.read_text().split()) if wanted_file.is_file() else None
	pieces = []
	for path in sorted((args.root / "calls").glob("hap*.summary.txt")):
		match = re.match(r"hap(.+)\.summary\.txt$", path.name)
		if match:
			if wanted is not None and match.group(1) not in wanted:
				continue
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
	destination = final / args.filename
	stage = destination.with_name(destination.name + ".part")
	hap[order].to_csv(stage, sep="\t", index=False, compression="gzip")
	stage.replace(destination)
	print(f"TRACE segments: {len(hap)} -> {destination}")


def trace_combine_cli():
	trace_combine_main()


def validate_map_cli():
    """Validate authoritative ARG export maps against every posterior tree."""
    import csv
    ap = argparse.ArgumentParser()
    ap.add_argument('--tree', required=True, type=Path)
    ap.add_argument('--sample-map', required=True, type=Path)
    a = ap.parse_args()
    ts = _common.load_ts(a.tree)
    with a.sample_map.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    nodes = [int(r['tree_node_id']) for r in rows]
    identities = [(r['sample'], int(r['haplotype'])) for r in rows]
    if len(nodes) != len(set(nodes)) or set(nodes) != set(map(int, ts.samples())) or len(identities) != len(set(identities)):
        raise SystemExit(f'ERROR: sample map is not a bijection of tree samples: {a.tree}')
    by_node = dict(zip(nodes, identities))
    counts = {}
    keys = ('variant_data_sample_id', 'sample', 'sample_id', 'name')
    for node_id in ts.samples():
        node = ts.node(node_id)
        metadata = _common.meta(node)
        if node.individual >= 0:
            metadata = {**_common.meta(ts.individual(node.individual)), **metadata}
        name = next((str(metadata[k]) for k in keys if metadata.get(k) is not None), None)
        if name is not None:
            counts[name] = counts.get(name, 0) + 1
            hap = int(metadata.get('haplotype', counts[name]))
            if by_node[int(node_id)] != (name, hap):
                raise SystemExit(f'ERROR: ARG metadata disagrees with sample map for node {node_id}: {a.tree}')
    print(f'TRACE tree/map check PASS: {a.tree} samples={len(nodes)}')


def normalize_map_cli():
    """Convert the common 3/4-column genetic map layouts to TRACE's format."""
    import gzip
    import math
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--chr', required=True)
    a = ap.parse_args()
    opener = gzip.open if a.input.suffix == '.gz' else open
    rows = []
    with opener(a.input, 'rt') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) not in (3, 4):
                raise SystemExit(f'ERROR: expected position/rate/cM or chrom/position/rate/cM: {a.input}')
            values = fields[-3:]
            try:
                pos, rate, cm = map(float, values)
            except ValueError:
                if not rows and any(x.lower() in ('position', 'pos', 'bp', 'rate', 'genetic_map', 'map', 'gen_dist', 'genetic_position') or 'cm' in x.lower() or 'rate' in x.lower() for x in values):
                    continue
                raise SystemExit(f'ERROR: nonnumeric genetic map row: {line.strip()}')
            if len(fields) == 4 and fields[0].removeprefix('chr').upper() != a.chr.upper():
                raise SystemExit(f'ERROR: genetic map chromosome does not match chr{a.chr}: {a.input}')
            if not all(map(math.isfinite, (pos, rate, cm))) or min(pos, rate, cm) < 0 or not pos.is_integer():
                raise SystemExit(f'ERROR: invalid genetic map position/rate/distance: {a.input}')
            if rows and (pos <= rows[-1][0] or cm < rows[-1][2]):
                raise SystemExit(f'ERROR: map requires strictly increasing positions and nondecreasing cM: {a.input}')
            rows.append((int(pos), rate, cm))
    if len(rows) < 2:
        raise SystemExit(f'ERROR: genetic map needs at least two positions: {a.input}')
    a.output.parent.mkdir(parents=True, exist_ok=True)
    text = 'chrom\tpos\trate\tgen_dist\n' + ''.join(f'chr{a.chr}\t{p}\t{r:.12g}\t{c:.12g}\n' for p,r,c in rows)
    if not a.output.is_file() or a.output.read_text() != text:
        a.output.write_text(text)


# Shared CPU/CUDA inference uses the pinned official 2026 TRACE initializer.
RUNTIME_VERSION='gu-trace-2026-v1'


def trace_infer_main():
    backend=os.environ.get('TRACE_BACKEND','cpu')
    if backend=='cuda':
        sys.path.insert(0,str(Path(__file__).with_name('trace_gpu')))
        from backend import install
        install()
    elif backend!='cpu':
        raise RuntimeError(f'Unknown TRACE backend: {backend}')
    print(f'[GU TRACE] runtime={RUNTIME_VERSION} backend={backend}',flush=True)
    if sys.argv[1:]==['--gpu-check']:
        if backend!='cuda': raise RuntimeError('GPU check requires TRACE_BACKEND=cuda')
        return
    from tracehmm.infer_cli import main as cli
    cli()


SUBCOMMANDS = {
	"infer": trace_infer_main,
	"output": trace_output_cli,
	"combine": trace_combine_cli,
	"validate-map": validate_map_cli,
	"normalize-map": normalize_map_cli,
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
