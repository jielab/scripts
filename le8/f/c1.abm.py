# 🚩 c1.abm
"""C1 Python: agent-based modeling (ABM), public annotations and aggregate figures.

Subcommands: annotations, figures. Other arguments run ABM reference/TabICLv2.
Scientific dependencies load after lightweight commands and memory supervision.
"""

from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
import argparse
import csv
import urllib.request
import ast
import copy
import gzip
import hashlib
import importlib.metadata
import inspect
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import warnings

HERE = Path(__file__).resolve().parent


def read_ids(path):
	with open(path, encoding="utf-8") as handle:
		return list(dict.fromkeys(x.strip() for x in handle if x.strip()))


def run_gprofiler(argv):
	if len(argv) != 3:
		raise SystemExit(
			"usage: c1.abm.py annotations --gprofiler SIGNIFICANT_IDS BACKGROUND_IDS OUTPUT_CSV"
		)
	query, background = read_ids(argv[0]), read_ids(argv[1])
	payload = {
		"organism": "hsapiens",
		"query": query,
		"sources": ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
		# Return ranked terms even when none passes FDR 0.05. The R plot marks
		# how many displayed terms are significant instead of producing a
		# visually blank figure for a biologically null result.
		"user_threshold": 1.0,
		"significance_threshold_method": "fdr",
		"domain_scope": "custom",
		"background": background,
		"no_evidences": True,
	}
	request = urllib.request.Request(
		"https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
		data=json.dumps(payload).encode("utf-8"),
		headers={"Content-Type": "application/json", "User-Agent": "LE8-C1/1.0"},
		method="POST",
	)
	with urllib.request.urlopen(request, timeout=180) as response:
		result = json.load(response).get("result", [])
	fields = [
		"source",
		"term_id",
		"term_name",
		"adjusted_p",
		"intersection_size",
		"term_size",
		"query_size",
	]
	with open(argv[2], "w", newline="", encoding="utf-8") as handle:
		writer = csv.DictWriter(handle, fieldnames=fields)
		writer.writeheader()
		for row in result:
			if row.get("p_value") is not None:
				writer.writerow(
					{
						"source": row.get("source", ""),
						"term_id": row.get("native", ""),
						"term_name": row.get("name", ""),
						"adjusted_p": row.get("p_value", ""),
						"intersection_size": row.get("intersection_size", 0),
						"term_size": row.get("term_size", ""),
						"query_size": row.get("query_size", len(query)),
					}
				)


def run_annotations(argv):
	import csv, hashlib, io, json, math, sys, time, os, tempfile
	from pathlib import Path
	import requests

	src, out, cache = map(Path, argv)
	cache.mkdir(parents=True, exist_ok=True)
	r = list(csv.DictReader(src.open()))
	background = {x["gene"].upper() for x in r}
	sig = {x["gene"].upper() for x in r if x["selected"].upper() == "TRUE"}
	terms = []
	tf = []
	ppi = []
	audit = []

	def atomic_text(path, text):
		fd, tmp = tempfile.mkstemp(dir=path.parent)
		try:
			with os.fdopen(fd, "w") as f:
				f.write(text)
			os.replace(tmp, path)
		finally:
			if os.path.exists(tmp):
				os.unlink(tmp)

	def write(name, rows, cols):
		with (out / name).open("w", newline="") as f:
			w = csv.DictWriter(f, fieldnames=cols)
			w.writeheader()
			w.writerows(rows)

	def tail(k, m, N, n):
		def lc(a, b):
			return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

		logs = [
			lc(m, j) + lc(N - m, n - j) - lc(N, n)
			for j in range(max(k, n - (N - m)), min(m, n) + 1)
		]
		if not logs:
			return 1.0
		z = max(logs)
		return min(1.0, math.exp(z) * sum(math.exp(x - z) for x in logs))

	def bh(rows):
		rows.sort(key=lambda x: x["p_raw"])
		best = 1.0
		for i in range(len(rows) - 1, -1, -1):
			best = min(best, rows[i]["p_raw"] * len(rows) / (i + 1))
			rows[i]["adjusted_p"] = best
		return rows

	def library(label, names):
		for name in names:
			p = cache / (name + ".gmt")
			try:
				if not p.exists():
					z = requests.get(
						"https://maayanlab.cloud/Enrichr/geneSetLibrary",
						params={"mode": "text", "libraryName": name},
						timeout=45,
					)
					z.raise_for_status()
					if "\t" not in z.text:
						continue
					atomic_text(p, z.text)
				lines = p.read_text().splitlines()
				rows = []
				for line in lines:
					c = line.split("\t")
					if label == "TF" and not c[0].lower().endswith(" human"):
						continue
					g = {x.split(",")[0].upper() for x in c[2:]} & background
					if len(g) < 3:
						continue
					hits = g & sig
					k = len(hits)
					# Retain every eligible term for multiplicity, including zero overlaps.
					rows.append(
						dict(
							source=label,
							term_name=c[0],
							p_raw=(
								tail(k, len(g), len(background), len(sig)) if k else 1.0
							),
							intersection_size=k,
							genes=";".join(sorted(hits)),
						)
					)
				if not rows:
					continue
				audit.append(dict(source=label, status="ok", detail=name))
				return bh(rows)
			except Exception as e:
				audit.append(dict(source=label, status="unavailable", detail=str(e)))
		return []

	def mgi_high_level():
		try:
			paths = []
			for name in ["HMD_HumanPhenotype.rpt", "VOC_MammalianPhenotype.rpt"]:
				p = cache / name
				if not p.exists():
					z = requests.get(
						"https://www.informatics.jax.org/downloads/reports/" + name,
						timeout=45,
					)
					z.raise_for_status()
					atomic_text(p, z.text)
				paths.append(p)
			names = {
				c[0]: c[1]
				for c in csv.reader(paths[1].open(), delimiter="\t")
				if len(c) > 1
			}
			groups = {}
			for c in csv.reader(paths[0].open(), delimiter="\t"):
				if len(c) < 5 or c[0].upper() not in background:
					continue
				for term in c[4].replace(",", " ").split():
					groups.setdefault(term, set()).add(c[0].upper())
			universe = set().union(*groups.values()) if groups else set()
			query = sig & universe
			rows = []
			if not universe or not query:
				return []
			for term, g in groups.items():
				hits = g & query
				rows.append(
					dict(
						source="MGI",
						term_name=names.get(term, term),
						p_raw=(
							tail(len(hits), len(g), len(universe), len(query))
							if hits
							else 1.0
						),
						intersection_size=len(hits),
						genes=";".join(sorted(hits)),
					)
				)
			audit.append(
				dict(
					source="MGI",
					status="ok",
					detail=f"MGI high-level phenotypes; annotated assayed background N={len(universe)}, selected N={len(query)}",
				)
			)
			return bh(rows)
		except Exception as e:
			audit.append(dict(source="MGI", status="unavailable", detail=str(e)))
			return []

	if sig:
		terms += mgi_high_level()
		for label, names in [
			("KEGG", ["KEGG_2021_Human"]),
			("TF", ["TRRUST_Transcription_Factors_2019"]),
		]:
			rows = library(label, names)
			terms += rows
			if label == "TF":
				for r in [x for x in rows if x["adjusted_p"] < 0.05][:15]:
					regulator = r["term_name"].split("_")[0].split(" ")[0]
					for g in r["genes"].split(";"):
						if g and regulator != g:
							tf.append({"from": regulator, "to": g})
		key = hashlib.sha256("\n".join(sorted(sig)).encode()).hexdigest()[:20]
		p = cache / ("string_physical_" + key + ".tsv")
		try:
			if not p.exists():
				z = requests.post(
					"https://version-12-0.string-db.org/api/tsv/network",
					data={
						"identifiers": "\r".join(sorted(sig)),
						"species": 9606,
						"network_type": "physical",
						"required_score": 400,
					},
					timeout=60,
				)
				z.raise_for_status()
				atomic_text(p, z.text)
			for r in csv.DictReader(io.StringIO(p.read_text()), delimiter="\t"):
				a, b = r.get("preferredName_A", ""), r.get("preferredName_B", "")
				if a.upper() in sig and b.upper() in sig:
					ppi.append({"from": a, "to": b, "score": r["score"]})
			audit.append(
				dict(source="STRING", status="ok", detail=f"{len(ppi)} physical edges")
			)
		except Exception as e:
			audit.append(dict(source="STRING", status="unavailable", detail=str(e)))
	write(
		"c1.mock_function_terms.csv",
		terms,
		["source", "term_name", "p_raw", "intersection_size", "genes", "adjusted_p"],
	)
	write("c1.mock_tf_edges.csv", tf, ["from", "to"])
	write("c1.mock_ppi_edges.csv", ppi, ["from", "to", "score"])
	write("c1.mock_annotation_status.csv", audit, ["source", "status", "detail"])
	print(
		"MOCK enrichment:",
		len(terms),
		"terms;",
		len(tf),
		"TF edges;",
		len(ppi),
		"physical edges",
		flush=True,
	)


def annotations_main(argv):
	if argv[:1] == ["--gprofiler"]:
		run_gprofiler(argv[1:])
	elif argv[:1] in (["--help"], ["-h"]):
		print(
			"Usage: c1.abm.py annotations GENES_CSV OUTPUT_DIR CACHE_DIR\n"
			"       c1.abm.py annotations --gprofiler SIGNIFICANT_IDS BACKGROUND_IDS OUTPUT_CSV"
		)
	else:
		run_annotations(argv)


def render_abm_figures(out, primary="abm_transformer", missing_only=True):
	import pandas as pd

	out = Path(out)
	frozen = out / "MODEL_FROZEN.json"
	if frozen.is_file():
		primary = json.loads(frozen.read_text()).get("primary", primary)
	import matplotlib

	matplotlib.use("Agg")
	import matplotlib.pyplot as plt

	written = []

	def needed(filename):
		return not missing_only or not (out / filename).is_file()

	def save(fig, filename):
		fig.tight_layout()
		fig.savefig(out / filename, dpi=180)
		plt.close(fig)
		written.append(out / filename)

	name = "Fig_model_comparison.png"
	if needed(name) and (out / "test_metrics.csv").is_file():
		result = pd.read_csv(out / "test_metrics.csv")
		if {"subset", "model", "AUC_IPCW"} <= set(result.columns):
			main = result[result.subset.eq("all")].sort_values("AUC_IPCW")
			if not main.empty:
				fig, ax = plt.subplots(figsize=(10, max(6, 0.27 * len(main))))
				colors = ["#c04a31" if s == primary else "#4a718e" for s in main.model]
				ax.barh(main.model, main.AUC_IPCW, color=colors)
				ax.axvline(0.5, color="grey", ls="--")
				ax.set_xlim(0, 1)
				ax.set_xlabel("Full test IPCW AUC")
				save(fig, name)

	name = "Fig_masked_reconstruction.png"
	masked = out / "masked_reconstruction_summary.csv"
	if not masked.is_file():
		masked = out / "masked_reconstruction.csv"
	if needed(name) and masked.is_file():
		result = pd.read_csv(masked)
		if {"method", "masked_RMSE"} <= set(result.columns) and not result.empty:
			data = result.groupby("method").masked_RMSE.mean().sort_values()
			fig, ax = plt.subplots(figsize=(8, 4))
			ax.barh(data.index, data.values, color="#4a718e")
			ax.set_xlabel("Mean person-level masked RMSE (lower is better)")
			save(fig, name)

	name = "Fig_coverage.png"
	if needed(name) and (out / "coverage_curve.csv").is_file():
		data = pd.read_csv(out / "coverage_curve.csv")
		if {"model", "coverage", "Brier_IPCW"} <= set(data.columns) and not data.empty:
			fig, ax = plt.subplots(figsize=(7, 4))
			for model in [primary, "elasticnet", "abm_clinical"]:
				sub = data[data.model.eq(model)]
				if not sub.empty:
					ax.plot(sub.coverage, sub.Brier_IPCW, "o-", label=model)
			ax.set_xlabel("Supported fraction; thresholds frozen on tune")
			ax.set_ylabel("IPCW Brier on identical people")
			if ax.lines:
				ax.legend(fontsize=8)
			save(fig, name)
	return written


def figures_main(argv):
	p = argparse.ArgumentParser(
		prog="c1.abm.py figures",
		description="Redraw ABM figures from saved aggregate CSVs",
	)
	p.add_argument("run_dir", type=Path)
	p.add_argument("--primary", default="abm_transformer")
	p.add_argument("--replace", action="store_true", help="Redraw existing figures too")
	a = p.parse_args(argv)
	for path in render_abm_figures(a.run_dir, a.primary, missing_only=not a.replace):
		print(path)


# Annotation and figure recovery do not import Torch/sklearn or launch training.
if __name__ == "__main__" and sys.argv[1:2] in (["annotations"], ["figures"]):
	action = annotations_main if sys.argv[1] == "annotations" else figures_main
	action(sys.argv[2:])
	raise SystemExit(0)


def run(kind=None, argv=None):
	p = argparse.ArgumentParser(add_help=False)
	p.add_argument(
		"--backend",
		choices=["reference", "tabicl", "tf", "both"],
		default=kind or "reference",
		help="ABM method; reference also includes Transformer models; tf is a legacy alias for tabicl",
	)
	p.add_argument(
		"--memory-limit-gb",
		type=int,
		default=int(os.getenv("ABM_MEMORY_CAP_GB", "24")),
	)
	p.add_argument(
		"--memory-swap-gb", type=int, default=int(os.getenv("ABM_SWAP_CAP_GB", "2"))
	)
	a, rest = p.parse_known_args(sys.argv[1:] if argv is None else argv)
	if a.memory_limit_gb < 0 or a.memory_swap_gb < 0:
		raise ValueError("Invalid process-tree memory cap")
	selected = "tabicl" if a.backend == "tf" else a.backend
	if selected == "both":
		if any(x == "--run-dir" or x.startswith("--run-dir=") for x in rest):
			raise ValueError(
				"For --backend both use --analysis-root; --run-dir identifies one backend run."
			)
		for backend in ["reference", "tabicl"]:
			code = run(
				argv=[
					"--backend",
					backend,
					"--memory-limit-gb",
					str(a.memory_limit_gb),
					"--memory-swap-gb",
					str(a.memory_swap_gb),
					*rest,
				]
			)
			if code:
				return code
		return 0
	kind = selected
	defaults = ["--device", "cuda", "--cores", "16"]
	if kind == "reference":
		defaults += ["--quality-teacher", "all"]
	if any(x in rest for x in ["--help", "-h"]):
		print(
			"Unified ABM: --backend reference|tabicl|both (reference includes Transformer models)",
			flush=True,
		)
	command = [
		sys.executable,
		str(HERE / "c1.abm.py"),
		"--_abm-worker",
		kind,
		*defaults,
		*rest,
	]
	# Help/dry run is entirely local and must not require GPU, data or systemd.
	if any(x in rest for x in ["--help", "-h", "--dry-run"]):
		if "--dry-run" in rest:
			print("[LE8 ABM] " + shlex.join(command))
			return 0
		return subprocess.call(command)
	env = dict(os.environ)
	# Actual fits/project/evaluate keep the process-tree guard (not per-child ulimit).
	if a.memory_limit_gb and env.get("LE8_MEMORY_SCOPE_ACTIVE") != "1":
		if not shutil.which("systemd-run"):
			raise RuntimeError(
				"systemd-run unavailable. Set up user systemd, or explicitly --memory-limit-gb 0 to opt out."
			)
		probe = subprocess.run(
			["systemctl", "--user", "is-system-running"], capture_output=True, text=True
		)
		if probe.stdout.strip() not in {"running", "degraded"}:
			raise RuntimeError(
				"User systemd is not ready for the requested process-tree guard. Explicit --memory-limit-gb 0 disables it."
			)
		command = [
			"systemd-run",
			"--user",
			"--scope",
			"--quiet",
			"--collect",
			"-p",
			f"MemoryMax={a.memory_limit_gb}G",
			"-p",
			f"MemorySwapMax={a.memory_swap_gb}G",
			"-p",
			"OOMPolicy=kill",
			"--",
			"env",
			"LE8_MEMORY_SCOPE_ACTIVE=1",
			*command,
		]
	print("[LE8 ABM] " + shlex.join(command), flush=True)
	return subprocess.call(command, env=env)


# Re-enter this file in a supervised worker. Help/dry-run needs no scientific stack.
if __name__ == "__main__" and "--_abm-worker" not in sys.argv[1:]:
	raise SystemExit(run())


# ABM
# ABM 5: Transformer-assisted individual reference copying and disease risk.


def reference_parser():
	p = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog="""Usage examples:
  python f/c1.abm.py --Y cvd_cad --biom prot --dry-run
  python f/c1.abm.py --Y cvd_cad --biom prot --preflight
  python f/c1.abm.py --Y cvd_cad --biom prot --cores 16
  python f/c1.abm.py evaluate --run-dir /mnt/d/analysis/le8/cvd_cad/prot/c1_correlate/abm_reference
  python f/c1.abm.py --demo --tree hist --max-samples 1600 --analysis-root /tmp/abm_check
""",
	)
	p.add_argument(
		"command", nargs="?", choices=["run", "evaluate", "project"], default="run"
	)
	p.add_argument("--Y", dest="trait", default="cvd_cad")
	p.add_argument("--biom", choices=["prot", "met"], default="prot")
	p.add_argument(
		"--ukb-phe", default=os.environ.get("UKB_PHE", "/mnt/d/data/ukb/phe")
	)
	p.add_argument("--phe-file")
	p.add_argument("--omics-file")
	p.add_argument("--input-source", choices=["raw", "cleaned"], default="raw")
	p.add_argument("--met-input", choices=["raw", "named"], default=None)
	p.add_argument("--met-map")
	p.add_argument("--id-col", default="eid")
	p.add_argument("--baseline-col", default="date_attend")
	p.add_argument("--diagnosis-col")
	p.add_argument("--death-col", default="date_death")
	p.add_argument("--lost-col", default="date_lost")
	p.add_argument("--end-date", default="2023-04-01")
	p.add_argument("--disease-evidence-col", default="")
	p.add_argument("--healthy-date-cols", default="")
	p.add_argument("--covariates", default="age,sex,tdi,PC1,PC2,center")
	p.add_argument("--residualize", default=None)
	p.add_argument("--categorical", default=None)
	p.add_argument("--group-col", default="")
	p.add_argument("--split-file", default="")
	p.add_argument("--module-file", default="")
	p.add_argument("--transform", choices=["none", "log1p"], default=None)
	p.add_argument("--exclude-features", default="")
	p.add_argument("--feature-missing", type=float, default=0.2)
	p.add_argument("--sample-missing", type=float, default=0.2)
	p.add_argument("--horizon", type=float, default=10)
	p.add_argument("--min-censor-survival", type=float, default=0.05)
	p.add_argument("--min-events", type=int, default=20)
	p.add_argument("--panel-size", type=int, default=100)
	p.add_argument("--panel-sizes", default="100,300,1000")
	p.add_argument("--match-k", default="5,10,20")
	p.add_argument("--all-k", default="30,100,300")
	p.add_argument("--folds", type=int, default=5)
	p.add_argument("--oof-repeats", type=int, default=3)
	p.add_argument("--teacher-c", type=float, default=0.01)
	p.add_argument(
		"--quality-teacher",
		choices=["elasticnet", "ensemble", "neural", "all"],
		default="ensemble",
	)
	p.add_argument("--quality-neural-epochs", type=int, default=20)
	p.add_argument("--quality-pretrain-epochs", type=int, default=5)
	p.add_argument("--quality-trees", type=int, default=100)
	p.add_argument("--c-grid", default="0.001,0.01,0.1,1")
	p.add_argument("--max-iter", type=int, default=3000)
	p.add_argument("--fit-stability", type=float, default=2 / 3)
	p.add_argument("--min-gain", type=float, default=0.0)
	p.add_argument("--dimensions", type=int, default=32)
	p.add_argument("--accept-quantile", type=float, default=0.95)
	p.add_argument("--min-match-ess", type=float, default=3.0)
	p.add_argument("--min-coverage", type=float, default=0.5)
	p.add_argument("--tree", choices=["lightgbm", "hist", "none"], default="lightgbm")
	p.add_argument("--tree-estimators", type=int, default=300)
	p.add_argument("--bootstrap", type=int, default=200)
	p.add_argument("--pair-caliper", type=float, default=0.01)
	p.add_argument("--max-pairs", type=int, default=100)
	p.add_argument("--explanation-samples", type=int, default=1000)
	p.add_argument("--mask-fraction", type=float, default=0.1)
	p.add_argument("--mask-repeats", type=int, default=3)
	p.add_argument(
		"--analysis-root", default=os.getenv("LE8_ANALYSIS_ROOT", "/mnt/d/analysis/le8")
	)
	p.add_argument("--run-dir")
	p.add_argument("--output", help="CSV output for project")
	p.add_argument("--seed", type=int, default=2026)
	p.add_argument("--cores", type=int, default=4)
	p.add_argument("--r-bin", default="Rscript")
	p.add_argument("--max-samples", type=int, default=0)
	p.add_argument("--demo-features", type=int, default=80)
	p.add_argument("--demo", action="store_true")
	p.add_argument("--train-only", action="store_true")
	p.add_argument("--dry-run", action="store_true")
	p.add_argument("--preflight", action="store_true")
	p.add_argument("--full-input-hash", action="store_true")
	p.add_argument("--replace", action="store_true")
	p.add_argument(
		"--experiments", default="transformer,mlp,direct,no_pretrain,uniform,metric"
	)
	p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
	p.add_argument("--check-device", action="store_true")
	p.add_argument("--tokens", type=int, default=48)
	p.add_argument("--width", type=int, default=64)
	p.add_argument("--heads", type=int, default=4)
	p.add_argument("--layers", type=int, default=2)
	p.add_argument("--dropout", type=float, default=0.1)
	p.add_argument("--batch-size", type=int, default=128)
	p.add_argument("--pretrain-epochs", type=int, default=20)
	p.add_argument("--epochs", type=int, default=60)
	p.add_argument("--patience", type=int, default=10)
	p.add_argument("--neural-k", type=int, default=64)
	p.add_argument("--learning-rate", type=float, default=0.0003)
	p.add_argument("--weight-decay", type=float, default=0.001)
	p.add_argument("--reconstruction-weight", type=float, default=0.2)
	p.add_argument("--direct-weight", type=float, default=0.3)
	p.add_argument("--train-prior", type=float, default=2.0)
	p.add_argument("--gradient-clip", type=float, default=5.0)
	p.add_argument("--temperatures", default="0.05,0.2,1,5")
	p.add_argument("--prior-grid", default="0,2,10")
	p.add_argument(
		"--attention-samples",
		type=int,
		default=128,
		help="Outcome-blind ID sample for head maps and mechanistic perturbations; 0 disables",
	)
	p.add_argument("--perturbation-samples", type=int, default=256)
	p.add_argument("--amp", action=argparse.BooleanOptionalAction, default=True)
	p.add_argument("--deterministic", action="store_true")
	p.add_argument(
		"--resume",
		action="store_true",
		help="Resume neural epochs after interruption; configuration and inputs must match",
	)
	p.add_argument(
		"--shuffle-development-outcomes",
		action="store_true",
		help="Negative control: permute time/event pairs independently within each non-test role",
	)
	return p


def configure(a):
	a.outcome_type = getattr(a, "outcome_type", "survival")
	for field, cast in [
		("panel_sizes", int),
		("match_k", int),
		("all_k", int),
		("c_grid", float),
	]:
		setattr(a, field, sorted(set(cast(v) for v in getattr(a, field).split(","))))
		if not getattr(a, field) or min(getattr(a, field)) <= 0:
			raise ValueError(f"Invalid {field}")
	a.panel_sizes = sorted(set(a.panel_sizes + [a.panel_size]))
	a.outcome_type, a.target_col = "survival", a.trait
	a.diagnosis_col = a.diagnosis_col or "fod_icd10_" + a.trait
	base = Path(a.ukb_phe)
	a.phe_file = a.phe_file or str(base / "Rdata/all.rds")
	a.omics_file = a.omics_file or str(
		base
		/ (
			f"Rdata/{a.biom}.rds"
			if a.input_source == "cleaned"
			else "rap/raw/prot.tab.gz"
			if a.biom == "prot"
			else "rap/met.tab.gz"
		)
	)
	a.met_map = a.met_map or str(base / "common/met.lst")
	a.met_input = a.met_input or (
		"raw" if a.biom == "met" and a.input_source == "raw" else "named"
	)
	a.residualize = f"{a.biom}.plate" if a.residualize is None else a.residualize
	a.categorical = (
		f"sex,center,{a.biom}.plate" if a.categorical is None else a.categorical
	)
	a.transform = a.transform or (
		"log1p"
		if a.biom == "met" and a.input_source == "raw" and not a.demo
		else "none"
	)
	if (
		a.horizon <= 0
		or a.panel_size < 3
		or a.dimensions < 1
		or a.folds < 2
		or a.oof_repeats < 2
	):
		raise ValueError(
			"Invalid horizon, panel size, dimensions, folds, or OOF repeats (minimum 2)"
		)
	if not 0 < a.accept_quantile <= 1 or not 0 < a.fit_stability <= 1:
		raise ValueError("Invalid acceptance/stability quantile")
	if not 0 < a.mask_fraction < 1 or a.mask_repeats < 1 or a.explanation_samples < 0:
		raise ValueError("Invalid masking sensitivity settings")
	for field in ["feature_missing", "sample_missing", "min_censor_survival"]:
		if not 0 < getattr(a, field) < 1:
			raise ValueError(f"Invalid {field}")
	a.experiments = [v.strip() for v in a.experiments.split(",") if v.strip()]
	if (
		"transformer" not in a.experiments
		or len(set(a.experiments)) != len(a.experiments)
		or not set(a.experiments)
		<= {"transformer", "mlp", "direct", "no_pretrain", "uniform", "metric"}
	):
		raise ValueError(
			"--experiments requires transformer; other choices: mlp,direct,no_pretrain,uniform,metric; no duplicates"
		)
	a.temperatures = sorted(set(float(v) for v in a.temperatures.split(",")))
	a.prior_grid = sorted(set(float(v) for v in a.prior_grid.split(",")))
	if (
		not a.temperatures
		or min(a.temperatures) <= 0
		or not a.prior_grid
		or min(a.prior_grid) < 0
	):
		raise ValueError("Invalid temperature/prior grid")
	if (
		a.width < 8
		or a.heads < 1
		or a.width % a.heads
		or a.tokens < 2
		or a.layers < 1
		or a.batch_size < 2
	):
		raise ValueError("Invalid Transformer dimensions; width must divide by heads")
	if (
		a.epochs < 1
		or a.pretrain_epochs < 0
		or a.patience < 1
		or a.neural_k < 2
		or not 0 <= a.dropout < 1
	):
		raise ValueError("Invalid neural training settings")
	if a.quality_neural_epochs < 1 or a.quality_pretrain_epochs < 0:
		raise ValueError("Invalid neural OOF budget")
	if (
		a.learning_rate <= 0
		or a.weight_decay < 0
		or a.reconstruction_weight < 0
		or a.direct_weight < 0
		or a.train_prior < 0
		or a.gradient_clip <= 0
	):
		raise ValueError("Invalid loss/optimizer settings")
	if (
		a.attention_samples < 0
		or a.perturbation_samples < 0
		or a.bootstrap < 0
		or a.min_events < 1
		or not 0 <= a.min_coverage <= 1
		or a.min_match_ess < 1
	):
		raise ValueError("Invalid evaluation settings")
	if a.resume and a.replace:
		raise ValueError("Choose --resume or --replace")
	return a


def output_directory(a, backend="reference"):
	return (
		Path(a.run_dir)
		if a.run_dir
		else Path(a.analysis_root)
		/ a.trait
		/ a.biom
		/ "c1_correlate"
		/ ("abm_" + backend)
	)


def abm_cache_dir(out, *parts):
	out = Path(out).resolve()
	if out.parent.name == "c1_correlate" and out.name in {"abm_reference", "abm_tabicl"}:
		project = out.parents[3]
		if project == Path(os.getenv("LE8_ANALYSIS_ROOT", str(project))).resolve():
			project = Path(os.getenv("LE8_PUBLISHED_ROOT", str(project))).resolve()
		relative = Path(out.parents[2].name) / out.parents[1].name / out.name
	else:
		project = out
		relative = Path("abm")
	key = hashlib.sha256(str(project).encode()).hexdigest()[:12]
	cache = Path("/tmp/le8-cache") / key / relative
	cache = cache.joinpath(*parts)
	cache.mkdir(parents=True, exist_ok=True)
	return cache


def reference_main():
	a = configure(reference_parser().parse_args())
	if a.cores < 1:
		raise ValueError("--cores must be positive")
	for name in [
		"OMP_NUM_THREADS",
		"OPENBLAS_NUM_THREADS",
		"MKL_NUM_THREADS",
		"NUMEXPR_NUM_THREADS",
	]:
		os.environ[name] = str(a.cores)
	# Imports occur after thread limits have been set.

	if a.check_device:
		import torch, json

		print(
			json.dumps(
				{
					"torch": torch.__version__,
					"cuda_available": torch.cuda.is_available(),
					"cuda_version": torch.version.cuda,
					"gpu": (
						torch.cuda.get_device_name(0)
						if torch.cuda.is_available()
						else None
					),
				},
				indent=2,
			)
		)
		return

	from threadpoolctl import threadpool_limits

	out = output_directory(a)
	if a.command == "project":
		if not a.run_dir or not a.output:
			raise ValueError("project requires --run-dir and --output")
		with threadpool_limits(a.cores):
			reference_project(
				out, a.phe_file, a.omics_file, a.output, a.r_bin, a.met_input, a.device
			)
		return
	if a.dry_run:
		import json

		print(json.dumps(dict(output=str(out), config=vars(a)), indent=2))
		return
	if a.preflight:
		import importlib.util

		paths = [] if a.demo else [a.phe_file, a.omics_file]
		if a.module_file:
			paths.append(a.module_file)
		if a.split_file:
			paths.append(a.split_file)
		if a.biom == "met" and a.met_input == "raw" and not a.demo:
			paths.append(a.met_map)
		missing = [v for v in paths if not Path(v).is_file()]
		if missing:
			raise FileNotFoundError("Missing inputs: " + ", ".join(missing))
		if any(
			str(v).lower().endswith(".rds") for v in paths
		) and not importlib.util.find_spec("pyreadr"):
			raise RuntimeError(
				"RDS needs pyreadr for lossless binary reading; install requirements.txt"
			)
		if a.tree == "lightgbm" and not importlib.util.find_spec("lightgbm"):
			raise RuntimeError(
				"LightGBM missing; install requirements.txt or explicitly choose --tree hist"
			)

		device_for(a.device)
		log(
			"DONE",
			"preflight",
			"paths and required reader/tree dependencies available; data values not yet read",
		)
		return
	if a.command == "evaluate":
		if not (out / "MODEL_FROZEN.json").is_file():
			raise ValueError("No frozen model found")
		with run_lock(out), threadpool_limits(a.cores):
			reference_evaluate(out)
		return
	if a.tree == "lightgbm":
		import importlib.util

		if not importlib.util.find_spec("lightgbm"):
			raise RuntimeError(
				"LightGBM missing; install requirements.txt or explicitly choose --tree hist before training"
			)
	out.mkdir(parents=True, exist_ok=True)
	with run_lock(out), threadpool_limits(a.cores):
		if not prepare_run_directory(
			out, resume=a.resume, replace=a.replace, train_only=a.train_only
		):
			return
		reference_train(a, out)
		dump(out / "TRAIN_DONE.json", dict(version=VERSION))
		if not a.train_only:
			log("START", "test_evaluation")
			reference_evaluate(out)


# Tf Download
# Pinned official TabICLv2 weights; download is explicit, never part of report/Shiny.
REPO = "jingang/TabICL"
MODEL_COMMIT = "4dcd344ece2c00be9e831fdd35bed57b5ad83e19"
FILENAME = "tabicl-classifier-v2-20260212.ckpt"
SHA256 = "bdc7dbd5e4ff21f8f0456fcf90c6b7cdf72dbea960f2d05b19bec19f9b3d4ed0"
DEFAULT_MODEL = "/mnt/f/AI_model/TabICL/" + FILENAME


def verify_model(path):
	path = Path(path)
	if not path.is_file() or path.stat().st_size < 1_000_000:
		raise FileNotFoundError(
			f"Missing checkpoint: {path}; use --download-model explicitly"
		)
	h = hashlib.sha256()
	with path.open("rb") as f:
		for block in iter(lambda: f.read(1024 * 1024), b""):
			h.update(block)
	if h.hexdigest() != SHA256:
		raise ValueError(
			"Checkpoint checksum does not match the pinned official TabICLv2 weights"
		)
	return dict(
		repo=REPO,
		revision=MODEL_COMMIT,
		file=str(path),
		sha256=h.hexdigest(),
		bytes=path.stat().st_size,
	)


def download_model(path=DEFAULT_MODEL):
	path = Path(path)
	if path.name != FILENAME:
		raise ValueError(f"Checkpoint filename must be {FILENAME}")
	try:
		return verify_model(path)
	except (ValueError, FileNotFoundError):
		pass
	from huggingface_hub import snapshot_download

	path.parent.mkdir(parents=True, exist_ok=True)
	snapshot_download(
		repo_id=REPO,
		revision=MODEL_COMMIT,
		local_dir=str(path.parent),
		allow_patterns=[FILENAME, "README.md"],
		token=os.environ.get("HF_TOKEN", False),
	)
	return verify_model(path)


# Tabicl Backend
# Embedded numerical Transformer entry; external weights/environments remain explicit.


def tabicl_parser():
	p = reference_parser()
	p.description = __doc__
	p.set_defaults(
		epochs=10,
		patience=3,
		learning_rate=1e-5,
		weight_decay=0.01,
		batch_size=64,
	)
	p.add_argument("--model-path", default=DEFAULT_MODEL)
	p.add_argument("--download-model", action="store_true")
	p.add_argument("--max-features", type=int, default=256)
	p.add_argument("--context-size", type=int, default=256)
	p.add_argument("--predict-batch-size", type=int, default=256)
	p.add_argument("--validation-samples", type=int, default=2048)
	p.add_argument("--with-clinical", action="store_true")
	p.add_argument("--unfreeze-encoder", action="store_true")
	return p


def tabicl_main():
	a = configure(tabicl_parser().parse_args())
	for k in [
		"OMP_NUM_THREADS",
		"OPENBLAS_NUM_THREADS",
		"MKL_NUM_THREADS",
		"NUMEXPR_NUM_THREADS",
	]:
		os.environ[k] = str(a.cores)
	if (
		a.max_features < 0
		or a.context_size < 4
		or a.predict_batch_size < 1
		or a.validation_samples < 0
	):
		raise ValueError("Invalid TF settings")
	if getattr(a, "resume", False):
		raise ValueError(
			"TF epoch resume is not implemented. Keep the old run and use a new --run-dir, or explicitly --replace it."
		)
	if a.dry_run:
		print(json.dumps(vars(a), indent=2))
		return
	if a.download_model:
		print(json.dumps(download_model(a.model_path), indent=2))
		return
	import torch

	from threadpoolctl import threadpool_limits

	if a.check_device:
		print(
			json.dumps(
				dict(torch=torch.__version__, cuda_available=torch.cuda.is_available())
			)
		)
		return
	if importlib.metadata.version("tabicl") != "2.2.0":
		raise RuntimeError("This adapter requires tabicl==2.2.0")
	a.device = device_for(a.device)
	out = output_directory(a, "tabicl")
	if a.command == "project":
		if not a.run_dir or not a.output:
			raise ValueError("project requires --run-dir and --output")
		with threadpool_limits(a.cores):
			tabicl_project(a, out)
		return
	if a.command == "evaluate":
		with run_lock(out), threadpool_limits(a.cores):
			tabicl_evaluate(out)
		return
	info = verify_model(a.model_path)
	if not a.demo and any(not Path(x).is_file() for x in [a.phe_file, a.omics_file]):
		raise FileNotFoundError("TF phenotype/omics input missing")
	if a.preflight:
		print("TF dependencies/weights/input paths checked; no model fitted")
		return
	out.mkdir(parents=True, exist_ok=True)
	with run_lock(out), threadpool_limits(a.cores):
		if not prepare_run_directory(out, replace=a.replace, train_only=a.train_only):
			return
		tabicl_train(a, out, info)
		if not a.train_only:
			tabicl_evaluate(out)


# Parse help and set thread limits before importing NumPy, Torch or sklearn.
_WORKER = None
if __name__ == "__main__":
	if len(sys.argv) < 3 or sys.argv[1] != "--_abm-worker":
		raise ValueError("Invalid internal worker invocation")
	_WORKER = sys.argv[2]
	if _WORKER not in {"reference", "tabicl"}:
		raise ValueError("Unknown ABM worker")
	del sys.argv[1:3]
	_args = (
		reference_parser() if _WORKER == "reference" else tabicl_parser()
	).parse_args()
	if _args.cores < 1:
		raise ValueError("--cores must be positive")
	for _key in (
		"OMP_NUM_THREADS",
		"OPENBLAS_NUM_THREADS",
		"MKL_NUM_THREADS",
		"NUMEXPR_NUM_THREADS",
	):
		os.environ[_key] = str(_args.cores)
	# Pickles created through the CLI use an importable, stable module name.
	sys.modules["c1_abm"] = sys.modules[__name__]


from scipy import sparse
from scipy.optimize import minimize
from scipy.special import logit, expit
from scipy.special import softmax
from scipy.special import softmax, logsumexp
from scipy.stats import spearmanr
from sklearn.cluster import MiniBatchKMeans
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import roc_curve
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.model_selection import train_test_split, GroupShuffleSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler, SplineTransformer
from torch import nn
from torch.nn import functional as F
import joblib
import numpy as np
import pandas as pd
import torch

# Common
# Small shared helpers. All learned objects are local, trusted research artifacts.

VERSION = "5.1.0"


def words(value):
	return [s.strip() for s in (value or "").split(",") if s.strip()]


def json_safe(value):
	if isinstance(value, dict):
		return {str(k): json_safe(v) for k, v in value.items()}
	if isinstance(value, (list, tuple, np.ndarray)):
		return [json_safe(v) for v in value]
	if isinstance(value, (np.integer, np.bool_)):
		return value.item()
	if isinstance(value, (float, np.floating)):
		return float(value) if np.isfinite(value) else None
	if isinstance(value, Path):
		return str(value)
	return value


def dump(path, value):
	path = Path(path)
	temp = path.with_suffix(path.suffix + ".tmp")
	temp.write_text(
		json.dumps(json_safe(value), ensure_ascii=False, indent=2, allow_nan=False)
		+ "\n",
		encoding="utf-8",
	)
	temp.replace(path)


def write_array_rds(path, arrays, r_bin="Rscript"):
	"""Preserve participant attention arrays in one named RDS, without text rounding."""
	path = Path(path)
	with tempfile.TemporaryDirectory(prefix="le8-array-rds-", dir="/tmp") as scratch:
		work = Path(scratch)
		specs = {}
		for i, (name, value) in enumerate(arrays.items()):
			value = np.asarray(value)
			entry = dict(shape=list(value.shape), dtype=str(value.dtype))
			if value.dtype.kind in "US":
				entry["values"] = value.astype(str).ravel(order="F").tolist()
			elif value.dtype.kind == "f":
				entry["file"] = str(work / f"{i}.bin")
				value.astype("<f8").ravel(order="F").tofile(entry["file"])
			else:
				raise ValueError(f"Unsupported participant array type: {value.dtype}")
			specs[name] = entry
		(work / "arrays.json").write_text(json.dumps(specs), encoding="utf-8")
		script = """args <- commandArgs(TRUE)
specs <- jsonlite::fromJSON(args[1], simplifyVector = FALSE)
value <- lapply(specs, function(entry) {
	dims <- as.integer(unlist(entry$shape))
	x <- if (is.null(entry$file)) unlist(entry$values, use.names = FALSE) else
		readBin(entry$file, what = 'double', n = prod(dims), size = 8L, endian = 'little')
	stopifnot(length(x) == prod(dims))
	dim(x) <- dims
	attr(x, 'numpy_dtype') <- entry$dtype
	x
})
saveRDS(value, args[2], compress = 'gzip')
stopifnot(identical(value, readRDS(args[2])))
"""
		(work / "save.R").write_text(script, encoding="utf-8")
		stored = work / "arrays.rds"
		subprocess.run([r_bin, str(work / "save.R"), str(work / "arrays.json"), str(stored)], check=True)
		path.parent.mkdir(parents=True, exist_ok=True)
		shutil.copyfile(stored, path)
		if hashlib.sha256(stored.read_bytes()).digest() != hashlib.sha256(path.read_bytes()).digest():
			raise OSError(f"RDS copy verification failed: {path}")
		path.chmod(0o600)


def digest(value):
	return hashlib.sha256(
		json.dumps(json_safe(value), sort_keys=True).encode()
	).hexdigest()


def log(action, step, detail=""):
	print(f"[ABM] {action} {step}" + (f" | {detail}" if detail else ""), flush=True)


@contextmanager
def stage_log(step):
	start = time.monotonic()
	log("START", step)
	yield
	log("DONE", step, f"{time.monotonic() - start:.1f}s")


def bh(p):
	p = np.asarray(p, float)
	q = np.full(len(p), np.nan)
	ok = np.flatnonzero(np.isfinite(p))
	ix = ok[np.argsort(p[ok])]
	q[ix] = np.minimum(
		1,
		np.minimum.accumulate((p[ix] * len(ix) / np.arange(1, len(ix) + 1))[::-1])[
			::-1
		],
	)
	return q


def fingerprints(paths, full=False):
	rows = []
	for path in sorted(set(str(Path(p).resolve()) for p in paths if p)):
		p = Path(path)
		if not p.is_file():
			raise FileNotFoundError(p)
		stat = p.stat()
		row = dict(path=path, bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
		if full or stat.st_size < 2_000_000:
			h = hashlib.sha256()
			with p.open("rb") as handle:
				for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
					h.update(block)
			row["sha256"] = h.hexdigest()
		rows.append(row)
	return rows


@contextmanager
def run_lock(root):
	path = root / ".lock"
	try:
		fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
	except FileExistsError as exc:
		raise RuntimeError(
			f"Run is locked: {path}. Inspect its PID before removing a stale lock."
		) from exc
	with os.fdopen(fd, "w") as f:
		f.write(str(os.getpid()))
	try:
		yield
	finally:
		path.unlink(missing_ok=True)


def run_complete(root, train_only=False):
	"""A frozen model alone is not a completed run: audits may still fail."""
	root = Path(root)
	required = ["manifest.json", "MODEL_FROZEN.json", "model_bundle.joblib"]
	marker = "TRAIN_DONE.json" if train_only else "DONE.json"
	if train_only and (root / "DONE.json").is_file():
		return run_complete(root)
	if not train_only:
		required += ["REPORT.md", "test_metrics.csv"]
	if any(
		not (root / name).is_file() or (root / name).stat().st_size == 0
		for name in required + [marker]
	):
		return False
	try:
		record = json.loads((root / marker).read_text())
		return isinstance(record, dict) and bool(record.get("version"))
	except (ValueError, OSError):
		return False


def prepare_run_directory(root, resume=False, replace=False, train_only=False):
	"""Called under run_lock; keep that lock held while clearing old outputs."""
	root = Path(root)
	entries = [p for p in root.iterdir() if p.name != ".lock"]
	if not entries:
		return True
	manifest = root / "manifest.json"
	if not manifest.is_file():
		raise ValueError(
			"Refusing to overwrite a nonempty directory without a ABM manifest"
		)
	try:
		record = json.loads(manifest.read_text())
	except (ValueError, OSError) as exc:
		raise ValueError(
			"Refusing to overwrite a directory with an unreadable ABM manifest"
		) from exc
	if (
		not isinstance(record, dict)
		or not {"version", "config", "signature"} <= record.keys()
	):
		raise ValueError(
			"Refusing to overwrite a directory without a valid ABM manifest"
		)
	if resume:
		return True
	if not replace and run_complete(root, train_only):
		log("SKIP", "completed", str(root))
		return False
	log("RESTART", "replace" if replace else "incomplete_or_failed", str(root))
	for path in entries:
		if path.is_dir() and not path.is_symlink():
			shutil.rmtree(path)
		else:
			path.unlink()
	return True


# Survival
# Fixed-horizon IPCW learning. Death is censored: this estimates net risk.
#
# Censoring distribution is estimated in the reference-building sample only.
# Weights require independent censoring; no clipping conceals lack of support.
#


class CensoringKM:
	def fit(self, p, horizon, min_g=0.05):
		t = p.time.to_numpy(float)
		e = p.event.to_numpy(bool)
		self.horizon = float(horizon)
		self.times, inverse, counts = np.unique(
			t, return_inverse=True, return_counts=True
		)
		at_risk = len(t) - np.r_[0, np.cumsum(counts)[:-1]]
		# With tied event/censor times, observed failures leave before censoring.
		failures = np.bincount(inverse, weights=e, minlength=len(counts))
		censored = counts - failures
		denom = at_risk - failures
		step = np.divide(censored, denom, out=np.zeros_like(failures), where=denom > 0)
		self.survival = np.cumprod(1 - step)
		self.min_g = min_g
		if np.sum(t >= horizon) < 10 or self.at([horizon])[0] < min_g:
			raise ValueError(
				"Horizon lacks censoring support in build sample; shorten --horizon"
			)
		return self

	def at(self, t, left=False):
		ix = (
			np.searchsorted(self.times, np.asarray(t), side="left" if left else "right")
			- 1
		)
		return np.where(ix >= 0, self.survival[np.maximum(ix, 0)], 1.0)

	def labels_weights(self, p):
		t, e = p.time.to_numpy(float), p.event.to_numpy(bool)
		case = e & (t <= self.horizon)
		control = t > self.horizon
		g = np.ones(len(t))
		g[case] = self.at(t[case], left=True)
		g[control] = self.at(np.full(control.sum(), self.horizon))
		known = case | control
		if np.any(g[known] < self.min_g):
			raise ValueError(
				"IPCW unsupported: censor survival below --min-censor-survival"
			)
		w = np.zeros(len(t))
		w[known] = 1 / g[known]
		return case.astype(int), w


def loss(y, p):
	p = np.clip(p, 1e-6, 1 - 1e-6)
	return -(y * np.log(p) + (1 - y) * np.log1p(-p))


def metrics(y, w, p):
	y, w, p = np.asarray(y), np.asarray(w), np.asarray(p)
	if not len(y) or not np.isfinite(p).all() or w.sum() <= 0:
		return dict(
			AUC_IPCW=np.nan,
			AUPRC_IPCW=np.nan,
			Brier_IPCW=np.nan,
			Brier_Hajek=np.nan,
			LogLoss_IPCW=np.nan,
			observed_IPCW=np.nan,
		)
	auc = ap = np.nan
	if len(np.unique(y[w > 0])) == 2:
		auc = roc_auc_score(y, p, sample_weight=w)
		ap = average_precision_score(y, p, sample_weight=w)
	return dict(
		AUC_IPCW=float(auc),
		AUPRC_IPCW=float(ap),
		Brier_IPCW=float(np.mean(w * (y - p) ** 2)),
		Brier_Hajek=float(np.average((y - p) ** 2, weights=w)),
		LogLoss_IPCW=float(np.mean(w * loss(y, p))),
		observed_IPCW=float(np.average(y, weights=w)),
	)


# Preprocess
# Training-only QC, clinical coding and molecular residualization.


def parsed_date(raw, name):
	if pd.api.types.is_numeric_dtype(raw) and raw.notna().any():
		raise ValueError(
			f"{name}: numeric dates are ambiguous; use ISO dates or R Date"
		)
	converted = pd.to_datetime(raw, errors="coerce")
	bad = raw.notna() & raw.astype(str).str.strip().ne("") & converted.isna()
	if bad.any():
		raise ValueError(f"Invalid dates in {name}")
	return converted


def outcomes(p, a):
	p = p.copy()
	if a.outcome_type == "quantitative":
		raw = p[a.target_col]
		value = pd.to_numeric(raw, errors="coerce")
		if (raw.notna() & value.isna()).any():
			raise ValueError(f"Non-numeric quantitative outcome: {a.target_col}")
		p["target"] = value
		p["eligible"] = np.isfinite(value)
		return p, dict(
			joined=len(p),
			eligible=int(p.eligible.sum()),
			missing_target=int((~p.eligible).sum()),
			outcome_type="quantitative",
		)
	baseline = parsed_date(p[a.baseline_col], a.baseline_col)
	diagnosis = parsed_date(p[a.diagnosis_col], a.diagnosis_col)
	censor = pd.concat(
		[
			parsed_date(p[a.death_col], a.death_col),
			parsed_date(p[a.lost_col], a.lost_col),
			pd.Series(pd.Timestamp(a.end_date), index=p.index),
		],
		axis=1,
	).min(axis=1)
	prevalent = diagnosis.notna() & (diagnosis <= baseline)
	unknown = pd.Series(False, index=p.index)
	if a.disease_evidence_col:
		evidence = pd.to_numeric(p[a.disease_evidence_col], errors="raise")
		unknown = evidence.gt(0) & diagnosis.isna()
	valid = baseline.notna() & (censor > baseline)
	healthy = pd.Series(True, index=p.index)
	for col in words(a.healthy_date_cols):
		d = parsed_date(p[col], col)
		healthy &= ~(d.notna() & (d <= baseline))
	event = (
		valid
		& ~prevalent
		& diagnosis.notna()
		& (diagnosis > baseline)
		& (diagnosis <= censor)
	)
	p["time"] = (censor.where(~event, diagnosis) - baseline).dt.days / 365.25
	p["event"] = event.astype(int)
	p["eligible"] = valid & ~prevalent & ~unknown & healthy
	return p, dict(
		joined=len(p),
		eligible=int(p.eligible.sum()),
		prevalent=int(prevalent.sum()),
		same_day=int((diagnosis == baseline).sum()),
		invalid_followup=int((~valid).sum()),
		unknown_diagnosis=int(unknown.sum()),
		other_baseline_disease=int((~healthy & ~prevalent).sum()),
		diagnosis_after_censor=int((diagnosis > censor).sum()),
		incident=int(p.loc[p.eligible, "event"].sum()),
		outcome_type="survival",
	)


class MetadataDesign:
	"""Normalize categorical types and fit all encodings on training rows only."""

	def __init__(self, cols, categorical, age_spline=True, sparse_output=False):
		self.cols = list(cols)
		self.categorical = [c for c in self.cols if c in categorical]
		self.age_spline = age_spline
		self.sparse_output = sparse_output

	def normalized(self, p):
		p = p[self.cols].copy()
		for col in self.cols:
			if col in self.categorical:
				p[col] = (
					p[col]
					.map(
						lambda v: (
							np.nan
							if pd.isna(v)
							else (
								str(int(v))
								if isinstance(v, (float, np.floating))
								and np.isfinite(v)
								and v.is_integer()
								else str(v)
							)
						)
					)
					.astype(object)
				)
			else:
				raw = p[col]
				p[col] = pd.to_numeric(raw, errors="coerce").replace(
					[np.inf, -np.inf], np.nan
				)
				if (raw.notna() & p[col].isna()).any():
					raise ValueError(
						f"Non-numeric covariate {col}; declare it categorical if appropriate"
					)
		return p

	def fit(self, p):
		p = self.normalized(p)
		if p.isna().all().any():
			raise ValueError(
				"All-missing training covariate: " + ",".join(p.columns[p.isna().all()])
			)
		numeric = [c for c in self.cols if c not in self.categorical]
		spline = ["age"] if self.age_spline and "age" in numeric else []
		numeric = [c for c in numeric if c not in spline]
		pieces = []
		if numeric:
			pieces.append(
				(
					"num",
					make_pipeline(SimpleImputer(strategy="median"), StandardScaler()),
					numeric,
				)
			)
		if spline:
			pieces.append(
				(
					"age",
					make_pipeline(
						SimpleImputer(strategy="median"),
						SplineTransformer(n_knots=4, degree=3, include_bias=False),
						StandardScaler(),
					),
					spline,
				)
			)
		if self.categorical:
			pieces.append(
				(
					"cat",
					make_pipeline(
						SimpleImputer(strategy="most_frequent"),
						OneHotEncoder(
							handle_unknown="ignore",
							drop="first",
							sparse_output=self.sparse_output,
						),
					),
					self.categorical,
				)
			)
		self.transformer = ColumnTransformer(
			pieces, sparse_threshold=1.0 if self.sparse_output else 0.0
		)
		self.transformer.fit(p)
		self.levels = {c: set(p[c].dropna()) for c in self.categorical}
		return self

	def transform(self, p):
		if not self.cols:
			return np.empty((len(p), 0), dtype="float32")
		result = self.transformer.transform(self.normalized(p))
		# getattr keeps previously frozen bundles without this attribute usable.
		if getattr(self, "sparse_output", False):
			return sparse.csr_matrix(result, dtype="float32")
		return np.asarray(result, dtype="float32")

	def unknown_categories(self, p):
		p = self.normalized(p)
		return {
			c: int((p[c].notna() & ~p[c].isin(self.levels[c])).sum())
			for c in self.categorical
		}


class MolecularPreprocessor:
	def __init__(
		self, max_missing=0.2, residual_cols=(), categorical=(), transform="none"
	):
		self.max_missing = max_missing
		self.residual_cols = list(residual_cols)
		self.categorical = list(categorical)
		self.transform_name = transform

	def scale_transform(self, x):
		x = np.asarray(x, dtype="float32")
		x = np.where(np.isfinite(x), x, np.nan)
		if self.transform_name == "log1p":
			if np.any(x < 0):
				raise ValueError(
					"log1p requires nonnegative measured abundances; use --transform none for pretransformed data"
				)
			x = np.log1p(x)
		return x

	def fit(self, x, p, allowed=None):
		x = self.scale_transform(x)
		keep = (np.mean(~np.isfinite(x), axis=0) < self.max_missing) & (
			np.nanstd(x, axis=0) > 1e-8
		)
		if allowed is not None:
			# Freeze the feature set used for sample QC; do not move its denominator.
			keep = np.asarray(allowed, dtype=bool).copy()
		self.keep = np.flatnonzero(keep)
		if len(self.keep) < 3:
			raise ValueError("Fewer than three molecular features pass training QC")
		z = x[:, self.keep]
		if np.any(np.all(~np.isfinite(z), axis=0)):
			raise ValueError(
				"A retained assay has no observed training values after sample QC"
			)
		self.lower, self.upper = np.nanquantile(z, [0.005, 0.995], axis=0)
		z = np.clip(z, self.lower, self.upper)
		self.median = np.nanmedian(z, axis=0)
		z = np.where(np.isfinite(z), z, self.median)
		self.adjust = None
		if self.residual_cols:
			self.design = MetadataDesign(self.residual_cols, self.categorical).fit(p)
			c = self.design.transform(p)
			self.adjust = Ridge(alpha=1, solver="cholesky").fit(c, z)
			z = z - self.adjust.predict(c)
		self.mean, self.scale = z.mean(axis=0), z.std(axis=0)
		self.scale[self.scale < 1e-6] = 1
		return self

	def transform(self, x, p, batch_size=1024):
		# Plate one-hot coding can have thousands of columns. Predict residuals
		# in bounded row batches instead of constructing the whole cohort's
		# dense design (19.4 GiB float64 for the current metabolomics cohort).
		if batch_size < 1 or len(x) != len(p):
			raise ValueError("Invalid preprocessing batch size or metadata row count")
		result = np.empty((len(x), len(self.keep)), dtype="float32")
		observed = np.empty_like(result, dtype=bool)
		for start in range(0, len(x), batch_size):
			stop = min(start + batch_size, len(x))
			z = self.scale_transform(x[start:stop])[:, self.keep]
			measured = np.isfinite(z)
			z = np.clip(np.where(measured, z, self.median), self.lower, self.upper)
			if self.adjust is not None:
				z -= self.adjust.predict(self.design.transform(p.iloc[start:stop]))
			result[start:stop] = np.clip((z - self.mean) / self.scale, -10, 10)
			observed[start:stop] = measured
		return result, observed


# Reference
# Classical PCA geometry and quality-weighted real-person coverage selection.


class MolecularGeometry:
	def fit(self, x, importance, dimensions, seed):
		d = min(dimensions, len(x) - 2, x.shape[1] - 1)
		self.unsupervised = PCA(d, svd_solver="randomized", random_state=seed).fit(x)
		importance = np.maximum(importance, 0)
		self.has_supervised_signal = bool(importance.max() > 0)
		self.weight = (
			np.sqrt(0.02 + 0.98 * importance / max(importance.max(), 1e-12))
			if self.has_supervised_signal
			else np.ones_like(importance)
		).astype("float32")
		self.supervised = PCA(d, svd_solver="randomized", random_state=seed + 1).fit(
			x * self.weight
		)
		self.scale_u = np.sqrt(
			np.maximum(self.unsupervised.explained_variance_, 1e-6)
		).astype("float32")
		self.scale_s = np.sqrt(
			np.maximum(self.supervised.explained_variance_, 1e-6)
		).astype("float32")
		self.normalizer = np.sqrt(2 * d)
		return self

	def transform(self, x):
		return (
			np.c_[
				self.unsupervised.transform(x) / self.scale_u,
				self.supervised.transform(x * self.weight) / self.scale_s,
			].astype("float32")
			/ self.normalizer
		)


def diverse(z, pool, size, quality, seed):
	"""Quality-weighted farthest-first real exemplars, not synthetic centroids."""
	pool = np.asarray(pool)
	size = min(int(size), len(pool))
	if size <= 0:
		return np.array([], int)
	rng = np.random.default_rng(seed)
	q = np.asarray(quality)[pool]
	q = 0.25 + 0.75 * pd.Series(q).rank(pct=True).to_numpy()
	# Start near a representative centre; cap leverage of very distant points.
	cloud = z[pool]
	central = np.sum((cloud - np.median(cloud, axis=0)) ** 2, axis=1)
	density_cap = np.quantile(central, 0.95) + 1e-8
	q *= np.minimum(1, density_cap / np.maximum(central, 1e-8))
	first = int(np.argmax(q / (1 + central)))
	chosen, distance = [], np.full(len(pool), np.inf)
	available = np.ones(len(pool), bool)
	current = first
	for _ in range(size):
		chosen.append(pool[current])
		available[current] = False
		delta = np.sum((cloud - cloud[current]) ** 2, axis=1)
		distance = np.minimum(distance, delta)
		score = distance * q
		score[~available] = -1
		current = int(np.argmax(score + rng.uniform(0, 1e-12, len(score))))
	return np.asarray(chosen, int)


# Attention
# Q/K/V attention with inspectable heads, and label-preserving cross attention.
#
# Self attention has learned Q,K,V and output projections, residuals, pre-LN,
# FFN and dropout. Cross attention uses learned Q,K but observed donor Y as V;
# convex head mixing preserves a literal decomposition into real people's labels.
#


class ContextBlock(nn.Module):
	def __init__(self, width, heads, dropout, mode="learned"):
		super().__init__()
		self.heads, self.width, self.mode = heads, width, mode
		self.norm1, self.norm2 = nn.LayerNorm(width), nn.LayerNorm(width)
		self.qkv = nn.Linear(width, 3 * width)
		self.out = nn.Linear(width, width)
		self.ffn = nn.Sequential(
			nn.Linear(width, 4 * width),
			nn.GELU(),
			nn.Dropout(dropout),
			nn.Linear(4 * width, width),
		)
		self.dropout = nn.Dropout(dropout)
		self.attention_dropout = dropout
		self.intervention = (
			None  # evaluation-only intervention; never a fitted parameter
		)

	def forward(self, x, capture=False):
		b, n, _ = x.shape
		q, k, v = (
			self.qkv(self.norm1(x))
			.view(b, n, 3, self.heads, self.width // self.heads)
			.permute(2, 0, 3, 1, 4)
		)
		mode = self.intervention or self.mode
		weights = None
		if mode == "learned" and not capture:
			context = F.scaled_dot_product_attention(
				q, k, v, dropout_p=self.attention_dropout if self.training else 0.0
			)
		else:
			if mode == "uniform":
				weights = x.new_full((b, self.heads, n, n), 1 / n)
			elif mode == "identity":
				weights = torch.eye(n, device=x.device, dtype=x.dtype)[
					None, None
				].expand(b, self.heads, -1, -1)
			else:
				weights = torch.softmax(
					q @ k.transpose(-2, -1) / math.sqrt(q.shape[-1]), dim=-1
				)
			context = (
				F.dropout(weights, p=self.attention_dropout, training=self.training) @ v
			)
		x = x + self.dropout(
			self.out(context.transpose(1, 2).reshape(b, n, self.width))
		)
		x = x + self.dropout(self.ffn(self.norm2(x)))
		return x, weights


class ContextStack(nn.Module):
	def __init__(self, width, heads, layers, dropout, mode="learned"):
		super().__init__()
		self.layers = nn.ModuleList(
			[ContextBlock(width, heads, dropout, mode) for _ in range(layers)]
		)
		self.norm = nn.LayerNorm(width)

	def forward(self, x, capture=False):
		maps = []
		for layer in self.layers:
			x, weights = layer(x, capture)
			if capture:
				maps.append(weights)
		return self.norm(x), torch.stack(maps, dim=1) if capture else None


def numpy_head_probabilities(logits):
	"""Softmax over donors; all forbidden returns zeros, not NaNs."""
	maximum = np.max(logits, axis=1, keepdims=True)
	maximum = np.where(np.isfinite(maximum), maximum, 0)
	values = np.exp(logits.astype("float64") - maximum)
	total = values.sum(1, keepdims=True)
	return np.divide(values, total, out=np.zeros_like(values), where=total > 0)


def numpy_mixture_scores(logits, gate):
	probability = np.sum(numpy_head_probabilities(logits) * gate[:, None, :], axis=-1)
	with np.errstate(divide="ignore"):
		return np.log(probability)


def torch_mixture_scores(logits, gate):
	return torch.logsumexp(
		torch.log_softmax(logits, dim=1) + gate.clamp_min(1e-30).log()[:, None, :],
		dim=-1,
	)


# Neural
# Hierarchical molecular Transformer and differentiable person-to-person copying.
#
# Outcome values never enter the encoder or query/key networks. During training,
# all donors from the query's inner fold (and therefore its family) are masked.
# Cached donor keys are refreshed each epoch, detached from the current gradient.
#


def device_for(request):
	if request == "auto":
		return "cuda" if torch.cuda.is_available() else "cpu"
	if request == "cuda" and not torch.cuda.is_available():
		raise RuntimeError(
			"CUDA requested but unavailable. Run --check-device or use --device cpu."
		)
	return request


def token_partition(x, features, count, seed, module_file=""):
	"""Every retained assay belongs to exactly one token. No Y is used."""
	labels = np.full(len(features), -1, int)
	names = []
	if module_file:
		tab = pd.read_csv(module_file, sep="\t", dtype=str)
		if not {"feature", "module"} <= set(tab):
			raise ValueError("--module-file requires feature and module TSV columns")
		# Deterministic first alphabetical membership for overlapping pathways.
		lookup = {f: i for i, f in enumerate(features)}
		for name, sub in tab.sort_values(["module", "feature"]).groupby("module"):
			ix = [
				lookup[f]
				for f in sub.feature.unique()
				if f in lookup and labels[lookup[f]] < 0
			]
			if ix:
				labels[ix] = len(names)
				names.append(str(name))
	left = np.flatnonzero(labels < 0)
	if len(left):
		k = min(max(1, count - len(names)), len(left))
		if k == 1:
			cluster = np.zeros(len(left), int)
		else:
			d = min(24, len(x) - 1, len(left))
			load = (
				PCA(d, svd_solver="randomized", random_state=seed)
				.fit(x[:, left])
				.components_.T
			)
			load /= np.maximum(np.linalg.norm(load, axis=1, keepdims=True), 1e-8)
			cluster = MiniBatchKMeans(
				k, random_state=seed, n_init=5, batch_size=512
			).fit_predict(load)
		for group in np.unique(cluster):
			labels[left[cluster == group]] = len(names)
			names.append(f"data_token_{len(names) + 1:03d}")
	return labels, names


class RowEncoder(nn.Module):
	def __init__(
		self,
		membership,
		width=64,
		heads=4,
		layers=2,
		dropout=0.1,
		architecture="transformer",
		cross_mode="qkv",
	):
		super().__init__()
		membership = np.asarray(membership, int)
		self.n_features, self.n_tokens = len(membership), int(membership.max()) + 1
		self.width, self.heads, self.architecture = width, heads, architecture
		self.cross_mode = cross_mode
		self.register_buffer("membership", torch.tensor(membership, dtype=torch.long))
		counts = np.bincount(membership, minlength=self.n_tokens).astype("float32")
		self.register_buffer("counts", torch.tensor(counts))
		# Learned feature identities distinguish every assay within its module.
		self.value_embedding = nn.Parameter(torch.randn(self.n_features, width) * 0.05)
		self.missing_embedding = nn.Parameter(
			torch.randn(self.n_features, width) * 0.02
		)
		self.token_embedding = nn.Parameter(torch.randn(self.n_tokens, width) * 0.02)
		self.cls = nn.Parameter(torch.randn(1, 1, width) * 0.02)
		self.input_norm = nn.LayerNorm(width)
		if architecture in ["transformer", "uniform"]:
			self.context = ContextStack(
				width,
				heads,
				layers,
				dropout,
				"uniform" if architecture == "uniform" else "learned",
			)
		else:
			self.context = nn.Sequential(
				nn.Linear(self.n_tokens * width, width * 4),
				nn.GELU(),
				nn.Dropout(dropout),
				nn.Linear(width * 4, width),
				nn.LayerNorm(width),
			)
		self.project = nn.Sequential(
			nn.Linear(width, width), nn.GELU(), nn.Linear(width, width)
		)
		self.direct = nn.Linear(width, 1)
		self.decoder_weight = nn.Parameter(torch.randn(self.n_features, width) * 0.05)
		self.decoder_bias = nn.Parameter(torch.zeros(self.n_features))
		self.gate = nn.Linear(width, heads)
		self.temperature = nn.Parameter(
			torch.full((heads,), -1.0 if cross_mode == "qkv" else -3.0)
		)
		self.query = nn.Linear(width, width, bias=False)
		self.key = nn.Linear(width, width, bias=False)
		# An asymmetric learned Q/K projection; sqrt(width) restores unit scale
		# after the normalized row embedding, before scaled dot product.
		with torch.no_grad():
			self.query.weight.copy_(torch.eye(width) + torch.randn(width, width) * 0.01)
			self.key.weight.copy_(torch.eye(width) + torch.randn(width, width) * 0.01)

	def forward(self, x, observed, decode=False, return_attention=False):
		# No feature value hidden by observed=False can enter the encoder.
		x = torch.where(observed, x, torch.zeros_like(x))
		token = x.new_zeros((len(x), self.n_tokens, self.width))
		token.index_add_(1, self.membership, x.unsqueeze(-1) * self.value_embedding)
		token.index_add_(
			1, self.membership, (~observed).unsqueeze(-1) * self.missing_embedding
		)
		token = token / self.counts.sqrt()[None, :, None] + self.token_embedding
		token = self.input_norm(token)
		maps = None
		if self.architecture in ["transformer", "uniform"]:
			context, maps = self.context(
				torch.cat([self.cls.expand(len(x), -1, -1), token], dim=1),
				return_attention,
			)
			row, contextual = context[:, 0], context[:, 1:]
		else:
			row = self.context(token.flatten(1))
			contextual = token + row[:, None, :]
		z = F.normalize(self.project(row), dim=-1)
		reconstruction = None
		if decode:
			reconstruction = (contextual[:, self.membership] * self.decoder_weight).sum(
				-1
			) + self.decoder_bias
		result = (z, self.direct(row).squeeze(-1), reconstruction)
		return (*result, maps) if return_attention else result

	def head_scores(self, q, bank):
		h, d = self.heads, self.width // self.heads
		scale = self.width**0.5
		qh = self.query(q * scale).reshape(-1, h, d)
		kh = self.key(bank * scale).reshape(-1, h, d)
		temperature = 0.02 + 1.98 * torch.sigmoid(self.temperature)
		return torch.einsum("bhd,rhd->brh", qh, kh) / (d**0.5 * temperature)

	def scores(self, q, bank):
		if self.cross_mode == "qkv":
			return torch_mixture_scores(
				self.head_scores(q, bank), torch.softmax(self.gate(q), -1)
			)
		h, d = self.heads, self.width // self.heads
		qh, bh = q.reshape(-1, h, d), bank.reshape(-1, h, d)
		distances = (
			qh.square().sum(-1)[:, None, :]
			+ bh.square().sum(-1)[None, :, :]
			- 2 * torch.einsum("bhd,rhd->brh", qh, bh)
		).clamp_min(0)
		gate = torch.softmax(self.gate(q), -1)
		temperature = 0.02 + 1.98 * torch.sigmoid(self.temperature)
		return -(distances * gate[:, None, :] / temperature).sum(-1)


def torch_copy(model, q, bank, by, bw, forbidden, k, prior, strength):
	if torch.any((~forbidden).sum(1) < k):
		raise ValueError("Insufficient independent neural donors")
	q, bank = q.float(), bank.float()
	if model.cross_mode == "qkv":
		logits = model.head_scores(q, bank).masked_fill(
			forbidden[:, :, None], -torch.inf
		)
		gate = torch.softmax(model.gate(q), -1)
		score = torch_mixture_scores(logits, gate)
		jj = torch.argsort(score, descending=True, stable=True)[:, :k]
		chosen = logits.gather(1, jj[:, :, None].expand(-1, -1, model.heads))
		# Label donors only; IPCW is a positive measure over keys, not part of Q/K.
		head_weights = torch.softmax(chosen + bw[jj].log()[:, :, None], dim=1)
		weights = (head_weights * gate[:, None, :]).sum(-1)
	else:
		score = model.scores(q, bank).masked_fill(forbidden, -torch.inf)
		jj = torch.argsort(score, descending=True, stable=True)[:, :k]
		weights = torch.softmax(score.gather(1, jj), dim=1) * bw[jj]
		weights = weights / weights.sum(1, keepdim=True)
	ess = weights.square().sum(1).reciprocal()
	p = (ess * (weights * by[jj]).sum(1) + strength * prior) / (ess + strength)
	return p, jj, weights


@torch.no_grad()
def encode(model, x, observed, device="cpu", batch_size=256, decode=False):
	model.to(device).eval()
	zs, ps, rs = [], [], []
	for start in range(0, len(x), batch_size):
		xx = torch.as_tensor(
			np.ascontiguousarray(x[start : start + batch_size]),
			dtype=torch.float32,
			device=device,
		)
		mm = torch.as_tensor(
			np.ascontiguousarray(observed[start : start + batch_size]),
			dtype=torch.bool,
			device=device,
		)
		z, p, r = model(xx, mm, decode)
		zs.append(z.cpu().numpy())
		ps.append(torch.sigmoid(p).cpu().numpy())
		if decode:
			rs.append(r.cpu().numpy())
	return (
		np.concatenate(zs),
		np.concatenate(ps),
		np.concatenate(rs) if decode else None,
	)


def train_encoder(
	x,
	observed,
	y,
	w,
	ids,
	groups,
	membership,
	a,
	output,
	architecture="transformer",
	objective="retrieval",
	pretrain=True,
	cross_mode="qkv",
):
	"""Only build rows enter here; an inner fold is reserved for early stopping.

	The OOF quality score used elsewhere is independently cross-fitted. This
	neural model's training predictions are never relabelled as OOF evidence.
	"""
	output = Path(output)
	output.mkdir(parents=True, exist_ok=True)
	dev = device_for(a.device)
	torch.set_num_threads(a.cores)
	torch.manual_seed(a.seed)
	np.random.seed(a.seed)
	if torch.cuda.is_available():
		torch.cuda.manual_seed_all(a.seed)
	torch.use_deterministic_algorithms(a.deterministic, warn_only=False)
	cv = (
		StratifiedGroupKFold(5, shuffle=True, random_state=a.seed + 43)
		if groups is not None
		else StratifiedKFold(5, shuffle=True, random_state=a.seed + 43)
	)
	strata = np.where(w > 0, y, 2)
	tr, va = next(cv.split(x, strata, groups))
	folds = np.full(len(x), -1, int)
	cv2 = (
		StratifiedGroupKFold(a.folds, shuffle=True, random_state=a.seed + 44)
		if groups is not None
		else StratifiedKFold(a.folds, shuffle=True, random_state=a.seed + 44)
	)
	for fold, (_, ix) in enumerate(
		cv2.split(x[tr], strata[tr], None if groups is None else groups[tr])
	):
		folds[tr[ix]] = fold
	pd.DataFrame(
		{
			"eid": ids,
			"inner_role": np.where(folds < 0, "early_stop", "train"),
			"donor_exclusion_fold": folds,
		}
	).to_csv(output / "inner_split.csv", index=False)
	donor = tr[w[tr] > 0]
	if min(np.bincount(y[donor], minlength=2)) < 5:
		raise ValueError("Neural inner training needs >=5 known outcomes in each class")
	k = min(a.neural_k, min(np.sum(folds[donor] != f) for f in np.unique(folds[tr])))
	if k < 2:
		raise ValueError("Too few fold-disjoint donors")
	model = RowEncoder(
		membership, a.width, a.heads, a.layers, a.dropout, architecture, cross_mode
	).to(dev)
	prior = float(np.average(y[donor], weights=w[donor]))
	with torch.no_grad():
		model.direct.bias.fill_(np.log(prior / (1 - prior)))
	optimizer = torch.optim.AdamW(
		model.parameters(), lr=a.learning_rate, weight_decay=a.weight_decay
	)
	amp = a.amp and dev == "cuda"
	scaler = torch.amp.GradScaler("cuda", enabled=amp)
	rng = np.random.default_rng(a.seed + 47)
	logs, completed, best, bad, best_state = [], {}, np.inf, 0, None
	checkpoint = output / "checkpoint.pt"
	ssl_state = None
	if a.resume and checkpoint.exists():
		saved = torch.load(checkpoint, map_location="cpu", weights_only=False)
		model.load_state_dict(saved["state"])
		optimizer.load_state_dict(saved["optimizer"])
		for state in optimizer.state.values():
			for key, value in state.items():
				if torch.is_tensor(value):
					state[key] = value.to(dev)
		scaler.load_state_dict(saved["scaler"])
		logs, completed, best, bad = (
			saved["logs"],
			saved["completed"],
			saved["best"],
			saved["bad"],
		)
		best_state, ssl_state = saved["best_state"], saved["ssl_state"]
		rng.bit_generator.state = saved["numpy_rng"]
		torch.set_rng_state(saved["torch_rng"])
		if dev == "cuda" and saved["cuda_rng"] is not None:
			torch.cuda.set_rng_state_all(saved["cuda_rng"])
		log("RESUME", output.name, str(completed))
	val_mask = observed[va] & (
		np.random.default_rng(a.seed + 49).random(observed[va].shape) < a.mask_fraction
	)
	by = torch.tensor(y[donor], dtype=torch.float32, device=dev)
	bw = torch.tensor(w[donor], dtype=torch.float32, device=dev)
	bank_fold = torch.tensor(folds[donor], device=dev)
	start_time = time.monotonic()
	for phase, epochs in [
		("pretrain", a.pretrain_epochs if pretrain else 0),
		("joint", a.epochs),
	]:
		if phase == "joint" and completed.get("joint", 0) == 0:
			best, bad, best_state = np.inf, 0, None
		if completed.get(phase, 0) >= epochs or (
			phase == "joint" and bad >= a.patience
		):
			continue
		for epoch in range(completed.get(phase, 0), epochs):
			epoch_started = time.monotonic()
			if dev == "cuda":
				torch.cuda.reset_peak_memory_stats()
			if phase == "joint":
				bz = torch.tensor(
					encode(model, x[donor], observed[donor], dev, a.batch_size)[0],
					device=dev,
				)
			model.train()
			totals = np.zeros(4)
			nb = 0
			order = rng.permutation(tr)
			for begin in range(0, len(order), a.batch_size):
				ix = order[begin : begin + a.batch_size]
				xx = torch.tensor(x[ix], dtype=torch.float32, device=dev)
				mm = torch.tensor(observed[ix], dtype=torch.bool, device=dev)
				hidden = mm & (torch.rand(mm.shape, device=dev) < a.mask_fraction)
				optimizer.zero_grad(set_to_none=True)
				with torch.autocast(device_type=dev, dtype=torch.float16, enabled=amp):
					z, direct, reconstruction = model(xx, mm & ~hidden, decode=True)
					rec = (
						(reconstruction.float() - xx).square()[hidden].mean()
						if hidden.any()
						else reconstruction.sum() * 0
					)
					direct_loss = rec * 0
					retrieval_loss = rec * 0
					total = (
						rec if phase == "pretrain" else a.reconstruction_weight * rec
					)
					if phase == "joint":
						yy = torch.tensor(y[ix], dtype=torch.float32, device=dev)
						ww = torch.tensor(w[ix], dtype=torch.float32, device=dev)
						direct_loss = (
							F.binary_cross_entropy_with_logits(
								direct.float(), yy, reduction="none"
							)
							* ww
						).sum() / ww.sum().clamp_min(1)
						if objective == "retrieval":
							forbidden = (
								torch.tensor(folds[ix], device=dev)[:, None]
								== bank_fold[None, :]
							)
							with torch.autocast(device_type=dev, enabled=False):
								borrowed, _, _ = torch_copy(
									model,
									z.float(),
									bz.float(),
									by.float(),
									bw.float(),
									forbidden,
									k,
									prior,
									a.train_prior,
								)
								retrieval_loss = (
									F.binary_cross_entropy(
										borrowed.float().clamp(1e-5, 1 - 1e-5),
										yy.float(),
										reduction="none",
									)
									* ww.float()
								).sum() / ww.sum().clamp_min(1)
							total = (
								total + retrieval_loss + a.direct_weight * direct_loss
							)
						else:
							total = total + direct_loss
				if not torch.isfinite(total):
					raise FloatingPointError(
						"Nonfinite neural loss; inspect input/learning rate"
					)
				scaler.scale(total).backward()
				scaler.unscale_(optimizer)
				nn.utils.clip_grad_norm_(model.parameters(), a.gradient_clip)
				scaler.step(optimizer)
				scaler.update()
				totals += [
					float(v.detach()) for v in [total, rec, direct_loss, retrieval_loss]
				]
				nb += 1
				if time.monotonic() - start_time > 55:
					log(
						"RUNNING", output.name, f"{phase} epoch={epoch + 1}, batch={nb}"
					)
					start_time = time.monotonic()
			model.eval()
			with torch.no_grad():
				zv, pv, rv = encode(
					model, x[va], observed[va] & ~val_mask, dev, a.batch_size, True
				)
				val_rec = float(np.mean((rv[val_mask] - x[va][val_mask]) ** 2))
				if phase == "pretrain":
					score = val_rec
					val_ll = np.nan
				else:
					if objective == "retrieval":
						bz = torch.tensor(
							encode(model, x[donor], observed[donor], dev, a.batch_size)[
								0
							],
							device=dev,
						)
						vals = []
						for begin in range(0, len(va), a.batch_size):
							q = torch.tensor(
								zv[begin : begin + a.batch_size], device=dev
							)
							pred, _, _ = torch_copy(
								model,
								q,
								bz,
								by,
								bw,
								torch.zeros(
									(len(q), len(donor)), dtype=torch.bool, device=dev
								),
								k,
								prior,
								a.train_prior,
							)
							vals.append(pred.cpu().numpy())
						pv = np.concatenate(vals)
					pv = np.clip(pv, 1e-6, 1 - 1e-6)
					val_ll = float(
						np.average(
							-(y[va] * np.log(pv) + (1 - y[va]) * np.log1p(-pv)),
							weights=w[va],
						)
					)
					score = val_ll
			improved = score < best - 1e-6
			if improved:
				best, bad = score, 0
				best_state = {
					k: v.detach().cpu().clone() for k, v in model.state_dict().items()
				}
			else:
				bad += 1
			completed[phase] = epoch + 1
			row = dict(
				phase=phase,
				epoch=epoch + 1,
				train_loss=totals[0] / nb,
				train_masked_mse=totals[1] / nb,
				train_direct_loss=totals[2] / nb,
				train_retrieval_loss=totals[3] / nb,
				early_stop_logloss=val_ll,
				early_stop_masked_mse=val_rec,
				best=improved,
				seconds=time.monotonic() - epoch_started,
				cuda_peak_allocated_mb=(
					torch.cuda.max_memory_allocated() / 2**20 if dev == "cuda" else None
				),
			)
			logs.append(row)
			pd.DataFrame(logs).to_csv(output / "learning_curve.csv", index=False)
			if phase == "pretrain" and epoch + 1 == epochs:
				# Start supervised training from the best masked-reconstruction epoch.
				model.load_state_dict(best_state)
				ssl_state = copy.deepcopy(best_state)
			saved = dict(
				state={k: v.detach().cpu() for k, v in model.state_dict().items()},
				optimizer=optimizer.state_dict(),
				scaler=scaler.state_dict(),
				logs=logs,
				completed=completed,
				best=best,
				bad=bad,
				best_state=best_state,
				ssl_state=ssl_state,
				numpy_rng=rng.bit_generator.state,
				torch_rng=torch.get_rng_state(),
				cuda_rng=torch.cuda.get_rng_state_all() if dev == "cuda" else None,
			)
			tmp = checkpoint.with_suffix(".tmp")
			torch.save(saved, tmp)
			tmp.replace(checkpoint)
			log(
				"EPOCH",
				output.name,
				f"{phase} {epoch + 1}/{epochs}; masked_MSE={val_rec:.4f}; logloss={val_ll:.4f}; seconds={row['seconds']:.1f}",
			)
			if phase == "joint" and bad >= a.patience:
				break
	if best_state is None:
		raise RuntimeError("No neural checkpoint selected")
	model.load_state_dict(best_state)
	model.cpu().eval()
	ssl = copy.deepcopy(model) if ssl_state is not None else None
	if ssl is not None:
		ssl.load_state_dict(ssl_state)
	dump(
		output / "training_summary.json",
		dict(
			device=dev,
			architecture=architecture,
			objective=objective,
			parameters=sum(p.numel() for p in model.parameters()),
			build=len(x),
			gradient_rows=len(tr),
			early_stop_rows=len(va),
			known_donors=len(donor),
			donor_fold_exclusion=True,
			bank_refresh="each_epoch; detached row embeddings; Q/K projections remain differentiable",
			cross_mode=cross_mode,
			completed_epochs=completed,
			best_early_stop_score=best,
			torch_version=torch.__version__,
			mixed_precision=amp,
		),
	)
	torch.save(model.state_dict(), output / "selected_weights.pt")
	return model, ssl


# Borrowing
# Auditable COPY and risk borrowing. Y enters values, never matching scores.


def stable_seed(eid, salt=0):
	return int.from_bytes(
		hashlib.blake2b((str(eid) + ":" + str(salt)).encode(), digest_size=8).digest(),
		"little",
	)


def select_references(z, quality, size, mode, seed, utility=None, balanced=False):
	"""Actual participants; fit scores are never interpreted as data validity."""
	known = quality.known_label.to_numpy(bool)
	y = quality.horizon_label.to_numpy(int)
	pool = np.flatnonzero(known)
	if mode == "all":
		return pool
	if mode == "topfit":
		return pool[
			np.argsort(quality.OOF_logloss.to_numpy()[pool], kind="stable")[:size]
		]
	if mode in ["reliable", "utility"]:
		pool = np.flatnonzero(known & quality.reliable_candidate.to_numpy(bool))
	size = min(size, len(pool))
	if size < 2:
		return np.array([], int)
	fraction = 0.5 if balanced else float(y[known].mean())
	counts = [size - int(round(size * fraction)), int(round(size * fraction))]
	if len(np.unique(y[pool])) == 2:
		counts[1] = max(1, min(size - 1, counts[1]))
		counts[0] = size - counts[1]
	counts = [min(counts[v], int(np.sum(y[pool] == v))) for v in [0, 1]]
	for v in [0, 1]:
		counts[v] += min(size - sum(counts), int(np.sum(y[pool] == v)) - counts[v])
	gain = quality.fit_gain.fillna(-1e6).to_numpy()
	result = []
	for label, count in enumerate(counts):
		candidates = pool[y[pool] == label]
		if not count:
			continue
		if mode == "random":
			take = np.random.default_rng(seed + label).choice(
				candidates, count, replace=False
			)
		elif mode == "stratified_fit":
			take = candidates[
				np.argsort(quality.OOF_logloss.to_numpy()[candidates], kind="stable")[
					:count
				]
			]
		else:
			score = np.ones(len(y)) if mode == "diversity" else gain.copy()
			if mode == "utility" and utility is not None:
				# Rank-based blend keeps scale/large-exposure donors from dominating.
				score = (
					pd.Series(score).rank(pct=True).to_numpy()
					+ 2 * pd.Series(utility).rank(pct=True).to_numpy()
				)
			take = diverse(z, candidates, count, score, seed + label)
		result.extend(take)
	return np.sort(result)


def risk_from_weights(weights, labels, prior, strength):
	sum_sq = np.sum(weights**2, axis=1)
	ess = np.divide(1.0, sum_sq, out=np.zeros_like(sum_sq), where=sum_sq > 0)
	prior_weight = np.divide(
		strength, ess + strength, out=np.ones_like(ess), where=(ess + strength) > 0
	)
	coefficients = weights * (1 - prior_weight[:, None])
	p = np.sum(coefficients * labels, axis=1) + prior_weight * prior
	return p, coefficients, prior_weight, ess


class ReferenceBank:
	def __init__(
		self,
		z,
		x,
		observed,
		ids,
		groups,
		y,
		w,
		indices,
		encoder=None,
		inclusion_correction=False,
		seed=2026,
	):
		self.indices = np.asarray(indices, int)
		self.z, self.x, self.observed = z[indices], x[indices], observed[indices]
		self.ids, self.groups = (
			ids[indices],
			(None if groups is None else groups[indices]),
		)
		self.y, self.w = y[indices].astype(float), w[indices].astype(float)
		self.prior = float(np.average(y, weights=w))
		self.seed = seed
		self.gate_weight = self.gate_bias = self.temperature = None
		self.cross_mode = "metric"
		self.query_weight = self.key_weight = None
		if encoder is not None:
			self.cross_mode = encoder.cross_mode
			self.query_weight = encoder.query.weight.detach().cpu().numpy().copy()
			self.key_weight = encoder.key.weight.detach().cpu().numpy().copy()
			self.gate_weight = encoder.gate.weight.detach().cpu().numpy().copy()
			self.gate_bias = encoder.gate.bias.detach().cpu().numpy().copy()
			t = encoder.temperature.detach().cpu().numpy()
			self.temperature = 0.02 + 1.98 / (1 + np.exp(-t))
		self.correction = np.ones(len(indices))
		if inclusion_correction:
			# Only a class sampling correction; does not correct within-class selection.
			for label in [0, 1]:
				ix = self.y == label
				if ix.any():
					self.correction[ix] = np.mean(y[w > 0] == label) / np.mean(ix)
		permutation = np.random.default_rng(seed + 841).permutation(len(indices))
		self.permuted_y, self.permuted_w = self.y[permutation], self.w[permutation]

	def head_scores(self, q):
		h = len(self.temperature)
		d = q.shape[1] // h
		scale = q.shape[1] ** 0.5
		qq = ((q * scale) @ self.query_weight.T).reshape(-1, h, d)
		kk = ((self.z * scale) @ self.key_weight.T).reshape(-1, h, d)
		return np.einsum("bhd,rhd->brh", qq, kk, optimize=True) / (
			d**0.5 * self.temperature
		)

	def head_gate(self, q):
		return softmax(
			(q @ self.gate_weight.T + self.gate_bias).astype("float64"), axis=1
		)

	def scores(self, q, kernel="attention"):
		if self.gate_weight is None or kernel == "euclidean":
			return -np.maximum(
				np.sum(q * q, 1)[:, None]
				+ np.sum(self.z * self.z, 1)[None, :]
				- 2 * q @ self.z.T,
				0,
			)
		if self.cross_mode == "qkv":
			return numpy_mixture_scores(self.head_scores(q), self.head_gate(q))
		h = len(self.temperature)
		d = q.shape[1] // h
		qq, zz = q.reshape(-1, h, d), self.z.reshape(-1, h, d)
		distances = np.maximum(
			np.sum(qq * qq, 2)[:, None, :]
			+ np.sum(zz * zz, 2)[None, :, :]
			- 2 * np.einsum("bhd,rhd->brh", qq, zz, optimize=True),
			0,
		)
		gate = softmax(q @ self.gate_weight.T + self.gate_bias, axis=1)
		return -np.sum(distances * gate[:, None, :] / self.temperature, axis=2)

	def match(
		self,
		q,
		ids,
		groups,
		k=20,
		temperature=1.0,
		strength=2.0,
		kernel="attention",
		random_candidates=False,
		permuted=False,
		literal=False,
		batch_size=256,
	):
		k = min(int(k), len(self.ids))
		if k < 1:
			raise ValueError("Empty reference panel")
		all_j, all_w, all_distance = [], [], []
		all_heads, all_gates = [], []
		multihead = self.cross_mode == "qkv" and kernel in ["attention", "equal_heads"]
		for begin in range(0, len(q), batch_size):
			query = q[begin : begin + batch_size]
			forbidden = ids[begin : begin + batch_size, None] == self.ids[None, :]
			if groups is not None and self.groups is not None:
				forbidden |= (
					groups[begin : begin + batch_size, None] == self.groups[None, :]
				)
			if multihead:
				logits = self.head_scores(query) / temperature
				logits[forbidden] = -np.inf
				gate = self.head_gate(query)
				if kernel == "equal_heads":
					gate[:] = 1 / gate.shape[1]
				score = numpy_mixture_scores(logits, gate)
			else:
				score = self.scores(query, kernel) / temperature
				score[forbidden] = -np.inf
			if random_candidates:
				jj = []
				for i, eid in enumerate(ids[begin : begin + batch_size]):
					valid = np.flatnonzero(~forbidden[i])
					selected = np.random.default_rng(
						stable_seed(eid, self.seed)
					).choice(valid, min(k, len(valid)), replace=False)
					padding = np.flatnonzero(forbidden[i])[: k - len(selected)]
					jj.append(np.r_[selected, padding])
				jj = np.asarray(jj)
			else:
				jj = np.argsort(-score, axis=1, kind="stable")[:, :k]
			values = np.take_along_axis(score, jj, axis=1)
			ww = self.permuted_w if permuted else self.w
			measure = ww[jj] * self.correction[jj]
			if multihead:
				chosen = np.take_along_axis(logits, jj[:, :, None], axis=1)
				with np.errstate(divide="ignore"):
					chosen = chosen + np.log(measure)[:, :, None]
				head_weights = numpy_head_probabilities(chosen)
				weights = np.sum(head_weights * gate[:, None, :], axis=2)
				all_heads.append(head_weights.astype("float32"))
				all_gates.append(gate)
				# A geometric support diagnostic independent of logit offsets.
				distance = np.sum((query - self.z[jj[:, 0]]) ** 2, axis=1)
				distance[weights.sum(1) == 0] = np.inf
			else:
				finite = np.isfinite(values)
				maximum = values.max(1, keepdims=True)
				maximum = np.where(np.isfinite(maximum), maximum, 0)
				weights = (
					np.where(finite, np.exp(values.astype("float64") - maximum), 0)
					* measure
				)
				total = weights.sum(1, keepdims=True)
				np.divide(weights, total, out=weights, where=total > 0)
				distance = -values.max(1) * temperature
			all_j.append(jj)
			all_w.append(weights)
			all_distance.append(distance)
		jj, weights = np.concatenate(all_j), np.concatenate(all_w)
		labels = (self.permuted_y if permuted else self.y)[jj]
		if literal:
			if k != 1:
				raise ValueError("Literal COPY needs k=1")
			strength = 0
			weights = (weights > 0).astype(float)  # literal binary copying is exact
		p, coefficient, pw, ess = risk_from_weights(
			weights, labels, self.prior, strength
		)
		if literal:
			p[ess == 0] = np.nan  # no donor's Y exists to COPY; do not invent one
		diagnostics = dict(
			jj=jj,
			weights=weights,
			coefficient=coefficient,
			prior_weight=pw,
			ESS=ess,
			nearest_distance=np.concatenate(all_distance),
			labels=labels,
			entropy=-np.sum(weights * np.log(np.maximum(weights, 1e-30)), axis=1),
			case_weight=np.sum(weights * labels, axis=1),
		)
		if multihead:
			diagnostics["head_weights"] = np.concatenate(all_heads)
			diagnostics["head_gate"] = np.concatenate(all_gates)
			diagnostics["head_risk"] = np.sum(
				diagnostics["head_weights"] * labels[:, :, None], axis=1
			)
		return p, diagnostics

	def reconstruct(self, detail):
		jj, ww = detail["jj"], detail["weights"]
		reconstructed = np.zeros((len(jj), self.x.shape[1]), dtype="float32")
		denominator = np.zeros_like(reconstructed)
		# A missing donor assay is never silently treated as a measured value.
		for rank in range(jj.shape[1]):
			weight = ww[:, rank : rank + 1] * self.observed[jj[:, rank]]
			reconstructed += weight * self.x[jj[:, rank]]
			denominator += weight
		np.divide(reconstructed, denominator, out=reconstructed, where=denominator > 0)
		return reconstructed, denominator > 0

	def describe(self, x, observed, p, detail):
		reconstruction, covered = self.reconstruct(detail)
		valid = observed & covered
		rmse = np.sqrt(
			np.sum(np.where(valid, (x - reconstruction) ** 2, 0), axis=1)
			/ np.maximum(valid.sum(1), 1)
		)
		ww, labels, jj = detail["weights"], detail["labels"], detail["jj"]
		variance = np.sum(ww * (labels - detail["case_weight"][:, None]) ** 2, axis=1)
		most_weighted = np.argmax(ww, axis=1)
		return pd.DataFrame(
			dict(
				nearest_reference=np.where(detail["ESS"] > 0, self.ids[jj[:, 0]], ""),
				highest_weight_reference=np.where(
					detail["ESS"] > 0,
					self.ids[jj[np.arange(len(jj)), most_weighted]],
					"",
				),
				nearest_distance=detail["nearest_distance"],
				reconstruction_RMSE=rmse,
				reference_ESS=detail["ESS"],
				attention_entropy=detail["entropy"],
				independent_references=np.sum(ww > 0, axis=1),
				normalized_entropy=detail["entropy"] / max(np.log(jj.shape[1]), 1),
				reference_disagreement=np.sqrt(variance),
				case_weight=detail["case_weight"],
				case_references=np.sum(labels * (ww > 0), axis=1),
				prior_weight=detail["prior_weight"],
				reconstructed_assay_fraction=covered.mean(1),
			)
		)


def tune_bank(
	bank,
	z,
	ids,
	groups,
	y,
	w,
	ks,
	temperatures,
	priors,
	name,
	kernel="attention",
	literal=False,
):
	rows, best = [], None
	for k in ks:
		if k > len(bank.ids):
			continue
		for temp in temperatures:
			_, cached = bank.match(
				z,
				ids,
				groups,
				k=k,
				temperature=temp,
				strength=0.0,
				kernel=kernel,
				literal=literal,
			)
			for strength in priors:
				spec = dict(
					k=k,
					temperature=temp,
					strength=strength,
					kernel=kernel,
					literal=literal,
				)
				p = risk_from_weights(
					cached["weights"], cached["labels"], bank.prior, strength
				)[0]
				ll = float(np.average(loss(y, p), weights=w))
				rows.append(dict(model=name, **spec, tune_logloss=ll))
				if best is None or ll < best[0]:
					best = ll, spec
	if best is None:
		raise ValueError(f"No usable neighbor setting for {name}")
	return best[1], rows


def reference_utility(bank, z, ids, groups, y, w, k=32):
	"""Tuning-set deletion utility, NOT an OOF subject-fit score or causality."""
	p, d = bank.match(z, ids, groups, k=min(k, len(bank.ids)), strength=2)
	summed, exposure = np.zeros(len(bank.ids)), np.zeros(len(bank.ids))
	ww, jj, labels = d["weights"], d["jj"], d["labels"]
	for r in range(jj.shape[1]):
		reduced = ww.copy()
		reduced[:, r] = 0
		reduced /= reduced.sum(1, keepdims=True)
		removed = risk_from_weights(reduced, labels, bank.prior, 2)[0]
		delta = w * (loss(y, removed) - loss(y, p))
		np.add.at(summed, jj[:, r], delta)
		np.add.at(exposure, jj[:, r], w)
	value = summed / (exposure + 20)  # shrink rarely used donor estimates toward zero
	return value, pd.DataFrame(
		dict(
			reference_eid=bank.ids,
			tuning_exposure=exposure,
			shrunk_deletion_utility=value,
			total_loss_difference=summed,
		)
	)


class SupportRule:
	"""X-only support plus a tuning-trained expected-error diagnostic.

	Expected error is evaluated on a separate audit set. It is not a posterior
	uncertainty interval, guarantee for a person, or evidence of cohort invalidity.
	"""

	columns = [
		"nearest_distance",
		"reconstruction_RMSE",
		"reference_ESS",
		"normalized_entropy",
		"prior_weight",
		"reconstructed_assay_fraction",
	]

	def fit(self, info, risk, y, w, quantile=0.95, min_ess=3):
		from sklearn.ensemble import HistGradientBoostingRegressor

		self.quantile, self.min_ess = quantile, min_ess
		self.thresholds = {
			k: float(info[k].quantile(quantile))
			for k in ["nearest_distance", "reconstruction_RMSE"]
		}
		known = w > 0
		self.error_model = HistGradientBoostingRegressor(
			max_iter=80,
			max_leaf_nodes=5,
			min_samples_leaf=max(10, min(100, int(known.sum() / 15))),
			l2_regularization=10,
			early_stopping=False,
			random_state=583,
		)
		self.error_model.fit(
			info[self.columns].to_numpy()[known],
			(y[known] - risk[known]) ** 2,
			sample_weight=w[known],
		)
		self.error_threshold = float(np.quantile(self.error_score(info), quantile))
		return self

	def error_score(self, info):
		return np.maximum(self.error_model.predict(info[self.columns].to_numpy()), 0)

	def apply(self, info, missing, unknown, max_missing=0.2, use_error=True):
		reasons = [[] for _ in range(len(info))]
		tests = {
			key: info[key].to_numpy() > threshold
			for key, threshold in self.thresholds.items()
		}
		tests.update(
			low_reference_ESS=info.reference_ESS.to_numpy() < self.min_ess,
			no_independent_reference=info.reference_ESS.to_numpy() == 0,
			excess_missingness=np.asarray(missing) > max_missing,
			unknown_category=unknown,
		)
		if use_error:
			tests["high_expected_error"] = self.error_score(info) > self.error_threshold
		for key, bad in tests.items():
			for j in np.flatnonzero(bad):
				reasons[j].append(key)
		return np.array([not v for v in reasons]), np.array(
			[";".join(v) for v in reasons]
		)


# Models
# Classical comparators and strictly out-of-fold prototype quality.


def fit_logistic(x, y, w, c, ratio, seed, max_iter=2000):
	keep = w > 0
	if np.bincount(y[keep], minlength=2).min() < 5:
		raise ValueError("Too few observed cases/controls for a logistic fit")
	kw = dict(
		C=float(c),
		solver="saga" if ratio else "lbfgs",
		max_iter=max_iter,
		tol=1e-4,
		random_state=seed,
	)
	# sklearn 1.8+ encodes the penalty by l1_ratio; earlier supported versions
	# require an explicit penalty. Do not emit thousands of deprecation warnings.
	if (
		inspect.signature(LogisticRegression).parameters["penalty"].default
		== "deprecated"
	):
		kw["l1_ratio"] = ratio
	else:
		kw["penalty"] = "elasticnet" if ratio else "l2"
		if ratio:
			kw["l1_ratio"] = ratio
	with warnings.catch_warnings(record=True) as caught:
		warnings.simplefilter("always", ConvergenceWarning)
		model = LogisticRegression(**kw).fit(
			np.asarray(x[keep], dtype=np.float64),
			y[keep],
			sample_weight=w[keep] / np.mean(w[keep]),
		)
	if any(issubclass(v.category, ConvergenceWarning) for v in caught):
		raise ValueError(f"Logistic model did not converge in {max_iter} iterations")
	return model


def tune_logistic(x, y, w, build, tune, name, a, ratio):
	best = None
	rows = []
	for c in a.c_grid:
		try:
			model = fit_logistic(
				x[build], y[build], w[build], c, ratio, a.seed, a.max_iter
			)
			pred = model.predict_proba(x[tune])[:, 1]
			value = np.average(loss(y[tune], pred), weights=w[tune])
			rows.append(
				dict(model=name, C=c, validation_logloss=value, status="completed")
			)
			if best is None or value < best[0]:
				best = (value, model)
		except ValueError as exc:
			rows.append(
				dict(model=name, C=c, validation_logloss=np.nan, status=str(exc))
			)
	if best is None:
		raise ValueError(
			f"All candidates failed: {name}; see --max-iter and event counts"
		)
	return best[1], rows


class RiskCalibrator:
	"""Monotone logistic recalibration; optional additive clinical log-odds.

	Coefficients are nonnegative. Constant raw scores remain intercept-only.
	Fit only on the dedicated calibration subset; never on test outcomes.
	"""

	def fit(self, raw, y, w, clinical=None):
		x = logit(np.clip(raw, 1e-5, 1 - 1e-5))[:, None]
		if clinical is not None:
			x = np.c_[x, logit(np.clip(clinical, 1e-5, 1 - 1e-5))]
		self.mean, self.sd = x.mean(0), x.std(0)
		self.sd[self.sd < 1e-8] = 1
		z = np.c_[np.ones(len(x)), (x - self.mean) / self.sd]
		wn = w / w.sum()

		def objective(b):
			eta = z @ b
			val = np.sum(wn * (np.logaddexp(0, eta) - y * eta)) + 0.0001 * np.sum(
				b[1:] ** 2
			)
			grad = z.T @ (wn * (expit(eta) - y)) + np.r_[0, 0.0002 * b[1:]]
			return val, grad

		prior = np.clip(np.average(y, weights=w), 1e-5, 1 - 1e-5)
		fit = minimize(
			objective,
			np.r_[logit(prior), np.ones(x.shape[1])],
			jac=True,
			method="L-BFGS-B",
			bounds=[(None, None)] + [(0, None)] * x.shape[1],
		)
		if not fit.success:
			raise ValueError("Risk recalibration failed: " + fit.message)
		self.coef = fit.x
		return self

	def predict(self, raw, clinical=None):
		x = logit(np.clip(raw, 1e-5, 1 - 1e-5))[:, None]
		if clinical is not None:
			x = np.c_[x, logit(np.clip(clinical, 1e-5, 1 - 1e-5))]
		return expit(self.coef[0] + ((x - self.mean) / self.sd) @ self.coef[1:])


def oof_quality(raw, people, features, a, out):
	n, nf = raw.shape
	repeat_predictions = np.full((a.oof_repeats, n), np.nan)
	repeat_gains = np.full_like(repeat_predictions, np.nan)
	linear_predictions = np.full_like(repeat_predictions, np.nan)
	tree_predictions = np.full_like(repeat_predictions, np.nan)
	neural_predictions = np.full_like(repeat_predictions, np.nan)
	coefs, fit_rows = [], []
	ids = people[a.id_col].to_numpy(str)
	groups = people[a.group_col].to_numpy(str) if a.group_col else None
	# Unknown horizon outcomes form their own fold stratum.
	y_all = (people.event.eq(1) & people.time.le(a.horizon)).to_numpy(int)
	known = people.time.gt(a.horizon).to_numpy() | y_all.astype(bool)
	strata = np.where(known, y_all, 2)
	for repeat in range(a.oof_repeats):
		cv = (
			StratifiedGroupKFold(a.folds, shuffle=True, random_state=a.seed + repeat)
			if groups is not None
			else StratifiedKFold(a.folds, shuffle=True, random_state=a.seed + repeat)
		)
		for fold, (tr, va) in enumerate(cv.split(raw, strata, groups)):
			prep = MolecularPreprocessor(
				a.feature_missing,
				words(a.residualize),
				words(a.categorical),
				a.transform,
			).fit(raw[tr], people.iloc[tr])
			xt, mt = prep.transform(raw[tr], people.iloc[tr])
			xv, mv = prep.transform(raw[va], people.iloc[va])
			train_ok = (
				np.mean(~np.isfinite(raw[tr][:, prep.keep]), axis=1) <= a.sample_missing
			)
			km = CensoringKM().fit(
				people.iloc[tr[train_ok]], a.horizon, a.min_censor_survival
			)
			yt, wt = km.labels_weights(people.iloc[tr[train_ok]])
			yv, wv = km.labels_weights(people.iloc[va])
			model = fit_logistic(
				xt[train_ok],
				yt,
				wt,
				a.teacher_c,
				0.5,
				a.seed + repeat * a.folds + fold,
				a.max_iter,
			)
			linear_pred = model.predict_proba(xv)[:, 1]
			linear_predictions[repeat, va] = linear_pred
			pred = linear_pred
			if a.quality_teacher in ["ensemble", "all"]:
				tree = HistGradientBoostingClassifier(
					max_iter=a.quality_trees,
					max_leaf_nodes=7,
					learning_rate=0.05,
					min_samples_leaf=30,
					l2_regularization=10,
					early_stopping=False,
					random_state=a.seed + repeat * a.folds + fold,
				)
				known_train = wt > 0
				tree.fit(
					xt[train_ok][known_train],
					yt[known_train],
					sample_weight=wt[known_train] / wt[known_train].mean(),
				)
				tree_pred = tree.predict_proba(xv)[:, 1]
				tree_predictions[repeat, va] = tree_pred
				pred = 0.5 * (linear_pred + tree_pred)
			if a.quality_teacher in ["neural", "all"]:
				import copy

				child = copy.copy(a)
				child.seed = a.seed + repeat * a.folds + fold
				child.epochs, child.pretrain_epochs = (
					a.quality_neural_epochs,
					a.quality_pretrain_epochs,
				)
				child.patience = min(a.patience, 5)
				xx, mm = xt[train_ok], mt[train_ok]
				ii = ids[tr[train_ok]]
				gg = None if groups is None else groups[tr[train_ok]]
				membership, _ = token_partition(
					xx,
					[features[j] for j in prep.keep],
					a.tokens,
					child.seed,
					a.module_file,
				)
				net, _ = train_encoder(
					xx,
					mm,
					yt,
					wt,
					ii,
					gg,
					membership,
					child,
					abm_cache_dir(out, "quality_neural", f"repeat_{repeat + 1}_fold_{fold + 1}"),
				)
				zt = encode(net, xx, mm, device_for(a.device), a.batch_size)[0]
				zv = encode(net, xv, mv, device_for(a.device), a.batch_size)[0]
				net.cpu()
				bank = ReferenceBank(
					zt,
					xx,
					mm,
					ii,
					gg,
					yt,
					wt,
					np.flatnonzero(wt > 0),
					net,
					seed=child.seed,
				)
				neural_pred = bank.match(
					zv,
					ids[va],
					None if groups is None else groups[va],
					k=min(a.neural_k, len(bank.ids)),
					strength=a.train_prior,
				)[0]
				neural_predictions[repeat, va] = neural_pred
				pred = (
					neural_pred
					if a.quality_teacher == "neural"
					else (linear_pred + tree_pred + neural_pred) / 3
				)
			prior = np.average(yt, weights=wt)
			gain = loss(yv, prior) - loss(yv, pred)
			valid = (wv > 0) & (1 - mv.mean(1) <= a.sample_missing)
			repeat_predictions[repeat, va] = pred
			repeat_gains[repeat, va[valid]] = gain[valid]
			beta = np.zeros(nf)
			beta[prep.keep] = model.coef_[0]
			coefs.append(beta)
			fit_rows.append(
				dict(
					repeat=repeat + 1,
					fold=fold + 1,
					train_n=len(tr),
					validation_n=len(va),
					retained_features=len(prep.keep),
					iterations=int(model.n_iter_[0]),
				)
			)
		log("DONE", "OOF quality", f"repeat={repeat + 1}/{a.oof_repeats}")
	finite = np.all(np.isfinite(repeat_gains), axis=0)
	mean_gain = np.full(n, np.nan)
	sd_gain = np.full(n, np.nan)
	mean_gain[finite] = repeat_gains[:, finite].mean(0)
	sd_gain[finite] = repeat_gains[:, finite].std(0)
	fraction = (repeat_gains > 0).mean(0)
	reliable = finite & (mean_gain > a.min_gain) & (fraction >= a.fit_stability)
	table = pd.DataFrame(
		{
			a.id_col: ids,
			"horizon_label": y_all,
			"known_label": known,
			"OOF_probability": repeat_predictions.mean(0),
			"OOF_probability_SD": repeat_predictions.std(0),
			"OOF_elasticnet_probability": linear_predictions.mean(0),
			"OOF_tree_probability": tree_predictions.mean(0),
			"OOF_neural_probability": neural_predictions.mean(0),
			"OOF_logloss": loss(y_all, repeat_predictions.mean(0)),
			"fit_gain": mean_gain,
			"fit_gain_SD": sd_gain,
			"positive_gain_fraction": fraction,
			"reliable_candidate": reliable,
		}
	)
	table.to_csv(out / "reference_candidates.csv", index=False)
	pd.DataFrame(fit_rows).to_csv(out / "oof_fits.csv", index=False)
	beta = np.asarray(coefs)
	magnitude = np.mean(np.abs(beta), axis=0)
	stability = np.abs(np.mean(np.sign(beta), axis=0))
	importance = magnitude * stability
	pd.DataFrame(
		dict(
			feature=features,
			mean_abs_beta=magnitude,
			sign_stability=stability,
			metric_importance=importance,
		)
	).to_csv(out / "metric_features.csv", index=False)
	return table, importance


# Mosaic
# Module-specific literal outcome borrowing as an interpretable sensitivity.


class ModuleMosaic:
	def fit(
		self, x, observed, ids, groups, y, w, indices, membership, names, seed, k=10
	):
		self.banks, self.columns, self.names = [], [], names
		self.k = min(k, len(indices))
		for module in range(len(names)):
			ix = np.flatnonzero(membership == module)
			self.columns.append(ix)
			self.banks.append(
				ReferenceBank(
					x[:, ix] / np.sqrt(len(ix)),
					x[:, ix],
					observed[:, ix],
					ids,
					groups,
					y,
					w,
					indices,
					seed=seed,
				)
			)
		self.weights = np.full(len(names), 1 / len(names))
		return self

	def components(self, x, ids, groups):
		values, references = [], []
		for ix, bank in zip(self.columns, self.banks):
			risk, detail = bank.match(
				x[:, ix] / np.sqrt(len(ix)),
				ids,
				groups,
				k=self.k,
				temperature=1.0,
				strength=2.0,
				kernel="euclidean",
			)
			values.append(risk)
			references.append(bank.ids[detail["jj"][:, 0]])
		return np.array(values).T, np.array(references).T

	def tune(self, values, y, w):
		weights = w / w.sum()

		def objective(beta):
			p = np.clip(values @ beta, 1e-5, 1 - 1e-5)
			loss = np.sum(
				weights * (-y * np.log(p) - (1 - y) * np.log1p(-p))
			) + 0.01 * np.sum(beta**2)
			grad = values.T @ (weights * (p - y) / (p * (1 - p))) + 0.02 * beta
			return loss, grad

		fit = minimize(
			objective,
			self.weights,
			jac=True,
			method="SLSQP",
			bounds=[(0, 1)] * len(self.weights),
			constraints=[
				{
					"type": "eq",
					"fun": lambda b: b.sum() - 1,
					"jac": lambda b: np.ones_like(b),
				}
			],
			options={"maxiter": 1000, "ftol": 1e-10},
		)
		if not fit.success:
			raise ValueError("Module mixture failed: " + fit.message)
		self.weights = fit.x / fit.x.sum()
		return self


# Io Data
# Read LE8-layout inputs without executing the LE8 analysis environment.


def validate_ids(frame, id_col):
	if id_col not in frame:
		raise ValueError(f"Missing ID column: {id_col}")
	if frame.columns.duplicated().any():
		raise ValueError("Duplicate column names")
	if frame[id_col].isna().any():
		raise ValueError("Missing person IDs")
	raw = frame[id_col]
	if pd.api.types.is_numeric_dtype(raw):
		num = raw.to_numpy(float)
		if not np.isfinite(num).all() or (num != np.floor(num)).any():
			raise ValueError(
				"Numeric IDs must be finite integers; supply strings otherwise"
			)
		frame[id_col] = raw.astype("int64").astype(str)
	else:
		frame[id_col] = raw.astype(str).str.strip()
	if frame[id_col].eq("").any() or frame[id_col].duplicated().any():
		raise ValueError(
			"Empty or duplicate person IDs; select one baseline visit first"
		)
	return frame


def read_table(path, id_col="eid", columns=None, r_bin="Rscript"):
	path = Path(path)
	if path.suffix.lower() == ".rds":
		if shutil.which(r_bin):
			# R selects columns before serialization, avoiding a second full RDS in Python.
			with tempfile.TemporaryDirectory(prefix="abm_rds_") as tmp:
				out = Path(tmp) / "selected.rds"
				cols = Path(tmp) / "columns.txt"
				cols.write_text("\n".join(columns or []))
				bridge = Path(tmp) / "export_rds.R"
				bridge.write_text(RDS_EXPORT_SCRIPT, encoding="utf-8")
				subprocess.run(
					[r_bin, str(bridge), str(path), str(out), str(cols)], check=True
				)
				import pyreadr

				frame = pyreadr.read_r(str(out))[None].reset_index(drop=True)
		else:
			import pyreadr

			objects = pyreadr.read_r(str(path))
			if len(objects) != 1 or not isinstance(
				next(iter(objects.values())), pd.DataFrame
			):
				raise ValueError(f"{path}: expected a single R data.frame")
			frame = next(iter(objects.values()))
	elif path.suffix.lower() == ".parquet":
		frame = pd.read_parquet(path, columns=columns)
	else:
		sep = "," if ".csv" in path.name.lower() else "\t"
		# Reject duplicate raw headers before pandas can silently suffix them.
		header = pd.read_csv(path, sep=sep, header=None, nrows=1).iloc[0].astype(str)
		if header.duplicated().any():
			raise ValueError(f"Duplicate header in {path}")
		frame = pd.read_csv(
			path, sep=sep, dtype={id_col: str}, usecols=columns, low_memory=False
		)
	if columns:
		missing = set(columns) - set(frame)
		if missing:
			raise ValueError(f"{path}: missing columns {sorted(missing)}")
		frame = frame[columns].copy()
	return validate_ids(frame, id_col)


def numeric(frame, columns, context):
	values = []
	for name in columns:
		raw = frame[name]
		val = pd.to_numeric(raw, errors="coerce")
		bad = raw.notna() & val.isna()
		if bad.any():
			raise ValueError(f"{context}: nonnumeric values in {name}")
		values.append(val.to_numpy(dtype="float32", na_value=np.nan))
	if not values:
		raise ValueError(f"{context}: no molecular columns")
	result = np.column_stack(values)
	result[~np.isfinite(result)] = np.nan
	return result


def safe_expression(expression, frame):
	"""Only arithmetic over explicit columns; never eval spreadsheet expressions."""

	def calculate(node):
		if isinstance(node, ast.Name):
			if node.id not in frame:
				raise KeyError(node.id)
			return pd.to_numeric(frame[node.id], errors="raise").to_numpy(float)
		if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
			return float(node.value)
		if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
			value = calculate(node.operand)
			return value if isinstance(node.op, ast.UAdd) else -value
		if isinstance(node, ast.BinOp):
			left, right = calculate(node.left), calculate(node.right)
			with np.errstate(divide="ignore", invalid="ignore"):
				if isinstance(node.op, ast.Div):
					return np.divide(left, right)
				if isinstance(node.op, ast.Mult):
					return left * right
				if isinstance(node.op, ast.Add):
					return left + right
				if isinstance(node.op, ast.Sub):
					return left - right
		raise ValueError(f"Unsupported metabolite expression: {expression}")

	return calculate(ast.parse(expression, mode="eval").body)


def map_metabolites(frame, mapping, id_col):
	# Match ukb/f/phe.R: baseline _i0 fields and common/met.lst.
	repeated = [c for c in frame if re.search(r"_i[1-9]\d*$", str(c))]
	frame = frame.drop(columns=repeated).copy()
	frame.columns = [re.sub(r"_i0$", "", str(c)) for c in frame]
	if frame.columns.duplicated().any():
		raise ValueError("Metabolite columns collide after removing _i0")
	mapping = Path(mapping)
	if not mapping.is_file():
		raise FileNotFoundError(f"Raw metabolite mapping required: {mapping}")
	spec = pd.read_csv(
		mapping,
		sep="\t",
		header=None,
		usecols=[0, 1],
		names=["expression", "feature"],
		dtype=str,
	)
	# Support both the original headerless mapping and the current UKB mapping.
	if len(spec) and tuple(spec.iloc[0]) == ("data_field", "met_name"):
		spec = spec.iloc[1:].copy()
	if spec.isna().any().any() or spec.feature.duplicated().any():
		raise ValueError("met.lst must contain unique, nonempty feature names")
	result = pd.DataFrame({id_col: frame[id_col]})
	audit = []
	for row in spec.itertuples():
		expr, feature = row.expression.strip(), row.feature.strip()
		try:
			values = (
				pd.to_numeric(frame[expr], errors="raise").to_numpy(float)
				if expr in frame
				else safe_expression(expr, frame)
			)
			values = np.broadcast_to(values, (len(frame),)).copy()
			invalid = ~np.isfinite(values)
			values[invalid] = np.nan
			result[feature] = values.astype("float32")
			audit.append(
				dict(
					feature=feature,
					expression=expr,
					status="mapped",
					nonfinite_to_missing=int(invalid.sum()),
				)
			)
		except KeyError as exc:
			audit.append(
				dict(
					feature=feature,
					expression=expr,
					status="missing_source",
					reason=str(exc),
				)
			)
	if result.shape[1] < 4:
		raise ValueError(
			"Fewer than three metabolites map; check raw baseline fields and met.lst"
		)
	return result, pd.DataFrame(audit)


def phenotype_columns(a):
	cols = [a.id_col] + words(a.covariates) + words(a.residualize)
	if a.outcome_type == "survival":
		cols += [a.baseline_col, a.diagnosis_col, a.death_col, a.lost_col]
		cols += words(a.healthy_date_cols)
		if a.disease_evidence_col:
			cols += [a.disease_evidence_col]
	else:
		cols += [a.target_col]
	if a.group_col:
		cols += [a.group_col]
	return list(dict.fromkeys(cols))


def generate_demo(a):
	rng = np.random.default_rng(a.seed)
	n = a.max_samples or 1500
	nf = a.demo_features
	state = rng.integers(0, 3, n)
	latent = rng.normal(size=(n, 6))
	latent[:, 0] += (state - 1) * 4
	x = (latent @ rng.normal(size=(6, nf)) + rng.normal(size=(n, nf))).astype("float32")
	age = rng.uniform(40, 75, n)
	sex = rng.integers(0, 2, n)
	plate = rng.integers(1, 6, n)
	x += plate[:, None] * 0.06
	x[rng.random(x.shape) < 0.03] = np.nan
	risk = 0.03 * (age - 55) + 0.35 * latent[:, 1] + 0.4 * (state - 1) * latent[:, 2]
	times = rng.exponential(18 * np.exp(-risk))
	censor = rng.uniform(7, 12, n)
	base = pd.Timestamp("2010-01-01")
	p = pd.DataFrame(
		{
			a.id_col: [f"demo_{i}" for i in range(n)],
			"age": age,
			"sex": sex,
			"tdi": rng.normal(size=n),
			"PC1": rng.normal(size=n),
			"PC2": rng.normal(size=n),
			"center": rng.choice(["A", "B", "C"], n),
			f"{a.biom}.plate": plate,
			a.baseline_col: "2010-01-01",
			a.death_col: pd.Series([None] * n, dtype=object),
			a.lost_col: (base + pd.to_timedelta(censor * 365.25, unit="D")).strftime(
				"%Y-%m-%d"
			),
			a.diagnosis_col: (
				base + pd.to_timedelta(np.minimum(times, 100) * 365.25, unit="D")
			).strftime("%Y-%m-%d"),
		}
	)
	if a.outcome_type == "quantitative":
		p[a.target_col] = (
			165
			+ 7 * sex
			+ 2.5 * latent[:, 1]
			+ 2 * (state - 1) * latent[:, 2]
			+ rng.normal(size=n)
		)
	if a.group_col:
		p[a.group_col] = np.arange(n) // 2
	features = [f"F{j:04d}" for j in range(nf)]
	omics = pd.DataFrame(x, columns=features)
	omics.insert(0, a.id_col, p[a.id_col])
	return p, omics


def prepare(a, out):
	mapping_audit = None
	if a.demo:
		p, omics = generate_demo(a)
	else:
		p = read_table(a.phe_file, a.id_col, phenotype_columns(a), a.r_bin)
		omics = read_table(a.omics_file, a.id_col, r_bin=a.r_bin)
		if a.biom == "met" and a.met_input == "raw":
			omics, mapping_audit = map_metabolites(omics, a.met_map, a.id_col)
		if a.biom == "prot":
			omics.columns = [c if c == a.id_col else str(c).upper() for c in omics]
			if omics.columns.duplicated().any():
				raise ValueError("Protein names collide after uppercasing")
	missing = set(phenotype_columns(a)) - set(p)
	if missing:
		raise ValueError(f"Missing phenotype columns: {sorted(missing)}")
	omics = omics[omics[a.id_col].isin(p[a.id_col])].copy()
	if a.max_samples and len(omics) > a.max_samples:
		omics = omics.sample(a.max_samples, random_state=a.seed)
	if omics.empty:
		raise ValueError("No matched phenotype/omics IDs")
	p = p.set_index(a.id_col).loc[omics[a.id_col]].reset_index()
	features = [
		c for c in omics if c != a.id_col and c not in words(a.exclude_features)
	]
	forbidden = set(phenotype_columns(a)) | {"event", "time", "target", "split", "Y"}
	collision = {c.casefold() for c in features} & {c.casefold() for c in forbidden}
	if collision:
		raise ValueError(
			"Metadata/outcome columns in omics matrix: " + ", ".join(sorted(collision))
		)
	x = numeric(omics, features, "omics")
	# Retain original missingness. Scale transforms are part of train-fitted preprocessing.
	p.to_csv(out / "phenotype.csv", index=False)
	np.save(out / "raw.npy", x)
	(out / "features.txt").write_text("\n".join(features) + "\n")
	if mapping_audit is not None:
		mapping_audit.to_csv(out / "metabolite_mapping.csv", index=False)
	upstream = "supplied matrix; upstream preprocessing not inferred"
	if str(a.omics_file).endswith(f"Rdata/{a.biom}.rds"):
		upstream = "LE8 cleaned RDS: upstream all-person QC/imputation may already have occurred"
	dump(
		out / "input_audit.json",
		dict(
			n=len(p),
			features=len(features),
			source="SYNTHETIC" if a.demo else upstream,
			met_input=a.met_input,
			missing_fraction=float(np.mean(~np.isfinite(x))),
		),
	)
	return p


# Evidence
# Individual evidence and falsifiable explanation checks; no causal claims.


def write_match_table(path, idcol, ids, bank, detail, extra=None):
	jj = detail["jj"]
	data = {
		idcol: np.repeat(ids, jj.shape[1]),
		"rank": np.tile(np.arange(jj.shape[1]) + 1, len(ids)),
		"reference_eid": bank.ids[jj.ravel()],
		"reference_observed_Y": detail["labels"].ravel(),
		"normalized_copy_weight": detail["weights"].ravel(),
		"risk_coefficient": detail["coefficient"].ravel(),
		"raw_risk_contribution": (detail["coefficient"] * detail["labels"]).ravel(),
		"reference_IPCW": bank.w[jj.ravel()],
	}
	if "head_weights" in detail:
		for h in range(detail["head_gate"].shape[1]):
			head = detail["head_weights"][:, :, h]
			gate = np.repeat(detail["head_gate"][:, h, None], jj.shape[1], axis=1)
			coefficient = (1 - detail["prior_weight"][:, None]) * gate * head
			data[f"head_{h + 1}_within_head_weight"] = head.ravel()
			data[f"head_{h + 1}_gate"] = gate.ravel()
			data[f"head_{h + 1}_raw_contribution"] = (
				coefficient * detail["labels"]
			).ravel()
	if extra:
		data.update({k: v.ravel() for k, v in extra.items()})
	table = pd.DataFrame(data)
	table.loc[table.normalized_copy_weight > 0].to_csv(
		path, index=False, compression="infer"
	)


def write_evidence(
	out,
	a,
	ids,
	x,
	observed,
	bank,
	d,
	calibrated,
	features,
	membership,
	token_names,
	calibrator,
	strength,
):
	raw = (
		np.sum(d["coefficient"] * d["labels"], axis=1) + d["prior_weight"] * bank.prior
	)
	ww, labels = d["weights"], d["labels"]
	deletion = np.zeros_like(ww)
	for r in range(ww.shape[1]):
		reduced = ww.copy()
		reduced[:, r] = 0
		sums = reduced.sum(1, keepdims=True)
		np.divide(reduced, sums, out=reduced, where=sums > 0)
		risk = risk_from_weights(reduced, labels, bank.prior, strength)[0]
		risk = np.where(sums[:, 0] > 0, risk, bank.prior)
		deletion[:, r] = calibrated - calibrator.predict(risk)
	flip_raw = d["coefficient"] * (1 - 2 * labels)
	write_match_table(
		out / "test_reference_matches.csv.gz",
		a.id_col,
		ids,
		bank,
		d,
		{
			"calibrated_delta_delete_fixed_neighbor": deletion,
			"raw_delta_flip_donor_Y": flip_raw,
		},
	)
	pd.DataFrame(
		{
			a.id_col: ids,
			"reference_contribution_sum": np.sum(d["coefficient"] * labels, axis=1),
			"prior_weight": d["prior_weight"],
			"build_prior": bank.prior,
			"prior_contribution": d["prior_weight"] * bank.prior,
			"reconstructed_raw_risk": raw,
			"reconstructed_calibrated_risk": calibrator.predict(raw),
			"calibrated_risk": calibrated,
			"max_absolute_reference_deletion_delta": np.abs(deletion).max(1),
		}
	).to_csv(out / "individual_risk_decomposition.csv", index=False)
	if not np.allclose(calibrator.predict(raw), calibrated, atol=1e-8):
		raise AssertionError(
			"Reference contributions fail to reconstruct deployed probability"
		)
	select = sorted(range(len(ids)), key=lambda i: stable_seed(ids[i], a.seed))[
		: a.explanation_samples
	]
	reconstruction, _ = (
		bank.reconstruct({k: v[select] for k, v in d.items()})
		if select
		else (None, None)
	)
	with gzip.open(
		out / "individual_explanations.jsonl.gz", "wt", encoding="utf-8"
	) as handle:
		for pos, i in enumerate(select):
			overlap = observed[i]
			mismatch = np.where(overlap, np.abs(x[i] - reconstruction[pos]), -np.inf)
			common = np.where(
				overlap & (x[i] * reconstruction[pos] > 0),
				np.minimum(np.abs(x[i]), np.abs(reconstruction[pos])),
				-np.inf,
			)

			def describe(score):
				ix = np.argsort(-score, kind="stable")[:8]
				return [
					dict(
						feature=features[j],
						token=token_names[membership[j]],
						query_standardized=float(x[i, j]),
						reference_reconstruction=float(reconstruction[pos, j]),
					)
					for j in ix
					if np.isfinite(score[j])
				]

			row = dict(
				eid=str(ids[i]),
				raw_borrowed_risk=float(raw[i]),
				calibrated_net_risk=float(calibrated[i]),
				prior_contribution=float(d["prior_weight"][i] * bank.prior),
				reference_ESS=float(d["ESS"][i]),
				shared_abnormal_assays=describe(common),
				largest_mismatches=describe(mismatch),
				interpretation="descriptive matching evidence; attention/contributions are not causal protein effects",
			)
			handle.write(json.dumps(row, ensure_ascii=False) + "\n")
	count = pd.Series(bank.ids[d["jj"][d["weights"] > 0]]).value_counts()
	utilization = pd.DataFrame({"reference_eid": bank.ids, "observed_Y": bank.y})
	utilization["test_matches"] = (
		utilization.reference_eid.map(count).fillna(0).astype(int)
	)
	utilization.to_csv(out / "reference_utilization.csv", index=False)


def main_prediction(bundle, x, observed, ids, groups, random=False):
	name = bundle["primary"]
	spec = bundle["specs"][name]
	model = bundle["encoders"][spec["encoder"]]
	cfg = bundle["config"]
	z, _, reconstructed = encode(
		model, x, observed, device_for(cfg["device"]), cfg["batch_size"], True
	)
	model.cpu()
	kw = {k: v for k, v in spec.items() if k not in ["encoder", "space"]}
	p, detail = bundle["banks"][name].match(
		z, ids, groups, **kw, random_candidates=random
	)
	return p, detail, reconstructed


def subset_indices(ids, count, seed):
	return np.array(
		sorted(range(len(ids)), key=lambda i: stable_seed(ids[i], seed))[:count], int
	)


def jaccard(a, b):
	return np.array([len(set(x) & set(y)) / len(set(x) | set(y)) for x, y in zip(a, b)])


def masked_validation(bundle, x, observed, ids, groups, out, a):
	selected = subset_indices(ids, a.explanation_samples, a.seed)
	x, observed, ids = x[selected], observed[selected], ids[selected]
	groups = None if groups is None else groups[selected]
	original, orig_detail, _ = main_prediction(bundle, x, observed, ids, groups)
	bank = bundle["banks"][bundle["primary"]]
	rows, feature_stats = [], {}
	for repeat in range(a.mask_repeats):
		hidden = np.zeros_like(observed)
		for i, eid in enumerate(ids):
			valid = np.flatnonzero(observed[i])
			n = max(1, int(round(len(valid) * a.mask_fraction)))
			if len(valid):
				hidden[
					i,
					np.random.default_rng(
						stable_seed(eid, a.seed + repeat + 981)
					).choice(valid, min(n, len(valid)), replace=False),
				] = True
		query = x.copy()
		query[hidden] = 0
		mask = observed & ~hidden
		pred, detail, decoder = main_prediction(bundle, query, mask, ids, groups)
		_, random_detail, _ = main_prediction(
			bundle, query, mask, ids, groups, random=True
		)
		copied, cover = bank.reconstruct(detail)
		random_copy, _ = bank.reconstruct(random_detail)
		pca = bundle["geometry"].unsupervised
		pca_copy = pca.inverse_transform(pca.transform(query))
		methods = {
			"reference_attention": copied,
			"decoder": decoder,
			"random_references": random_copy,
			"pca_reconstruction": pca_copy,
			"zero_build_mean": np.zeros_like(x),
		}
		nj = jaccard(orig_detail["jj"], detail["jj"])
		for name, estimate in methods.items():
			square = np.where(hidden, (estimate - x) ** 2, 0)
			n = hidden.sum(1)
			for i, eid in enumerate(ids):
				rows.append(
					dict(
						eid=eid,
						repeat=repeat + 1,
						method=name,
						masked_count=int(n[i]),
						masked_RMSE=float(np.sqrt(square[i].sum() / max(n[i], 1))),
						match_Jaccard=(
							float(nj[i]) if name == "reference_attention" else np.nan
						),
						absolute_raw_risk_change=(
							float(abs(pred[i] - original[i]))
							if name == "reference_attention"
							else np.nan
						),
						copied_assay_coverage=(
							float(cover[i, hidden[i]].mean())
							if name == "reference_attention"
							else np.nan
						),
					)
				)
			stats = feature_stats.setdefault(name, np.zeros((4, x.shape[1])))
			stats[0] += hidden.sum(0)
			stats[1] += np.where(hidden, x, 0).sum(0)
			stats[2] += np.where(hidden, x * x, 0).sum(0)
			stats[3] += square.sum(0)
		log(
			"DONE",
			"masked_validation",
			f"repeat={repeat + 1}/{a.mask_repeats}, n={len(ids)}",
		)
	person_metrics = pd.DataFrame(rows)
	person_metrics.to_csv(out / "masked_reconstruction.csv", index=False)
	person_metrics.groupby("method", as_index=False)["masked_RMSE"].mean().to_csv(
		out / "masked_reconstruction_summary.csv", index=False
	)
	table = []
	for name, stats in feature_stats.items():
		n, total, total_sq, sse = stats
		ss = total_sq - total**2 / np.maximum(n, 1)
		r2 = np.full(len(n), np.nan)
		ok = (n >= 3) & (ss > 1e-8)
		r2[ok] = 1 - sse[ok] / ss[ok]
		for j, feature in enumerate(bundle["retained_features"]):
			table.append(
				dict(
					method=name,
					feature=feature,
					masked_count=int(n[j]),
					masked_R2=r2[j],
					masked_RMSE=np.sqrt(sse[j] / n[j]) if n[j] else np.nan,
				)
			)
	pd.DataFrame(table).to_csv(out / "masked_feature_metrics.csv", index=False)


def module_perturbations(bundle, x, observed, ids, groups, out, a):
	ix = subset_indices(
		ids, min(a.explanation_samples, a.perturbation_samples), a.seed + 3
	)
	if not len(ix):
		return
	x, observed, ids = x[ix], observed[ix], ids[ix]
	groups = None if groups is None else groups[ix]
	baseline, detail, _ = main_prediction(bundle, x, observed, ids, groups)
	cal = bundle["calibrators"][bundle["primary"]]
	original = cal.predict(baseline)
	rows = []
	for module, name in enumerate(bundle["token_names"]):
		query, mask = x.copy(), observed.copy()
		mask[:, bundle["membership"] == module] = False
		query[:, bundle["membership"] == module] = 0
		changed, other, _ = main_prediction(bundle, query, mask, ids, groups)
		calibrated_changed = cal.predict(changed)
		overlap = jaccard(detail["jj"], other["jj"])
		for i, eid in enumerate(ids):
			rows.append(
				dict(
					eid=eid,
					token=name,
					original_risk=original[i],
					risk_after_module_mask=calibrated_changed[i],
					delta_original_minus_masked=original[i] - calibrated_changed[i],
					match_Jaccard=overlap[i],
					interpretation="missing-module sensitivity; not a treatment effect",
				)
			)
		log("DONE", "module_perturbation", name)
	pd.DataFrame(rows).to_csv(
		out / "module_perturbations.csv.gz", index=False, compression="gzip"
	)


def same_risk_pairs(ids, risk, profile, names, caliper=0.01, max_pairs=100):
	order = np.argsort(risk, kind="stable")
	used, rows = set(), []
	# Deterministic, no outcome inspection; cap candidates to avoid an N^2 matrix.
	for pos in np.linspace(
		0, max(0, len(order) - 1), min(len(order), max_pairs * 20), dtype=int
	):
		i = order[pos]
		if i in used:
			continue
		lo = np.searchsorted(risk[order], risk[i] - caliper)
		hi = np.searchsorted(risk[order], risk[i] + caliper, side="right")
		candidate = order[lo:hi]
		candidate = np.array([j for j in candidate if j != i and j not in used], int)
		if not len(candidate):
			continue
		if len(candidate) > 500:
			candidate = candidate[np.linspace(0, len(candidate) - 1, 500, dtype=int)]
		distance = np.linalg.norm(profile[candidate] - profile[i], axis=1)
		j = candidate[np.argmax(distance)]
		used.update([i, j])
		contrast = np.argsort(-np.abs(profile[i] - profile[j]))[:3]
		rows.append(
			dict(
				eid_A=ids[i],
				eid_B=ids[j],
				predicted_risk_A=risk[i],
				predicted_risk_B=risk[j],
				profile_distance=float(np.max(distance)),
				largest_contrasts=";".join(names[k] for k in contrast),
			)
		)
		if len(rows) >= max_pairs:
			break
	return pd.DataFrame(
		rows,
		columns=[
			"eid_A",
			"eid_B",
			"predicted_risk_A",
			"predicted_risk_B",
			"profile_distance",
			"largest_contrasts",
		],
	)


# Attention Audit
# Export actual attention tensors and perturbation evidence, without using test Y.


def rollout(maps):
	"""Mean-head + residual rollout; omits FFN/LN/value mixing, NOT attribution."""
	b, layers, heads, tokens, _ = maps.shape
	joint = np.broadcast_to(np.eye(tokens), (b, tokens, tokens)).copy()
	for level in range(layers):
		flow = maps[:, level].mean(1) + np.eye(tokens)[None]
		flow /= flow.sum(-1, keepdims=True)
		joint = flow @ joint
	importance = joint[:, 0, 1:]
	return importance / np.maximum(importance.sum(1, keepdims=True), 1e-30)


def audit_attention(bundle, x, observed, ids, groups, out, a):
	indices = subset_indices(ids, a.attention_samples, a.seed + 3)
	if not len(indices):
		return
	x, observed, ids = x[indices], observed[indices], ids[indices]
	groups = None if groups is None else groups[indices]
	out = Path(out) / "attention"
	out.mkdir(exist_ok=True)
	name = bundle["primary"]
	spec = bundle["specs"][name]
	kw = {k: v for k, v in spec.items() if k not in ["encoder", "space"]}
	model = bundle["encoders"][spec["encoder"]]
	bank = bundle["banks"][name]
	cal = bundle["calibrators"][name]
	dev = device_for(a.device)
	model.to(dev).eval()
	encoded, maps = [], []
	with torch.no_grad():
		for begin in range(0, len(ids), a.batch_size):
			z, _, _, att = model(
				torch.as_tensor(x[begin : begin + a.batch_size], device=dev),
				torch.as_tensor(observed[begin : begin + a.batch_size], device=dev),
				return_attention=True,
			)
			encoded.append(z.cpu().numpy())
			maps.append(att.cpu().numpy())
	model.cpu()
	z, maps = np.concatenate(encoded), np.concatenate(maps)
	if not np.allclose(maps.sum(-1), 1, atol=2e-6):
		raise AssertionError("Attention rows do not sum to one")
	# Capturing the maps must preserve the deployed model output.
	normal_z = encode(model, x, observed, dev, a.batch_size)[0]
	model.cpu()
	np.testing.assert_allclose(z, normal_z, rtol=2e-5, atol=2e-6)
	# Use normal execution for every risk; captured maps are evidence only.
	raw, detail = bank.match(normal_z, ids, groups, **kw)
	baseline = cal.predict(raw)
	tokens = ["CLS"] + list(bundle["token_names"])
	importance = rollout(maps)
	write_array_rds(
		out / "self_attention.rds",
		dict(
			eid=ids.astype(str),
			token=np.array(tokens),
			attention=maps,
			rollout=importance,
			axis_order=np.array(["person", "layer", "head", "query_token", "key_token"]),
		),
		a.r_bin,
	)
	table = []
	for i, eid in enumerate(ids):
		for level in range(maps.shape[1]):
			for head in range(maps.shape[2]):
				matrix = maps[i, level, head]
				for query in range(len(tokens)):
					for key in np.argsort(-matrix[query], kind="stable")[:3]:
						table.append(
							(
								eid,
								level + 1,
								head + 1,
								tokens[query],
								tokens[key],
								matrix[query, key],
							)
						)
	pd.DataFrame(
		table,
		columns=[
			"eid",
			"layer",
			"head",
			"query_token",
			"key_token",
			"attention_weight",
		],
	).to_csv(out / "self_attention_top_edges.csv.gz", index=False)
	summary = []
	for i, eid in enumerate(ids):
		for level in range(maps.shape[1]):
			for head in range(maps.shape[2]):
				matrix = maps[i, level, head]
				summary.append(
					dict(
						eid=eid,
						layer=level + 1,
						head=head + 1,
						mean_query_entropy=float(
							-np.sum(
								matrix * np.log(np.maximum(matrix, 1e-30)), -1
							).mean()
						),
						mean_self_weight=float(np.trace(matrix) / len(matrix)),
						across_query_weight_SD=float(matrix.std(0).mean()),
					)
				)
	pd.DataFrame(summary).to_csv(out / "self_attention_diagnostics.csv", index=False)
	pd.DataFrame(importance, index=ids, columns=tokens[1:]).rename_axis("eid").to_csv(
		out / "rollout_heuristic.csv"
	)
	experiments = []

	def record(label, changed, changed_detail):
		for i, eid in enumerate(ids):
			experiments.append(
				dict(
					eid=eid,
					intervention=label,
					original_raw=raw[i],
					changed_raw=changed[i],
					original_calibrated=baseline[i],
					changed_calibrated=cal.predict(changed)[i],
					absolute_calibrated_change=abs(
						baseline[i] - cal.predict(changed)[i]
					),
					selected_reference_Jaccard=len(
						set(detail["jj"][i]) & set(changed_detail["jj"][i])
					)
					/ len(set(detail["jj"][i]) | set(changed_detail["jj"][i])),
				)
			)

	# Real head weights reconstruct the marginal donor weight, then raw risk.
	if "head_weights" in detail:
		heads, gate = detail["head_weights"], detail["head_gate"]
		np.testing.assert_allclose(
			(heads * gate[:, None, :]).sum(-1), detail["weights"], atol=2e-7
		)
		rows = []
		for i, eid in enumerate(ids):
			for rank, j in enumerate(detail["jj"][i]):
				for h in range(gate.shape[1]):
					coefficient = (
						(1 - detail["prior_weight"][i]) * gate[i, h] * heads[i, rank, h]
					)
					rows.append(
						(
							eid,
							h + 1,
							rank + 1,
							bank.ids[j],
							detail["labels"][i, rank],
							gate[i, h],
							heads[i, rank, h],
							coefficient,
							coefficient * detail["labels"][i, rank],
						)
					)
		pd.DataFrame(
			rows,
			columns=[
				"eid",
				"head",
				"reference_rank",
				"reference_eid",
				"observed_Y",
				"head_gate",
				"within_head_copy_weight",
				"risk_coefficient",
				"raw_risk_contribution",
			],
		).to_csv(out / "reference_head_contributions.csv.gz", index=False)
		head_table = []
		for i, eid in enumerate(ids):
			for h in range(gate.shape[1]):
				head_table.append(
					dict(
						eid=eid,
						head=h + 1,
						gate=gate[i, h],
						borrowed_head_risk=detail["head_risk"][i, h],
						effective_donors=1 / max(np.sum(heads[i, :, h] ** 2), 1e-30),
						prior_weight=detail["prior_weight"][i],
					)
				)
		pd.DataFrame(head_table).to_csv(out / "reference_heads.csv", index=False)
		if gate.shape[1] > 1:
			for h in range(gate.shape[1]):
				reduced = gate.copy()
				reduced[:, h] = 0
				reduced /= reduced.sum(1, keepdims=True)
				weights = (heads * reduced[:, None, :]).sum(-1)
				changed = risk_from_weights(
					weights, detail["labels"], bank.prior, kw["strength"]
				)[0]
				record(f"drop_cross_head_{h + 1}_fixed_neighbors", changed, detail)
	# Uniform/identity attention is an intervention, not an independently retrained comparator.
	old = [layer.intervention for layer in model.context.layers]
	try:
		for mode in ["uniform", "identity"]:
			for layer in model.context.layers:
				layer.intervention = mode
			changed_z = encode(model, x, observed, dev, a.batch_size)[0]
			model.cpu()
			changed, other = bank.match(changed_z, ids, groups, **kw)
			record(mode + "_query_only_frozen_bank", changed, other)
			changed_bank = copy.copy(bank)
			changed_bank.z = encode(model, bank.x, bank.observed, dev, a.batch_size)[0]
			model.cpu()
			changed, other = changed_bank.match(changed_z, ids, groups, **kw)
			record(mode + "_query_and_bank", changed, other)
	finally:
		for layer, value in zip(model.context.layers, old):
			layer.intervention = value
		model.cpu()
	# Does high attention predict sensitivity better than masking random modules?
	n_mask = min(2, len(tokens) - 1)
	sets = {
		"top_rollout": np.argsort(-importance, axis=1)[:, :n_mask],
		"bottom_rollout": np.argsort(importance, axis=1)[:, :n_mask],
		"random_modules": np.array(
			[
				np.random.default_rng(stable_seed(eid, a.seed + 931)).choice(
					len(tokens) - 1, n_mask, replace=False
				)
				for eid in ids
			]
		),
	}
	for label, selected in sets.items():
		mask = observed.copy()
		for i in range(len(ids)):
			mask[i, np.isin(bundle["membership"], selected[i])] = False
		zz = encode(model, x, mask, dev, a.batch_size)[0]
		model.cpu()
		changed, other = bank.match(zz, ids, groups, **kw)
		record("mask_" + label, changed, other)
	frame = pd.DataFrame(experiments)
	frame.to_csv(out / "attention_interventions.csv.gz", index=False)
	frame.groupby("intervention").agg(
		n=("eid", "size"),
		mean_absolute_change=("absolute_calibrated_change", "mean"),
		median_absolute_change=("absolute_calibrated_change", "median"),
		mean_match_Jaccard=("selected_reference_Jaccard", "mean"),
	).to_csv(out / "intervention_summary.csv")
	# All-module perturbations from the existing module audit provide a stronger
	# person-wise correspondence check; no manufactured correlation if constant.
	module_file = Path(out).parent / "module_perturbations.csv.gz"
	correspondence = []
	if module_file.exists():
		module = pd.read_csv(module_file, dtype={"eid": str})
		lookup = {str(eid): i for i, eid in enumerate(ids)}
		for eid, sub in module.groupby("eid"):
			if eid not in lookup:
				continue
			ordered = sub.set_index("token").reindex(tokens[1:])
			sensitivity = ordered.delta_original_minus_masked.abs().to_numpy()
			score = importance[lookup[eid]]
			ok = np.isfinite(sensitivity)
			corr = (
				spearmanr(score[ok], sensitivity[ok]).statistic
				if ok.sum() > 2
				and np.std(sensitivity[ok]) > 1e-12
				and np.std(score[ok]) > 1e-12
				else np.nan
			)
			correspondence.append(
				dict(eid=eid, rollout_vs_masking_spearman=corr, n_tokens=int(ok.sum()))
			)
	pd.DataFrame(
		correspondence, columns=["eid", "rollout_vs_masking_spearman", "n_tokens"]
	).to_csv(out / "attention_vs_sensitivity.csv", index=False)
	dump(
		out / "audit.json",
		dict(
			n_people=len(ids),
			selection="outcome-blind stable ID hash",
			self_attention_shape=maps.shape,
			axis_order=["person", "layer", "head", "query_token", "key_token"],
			maps="actual pre-dropout Q/K softmax in eval mode; exported for every layer/head",
			cross_values="observed build donor outcomes, IPCW-adjusted within each head",
			caveats=[
				"Attention and rollout are not causal effects or complete feature attributions.",
				"Interventions keep fitted parameters/calibration frozen; may move inputs out of distribution.",
				"Module count is matched in masking controls, not number of assays; modules differ in size.",
				"Retrained uniform/metric controls in main comparison answer different questions from these interventions.",
			],
		),
	)
	log(
		"DONE",
		"attention_audit",
		f"n={len(ids)}, layers={maps.shape[1]}, heads={maps.shape[2]}",
	)


# Evaluation
# Evaluation starts only from frozen predictions and model metadata.


def paired_intervals(y, w, prediction, masks, groups, primary, count, seed):
	contrasts = [
		(primary, name)
		for name in [
			"elasticnet",
			"protein_elasticnet",
			"transformer_random_panel",
			"transformer_diversity_panel",
			"transformer_fullbank",
			"pca_same_panel",
			"embedding_knn_same_panel",
			"mlp_same_panel",
			"uniform_same_panel",
			"metric_same_panel",
			"equal_heads_same_panel",
			"no_pretrain_same_panel",
			"ssl_same_panel",
			"permuted_values_same_panel",
			"random_donors_same_panel",
			"copy1_reliable_panel",
		]
		if name in prediction
	]
	contrasts += [("abm_clinical", "clinical"), ("abm_clinical", "elasticnet")]
	contrasts += (
		[("transformer_fullbank", "mlp_fullbank")]
		if "mlp_fullbank" in prediction
		else []
	)
	rng = np.random.default_rng(seed)
	rows = []
	for subset, mask in masks.items():
		idx = np.flatnonzero(mask)
		if len(idx) < 30 or len(np.unique(y[idx][w[idx] > 0])) < 2:
			continue
		units = (
			[idx[groups[idx] == g] for g in np.unique(groups[idx])]
			if groups is not None
			else None
		)
		names = sorted(set(v for pair in contrasts for v in pair))
		point = {name: metrics(y[idx], w[idx], prediction[name][idx]) for name in names}
		draws = {
			(m, b, key): []
			for m, b in contrasts
			for key in ["AUC_IPCW", "Brier_IPCW", "LogLoss_IPCW"]
		}
		for _ in range(count):
			take = (
				np.concatenate(
					[units[j] for j in rng.integers(0, len(units), len(units))]
				)
				if units is not None
				else rng.choice(idx, len(idx), replace=True)
			)
			value = {
				name: metrics(y[take], w[take], prediction[name][take])
				for name in names
			}
			for (m, b, key), vals in draws.items():
				delta = value[m][key] - value[b][key]
				if np.isfinite(delta):
					vals.append(delta)
		for (m, b, key), vals in draws.items():
			ci = (
				np.quantile(vals, [0.025, 0.975])
				if len(vals) >= 20
				else [np.nan, np.nan]
			)
			rows.append(
				dict(
					subset=subset,
					model=m,
					reference=b,
					metric=key,
					delta=point[m][key] - point[b][key],
					lower=ci[0],
					upper=ci[1],
					replicates=len(vals),
					uncertainty="conditional_on_fitted_model; exploratory_unadjusted_intervals",
				)
			)
	return pd.DataFrame(rows)


def evaluate_frozen(out):
	if not (out / "MODEL_FROZEN.json").exists():
		raise ValueError("Missing frozen-model record")
	bundle = joblib.load(out / "model_bundle.joblib")
	cfg, primary = bundle["config"], bundle["primary"]
	people = pd.read_csv(out / "test_outcomes.csv", dtype={cfg["id_col"]: str})
	info = pd.read_csv(out / "test_individuals.csv", dtype={cfg["id_col"]: str})
	if not people[cfg["id_col"]].equals(info[cfg["id_col"]]):
		raise ValueError("Test outcome/prediction ID mismatch")
	names = json.loads((out / "prediction_columns.json").read_text())
	predictions = {name: info[name].to_numpy() for name in names}
	y, w = bundle["censoring"].labels_weights(people)
	accepted = info.supported_match.to_numpy(bool)
	masks = {"all": np.ones(len(y), bool), "supported": accepted, "rejected": ~accepted}
	rows = []
	for subset, mask in masks.items():
		for name, prob in predictions.items():
			rows.append(
				dict(
					model=name,
					subset=subset,
					n=int(mask.sum()),
					known=int((w[mask] > 0).sum()),
					cases=int(y[mask].sum()),
					coverage=float(mask.mean()),
					**metrics(y[mask], w[mask], prob[mask]),
				)
			)
	result = pd.DataFrame(rows)
	result.to_csv(out / "test_metrics.csv", index=False)
	groups = people[cfg["group_col"]].to_numpy(str) if cfg["group_col"] else None
	paired_intervals(
		y,
		w,
		predictions,
		{key: masks[key] for key in ["all", "supported", "rejected"]},
		groups,
		primary,
		cfg["bootstrap"],
		cfg["seed"] + 823,
	).to_csv(out / "paired_contrasts.csv", index=False)
	curves = []
	for quantile, rule in bundle["coverage_rules"].items():
		mask, _ = rule.apply(
			info,
			info.missing_fraction.to_numpy(),
			info.unknown_category.to_numpy(bool),
			cfg["sample_missing"],
		)
		for name, prob in predictions.items():
			curves.append(
				dict(
					quantile=quantile,
					model=name,
					coverage=float(mask.mean()),
					n=int(mask.sum()),
					known=int((w[mask] > 0).sum()),
					cases=int(y[mask].sum()),
					threshold_source="frozen on tune; test outcomes not used for support selection",
					**metrics(y[mask], w[mask], prob[mask]),
				)
			)
	pd.DataFrame(curves).to_csv(out / "coverage_curve.csv", index=False)
	# Compare expected-error scores with actual squared error without claiming
	# that lower outcome prevalence alone proves better matching.
	errbins = pd.qcut(info.estimated_squared_error, 5, labels=False, duplicates="drop")
	error_rows = []
	for b in sorted(errbins.dropna().unique()):
		mask = errbins.eq(b).to_numpy()
		error_rows.append(
			dict(
				bin=int(b) + 1,
				n=int(mask.sum()),
				expected_error=float(info.loc[mask, "estimated_squared_error"].mean()),
				observed_error=float(
					np.average(
						(y[mask] - predictions[primary][mask]) ** 2, weights=w[mask]
					)
				),
				event_rate_IPCW=float(np.average(y[mask], weights=w[mask])),
				brier_gain_vs_elasticnet=float(
					np.mean(
						w[mask]
						* (
							(y[mask] - predictions["elasticnet"][mask]) ** 2
							- (y[mask] - predictions[primary][mask]) ** 2
						)
					)
				),
			)
		)
	pd.DataFrame(error_rows).to_csv(out / "support_error_audit.csv", index=False)
	calibration = []
	for name, prob in predictions.items():
		bins = pd.qcut(pd.Series(prob), 5, labels=False, duplicates="drop")
		for b in sorted(bins.dropna().unique()):
			mask = bins.eq(b).to_numpy()
			calibration.append(
				dict(
					model=name,
					bin=int(b) + 1,
					n=int(mask.sum()),
					predicted=float(prob[mask].mean()),
					observed_IPCW=(
						float(np.average(y[mask], weights=w[mask]))
						if w[mask].sum() > 0
						else np.nan
					),
				)
			)
	pd.DataFrame(calibration).to_csv(out / "test_calibration.csv", index=False)
	demographic = pd.read_csv(out / "test_demographics.csv")
	if "age" in demographic:
		demographic["age_group"] = pd.cut(
			demographic.age, [0, 50, 60, 70, 150], right=False
		).astype(str)
	subgroup = []
	for column in [c for c in ["sex", "center", "age_group"] if c in demographic]:
		for level, sub in demographic.groupby(column, dropna=False):
			ix = sub.index.to_numpy()
			for name in [primary, "elasticnet", "abm_clinical"]:
				subgroup.append(
					dict(
						variable=column,
						level=level,
						model=name,
						n=len(ix),
						cases=int(y[ix].sum()),
						supported_fraction=float(accepted[ix].mean()),
						**metrics(y[ix], w[ix], predictions[name][ix]),
					)
				)
	pd.DataFrame(subgroup).to_csv(out / "subgroup_metrics.csv", index=False)
	raw_rows = []
	for name in names:
		col = name + "_raw"
		if col in info:
			raw_rows.append(dict(model=name, **metrics(y, w, info[col].to_numpy())))
	pd.DataFrame(raw_rows).to_csv(out / "raw_test_metrics.csv", index=False)
	compare = result[result.subset.eq("all")].copy()
	for name, bank in bundle["banks"].items():
		ix = compare.model.eq(name)
		compare.loc[ix, "reference_count"] = len(bank.ids)
		compare.loc[ix, "reference_case_fraction"] = bank.y.mean()
	compare.to_csv(out / "approach_comparison.csv", index=False)
	write_report(out, bundle, result, info)
	figures(out, bundle, result)
	dump(
		out / "DONE.json",
		dict(
			version=VERSION,
			n_test=len(y),
			cases_by_horizon=int(y.sum()),
			supported_fraction=float(accepted.mean()),
			panel_status=bundle["readiness"]["status"],
			results_synthetic=cfg["demo"],
		),
	)
	log(
		"DONE",
		"evaluation",
		f"N={len(y)}; horizon events={y.sum()}; support={accepted.mean():.1%}",
	)


def write_report(out, bundle, result, info):
	primary = bundle["primary"]
	all_rows = result[result.subset.eq("all")].sort_values("AUC_IPCW", ascending=False)
	text = [
		"# ABM 5 run report",
		"",
		(
			"SYNTHETIC demonstration; not UKB evidence."
			if bundle["config"]["demo"]
			else "Research results on the supplied input cohort."
		),
		"",
		f"Prespecified primary: `{primary}`. Reference readiness: `{bundle['readiness']['status']}`.",
		f"Readiness reasons: {bundle['readiness']['reasons']}",
		"",
		"| Model | AUC IPCW | Brier IPCW |",
		"|---|---:|---:|",
	]
	for row in all_rows.itertuples():
		text.append(f"| {row.model} | {row.AUC_IPCW:.4f} | {row.Brier_IPCW:.5f} |")
	text += [
		"",
		"## Questions to answer before claiming a useful individual reference model",
		"",
		"1. Does the learned encoder + copying beat same-reference PCA, Euclidean copying, MLP, and random reference controls?",
		"2. Does a 100-person panel preserve useful signal from the full bank? Does class stratification prevent constant COPY?",
		"3. Do reference label permutation and fixed-neighbour deletion measurably change predictions? A head-only predictor is not evidence of useful outcome borrowing.",
		"4. Does masked reconstruction beat a zero build-mean profile and PCA reconstruction on exactly the same masked entries? Stable but constant predictions are not informative explanations.",
		"5. Does support improve error relative to a comparator, or mostly reject higher-risk people? Read subgroup_metrics and support_error_audit together.",
		"6. Are same-risk molecular contrasts reproducible across seeds/splits/residualization, and externally validated? Token names are not biological pathway claims.",
		"",
		"## Individual interpretation",
		"",
		"The uncalibrated estimate is exactly the sum of nonnegative reference outcome contributions and an explicit build-prior contribution. Monotone calibration follows this sum. COPY1 outputs an observed binary outcome, not a person's certain future or a probability estimate.",
		"Deleting a reference holds the retrieved neighbour set fixed and renormalizes the remaining weights. Module masking changes information availability; neither operation is a causal intervention.",
		"Reference identities alone do not establish explanatory validity. Use the recorded contribution identity, deletion sensitivity, missing-protein reconstruction and controls.",
		"",
		"## Validation limits",
		"",
		"The outer test half is not used for fitting, model selection, calibration, panel utility or support thresholds. Calibration is split into fitting and independent audit halves. Inner neural validation is only an early-stopping device; OOF reference-fit scores come from separate classical cross-fitting.",
		"Report all prespecified models. Test-set ranking is exploratory: selecting a winner after seeing this report needs a new untouched cohort. Bootstrap intervals condition on the fitted model and are unadjusted for multiple contrasts; use the repeat runner for training instability.",
		"A readiness pass is an audit point-estimate screen, not proof of individual accuracy. A failed panel does not invalidate the source cohort. Death is censored, so the endpoint is fixed-horizon net risk, not a competing-risk cumulative incidence. IPCW uses a marginal build-set censoring model and assumes independent censoring.",
		"",
	]
	(out / "REPORT.md").write_text("\n".join(text), encoding="utf-8")


def figures(out, bundle, result):
	import matplotlib

	matplotlib.use("Agg")
	import matplotlib.pyplot as plt

	main = result[result.subset.eq("all")].sort_values("AUC_IPCW")
	fig, ax = plt.subplots(figsize=(10, max(6, 0.27 * len(main))))
	colors = ["#c04a31" if s == bundle["primary"] else "#4a718e" for s in main.model]
	ax.barh(main.model, main.AUC_IPCW, color=colors)
	ax.axvline(0.5, color="grey", ls="--")
	ax.set_xlim(0, 1)
	ax.set_xlabel("Full test IPCW AUC")
	fig.tight_layout()
	fig.savefig(out / "Fig_model_comparison.png", dpi=180)
	plt.close(fig)
	if (out / "masked_reconstruction.csv").exists():
		dat = (
			pd.read_csv(out / "masked_reconstruction.csv")
			.groupby("method")
			.masked_RMSE.mean()
			.sort_values()
		)
		fig, ax = plt.subplots(figsize=(8, 4))
		ax.barh(dat.index, dat.values, color="#4a718e")
		ax.set_xlabel("Mean person-level masked RMSE (lower is better)")
		fig.tight_layout()
		fig.savefig(out / "Fig_masked_reconstruction.png", dpi=180)
		plt.close(fig)
	dat = pd.read_csv(out / "coverage_curve.csv")
	fig, ax = plt.subplots(figsize=(7, 4))
	for name in [bundle["primary"], "elasticnet", "abm_clinical"]:
		sub = dat[dat.model.eq(name)]
		ax.plot(sub.coverage, sub.Brier_IPCW, "o-", label=name)
	ax.set_xlabel("Supported fraction; thresholds frozen on tune")
	ax.set_ylabel("IPCW Brier on identical people")
	ax.legend(fontsize=8)
	fig.tight_layout()
	fig.savefig(out / "Fig_coverage.png", dpi=180)
	plt.close(fig)


# Pipeline
# Train / freeze / evaluate / project with an untouched outer half-cohort.


REFERENCE_PRIMARY = "abm_transformer"


def split_people(p, a):
	if a.split_file:
		table = read_table(a.split_file, a.id_col)
		if "split" not in table or not set(table.split) <= {
			"build",
			"tune",
			"calibration",
			"test",
		}:
			raise ValueError("Split file requires build/tune/calibration/test labels")
		part = table.set_index(a.id_col).reindex(p[a.id_col]).split
		if part.isna().any():
			raise ValueError("Split file must cover all eligible participants")
		part = part.to_numpy(str)
	else:

		def divide(ix, fraction, seed):
			if a.group_col:
				one, two = next(
					GroupShuffleSplit(1, test_size=fraction, random_state=seed).split(
						ix, groups=p.iloc[ix][a.group_col].to_numpy(str)
					)
				)
				return ix[one], ix[two]
			labels = (p.iloc[ix].event.eq(1) & p.iloc[ix].time.le(a.horizon)).to_numpy(
				int
			)
			return train_test_split(
				ix, test_size=fraction, random_state=seed, stratify=labels
			)

		develop, test = divide(np.arange(len(p)), 0.5, a.seed)
		bt, cal = divide(develop, 0.2, a.seed + 1)
		build, tune = divide(bt, 0.25, a.seed + 2)
		part = np.full(len(p), "test", object)
		part[build], part[tune], part[cal] = "build", "tune", "calibration"
	if set(part) != {"build", "tune", "calibration", "test"}:
		raise ValueError("All four outer partitions are required")
	if a.group_col:
		if (
			p[a.group_col].isna().any()
			or p[a.group_col].astype(str).str.strip().eq("").any()
		):
			raise ValueError("Complete family component IDs required")
		if (
			pd.DataFrame({"group": p[a.group_col].to_numpy(str), "split": part})
			.groupby("group")
			.split.nunique()
			.max()
			> 1
		):
			raise ValueError("Family components cross split boundaries")
	return part


def calibration_halves(p, mask, idcol, groupcol, seed):
	# Outcome-blind, stable under row order. Families stay together.
	unit = p[groupcol if groupcol else idcol].astype(str).to_numpy()
	levels = sorted(set(unit[mask]), key=lambda s: stable_seed(s, seed + 832))
	if len(levels) < 2:
		raise ValueError("Need two independent calibration units")
	audit_units = set(levels[: len(levels) // 2])
	audit = mask & np.array([s in audit_units for s in unit])
	return mask & ~audit, audit


def unknown_rows(p, designs):
	unknown = np.zeros(len(p), bool)
	for design in designs:
		if design is None:
			continue
		norm = design.normalized(p)
		for col, levels in design.levels.items():
			unknown |= (norm[col].notna() & ~norm[col].isin(levels)).to_numpy()
	return unknown


def runtime_manifest(a):
	deps = {}
	for name in [
		"numpy",
		"pandas",
		"scipy",
		"scikit-learn",
		"joblib",
		"lightgbm",
		"torch",
		"pyreadr",
	]:
		try:
			deps[name] = importlib.metadata.version(name)
		except importlib.metadata.PackageNotFoundError:
			deps[name] = "not_installed"
	paths = [] if a.demo else [a.phe_file, a.omics_file]
	paths += [s for s in [a.split_file, a.module_file] if s]
	if a.biom == "met" and a.met_input == "raw" and not a.demo:
		paths.append(a.met_map)
	config = {
		k: v for k, v in vars(a).items() if k not in ["resume", "replace", "train_only"]
	}
	data = dict(
		version=VERSION,
		config=config,
		dependencies=deps,
		inputs=fingerprints(paths, a.full_input_hash),
		code=fingerprints([Path(__file__)], True),
	)
	data["signature"] = digest(data)
	return data


def reference_train(a, out):
	out = Path(out)
	manifest = runtime_manifest(a)
	if a.resume and (out / "manifest.json").exists():
		old = json.loads((out / "manifest.json").read_text())
		if manifest["signature"] != old["signature"]:
			raise ValueError(
				"Resume requires identical code, data fingerprints, dependencies and configuration"
			)
	dump(out / "manifest.json", manifest)
	prepared = abm_cache_dir(out, "input")
	prepared.mkdir(exist_ok=True)
	log("START", "prepare")
	p = prepare(a, prepared)
	raw = np.load(prepared / "raw.npy")
	features = (prepared / "features.txt").read_text().splitlines()
	p, audit = outcomes(p, a)
	eligible = p.eligible.to_numpy(bool)
	raw, p = raw[eligible], p.loc[eligible].reset_index(drop=True)
	p["split"] = split_people(p, a)
	p[[a.id_col, "split"]].to_csv(out / "split_before_qc.csv", index=False)
	build = p.split.eq("build").to_numpy()
	proto = MolecularPreprocessor(
		a.feature_missing, words(a.residualize), words(a.categorical), a.transform
	)
	scaled = proto.scale_transform(raw[build])
	keep = (np.mean(~np.isfinite(scaled), axis=0) < a.feature_missing) & (
		np.nanstd(scaled, axis=0) > 1e-8
	)
	if keep.sum() < 3:
		raise ValueError("Fewer than three usable assays")
	missing = np.mean(~np.isfinite(proto.scale_transform(raw[:, keep])), axis=1)
	qc = missing <= a.sample_missing
	pd.DataFrame(
		{
			a.id_col: p[a.id_col],
			"split": p.split,
			"missing_fraction": missing,
			"included": qc,
		}
	).to_csv(out / "sample_qc.csv", index=False)
	raw, p, missing = raw[qc], p.loc[qc].reset_index(drop=True), missing[qc]
	build, tune, cal, test = [
		p.split.eq(s).to_numpy() for s in ["build", "tune", "calibration", "test"]
	]
	calfit, cala = calibration_halves(p, cal, a.id_col, a.group_col, a.seed)
	if min(build.sum(), tune.sum(), calfit.sum(), cala.sum(), test.sum()) < 30:
		raise ValueError("Too few people after QC / calibration-audit separation")
	p["role"] = p.split
	p.loc[calfit, "role"], p.loc[cala, "role"] = "calibration_fit", "calibration_audit"
	p[[a.id_col, "split", "role"]].to_csv(out / "split.csv", index=False)
	if a.shuffle_development_outcomes:
		rng = np.random.default_rng(a.seed + 912)
		for role in ["build", "tune", "calibration_fit", "calibration_audit"]:
			ix = np.flatnonzero(p.role.eq(role))
			permuted = rng.permutation(ix)
			p.loc[ix, ["time", "event"]] = p.loc[permuted, ["time", "event"]].to_numpy()
		log(
			"CONTROL",
			"development_outcomes",
			"Joint time/event pairs permuted within each development role; not a permutation p-value",
		)
	audit.update(
		n_after_qc=len(p),
		sample_missing_exclusions=int((~qc).sum()),
		retained_features=int(keep.sum()),
		group_split=bool(a.group_col),
		role_counts=p.role.value_counts().to_dict(),
		test_Y_not_used_after_initial_stratification=True,
		development_outcomes_shuffled=a.shuffle_development_outcomes,
	)
	dump(out / "cohort_audit.json", audit)
	prep = proto.fit(raw[build], p.loc[build], allowed=keep)
	x, observed = prep.transform(raw, p)
	clinical = MetadataDesign(words(a.covariates), words(a.categorical)).fit(
		p.loc[build]
	)
	c = clinical.transform(p)
	if c.shape[1] == 0:
		raise ValueError("At least one clinical comparator covariate required")
	tech = MetadataDesign(
		list(dict.fromkeys(words(a.covariates) + words(a.residualize))),
		words(a.categorical),
	).fit(p.loc[build])
	unknown = unknown_rows(p, [clinical, getattr(prep, "design", None)])
	km = CensoringKM().fit(p.loc[build], a.horizon, a.min_censor_survival)
	y, w = np.zeros(len(p), int), np.zeros(len(p))
	y[~test], w[~test] = km.labels_weights(p.loc[~test])
	for name, mask in [
		("build", build),
		("tune", tune),
		("calibration_fit", calfit),
		("calibration_audit", cala),
	]:
		if (
			np.sum(y[mask]) < a.min_events
			or np.sum((w[mask] > 0) & (y[mask] == 0)) < a.min_events
		):
			raise ValueError(
				f"Too few known cases/controls in {name}; change horizon/data or --min-events for a demo"
			)
	ids = p[a.id_col].to_numpy(str)
	groups = p[a.group_col].to_numpy(str) if a.group_col else None
	bg = None if groups is None else groups[build]
	tg = None if groups is None else groups[tune]
	f_names = [features[j] for j in prep.keep]
	log(
		"DONE",
		"prepare",
		f"N={len(p)}, assays={x.shape[1]}, build={build.sum()}, test={test.sum()}",
	)
	quality, importance = oof_quality(
		raw[build], p.loc[build].reset_index(drop=True), features, a, out
	)
	geometry = MolecularGeometry().fit(
		x[build], importance[prep.keep], a.dimensions, a.seed
	)
	pca_z = (
		geometry.unsupervised.transform(x)
		/ geometry.scale_u
		/ np.sqrt(len(geometry.scale_u))
	)
	membership, token_names = token_partition(
		x[build], f_names, a.tokens, a.seed, a.module_file
	)
	pd.DataFrame(
		{"feature": f_names, "token": [token_names[j] for j in membership]}
	).to_csv(out / "token_membership.csv", index=False)
	models, rawpred, model_tuning = {}, {}, []
	blocks = {
		"clinical": c,
		"elasticnet": np.c_[c, x],
		"protein_elasticnet": x,
		"clinical_pca": np.c_[c, pca_z],
		"clinical_technical": np.c_[tech.transform(p), 1 - observed.mean(1)],
	}
	log("START", "classical_comparators")
	for name, block in blocks.items():
		model, rows = tune_logistic(
			block, y, w, build, tune, name, a, 0.5 if "elasticnet" in name else 0
		)
		models[name] = model
		model_tuning.extend(rows)
		rawpred[name] = model.predict_proba(block)[:, 1]
	if a.tree != "none":
		best = None
		for leaves in [7, 15, 31]:
			if a.tree == "lightgbm":
				from lightgbm import LGBMClassifier

				obj = LGBMClassifier(
					n_estimators=a.tree_estimators,
					num_leaves=leaves,
					learning_rate=0.03,
					min_child_samples=40,
					reg_lambda=10,
					random_state=a.seed,
					n_jobs=a.cores,
					verbosity=-1,
					deterministic=True,
					force_col_wise=True,
				)
			else:
				from sklearn.ensemble import HistGradientBoostingClassifier

				obj = HistGradientBoostingClassifier(
					max_iter=a.tree_estimators,
					max_leaf_nodes=leaves,
					learning_rate=0.03,
					min_samples_leaf=40,
					l2_regularization=10,
					early_stopping=False,
					random_state=a.seed,
				)
			known = build & (w > 0)
			obj.fit(
				blocks["elasticnet"][known],
				y[known],
				sample_weight=w[known] / w[known].mean(),
			)
			pred = obj.predict_proba(blocks["elasticnet"][tune])[:, 1]
			score = float(np.average(loss(y[tune], pred), weights=w[tune]))
			model_tuning.append(
				dict(model=a.tree, leaves=leaves, validation_logloss=score)
			)
			if best is None or score < best[0]:
				best = score, obj
		models[a.tree] = best[1]
		rawpred[a.tree] = best[1].predict_proba(blocks["elasticnet"])[:, 1]
	pd.DataFrame(model_tuning).to_csv(out / "comparator_tuning.csv", index=False)
	# Independent models answer whether attention, pretraining, and retrieval help.
	encoders, embeddings = {}, {}
	for name in a.experiments:
		objective = "direct" if name == "direct" else "retrieval"
		architecture = name if name in ["mlp", "uniform"] else "transformer"
		log("START", "neural_" + name)
		obj, ssl = train_encoder(
			x[build],
			observed[build],
			y[build],
			w[build],
			ids[build],
			bg,
			membership,
			a,
			abm_cache_dir(out, "neural", name),
			architecture,
			objective,
			pretrain=name != "no_pretrain",
			cross_mode="metric" if name == "metric" else "qkv",
		)
		encoders[name] = obj
		z, direct, _ = encode(obj, x, observed, device_for(a.device), a.batch_size)
		obj.cpu()
		embeddings[name] = z
		rawpred[name + "_direct_head"] = direct
		if name == "transformer" and ssl is not None:
			encoders["ssl"] = ssl
			embeddings["ssl"] = encode(
				ssl, x, observed, device_for(a.device), a.batch_size
			)[0]
			ssl.cpu()
		log("DONE", "neural_" + name)
	primary_z = embeddings["transformer"]
	embedding_audit = []
	for name, zz in embeddings.items():
		eigen = np.maximum(np.linalg.eigvalsh(np.cov(zz[build].T)), 0)
		embedding_audit.append(
			dict(
				model=name,
				mean_coordinate_SD=float(zz[build].std(0).mean()),
				effective_rank=float(eigen.sum() ** 2 / max(np.sum(eigen**2), 1e-30)),
				dimensions=zz.shape[1],
				warning="representation_collapse" if eigen.sum() < 1e-6 else "",
			)
		)
	pd.DataFrame(embedding_audit).to_csv(out / "embedding_diagnostics.csv", index=False)
	xb, ob, yb, wb, bid = x[build], observed[build], y[build], w[build], ids[build]
	known_indices = np.flatnonzero(wb > 0)
	primary_anchor = select_references(
		primary_z[build], quality, a.panel_size, "reliable", a.seed
	)
	qualified = (
		len(primary_anchor) == a.panel_size and len(np.unique(yb[primary_anchor])) == 2
	)
	fallback = len(primary_anchor) < 2
	if fallback:
		primary_anchor = select_references(
			primary_z[build], quality, a.panel_size, "diversity", a.seed
		)
	banks, specs, tuning, direct_map = {}, {}, [], {}

	def add_bank(
		name,
		encoder_name,
		anchor,
		kernel="attention",
		literal=False,
		balanced=False,
		z_override=None,
		full=False,
	):
		if len(anchor) < 2:
			log("SKIP", name, "fewer than 2 reference candidates")
			return
		zz = embeddings[encoder_name] if z_override is None else z_override
		encoder = (
			encoders.get(encoder_name)
			if kernel in ["attention", "equal_heads"]
			else None
		)
		bank = ReferenceBank(
			zz[build],
			xb,
			ob,
			bid,
			bg,
			yb,
			wb,
			anchor,
			encoder,
			inclusion_correction=balanced,
			seed=a.seed,
		)
		ks = [1] if literal else a.all_k if full else a.match_k
		ks = sorted(set(min(k, len(anchor)) for k in ks))
		# Scalar temperature candidates matter: unlike v4, the bandwidth does not
		# force every chosen donor into a narrow, nearly uniform weight interval.
		chosen, rows = tune_bank(
			bank,
			zz[tune],
			ids[tune],
			tg,
			y[tune],
			w[tune],
			ks,
			[1.0] if literal else a.temperatures,
			[0.0] if literal else a.prior_grid,
			name,
			kernel,
			literal,
		)
		specs[name] = dict(
			encoder=encoder_name,
			**chosen,
			space=(
				"full_proteome"
				if z_override is x
				else "pca"
				if z_override is pca_z
				else "neural"
			),
		)
		banks[name] = bank
		tuning.extend(rows)
		log(
			"TUNED",
			name,
			f"references={len(bank.ids)}; k={chosen['k']}; temperature={chosen['temperature']}; prior={chosen['strength']}",
		)
		og = None if groups is None else groups[~build]
		pred, _ = bank.match(zz[~build], ids[~build], og, **chosen)
		rawpred[name] = np.full(len(p), np.nan)
		rawpred[name][~build] = pred
		pd.DataFrame(
			{
				a.id_col: bank.ids,
				"horizon_label": bank.y,
				"fit_gain": quality.iloc[anchor].fit_gain.to_numpy(),
				"OOF_probability": quality.iloc[anchor].OOF_probability.to_numpy(),
				"class_sampling_correction": bank.correction,
			}
		).to_csv(out / f"panel_{name}.csv", index=False)

	log("START", "reference_panels")
	add_bank(REFERENCE_PRIMARY, "transformer", primary_anchor)
	for mode in ["topfit", "stratified_fit", "random", "diversity"]:
		anchor = select_references(
			primary_z[build], quality, a.panel_size, mode, a.seed
		)
		if mode in ["topfit", "stratified_fit"]:
			add_bank("copy1_" + mode, "transformer", anchor, literal=True)
		add_bank("transformer_" + mode + "_panel", "transformer", anchor)
	topfit = select_references(
		primary_z[build], quality, a.panel_size, "topfit", a.seed
	)
	add_bank(
		"copy1_fullproteome",
		"none",
		topfit,
		kernel="euclidean",
		literal=True,
		z_override=x,
	)
	add_bank("copy1_reliable_panel", "transformer", primary_anchor, literal=True)
	balanced = select_references(
		primary_z[build], quality, a.panel_size, "reliable", a.seed, balanced=True
	)
	add_bank("transformer_balanced_panel", "transformer", balanced, balanced=True)
	for size in a.panel_sizes:
		if size != a.panel_size:
			add_bank(
				f"transformer_reliable{size}",
				"transformer",
				select_references(primary_z[build], quality, size, "reliable", a.seed),
			)
	add_bank("transformer_fullbank", "transformer", known_indices, full=True)
	full = banks["transformer_fullbank"]
	utility, utility_table = reference_utility(
		full, primary_z[tune], ids[tune], tg, y[tune], w[tune]
	)
	utility_table.to_csv(out / "tuning_reference_utility.csv", index=False)
	all_utility = np.zeros(build.sum())
	all_utility[full.indices] = utility
	utility_anchor = select_references(
		primary_z[build], quality, a.panel_size, "utility", a.seed, all_utility
	)
	add_bank("transformer_utility_panel", "transformer", utility_anchor)
	add_bank(
		"equal_heads_same_panel", "transformer", primary_anchor, kernel="equal_heads"
	)
	add_bank(
		"embedding_knn_same_panel", "transformer", primary_anchor, kernel="euclidean"
	)
	add_bank(
		"pca_same_panel", "none", primary_anchor, kernel="euclidean", z_override=pca_z
	)
	add_bank(
		"pca_fullbank",
		"none",
		known_indices,
		kernel="euclidean",
		z_override=pca_z,
		full=True,
	)
	for name in ["mlp", "no_pretrain", "ssl", "uniform", "metric"]:
		if name in encoders:
			add_bank(
				name + "_same_panel",
				name,
				primary_anchor,
				kernel="euclidean" if name == "ssl" else "attention",
			)
			if name == "mlp":
				add_bank("mlp_fullbank", name, known_indices, full=True)
	primary = banks[REFERENCE_PRIMARY]
	primary_spec = specs[REFERENCE_PRIMARY]
	match_kw = {k: v for k, v in primary_spec.items() if k not in ["encoder", "space"]}
	og = None if groups is None else groups[~build]
	for name, change in [
		("random_donors_same_panel", {"random_candidates": True}),
		("permuted_values_same_panel", {"permuted": True}),
	]:
		specs[name] = dict(primary_spec, **change)
		banks[name] = primary
		rawpred[name] = np.full(len(p), np.nan)
		rawpred[name][~build] = primary.match(
			primary_z[~build], ids[~build], og, **match_kw, **change
		)[0]
	primary_raw, detail = primary.match(primary_z[~build], ids[~build], og, **match_kw)
	info = primary.describe(x[~build], observed[~build], primary_raw, detail)
	info.insert(0, a.id_col, ids[~build])
	info["role"] = p.role.to_numpy()[~build]
	log("START", "module_mosaic")
	mosaic = ModuleMosaic().fit(
		xb, ob, bid, bg, yb, wb, primary_anchor, membership, token_names, a.seed
	)
	module_risk, module_reference = mosaic.components(x[~build], ids[~build], og)
	mosaic.tune(module_risk[tune[~build]], y[tune], w[tune])
	for name, risk in [
		("mosaic_equal", module_risk.mean(1)),
		("mosaic_weighted", module_risk @ mosaic.weights),
	]:
		rawpred[name] = np.full(len(p), np.nan)
		rawpred[name][~build] = risk
	pd.DataFrame({"token": token_names, "mixture_weight": mosaic.weights}).to_csv(
		out / "mosaic_weights.csv", index=False
	)
	pd.DataFrame(tuning).to_csv(out / "panel_tuning.csv", index=False)
	dump(out / "panel_specs.json", specs)
	# Tune set selects candidate hyperparameters. Calibration and auditing are disjoint.
	calibrators, predictions = {}, {}
	for name, risk in rawpred.items():
		if name.startswith("copy1_"):
			predictions[name] = risk.copy()
		else:
			calibrators[name] = RiskCalibrator().fit(risk[calfit], y[calfit], w[calfit])
			predictions[name] = calibrators[name].predict(risk)
	combinations = {
		"abm_clinical": (REFERENCE_PRIMARY, "clinical"),
		"abm_elasticnet_hybrid": (REFERENCE_PRIMARY, "elasticnet"),
		"mosaic_clinical": ("mosaic_weighted", "clinical"),
	}
	for name, (source, covariate) in combinations.items():
		obj = RiskCalibrator().fit(
			rawpred[source][calfit], y[calfit], w[calfit], rawpred[covariate][calfit]
		)
		calibrators[name] = obj
		predictions[name] = obj.predict(rawpred[source], rawpred[covariate])
	dump(
		out / "calibration_parameters.json",
		{
			name: dict(coefficients=obj.coef, logit_mean=obj.mean, logit_scale=obj.sd)
			for name, obj in calibrators.items()
		},
	)
	rule = SupportRule().fit(
		info.loc[tune[~build]],
		predictions[REFERENCE_PRIMARY][tune],
		y[tune],
		w[tune],
		a.accept_quantile,
		a.min_match_ess,
	)
	supported, reasons = rule.apply(
		info, missing[~build], unknown[~build], a.sample_missing
	)
	info["supported_match"], info["rejection_reason"] = supported, reasons
	info["estimated_squared_error"] = rule.error_score(info)
	info["missing_fraction"], info["unknown_category"] = (
		missing[~build],
		unknown[~build],
	)
	dev_rows = []
	for subset, mask in [("tune", tune), ("calibration_audit", cala)]:
		for name, risk in predictions.items():
			dev_rows.append(
				dict(subset=subset, model=name, **metrics(y[mask], w[mask], risk[mask]))
			)
	dev_rows = pd.DataFrame(dev_rows)
	dev_rows.to_csv(out / "development_metrics.csv", index=False)
	primary_met = metrics(y[cala], w[cala], predictions[REFERENCE_PRIMARY][cala])
	constant = float(np.average(y[calfit], weights=w[calfit]))
	null = metrics(y[cala], w[cala], np.full(cala.sum(), constant))
	ready_reasons = []
	if not qualified:
		ready_reasons.append("requested_reliable_panel_incomplete_or_one_class")
	if supported[cala[~build]].mean() < a.min_coverage:
		ready_reasons.append("low_audit_coverage")
	if primary_met["Brier_IPCW"] >= null["Brier_IPCW"]:
		ready_reasons.append("no_independent_audit_Brier_gain")
	if primary_met["LogLoss_IPCW"] >= null["LogLoss_IPCW"]:
		ready_reasons.append("no_independent_audit_logloss_gain")
	if not np.isfinite(primary_met["AUC_IPCW"]) or primary_met["AUC_IPCW"] <= 0.5:
		ready_reasons.append("no_independent_audit_discrimination")
	readiness = dict(
		status="provisional_pass" if not ready_reasons else "not_ready",
		reasons=ready_reasons,
		reliable_candidates=int(quality.reliable_candidate.sum()),
		requested_size=a.panel_size,
		actual_size=len(primary.ids),
		primary_model=REFERENCE_PRIMARY,
		fallback_diversity=fallback,
		audit_metrics=primary_met,
		audit_null=null,
		audit_coverage=float(supported[cala[~build]].mean()),
		interpretation="independent audit point-estimate screen; not an individual accuracy guarantee or cohort validity test",
	)
	dump(out / "reference_readiness.json", readiness)
	dump(
		out / "support_thresholds.json",
		dict(
			thresholds=rule.thresholds,
			min_ess=rule.min_ess,
			error_threshold=rule.error_threshold,
			source="tune; error model uses tune Y; audit is separate",
		),
	)
	for name, value in predictions.items():
		info[name] = value[~build]
	raw_frame = pd.DataFrame(
		{name + "_raw": value[~build] for name, value in rawpred.items()}
	)
	info = pd.concat([info.reset_index(drop=True), raw_frame], axis=1)
	info["raw_borrowed_risk"] = primary_raw
	info["reference_ready"] = readiness["status"] == "provisional_pass"
	info["prediction_released"] = info.supported_match & info.reference_ready
	info["released_net_risk"] = np.where(
		info.prediction_released, info[REFERENCE_PRIMARY], np.nan
	)
	info.to_csv(
		out / "development_and_test_individuals.csv.gz", index=False, compression="gzip"
	)
	coverage_rules = {}
	for quantile in [0.5, 0.75, 0.9, 0.95, 0.99, 1.0]:
		alternate = copy.copy(rule)
		alternate.thresholds = {
			key: float(info.loc[tune[~build], key].quantile(quantile))
			for key in rule.thresholds
		}
		alternate.error_threshold = float(
			np.quantile(rule.error_score(info.loc[tune[~build]]), quantile)
		)
		coverage_rules[quantile] = alternate
	profile_build = np.column_stack(
		[xb[:, membership == j].mean(1) for j in range(len(token_names))]
	)
	profile_center, profile_scale = (
		profile_build.mean(0),
		np.maximum(profile_build.std(0), 1e-6),
	)
	bundle = dict(
		version=VERSION,
		config=vars(a),
		features=features,
		retained_features=f_names,
		prep=prep,
		clinical=clinical,
		technical=tech,
		geometry=geometry,
		models=models,
		encoders=encoders,
		banks=banks,
		specs=specs,
		calibrators=calibrators,
		combinations=combinations,
		rule=rule,
		readiness=readiness,
		censoring=km,
		membership=membership,
		token_names=token_names,
		mosaic=mosaic,
		constant_risk=constant,
		primary=REFERENCE_PRIMARY,
		coverage_rules=coverage_rules,
		profile_center=profile_center,
		profile_scale=profile_scale,
	)
	joblib.dump(bundle, out / "model_bundle.joblib", compress=3)
	dump(
		out / "MODEL_FROZEN.json",
		dict(
			version=VERSION,
			primary=REFERENCE_PRIMARY,
			test_Y_used_for_fitting_or_selection=False,
			test_Y_used_for_initial_stratification=not bool(a.split_file),
			independent_calibration_audit=True,
			artifact=fingerprints([out / "model_bundle.joblib"], True),
		),
	)
	testrel = test[~build]
	info.loc[testrel].reset_index(drop=True).to_csv(
		out / "test_individuals.csv", index=False
	)
	p.loc[
		test, [a.id_col, "time", "event"] + ([a.group_col] if a.group_col else [])
	].to_csv(out / "test_outcomes.csv", index=False)
	p.loc[test, [a.id_col] + [v for v in ["age", "sex", "center"] if v in p]].to_csv(
		out / "test_demographics.csv", index=False
	)
	dump(out / "prediction_columns.json", list(predictions))
	# Assay means within train-defined tokens describe profiles; they are not pathways.
	profile = (
		np.column_stack(
			[x[test][:, membership == j].mean(1) for j in range(len(token_names))]
		)
		- profile_center
	) / profile_scale
	profiles = pd.DataFrame(profile, columns=token_names)
	profiles.insert(0, a.id_col, ids[test])
	profiles.to_csv(out / "molecular_profiles.csv", index=False)

	testdetail = {k: v[testrel] for k, v in detail.items()}
	write_evidence(
		out,
		a,
		ids[test],
		x[test],
		observed[test],
		primary,
		testdetail,
		predictions[REFERENCE_PRIMARY][test],
		f_names,
		membership,
		token_names,
		calibrators[REFERENCE_PRIMARY],
		match_kw["strength"],
	)
	same_risk_pairs(
		ids[test],
		predictions["elasticnet"][test],
		profile,
		token_names,
		a.pair_caliper,
		a.max_pairs,
	).to_csv(out / "same_risk_pairs.csv", index=False)
	pd.DataFrame(
		{
			a.id_col: np.repeat(ids[test], len(token_names)),
			"token": np.tile(token_names, test.sum()),
			"reference_eid": module_reference[testrel].ravel(),
			"raw_module_risk": module_risk[testrel].ravel(),
			"mixture_weight": np.tile(mosaic.weights, test.sum()),
		}
	).to_csv(out / "test_mosaic_matches.csv.gz", index=False, compression="gzip")
	if a.explanation_samples:
		log("START", "masked_validation")
		masked_validation(
			bundle,
			x[test],
			observed[test],
			ids[test],
			None if groups is None else groups[test],
			out,
			a,
		)
		module_perturbations(
			bundle,
			x[test],
			observed[test],
			ids[test],
			None if groups is None else groups[test],
			out,
			a,
		)
	if a.attention_samples:
		log("START", "attention_audit")
		audit_attention(
			bundle,
			x[test],
			observed[test],
			ids[test],
			None if groups is None else groups[test],
			out,
			a,
		)
	log("DONE", "model_frozen", readiness["status"])
	return bundle


def predict_bundle(bundle, x, observed, p, device="cpu", details=True):
	cfg = bundle["config"]
	ids = p[cfg["id_col"]].to_numpy(str)
	groups = p[cfg["group_col"]].to_numpy(str) if cfg["group_col"] else None
	if groups is not None and p[cfg["group_col"]].isna().any():
		raise ValueError("Missing query family component")
	c = bundle["clinical"].transform(p)
	geometry = bundle["geometry"]
	pca = (
		geometry.unsupervised.transform(x)
		/ geometry.scale_u
		/ np.sqrt(len(geometry.scale_u))
	)
	blocks = {
		"clinical": c,
		"elasticnet": np.c_[c, x],
		"protein_elasticnet": x,
		"clinical_pca": np.c_[c, pca],
		"clinical_technical": np.c_[
			bundle["technical"].transform(p), 1 - observed.mean(1)
		],
	}
	raw = {
		name: obj.predict_proba(blocks.get(name, blocks["elasticnet"]))[:, 1]
		for name, obj in bundle["models"].items()
	}
	embeddings = {}
	for name, obj in bundle["encoders"].items():
		embeddings[name], head, _ = encode(obj, x, observed, device, cfg["batch_size"])
		obj.cpu()
		if name != "ssl":
			raw[name + "_direct_head"] = head
	detail = None
	for name, bank in bundle["banks"].items():
		spec = bundle["specs"][name]
		z = (
			x
			if spec["space"] == "full_proteome"
			else pca
			if spec["space"] == "pca"
			else embeddings[spec["encoder"]]
		)
		kw = {k: v for k, v in spec.items() if k not in ["encoder", "space"]}
		raw[name], result = bank.match(z, ids, groups, **kw)
		if name == REFERENCE_PRIMARY:
			detail = result
	component, ref = bundle["mosaic"].components(x, ids, groups)
	raw["mosaic_equal"], raw["mosaic_weighted"] = (
		component.mean(1),
		component @ bundle["mosaic"].weights,
	)
	predictions = {
		name: (
			bundle["calibrators"][name].predict(v)
			if name in bundle["calibrators"]
			else v
		)
		for name, v in raw.items()
	}
	for name, (source, covariate) in bundle["combinations"].items():
		predictions[name] = bundle["calibrators"][name].predict(
			raw[source], raw[covariate]
		)
	info = bundle["banks"][REFERENCE_PRIMARY].describe(
		x, observed, raw[REFERENCE_PRIMARY], detail
	)
	unknown = unknown_rows(
		p, [bundle["clinical"], getattr(bundle["prep"], "design", None)]
	)
	supported, reason = bundle["rule"].apply(
		info, 1 - observed.mean(1), unknown, cfg["sample_missing"]
	)
	info.insert(0, cfg["id_col"], ids)
	for name, v in predictions.items():
		info[name] = v
	info["supported_match"], info["rejection_reason"] = supported, reason
	info["reference_ready"] = bundle["readiness"]["status"] == "provisional_pass"
	info["prediction_released"] = supported & info.reference_ready
	info["released_net_risk"] = np.where(
		info.prediction_released, info[REFERENCE_PRIMARY], np.nan
	)
	info["raw_borrowed_risk"] = raw[REFERENCE_PRIMARY]
	info["estimated_squared_error"] = bundle["rule"].error_score(info)
	info["missing_fraction"], info["unknown_category"] = 1 - observed.mean(1), unknown
	return info, detail


def reference_project(
	run_dir,
	phe_path,
	omics_path,
	output,
	r_bin="Rscript",
	met_input="named",
	device="cpu",
):
	bundle = joblib.load(Path(run_dir) / "model_bundle.joblib")
	cfg = bundle["config"]
	idcol = cfg["id_col"]
	cols = list(
		dict.fromkeys(
			[idcol]
			+ words(cfg["covariates"])
			+ words(cfg["residualize"])
			+ ([cfg["group_col"]] if cfg["group_col"] else [])
		)
	)
	p = read_table(phe_path, idcol, cols, r_bin)
	omics = read_table(omics_path, idcol, r_bin=r_bin)
	if cfg["biom"] == "prot":
		omics.columns = [c if c == idcol else str(c).upper() for c in omics.columns]
	elif met_input == "raw":
		omics, _ = map_metabolites(omics, cfg["met_map"], idcol)
	if omics.columns.duplicated().any():
		raise ValueError("Assay names collide after normalization")
	if not set(omics[idcol]) <= set(p[idcol]):
		raise ValueError("Query IDs lack baseline metadata")
	p = p.set_index(idcol).loc[omics[idcol]].reset_index()
	# Entirely absent assays can be imputed, but their missingness affects support.
	absent = set(bundle["features"]) - set(omics)
	omics = omics.reindex(columns=[idcol] + bundle["features"])
	raw = numeric(omics, bundle["features"], "projection")
	x, observed = bundle["prep"].transform(raw, p)
	info, detail = predict_bundle(bundle, x, observed, p, device_for(device))
	output = Path(output)
	output.parent.mkdir(parents=True, exist_ok=True)
	info.to_csv(output, index=False)

	bank = bundle["banks"][REFERENCE_PRIMARY]
	write_match_table(
		output.with_name(output.stem + "_references.csv.gz"),
		cfg["id_col"],
		p[idcol].to_numpy(str),
		bank,
		detail,
	)
	dump(
		output.with_name(output.stem + "_audit.json"),
		dict(
			absent_assays=sorted(absent),
			n=len(p),
			outcome_columns_required=False,
			query_rows_attend_only_to_build_references=True,
			released=int(info.prediction_released.sum()),
		),
	)
	log("DONE", "project", f"N={len(p)}; released={info.prediction_released.sum()}")


def reference_evaluate(out):

	return evaluate_frozen(Path(out))


# Tf Model
# TabICLv2 fine-tuning: exact IPCW query loss and family-disjoint context folds.
# Adapted from the jielab/scripts TF model at commit 7be82bbc249187ff28b2d5fe4f4d61224e09b93b.
# Private numerical API pinned to tabicl==2.2.0. This actually optimizes parameters;
# sklearn fit alone only prepares in-context inference and is not fine-tuning.
#


def select_features(x, y, w, maximum):
	w = np.asarray(w, float)
	if w.sum() <= 0:
		raise ValueError("No positive build weights")
	w = w / w.sum()
	mean_y = w @ y
	mean_x = w @ x
	cov = (w * (y - mean_y)) @ x
	var = np.maximum(w @ np.square(x.astype(float)) - mean_x**2, 1e-12)
	score = np.abs(cov) / np.sqrt(var * max(mean_y * (1 - mean_y), 1e-12))
	n = x.shape[1] if maximum == 0 else min(maximum, x.shape[1])
	return np.sort(np.argsort(-score, kind="stable")[:n]), score


def context_indices(pool, y, w, size, rng):
	pool = np.asarray(pool, int)
	size = min(size, len(pool))
	cases, controls = pool[y[pool] == 1], pool[y[pool] == 0]
	if min(len(cases), len(controls)) < 2 or size < 4:
		raise ValueError("Context requires at least two known cases and controls")
	nc = int(round(size * w[cases].sum() / w[pool].sum()))
	nc = min(len(cases), size - 2, max(2, size - len(controls), nc))
	result = []
	for part, count in [(cases, nc), (controls, size - nc)]:
		result.extend(rng.choice(part, count, replace=False, p=w[part] / w[part].sum()))
	return rng.permutation(result)


def fold_labels(y, groups, count, seed):
	cv = (
		StratifiedKFold(count, shuffle=True, random_state=seed)
		if groups is None
		else StratifiedGroupKFold(count, shuffle=True, random_state=seed)
	)
	folds = np.empty(len(y), int)
	for f, (_, ix) in enumerate(cv.split(np.zeros(len(y)), y, groups)):
		folds[ix] = f
	if (
		groups is not None
		and pd.DataFrame({"group": groups, "fold": folds})
		.groupby("group")
		.fold.nunique()
		.max()
		> 1
	):
		raise ValueError("Family leakage in context folds")
	return folds


def _predictor(model, config, cx, cy, a):
	from tabicl import TabICLClassifier

	estimator = TabICLClassifier(
		n_estimators=1,
		norm_methods=["none"],
		feat_shuffle_method="none",
		class_shuffle_method="none",
		softmax_temperature=1.0,
		device=a.device,
		use_amp=a.amp and a.device == "cuda",
		use_fa3=False,
		n_jobs=a.cores,
		kv_cache=True,
		allow_auto_download=False,
		random_state=a.seed,
	)
	estimator.model_, estimator.model_config_, estimator.model_path_ = (
		model,
		config,
		a.model_path,
	)
	estimator._load_model = lambda: None
	model.eval()
	return estimator.fit(cx, cy)


def predict(estimator, x, batch_size):
	rows = []
	for start in range(0, len(x), batch_size):
		pred = estimator.predict_proba(x[start : start + batch_size])[:, 1]
		if not np.isfinite(pred).all():
			raise FloatingPointError("Nonfinite Transformer predictions")
		rows.append(pred)
	return np.concatenate(rows) if rows else np.empty(0)


def train_transformer(x, y, w, ids, groups, xv, yv, wv, a, out):
	from tabicl._model import TabICL
	from tabicl._sklearn.preprocessing import EnsembleGenerator

	torch.manual_seed(a.seed)
	if a.device == "cuda":
		torch.cuda.manual_seed_all(a.seed)
	checkpoint = torch.load(a.model_path, map_location="cpu", weights_only=True)
	config = checkpoint["config"]
	model = TabICL(**config)
	model.load_state_dict(checkpoint["state_dict"])
	if not a.unfreeze_encoder:
		for m in [model.col_embedder, model.row_interactor]:
			m.requires_grad_(False)
	model.to(a.device)
	folds = fold_labels(y, groups, a.folds, a.seed + 41)
	pd.DataFrame({"eid": ids, "context_exclusion_fold": folds}).to_csv(
		out / "inner_folds.csv", index=False
	)
	context = context_indices(
		np.arange(len(y)), y, w, a.context_size, np.random.default_rng(a.seed + 42)
	)
	pd.DataFrame({"eid": ids[context], "Y": y[context], "IPCW": w[context]}).to_csv(
		out / "context.csv", index=False
	)
	optimizer = torch.optim.AdamW(
		[p for p in model.parameters() if p.requires_grad],
		lr=a.learning_rate,
		weight_decay=a.weight_decay,
	)
	scaler = torch.amp.GradScaler("cuda", enabled=a.amp and a.device == "cuda")
	history = []
	best = np.inf
	bad = 0
	selected_epoch = -1
	best_file = abm_cache_dir(out) / "selected_weights.ckpt"
	for epoch in range(a.epochs + 1):
		start_time = time.monotonic()
		numerator = denominator = 0.0
		steps = 0
		if epoch:
			rng = np.random.default_rng(a.seed + epoch)
			model.train()
			for fold in np.unique(folds):
				pool = np.flatnonzero(folds != fold)
				ctx = context_indices(pool, y, w, a.context_size, rng)
				query = rng.permutation(np.flatnonzero(folds == fold))
				gen = EnsembleGenerator(
					classification=True,
					n_estimators=1,
					norm_methods=["none"],
					feat_shuffle_method="none",
					class_shuffle_method="none",
					random_state=a.seed,
				)
				gen.fit(x[ctx], y[ctx])
				for begin in range(0, len(query), a.batch_size):
					ix = query[begin : begin + a.batch_size]
					xx, yy = next(iter(gen.transform(x[ix], mode="both").values()))
					xx = torch.as_tensor(xx, dtype=torch.float32, device=a.device)
					yy = torch.as_tensor(yy, dtype=torch.float32, device=a.device)
					target = torch.as_tensor(y[ix], dtype=torch.long, device=a.device)
					weight = torch.as_tensor(
						w[ix], dtype=torch.float32, device=a.device
					)
					optimizer.zero_grad(set_to_none=True)
					with torch.autocast(
						device_type=a.device,
						dtype=torch.float16,
						enabled=a.amp and a.device == "cuda",
					):
						logits = model(xx, yy)[0, :, :2]
						value = (
							F.cross_entropy(logits.float(), target, reduction="none")
							* weight
						).sum() / weight.sum()
					if not torch.isfinite(value):
						raise FloatingPointError("Nonfinite IPCW training loss")
					scaler.scale(value).backward()
					scaler.unscale_(optimizer)
					torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
					scaler.step(optimizer)
					scaler.update()
					numerator += value.item() * float(weight.sum())
					denominator += float(weight.sum())
					steps += 1
		estimator = _predictor(model, config, x[context], y[context], a)
		valid = predict(estimator, xv, a.predict_batch_size)
		score = float(np.average(loss(yv, valid), weights=wv))
		if not np.isfinite(score):
			raise FloatingPointError("Nonfinite TF validation loss")
		improved = score < best - 1e-6
		if improved:
			best, bad, selected_epoch = score, 0, epoch
			state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
			temp = best_file.with_suffix(".tmp")
			torch.save(dict(config=config, state_dict=state), temp)
			temp.replace(best_file)
		else:
			bad += 1
		history.append(
			dict(
				epoch=epoch,
				steps=steps,
				train_IPCW_logloss=numerator / denominator if denominator else None,
				tune_IPCW_logloss=score,
				selected=improved,
				seconds=time.monotonic() - start_time,
			)
		)
		pd.DataFrame(history).to_csv(out / "learning_curve.csv", index=False)
		log(
			"EPOCH",
			"TF",
			f"{epoch}/{a.epochs}; tune={score:.5f}; best_epoch={selected_epoch}",
		)
		del estimator
		model.clear_cache()
		if epoch and bad >= a.patience:
			break
	selected = torch.load(best_file, map_location="cpu", weights_only=True)
	model.load_state_dict(selected["state_dict"])
	changes = sum(
		not torch.equal(v, checkpoint["state_dict"][k])
		for k, v in selected["state_dict"].items()
	)
	dump(
		out / "training_summary.json",
		dict(
			gradient_steps=sum(r["steps"] for r in history),
			best_epoch=selected_epoch,
			selected_tensors_changed=changes,
			loss="IPCW weighted cross entropy",
			context_rows=len(context),
			family_context_exclusion=groups is not None,
			pretrained_baseline_epoch=0,
			context_weighting="weighted sampling without replacement; approximate, not exact IPCW attention",
		),
	)
	return _predictor(model, config, x[context], y[context], a), selected, context


def restore_predictor(bundle, a):
	from tabicl._model import TabICL

	cp = bundle["transformer"]
	model = TabICL(**cp["config"])
	model.load_state_dict(cp["state_dict"])
	model.to(a.device)
	return _predictor(model, cp["config"], bundle["context_x"], bundle["context_y"], a)


# Tf Pipeline
# Supervised numerical Transformer workflow; test outcomes stay sealed during fitting.
# Adapted from the currently inspected abm_TF protocol. Inference is not literal COPY1.
#


TABICL_VERSION = "TF-LE8-2026-10-02"
TABICL_PRIMARY = "tabicl_finetuned"


def manifest(a, out, model_info):
	paths = [] if a.demo else [a.phe_file, a.omics_file]
	paths += [
		v
		for v in [
			a.split_file,
			(
				a.met_map
				if a.biom == "met" and a.met_input == "raw" and not a.demo
				else ""
			),
		]
		if v
	]
	value = dict(
		version=TABICL_VERSION,
		config=vars(a),
		pretrained=model_info,
		dependencies={
			name: importlib.metadata.version(name)
			for name in ["torch", "tabicl", "scikit-learn", "numpy", "pandas"]
		},
		inputs=fingerprints(paths, a.full_input_hash),
		code=fingerprints([Path(__file__)], True),
	)
	value["signature"] = digest(value)
	dump(out / "manifest.json", value)


def cohort(a, out):
	folder = abm_cache_dir(out, "input")
	folder.mkdir(exist_ok=True)
	p = prepare(a, folder)
	raw = np.load(folder / "raw.npy")
	features = (folder / "features.txt").read_text().splitlines()
	p, audit = outcomes(p, a)
	eligible = p.eligible.to_numpy(bool)
	raw, p = raw[eligible], p.loc[eligible].reset_index(drop=True)
	p["split"] = split_people(p, a)
	prep = MolecularPreprocessor(
		a.feature_missing, words(a.residualize), words(a.categorical), a.transform
	)
	build = p.split.eq("build").to_numpy()
	scaled = prep.scale_transform(raw[build])
	keep = (np.mean(~np.isfinite(scaled), axis=0) < a.feature_missing) & (
		np.nanstd(scaled, axis=0) > 1e-8
	)
	if keep.sum() < 3:
		raise ValueError("Fewer than three assays pass build QC")
	missing = np.mean(~np.isfinite(prep.scale_transform(raw[:, keep])), axis=1)
	qc = missing <= a.sample_missing
	pd.DataFrame(
		{
			a.id_col: p[a.id_col],
			"split": p.split,
			"missing_fraction": missing,
			"included": qc,
		}
	).to_csv(out / "sample_qc.csv", index=False)
	raw, p = raw[qc], p.loc[qc].reset_index(drop=True)
	masks = {
		name: p.split.eq(name).to_numpy()
		for name in ["build", "tune", "calibration", "test"]
	}
	if min(m.sum() for m in masks.values()) < 30:
		raise ValueError("Each partition requires at least 30 participants after QC")
	prep.fit(raw[masks["build"]], p.loc[masks["build"]], allowed=keep)
	x, observed = prep.transform(raw, p)
	clinical = MetadataDesign(words(a.covariates), words(a.categorical)).fit(
		p.loc[masks["build"]]
	)
	c = clinical.transform(p)
	if c.shape[1] == 0:
		raise ValueError("Clinical covariates are needed for the comparator")
	km = CensoringKM().fit(p.loc[masks["build"]], a.horizon, a.min_censor_survival)
	y, w = np.zeros(len(p), int), np.zeros(len(p))
	dev = ~masks["test"]
	y[dev], w[dev] = km.labels_weights(p.loc[dev])
	for name in ["build", "tune", "calibration"]:
		known = masks[name] & (w > 0)
		if min(np.sum(y[known] == label) for label in [0, 1]) < a.min_events:
			raise ValueError(f"Too few known outcomes in {name}")
	p[[a.id_col, "split"]].to_csv(out / "split.csv", index=False)
	audit.update(
		after_qc=len(p),
		retained_features=len(prep.keep),
		role_counts=p.split.value_counts().to_dict(),
		family_split=bool(a.group_col),
		horizon=a.horizon,
	)
	dump(out / "cohort_audit.json", audit)
	return p, x, c, observed, y, w, masks, prep, clinical, km, features


def decision_threshold(y, w, risk):
	valid = w > 0
	fpr, tpr, thresholds = roc_curve(y[valid], risk[valid], sample_weight=w[valid])
	ok = np.isfinite(thresholds) & (thresholds >= 0) & (thresholds <= 1)
	if not ok.any():
		return 0.5
	return float(thresholds[ok][np.argmax((tpr - fpr)[ok])])


def tabicl_train(a, out, model_info):
	manifest(a, out, model_info)
	p, x, c, observed, y, w, masks, prep, clinical, km, features = cohort(a, out)
	build, tune, cal, test = [
		masks[s] for s in ["build", "tune", "calibration", "test"]
	]
	selected, importance = select_features(x[build], y[build], w[build], a.max_features)
	names = np.asarray(features)[prep.keep]
	pd.DataFrame(
		{
			"feature": names,
			"build_weighted_association": importance,
			"selected": np.isin(np.arange(len(names)), selected),
		}
	).to_csv(out / "feature_selection.csv", index=False)
	z = x[:, selected]
	if a.with_clinical:
		z = np.c_[z, c]
	known = build & (w > 0)
	val = np.flatnonzero(tune & (w > 0))
	if min(np.bincount(y[known], minlength=2)) < a.folds:
		raise ValueError("Insufficient known outcomes for context folds")
	if a.validation_samples and len(val) > a.validation_samples:
		from sklearn.model_selection import train_test_split

		val, _ = train_test_split(
			val,
			train_size=a.validation_samples,
			stratify=y[val],
			random_state=a.seed + 53,
		)
	ids = p[a.id_col].to_numpy(str)
	groups = p[a.group_col].to_numpy(str) if a.group_col else None
	p.iloc[val][[a.id_col]].to_csv(out / "early_stop_ids.csv", index=False)
	estimator, checkpoint, ctx = train_transformer(
		z[known],
		y[known],
		w[known],
		ids[known],
		None if groups is None else groups[known],
		z[val],
		y[val],
		w[val],
		a,
		out,
	)
	holdout = ~build
	raw = {TABICL_PRIMARY: predict(estimator, z[holdout], a.predict_batch_size)}
	estimator.model_.cpu()
	del estimator
	models = {}
	tuning = []
	for name, block in [("clinical", c), ("elasticnet", z)]:
		obj, rows = tune_logistic(
			block, y, w, build, tune, name, a, 0.5 if name == "elasticnet" else 0
		)
		models[name] = obj
		tuning.extend(rows)
		raw[name] = obj.predict_proba(block[holdout])[:, 1]
	pd.DataFrame(tuning).to_csv(out / "comparator_tuning.csv", index=False)
	raw["constant"] = np.full(holdout.sum(), np.average(y[build], weights=w[build]))
	calrel, tunerel, testrel = cal[holdout], tune[holdout], test[holdout]
	calibrators = {}
	thresholds = {}
	probabilities = {}
	for name, risk in raw.items():
		calibrator = RiskCalibrator().fit(risk[calrel], y[cal], w[cal])
		calibrators[name] = calibrator
		probabilities[name] = calibrator.predict(risk)
		th = decision_threshold(y[tune], w[tune], risk[tunerel])
		thresholds[name] = float(calibrator.predict(np.array([th]))[0])
	info = pd.DataFrame(
		{a.id_col: ids[test], "missing_fraction": 1 - observed[test].mean(1)}
	)
	for name in raw:
		info[name + "_raw"] = raw[name][testrel]
		info[name] = probabilities[name][testrel]
		info[name + "_predicted_Y"] = (info[name] >= thresholds[name]).astype(int)
	bundle = dict(
		version=TABICL_VERSION,
		config=vars(a),
		primary=TABICL_PRIMARY,
		prep=prep,
		clinical=clinical,
		censoring=km,
		features=features,
		selected=selected,
		transformer=checkpoint,
		context_x=z[known][ctx],
		context_y=y[known][ctx],
		context_ids=ids[known][ctx],
		context_groups=None if groups is None else groups[known][ctx],
		models=models,
		calibrators=calibrators,
		thresholds=thresholds,
		constant=float(raw["constant"][0]),
		pretrained=model_info,
	)
	joblib.dump(bundle, out / "model_bundle.joblib", compress=3)
	info.to_csv(out / "test_predictions.csv", index=False)
	p.loc[
		test, [a.id_col, "time", "event"] + ([a.group_col] if a.group_col else [])
	].to_csv(out / "test_outcomes.csv", index=False)
	dump(
		out / "thresholds.json",
		dict(
			method="Tune-only Youden J mapped through independent calibration fit",
			thresholds=thresholds,
		),
	)
	dump(
		out / "MODEL_FROZEN.json",
		dict(
			version=TABICL_VERSION,
			test_outcomes_used_for_training=False,
			artifact=fingerprints([out / "model_bundle.joblib"], True),
		),
	)
	dump(out / "TRAIN_DONE.json", dict(version=TABICL_VERSION))


def tabicl_evaluate(out):
	if not (out / "MODEL_FROZEN.json").is_file():
		raise ValueError("No frozen model")
	b = joblib.load(out / "model_bundle.joblib")
	a = SimpleNamespace(**b["config"])
	pred = pd.read_csv(out / "test_predictions.csv", dtype={a.id_col: str})
	truth = pd.read_csv(out / "test_outcomes.csv", dtype={a.id_col: str})
	if pred[a.id_col].duplicated().any() or truth[a.id_col].duplicated().any():
		raise ValueError("Duplicate evaluation IDs")
	if (
		pred[a.id_col].isna().any()
		or truth[a.id_col].isna().any()
		or set(pred[a.id_col]) != set(truth[a.id_col])
	):
		raise ValueError(
			"Evaluation predictions and outcomes must have exactly the same nonmissing IDs"
		)
	d = truth.merge(pred, on=a.id_col, validate="one_to_one")
	y, w = b["censoring"].labels_weights(d)
	rows = []
	cal = []
	for name in b["calibrators"]:
		rows.append(
			dict(
				model=name,
				N=len(d),
				known_N=int((w > 0).sum()),
				horizon=a.horizon,
				**metrics(y, w, d[name].to_numpy()),
			)
		)
		# Equal-sized risk groups, including all test rows; descriptive IPCW calibration.
		g = pd.qcut(
			pd.Series(d[name]).rank(method="first"), q=min(10, len(d)), labels=False
		)
		for group in sorted(g.unique()):
			ix = g.eq(group).to_numpy()
			wt = w[ix]
			cal.append(
				dict(
					model=name,
					group=int(group) + 1,
					N=int(ix.sum()),
					predicted=float(d.loc[ix, name].mean()),
					observed_IPCW=(
						float(np.average(y[ix], weights=wt)) if wt.sum() > 0 else np.nan
					),
				)
			)
	result = pd.DataFrame(rows)
	result.to_csv(out / "test_metrics.csv", index=False)
	pd.DataFrame(cal).to_csv(out / "test_calibration.csv", index=False)
	text = [
		"# ABM TF / LE8",
		"",
		f"Outcome {a.trait}; layer {a.biom}; horizon {a.horizon:g} years.",
		"",
		result.to_markdown(index=False),
		"",
		"All preprocessing, selection, context, weights and calibration were frozen before test scoring.",
		"TabICLv2 predicts from labeled context and a learned head; it does not literally copy a donor outcome.",
		"Death is censored: net risk under independent-censoring assumptions, not a competing-risk cumulative incidence.",
		"Each method is shown. Test ranking must not be used to choose a new winning model.",
	]
	(out / "REPORT.md").write_text("\n".join(text) + "\n")
	dump(out / "DONE.json", dict(version=TABICL_VERSION))
	return result


def tabicl_project(a, out):
	b = joblib.load(out / "model_bundle.joblib")
	cfg = b["config"]
	idcol = cfg["id_col"]
	cols = list(
		dict.fromkeys(
			[idcol]
			+ words(cfg["covariates"])
			+ words(cfg["residualize"])
			+ ([cfg["group_col"]] if cfg["group_col"] else [])
		)
	)
	p = read_table(a.phe_file, idcol, cols, a.r_bin)
	m = read_table(a.omics_file, idcol, r_bin=a.r_bin)
	if cfg["biom"] == "prot":
		m.columns = [c if c == idcol else str(c).upper() for c in m]
	elif a.met_input == "raw":
		m, _ = map_metabolites(m, cfg["met_map"], idcol)
	if not set(m[idcol]) <= set(p[idcol]):
		raise ValueError("Missing query metadata")
	p = p.set_index(idcol).loc[m[idcol]].reset_index()
	ids = p[idcol].to_numpy(str)
	if set(ids) & set(b["context_ids"]):
		raise ValueError(
			"Query contains labeled context participants; use a truly held-out query"
		)
	if cfg["group_col"] and b["context_groups"] is not None:
		if set(p[cfg["group_col"]].astype(str)) & set(b["context_groups"]):
			raise ValueError("Query families overlap labeled context")
	m = m.reindex(columns=[idcol] + b["features"])
	raw = numeric(m, b["features"], "TF projection")
	x, observed = b["prep"].transform(raw, p)
	z = x[:, b["selected"]]
	c = b["clinical"].transform(p)
	if cfg["with_clinical"]:
		z = np.c_[z, c]
	runtime = SimpleNamespace(
		**{
			**cfg,
			"device": a.device,
			"cores": a.cores,
			"amp": a.amp,
			"model_path": a.model_path,
		}
	)
	estimator = restore_predictor(b, runtime)
	r = predict(estimator, z, a.predict_batch_size)
	values = b["calibrators"][TABICL_PRIMARY].predict(r)
	missing = 1 - observed.mean(1)
	unknown = unknown_rows(p, [b["clinical"], getattr(b["prep"], "design", None)])
	release = (missing <= cfg["sample_missing"]) & ~unknown
	result = pd.DataFrame(
		{
			idcol: ids,
			"net_risk": values,
			"predicted_Y": (values >= b["thresholds"][TABICL_PRIMARY]).astype(int),
			"missing_fraction": missing,
			"unknown_category": unknown,
			"prediction_released": release,
			"released_net_risk": np.where(release, values, np.nan),
		}
	)
	Path(a.output).parent.mkdir(parents=True, exist_ok=True)
	result.to_csv(a.output, index=False)
	return result


# Lossless RDS column-selection bridge, written only to a temporary directory.
RDS_EXPORT_SCRIPT = '# Read-only bridge. No sourcing of phenotype pipelines or installing R packages.\nargs <- commandArgs(trailingOnly = TRUE)\nstopifnot(length(args) == 3L)\nx <- readRDS(args[1])\nif (!is.data.frame(x)) stop("Expected one data.frame in RDS")\ncols <- readLines(args[3], warn = FALSE)\nif (length(cols)) {\n  absent <- setdiff(cols, names(x))\n  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))\n  x <- x[, cols, drop = FALSE]\n}\nfor (nm in names(x)) {\n  if (inherits(x[[nm]], "Date") || inherits(x[[nm]], "POSIXt"))\n    x[[nm]] <- format(x[[nm]], "%Y-%m-%d")\n  if (inherits(x[[nm]], "integer64")) x[[nm]] <- as.character(x[[nm]])\n}\n# Binary transfer preserves subnormal numeric metadata used as category codes.\nsaveRDS(x, args[2], compress = FALSE)\n\n'


# 🚩 Canonical model serialization
# A single module identity is used for newly fitted models, including CLI runs.
sys.modules["c1_abm"] = sys.modules[__name__]
for _class_name in (
	"ContextBlock",
	"ContextStack",
	"ReferenceBank",
	"SupportRule",
	"RiskCalibrator",
	"ModuleMosaic",
	"RowEncoder",
	"MetadataDesign",
	"MolecularPreprocessor",
	"MolecularGeometry",
	"CensoringKM",
):
	globals()[_class_name].__module__ = "c1_abm"


if __name__ == "__main__":
	try:
		reference_main() if _WORKER == "reference" else tabicl_main()
	except Exception as exc:
		print(f"[ABM] ERROR {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
		raise
