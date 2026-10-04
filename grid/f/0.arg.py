#!/usr/bin/env python3
"""0.arg.py: build-asmc-memory, build-argneedle-map, validate-haps, make-arg-keep, run-argneedle-advanced, argn-to-trees, make-sample-map, make-anchors, arg-features, arg-affinity. Use a subcommand followed by --help."""

from __future__ import annotations


# 🚩 argneedle-memory: argneedle_memory
"""Memory-bounded adapter for ARG-Needle's hashed threading path.

Keep the installed package untouched. Transform the dense-posterior operations
and native batch size, refusing unknown upstream code rather than falling back.
For each site, smoothing reads only the posterior of the selected cousin.
Across a smooth interval the selected cousin is constant, so storing those
selected values gives the same float64 means without a samples-by-sites matrix.
Hashing, ASMC pair requests, random draws, MAP boundaries and threading stay
upstream. Native ASMC workspaces use batches of 16 instead of the default 64.
"""
import inspect
import logging
import textwrap
import importlib.util
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import sysconfig


def backend_path():
	return (
		Path.home()
		/ ".cache/gu/asmc-memory/v1.4.0-p2"
		/ sys.implementation.cache_tag
		/ ("asmc_python_bindings" + sysconfig.get_config_var("EXT_SUFFIX"))
	)


def _load_backend():
	name = "asmc.asmc.asmc_python_bindings"
	if name in sys.modules:
		if getattr(sys.modules[name], "_gu_bounded_workspace", 0) == 2:
			return
		raise RuntimeError("Load argneedle_memory before importing arg_needle or asmc.asmc")
	path = backend_path()
	if not path.is_file():
		raise RuntimeError(
			"Missing bounded ASMC backend. Run grid Python on "
			f"{Path(__file__).with_name('0.arg.py')} "
			"--source /path/to/ASMC-1.4.0"
		)
	spec = importlib.util.spec_from_file_location(name, path)
	module = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(module)
	if getattr(module, "_gu_bounded_workspace", 0) != 2:
		raise RuntimeError("ASMC backend does not contain the workspace fix")
	sys.modules[name] = module


def _rewrite(function, replacements):
	source = textwrap.dedent(inspect.getsource(function))
	for old, new in replacements:
		if source.count(old) != 1:
			raise RuntimeError(
				f"Unsupported ARG-Needle implementation in {function.__name__}: expected exactly one {old!r}"
			)
		source = source.replace(old, new, 1)
	namespace = {}
	exec(compile(source, inspect.getfile(function) + "[bounded-memory]", "exec"), function.__globals__, namespace)
	return namespace[function.__name__]


def install():
	_load_backend()
	import arg_needle.inference as inference
	import arg_needle.decoders as decoders
	from arg_needle.decoders import ASMCDecoder

	if getattr(inference.thread_samples, "_bounded_memory", False):
		return
	thread = _rewrite(
		inference.thread_samples,
		[
			(
				"np.full((start_thread_id + num_next_samples - 1, len(posterior_phys_pos)), np.nan)",
				"np.full(len(posterior_phys_pos), np.nan)",
			),
			(
				"np.mean(tmrca_mean[indices[begin], begin:end])",
				"np.mean(tmrca_mean[begin:end] if hash_topk > 0 else tmrca_mean[indices[begin], begin:end])",
			),
		],
	)
	decode = _rewrite(
		ASMCDecoder.compute_with_hashing,
		[
			(
				"tmrca_mean[other_ids, from_pos:to_pos] = batch_mean\n            foo = np.argmin(batch_mean, axis=0)",
				"foo = np.argmin(batch_mean, axis=0)\n"
				"            tmrca_mean[from_pos:to_pos] = "
				"batch_mean[foo, np.arange(to_pos - from_pos)]",
			),
		],
	)
	make_decoder = _rewrite(
		decoders.make_asmc_decoder,
		[
			# Keep all 64 candidates. Only the internal SIMD batch/workspace shrinks.
			("asmc_obj = ASMC(params)", "params.batchSize = 16\n    asmc_obj = ASMC(params)"),
		],
	)
	# Commit together only after all upstream functions pass compatibility checks.
	thread._bounded_memory = True
	inference.thread_samples = thread
	ASMCDecoder.compute_with_hashing = decode
	decoders.make_asmc_decoder = make_decoder
	inference.make_asmc_decoder = make_decoder
	logging.info(
		"Bounded-memory threading enabled: one float64 posterior per site; "
		"ASMC batches <=16; native scratch pages discarded after each window"
	)


# 🚩 build-asmc-memory: build_asmc_memory
"""Build an isolated ASMC 1.4.0 backend with discardable HMM workspace pages.

Usage: grid-python 0.arg.py --source /path/to/ASMC-1.4.0
Requires CMake, C++17, Boost and zlib. CMake fetches the pinned upstream deps;
repeat --cmake-arg=... to supply local FETCHCONTENT_SOURCE_DIR_* overrides.
The installed ASMC package is never modified.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys


PATCH = r"""
// GU: results have been copied out before this is called. These two workspaces
// are scratch, fully rewritten on the next forward/backward pass. Discard only
// whole pages strictly inside each Eigen allocation (never allocator metadata).
void HMM::releaseWorkspacePages()
{
  const auto page = static_cast<std::uintptr_t>(::sysconf(_SC_PAGESIZE));
  if (page == 0 || page > 1024 * 1024) {
    throw std::runtime_error("Cannot determine workspace page size");
  }
  for (auto* buffer : {&m_alphaBuffer, &m_betaBuffer}) {
    const auto address = reinterpret_cast<std::uintptr_t>(buffer->data());
    const auto begin = ((address + page - 1) / page) * page;
    const auto end = ((address + buffer->size() * sizeof(float)) / page) * page;
    if (end > begin && ::madvise(reinterpret_cast<void*>(begin), end - begin, MADV_DONTNEED) != 0) {
      throw std::runtime_error("Failed to discard ASMC scratch workspace pages");
    }
  }
}
"""


def build_asmc_memory_main():
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("--source", type=Path, required=True)
	p.add_argument("--jobs", type=int, default=2)
	p.add_argument("--cmake-arg", action="append", default=[])
	a = p.parse_args()
	if sys.platform != "linux":
		p.error("This backend requires Linux madvise (WSL is supported)")
	if "project(asmc LANGUAGES CXX VERSION 1.4.0)" not in (a.source / "CMakeLists.txt").read_text():
		p.error("Expected upstream ASMC v1.4.0 source")
	target = backend_path()
	work = target.parent / "build-source"
	work.mkdir(parents=True, exist_ok=True)
	shutil.copytree(a.source, work, dirs_exist_ok=True)

	def change(name, old, new):
		path = work / "src" / name
		source = path.read_text()
		if source.count(old) != 1:
			raise RuntimeError(f"Unexpected upstream source: {name}: {old!r}")
		path.write_text(source.replace(old, new, 1))

	change(
		"HMM.hpp",
		"DecodePairsReturnStruct& getDecodePairsReturnStruct();",
		"DecodePairsReturnStruct& getDecodePairsReturnStruct();\n  void releaseWorkspacePages();",
	)
	change(
		"HMM.cpp",
		'#include "HMM.hpp"',
		'#include "HMM.hpp"\n#include <cstdint>\n#include <sys/mman.h>\n#include <unistd.h>',
	)
	change("HMM.cpp", "void HMM::finishDecoding()", PATCH + "\nvoid HMM::finishDecoding()")
	change(
		"ASMC.cpp",
		"mHmm.getDecodePairsReturnStruct().finaliseCalculations();",
		"mHmm.getDecodePairsReturnStruct().finaliseCalculations();\n  mHmm.releaseWorkspacePages();",
	)
	# chr1: 2,054,102 sites * 69 states * 16 lanes exceeds signed int32.
	simd_path = work / "src/Simd.cpp"
	simd = simd_path.read_text()
	for old, new, count in [
		("for (int pos = from; pos < to; ++pos)", "for (Eigen::Index pos = from; pos < to; ++pos)", 4),
		("const int ind = ", "const Eigen::Index ind = ", 3),
	]:
		if simd.count(old) != count:
			raise RuntimeError("Unexpected SIMD index implementation")
		simd = simd.replace(old, new)
	simd_path.write_text(simd)
	change(
		"HMM.cpp",
		"m_alphaBuffer.resize(sequenceLength * states * m_batchSize);",
		"m_alphaBuffer.resize(static_cast<Eigen::Index>(sequenceLength) * states * m_batchSize);",
	)
	change(
		"HMM.cpp",
		"m_betaBuffer.resize(sequenceLength * states * m_batchSize);",
		"m_betaBuffer.resize(static_cast<Eigen::Index>(sequenceLength) * states * m_batchSize);",
	)
	change(
		"pybind.cpp",
		"PYBIND11_MODULE(asmc_python_bindings, m)\n{",
		'PYBIND11_MODULE(asmc_python_bindings, m)\n{\n  m.attr("_gu_bounded_workspace") = 2;',
	)
	build = target.parent / "build"
	subprocess.run(
		[
			"cmake",
			"-S",
			str(work),
			"-B",
			str(build),
			"-DASMC_PYTHON_BINDINGS=ON",
			"-DASMC_TESTING=OFF",
			"-DCMAKE_BUILD_TYPE=Release",
			"-DPYTHON_EXECUTABLE=" + sys.executable,
			"-DPython_EXECUTABLE=" + sys.executable,
			"-DCMAKE_PREFIX_PATH=" + sys.prefix,
			*a.cmake_arg,
		],
		check=True,
	)
	subprocess.run(["cmake", "--build", str(build), "--target", "asmc_python_bindings", "-j", str(a.jobs)], check=True)
	(artifact,) = build.glob("asmc_python_bindings*.so")
	shutil.copy2(artifact, target.with_suffix(".next"))
	target.with_suffix(".next").replace(target)
	(target.parent / "build.json").write_text(
		json.dumps(
			{
				"upstream": "https://github.com/PalamaraLab/ASMC/tree/v1.4.0",
				"patch": 2,
				"python": sys.version,
				"artifact": str(target),
			},
			indent=2,
		)
		+ "\n"
	)
	print(target)


def build_asmc_memory_cli():
	build_asmc_memory_main()


# 🚩 build-argneedle-map: build_argneedle_map
"""Interpolate a GRCh37 cumulative genetic map onto every Oxford HAPS variant."""
import argparse, gzip, re
from pathlib import Path
import numpy as np


def map_open(path):
	return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def numeric(x):
	try:
		return float(x)
	except:
		return None


def read_source(path):
	rows = []
	header = None
	with map_open(path) as h:
		for line in h:
			if not line.strip() or line.lstrip().startswith("#"):
				continue
			z = re.split(r"\s+", line.strip())
			if header is None and any(
				x.lower() in {"position", "pos", "bp", "base_pair_location"} or "position" in x.lower() for x in z
			):
				header = [x.lower() for x in z]
				continue
			rows.append(z)
	if not rows:
		raise SystemExit(f"No map rows in {path}")
	if header:

		def idx(pred):
			return next((i for i, x in enumerate(header) if pred(x)), None)

		ip = idx(lambda x: x in {"position", "pos", "bp", "base_pair_location"} or "position" in x)
		ic = idx(
			lambda x: ("cm" in x and "rate" not in x) or x in {"map", "genetic_map"} or ("map" in x and "rate" not in x)
		)
		ir = idx(lambda x: "rate" in x)
	else:
		ip = ic = ir = None
	vals = []
	for z in rows:
		nums = [numeric(x) for x in z]
		if ip is not None and ip < len(nums) and nums[ip] is not None:
			bp = nums[ip]
		elif len(nums) >= 4 and nums[-1] is not None:
			bp = nums[-1]
		else:
			bp = next((x for x in nums if x is not None and x >= 1), None)
		if ic is not None and ic < len(nums):
			cm = nums[ic]
		elif len(nums) >= 4:
			cm = nums[-2]
		elif len(nums) >= 3:
			cm = nums[-1]
		else:
			cm = None
		rate = nums[ir] if ir is not None and ir < len(nums) else None
		if bp is not None:
			vals.append((float(bp), None if cm is None else float(cm), rate))
	vals = sorted({int(bp): (bp, cm, rate) for bp, cm, rate in vals}.values())
	bp = np.array([x[0] for x in vals], float)
	cm = np.array([np.nan if x[1] is None else x[1] for x in vals], float)
	if np.isfinite(cm).sum() < 2:
		rate = np.array([np.nan if x[2] is None else x[2] for x in vals], float)
		if np.isfinite(rate).sum() < 2:
			raise SystemExit("Map needs cumulative cM or recombination rate")
		rate = np.interp(bp, bp[np.isfinite(rate)], rate[np.isfinite(rate)])
		cm = np.zeros(len(bp))
		cm[1:] = np.cumsum((bp[1:] - bp[:-1]) * (rate[1:] + rate[:-1]) / 2 / 1e6)
	else:
		cm = np.interp(bp, bp[np.isfinite(cm)], cm[np.isfinite(cm)])
	if np.any(np.diff(bp) < 0) or np.any(np.diff(cm) < -1e-6):
		raise SystemExit("Source map is not monotonic")
	return bp, cm


def haps_variants(path):
	with map_open(path) as h:
		for line in h:
			if not line.strip():
				continue
			z = line.split()
			# Oxford/PLINK HAPS: chr ID BP A0 A1 ...; accept chr ID rsID BP A0 A1 ...
			if len(z) < 6:
				raise SystemExit("Malformed HAPS")
			if numeric(z[2]) is not None:
				chrom, sid, pos = z[0], z[1], int(float(z[2]))
			elif numeric(z[3]) is not None:
				chrom, sid, pos = z[0], z[2], int(float(z[3]))
			else:
				raise SystemExit(f"Cannot identify HAPS position: {' '.join(z[:6])}")
			yield chrom, sid, pos


def build_argneedle_map_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--source", required=True)
	ap.add_argument("--haps", required=True)
	ap.add_argument("--out", required=True)
	ap.add_argument("--chr", required=True)
	a = ap.parse_args()
	bp, cm = read_source(a.source)
	Path(a.out).parent.mkdir(parents=True, exist_ok=True)
	n = 0
	prev = -1
	prev_g = -float("inf")
	adjusted = 0
	with open(a.out, "w") as o:
		for chrom, sid, pos in haps_variants(a.haps):
			if pos <= prev:
				raise SystemExit("HAPS positions must be unique and strictly increasing")
			prev = pos
			g = round(float(np.interp(pos, bp, cm, left=cm[0], right=cm[-1])), 9)
			# Map plateaus are common. A 1e-9 cM numerical tie-break satisfies ASMC.
			if g <= prev_g:
				g = round(prev_g + 1e-9, 9)
				adjusted += 1
			prev_g = g
			o.write(f"{a.chr}\t{sid}\t{g:.9f}\t{pos}\n")
			n += 1
	if n < 2:
		raise SystemExit("Too few mapped variants")
	print(f"variants={n} source_points={len(bp)} plateau_tiebreaks={adjusted} map={a.out}")


def build_argneedle_map_cli():
	build_argneedle_map_main()


# 🚩 validate-haps: validate_haps
import argparse, gzip, json
from pathlib import Path


def haps_open(path):
	return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def validate_haps_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--haps", required=True)
	ap.add_argument("--sample", required=True)
	ap.add_argument("--out", required=True)
	ap.add_argument("--keep")
	ap.add_argument("--map")
	a = ap.parse_args()
	# Validate every genotype.
	sl = Path(a.sample).read_text().splitlines()
	n = max(0, len([x for x in sl[2:] if x.strip()]))
	exp = 2 * n
	nv = 0
	prev = -1
	errors = []
	first = None
	last = None
	ids = [x.split()[1] for x in sl[2:] if x.strip()]
	if len(set(ids)) != n:
		errors.append("duplicate sample IID")
	if a.keep:
		wanted = [x.split()[-1] for x in Path(a.keep).read_text().splitlines() if x.strip() and not x.startswith("#")]
		if ids != wanted:
			errors.append("sample identities/order differ from selected panel")
	maps = open(a.map) if a.map else None
	prev_cm = -float("inf")
	with haps_open(a.haps) as h:
		for line in h:
			if not line.strip():
				continue
			z = line.split()
			nv += 1
			# PLINK/Oxford HAPS has five metadata columns.
			meta = 5 if len(z) - 5 == exp else None
			if meta is None:
				errors.append(f"line {nv}: columns={len(z)} not metadata+2N ({exp})")
				if len(errors) > 10:
					break
				continue
			try:
				pos = int(float(z[2] if meta == 5 else z[3]))
			except Exception:
				errors.append(f"line {nv}: bad position")
				continue
			if pos <= prev:
				errors.append(f"line {nv}: positions not strictly increasing {prev}>={pos}")
			prev = pos
			first = pos if first is None else first
			last = pos
			bad = {x for x in z[meta:] if x not in {"0", "1"}}
			if bad:
				errors.append(f"line {nv}: invalid hap values {sorted(bad)[:5]}")
			if len(z[3]) != 1 or len(z[4]) != 1 or z[3] not in "ACGT" or z[4] not in "ACGT" or z[3] == z[4]:
				errors.append(f"line {nv}: invalid alleles")
			if maps:
				m = maps.readline().split()
				if (
					len(m) != 4
					or m[0].removeprefix("chr") != z[0].removeprefix("chr")
					or m[1] != z[1]
					or int(m[3]) != pos
				):
					errors.append(f"line {nv}: map/HAPS mismatch")
				elif float(m[2]) <= prev_cm:
					errors.append(f"line {nv}: genetic positions not strictly increasing")
				else:
					prev_cm = float(m[2])
			if len(errors) > 10:
				break
	if maps:
		if maps.readline():
			errors.append("extra map records")
		maps.close()
	q = {
		"individuals": n,
		"haplotypes": exp,
		"variants": nv,
		"first_position": first,
		"last_position": last,
		"errors": errors,
	}
	Path(a.out).parent.mkdir(parents=True, exist_ok=True)
	Path(a.out).write_text(json.dumps(q, indent=2) + "\n")
	if errors or nv == 0 or n == 0:
		raise SystemExit("HAPS validation failed: " + "; ".join(errors[:3]))
	if a.keep:
		# Inference identifies samples by IID; Oxford requires identical ID_1/ID_2.
		# Only normalize after proving the IID order matches the selected panel.
		rows = [x.split() for x in sl[2:] if x.strip()]
		if any(r[0] != r[1] for r in rows):
			for r in rows:
				r[0] = r[1]
			Path(a.sample).write_text("\n".join(sl[:2]) + "\n" + "".join(" ".join(r) + "\n" for r in rows))
	print(json.dumps(q))


def validate_haps_cli():
	validate_haps_main()


# 🚩 make-arg-keep: make_arg_keep
"""Create a deterministic, ancestry-balanced UKB ARG sample order."""
import argparse, gzip
from pathlib import Path
import numpy as np
import pandas as pd


def read_table(path: str) -> pd.DataFrame:
	return pd.read_csv(path, sep=None, engine="python", compression="infer", dtype=str)


def read_sample(path: str) -> pd.DataFrame:
	lines = [x for x in Path(path).read_text().splitlines() if x.strip()]
	if len(lines) < 2:
		raise SystemExit(f"Bad Oxford sample/PSAM file: {path}")
	h = lines[0].lstrip("#").split()
	# Oxford SAMPLE has a type row after the header; PLINK PSAM does not.
	start = 2 if lines[1].split() and all(x in {"0", "D", "B", "C", "P"} for x in lines[1].split()) else 1
	rows = [x.split() for x in lines[start:] if x.strip()]
	if not rows:
		raise SystemExit(f"No samples in {path}")
	d = pd.DataFrame(rows, columns=h[: len(rows[0])])
	lower = {x.lower(): x for x in d.columns}
	col = (
		lower.get("id_2")
		or lower.get("iid")
		or lower.get("id")
		or lower.get("id_1")
		or d.columns[min(1, len(d.columns) - 1)]
	)
	out = pd.DataFrame({"eid": d[col].astype(str)})
	if out.eid.duplicated().any():
		raise SystemExit("Duplicate sample IID is ambiguous")
	fid = lower.get("fid") or lower.get("id_1")
	out["fid"] = d[fid].astype(str) if fid else "0"
	out = out[~out.eid.isin(["0", "NA", "nan", ""])]
	out["source_order"] = np.arange(len(out), dtype=np.int64)
	return out


def norm_pop(x: pd.Series) -> pd.Series:
	z = (
		x.fillna("")
		.astype(str)
		.str.upper()
		.str.strip()
		.str.replace("_", " ", regex=False)
		.str.replace("-", " ", regex=False)
	)
	out = pd.Series("OTH", index=x.index, dtype="object")
	rules = {
		"AFR": r"^(?:AFR|AFRICAN|BLACK)",
		"EAS": r"^(?:EAS|EAST ASIAN|CHINESE|JAPANESE|KOREAN)",
		"EUR": r"^(?:EUR|EUROPEAN|WHITE)",
		"SAS": r"^(?:SAS|SOUTH ASIAN|INDIAN|PAKISTANI|BANGLADESHI)",
		"AMR": r"^(?:AMR|HISPANIC|LATINO|ADMIXED AMERICAN)",
	}
	for p, pat in rules.items():
		out[z.str.contains(pat, regex=True, na=False)] = p
	return out


def make_arg_keep_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--sample", required=True)
	ap.add_argument("--ancestry", required=True)
	ap.add_argument("--out-keep", required=True)
	ap.add_argument("--out-panel", required=True)
	ap.add_argument("--max-individuals", type=int, default=20000)
	ap.add_argument("--full", action="store_true")
	ap.add_argument("--anchors-per-pop", type=int, default=1000)
	ap.add_argument("--prob-min", type=float, default=0.999)
	ap.add_argument("--seed", type=int, default=20260904)
	ap.add_argument("--custom-keep", default="")
	a = ap.parse_args()
	rng = np.random.default_rng(a.seed)
	s = read_sample(a.sample)
	anc = read_table(a.ancestry)
	idc = next((c for c in ["eid", "IID", "#IID", "sample", "ID_2"] if c in anc.columns), None)
	if idc is None:
		raise SystemExit(f"No ID column in {a.ancestry}: {list(anc.columns)}")
	ac = next(
		(
			c
			for c in ["ancestry", "genetic_ancestry", "predicted_ancestry", "super_pop", "population"]
			if c in anc.columns
		),
		None,
	)
	if ac is None:
		raise SystemExit(f"No ancestry column in {a.ancestry}: {list(anc.columns)}")
	pc = next(
		(
			c
			for c in ["ancestry_prob", "ancestry_probability", "probability", "max_posterior", "posterior_prob"]
			if c in anc.columns
		),
		None,
	)
	anc = anc.rename(columns={idc: "eid"})
	anc["eid"] = anc.eid.astype(str)
	anc["pop"] = norm_pop(anc[ac])
	anc["prob"] = pd.to_numeric(anc[pc], errors="coerce") if pc else np.nan
	# If no single probability column exists, use posterior of the assigned group.
	if not np.isfinite(anc["prob"]).any():
		vals = []
		for _, r in anc.iterrows():
			c = f"posterior_{r['pop']}"
			vals.append(pd.to_numeric(r.get(c, np.nan), errors="coerce"))
		anc["prob"] = vals
	if anc.eid.duplicated().any():
		raise SystemExit("Duplicate sample IDs in ancestry table")
	d = s.merge(anc[["eid", "pop", "prob"]], on="eid", how="left", validate="one_to_one")
	d["pop"] = d["pop"].fillna("OTH")
	if a.custom_keep:
		k = pd.read_csv(a.custom_keep, sep=None, engine="python", comment="#", header=None, dtype=str)
		wanted = k.iloc[:, 1 if k.shape[1] > 1 else 0].astype(str).tolist()
		if len(set(wanted)) != len(wanted) or not set(wanted) <= set(d.eid):
			raise SystemExit("Custom keep contains duplicate or unknown sample IDs")
		rank = {x: i for i, x in enumerate(wanted)}
		d = d[d.eid.isin(rank)].copy()
		d["rank"] = d.eid.map(rank)
		d = d.sort_values("rank")
	else:
		from itertools import zip_longest

		anchor_lists = []
		for p in ["AFR", "EAS", "EUR", "SAS"]:
			x = d[(d["pop"] == p) & ((d["prob"] >= a.prob_min) | d["prob"].isna())].copy()
			x = x.sort_values(["prob", "source_order"], ascending=[False, True])
			anchor_lists.append(x.head(a.anchors_per_pop).eid.tolist())
		# Round-robin ordering ensures the initial ARG scaffold contains all groups.
		selected = [eid for row in zip_longest(*anchor_lists) for eid in row if eid is not None]
		selected = list(dict.fromkeys(selected))
		if not a.full and a.max_individuals > 0:
			selected = selected[: a.max_individuals]
		remaining = d[~d.eid.isin(selected)].copy()
		if a.full or a.max_individuals <= 0:
			fill = remaining.sort_values("source_order").eid.tolist()
		else:
			n = max(0, a.max_individuals - len(selected))
			if n > len(remaining):
				n = len(remaining)
			# Stratified random fill so the pilot is not almost entirely EUR.
			fill = []
			groups = [g for g in ["AFR", "EAS", "EUR", "SAS", "AMR", "OTH"] if (remaining["pop"] == g).any()]
			if groups and n:
				alloc = {g: min(len(remaining[remaining["pop"] == g]), n // len(groups)) for g in groups}
				for g, m in alloc.items():
					ids = remaining.loc[remaining["pop"] == g, "eid"].to_numpy()
					if m:
						fill.extend(rng.choice(ids, size=m, replace=False).tolist())
				remn = n - len(fill)
				pool = remaining[~remaining.eid.isin(fill)].eid.to_numpy()
				if remn:
					fill.extend(rng.choice(pool, size=min(remn, len(pool)), replace=False).tolist())
		selected = selected + fill
		order = {x: i for i, x in enumerate(selected)}
		d = d[d.eid.isin(order)].copy()
		d["rank"] = d.eid.map(order)
		d = d.sort_values("rank")
	if not len(d):
		raise SystemExit("No ARG samples selected")
	Path(a.out_keep).parent.mkdir(parents=True, exist_ok=True)
	with open(a.out_keep, "w") as h:
		h.write("#FID\tIID\n")
		h.writelines(f"{r.fid}\t{r.eid}\n" for r in d.itertuples())
	d[["eid", "pop", "prob", "source_order"]].to_csv(a.out_panel, sep="\t", index=False)
	print(f"selected_individuals={len(d)} pop_counts={d['pop'].value_counts().to_dict()}")


def make_arg_keep_cli():
	make_arg_keep_main()


# 🚩 run-argneedle-advanced: run_argneedle_advanced
"""Request-bound, validated ARG-Needle checkpoints and visible progress."""
import argparse, fcntl, hashlib, json, os, shlex, signal, subprocess, sys, time
from pathlib import Path


def locate(home):
	paths = list(Path(home).glob("**/infer_args_advanced.py")) if home else []
	try:
		import importlib.util

		s = importlib.util.find_spec("arg_needle.scripts.infer_args_advanced")
		if s and s.origin:
			paths.append(Path(s.origin))
	except (ImportError, AttributeError, ValueError):
		pass
	return next((p.resolve() for p in paths if p.is_file()), None)


def stamp(p):
	p = Path(p).resolve(strict=True)
	s = p.stat()
	return [str(p), s.st_size, s.st_mtime_ns]


def run_argneedle_advanced_main():
	def stopped(signum, frame):
		raise SystemExit(128 + signum)

	signal.signal(signal.SIGTERM, stopped)
	p = argparse.ArgumentParser()
	p.add_argument("--home", default="")
	for k in ("haps", "map", "out", "log"):
		p.add_argument("--" + k, required=True)
	p.add_argument("--chr", type=int, required=True)
	p.add_argument("--seed-haplotypes", type=int, required=True)
	p.add_argument("--threads", type=int, default=8)
	p.add_argument("--normalize", type=int, choices=(0, 1), default=1)
	p.add_argument("--random-seed", type=int, default=20260904)
	p.add_argument("--mode", choices=("array", "sequence"), default="array")
	p.add_argument("--normalize-demography", default="")
	p.add_argument("--replace", action="store_true")
	a = p.parse_args()
	script = locate(a.home)
	if script is None:
		raise SystemExit("Cannot find infer_args_advanced.py")
	env = os.environ.copy()
	if a.home:
		env["PYTHONPATH"] = a.home + os.pathsep + env.get("PYTHONPATH", "")
	env["PYTHONUNBUFFERED"] = "1"
	env["PYTHONHASHSEED"] = str(a.random_seed)
	helptext = subprocess.check_output([sys.executable, str(script), "--help"], env=env, text=True)
	prefix = Path(a.out)
	prefix.parent.mkdir(parents=True, exist_ok=True)
	Path(a.log).parent.mkdir(parents=True, exist_ok=True)
	lock = open(str(prefix) + ".lock", "a")
	try:
		fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
	except BlockingIOError:
		raise SystemExit(f"Another inference uses {prefix}")
	sample = a.haps.removesuffix(".gz").removesuffix(".haps").removesuffix(".hap") + ".sample"
	n = 2 * len(Path(sample).read_text().splitlines()[2:])
	if not 2 <= a.seed_haplotypes <= n:
		raise SystemExit("Invalid scaffold size")
	import importlib.metadata
	import arg_needle

	package = Path(arg_needle.__file__).parent
	memory_adapter = Path(__file__).resolve().with_name("0.arg.py")
	backend = backend_path()
	if not backend.is_file():
		raise SystemExit(
			f"Missing bounded ASMC backend: {backend}; run f/0.arg.py build-asmc-memory --source /path/to/ASMC-1.4.0"
		)
	dependencies = [
		memory_adapter,
		backend,
		package / "inference.py",
		package / "decoders.py",
		package / "resources/30-100-2000_CEU.decodingQuantities.gz",
		package / "resources/CEU.demo",
	]
	request = {
		"schema": 3,
		"inputs": [stamp(x) for x in (a.haps, a.map, sample, script, __file__, *dependencies) if Path(x).is_file()],
		"versions": {name: importlib.metadata.version(name) for name in ("arg-needle", "arg-needle-lib")},
		"mode": a.mode,
		"scaffold": a.seed_haplotypes,
		"seed": a.random_seed,
		"normalize": a.normalize,
		"chr": a.chr,
		"trim": 0,
		"demography": stamp(a.normalize_demography) if a.normalize_demography else "default CEU",
	}
	rid = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()
	paths = {1: Path(str(prefix) + ".step1.argn"), 2: Path(str(prefix) + ".step2.argn"), 3: Path(str(prefix) + ".argn")}
	bootstrap = "import runpy,sys,random,numpy as np; adapter=sys.argv.pop(1); runpy.run_path(adapter)['install'](); seed=int(sys.argv.pop(1)); random.seed(seed); np.random.seed(seed); script=sys.argv.pop(1); sys.argv[0]=script; runpy.run_path(script,run_name='__main__')"
	dirty = a.replace
	with open(a.log, "a", buffering=1) as log:
		th = next((x for x in ("--num_threads", "--threads", "--n_threads") if x in helptext), None)
		if th is None:
			msg = "ARG-Needle inference is single-threaded in this version; --threads controls preparation; --jobs controls chromosome concurrency."
			print(msg, flush=True)
			log.write(msg + "\n")
		for step in (1, 2, 3):
			output = paths[step]
			marker = Path(str(prefix) + f".step{step}.done")
			valid = False
			if not dirty and marker.is_file() and output.is_file():
				try:
					valid = (
						json.loads(marker.read_text()) == {"request": rid, "output": stamp(output)}
						and output.stat().st_size > 0
					)
				except (ValueError, OSError):
					pass
			if valid:
				print(f"SKIP verified chr{a.chr} step={step}", flush=True)
				continue
			dirty = True
			for ds in range(step, 4):
				Path(str(prefix) + f".step{ds}.done").unlink(missing_ok=True)
				paths[ds].unlink(missing_ok=True)
			cmd = [
				sys.executable,
				"-u",
				"-c",
				bootstrap,
				str(memory_adapter),
				str(a.random_seed),
				str(script),
				"--hap_gz",
				a.haps,
				"--map",
				a.map,
				"--out",
				str(prefix),
				"--mode",
				a.mode,
				"--step",
				str(step),
				"--chromosome",
				str(a.chr),
				"--verbose",
				"1",
				"--trim_num_snps",
				"0",
			]
			if step == 1:
				cmd += ["--num_snp_samples" if a.mode == "array" else "--num_sequence_samples", str(a.seed_haplotypes)]
			if step == 3:
				cmd += ["--normalize", str(a.normalize)]
			if a.normalize_demography:
				cmd += ["--normalize_demography", a.normalize_demography]
			if th:
				cmd += [th, str(a.threads)]
			log.write("COMMAND " + shlex.join(cmd) + "\n")
			print(f"INFER chr{a.chr} step={step}/3 haplotypes={n} scaffold={a.seed_haplotypes} log={a.log}", flush=True)
			proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
			try:
				code = proc.wait()
				if code:
					raise subprocess.CalledProcessError(code, cmd)
			finally:
				if proc.poll() is None:
					proc.terminate()
					try:
						proc.wait(timeout=10)
					except subprocess.TimeoutExpired:
						proc.kill()
						proc.wait()
			if not output.is_file() or output.stat().st_size == 0:
				raise SystemExit(f"Step {step} did not create {output}")
			import arg_needle_lib

			arg = arg_needle_lib.deserialize_arg(str(output))
			if arg.num_samples() != (a.seed_haplotypes if step == 1 else n):
				raise SystemExit(f"Step {step}: wrong sample count")
			del arg
			tmp = Path(str(marker) + ".next")
			tmp.write_text(json.dumps({"request": rid, "output": stamp(output)}) + "\n")
			tmp.replace(marker)
	print(paths[3], flush=True)


def run_argneedle_advanced_cli():
	run_argneedle_advanced_main()


# 🚩 argn-to-trees: argn_to_trees
"""Convert an ARG-Needle .argn object to a mutation-bearing tskit file."""
import argparse, gzip, importlib, json, sys
from pathlib import Path


def is_ts(x):
	return hasattr(x, "dump") and hasattr(x, "num_trees") and hasattr(x, "samples")


def unwrap(x):
	if is_ts(x):
		return x
	if isinstance(x, (tuple, list)):
		for y in x:
			z = unwrap(y)
			if z is not None:
				return z
	if isinstance(x, dict):
		for y in x.values():
			z = unwrap(y)
			if z is not None:
				return z
	return None


def add_haps_mutations(ts, path):
	import tskit

	"""Place HAPS alleles parsimoniously when arg_to_tskit omits mutations."""
	if not path or ts.num_sites:
		return ts, 0
	metadata = ts.metadata if isinstance(ts.metadata, dict) else {}
	offset = int(metadata.get("offset", 0))
	samples = ts.samples()
	tables = ts.dump_tables()
	tables.sites.metadata_schema = tskit.MetadataSchema.permissive_json()
	added = 0
	opener = gzip.open if str(path).endswith(".gz") else open
	with opener(path, "rt") as src:
		for line in src:
			z = line.split()
			if len(z) < 6:
				continue
			try:
				bp = int(float(z[2]))
			except ValueError:
				continue
			pos = bp - offset
			if not 0 <= pos < ts.sequence_length:
				continue
			genotypes = [int(x) if x in ("0", "1") else -1 for x in z[5:]]
			if len(genotypes) != len(samples):
				raise SystemExit(
					f"HAPS sample count mismatch at {z[1]}: "
					f"{len(genotypes)} haplotypes versus {len(samples)} ARG samples"
				)
			ancestral, mutations = ts.at(pos).map_mutations(genotypes, alleles=[z[3], z[4]])
			site_id = tables.sites.add_row(
				position=pos,
				ancestral_state=ancestral,
				metadata={"bp": bp, "id": z[1], "allele0": z[3], "allele1": z[4]},
			)
			mutation_base = tables.mutations.num_rows
			for mutation in mutations:
				parent = tskit.NULL if mutation.parent == tskit.NULL else mutation_base + mutation.parent
				tables.mutations.add_row(
					site=site_id,
					node=mutation.node,
					derived_state=mutation.derived_state,
					parent=parent,
				)
			added += 1
	tables.sort()
	return tables.tree_sequence(), added


def argn_to_trees_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--argn", required=True)
	ap.add_argument("--haps", default="")
	ap.add_argument("--out", required=True)
	ap.add_argument("--home", default="")
	a = ap.parse_args()
	if a.home:
		sys.path.insert(0, a.home)
	import arg_needle_lib

	arg = arg_needle_lib.deserialize_arg(a.argn)
	# Newer builds may expose a method directly.
	for name in ("to_tskit", "to_tree_sequence", "as_tskit"):
		fn = getattr(arg, name, None)
		if callable(fn):
			try:
				ts = unwrap(fn())
				if ts is not None:
					break
			except Exception:
				pass
	else:
		ts = None
	modules = [arg_needle_lib]
	for mn in ("arg_needle_lib.convert", "arg_needle_lib.utils"):
		try:
			modules.append(importlib.import_module(mn))
		except Exception:
			pass
	names = (
		"arg_to_ts",
		"arg_to_tree_sequence",
		"arg_to_tskit",
		"convert_arg_to_ts",
		"convert_to_tskit",
		"convert_arg",
	)
	errors = []
	if ts is None:
		for m in modules:
			for name in names:
				fn = getattr(m, name, None)
				if not callable(fn):
					continue
				attempts = [(arg,), (arg, None), (arg, []), (arg, None, None), (arg, [], [])]
				for aa in attempts:
					try:
						z = unwrap(fn(*aa))
						if z is not None:
							ts = z
							break
					except Exception as e:
						errors.append(f"{m.__name__}.{name}{len(aa)}:{type(e).__name__}:{e}")
				if ts is not None:
					break
			if ts is not None:
				break
	if ts is None:
		avail = []
		for m in modules:
			avail += [
				f"{m.__name__}.{x}"
				for x in dir(m)
				if "ts" in x.lower() or "tree" in x.lower() or "convert" in x.lower()
			]
		raise SystemExit(
			"No compatible ARG-to-tskit converter found. Available="
			+ ",".join(avail[:50])
			+ " Errors="
			+ " | ".join(errors[:10])
		)
	if ts.num_samples <= 0 or ts.num_trees <= 0:
		raise SystemExit("Converted tree sequence is empty")
	ts, added_sites = add_haps_mutations(ts, a.haps)
	tables = ts.dump_tables()
	tables.time_units = "generations"
	tables.provenances.add_row(
		record=json.dumps(
			{
				"schema_version": "1.0.0",
				"software": {"name": "refGen.argn_to_trees", "version": "1"},
				"parameters": {
					"source_method": "needle",
					"node_time_units": "generations",
					"ancestral_state_source": "tree_parsimony_not_external_AA",
					"mutations": "HAPS alleles mapped by parsimony",
				},
			},
			separators=(",", ":"),
		)
	)
	ts = tables.tree_sequence()
	Path(a.out).parent.mkdir(parents=True, exist_ok=True)
	temp = Path(a.out + ".next")
	ts.dump(temp)
	temp.replace(a.out)
	print(
		f"trees={ts.num_trees} nodes={ts.num_nodes} samples={ts.num_samples} sites={ts.num_sites} mutations={ts.num_mutations} added_sites={added_sites} sequence_length={ts.sequence_length}"
	)


def argn_to_trees_cli():
	argn_to_trees_main()


# 🚩 make-sample-map: make_sample_map
import argparse
from pathlib import Path
import pandas as pd


def sample_ids(path):
	lines = Path(path).read_text().splitlines()
	h = lines[0].split()
	rows = [x.split() for x in lines[2:] if x.strip()]
	ix = {x.lower(): i for i, x in enumerate(h)}
	j = ix.get("id_2", ix.get("iid", ix.get("id", ix.get("id_1", min(1, len(h) - 1)))))
	return [r[j] for r in rows]


def make_sample_map_main():
	import tskit

	ap = argparse.ArgumentParser()
	ap.add_argument("--trees", required=True)
	ap.add_argument("--sample", required=True)
	ap.add_argument("--out", required=True)
	a = ap.parse_args()
	ids = sample_ids(a.sample)
	ts = tskit.load(a.trees)
	nodes = list(map(int, ts.samples()))
	if len(nodes) != 2 * len(ids):
		raise SystemExit(f"tree samples={len(nodes)} != 2*individuals={2 * len(ids)}")
	out = []
	for i, eid in enumerate(ids):
		out.append((eid, 0, nodes[2 * i]))
		out.append((eid, 1, nodes[2 * i + 1]))
	pd.DataFrame(out, columns=["eid", "hap", "sample_node"]).to_csv(a.out, sep="\t", index=False)
	print(f"individuals={len(ids)} sample_nodes={len(nodes)}")


def make_sample_map_cli():
	make_sample_map_main()


# 🚩 make-anchors: make_anchors
import argparse
import pandas as pd


def make_anchors_main():
	ap = argparse.ArgumentParser()
	ap.add_argument("--panel", required=True)
	ap.add_argument("--sample-map", required=True)
	ap.add_argument("--out", required=True)
	ap.add_argument("--per-pop", type=int, default=1000)
	ap.add_argument("--prob-min", type=float, default=0.999)
	a = ap.parse_args()
	p = pd.read_csv(a.panel, sep="\t", dtype={"eid": str})
	m = pd.read_csv(a.sample_map, sep="\t", dtype={"eid": str})
	p["pop"] = p["pop"].astype(str).str.upper()
	p["prob"] = pd.to_numeric(p.get("prob"), errors="coerce")
	keep = []
	for pop in ["AFR", "EAS", "EUR", "SAS"]:
		x = (
			p[(p["pop"] == pop) & ((p["prob"] >= a.prob_min) | p["prob"].isna())]
			.sort_values(["prob", "source_order"], ascending=[False, True])
			.head(a.per_pop)
		)
		if len(x) < 10:
			raise SystemExit(f"Only {len(x)} anchors for {pop}")
		y = m.merge(x[["eid", "prob"]], on="eid", how="inner")
		y["pop"] = pop
		keep.append(y)
	z = pd.concat(keep, ignore_index=True)
	z[["eid", "hap", "sample_node", "pop", "prob"]].to_csv(a.out, sep="\t", index=False)
	print(z.groupby("pop").size().to_dict())


def make_anchors_cli():
	make_anchors_main()


# 🚩 arg-features: arg_features
"""Extract mutation age and population-level local genealogy from a tskit ARG."""
import argparse, gzip, math
from pathlib import Path
import numpy as np, pandas as pd

POPS = ["AFR", "EAS", "EUR", "SAS"]
PAIRS = [("AFR", "EAS"), ("AFR", "EUR"), ("AFR", "SAS"), ("EAS", "EUR"), ("EAS", "SAS"), ("EUR", "SAS")]


def features_open(path):
	return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_haps(path):
	d = {}
	with features_open(path) as h:
		for line in h:
			if not line.strip():
				continue
			z = line.split()
			if len(z) < 6:
				continue
			if z[2].replace(".", "", 1).isdigit():
				chrom, sid, bp, a0, a1 = z[0], z[1], int(float(z[2])), z[3], z[4]
			elif z[3].replace(".", "", 1).isdigit():
				chrom, sid, bp, a0, a1 = z[0], z[2], int(float(z[3])), z[4], z[5]
			else:
				continue
			d.setdefault(bp, (sid, a0, a1))
	return d


def stat_arrays(ts, sets, windows):
	pairs = [(POPS.index(a), POPS.index(b)) for a, b in PAIRS]
	try:
		div = ts.divergence(sample_sets=sets, indexes=pairs, windows=windows, mode="branch", span_normalise=True)
	except TypeError:
		div = ts.divergence(sets, indexes=pairs, windows=windows, mode="branch")
	try:
		within = ts.diversity(sample_sets=sets, windows=windows, mode="branch", span_normalise=True)
	except TypeError:
		within = ts.diversity(sets, windows=windows, mode="branch")
	div = np.asarray(div)
	within = np.asarray(within)
	if div.ndim == 1:
		div = div[:, None]
	if within.ndim == 1:
		within = within[:, None]
	return div, within


def root_by_window(ts, windows):
	n = len(windows) - 1
	acc = np.zeros(n)
	span = np.zeros(n)
	for tree in ts.trees():
		l, r = tree.interval
		roots = list(tree.roots)
		rt = max((ts.node(x).time for x in roots), default=np.nan)
		if not np.isfinite(rt):
			continue
		a = max(0, int(np.searchsorted(windows, l, side="right") - 1))
		b = min(n - 1, int(np.searchsorted(windows, r, side="left")))
		for w in range(a, b + 1):
			ov = max(0, min(r, windows[w + 1]) - max(l, windows[w]))
			if ov:
				acc[w] += ov * rt
				span[w] += ov
	return np.divide(acc, span, out=np.full(n, np.nan), where=span > 0)


def mutation_info(ts, site):
	import tskit

	tree = ts.at(site.position)
	ages = []
	lengths = []
	for m in site.mutations:
		node = int(m.node)
		nt = ts.node(node).time
		par = tree.parent(node)
		pt = ts.node(par).time if par != tskit.NULL else np.nan
		mt = getattr(m, "time", np.nan)
		if mt is not None and np.isfinite(mt) and mt != tskit.UNKNOWN_TIME:
			age = float(mt)
		elif np.isfinite(pt):
			age = float((nt + pt) / 2)
		else:
			age = float(nt)
		ages.append(age)
		lengths.append(float(pt - nt) if np.isfinite(pt) else np.nan)
	return (
		max(ages) if ages else np.nan,
		max([x for x in lengths if np.isfinite(x)], default=np.nan),
		len(site.mutations),
	)


def entropy(v):
	x = np.asarray(v, float)
	x = x[x > 0]
	if len(x) <= 1:
		return 0.0
	p = x / x.sum()
	return float(-(p * np.log(p)).sum() / math.log(len(POPS)))


def arg_features_main():
	import tskit

	ap = argparse.ArgumentParser()
	ap.add_argument("--trees", required=True)
	ap.add_argument("--haps", required=True)
	ap.add_argument("--anchors", required=True)
	ap.add_argument("--out", required=True)
	ap.add_argument("--chr", required=True)
	ap.add_argument("--window-bp", type=int, default=1000000)
	a = ap.parse_args()
	ts = tskit.load(a.trees)
	hv = read_haps(a.haps)
	anc = pd.read_csv(a.anchors, sep="\t")
	sets = [anc.loc[anc["pop"] == p, "sample_node"].astype(int).tolist() for p in POPS]
	if any(len(x) < 2 for x in sets):
		raise SystemExit("Each population needs >=2 anchor haplotypes")
	L = float(ts.sequence_length)
	windows = np.arange(0, L, a.window_bp, dtype=float)
	if len(windows) == 0 or windows[0] != 0:
		windows = np.insert(windows, 0, 0.0)
	if windows[-1] != L:
		windows = np.append(windows, L)
	offset = int(ts.metadata.get("offset", 0)) if isinstance(ts.metadata, dict) else 0
	div, within = stat_arrays(ts, sets, windows)
	root = root_by_window(ts, windows)
	all_nodes = [x for s in sets for x in s]
	slices = []
	k = 0
	for s in sets:
		slices.append(slice(k, k + len(s)))
		k += len(s)
	sample_nodes = np.asarray(ts.samples(), dtype=int)
	node_to_index = {n: i for i, n in enumerate(sample_nodes)}
	anchor_sample_idx = np.array([node_to_index[n] for n in all_nodes], dtype=int)
	Path(a.out).parent.mkdir(parents=True, exist_ok=True)
	cols = [
		"chr",
		"bp",
		"SNP",
		"allele0",
		"allele1",
		"arg_age_gen",
		"mutation_branch_gen",
		"n_mutations",
		"root_time_gen",
		"carrier_frequency",
		"lineage_breadth",
		"lineage_entropy",
	]
	cols += [f"af_{p}" for p in POPS] + [f"div_{x}_{y}" for x, y in PAIRS] + [f"within_{p}" for p in POPS]
	n = 0
	missing_id = 0
	with gzip.open(a.out, "wt") as o:
		o.write("\t".join(cols) + "\n")
		for var in ts.variants(samples=sample_nodes, isolated_as_missing=False):
			site = var.site
			bp = int(round(site.position)) + offset
			info = hv.get(bp)
			if info is None:
				# ARG conversions can use 0-based site coordinates while HAPS is 1-based.
				info = hv.get(bp + 1) or hv.get(bp - 1)
			if info is None:
				sid, a0, a1 = f"{a.chr}:{bp}", "NA", "NA"
				missing_id += 1
			else:
				sid, a0, a1 = info
			wi = min(len(windows) - 2, max(0, int(np.searchsorted(windows, site.position, side="right") - 1)))
			age, bl, nmut = mutation_info(ts, site)
			g = np.asarray(var.genotypes)[anchor_sample_idx]
			counts = []
			af = []
			for sl in slices:
				z = g[sl]
				z = z[z >= 0]
				c = float(np.sum(z > 0))
				counts.append(c)
				af.append(c / max(1, len(z)))
			breadth = sum(x > 0 for x in counts)
			freq = float(np.mean(g > 0))
			row = [a.chr, bp, sid, a0, a1, age, bl, nmut, root[wi], freq, breadth, entropy(counts)]
			row += af + list(np.asarray(div[wi]).ravel()) + list(np.asarray(within[wi]).ravel())
			o.write("\t".join("NA" if (isinstance(x, float) and not np.isfinite(x)) else str(x) for x in row) + "\n")
			n += 1
	print(f"sites={n} missing_haps_id={missing_id} windows={len(windows) - 1}")


def arg_features_cli():
	arg_features_main()


# 🚩 arg-affinity: arg_affinity
"""Optional chromosome-averaged genealogical nearest-neighbour affinities."""
import argparse, gzip
from pathlib import Path
import numpy as np, pandas as pd

POPS = ["AFR", "EAS", "EUR", "SAS"]


def arg_affinity_main():
	import tskit

	ap = argparse.ArgumentParser()
	ap.add_argument("--trees", required=True)
	ap.add_argument("--sample-map", required=True)
	ap.add_argument("--anchors", required=True)
	ap.add_argument("--out", required=True)
	ap.add_argument("--batch", type=int, default=10000)
	a = ap.parse_args()
	ts = tskit.load(a.trees)
	m = pd.read_csv(a.sample_map, sep="\t", dtype={"eid": str})
	anc = pd.read_csv(a.anchors, sep="\t")
	refs = [anc.loc[anc["pop"] == p, "sample_node"].astype(int).tolist() for p in POPS]
	focal = m.sample_node.astype(int).to_numpy()
	Path(a.out).parent.mkdir(parents=True, exist_ok=True)
	with gzip.open(a.out, "wt") as o:
		o.write("eid\thap\tq_AFR\tq_EAS\tq_EUR\tq_SAS\n")
		for st in range(0, len(focal), a.batch):
			en = min(len(focal), st + a.batch)
			q = np.asarray(ts.genealogical_nearest_neighbours(focal[st:en], refs))
			q = np.nan_to_num(q, nan=0.0)
			den = q.sum(1)
			q = np.divide(q, den[:, None], out=np.full_like(q, 0.25), where=den[:, None] > 0)
			for (_, r), v in zip(m.iloc[st:en].iterrows(), q):
				o.write(f"{r.eid}\t{r.hap}\t" + "\t".join(map(str, v)) + "\n")
	print(f"haplotypes={len(focal)}")


def arg_affinity_cli():
	arg_affinity_main()


COMMANDS = {
	"build-asmc-memory": build_asmc_memory_cli,
	"build-argneedle-map": build_argneedle_map_cli,
	"validate-haps": validate_haps_cli,
	"make-arg-keep": make_arg_keep_cli,
	"run-argneedle-advanced": run_argneedle_advanced_cli,
	"argn-to-trees": argn_to_trees_cli,
	"make-sample-map": make_sample_map_cli,
	"make-anchors": make_anchors_cli,
	"arg-features": arg_features_cli,
	"arg-affinity": arg_affinity_cli,
}


def main():
	import sys

	if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help", "help"):
		print(__doc__)
		print("Commands: " + ", ".join(COMMANDS))
		return
	command = sys.argv.pop(1)
	if command not in COMMANDS:
		raise SystemExit("Unknown command: " + command)
	sys.argv[0] += " " + command
	COMMANDS[command]()


if __name__ == "__main__":
	main()
