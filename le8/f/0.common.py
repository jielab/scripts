#!/usr/bin/env python3


# 🚩 0.common
# Shared LE8 execution dispatcher and aggregate evidence catalogue.
#
# Use ``dispatch`` (the default) for pipeline modules and ``index`` for the catalogue.
#
# The index reads aggregate results without fitting models; dispatch runs the selected
# pipeline stages. Measured-significant proteins are not an eligibility gate for the
# genetic catalogue. Records retain their comparison scope, source row and hash.
"""Shared LE8 execution dispatcher and aggregate evidence catalogue.

Use ``dispatch`` (the default) for pipeline modules and ``index`` for the catalogue.

The index reads aggregate results without fitting models; dispatch runs the selected
pipeline stages. Measured-significant proteins are not an eligibility gate for the
genetic catalogue. Records retain their comparison scope, source row and hash.
"""

from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
import shutil
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from contextlib import contextmanager, ExitStack


def _load_index_dependencies():
	global np, pd
	import numpy as np
	import pandas as pd


if __name__ != "__main__":
	_load_index_dependencies()


VERSION = "2026-10-03.1"
MODULES = [
	"c1_correlate",
	"c2_cause",
	"c3_coloc",
	"c4_connect",
	"c5_cellulation",
	"final",
]
PRIVATE_PARTS = re.compile(r"^(?:_|\.|private|cache|input|neural|checkpoints)", re.I)
PRIVATE_NAMES = re.compile(
	r"(?:^|[._])(roles|predictions|outcomes|individual|participants|context|split|matches|early_stop_ids|reference_candidates)(?:[._]|$)",
	re.I,
)
PRIVATE_COLUMNS = {
	"eid",
	"iid",
	"fid",
	"participant_id",
	"individual_id",
	"subject_id",
	"patient_id",
	"query_id",
	"target_id",
	"reference_id",
	"reference_eid",
	"target_eid",
	"query_eid",
	"eid_a",
	"eid_b",
	"donor_id",
	"date_birth",
	"birth_date",
}


def has_private_columns(columns):
	names = [str(name).lower() for name in columns]
	# Aggregate model tables may identify the comparator with reference_id.
	if {"run_id", "model_id", "primary_id", "reference_id", "n", "uno_c_horizon"} <= set(names):
		names = [name for name in names if name != "reference_id"]
	private = {re.sub(r"[^a-z0-9]", "", name) for name in PRIVATE_COLUMNS}
	private.update({"id", "id1", "id2", "sampleid", "ukbid", "ukbeid"})
	return any(
		re.sub(r"[^a-z0-9]", "", name) in private or re.search(r"(^|[_. #])(eid|iid|fid)($|[_. ])", name)
		for name in names
	)


FEATURE_COLUMNS = (
	"feature",
	"assay",
	"exposure",
	"protein",
	"Protein",
	"metabolite",
	"Metabolite",
	"term",
	"gene",
)
# Fields which identify an estimate, not an individual.
META = [
	"Y",
	"layer",
	"feature",
	"scope",
	"adjustment",
	"model",
	"series",
	"kind",
	"N",
	"events",
	"units",
	"source_id",
	"source_file",
	"source_row",
	"sample_hash",
]


def final_dir(root, trait, layer=None, kind=None):
	path = Path(root) / "final" / trait
	if layer is not None:
		path /= layer
	if kind is not None:
		path /= kind
	return path


def result_file(root, trait, layer, module, filename):
	if module == "c1_prediction":
		return (
			Path(root) / trait / layer / "c1_correlate" / "prediction" / re.sub(r"^final[._]", "", filename)
		)
	if module == "c4_panel_validation":
		module = "c4_connect"
	return Path(root) / trait / layer / module / filename


def digest(path: Path) -> str:
	h = hashlib.sha256()
	with path.open("rb") as f:
		for block in iter(lambda: f.read(1024 * 1024), b""):
			h.update(block)
	return h.hexdigest()


def bh(values, family_n=None):
	p = np.asarray(values, dtype=float)
	q = np.full(len(p), np.nan)
	ok = np.flatnonzero(np.isfinite(p) & (p >= 0) & (p <= 1))
	ix = ok[np.argsort(p[ok], kind="stable")]
	n = max(len(p), int(family_n or len(p)))
	if len(ix):
		q[ix] = np.minimum(
			1,
			np.minimum.accumulate((p[ix] * n / np.arange(1, len(ix) + 1))[::-1])[::-1],
		)
	return q


def normalize_feature(value):
	if pd.isna(value):
		return ""
	v = str(value).strip().strip("`")
	v = re.sub(r"^(?:prot|met)__", "", v)
	v = re.sub(r"\.pgs$", "", v, flags=re.I)
	# No mapping of NTPROBNP to NPPB, nor of percentages to concentrations.
	return "Lactate" if v.casefold() in {"lactat", "lactate"} else v


def series(d, names, default=float("nan")):
	for name in names:
		if name in d:
			return d[name]
	return pd.Series(default, index=d.index)


def numeric(d, names, default=float("nan")):
	return pd.to_numeric(series(d, names, default), errors="coerce")


def safe_path(path: Path, root: Path) -> bool:
	try:
		rel = path.resolve().relative_to(root.resolve())
	except ValueError:
		return False
	return not any(PRIVATE_PARTS.search(x) for x in rel.parts[:-1])


def source_module(path):
	parts = path.parts
	if "c5_cellulation" in parts or re.search(r"^c5\.(?:(?:systematic[.])?cell[._]|CIGMA)", path.name):
		return "c5_cellulation"
	if any(x in parts for x in ["final_prediction", "final"]):
		return "final"
	if "c4_connect" in parts or path.name.startswith("c4."):
		return "c4_connect"
	for module in MODULES[:3]:
		if module in parts:
			return module
	if re.match(r"^final[.](proxy_|domain_support)", path.name):
		return "c4_connect"
	if path.name.startswith("final."):
		return "final"
	if "abm" in str(path).lower():
		return "c1_correlate"
	return "final"


def aggregate_allowed(path, root):
	if not path.is_file() or not safe_path(path, root):
		return False
	# These complete regional arrays remain in the reusable C3 fit. Publication
	# exports only the actual plotted loci/variants to aggregate workbooks.
	if path.name in {"c3.regional_rows.csv", "c3.variant_posteriors.csv"}:
		return False
	if PRIVATE_NAMES.search(path.name):
		return False
	if re.search(
		r"(?:manifest|provenance|output_index|omission_audit|genetic_score_weights)",
		path.name,
		re.I,
	):
		return False
	if not re.search(r"\.(csv|tsv)(\.gz)?$", path.name, re.I):
		return False
	return bool(
		re.match(
			r"^(?:c[0-5][._]|pwas_|mwas_|Fig|final[._]|test_metrics|raw_test_metrics|approach_comparison|support_error_audit|coverage_curve|paired_contrasts|subgroup_metrics|masked_feature_metrics|learning_curve|reference_readiness|cohort_audit)",
			path.name,
			re.I,
		)
	)


def atomic_csv(d, path):
	path.parent.mkdir(parents=True, exist_ok=True)
	temp = path.with_name("." + path.name + ".tmp")
	# Empty tables still have an explicit schema, never a one-byte fake result.
	d.to_csv(temp, index=False)
	temp.replace(path)


class Index:
	def __init__(self, root, out, traits, layers, max_mb=128):
		self.root = root.resolve()
		self.out = out.resolve()
		self.traits = traits
		self.layers = layers
		if self.out == self.root:
			raise ValueError("Index destination must not replace the analysis root")
		self.out.mkdir(parents=True, exist_ok=True)
		self.registry = []
		self.errors = []
		self.files = {}
		self.candidates = []
		self.effects = []
		self.loci = []
		self.proxy = []
		self.counts = []
		self.max_bytes = max_mb * 1024 * 1024
		self.origins = []

	def register_only(self, p, Y, layer):
		if p.stat().st_size > self.max_bytes:
			return
		sep = "\t" if ".tsv" in p.name else ","
		try:
			h = pd.read_csv(p, sep=sep, nrows=0)
			if has_private_columns(h.columns):
				return
			rel = str(p.relative_to(self.root))
			sid = hashlib.sha256(rel.encode()).hexdigest()[:16]
			if sid in self.files:
				return
			self.registry.append(
				dict(
					source_id=sid,
					Y=Y,
					layer=layer,
					module=source_module(p),
					role="browse_aggregate",
					path=rel,
					rows=np.nan,
					columns=len(h.columns),
					bytes=p.stat().st_size,
					sha256=digest(p),
					modified_utc=datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat(),
				)
			)
			self.files[sid] = p
		except (OSError, ValueError, pd.errors.ParserError) as exc:
			self.errors.append(
				dict(
					Y=Y,
					layer=layer,
					role="browse_aggregate",
					status="unreadable",
					detail=f"{p}: {exc}",
				)
			)

	def read(self, paths, Y, layer, role, required=False):
		available = [p for p in paths if p.is_file()]
		if not available:
			if required:
				self.errors.append(
					dict(
						Y=Y,
						layer=layer,
						role=role,
						status="unavailable",
						detail="No recognized aggregate table",
					)
				)
			return pd.DataFrame()
		p = available[0]
		# All canonical input paths are constructed internally. Never accept a URL as a path.
		if not safe_path(p, self.root):
			self.errors.append(dict(Y=Y, layer=layer, role=role, status="blocked", detail=str(p)))
			return pd.DataFrame()
		if p.stat().st_size > self.max_bytes:
			self.errors.append(dict(Y=Y, layer=layer, role=role, status="too_large", detail=str(p)))
			return pd.DataFrame()
		sep = "\t" if ".tsv" in p.name else ","
		try:
			head = pd.read_csv(p, sep=sep, nrows=0)
			if has_private_columns(head.columns):
				raise ValueError("Participant-level identifier column; not exposed in Shiny")
			d = pd.read_csv(p, sep=sep, low_memory=False)
		except Exception as exc:
			self.errors.append(
				dict(
					Y=Y,
					layer=layer,
					role=role,
					status="unreadable_or_private",
					detail=f"{p}: {exc}",
				)
			)
			return pd.DataFrame()
		rel = str(p.relative_to(self.root))
		h = digest(p)
		sid = hashlib.sha256(rel.encode()).hexdigest()[:16]
		if sid not in self.files:
			row = dict(
				source_id=sid,
				Y=Y,
				layer=layer,
				module=source_module(p),
				role=role,
				path=rel,
				rows=len(d),
				columns=len(d.columns),
				bytes=p.stat().st_size,
				sha256=h,
				modified_utc=datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat(),
			)
			self.registry.append(row)
			self.files[sid] = p
		if len(available) > 1 and any(digest(q) != h for q in available[1:]):
			self.errors.append(
				dict(
					Y=Y,
					layer=layer,
					role=role,
					status="alternative_conflict",
					detail="Preferred " + rel + "; alternatives differ, not combined",
				)
			)
		d["_source_id"] = sid
		d["_source_file"] = rel
		d["_source_row"] = np.arange(2, len(d) + 2)
		return d

	def paths(self, Y, layer, module, name):
		b = self.root / Y / layer
		if module == "c1_prediction":
			return [result_file(self.root, Y, layer, module, name)]
		return [result_file(self.root, Y, layer, module, name)]

	def metadata(self, d, Y, layer, feature=None):
		z = pd.DataFrame(index=d.index)
		z["Y"] = Y
		z["layer"] = layer
		z["feature"] = (feature if feature is not None else series(d, FEATURE_COLUMNS, "")).map(
			normalize_feature
		)
		for c, opts, default in [
			("scope", ["scope"], "existing"),
			("adjustment", ["adjustment"], "unspecified"),
			("model", ["model"], "separate"),
			("N", ["N", "N_total"], np.nan),
			("events", ["events", "N_event"], np.nan),
			("sample_hash", ["sample_hash"], ""),
		]:
			z[c] = series(d, opts, default)
		for c in ["source_id", "source_file", "source_row"]:
			z[c] = d["_" + c]
		return z

	def add_effect(
		self,
		d,
		Y,
		layer,
		kind,
		label,
		beta,
		se,
		p,
		q=None,
		units="Log HR per own SD",
		feature=None,
	):
		if d.empty:
			return
		z = self.metadata(d, Y, layer, feature)
		z["kind"] = kind
		z["series"] = label if isinstance(label, str) else label
		z["beta"] = beta
		z["se"] = se
		z["p"] = p
		z["FDR"] = np.nan if q is None else q
		z["lo"] = z.beta - 1.96 * z.se
		z["hi"] = z.beta + 1.96 * z.se
		z["units"] = units
		# Keep a reported coefficient even when its CI is missing; do not manufacture precision.
		self.effects.append(z.loc[z.feature.ne("") & np.isfinite(z.beta)])

	def run_layer(self, Y, layer):
		base = self.root / Y / layer
		prefix = "pwas" if layer == "prot" else "mwas"
		obs = self.read(
			self.paths(Y, layer, "c1_correlate", prefix + "_incident_adj2.csv"),
			Y,
			layer,
			"incident",
			True,
		)
		paired = self.read(
			self.paths(Y, layer, "c1_correlate", "c1.pgs_focus.paired.csv")
			+ [final_dir(self.root, Y) / f"Fig1.PGS_integrated.{layer}_paired.csv"],
			Y,
			layer,
			"matched_PGS",
			True,
		)
		master = pd.DataFrame()
		if len(paired) and {"feature", "measured_p", "pgs_p"} <= set(paired):
			master = self.metadata(paired, Y, layer)
			for new, old in [
				("measured_beta", "measured_beta"),
				("measured_se", "measured_se"),
				("measured_p", "measured_p"),
				("measured_FDR", "measured_FDR"),
				("pgs_beta", "pgs_beta"),
				("pgs_se", "pgs_se"),
				("pgs_p", "pgs_p"),
				("pgs_FDR", "pgs_FDR"),
				("pgs_joint_p", "pgs_joint_p"),
				("pgs_joint_FDR", "pgs_joint_FDR"),
				("pgs_measured_r", "correlation"),
			]:
				master[new] = numeric(paired, [old])
			# Recompute a missing FDR only across the ENTIRE original scope/adjustment family.
			for pcol, qcol in [("measured_p", "measured_FDR"), ("pgs_p", "pgs_FDR")]:
				if qcol not in paired:
					for _, ix in master.groupby(["scope", "adjustment"], dropna=False).groups.items():
						master.loc[ix, qcol] = bh(master.loc[ix, pcol])
			master["same_people"] = series(paired, ["comparison"], "").str.contains(
				"identical people", case=False, na=False
			)
			master["comparison"] = series(paired, ["comparison"], "Unverified matching")
			master["genetic_evidence"] = "PGS association; discovery overlap must be checked separately"
			self.add_effect(
				paired,
				Y,
				layer,
				"measured_vs_PGS",
				"Measured",
				master.measured_beta,
				master.measured_se,
				master.measured_p,
				master.measured_FDR,
			)
			self.add_effect(
				paired,
				Y,
				layer,
				"measured_vs_PGS",
				"PGS",
				master.pgs_beta,
				master.pgs_se,
				master.pgs_p,
				master.pgs_FDR,
			)
		if len(obs) and {"term", "p.value"} <= set(obs):
			self.add_effect(
				obs,
				Y,
				layer,
				"incident",
				"Measured (adjusted)",
				numeric(obs, ["beta"]),
				numeric(obs, ["std.error"]),
				numeric(obs, ["p.value"]),
				numeric(obs, ["FDR"]),
			)
			# Retain all measured-only assays, not only genetic matches.
			seen = set(master.feature) if len(master) else set()
			extra = obs[~obs.term.map(normalize_feature).isin(seen)].copy()
			if len(extra):
				m = self.metadata(extra, Y, layer, extra.term)
				for k, cols in [
					("measured_beta", ["beta"]),
					("measured_se", ["std.error"]),
					("measured_p", ["p.value"]),
					("measured_FDR", ["FDR"]),
				]:
					m[k] = numeric(extra, cols)
				m["same_people"] = False
				m["comparison"] = "No full matched-PGS table available for this assay"
				m["genetic_evidence"] = "Unavailable, not negative"
				master = pd.concat([master, m], ignore_index=True)
		# Single-biomarker anchors are a fallback display only. They are NOT an all-feature screening family.
		anchors = self.read(
			self.paths(Y, layer, "c1_correlate", "c1.paired_pgs_measured.csv"),
			Y,
			layer,
			"paired_anchors",
		)
		if len(anchors) and {"feature", "component", "beta", "std.error"} <= set(anchors):
			self.add_effect(
				anchors,
				Y,
				layer,
				"paired_anchors",
				anchors.component.astype(str),
				numeric(anchors, ["beta"]),
				numeric(anchors, ["std.error"]),
				numeric(anchors, ["p.value"]),
				numeric(anchors, ["FDR"]),
			)
		comp = self.read(
			self.paths(Y, layer, "c1_correlate", "c1.pgs_focus.components.csv")
			+ [final_dir(self.root, Y) / f"Fig1.PGS_integrated.{layer}_components.csv"],
			Y,
			layer,
			"components",
		)
		if len(comp) and {"feature", "term", "beta", "se"} <= set(comp):
			self.add_effect(
				comp,
				Y,
				layer,
				"G_R",
				comp.term.astype(str),
				numeric(comp, ["beta"]),
				numeric(comp, ["se"]),
				numeric(comp, ["p"]),
				numeric(comp, ["FDR"]),
				units="Log HR per whole-biomarker SD; calibrated components",
			)
		con = self.read(
			self.paths(Y, layer, "c1_correlate", "c1.pgs_focus.contrasts.csv")
			+ [final_dir(self.root, Y) / f"Fig1.PGS_integrated.{layer}_contrasts.csv"],
			Y,
			layer,
			"component_contrasts",
		)
		if len(con) and {"feature", "beta_difference", "se_difference"} <= set(con):
			self.add_effect(
				con,
				Y,
				layer,
				"G_minus_R",
				"G - R",
				numeric(con, ["beta_difference"]),
				numeric(con, ["se_difference"]),
				numeric(con, ["p"]),
				numeric(con, ["FDR"]),
				units="Difference in calibrated component log HR; conditional CI",
			)
		calibration = self.read(
			self.paths(Y, layer, "c1_correlate", "c1.pgs_focus.calibration.csv")
			+ [final_dir(self.root, Y) / f"Fig1.PGS_integrated.{layer}_calibration.csv"],
			Y,
			layer,
			"PGS_calibration",
		)
		# Temporal plots keep windows/landmarks separate; this is not repeated sampling.
		for role, file in [
			("landmark", prefix + "_incident_landmark_adj2.csv"),
			("event_window", prefix + "_diagnosis_window_riskset_adj2.csv"),
		]:
			t = self.read(self.paths(Y, layer, "c1_correlate", file), Y, layer, role)
			if len(t) and "beta" in t:
				ss = series(
					t,
					["landmark_years", "landmark", "window", "period", "interval"],
					"unspecified",
				).astype(str)
				if {"window_lo", "window_hi"} <= set(t):
					ss = (
						series(t, ["side"], "").astype(str)
						+ " | "
						+ t.window_lo.astype(str)
						+ " to "
						+ t.window_hi.astype(str)
						+ " years"
					)
				units = series(t, ["effect_measure"], "Log HR per SD; landmark-specific risk set")
				if role == "event_window" and "side" in t:
					prev = t.side.astype(str).str.contains("Pre-baseline", case=False, na=False)
					for mask, kind in [
						(~prev, "event_window"),
						(prev, "prevalent_window"),
					]:
						u = t.loc[mask]
						self.add_effect(
							u,
							Y,
							layer,
							kind,
							ss.loc[mask],
							numeric(u, ["beta"]),
							numeric(u, ["std.error", "se"]),
							numeric(u, ["p.value", "p"]),
							numeric(u, ["FDR"]),
							units=units.loc[mask],
						)
				else:
					self.add_effect(
						t,
						Y,
						layer,
						role,
						ss,
						numeric(t, ["beta"]),
						numeric(t, ["std.error", "se"]),
						numeric(t, ["p.value", "p"]),
						numeric(t, ["FDR"]),
						units=units,
					)
		co = self.read(
			self.paths(Y, layer, "c3_coloc", "c3.coloc_summary.csv")
			+ self.paths(Y, layer, "c3_coloc", "c3.credible_set_by_locus.csv"),
			Y,
			layer,
			"coloc",
		)
		cm = pd.DataFrame()
		if len(co) and "feature" in co:
			zz = self.metadata(co, Y, layer)
			zz["locus"] = series(co, ["locus", "region"], "").astype(str)
			zz["PP_H4"] = numeric(co, ["PP.H4", "PP.H4_p12_default"])
			zz["PP_H4_robust_min"] = numeric(co, ["PP.H4_robust_min", "PP.H4_p12_conservative"])
			zz["credible_set_n"] = numeric(co, ["credible_set_n"])
			zz["n_snps"] = numeric(co, ["n_snps", "nsnps"])
			zz["status"] = series(co, ["status"], "reported")
			zz["locus_class"] = series(co, ["locus_class", "analysis"], "unverified")
			zz["interpretation"] = (
				"Shared regional association is not causal direction; H4 SNP posterior is not trait fine-mapping"
			)
			self.loci.append(zz)
			# Display ranking only; no combined P or pooled independent-locus count.
			zz["_display_rank"] = zz.PP_H4_robust_min.fillna(zz.PP_H4)
			zz = zz[zz._display_rank.notna()]
			ix = zz.groupby("feature")._display_rank.idxmax().dropna().astype(int)
			cm = zz.loc[ix, ["feature", "PP_H4", "PP_H4_robust_min", "locus", "credible_set_n"]].rename(
				columns={"locus": "best_reported_coloc_region"}
			)
		mr = self.read(self.paths(Y, layer, "c2_cause", "c2.MR_all.csv"), Y, layer, "MR")
		mm = pd.DataFrame()
		if len(mr) and {"exposure", "pval"} <= set(mr):
			ana = series(mr, ["analysis", "instrument_class"], "").astype(str)
			cls = np.where(
				ana.str.contains(r"cis|local", case=False, regex=True)
				& ~ana.str.contains(r"trans|distal", case=False, regex=True),
				"cis_local",
				"other",
			)
			self.add_effect(
				mr,
				Y,
				layer,
				"MR",
				ana + " / " + series(mr, ["method"], "").astype(str),
				numeric(mr, ["b", "beta"]),
				numeric(mr, ["se"]),
				numeric(mr, ["pval"]),
				numeric(mr, ["FDR_all", "FDR"]),
				units="MR effect in original exposure units; do not compare magnitude to PGS association",
				feature=mr.exposure,
			)
			mr = mr.assign(
				_feature=mr.exposure.map(normalize_feature),
				_q=numeric(mr, ["FDR_all"]),
				_class=cls,
			)
			cis = mr[mr._class.eq("cis_local") & mr._q.notna()]
			if len(cis):
				mm = (
					cis.sort_values("_q")
					.drop_duplicates("_feature")[["_feature", "_q"]]
					.rename(
						columns={
							"_feature": "feature",
							"_q": "best_reported_cis_MR_FDR_all",
						}
					)
				)
		reverse = self.read(
			self.paths(Y, layer, "c2_cause", "c2.reverse_MR_all.csv"),
			Y,
			layer,
			"reverse_MR",
		)
		if len(reverse) and {"feature", "b", "se", "pval"} <= set(reverse):
			self.add_effect(
				reverse,
				Y,
				layer,
				"MR_reverse",
				"Disease liability → biomarker",
				numeric(reverse, ["b"]),
				numeric(reverse, ["se"]),
				numeric(reverse, ["pval"]),
				numeric(reverse, ["FDR_reverse"]),
				units="Reverse MR: disease liability, not the effect of diagnosed disease or treatment",
			)
		if master.empty and not cm.empty:
			master = cm[["feature"]].copy()
			master["Y"] = Y
			master["layer"] = layer
			master["scope"] = "coloc_only"
			master["adjustment"] = "unavailable"
			master["same_people"] = False
		if master.empty:
			return
		if len(con) and {"feature", "scope", "beta_difference", "p", "FDR"} <= set(con):
			g = con.copy()
			g["feature"] = g.feature.map(normalize_feature)
			counts = g.groupby(["feature", "scope"]).size()
			good = counts[counts.eq(1)].reset_index()[["feature", "scope"]]
			g = g.merge(good, on=["feature", "scope"])
			g = g[["feature", "scope", "beta_difference", "p", "FDR", "_source_file"]].rename(
				columns={
					"beta_difference": "GR_beta_difference",
					"p": "GR_difference_p",
					"FDR": "GR_difference_FDR",
					"_source_file": "GR_source_file",
				}
			)
			master = master.merge(g, on=["feature", "scope"], how="left", validate="many_to_one")
		if len(cm):
			master = master.merge(cm, on="feature", how="left", validate="many_to_one")
		if len(mm):
			master = master.merge(mm, on="feature", how="left", validate="many_to_one")
		for c in [
			"measured_p",
			"measured_FDR",
			"pgs_p",
			"pgs_FDR",
			"GR_beta_difference",
			"GR_difference_p",
			"GR_difference_FDR",
			"best_reported_cis_MR_FDR_all",
			"PP_H4_robust_min",
		]:
			if c not in master:
				master[c] = np.nan
		master["measured_P_ge_005_PGS_FDR_lt_005"] = (
			master.same_people.fillna(False) & master.measured_p.ge(0.05) & master.pgs_FDR.lt(0.05)
		)
		master["measured_FDR_ge_005_PGS_FDR_lt_005"] = (
			master.same_people.fillna(False) & master.measured_FDR.ge(0.05) & master.pgs_FDR.lt(0.05)
		)
		master["measured_P_ge_005_cis_MR_FDR_lt_005"] = master.measured_p.ge(
			0.05
		) & master.best_reported_cis_MR_FDR_all.lt(0.05)
		om = master.measured_FDR.lt(0.05)
		pg = master.pgs_FDR.lt(0.05)
		enough = master.measured_FDR.notna() & master.pgs_FDR.notna()
		master["evidence_pattern"] = np.select(
			[
				enough & om & pg,
				enough & ~om & pg,
				enough & om & ~pg,
				enough & ~om & ~pg,
			],
			["Both supported", "PGS only", "Measured only", "Neither FDR-supported"],
			default="Unavailable/incomplete",
		)
		master["PGS_association_is_MR"] = False
		master["coloc_MR_same_locus_confirmed"] = (
			"Not inferred by min-P/max-H4 joins; inspect exact locus/instrument table"
		)
		master["causal_status"] = "Not classified from observational/PGS significance alone"
		master["candidate_id"] = [
			hashlib.sha256("|".join(map(str, v)).encode()).hexdigest()[:20]
			for v in master[["Y", "layer", "feature", "scope", "adjustment"]].itertuples(
				index=False, name=None
			)
		]
		if master.candidate_id.duplicated().any():
			self.errors.append(
				dict(
					Y=Y,
					layer=layer,
					role="candidates",
					status="duplicate_scope",
					detail="Duplicate candidate comparison keys; records retained with source-row suffix",
				)
			)
			master["candidate_id"] = master.candidate_id + "_" + master.index.astype(str)
		self.candidates.append(master)
		for (scope, adj), g in master.groupby(["scope", "adjustment"], dropna=False):
			self.counts.append(
				dict(
					Y=Y,
					layer=layer,
					scope=scope,
					adjustment=adj,
					assays=g.feature.nunique(),
					matched_PGS_assays=int(g.same_people.fillna(False).sum()),
					measured_nominal_null_PGS_supported=int(g.measured_P_ge_005_PGS_FDR_lt_005.sum()),
					measured_FDR_null_PGS_supported=int(g.measured_FDR_ge_005_PGS_FDR_lt_005.sum()),
					measured_nominal_null_cis_MR_supported=int(g.measured_P_ge_005_cis_MR_FDR_lt_005.sum()),
				)
			)

	def run_proxy(self, Y, layer):
		d = self.read(
			self.paths(Y, layer, "c4_connect", "c4.focus.proxy_accuracy.csv"),
			Y,
			layer,
			"proxy_accuracy",
		)
		if not {"model", "component", "R2_basic_omics", "R2_omics", "N"} <= set(d):
			return
		d["budget"] = pd.to_numeric(d.model.str.extract(r"_(\d+)$")[0], errors="coerce")
		ns = d[d.model.str.match(r"^NS_\d+$")]
		ys = d[d.model.str.match(r"^YS(?:plus|balanced)?_(?:Yin|YinYang)_\d+$")]
		cols = ["component", "budget", "N", "R2_basic_omics", "R2_omics"]
		z = ys.merge(
			ns[cols],
			on=["component", "budget"],
			suffixes=("", "_NS"),
			validate="many_to_one",
		)
		z["Y"] = Y
		z["layer"] = layer
		z["delta_R2_vs_NS"] = z.R2_basic_omics - z.R2_basic_omics_NS
		z["delta_R2_omics_vs_NS"] = z.R2_omics - z.R2_omics_NS
		z["same_N"] = z.N.eq(z.N_NS)
		z.loc[~z.same_N, ["delta_R2_vs_NS", "delta_R2_omics_vs_NS"]] = np.nan
		z["uncertainty"] = (
			"Point contrast only; paired CI requires prediction-level paired bootstrap, not subtraction of marginal CIs"
		)
		z["interpretation"] = (
			"Held-out LE8 reconstruction, not intervention responsiveness; all eight targets retained"
		)
		ci = self.read(
			self.paths(Y, layer, "c4_connect", "c4.explain.contrasts.csv"),
			Y,
			layer,
			"paired_proxy_uncertainty",
		)
		if len(ci) and {
			"model",
			"component",
			"budget",
			"N",
			"delta_R2_vs_NS",
			"lo",
			"hi",
		} <= set(ci):
			use = [
				"model",
				"component",
				"budget",
				"N",
				"delta_R2_vs_NS",
				"lo",
				"hi",
			] + [c for c in ["p", "FDR_all_proxy_contrasts", "valid_boot"] if c in ci]
			ci = ci[use].rename(columns={"N": "N_bootstrap", "delta_R2_vs_NS": "delta_R2_refit"})
			z = z.merge(
				ci,
				on=["model", "component", "budget"],
				how="left",
				validate="many_to_one",
			)
			valid = z.N.eq(z.N_bootstrap) & np.isclose(
				z.delta_R2_vs_NS, z.delta_R2_refit, rtol=1e-5, atol=1e-7
			)
			z.loc[
				~valid,
				[c for c in ["lo", "hi", "p", "FDR_all_proxy_contrasts"] if c in z],
			] = np.nan
			z.loc[valid, "uncertainty"] = (
				"Paired bootstrap of reconstructed original frozen-panel predictions; training and selection uncertainty excluded"
			)
			if (z.N_bootstrap.notna() & ~valid).any():
				self.errors.append(
					dict(
						Y=Y,
						layer=layer,
						role="proxy_uncertainty",
						status="refit_mismatch",
						detail="Paired refit does not reproduce source N/point contrast; CI withheld.",
					)
				)
		self.proxy.append(z)

	def catalogue(self):
		known = {str(p.resolve()) for p in self.files.values()}
		gallery = []
		hashes = set()
		scopes = [(Y, self.root / Y) for Y in self.traits]
		scopes += [(Y, final_dir(self.root, Y)) for Y in self.traits]
		scopes.append(("all", self.root / "final"))
		for Y, yd in scopes:
			if not yd.is_dir():
				continue
			for directory, dirs, filenames in os.walk(yd, followlinks=False):
				dirs[:] = [
					name
					for name in dirs
					if not PRIVATE_PARTS.search(name) and not (Y == "all" and Path(directory) == yd)
				]
				for filename in filenames:
					p = Path(directory) / filename
					if not re.search(r"\.(png|jpg|jpeg|pdf|csv|tsv)(\.gz)?$", filename, re.I):
						continue
					if not p.is_file() or not safe_path(p, self.root):
						continue
					rel = str(p.relative_to(self.root))
					parts = p.relative_to(yd).parts
					if Y != "all" and parts[0] in {"prot", "met"} and parts[0] not in self.layers:
						continue
					layer = "all" if Y == "all" else parts[0] if parts[0] in {"prot", "met"} else "joint"
					# Historical root-level cell tables annotate protein assays only.
					if layer == "joint" and source_module(p) == "c5_cellulation":
						layer = "prot"
					if layer in {"prot", "met"} and layer not in self.layers:
						continue
					if p.suffix.lower() in {
						".png",
						".jpg",
						".jpeg",
						".pdf",
					} and not PRIVATE_NAMES.search(p.name):
						if p.stat().st_size > self.max_bytes:
							continue
						if not re.search(r"(Fig|figure|enrichment)", p.name, re.I):
							continue
						h = digest(p)
						key = (Y, layer, h)
						if key in hashes:
							continue
						hashes.add(key)
						gallery.append(
							dict(
								Y=Y,
								layer=layer,
								module=source_module(p),
								path=rel,
								name=p.name,
								sha256=h,
								bytes=p.stat().st_size,
							)
						)
					elif str(p.resolve()) not in known and aggregate_allowed(p, self.root):
						self.register_only(p, Y, layer)
		atomic_csv(
			pd.DataFrame(
				gallery,
				columns=["Y", "layer", "module", "path", "name", "sha256", "bytes"],
			),
			self.out / "figures.csv",
		)

	def finish(self):
		schemas = {
			"candidates": [
				"candidate_id",
				"Y",
				"layer",
				"feature",
				"scope",
				"adjustment",
				"same_people",
				"measured_p",
				"measured_FDR",
				"pgs_p",
				"pgs_FDR",
				"measured_P_ge_005_PGS_FDR_lt_005",
			],
			"effects": META + ["beta", "se", "lo", "hi", "p", "FDR"],
			"loci": META + ["locus", "PP_H4", "PP_H4_robust_min"],
			"proxy_comparisons": [
				"Y",
				"layer",
				"model",
				"component",
				"budget",
				"delta_R2_vs_NS",
				"same_N",
			],
			"discovery_counts": [
				"Y",
				"layer",
				"scope",
				"adjustment",
				"assays",
				"matched_PGS_assays",
				"measured_nominal_null_PGS_supported",
			],
		}
		for name, frames in [
			("candidates", self.candidates),
			("effects", self.effects),
			("loci", self.loci),
			("proxy_comparisons", self.proxy),
			("discovery_counts", [pd.DataFrame(self.counts)]),
		]:
			d = (
				pd.concat(frames, ignore_index=True)
				if frames and any(len(x) for x in frames)
				else pd.DataFrame(columns=schemas[name])
			)
			if name == "effects":
				d = d.drop_duplicates(
					[
						"Y",
						"layer",
						"feature",
						"kind",
						"model",
						"series",
						"scope",
						"adjustment",
						"source_id",
						"source_row",
					]
				)
			atomic_csv(d, self.out / (name + ".csv"))
		atomic_csv(
			pd.DataFrame(
				self.registry,
				columns=[
					"source_id",
					"Y",
					"layer",
					"module",
					"role",
					"path",
					"rows",
					"columns",
					"bytes",
					"sha256",
					"modified_utc",
				],
			),
			self.out / "tables.csv",
		)
		atomic_csv(
			pd.DataFrame(self.errors, columns=["Y", "layer", "role", "status", "detail"]),
			self.out / "status.csv",
		)
		print(
			json.dumps(
				dict(
					index=str(self.out),
					candidates=sum(map(len, self.candidates)),
					tables=len(self.registry),
					issues=len(self.errors),
				)
			)
		)


def index_main(argv=None):
	_load_index_dependencies()
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("--analysis-root", type=Path, required=True)
	p.add_argument("--out", type=Path)
	p.add_argument("--Y", default="cvd_cad,ra")
	p.add_argument("--biom", default="prot,met")
	p.add_argument("--max-table-mb", type=int, default=128)
	p.add_argument(
		"--question-figures",
		action="store_true",
		help="Render question-led aggregate report figures",
	)
	a = p.parse_args(argv)
	traits = a.Y.split(",")
	layers = a.biom.split(",")
	if any(not re.fullmatch(r"[A-Za-z0-9_]+", x) for x in traits) or not set(layers) <= {"prot", "met"}:
		p.error("Invalid outcome/layer")
	if a.max_table_mb < 1:
		p.error("--max-table-mb must be positive")
	idx = Index(
		a.analysis_root,
		a.out or Path(os.getenv("LE8_BROWSER_INDEX", str(a.analysis_root / "shiny"))),
		traits,
		layers,
		a.max_table_mb,
	)
	for Y in traits:
		for layer in layers:
			idx.run_layer(Y, layer)
			idx.run_proxy(Y, layer)
	from final import build_questions, question_figures

	build_questions(
		idx.root,
		traits,
		layers,
		idx.out,
		pd.concat(idx.candidates, ignore_index=True) if idx.candidates else pd.DataFrame(),
		pd.concat(idx.proxy, ignore_index=True) if idx.proxy else pd.DataFrame(),
		idx.max_bytes,
	)
	if a.question_figures:
		question_figures(idx.root / "final", traits, layers)
	idx.catalogue()
	idx.finish()


# Cross-module execution and public CLI dispatch.
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NATIVE = {
	"c1_correlate",
	"c2_cause",
	"c3_coloc",
	"pgs_focus",
	"c4_connect",
}
ORDER = [
	"c1_correlate",
	"c1_abm",
	"c2_cause",
	"c3_coloc",
	"pgs_focus",
	"c4_connect",
	"c5_cellulation",
	"final",
	"shiny",
	"share",
]


def canonical_modules(csv):
	aliases = {
		"c1_connect": "c1_correlate",
		"c4_focus": "c4_connect",
		"c4_panel_validation": "c4_connect",
		"c4_explain": "c4_connect",
		"s1_interact": "c4_connect",
		"s1_interaction": "c4_connect",
		"s2_nonlin": "c4_connect",
	}
	return list(dict.fromkeys(aliases.get(x.strip(), x.strip()) for x in csv.split(",")))


def call(cmd, env, dry=False):
	print("[LE8] " + shlex.join(map(str, cmd)), flush=True)
	if not dry:
		subprocess.run(list(map(str, cmd)), env=env, check=True)


# Shared public configuration; explicit CLI > environment > fixed default.
def shared_analysis_settings(values, environ=None):
	env = os.environ if environ is None else environ
	spec = {
		"group_file": ("LE8_GROUP_FILE", ""),
		"group_col": ("LE8_GROUP_COLUMN", env.get("PGS_GROUP_COLUMN", "")),
		"end_date": ("DATE_FOLLOW_END", "2023-04-01"),
		"diagnosis_col": ("LE8_Y_DATE", ""),
		"outer_roster": ("LE8_OUTER_ROSTER", ""),
		"shared_covariates": ("LE8_VARS_ADJ", ""),
	}
	resolved, sources = {}, {}
	for key, (name, default) in spec.items():
		explicit = getattr(values, key, None)
		value = str(explicit) if explicit is not None else env.get(name, default)
		resolved[key] = value
		sources[key] = dict(
			source="cli" if explicit is not None else name if name in env else "default", value=value
		)
		if explicit is not None and name in env and value != env[name]:
			sources[key]["overridden_environment"] = env[name]
			print(f"[LE8] {key}: explicit CLI overrides {name}", file=sys.stderr)
	if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", resolved["end_date"]):
		raise ValueError("Administrative end date must be YYYY-MM-DD")
	datetime.strptime(resolved["end_date"], "%Y-%m-%d")
	return resolved, sources, {spec[k][0]: v for k, v in resolved.items()}


def parser():
	p = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog="""Examples:
  ./le8.sh final,shiny --Y cvd_cad,ra --biom prot,met
  ./le8.sh shiny --Y cvd_cad,ra --port 3839
  ./le8.sh final --index-only --Y cvd_cad,ra
  ./le8.sh c3_coloc,c4_connect,c5_cellulation --Y cvd_cad --biom prot
  C4_EXPLAIN_BOOT=200 ./le8.sh c4_connect --Y cvd_cad --biom prot
  ./le8.sh c1_abm --Y cvd_cad --biom prot
  ./le8.sh final --fit-joint --Y cvd_cad,ra --biom prot,met --replace TRUE

Modules: C1 correlation (c1_connect is accepted as an alias), C2 cause, C3 coloc,
C4 connect (including interactions, nonlinearity, panel validation and reconstruction bootstrap), C5 cellulation, final, shiny.
--abm-args is parsed as arguments, NEVER evaluated by a shell.
""",
	)
	p.add_argument(
		"modules",
		nargs="?",
		default="c1_correlate,c1_abm,c2_cause,c3_coloc,c4_connect,c5_cellulation",
	)
	p.add_argument("--Y", "--trait", "-Y", default=os.getenv("Y", "cvd_cad,ra"))
	p.add_argument("--biom", "-b", default=os.getenv("BIOM", "prot,met"))
	p.add_argument(
		"--analysis-root",
		type=Path,
		default=Path(os.getenv("LE8_ANALYSIS_ROOT", "/mnt/d/analysis/le8")),
	)
	p.add_argument(
		"--out",
		type=Path,
		help="Final figure destination; default <analysis-root>/final",
	)
	p.add_argument("--share-out", type=Path, help="Portable Shiny ZIP outside the analysis root")
	p.add_argument("--ukb-phe", type=Path, default=None)
	p.add_argument("--cores", type=int, default=None)
	for name, variable, default in [
		("pwas-workers", "LE8_PWAS_WORKERS", "auto"),
		("pgs-workers", "LE8_PGS_WORKERS", "auto"),
	]:
		p.add_argument("--" + name, default=os.getenv(variable, default))
	p.add_argument("--mrlink2-workers", type=int, default=int(os.getenv("MRLINK2_WORKERS", "4")))
	p.add_argument(
		"--eigh-backend", choices=["auto", "cpu", "cuda"], default=os.getenv("MRLINK2_EIGH_BACKEND", "auto")
	)
	p.add_argument("--mrlink2-inner-threads", type=int, default=int(os.getenv("MRLINK2_INNER_THREADS", "1")))
	p.add_argument("--cpu-task-memory-gib", type=float, default=None)
	p.add_argument("--gpu-jobs", type=int, default=int(os.getenv("LE8_GPU_JOBS", "1")))
	p.add_argument("--seed", type=int, default=int(os.getenv("SEED", "2026")))
	p.add_argument("--group-file", type=Path)
	p.add_argument("--group-col", "--group-column", dest="group_col")
	p.add_argument("--end-date")
	p.add_argument("--Y-date", "--diagnosis-col", dest="diagnosis_col")
	p.add_argument("--outer-roster", type=Path, help="Shared eid,role CSV/TSV; roles training/test")
	p.add_argument("--vars.adj", dest="shared_covariates")
	p.add_argument("--r-bin", default=os.getenv("R_BIN", "Rscript"))
	p.add_argument("--replace", choices=["TRUE", "FALSE"], default="FALSE")
	p.add_argument("--memory-limit-gb", type=int, default=None)
	p.add_argument("--memory-swap-gb", type=int, default=None)
	abm = p.add_mutually_exclusive_group()
	abm.add_argument("--run-abm", dest="run_abm", action="store_true", default=None)
	abm.add_argument("--skip-abm", dest="run_abm", action="store_false")
	for flag in [
		"fit-reference",
		"fit-joint",
		"fit-genetic",
		"details",
		"strict",
		"dry-run",
		"preflight",
		"index-only",
		"prepare-only",
		"no-reindex",
		"shiny-review",
	]:
		p.add_argument("--" + flag, action="store_true")
	for name in [
		"atlas",
		"universe",
		"panels",
		"contrasts",
		"cigma-manifest",
		"cigma-results",
		"cigma-cells",
	]:
		p.add_argument("--" + name, type=Path)
	p.add_argument("--matched-draws", type=int, default=int(os.getenv("C5_MATCHED_DRAWS", "500")))
	p.add_argument("--allow-untested-cigma", action="store_true")
	p.add_argument("--no-plots", action="store_true", help="C5 only: omit annotation plots")
	p.add_argument("--max-table-mb", type=int, default=128)
	p.add_argument("--port", type=int, default=int(os.getenv("LE8_SHINY_PORT", "3839")))
	p.add_argument("--host", default=os.getenv("LE8_SHINY_HOST", "127.0.0.1"))
	return p


def dispatch_main(argv=None):
	p = parser()
	a, extra = p.parse_known_args(argv)
	if a.cigma_results is None and a.cigma_manifest is None and os.getenv("C5_CIGMA_RESULTS"):
		a.cigma_results = Path(os.environ["C5_CIGMA_RESULTS"])
	try:
		shared, config_sources, shared_env = shared_analysis_settings(a)
	except ValueError as exc:
		p.error(str(exc))
	if a.matched_draws != 0 and a.matched_draws < 100:
		p.error("--matched-draws must be 0 or >=100")
	if (a.universe is None) != (a.panels is None):
		p.error("Provide --universe and --panels together")
	if a.cigma_manifest and a.cigma_results:
		p.error("Choose --cigma-manifest OR --cigma-results")
	if a.cigma_cells and not a.cigma_results:
		p.error("--cigma-cells requires --cigma-results")
	traits = a.Y.split(",")
	layers = a.biom.split(",")
	if any(not re.fullmatch(r"[A-Za-z0-9_]+", x) for x in traits) or len(traits) != len(set(traits)):
		p.error("Use unique, comma-separated outcome identifiers.")
	if not layers or not set(layers) <= {"prot", "met"} or len(layers) != len(set(layers)):
		p.error("Invalid --biom")
	if (a.cores is not None and a.cores < 1) or a.max_table_mb < 1 or not 1024 <= a.port <= 65535:
		p.error("Invalid resource/port setting")
	if any(x is not None and x < 0 for x in [a.memory_limit_gb, a.memory_swap_gb]):
		p.error("Memory settings must be nonnegative")
	if a.host not in {"127.0.0.1", "localhost", "::1"} and os.getenv("LE8_SHINY_ALLOW_NETWORK") != "TRUE":
		p.error("Network binding needs LE8_SHINY_ALLOW_NETWORK=TRUE and your own authentication/proxy.")
	requested = canonical_modules(a.modules)
	if any(x not in ORDER for x in requested):
		p.error("Unknown module; see --help.")
	if a.run_abm and "c1_abm" not in requested:
		requested.append("c1_abm")
	if a.run_abm is False:
		requested = [x for x in requested if x != "c1_abm"]
	modules = [x for x in ORDER if x in requested]
	if a.fit_reference:
		p.error(
			"--fit-reference moved to c1_correlate; its basic prediction benchmark now runs automatically"
		)
	fitted = a.fit_joint or a.fit_genetic or a.details
	if fitted and "final" not in modules:
		p.error("--fit-* and --details require final")
	if extra and not (set(modules) & NATIVE or fitted or modules == ["c1_abm"]):
		p.error("Unrecognized options: " + " ".join(extra))
	if a.index_only and fitted:
		p.error("--index-only cannot be combined with model fitting")
	env = dict(
		os.environ,
		LE8_ANALYSIS_ROOT=str(a.analysis_root.resolve()),
		LE8_FDIR=str(HERE),
		DIRSCRIPT=str(ROOT),
		SCRIPT_LOG_DIR=str(a.analysis_root.resolve() / "logs"),
		BIOM=a.biom,
		PROT_DO=str("prot" in layers).upper(),
		MET_DO=str("met" in layers).upper(),
		SEED=str(a.seed),
		LE8_REPLACE=a.replace,
		R_BIN=a.r_bin,
		LE8_SHINY_PORT=str(a.port),
		LE8_SHINY_HOST=a.host,
		LE8_SHINY_MAX_TABLE_MB=str(a.max_table_mb),
	)
	capacity = int(os.getenv("LE8_TOTAL_CORES", "16"))
	if capacity < 1:
		raise ValueError("LE8_TOTAL_CORES must be positive")
	if a.cores is None:
		a.cores = min(capacity, len(os.sched_getaffinity(0)))
	if a.cores > capacity:
		raise ValueError("--cores exceeds LE8_TOTAL_CORES shared pool")
	for name in ("pwas_workers", "pgs_workers"):
		v = getattr(a, name)
		if v != "auto" and (not str(v).isdigit() or int(v) < 1):
			raise ValueError("Invalid " + name)
	if min(a.mrlink2_workers, a.mrlink2_inner_threads) < 1 or a.gpu_jobs != 1:
		raise ValueError("Positive MR workers/threads and --gpu-jobs 1 required")
	if a.mrlink2_inner_threads > a.cores:
		raise ValueError("MR inner threads exceed allocated cores")
	if a.cpu_task_memory_gib is not None and (
		not math.isfinite(a.cpu_task_memory_gib) or a.cpu_task_memory_gib <= 0
	):
		raise ValueError("Invalid task memory")
	env.update(
		LE8_TOTAL_CORES=str(capacity),
		LE8_RESOURCE_MEMORY_GIB=os.getenv("LE8_RESOURCE_MEMORY_GIB", "32"),
		LE8_POOL_CPUS=os.getenv(
			"LE8_POOL_CPUS", ",".join(map(str, sorted(os.sched_getaffinity(0))[:capacity]))
		),
		LE8_PWAS_WORKERS=str(a.pwas_workers),
		LE8_PGS_WORKERS=str(a.pgs_workers),
		MRLINK2_WORKERS=str(min(a.mrlink2_workers, a.cores // a.mrlink2_inner_threads)),
		MRLINK2_INNER_THREADS=str(a.mrlink2_inner_threads),
		MRLINK2_EIGH_BACKEND=a.eigh_backend,
		LE8_GPU_JOBS=str(a.gpu_jobs),
	)
	if a.cpu_task_memory_gib is not None:
		env["LE8_CPU_TASK_MEMORY_GIB"] = str(a.cpu_task_memory_gib)
	env.update(shared_env)
	env["LE8_SHARED_CONFIG_SOURCES"] = json.dumps(config_sources, sort_keys=True)
	if a.ukb_phe:
		env["UKB_PHE"] = str(a.ukb_phe)
	if a.cores:
		env["N_CORES"] = str(a.cores)
	base = ["--Y", a.Y, "--biom", a.biom, "--analysis-root", str(a.analysis_root)]
	native = ["--seed", str(a.seed), "--replace", a.replace, "--r-bin", a.r_bin]
	if a.ukb_phe:
		native += ["--ukb-phe", str(a.ukb_phe)]
	if a.cores:
		native += ["--cores", str(a.cores)]
	if a.memory_limit_gb is not None:
		native += ["--memory-limit-gb", str(a.memory_limit_gb)]
	if a.memory_swap_gb is not None:
		native += ["--memory-swap-gb", str(a.memory_swap_gb)]
	if a.preflight:
		native += ["--preflight"]
	engine = HERE / "0.engine.sh"
	indexed = False

	def report_base():
		# Analysis scope remains explicit; aggregate reports include already completed scopes.
		available = [
			(d.name, l)
			for d in (a.analysis_root.iterdir() if a.analysis_root.is_dir() else [])
			if d.is_dir() and re.fullmatch(r"[A-Za-z0-9_]+", d.name)
			for l in ("prot", "met")
			if (d / l).is_dir()
		]
		ts = sorted(set(traits) | {t for t, _ in available})
		ls = sorted(set(layers) | {l for _, l in available})
		return ["--Y", ",".join(ts), "--biom", ",".join(ls), "--analysis-root", str(a.analysis_root)]

	def index():
		nonlocal indexed
		cmd = [
			sys.executable,
			HERE / "0.common.py",
			"index",
			*report_base(),
			"--max-table-mb",
			a.max_table_mb,
		]
		if "final" in modules and not a.index_only:
			cmd += ["--question-figures"]
		call(cmd, env, a.dry_run)
		indexed = True

	def phase(command, label, cores=None, gpu=False, memory=None):
		if a.preflight:
			return call(command, env, a.dry_run)
		budget = float(env["LE8_RESOURCE_MEMORY_GIB"])
		claim = memory if memory is not None else min(budget, float(a.memory_limit_gb or 32))
		wrapped = [
			sys.executable,
			HERE / "0.resources.py",
			"phase",
			"--label",
			label,
			"--cores",
			str(cores or a.cores),
			"--memory-gib",
			str(claim),
			*(["--gpu"] if gpu else []),
			"--",
			*command,
		]
		call(wrapped, env, a.dry_run)

	def performance(stage):
		if not a.preflight and not a.index_only:
			call(
				[
					sys.executable,
					HERE / "final.py",
					"performance",
					*(report_base() if stage == "final" else base),
					"--stage",
					stage,
					*(
						["--skip-completed"]
						if stage in {"c1_correlate", "c1_abm"} and a.replace != "TRUE"
						else []
					),
				],
				env,
				a.dry_run,
			)

	for module in modules:
		if module == "share":
			if modules != ["share"]:
				p.error("Run share separately after preparing final/shiny")
			output = a.share_out or a.analysis_root.parent / "le8-share.zip"
			if not a.dry_run:
				share_viewer(a.analysis_root, output, a.r_bin)
		elif module in NATIVE:
			# Native entry points may execute prerequisites and in-process R scans.
			# Reserve the full phase budget; MR subprocesses subdivide this lease.
			phase(
				["bash", engine, module, *base, *native, *extra],
				module,
				gpu=module == "c3_coloc" and env.get("RUN_GPU_COLOC", "TRUE") == "TRUE",
			)
		elif module == "c1_abm":
			call([sys.executable, HERE / "0.common.py", "pipeline", "c1_abm", *argv[1:]], env, a.dry_run)
		elif module == "c5_cellulation":
			cmd = [
				sys.executable,
				HERE / "c5.cellulation.py",
				*base,
				"--seed",
				str(a.seed),
				"--matched-draws",
				str(a.matched_draws),
			]
			for flag in ["preflight", "allow_untested_cigma", "no_plots"]:
				if getattr(a, flag):
					cmd += ["--" + flag.replace("_", "-")]
			if a.replace == "TRUE":
				cmd += ["--replace"]
			for key in [
				"atlas",
				"universe",
				"panels",
				"contrasts",
				"cigma_manifest",
				"cigma_results",
				"cigma_cells",
			]:
				if getattr(a, key) is not None:
					cmd += ["--" + key.replace("_", "-"), str(getattr(a, key))]
			for trait in traits:
				command = list(cmd)
				command[command.index("--Y") + 1] = trait
				call(command, env, a.dry_run)
		elif module == "final":
			fit_options = dict(
				traits=traits,
				layers=layers,
				seed=a.seed,
				native_options=extra,
				settings={k: v for k, v in os.environ.items() if k.startswith("FINAL_")},
			)
			if a.fit_joint or a.fit_genetic:
				ev = dict(
					env,
					FINAL_RUN_JOINT=str(a.fit_joint or a.fit_genetic).upper(),
					LE8_FINAL_REQUEST=json.dumps(fit_options, sort_keys=True),
				)
				call(
					["bash", engine, "final_prediction", *base, *native, *extra],
					ev,
					a.dry_run,
				)
			if a.details:
				call(["bash", engine, "final", *base, *native, *extra], env, a.dry_run)
			if not a.index_only and not a.preflight:
				cmd = [sys.executable, HERE / "final.py", "report", *report_base()]
				for key in ["out"]:
					if getattr(a, key) is not None:
						cmd += ["--" + key.replace("_", "-"), str(getattr(a, key))]
				if a.strict:
					cmd += ["--strict"]
				call(cmd, env, a.dry_run)
			performance("final")
			if not a.preflight:
				index()
		elif module == "shiny":
			if not a.preflight and not a.no_reindex and not indexed:
				index()
			if a.prepare_only:
				continue
			if a.preflight:
				call(
					[
						a.r_bin,
						"-e",
						"p<-c('shiny','DT','data.table','ggplot2','digest');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste('Missing',paste(m,collapse=', ')))",
					],
					env,
					a.dry_run,
				)
				continue
			print(
				f"[LE8] Shiny: http://{a.host}:{a.port} (foreground; Ctrl+C stops the server)",
				flush=True,
			)
			call(
				[
					a.r_bin,
					ROOT / "shiny/app.R",
					*(["--review"] if a.shiny_review else []),
				],
				env,
				a.dry_run,
			)
		if module in NATIVE | {"c5_cellulation"} and module != "pgs_focus":
			ready = module != "c1_correlate" or os.getenv("LE8_C1_STAGE", "c1_prediction") == "c1_prediction"
			ready &= module != "c4_connect" or os.getenv("LE8_C4_STAGE", "c4_explain") == "c4_explain"
			if ready:
				performance(module)


# 🚩 Temporary execution and consolidated publication
# Existing analytical modules exchange delimited files only in this private
# /tmp workspace. The published tree contains workbooks, readable configuration and named RDS data.
def run_table_storage(mode, work, r_bin, env):
	log = work.parent / (mode.removeprefix("--") + ".log")
	try:
		with log.open("w") as stream:
			subprocess.run(
				[r_bin, str(HERE / "0.common.R"), mode, str(work)],
				env=env,
				stdout=stream,
				stderr=subprocess.STDOUT,
				check=True,
			)
	except subprocess.CalledProcessError:
		print(log.read_text()[-4000:], file=sys.stderr)
		raise


@contextmanager
def analysis_lock(path, shared=False, waiting="Waiting for another LE8 run"):
	import fcntl

	with open(path, "a") as lock:
		mode = fcntl.LOCK_SH if shared else fcntl.LOCK_EX
		try:
			fcntl.flock(lock, mode | fcntl.LOCK_NB)
		except BlockingIOError:
			print(f"[LE8] {waiting}", flush=True)
			fcntl.flock(lock, mode)
		try:
			yield
		finally:
			fcntl.flock(lock, fcntl.LOCK_UN)


def owns_scope(relative, scopes):
	return scopes is None or (len(relative.parts) >= 2 and tuple(relative.parts[:2]) in scopes)


# 🚩 Reuse independently validated C1 outputs before opening a table workspace


def c1_result_signature(root, traits, layers):
	files = []
	for trait in traits:
		for layer in layers:
			directory = root / trait / layer / "c1_correlate"
			if not directory.is_dir():
				return None
			for path in sorted(directory.rglob("*")):
				if path.is_file():
					info = path.stat()
					link = path.lstat()
					files.append(
						[
							str(path.relative_to(root)),
							info.st_size,
							info.st_mtime_ns,
							info.st_ctime_ns,
							link.st_mtime_ns,
							str(path.resolve()),
						]
					)
	return dict(
		schema=1,
		root=str(root),
		traits=traits,
		layers=layers,
		files=files,
		validator=hashlib.sha256((HERE / "0.common.R").read_bytes()).hexdigest(),
	)


def skip_completed_c1(parsed):
	"""Cache successful native validation in /tmp; changed files require full checks."""
	if parsed.replace == "TRUE" or parsed.details or parsed.run_abm or os.getenv("LE8_C1_STAGE"):
		return False
	root = parsed.analysis_root.resolve()
	traits, layers = parsed.Y.split(","), parsed.biom.split(",")
	if any(not re.fullmatch(r"[A-Za-z0-9_]+", trait) for trait in traits) or not set(layers) <= {
		"prot",
		"met",
	}:
		raise ValueError("Invalid outcome/layer")
	key = hashlib.sha256(json.dumps([str(root), traits, layers]).encode()).hexdigest()[:24]
	cache = Path("/tmp/le8-completed")
	cache.mkdir(exist_ok=True)
	receipt = cache / (key + ".json")
	lock_name = hashlib.sha256(str(root).encode()).hexdigest()[:20]
	with ExitStack() as locks:
		locks.enter_context(analysis_lock(Path("/tmp") / f"le8-{lock_name}.lock", shared=True))
		for trait in sorted(set(traits)):
			for layer in sorted(set(layers)):
				scope = hashlib.sha256(f"{root}/{trait}/{layer}".encode()).hexdigest()[:20]
				locks.enter_context(analysis_lock(Path("/tmp") / f"le8-scope-{scope}.lock"))
		before = c1_result_signature(root, traits, layers)
		if before is None:
			return False
		try:
			saved = json.loads(receipt.read_text())
		except (OSError, ValueError):
			saved = {}
		if saved.get("signature") != before:
			with tempfile.TemporaryDirectory(prefix="le8-check-completed-", dir="/tmp") as temporary:
				audit = Path(temporary) / "results.csv"
				result = subprocess.run(
					[
						parsed.r_bin,
						str(HERE / "0.common.R"),
						"--check-completed",
						str(root),
						parsed.Y,
						parsed.biom,
						"c1_correlate,c1_prediction",
						str(audit),
					],
					text=True,
					capture_output=True,
					check=True,
				)
				if result.stderr:
					print(result.stderr, end="", file=sys.stderr)
				with audit.open() as handle:
					rows = list(csv.DictReader(handle))
				if c1_result_signature(root, traits, layers) != before:
					raise RuntimeError(
						"C1 outputs changed during validation; retry after the writer finishes"
					)
				saved = dict(signature=before, audit=rows)
				pending = receipt.with_suffix(f".{os.getpid()}.part")
				try:
					pending.write_text(json.dumps(saved))
					os.replace(pending, receipt)
				finally:
					pending.unlink(missing_ok=True)
		validated = {
			(row["trait"], row["layer"], row["module"])
			for row in saved.get("audit", [])
			if row["status"] == "completed"
		}
		expected = {
			(trait, layer, module)
			for trait in traits
			for layer in layers
			for module in ("c1_correlate", "c1_prediction")
		}
		if not expected <= validated:
			return False
		# The report keeps its original comparison scope unless explicitly replaced.
		from final import performance_complete

		if not all(
			performance_complete(root / trait / layer / "c1_correlate/performance")
			for trait in traits
			for layer in layers
		):
			return False
		for trait, layer, module in sorted(expected):
			print(
				f"[LE8] SKIP {module} Y={trait} biom={layer}: validated results unchanged; no table expansion or fitting",
				flush=True,
			)
		return True


def _table_runtime(argv):
	def option(name, default=None):
		for i, value in enumerate(argv):
			if value == name and i + 1 < len(argv):
				return argv[i + 1]
			if value.startswith(name + "="):
				return value.split("=", 1)[1]
		return default

	if os.getenv("LE8_TABLE_WORKSPACE") == "1" or any(
		x in argv for x in ["--help", "-h", "--dry-run", "--preflight"]
	):
		return False
	root = Path(option("--analysis-root", os.getenv("LE8_ANALYSIS_ROOT", "/mnt/d/analysis/le8"))).resolve()
	r_bin = option("--r-bin", os.getenv("R_BIN", "Rscript"))
	arguments = argv[1:] if argv[:1] == ["dispatch"] else argv
	if argv[:1] == ["index"]:
		module_arg = "index"
		traits = set(option("--Y", "cvd_cad,ra").split(","))
		layers = set(option("--biom", "prot,met").split(","))
		train_abm = False
	else:
		parsed, _ = parser().parse_known_args(arguments)
		module_arg = parsed.modules
		traits = set(parsed.Y.split(","))
		layers = set(parsed.biom.split(","))
		train_abm = parsed.run_abm is not False and (parsed.run_abm or "c1_abm" in module_arg.split(","))
	if (
		not traits
		or any(not re.fullmatch(r"[A-Za-z0-9_]+", y) for y in traits)
		or not layers
		or not layers <= {"prot", "met"}
	):
		raise ValueError("Invalid outcome/layer")
	modules = set(canonical_modules(module_arg))
	if modules == {"c1_correlate"} and not train_abm and skip_completed_c1(parsed):
		return True
	if modules == {"c4_connect"} and not os.getenv("LE8_C4_STAGE"):
		# One public C4 entry, with a committed transaction for each substage.
		# Later validation/bootstrap failures must retain earlier fitted results.
		try:
			for stage in ("c4_connect", "c4_panel_validation", "c4_explain"):
				print(f"[LE8] C4 stage: {stage}; publish before the next stage", flush=True)
				os.environ["LE8_C4_STAGE"] = stage
				_table_runtime(argv)
		finally:
			os.environ.pop("LE8_C4_STAGE", None)
		return True
	scopes = (
		None if modules & {"final", "shiny", "index", "share"} else {(y, b) for y in traits for b in layers}
	)
	if (modules == {"c1_abm"} or (modules == {"c1_correlate"} and not train_abm)) and len(scopes) > 1:
		# A fresh R process releases the preceding layer's full-cohort PGS heap.
		# Publish each successful scope before starting another so its numerical
		# caches and figures survive a later layer/trait failure.
		for y in dict.fromkeys(parsed.Y.split(",")):
			for b in dict.fromkeys(parsed.biom.split(",")):
				print(
					f"[LE8] {module_arg} transaction: Y={y} biom={b}; publish before the next scope",
					flush=True,
				)
				_table_runtime([*argv, "--Y", y, "--biom", b])
		return True
	if modules == {"c1_correlate"} and not os.getenv("LE8_C1_STAGE"):
		# Commit association results before prediction, so an interrupted benchmark
		# resumes without discarding a successful PWAS/MWAS run.
		try:
			for stage in ("c1_correlate", "c1_prediction"):
				print(f"[LE8] C1 stage: {stage}; publish before the next stage", flush=True)
				os.environ["LE8_C1_STAGE"] = stage
				_table_runtime(argv)
		finally:
			os.environ.pop("LE8_C1_STAGE", None)
		return True
	launch_shiny = "shiny" in module_arg.split(",") and "--prepare-only" not in argv
	if module_arg == "shiny" and "--no-reindex" in argv:
		return False
	if module_arg == "share":
		return False
	root.mkdir(parents=True, exist_ok=True)
	lock_name = hashlib.sha256(str(root).encode()).hexdigest()[:20]
	publication_lock = Path("/tmp") / f"le8-{lock_name}-publish.lock"
	locks = ExitStack()
	try:
		locks.enter_context(
			analysis_lock(
				Path("/tmp") / f"le8-{lock_name}.lock",
				shared=scopes is not None,
				waiting=f"Waiting for analysis/report transaction at {root}.",
			)
		)
		for y, b in sorted(scopes or []):
			key = hashlib.sha256(f"{root}/{y}/{b}".encode()).hexdigest()[:20]
			locks.enter_context(
				analysis_lock(
					Path("/tmp") / f"le8-scope-{key}.lock",
					waiting=f"{y}/{b} is already running at {root}; waiting for that scope.",
				)
			)
		return _run_table_workspace(
			argv, root, r_bin, module_arg, launch_shiny, train_abm, scopes, publication_lock, locks, option
		)
	finally:
		locks.close()


def _run_table_workspace(
	argv, root, r_bin, module_arg, launch_shiny, train_abm, scopes, publication_lock, locks, option
):
	import fcntl

	workspace = Path(tempfile.mkdtemp(prefix="le8-run-", dir="/tmp"))
	work = workspace / "results"
	work.mkdir()
	initial = {}
	try:
		print(f"[LE8] Temporary tables and execution files: {workspace}", flush=True)
		with analysis_lock(publication_lock, waiting="Waiting to read a consistent result snapshot."):
			for path in root.rglob("*"):
				if not path.is_file() or any(
					x
					in {
						"_history",
						"_source_figures",
						"_previous",
						"__pycache__",
						".ruff_cache",
						".pgs_focus_cache",
						"prepared",
						"neural",
						"quality_neural",
						"checkpoints",
						"le8_annotations",
						"logs",
					}
					for x in path.relative_to(root).parts
				):
					continue
				relative = path.relative_to(root)
				parts = relative.parts
				# twin publishes its own transactional directories and LOCAL_ONLY files.
				# The native workbook exporter must never repack these artifacts.
				if any(
					x
					in {
						"c1_abm",
						"abm_twin",
						"molecular_deviation",
						"baseline_comparison",
						"evidence",
						"_pgs_fullprecision",
					}
					or x.startswith("_abm_staging_")
					for x in parts
				):
					continue
				if not owns_scope(relative, scopes):
					continue
				stat = path.stat()
				initial[str(relative)] = (stat.st_size, stat.st_mtime_ns)
				target = work / relative
				target.parent.mkdir(parents=True, exist_ok=True)
				# An independent copy protects published fits on failed or interrupted runs.
				shutil.copy2(path, target)
		env = dict(
			os.environ,
			LE8_TABLE_WORKSPACE="1",
			LE8_ANALYSIS_ROOT=str(work),
			LE8_PUBLISHED_ROOT=str(root),
			LE8_TABLE_ABM_PRIVATE=str(train_abm).upper(),
			TMPDIR="/tmp",
			TMP="/tmp",
			TEMP="/tmp",
			PYTHONDONTWRITEBYTECODE="1",
			PYTHONPYCACHEPREFIX="/tmp/python-cache",
		)
		run_table_storage("--tables-restore", work, r_bin, env)
		new_args = []
		skip = False

		def map_path(value):
			if value.startswith("--") and "=" in value:
				key, path = value.split("=", 1)
				return key + "=" + map_path(path)
			return (
				str(work) + value[len(str(root)) :]
				if value == str(root) or value.startswith(str(root) + "/")
				else value
			)

		for i, value in enumerate(argv):
			if skip:
				skip = False
				continue
			if value == "--analysis-root":
				new_args += [value, str(work)]
				skip = True
			elif value.startswith("--analysis-root="):
				new_args.append("--analysis-root=" + str(work))
			else:
				new_args.append(map_path(value))
		if not any(x == "--analysis-root" or x.startswith("--analysis-root=") for x in new_args):
			new_args += ["--analysis-root", str(work)]
		if launch_shiny:
			new_args.append("--prepare-only")
		command = [sys.executable, str(HERE / "0.common.py"), *new_args]
		subprocess.run(command, env=env, check=True)
		with analysis_lock(publication_lock, waiting="Waiting to publish completed results."):
			publish_workspace(work, root, initial, r_bin, env, scopes=scopes)
	except BaseException:
		print(
			f"[LE8] Run failed; published results were retained. Diagnostic workspace: {workspace}",
			file=sys.stderr,
		)
		raise
	else:
		if (work / "logs").is_dir():
			logs = Path("/tmp/le8-logs") / workspace.name
			logs.parent.mkdir(exist_ok=True)
			shutil.move(str(work / "logs"), str(logs))
			print(f"[LE8] Logs: {logs}", flush=True)
		shutil.rmtree(workspace)
	finally:
		locks.close()
	if launch_shiny:
		shiny_env = dict(
			os.environ,
			LE8_ANALYSIS_ROOT=str(root),
			TMPDIR="/tmp",
			LE8_SHINY_PORT=option("--port", "3839"),
			LE8_SHINY_HOST=option("--host", "127.0.0.1"),
			LE8_SHINY_MAX_TABLE_MB=option("--max-table-mb", "128"),
		)
		command = [r_bin, str(ROOT / "shiny/app.R")]
		if "--shiny-review" in argv:
			command.append("--review")
		viewer_key = hashlib.sha256((str(root) + shiny_env["LE8_SHINY_PORT"]).encode()).hexdigest()[:20]
		viewer_lock = open(Path("/tmp") / f"le8-shiny-{viewer_key}.lock", "a")
		try:
			fcntl.flock(viewer_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
		except BlockingIOError:
			print(
				"[LE8] Shiny is already serving this result directory; the updated index is available on reload.",
				flush=True,
			)
		else:
			subprocess.run(command, env=shiny_env, check=True)
		finally:
			viewer_lock.close()
	return True


def publish_workspace(work, root, initial, r_bin, env, scopes=None):
	# Report paths refer to the published location, not an expired scratch folder.
	for path in work.rglob("*"):
		if path.is_file() and path.suffix in {".json", ".md", ".html", ".log"}:
			text = path.read_text(errors="strict")
			if str(work) in text:
				path.write_text(text.replace(str(work), str(root)))
	# Workbooks contain reviewable results, including IDs; RDS is reserved for reusable raw/model objects.
	# Rendering and Shiny used the scratch CSV exchanges.
	print("[LE8] Consolidating result tables into workbooks...", flush=True)
	run_table_storage("--tables-pack", work, r_bin, env)
	published = set()
	updates = []
	for path in work.rglob("*"):
		if not path.is_file() or any(
			x
			in {
				"_history",
				"_source_figures",
				"_previous",
				"__pycache__",
				".ruff_cache",
				".pgs_focus_cache",
				"prepared",
				"neural",
				"quality_neural",
				"checkpoints",
				"le8_annotations",
				"logs",
			}
			for x in path.relative_to(work).parts
		):
			continue
		relative = path.relative_to(work)
		if not owns_scope(relative, scopes):
			continue
		# Workbooks preserve scratch table bytes for reuse and source verification.
		if path.name.endswith((".csv", ".csv.gz", ".tsv", ".tsv.gz", ".tmp", ".pyc", ".log")):
			continue
		published.add(str(relative))
		stat = path.stat()
		if initial.get(str(relative)) == (stat.st_size, stat.st_mtime_ns):
			continue
		target = root / relative
		updates.append((path, target, relative))
	obsolete = [
		root / relative
		for relative in set(initial) - published
		if (root / relative).is_file()
		and relative.endswith(
			(
				".rds",
				".csv",
				".csv.gz",
				".tsv",
				".tsv.gz",
				".xlsx",
				".png",
				".pdf",
				".json",
				".html",
				".md",
			)
		)
	]
	previous = work.parent / "publish-before"
	previous.mkdir()
	for target in [target for _, target, _ in updates] + obsolete:
		if target.is_file():
			backup = previous / target.relative_to(root)
			backup.parent.mkdir(parents=True, exist_ok=True)
			shutil.copy2(target, backup)
	touched = []
	try:
		for source, target, relative in updates:
			target.parent.mkdir(parents=True, exist_ok=True)
			touched.append(target)
			shutil.copy2(source, target)
		for target in obsolete:
			touched.append(target)
			target.unlink()
	except BaseException:
		for target in reversed(touched):
			backup = previous / target.relative_to(root)
			if backup.is_file():
				shutil.copy2(backup, target)
			else:
				target.unlink(missing_ok=True)
		raise
	for directory in sorted(
		(p for p in root.rglob("*") if p.is_dir() and owns_scope(p.relative_to(root), scopes)),
		key=lambda p: len(p.parts),
		reverse=True,
	):
		if not any(directory.iterdir()):
			directory.rmdir()
	print(
		"[LE8] Published consolidated workbooks, figures and reusable RDS results.",
		flush=True,
	)


# 🚩 Portable Shiny package


def share_viewer(root, output, r_bin):
	from final import share_viewer as export

	return export(root, output, r_bin)


def main(argv=None):
	argv = list(sys.argv[1:] if argv is None else argv)
	if argv[:1] == ["pipeline"]:
		return pipeline_main(argv[1:])
	if argv[:1] == ["run-state"]:
		return state_main(argv[1:])
	public = argv[1:] if argv[:1] == ["dispatch"] else argv
	if (
		public
		and not public[0].startswith("-")
		and set(canonical_modules(public[0])) & {"c1_abm", "shiny", "share"}
	):
		return subprocess.run([sys.executable, HERE / "0.common.py", "pipeline", *public], check=True)
	if (
		os.getenv("LE8_TABLE_WORKSPACE") != "1"
		and not any(v in argv for v in ("--help", "-h", "--dry-run", "--preflight"))
		and argv[:1] != ["index"]
	):
		args = argv[1:] if argv[:1] == ["dispatch"] else argv
		parsed, unknown = parser().parse_known_args(args)
		mods = canonical_modules(parsed.modules)
		if any(v not in ORDER for v in mods):
			raise ValueError("Unknown module; see --help.")
		if parsed.run_abm is False:
			mods = [v for v in mods if v != "c1_abm"]
		if parsed.run_abm and "c1_abm" not in mods:
			mods.append("c1_abm")
		mods = [v for v in ORDER if v in mods]
		if len(mods) > 1 and not set(mods) & {"final", "shiny", "share"}:
			options = args[1:] if args and not args[0].startswith("-") else args
			options = [v for v in options if v not in ("--run-abm", "--skip-abm")]
			abm_options = []
			actions = parser()._option_string_actions
			i = 0
			while i < len(options):
				token = options[i]
				key = token.split("=", 1)[0]
				action = actions.get(key)
				if action is None:
					i += 1
					continue
				abm_options.append(token)
				i += 1
				if "=" not in token and action.nargs != 0:
					count = 1 if action.nargs is None else action.nargs
					if not isinstance(count, int):
						raise ValueError("Unsupported variable-length public option: " + key)
					abm_options.extend(options[i : i + count])
					i += count
			if "c1_abm" in mods:
				dispatch_main(["c1_abm", *abm_options, "--preflight"])
			if "c5_cellulation" in mods:
				dispatch_main(["c5_cellulation", *abm_options, "--preflight"])
			previous_state = os.environ.get("LE8_RUN_STATE_DIR")
			with tempfile.TemporaryDirectory(prefix="le8-run-state-", dir="/tmp") as state:
				os.environ["LE8_RUN_STATE_DIR"] = state
				try:
					for module in mods:
						print(
							"[LE8] Module transaction: "
							+ module
							+ "; successful results publish before the next module",
							flush=True,
						)
						main(["dispatch", module, *(abm_options if module == "c1_abm" else options)])
				finally:
					if previous_state is None:
						os.environ.pop("LE8_RUN_STATE_DIR", None)
					else:
						os.environ["LE8_RUN_STATE_DIR"] = previous_state
			return
	if _table_runtime(argv):
		return
	if argv and argv[0] == "index":
		return index_main(argv[1:])
	if argv and argv[0] == "dispatch":
		argv = argv[1:]
	return dispatch_main(argv)


import runpy

pipeline_ROOT = Path(__file__).resolve().parents[1]
pipeline_HERE = pipeline_ROOT / "f"


def pipeline_parse(argv):
	options = argparse.ArgumentParser(add_help=False)
	options.add_argument("--abm-config", default=str(pipeline_ROOT / "config/abm.json"))
	options.add_argument("--deviation-config")
	options.add_argument("--objective-file")
	options.add_argument("--pgs-root")
	options.add_argument("--pgs-mapping")
	options.add_argument("--remerge-pgs", action="store_true")
	options.add_argument("--abm-host-args", default="")
	research, rest = options.parse_known_args(argv)
	common = globals()
	args, unknown = common["parser"]().parse_known_args(rest)
	modules = common["canonical_modules"](args.modules)
	if set(modules) - set(common["ORDER"]):
		raise ValueError("Unknown module: " + ",".join(sorted(set(modules) - set(common["ORDER"]))))
	if args.run_abm:
		modules.append("c1_abm")
	if args.run_abm is False:
		modules = [m for m in modules if m != "c1_abm"]
	modules = [m for m in common["ORDER"] if m in modules]
	if not args.Y or any(not re.fullmatch(r"[A-Za-z0-9_]+", y) for y in args.Y.split(",")):
		raise ValueError("Invalid outcome names")
	if not set(args.biom.split(",")) <= {"prot", "met"}:
		raise ValueError("Invalid omic layer")
	if args.cores is not None and args.cores < 1:
		raise ValueError("--cores must be positive")
	if args.host not in {"127.0.0.1", "localhost", "::1"} and os.getenv("LE8_SHINY_ALLOW_NETWORK") != "TRUE":
		raise ValueError("Network binding requires LE8_SHINY_ALLOW_NETWORK=TRUE")
	if not 1024 <= args.port <= 65535 or args.max_table_mb < 1:
		raise ValueError("Invalid port or table limit")
	# The shell supplies a positional module list. Direct Python use also works.
	tail = rest[1:] if rest and not rest[0].startswith("-") else rest
	tail = [v for v in tail if v not in ("--run-abm", "--skip-abm")]
	return research, args, modules, tail, common, unknown


def pipeline_worker_command(research, a, trait, layer):
	cmd = [
		os.getenv("ABM_PYTHON") or sys.executable,
		str(pipeline_HERE / "0.host.py"),
		"--Y",
		trait,
		"--biom",
		layer,
		"--analysis-root",
		str(a.analysis_root),
		"--config",
		research.abm_config,
		"--seed",
		str(a.seed),
		"--cores",
		str(a.cores or 4),
		"--r-bin",
		a.r_bin,
		"--replace",
		a.replace,
		"--host-args=" + research.abm_host_args,
	]
	for key in [
		"ukb_phe",
		"group_file",
		"group_col",
		"outer_roster",
		"end_date",
		"diagnosis_col",
		"shared_covariates",
	]:
		value = getattr(a, key)
		if value is not None:
			cmd += ["--" + key.replace("_", "-"), str(value)]
	for key in ["deviation_config", "objective_file", "pgs_root", "pgs_mapping"]:
		value = getattr(research, key)
		if value:
			cmd += ["--" + key.replace("_", "-"), value.replace("{biom}", layer)]
	if research.remerge_pgs:
		cmd += ["--remerge-pgs"]
	return cmd


def pipeline_main(argv=None):
	research, a, modules, tail, common, unknown = pipeline_parse(list(sys.argv[1:] if argv is None else argv))

	def call(cmd, dry=None):
		cmd = list(map(str, cmd))
		print("[LE8 research] " + shlex.join(cmd), flush=True)
		if not (a.dry_run if dry is None else dry):
			try:
				subprocess.run(cmd, check=True)
			except subprocess.CalledProcessError as exc:
				# systemd-run can report SIGKILL as -9, which a Python wrapper
				# exposes as 247. Keep a useful cause ahead of the long command.
				if exc.returncode in (-9, 137, 247):
					raise RuntimeError(
						f"Child process exited with code {exc.returncode} (possible SIGKILL/OOM); "
						"check the memory cap and kernel log. "
						"Earlier published results are retained; resume without --replace TRUE."
					) from exc
				raise

	if any(
		t.startswith(
			(
				"--legacy",
				"--abm-backend",
				"--abm-engine",
				"--abm-args",
				"--s7-",
				"--abm-tf-root",
				"--abm-root",
			)
		)
		for t in unknown
	):
		raise ValueError("Obsolete ABM options are unsupported; use --abm-config and --abm-host-args")
	if unknown and not (
		set(modules) & common["NATIVE"] or a.fit_joint or a.fit_reference or a.fit_genetic or a.details
	):
		raise ValueError("Unrecognized options: " + " ".join(unknown))
	# Fail before a native C1 scan if research configuration or files are unavailable.
	if "c1_abm" in modules:
		for trait in a.Y.split(","):
			for layer in a.biom.split(","):
				call([*pipeline_worker_command(research, a, trait, layer), "--preflight"])
	for module in modules:
		if module == "share":
			if modules != ["share"]:
				raise ValueError("Run share separately after preparing final/shiny")
			if a.preflight:
				call(
					[
						a.r_bin,
						"-e",
						"p<-c('shiny','bslib','DT','DBI','RSQLite','ggplot2','plotly','openxlsx','jsonlite');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste(m,collapse=', '))",
					]
				)
			else:
				call(
					[
						sys.executable,
						pipeline_HERE / "final.py",
						"share",
						"--analysis-root",
						a.analysis_root,
						"--out",
						a.share_out or a.analysis_root.parent / "le8-share.zip",
						"--r-bin",
						a.r_bin,
						"--max-table-mb",
						a.max_table_mb,
					]
				)
		elif module == "c1_abm":
			if a.preflight:
				continue
			for trait in a.Y.split(","):
				for layer in a.biom.split(","):
					worker = pipeline_worker_command(research, a, trait, layer)
					budget = a.memory_limit_gb if a.memory_limit_gb is not None else 32
					if budget < 0:
						raise ValueError("Invalid memory budget")
					# Same cgroup guard and CPU lease as the existing dispatcher.
					if budget:
						if os.getenv("LE8_MEMORY_SCOPE_ACTIVE") != "1":
							worker = [
								"systemd-run",
								"--user",
								"--scope",
								"--quiet",
								"--collect",
								"-p",
								f"MemoryMax={budget}G",
								"-p",
								f"MemorySwapMax={a.memory_swap_gb if a.memory_swap_gb is not None else 4}G",
								"-p",
								"OOMPolicy=kill",
								"--",
								"env",
								"LE8_MEMORY_SCOPE_ACTIVE=1",
								*worker,
							]
						worker = [
							sys.executable,
							pipeline_ROOT / "f/0.resources.py",
							"phase",
							"--label",
							"c1_abm",
							"--cores",
							str(a.cores or 4),
							"--memory-gib",
							str(budget),
							"--",
							*worker,
						]
					call(worker)
		elif module in {"final", "shiny"}:
			if module == "final":
				call([sys.executable, pipeline_ROOT / "f/0.common.py", "dispatch", "final", *tail])
			if not a.preflight and (module == "final" or not a.no_reindex):
				call(
					[
						sys.executable,
						pipeline_HERE / "final.py",
						"evidence",
						"--analysis-root",
						a.analysis_root,
						"--Y",
						a.Y,
						"--biom",
						a.biom,
						"--max-table-mb",
						a.max_table_mb,
					]
				)
			if module == "shiny":
				if a.preflight:
					call(
						[
							a.r_bin,
							"-e",
							"p<-c('shiny','bslib','DT','DBI','RSQLite','ggplot2','plotly');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste(m,collapse=', '))",
						]
					)
				elif not a.prepare_only:
					call(
						[
							a.r_bin,
							pipeline_ROOT / "shiny/app.R",
							"--db",
							a.analysis_root / "evidence/evidence.sqlite",
							"--port",
							a.port,
							"--host",
							a.host,
							*(["--review"] if a.shiny_review else []),
						]
					)
		else:
			call([sys.executable, pipeline_ROOT / "f/0.common.py", "dispatch", module, *tail])


state_RESULTS = {
	"c1_correlate": ["c1_correlate/c1.res.rds"],
	"c1_prediction": [
		"c1_correlate/prediction/res.rds",
		"c1_correlate/prediction/_LOCAL_ONLY/benchmark.json.gz",
	],
	"c2_cause": ["c2_cause/c2.res.rds"],
	"c3_coloc": ["c3_coloc/c3.res.rds"],
	"c4_connect": [
		"c4_connect/" + name
		for name in ("c4.res.rds", "c4.interactions.res.rds", "c4.nonlin.res.rds", "c4.penalty.res.rds")
	],
	"c4_panel_validation": ["c4_connect/c4.validation.rds"],
	"c4_explain": ["c4_connect/c4.explain.rds"],
}


def state_digest(path):
	with path.open("rb") as stream:
		return hashlib.file_digest(stream, "sha256").hexdigest()


def state_snapshot(root, trait, layer, module):
	paths = [root / trait / layer / name for name in state_RESULTS[module]]
	if any(not p.is_file() or not p.stat().st_size for p in paths):
		raise ValueError(f"Missing result for {trait}/{layer}/{module}")
	return {str(p.relative_to(root)): state_digest(p) for p in paths}


def state_receipts(action, root, traits, layers, modules):
	directory = os.getenv("LE8_RUN_STATE_DIR")
	if not directory:
		return []
	directory = Path(directory)
	identity = str(Path(os.getenv("LE8_PUBLISHED_ROOT", root)).resolve())
	found = []
	for trait in traits:
		for layer in layers:
			for module in modules:
				if module not in state_RESULTS:
					continue
				key = f"{trait}|{module}|{layer}"
				file = directory / (hashlib.sha256((identity + key).encode()).hexdigest() + ".json")
				if action == "read" and not file.is_file():
					continue
				try:
					current = dict(root=identity, key=key, outputs=state_snapshot(root, trait, layer, module))
					if action == "record":
						temporary = file.with_suffix(".tmp")
						temporary.write_text(json.dumps(current, sort_keys=True))
						temporary.replace(file)
					elif json.loads(file.read_text()) == current:
						found.append(key)
				except (ValueError, OSError):
					if action == "record":
						raise
					# Invalid receipts require ordinary output-completeness checks.
	return found


def state_main(argv=None):
	p = argparse.ArgumentParser(description=__doc__)
	p.add_argument("action", choices=["read", "record"])
	p.add_argument("root", type=Path)
	p.add_argument("traits")
	p.add_argument("layers")
	p.add_argument("modules")
	a = p.parse_args(argv)
	for key in state_receipts(
		a.action, a.root, a.traits.split(","), a.layers.split(","), a.modules.split(",")
	):
		print(key)


# 🚩 Shared table and figure I/O
import io, zipfile
from xml.etree import ElementTree as ET

NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
PRIVATE = re.compile(
	r"(^|[._ #])(eid|iid|fid|participant_id|patient_id|donor_id|query_id|subject_id)($|[._ ])", re.I
)


def sha256(path):
	h = hashlib.sha256()
	with open(path, "rb") as f:
		for b in iter(lambda: f.read(1024 * 1024), b""):
			h.update(b)
	return h.hexdigest()


def read_table(path, **kwargs):
	p = Path(path)
	if p.suffix.lower() == ".rds":
		import pyreadr

		out = pyreadr.read_r(str(p))
		if len(out) != 1:
			raise ValueError("RDS must contain exactly one data.frame")
		return next(iter(out.values()))
	if p.suffix.lower() == ".parquet":
		return pd.read_parquet(p, **kwargs)
	if p.suffix.lower() == ".xlsx":
		sheet = kwargs.pop("sheet", None)
		frames = read_xlsx(p)
		if sheet is None and len(frames) != 1:
			raise ValueError("Select an explicit XLSX sheet")
		return frames[sheet] if sheet else next(iter(frames.values()))
	return pd.read_csv(p, sep="," if p.name.endswith((".csv", ".csv.gz")) else "\t", **kwargs)


def id_table(path, id_col="eid"):
	# Read IDs as text so leading zeros cannot be silently lost in delimited inputs.
	d = (
		read_table(path, dtype={id_col: str})
		if not str(path).endswith((".rds", ".parquet", ".xlsx"))
		else read_table(path)
	)
	if id_col not in d:
		raise ValueError(f"{path}: missing ID {id_col}")
	if d[id_col].isna().any():
		raise ValueError("Missing IDs")

	if pd.api.types.is_float_dtype(d[id_col]) and np.equal(d[id_col], np.floor(d[id_col])).all():
		d[id_col] = d[id_col].astype("int64")
	d[id_col] = d[id_col].astype(str)
	if d[id_col].duplicated().any() or d[id_col].eq("").any():
		raise ValueError("Duplicate/empty IDs")
	return d.set_index(id_col, drop=False)


def win_to_local(path):
	s = str(path)
	if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", s):
		s = "/mnt/" + s[0].lower() + "/" + s[3:].replace("\\", "/")
	return Path(s)


def score_source(root, layer):
	p = win_to_local(root)
	if p.is_file():
		return p
	path = p / f"{layer}.pgs.gz"
	if path.exists():
		return path
	candidates = list(p.glob("*.pgs.gz"))
	if len(candidates) == 1:
		return candidates[0]
	raise ValueError(
		f"Cannot unambiguously find {layer}.pgs.gz in {p}; run format.sh pgs step2 or supply an exact file"
	)


def load_pgs(root, layer, ids, measured_features, mapping=None):
	"""Reads format.sh's merged eid + <assay>.SCORE_SUM; row-order-independent.

    A per-trait #IID/ALLELE_CT/SCORE_SUM file is also accepted. SCORE_AVG is NOT
    silently substituted. Constant/header-only scores are unavailable, never zeros.
    A wide merged score cannot certify individual variant-pattern QC.
    """
	path = score_source(root, layer)
	opener = gzip.open if path.name.endswith(".gz") else open
	with opener(path, "rt", encoding="utf-8") as h:
		head = next(csv.reader(h, delimiter="\t"))
	if len(set(head)) != len(head):
		raise ValueError("Duplicate raw PGS column names")
	id_candidates = [c for c in ("eid", "#IID", "IID") if c in head]
	if len(id_candidates) != 1:
		raise ValueError("PGS needs exactly one eid/#IID/IID column")
	id_col = id_candidates[0]
	amap = {}
	for c in head:
		if c.endswith(".SCORE_SUM"):
			amap[c] = c[: -len(".SCORE_SUM")]
	if mapping is not None:
		mp = read_table(mapping)
		if not {"pgs_column", "feature"} <= set(mp):
			raise ValueError("Mapping needs pgs_column,feature")
		if mp.pgs_column.duplicated().any() or mp.feature.duplicated().any():
			raise ValueError("Ambiguous mapping")
		amap.update(dict(zip(mp.pgs_column, mp.feature)))
	amap = {c: f for c, f in amap.items() if c in head and f in measured_features}
	if len(set(amap.values())) != len(amap):
		raise ValueError("Multiple PGSs map to one assay; use an explicit source-specific run")
	if not amap:
		raise ValueError("No exact PGS/assay matches; do not use fuzzy gene/assay matching")
	features = [f for f in measured_features if f in set(amap.values())]
	cols = [next(c for c, f2 in amap.items() if f2 == f) for f in features]
	# Keep only matched columns, and align by ID, not position.
	requested = set(map(str, ids))
	chunks = []
	seen = set()
	for chunk in pd.read_csv(path, sep="\t", usecols=[id_col] + cols, dtype={id_col: str}, chunksize=20000):
		current = chunk[id_col]
		if current.isna().any() or current.duplicated().any() or set(current) & seen:
			raise ValueError("PGS duplicate/missing IDs")
		seen.update(current)
		chunks.append(chunk[current.isin(requested)])
	d = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=[id_col] + cols)
	d = d.set_index(id_col).reindex(list(map(str, ids)))
	g = d[cols].apply(pd.to_numeric, errors="raise").to_numpy(float)
	if np.isinf(g).any():
		raise ValueError("Infinite PGS value")
	audit = []
	for j, f in enumerate(features):
		x = g[:, j]
		v = x[np.isfinite(x)]
		u = len(np.unique(v))
		audit.append(
			dict(
				feature=f,
				pgs_column=cols[j],
				source=str(path),
				source_sha256=sha256(path) if j == 0 else "",
				matched_N=len(v),
				unique_scores=u,
				rounded_3dp_fraction=(
					float(np.mean(np.abs(v * 1000 - np.round(v * 1000)) < 1e-7)) if len(v) else np.nan
				),
				status="available" if u > 1 else "constant_or_unavailable",
				variant_pattern_qc="not_inferred_from_scores",
			)
		)
	return features, g, pd.DataFrame(audit)


def merge_native_scores(manifest, output):
	"""Re-merge existing per-trait scores at full precision; never regenerate GWAS.
    Fails on missing people rather than inventing zero scores. Output is create-only.
    """
	manifest = win_to_local(manifest)
	output = win_to_local(output)
	if output.exists():
		raise FileExistsError(output)
	m = read_table(manifest)
	if not {"GWAS", "STATUS", "PGS"} <= set(m):
		raise ValueError("Need format.sh <layer>.pgs.files.tsv manifest")
	incomplete = m[~m.STATUS.isin(["included", "empty_no_matched_variants"])]
	if len(incomplete):
		raise ValueError("Manifest contains incomplete traits")
	merged = None
	seen = set()
	for r in m.itertuples(index=False):
		if r.STATUS != "included":
			continue
		if r.GWAS in seen:
			raise ValueError("Duplicate trait in PGS manifest")
		seen.add(r.GWAS)
		p = win_to_local(r.PGS)
		if not p.is_absolute():
			p = manifest.parent / p
		d = pd.read_csv(p, sep="\t", dtype={"#IID": str})
		col = f"{r.GWAS}.SCORE_SUM"
		if not {"#IID", col} <= set(d) or d.empty or d["#IID"].isna().any() or d["#IID"].duplicated().any():
			raise ValueError(f"Invalid native score: {p}")
		s = pd.to_numeric(d.set_index("#IID")[col], errors="raise")
		if not np.isfinite(s).all():
			raise ValueError("Nonfinite score")
		if merged is None:
			merged = s.to_frame()
		else:
			if set(merged.index) != set(s.index):
				raise ValueError("Different PGS participant sets")
			merged[col] = s.reindex(merged.index)
	if merged is None:
		raise ValueError("No nonempty PGSs")
	output.parent.mkdir(parents=True, exist_ok=True)
	merged.index.name = "eid"
	merged.to_csv(output, sep="\t", float_format="%.17g", compression="infer")
	return dict(
		output=str(output),
		sha256=sha256(output),
		assays=merged.shape[1],
		N=len(merged),
		precision="17 significant digits; cannot recover previously rounded native inputs",
	)


def _letter(n):
	s = ""
	while n:
		n, r = divmod(n - 1, 26)
		s = chr(65 + r) + s
	return s


def write_xlsx(path, tables):
	"""Minimal standards-compliant OOXML writer for the user's portable code.
    All strings are inline strings, not formulas. Clinical outputs remain numeric.
    """
	path = Path(path)
	path.parent.mkdir(parents=True, exist_ok=True)
	from xml.sax.saxutils import escape, quoteattr

	ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
	rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
	sheets = []
	payload = {}
	used = set()
	for i, (name, df) in enumerate(tables.items(), 1):
		name = re.sub(r"[\[\]:*?/\\]", "_", str(name))[:31] or f"Sheet{i}"
		if name in used:
			name = name[:25] + "_" + str(i)
		used.add(name)
		if len(df) > 1048575:
			raise ValueError("Sheet exceeds Excel row limit; use table store instead")
		rows = []
		for n, record in enumerate(
			[list(map(str, df.columns))] + list(df.itertuples(index=False, name=None)), 1
		):
			cells = []
			for j, v in enumerate(record, 1):
				ref = _letter(j) + str(n)
				style = ' s="1"' if n == 1 else ""
				if v is None or (isinstance(v, (float, np.floating)) and not np.isfinite(v)):
					continue
				if isinstance(v, (bool, np.bool_)):
					c = f'<c r="{ref}" t="b"{style}><v>{int(v)}</v></c>'
				elif isinstance(v, (int, float, np.number)):
					c = f'<c r="{ref}"{style}><v>{v}</v></c>'
				else:
					c = f'<c r="{ref}" t="inlineStr"{style}><is><t xml:space="preserve">{escape(str(v))}</t></is></c>'
				cells.append(c)
			rows.append(f'<row r="{n}">' + "".join(cells) + "</row>")
		cols = f'<cols><col min="1" max="{max(1,len(df.columns))}" width="20" customWidth="1"/></cols>'
		payload[f"xl/worksheets/sheet{i}.xml"] = (
			f'<worksheet xmlns="{ns}"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" state="frozen"/></sheetView></sheetViews>{cols}<sheetData>'
			+ "".join(rows)
			+ "</sheetData></worksheet>"
		)
		sheets.append(f'<sheet name={quoteattr(name)} sheetId="{i}" r:id="rId{i}"/>')
	payload["xl/workbook.xml"] = (
		f'<workbook xmlns="{ns}" xmlns:r="{rel}"><sheets>' + "".join(sheets) + "</sheets></workbook>"
	)
	rns = "http://schemas.openxmlformats.org/package/2006/relationships"
	payload["xl/_rels/workbook.xml.rels"] = (
		f'<Relationships xmlns="{rns}">'
		+ "".join(
			f'<Relationship Id="rId{i}" Type="{rel}/worksheet" Target="worksheets/sheet{i}.xml"/>'
			for i in range(1, len(sheets) + 1)
		)
		+ f'<Relationship Id="style" Type="{rel}/styles" Target="styles.xml"/></Relationships>'
	)
	payload["_rels/.rels"] = (
		f'<Relationships xmlns="{rns}"><Relationship Id="rId1" Type="{rel}/officeDocument" Target="xl/workbook.xml"/></Relationships>'
	)
	payload["xl/styles.xml"] = (
		f'<styleSheet xmlns="{ns}"><fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts><fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills><borders count="1"><border/></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs></styleSheet>'
	)
	types = [
		'<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
		'<Default Extension="xml" ContentType="application/xml"/>',
		'<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
		'<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
	]
	types += [
		f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
		for i in range(1, len(sheets) + 1)
	]
	payload["[Content_Types].xml"] = (
		'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
		+ "".join(types)
		+ "</Types>"
	)
	with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
		for n, v in payload.items():
			z.writestr(n, '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' + v)


def read_xlsx(path, max_bytes=256 * 1024**2):
	"""Read ordinary XLSX worksheet VALUES without executing formulas/macros.
    Raw embedded exports are indexed separately by evidence for full precision.
    """
	with zipfile.ZipFile(path) as z:
		if sum(x.file_size for x in z.infolist()) > max_bytes:
			raise ValueError("Workbook exceeds expanded size cap")
		shared = []
		if "xl/sharedStrings.xml" in z.namelist():
			root = ET.fromstring(z.read("xl/sharedStrings.xml"))
			shared = ["".join(e.itertext()) for e in root.findall("s:si", NS)]
		rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
		rm = {e.attrib["Id"]: e.attrib["Target"] for e in rels}
		wb = ET.fromstring(z.read("xl/workbook.xml"))
		out = {}
		for sh in wb.findall("s:sheets/s:sheet", NS):
			rid = sh.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
			target = rm[rid]
			target = target.lstrip("/") if target.startswith("/") else "xl/" + target
			if ".." in Path(target).parts:
				raise ValueError("Unsafe worksheet path")
			root = ET.fromstring(z.read(target))
			rows = []
			for row in root.findall("s:sheetData/s:row", NS):
				vals = []
				for c in row.findall("s:c", NS):
					ref = re.match(r"[A-Z]+", c.attrib["r"]).group()
					j = 0
					for letter in ref:
						j = j * 26 + ord(letter) - 64
					while len(vals) < j:
						vals.append(None)
					typ = c.attrib.get("t")
					v = c.find("s:v", NS)
					if typ == "inlineStr":
						value = "".join(c.find("s:is", NS).itertext())
					elif v is None:
						value = None
					elif typ == "s":
						value = shared[int(v.text)]
					elif typ == "b":
						value = v.text == "1"
					elif typ in ("str", "e"):
						value = v.text
					else:
						try:
							value = float(v.text)
						except (TypeError, ValueError):
							value = v.text
					vals[j - 1] = value
				rows.append(vals)
			if not rows:
				out[sh.attrib["name"]] = pd.DataFrame()
				continue
			width = max(map(len, rows))
			header = rows[0] + [None] * (width - len(rows[0]))
			seen = {}
			names = []
			for j, h in enumerate(header):
				h = str(h) if h is not None else f"column_{j+1}"
				seen[h] = seen.get(h, 0) + 1
				names.append(h if seen[h] == 1 else f"{h}__{seen[h]}")
			out[sh.attrib["name"]] = pd.DataFrame(
				[r + [None] * (width - len(r)) for r in rows[1:]], columns=names
			)
		return out


def write_tables(out, tables, module, layer, trait="cvd_cad", run_id="", workbook="results.xlsx"):
	"""One workbook + machine-readable aggregate TSVs; no participant identifiers."""
	out = Path(out)
	out.mkdir(parents=True, exist_ok=True)
	registry = []
	safe = {}
	for name, d in tables.items():
		d = pd.DataFrame(d)
		if not len(d.columns):
			d = pd.DataFrame(columns=["status"])
		if any(PRIVATE.search(c) for c in d.columns):
			raise ValueError(f"Private column in public table {name}")
		stem = re.sub("[^a-zA-Z0-9_.-]", "_", name)
		p = out / (stem + ".tsv.gz")
		d.to_csv(p, sep="\t", index=False)
		registry.append(
			dict(
				name=name,
				path=p.name,
				module=module,
				layer=layer,
				trait=trait,
				run_id=run_id,
				rows=len(d),
				sha256=sha256(p),
				status="available" if len(d) else "empty",
			)
		)
		if len(d) <= 200000 and len(d) * len(d.columns) <= 2000000:
			safe[name] = d
	if safe:
		write_xlsx(out / workbook, safe)
	(out / "table_manifest.json").write_text(json.dumps(registry, ensure_ascii=False, indent=2))
	return registry


# 🚩 Aggregate figure output
def save_figure(fig, path):
	fig.tight_layout()
	fig.savefig(path, dpi=160, bbox_inches="tight")
	import matplotlib.pyplot as plt

	plt.close(fig)


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		sys.exit(130)
	except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
		print("ERROR: " + str(exc), file=sys.stderr)
		sys.exit(2)
