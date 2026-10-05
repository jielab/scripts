#!/usr/bin/env python3


# 🚩 0.common
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
	# S7 aggregate tables identify the comparator model with reference_id.
	if {'run_id', 'model_id', 'primary_id', 'reference_id', 'n', 'uno_c_horizon'} <= set(names):
		names = [name for name in names if name != 'reference_id']
	private = {re.sub(r"[^a-z0-9]", "", name) for name in PRIVATE_COLUMNS}
	private.update({"id", "id1", "id2", "sampleid", "ukbid", "ukbeid"})
	return any(
		re.sub(r"[^a-z0-9]", "", name) in private
		or re.search(r"(^|[_. #])(eid|iid|fid)($|[_. ])", name)
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
	if module == "final_prediction":
		return final_dir(root, trait, layer, "prediction") / re.sub(
			r"^final[._]", "", filename
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
	if "c5_cellulation" in parts or re.search(
		r"^c5\.(?:(?:systematic[.])?cell[._]|CIGMA)", path.name
	):
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
					modified_utc=datetime.fromtimestamp(
						p.stat().st_mtime, timezone.utc
					).isoformat(),
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
			self.errors.append(
				dict(Y=Y, layer=layer, role=role, status="blocked", detail=str(p))
			)
			return pd.DataFrame()
		if p.stat().st_size > self.max_bytes:
			self.errors.append(
				dict(Y=Y, layer=layer, role=role, status="too_large", detail=str(p))
			)
			return pd.DataFrame()
		sep = "\t" if ".tsv" in p.name else ","
		try:
			head = pd.read_csv(p, sep=sep, nrows=0)
			if has_private_columns(head.columns):
				raise ValueError(
					"Participant-level identifier column; not exposed in Shiny"
				)
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
				modified_utc=datetime.fromtimestamp(
					p.stat().st_mtime, timezone.utc
				).isoformat(),
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
		if module == "final_prediction":
			return [result_file(self.root, Y, layer, module, name)]
		return [result_file(self.root, Y, layer, module, name)]

	def metadata(self, d, Y, layer, feature=None):
		z = pd.DataFrame(index=d.index)
		z["Y"] = Y
		z["layer"] = layer
		z["feature"] = (
			feature if feature is not None else series(d, FEATURE_COLUMNS, "")
		).map(normalize_feature)
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
					for _, ix in master.groupby(
						["scope", "adjustment"], dropna=False
					).groups.items():
						master.loc[ix, qcol] = bh(master.loc[ix, pcol])
			master["same_people"] = series(paired, ["comparison"], "").str.contains(
				"identical people", case=False, na=False
			)
			master["comparison"] = series(paired, ["comparison"], "Unverified matching")
			master["genetic_evidence"] = (
				"PGS association; discovery overlap must be checked separately"
			)
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
		if len(anchors) and {"feature", "component", "beta", "std.error"} <= set(
			anchors
		):
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
			+ [
				final_dir(self.root, Y) / f"Fig1.PGS_integrated.{layer}_calibration.csv"
			],
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
				units = series(
					t, ["effect_measure"], "Log HR per SD; landmark-specific risk set"
				)
				if role == "event_window" and "side" in t:
					prev = t.side.astype(str).str.contains(
						"Pre-baseline", case=False, na=False
					)
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
			zz["PP_H4_robust_min"] = numeric(
				co, ["PP.H4_robust_min", "PP.H4_p12_conservative"]
			)
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
			cm = zz.loc[
				ix, ["feature", "PP_H4", "PP_H4_robust_min", "locus", "credible_set_n"]
			].rename(columns={"locus": "best_reported_coloc_region"})
		mr = self.read(
			self.paths(Y, layer, "c2_cause", "c2.MR_all.csv"), Y, layer, "MR"
		)
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
			g = g[
				["feature", "scope", "beta_difference", "p", "FDR", "_source_file"]
			].rename(
				columns={
					"beta_difference": "GR_beta_difference",
					"p": "GR_difference_p",
					"FDR": "GR_difference_FDR",
					"_source_file": "GR_source_file",
				}
			)
			master = master.merge(
				g, on=["feature", "scope"], how="left", validate="many_to_one"
			)
		if len(cm):
			master = master.merge(cm, on="feature", how="left", validate="many_to_one")
		if len(mm):
			master = master.merge(mm, on="feature", how="left", validate="many_to_one")
		for c in [
			"measured_p",
			"measured_FDR",
			"pgs_p",
			"pgs_FDR",
			"GR_beta_difference", "GR_difference_p", "GR_difference_FDR",
			"best_reported_cis_MR_FDR_all",
			"PP_H4_robust_min",
		]:
			if c not in master:
				master[c] = np.nan
		master["measured_P_ge_005_PGS_FDR_lt_005"] = (
			master.same_people.fillna(False)
			& master.measured_p.ge(0.05)
			& master.pgs_FDR.lt(0.05)
		)
		master["measured_FDR_ge_005_PGS_FDR_lt_005"] = (
			master.same_people.fillna(False)
			& master.measured_FDR.ge(0.05)
			& master.pgs_FDR.lt(0.05)
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
		master["causal_status"] = (
			"Not classified from observational/PGS significance alone"
		)
		master["candidate_id"] = [
			hashlib.sha256("|".join(map(str, v)).encode()).hexdigest()[:20]
			for v in master[
				["Y", "layer", "feature", "scope", "adjustment"]
			].itertuples(index=False, name=None)
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
			master["candidate_id"] = (
				master.candidate_id + "_" + master.index.astype(str)
			)
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
					measured_nominal_null_PGS_supported=int(
						g.measured_P_ge_005_PGS_FDR_lt_005.sum()
					),
					measured_FDR_null_PGS_supported=int(
						g.measured_FDR_ge_005_PGS_FDR_lt_005.sum()
					),
					measured_nominal_null_cis_MR_supported=int(
						g.measured_P_ge_005_cis_MR_FDR_lt_005.sum()
					),
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
			ci = ci[use].rename(
				columns={"N": "N_bootstrap", "delta_R2_vs_NS": "delta_R2_refit"}
			)
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
					if not PRIVATE_PARTS.search(name)
					and not (Y == "all" and Path(directory) == yd)
				]
				for filename in filenames:
					p = Path(directory) / filename
					if not re.search(
						r"\.(png|jpg|jpeg|pdf|csv|tsv)(\.gz)?$", filename, re.I
					):
						continue
					if not p.is_file() or not safe_path(p, self.root):
						continue
					rel = str(p.relative_to(self.root))
					parts = p.relative_to(yd).parts
					if (
						Y != "all"
						and parts[0] in {"prot", "met"}
						and parts[0] not in self.layers
					):
						continue
					layer = (
						"all"
						if Y == "all"
						else parts[0]
						if parts[0] in {"prot", "met"}
						else "joint"
					)
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
					elif str(p.resolve()) not in known and aggregate_allowed(
						p, self.root
					):
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
			pd.DataFrame(
				self.errors, columns=["Y", "layer", "role", "status", "detail"]
			),
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


ABM_BACKENDS = {
	"reference": ("reference", "v5_attention_allteachers"),
	"tf": ("tabicl", "tabiclv2_finetuned"),
	"tabicl": ("tabicl", "tabiclv2_finetuned"),
	"selective_attention": ("selective_attention", "selective_attention"),
}


def result_candidates(root, trait, layer, kind, external_root=None):
	root = Path(root)
	backend, legacy_run = ABM_BACKENDS[kind]
	base = root / trait / layer / "c1_correlate" / ("abm_" + backend)
	candidates = [base]
	if external_root is not None:
		candidates.insert(0, Path(external_root) / trait / layer / legacy_run)
	return candidates


ABM_TABLES = [
	"test_metrics.csv",
	"raw_test_metrics.csv",
	"approach_comparison.csv",
	"paired_contrasts.csv",
	"support_error_audit.csv",
	"subgroup_metrics.csv",
	"coverage_curve.csv",
	"masked_feature_metrics.csv",
	"masked_reconstruction_summary.csv",
	"learning_curve.csv",
	"embedding_diagnostics.csv",
]

ABM_FIGURES = [
	"Fig_coverage.png",
	"Fig_masked_reconstruction.png",
	"Fig_model_comparison.png",
]


def import_abm_results(
	root, traits, layers, reference_root=None, tf_root=None, render_figures=False
):
	# Existing local results are read in place. Only explicit external results are copied.
	for trait in traits:
		for layer in layers:
			for kind, external in [("reference", reference_root), ("tf", tf_root)]:
				backend = "tabicl" if kind == "tf" else kind
				target = root / trait / layer / "c1_correlate" / ("abm_" + backend)
				source = next(
					(
						p
						for p in result_candidates(root, trait, layer, kind, external)
						if any(
							(p / name).is_file() for name in ABM_TABLES + ABM_FIGURES
						)
					),
					None,
				)
				if source is None:
					continue
				if source.resolve() != target.resolve():
					target.mkdir(parents=True, exist_ok=True)
					for name in ABM_TABLES + ABM_FIGURES + ["cohort_audit.json"]:
						path = source / name
						if not path.is_file():
							continue
						if path.stat().st_size > 128 * 1024**2:
							raise ValueError(
								f"External ABM aggregate is too large: {path}"
							)
						if name.endswith(".csv") and has_private_columns(
							pd.read_csv(path, nrows=0).columns
						):
							raise ValueError(
								f"Participant data cannot be imported as an aggregate: {path}"
							)
						shutil.copy2(path, target / name)
				if kind == "reference" and render_figures:
					subprocess.run(
						[
							sys.executable,
							str(Path(__file__).with_name("c1.abm.py")),
							"figures",
							str(target),
						],
						check=True,
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
	p.add_argument(
		"--import-abm",
		action="store_true",
		help="Import aggregate ABM outputs before indexing",
	)
	p.add_argument("--abm-root", type=Path)
	p.add_argument("--abm-tf-root", type=Path)
	a = p.parse_args(argv)
	traits = a.Y.split(",")
	layers = a.biom.split(",")
	if any(not re.fullmatch(r"[A-Za-z0-9_]+", x) for x in traits) or not set(
		layers
	) <= {"prot", "met"}:
		p.error("Invalid outcome/layer")
	if a.max_table_mb < 1:
		p.error("--max-table-mb must be positive")
	if a.import_abm:
		import_abm_results(
			a.analysis_root,
			traits,
			layers,
			a.abm_root,
			a.abm_tf_root,
			render_figures=a.question_figures,
		)
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
		pd.concat(idx.candidates, ignore_index=True)
		if idx.candidates
		else pd.DataFrame(),
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
	"c4_panel_validation",
}
ORDER = [
	"c1_correlate",
	"c1_abm",
	"c2_cause",
	"c3_coloc",
	"pgs_focus",
	"c4_connect",
	"c4_panel_validation",
	"c4_explain",
	"c5_cellulation",
	"final",
	"shiny",
	"share",
]


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
		sources[key] = dict(source="cli" if explicit is not None else name if name in env else "default", value=value)
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
  ./le8.sh c4_explain --Y cvd_cad --biom prot
  ./le8.sh c1_abm --Y cvd_cad --biom prot
  ./le8.sh final --fit-joint --Y cvd_cad,ra --biom prot,met --replace TRUE

Modules: C1 correlation (c1_connect is accepted as an alias), C2 cause, C3 coloc,
C4 connect (including interactions and nonlinearity)/panel validation/explain, C5 cellulation, final, shiny. 
--abm-args is parsed as arguments, NEVER evaluated by a shell.
""",
	)
	p.add_argument(
		"modules",
		nargs="?",
		default="c1_correlate,c1_abm,c2_cause,c3_coloc,c4_connect,c4_panel_validation,c5_cellulation",
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
	p.add_argument(
		"--share-out", type=Path, help="Portable Shiny ZIP outside the analysis root"
	)
	p.add_argument("--ukb-phe", type=Path, default=None)
	p.add_argument("--cores", type=int, default=None)
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
	p.add_argument(
		"--abm-root",
		type=Path,
		help="Read existing external reference results; not a code dependency",
	)
	p.add_argument(
		"--abm-tf-root",
		type=Path,
		help="Read existing external TF results; not a code dependency",
	)
	p.add_argument(
		"--abm-backend",
		"--abm-engine",
		dest="abm_backend",
		choices=["reference", "tabicl", "tf", "both"],
		default="reference",
		help="Unified ABM backends; reference defaults to selective_attention on CUDA, tf is an alias for tabicl",
	)
	p.add_argument(
		"--abm-args",
		default="",
		help='Additional explicit ABM arguments, e.g. "--epochs 10 --device cuda"',
	)
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
	for name in ["atlas", "universe", "panels", "contrasts", "cigma-manifest", "cigma-results", "cigma-cells"]:
		p.add_argument("--" + name, type=Path)
	p.add_argument("--matched-draws", type=int, default=0)
	p.add_argument("--allow-untested-cigma", action="store_true")
	p.add_argument("--no-plots", action="store_true", help="C5 only: omit annotation plots")
	p.add_argument("--max-table-mb", type=int, default=128)
	p.add_argument("--port", type=int, default=int(os.getenv("LE8_SHINY_PORT", "3839")))
	p.add_argument("--host", default=os.getenv("LE8_SHINY_HOST", "127.0.0.1"))
	return p


def dispatch_main(argv=None):
	p = parser()
	a, extra = p.parse_known_args(argv)
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
	if any(not re.fullmatch(r"[A-Za-z0-9_]+", x) for x in traits) or len(traits) != len(
		set(traits)
	):
		p.error("Use unique, comma-separated outcome identifiers.")
	if (
		not layers
		or not set(layers) <= {"prot", "met"}
		or len(layers) != len(set(layers))
	):
		p.error("Invalid --biom")
	if (
		(a.cores is not None and a.cores < 1)
		or a.max_table_mb < 1
		or not 1024 <= a.port <= 65535
	):
		p.error("Invalid resource/port setting")
	if any(x is not None and x < 0 for x in [a.memory_limit_gb, a.memory_swap_gb]):
		p.error("Memory settings must be nonnegative")
	if (
		a.host not in {"127.0.0.1", "localhost", "::1"}
		and os.getenv("LE8_SHINY_ALLOW_NETWORK") != "TRUE"
	):
		p.error(
			"Network binding needs LE8_SHINY_ALLOW_NETWORK=TRUE and your own authentication/proxy."
		)
	aliases = {
		"c1_connect": "c1_correlate",
		"c4_focus": "c4_panel_validation",
		"s1_interact": "c4_connect",
		"s1_interaction": "c4_connect",
		"s2_nonlin": "c4_connect",
	}
	requested = [aliases.get(x, x) for x in a.modules.split(",")]
	if any(x not in ORDER for x in requested):
		p.error("Unknown module; see --help.")
	if a.run_abm and "c1_abm" not in requested:
		requested.append("c1_abm")
	if a.run_abm is False:
		requested = [x for x in requested if x != "c1_abm"]
	modules = [x for x in ORDER if x in requested]
	fitted = a.fit_reference or a.fit_joint or a.fit_genetic or a.details
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
	env.update(shared_env)
	env["LE8_SHARED_CONFIG_SOURCES"] = json.dumps(config_sources, sort_keys=True)
	if extra or a.ukb_phe:
		env["LE8_RESUME_COMPLETED"] = "FALSE"
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
		available = [(d.name, l) for d in a.analysis_root.iterdir() if d.is_dir() and re.fullmatch(r"[A-Za-z0-9_]+", d.name)
			for l in ("prot", "met") if (d / l).is_dir()]
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
			"--import-abm",
			"--max-table-mb",
			a.max_table_mb,
		]
		for key in ["abm_root", "abm_tf_root"]:
			if getattr(a, key) is not None:
				cmd += ["--" + key.replace("_", "-"), str(getattr(a, key))]
		if "final" in modules and not a.index_only:
			cmd += ["--question-figures"]
		call(cmd, env, a.dry_run)
		indexed = True

	abm_commands = []
	if "c1_abm" in modules:
		py = os.getenv("ABM_PYTHON") or sys.executable
		more = shlex.split(a.abm_args) + (extra if modules == ["c1_abm"] else [])
		mode_parser = argparse.ArgumentParser(add_help=False)
		mode_parser.add_argument("--abm-design")
		mode, _ = mode_parser.parse_known_args(more)
		if mode.abm_design == "selective_attention" and a.abm_backend != "reference":
			raise ValueError("selective_attention requires --abm-backend reference; tabicl/both are incompatible")
		for Y in traits:
			for layer in layers:
				selected = "tabicl" if a.abm_backend == "tf" else a.abm_backend
				for kind in (
					["reference", "tabicl"] if selected == "both" else [selected]
				):
					worker = HERE / "c1.abm.py"
					args = [
						"--backend",
						kind,
						"--device",
						"cuda",
						"--abm-design",
						"selective_attention" if kind == "reference" else "selective",
						"--Y",
						Y,
						"--biom",
						layer,
						"--analysis-root",
						str(a.analysis_root),
						"--seed",
						str(a.seed),
					]
					for key in ["group_file", "group_col", "end_date", "diagnosis_col", "outer_roster"]:
						if shared[key]:
							args += ["--" + key.replace("_", "-"), shared[key]]
					if shared["shared_covariates"]:
						args += ["--covariates", shared["shared_covariates"]]
					if a.ukb_phe:
						args += ["--ukb-phe", str(a.ukb_phe)]
					if a.cores:
						args += ["--cores", str(a.cores)]
					if a.replace == "TRUE":
						args += ["--replace"]
					for key in ["memory_limit_gb", "memory_swap_gb"]:
						v = getattr(a, key)
						if v is not None:
							args += ["--" + key.replace("_", "-"), str(v)]
					abm_commands.append([py, worker, *args, "--r-bin", a.r_bin, *more])
		# Fail before costly native C1 scans, using the selected ABM interpreter.
		# Explicit evaluate/project/device/download actions must execute only once.
		management = {"evaluate", "project", "--check-device", "--download-model", "--dry-run", "--help", "-h"}
		if a.preflight or not management.intersection(more):
			for cmd in abm_commands:
				call([*cmd, "--preflight"], env, a.dry_run)

	for module in modules:
		if module == "share":
			if modules != ["share"]:
				p.error("Run share separately after preparing final/shiny")
			output = a.share_out or a.analysis_root.parent / "le8-share.zip"
			if not a.dry_run:
				share_viewer(a.analysis_root, output, a.r_bin)
		elif module in NATIVE:
			call(["bash", engine, module, *base, *native, *extra], env, a.dry_run)
		elif module == "c1_abm":
			if not a.preflight:
				for cmd in abm_commands:
					call(cmd, env, a.dry_run)
		elif module == "c4_explain":
			if a.preflight:
				call(
					[
						a.r_bin,
						"-e",
						"p<-c('data.table','dplyr','survival','glmnet');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste('Missing',paste(m,collapse=', ')))",
					],
					env,
					a.dry_run,
				)
				continue
			for Y in traits:
				call(
					[a.r_bin, HERE / "c4.connect.R"],
					dict(env, Y=Y, LE8_PANEL_ACTION="reconstruct"),
					a.dry_run,
				)
		elif module == "c5_cellulation":
			cmd = [sys.executable, HERE / "c5.cellulation.py", *base,
				"--seed", str(a.seed), "--matched-draws", str(a.matched_draws)]
			for flag in ["preflight", "allow_untested_cigma", "no_plots"]:
				if getattr(a, flag):
					cmd += ["--" + flag.replace("_", "-")]
			if a.replace == "TRUE":
				cmd += ["--replace"]
			for key in ["atlas", "universe", "panels", "contrasts", "cigma_manifest", "cigma_results", "cigma_cells"]:
				if getattr(a, key) is not None:
					cmd += ["--" + key.replace("_", "-"), str(getattr(a, key))]
			call(cmd, env, a.dry_run)
		elif module == "final":
			fit_options = dict(
				traits=traits,
				layers=layers,
				seed=a.seed,
				native_options=extra,
				settings={
					k: v for k, v in os.environ.items() if k.startswith("FINAL_")
				},
			)
			if a.fit_reference or a.fit_joint or a.fit_genetic:
				ev = dict(
					env,
					FINAL_RUN_REFERENCE=str(a.fit_reference).upper(),
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
				for key in ["out", "abm_root", "abm_tf_root"]:
					if getattr(a, key) is not None:
						cmd += ["--" + key.replace("_", "-"), str(getattr(a, key))]
				if a.strict:
					cmd += ["--strict"]
				call(cmd, env, a.dry_run)
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
	root = Path(
		option("--analysis-root", os.getenv("LE8_ANALYSIS_ROOT", "/mnt/d/analysis/le8"))
	).resolve()
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
	if not traits or any(not re.fullmatch(r"[A-Za-z0-9_]+", y) for y in traits) or not layers or not layers <= {"prot", "met"}:
		raise ValueError("Invalid outcome/layer")
	modules = set(module_arg.split(","))
	scopes = None if modules & {"final", "shiny", "index", "share"} else {(y, b) for y in traits for b in layers}
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
		locks.enter_context(analysis_lock(Path("/tmp") / f"le8-{lock_name}.lock", shared=scopes is not None,
			waiting=f"Waiting for analysis/report transaction at {root}."))
		for y, b in sorted(scopes or []):
			key = hashlib.sha256(f"{root}/{y}/{b}".encode()).hexdigest()[:20]
			locks.enter_context(analysis_lock(Path("/tmp") / f"le8-scope-{key}.lock",
				waiting=f"{y}/{b} is already running at {root}; waiting for that scope."))
		return _run_table_workspace(argv, root, r_bin, module_arg, launch_shiny, train_abm, scopes, publication_lock, locks, option)
	finally:
		locks.close()


def _run_table_workspace(argv, root, r_bin, module_arg, launch_shiny, train_abm, scopes, publication_lock, locks, option):
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
						"le8_annotations",
						"logs",
					}
					for x in path.relative_to(root).parts
				):
					continue
				relative = path.relative_to(root)
				parts = relative.parts
				if not owns_scope(relative, scopes):
					continue
				# Aggregate-only ABM reports never need checkpoints or participant matrices.
				if not train_abm and any(
					x in {"abm_reference", "abm_tabicl", "abm_selective_attention"} for x in parts
				):
					if (
						path.suffix == ".rds"
						or path.suffix in {".npy", ".npz", ".pt", ".joblib", ".ckpt"}
						or any(
							x in {"input", "neural", "quality_neural", "checkpoints"}
							for x in parts
						)
					):
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
		for path in work.rglob("c1.abm_import_manifest.csv"):
			path.unlink()
		new_args = []
		skip = False
		def map_path(value):
			if value.startswith('--') and '=' in value:
				key, path = value.split('=', 1)
				return key + '=' + map_path(path)
			return str(work) + value[len(str(root)):] if value == str(root) or value.startswith(str(root) + '/') else value
		for i, value in enumerate(argv):
			if skip:
				skip = False
				continue
			if value == "--analysis-root":
				new_args += [value, str(work)]
				skip = True
			elif value.startswith("--analysis-root="):
				new_args.append("--analysis-root=" + str(work))
			elif (i and argv[i - 1] == '--abm-args') or value.startswith('--abm-args='):
				# A sibling such as results-external is not inside results.
				prefix, content = ('--abm-args=', value.split('=', 1)[1]) if value.startswith('--abm-args=') else ('', value)
				new_args.append(prefix + shlex.join([map_path(word) for word in shlex.split(content)]))
			else:
				new_args.append(map_path(value))
		if not any(
			x == "--analysis-root" or x.startswith("--analysis-root=") for x in new_args
		):
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
			print("[LE8] Shiny is already serving this result directory; the updated index is available on reload.", flush=True)
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
	# Workbooks retain verified aggregate exports; participant tables use named RDS.
	# Rendering and Shiny used the scratch CSV exchanges.
	print("[LE8] Consolidating result tables into workbooks and RDS...", flush=True)
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
		if path.name.endswith(
			(".csv", ".csv.gz", ".tsv", ".tsv.gz", ".tmp", ".pyc", ".log")
		):
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
	import zipfile

	root = root.resolve()
	output = output.resolve()
	if output.suffix.lower() != ".zip" or output.is_relative_to(root):
		raise ValueError("Use a ZIP path outside the analysis result directory")
	with tempfile.TemporaryDirectory(prefix="le8-share-", dir="/tmp") as scratch:
		bundle = Path(scratch) / "le8-share"
		(bundle / "shiny").mkdir(parents=True)
		(bundle / "f").mkdir()
		shutil.copy2(ROOT / "shiny/app.R", bundle / "shiny/app.R")
		shutil.copy2(HERE / "0.common.R", bundle / "f/0.common.R")
		subprocess.run(
			[
				r_bin,
				str(HERE / "0.common.R"),
				"--tables-share",
				str(root),
				str(bundle / "results"),
			],
			check=True,
		)
		(bundle / "README.md").write_text(
			"# LE8 结果查看包\n\n"
			"包含已登记的汇总表、PNG 和 Shiny 查看代码，不含个体 RDS、原始 UKB 数据或模型权重。\n\n"
			"1. 安装 R，并在 R 控制台安装依赖：\n\n"
			"```r\ninstall.packages(c('shiny', 'DT', 'data.table', 'ggplot2', 'digest', 'openxlsx', 'jsonlite'))\n```\n\n"
			"2. 解压整个文件夹。在解压目录运行：\n\n"
			"```sh\nRscript shiny/app.R\n```\n\n"
			"3. 浏览器打开 http://127.0.0.1:3839 。无需 Python、WSL 或原作者的数据路径。\n\n"
			"Windows 的 RStudio 也可把工作目录切到解压后的 shiny 文件夹，然后运行 `source('app.R')`。\n\n"
			"保留 results 内的相对目录结构；只改最外层目录的位置。工作簿包含精确来源数据，供查看器读取，"
			"请勿在 Excel 中覆盖保存。`--port 3840` 可更改端口。\n",
			encoding="utf-8",
		)
		check_env = dict(os.environ)
		for name in ["LE8_ANALYSIS_ROOT", "LE8_BROWSER_INDEX"]:
			check_env.pop(name, None)
		subprocess.run(
			[r_bin, str(bundle / "shiny/app.R"), "--review"], env=check_env, check=True
		)
		archive = Path(scratch) / "le8-share.zip"
		with zipfile.ZipFile(
			archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1
		) as stream:
			for path in sorted(bundle.rglob("*")):
				if path.is_file():
					if path.suffix.lower() not in {
						".r",
						".md",
						".xlsx",
						".png",
						".jpg",
						".jpeg",
					}:
						raise ValueError(f"Unexpected share artifact: {path}")
					stream.write(path, path.relative_to(bundle.parent))
		with zipfile.ZipFile(archive) as stream:
			if stream.testzip() is not None:
				raise ValueError("ZIP verification failed")
		output.parent.mkdir(parents=True, exist_ok=True)
		shutil.copy2(archive, output)
	print(f"[LE8] Share package: {output} ({output.stat().st_size / 1024**2:.1f} MiB)")


def main(argv=None):
	argv = list(sys.argv[1:] if argv is None else argv)
	if _table_runtime(argv):
		return
	if argv and argv[0] == "index":
		return index_main(argv[1:])
	if argv and argv[0] == "dispatch":
		argv = argv[1:]
	return dispatch_main(argv)


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		sys.exit(130)
	except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
		print("ERROR: " + str(exc), file=sys.stderr)
		sys.exit(2)
