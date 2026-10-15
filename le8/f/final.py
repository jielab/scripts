#!/usr/bin/env python3


# 🚩 final
"""Assemble LE8 final reports, systematic tables and question-led summaries.

Read existing aggregate results only; retain missing evidence and provenance.
Use ``report`` for publication figures and ``tables`` for a trait evidence atlas.
The shared index calls ``build_questions`` and ``question_figures`` directly.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
import hashlib
import html
from importlib.util import module_from_spec, spec_from_file_location
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile
import sys
import textwrap
from typing import Callable

import numpy as np
import pandas as pd


def load_module(filename):
	"""Load a dotted script filename under a stable Python module identifier."""
	path = Path(__file__).with_name(filename)
	name = path.stem.replace(".", "_")
	if name not in sys.modules:
		spec = spec_from_file_location(name, path)
		module = module_from_spec(spec)
		sys.modules[name] = module
		try:
			spec.loader.exec_module(module)
		except BaseException:
			sys.modules.pop(name, None)
			raise
	return sys.modules[name]


_index = load_module("0.common.py")
final_dir = _index.final_dir
result_file = _index.result_file
from abm_summary import table_paths as abm_table_paths, source_root as abm_source_root, summarize as twin_summary, BOUNDARY as ABM_BOUNDARY
atomic_csv = _index.atomic_csv
digest = _index.digest
safe_path = _index.safe_path
PRIVATE_COLUMNS = _index.PRIVATE_COLUMNS

MODULE_ORDER = ("C1", "C2", "C4", "C3", "C5")
MODULE_TITLES = dict(C1="Individual predictability", C2="Genetic vs acquired/state",
	C4="LE8-supervised explainability", C3="Locus / colocalization support", C5="Cell-level context")
MET_TABLES = {
	"met_pgs_pairs": "c1.pgs_actual_concordance.csv",
	"met_anchor": "c1.L_VLDL_TG_pct_deep_dive.csv",
	"met_pgs_conditional": "c1.L_VLDL_TG_pct_pgs_conditional.csv",
	"met_measured_conditional": "c1.L_VLDL_TG_pct_measured_conditional.csv",
	"met_trajectory": "c1.L_VLDL_TG_pct_riskset_trajectory.csv",
}
MET_ANALYSES = ("Incident", "Attained-age delayed entry", "Baseline prevalent")
MET_BOUNDARY = (
	"PGS represents inherited propensity. High prospective concordance with low prevalent concordance "
	"is compatible with disease-state/treatment/survivor remodeling; it does not prove reverse causality. "
	"Behavioral LE4 is the primary lifestyle adjustment; full LE8 is a sensitivity analysis because "
	"BP/glucose/lipids/BMI can overlap with metabolic intermediates."
)
ANCHOR_BOUNDARY = (
	"L_VLDL_TG.pct is a prespecified illustrative compositional percentage trait. A negative inherited "
	"association does not mean higher absolute VLDL triglyceride burden is protective. Compare with "
	"L_VLDL_TG, Total_TG, ApoB and VLDL_size where available. Conditional PGS analyses are "
	"genetic-architecture sensitivity analyses, not multivariable MR."
)


def finite_column(d, name):
	return pd.to_numeric(d.get(name, pd.Series(np.nan, index=d.index)), errors="coerce").replace([np.inf, -np.inf], np.nan)




def metabolite_concordance(pairs):
	"""Pairwise finite effects, not significance-selected rows; recomputed on every build."""
	if "layer" in pairs:
		pairs = pairs[pairs.layer.eq("met")]
	rows = []
	for analysis in MET_ANALYSES:
		z = pairs[pairs.analysis.eq(analysis)].copy() if "analysis" in pairs else pd.DataFrame()
		x, y = finite_column(z, "pgs_full_beta"), finite_column(z, "observed_beta")
		valid = x.notna() & y.notna()
		x, y = x[valid], y[valid]
		rho = x.rank().corr(y.rank()) if len(x) > 1 and x.nunique() > 1 and y.nunique() > 1 else np.nan
		row = dict(analysis=analysis, N_pairs=len(x), spearman_rho=rho,
			same_direction_fraction=float((np.sign(x) == np.sign(y)).mean()) if len(x) else np.nan,
			status="exploratory" if np.isfinite(rho) else "insufficient evidence" if len(x) else "unavailable")
		for c in ("Y", "layer", "source_file", "source_sha256"):
			if c in z:
				row[c] = ";".join(z[c].dropna().astype(str).unique())
		row["source_rows"] = ";".join(z.loc[valid, "source_row"].astype(str)) if "source_row" in z else ""
		rows.append(row)
	return pd.DataFrame(rows)


def c4_endpoint_summary(contrasts, proxy):
	ct = contrasts.copy()
	for col, value in (("stratum", "All"), ("landmark", 0)):
		ct = ct.loc[ct[col].eq(value)] if col in ct else ct.iloc[:0]
	ct = ct.loc[ct.reference.str.match(r"^NS_\d+$", na=False)] if "reference" in ct else ct.iloc[:0]
	ct = ct.loc[ct.comparison_valid.eq(True)] if "comparison_valid" in ct else ct.iloc[:0]
	ct = ct.loc[finite_column(ct, "delta_AUC").notna()]
	pos = int((finite_column(ct, "delta_lo") > 0).sum())
	prediction = "insufficient evidence" if ct.empty else "supported" if pos == len(ct) else "no universal superiority"
	pr = proxy.loc[finite_column(proxy, "delta_R2_vs_NS").notna()]
	confirmed = int(((finite_column(pr, "lo") > 0) & (finite_column(pr, "FDR_all_proxy_contrasts") < .05)).sum())
	explain = "insufficient evidence" if pr.empty else "supported" if confirmed else "exploratory" if (finite_column(pr, "delta_R2_vs_NS") > 0).any() else "not supported"
	result = ("No universal prediction superiority was observed. " if prediction == "no universal superiority" else "")
	result += f"Prediction: {prediction}; {pos}/{len(ct)} same-budget paired ΔAUC intervals above zero. "
	result += f"Explainability: {explain}; {int((finite_column(pr, 'delta_R2_vs_NS') > 0).sum())}/{len(pr)} held-out ΔR²>0, {confirmed} with interval>0 and FDR<0.05."
	return dict(prediction_status=prediction, explainability_status=explain, key_result=result)


# Publication report and figure assembly.
ROLES = {
	**{name: ("c1_correlate", filename) for name, filename in MET_TABLES.items()},
	"incident": ("c1_correlate", "pwas_incident_adj2.csv"),
	"prevalent": ("c1_correlate", "pwas_prevalent_adj2.csv"),
	"cohort": ("c1_correlate", "c1.cohort.csv"),
	"temporal": ("c1_correlate", "c1.directionality_triage.csv"),
	"paired": ("c1_correlate", "c1.paired_pgs_measured.csv"),
	"enrichment": ("c1_correlate", "c1.enrichment_incident_sig.csv"),
	"mr": ("c2_cause", "c2.MR_all.csv"),
	"mr_support": ("c2_cause", "c2.mr_incident_prevalent_support.csv"),
	"mr_comparison": ("c2_cause", "c2.mr_incident_prevalent.csv"),
	"mr_candidates": ("c2_cause", "c2.top_candidates.csv"),
	"mr_scope": ("c2_cause", "c2.method_scope_audit.csv"),
	"coloc_audit": ("c3_coloc", "c3.credible_set_audit.csv"),
	"coloc_by_locus": ("c3_coloc", "c3.credible_set_by_locus.csv"),
	"mr_grades": ("c2_cause", "c2.evidence_grades.csv"),
	"dandelion": ("c2_cause", "c2.dandelion_input_audit.csv"),
	"coloc": ("c3_coloc", "c3.coloc_summary.csv"),
	"membership": ("c4_connect", "c4.proxy_membership_YS_YSP_NS.csv"),
	"mediation": ("c4_connect", "c4.mediation_all.csv"),
	"pgs_bridges": ("c4_connect", "c4.matched_PGS_bridges.csv"),
	"focus_metrics": ("c4_connect", "c4.focus.metrics.csv"),
	"focus_contrasts": ("c4_connect", "c4.focus.contrasts.csv"),
	"focus_pillars": ("c4_connect", "c4.focus.pillar_counts.csv"),
	"focus_proxy": ("c4_connect", "c4.focus.proxy_accuracy.csv"),
	"deployed_concept_fidelity": ("c4_connect", "c4.focus.deployed_concept_fidelity.csv"),
	"explain_contrasts": ("c4_connect", "c4.explain.contrasts.csv"),
	"focus_members": ("c4_connect", "c4.focus.panel_members.csv"),
	"focus_coefficients": ("c4_connect", "c4.focus.model_coefficients.csv"),
	"focus_design": ("c4_connect", "c4.focus.design.csv"),
	"focus_diagnostics": ("c4_connect", "c4.focus.fit_diagnostics.csv"),
	"focus_calibration": ("c4_connect", "c4.focus.calibration.csv"),
	"focus_heterogeneity": ("c4_connect", "c4.focus.heterogeneity.csv"),
	"performance": ("c1_prediction", "prediction_summary.csv"),
	"budget": ("c1_prediction", "review_budget_metrics.csv"),
	"paired_delta": ("c1_prediction", "review_paired_delta_CI.csv"),
	"prediction_design": ("c1_prediction", "review_design.csv"),
	"calibration": ("c1_prediction", "review_calibration.csv"),
	"leadtime": ("c1_prediction", "leadtime_discrimination.csv"),
	"cell_status": ("c5_cellulation", "c5.cellulation_status.csv"),
	"cell_enrichment": ("c5_cellulation", "c5.cell.enrichment.csv"),
	"cell_coverage": ("c5_cellulation", "c5.cell.coverage.csv"),
	"cigma": ("c5_cellulation", "c5.CIGMA_annotation.csv"),
	"cell_contrasts": ("c5_cellulation", "c5.cell.panel_contrasts.csv"),
	"cell_evidence": ("c5_cellulation", "c5.evidence_long.csv"),
	"concept_coefficients": ("c4_connect", "c4.focus.concept_coefficients.csv"),
	"concept_status": ("c4_connect", "c4.focus.concept_status.csv"),
	"mr_lolo": ("c2_cause", "c2.dandelion_leave_locus_out.csv"),
	"coloc_same_locus": ("c3_coloc", "c3.same_locus_evidence.csv"),
	"susie_pairs": ("c3_coloc", "c3.susie_signal_pairs.csv"),
	"mr_signal_evidence": ("c3_coloc", "c3.same_locus_signal_evidence.csv"),
	"susie_diagnostics": ("c3_coloc", "c3.susie_diagnostics.csv"),
}


def has(d: pd.DataFrame, *cols: str) -> bool:
	return not d.empty and set(cols) <= set(d.columns)


def num(d: pd.DataFrame, name: str) -> pd.Series:
	return pd.to_numeric(d[name], errors="coerce")


def display_budget(d: pd.DataFrame) -> int:
	requested = int(os.environ.get("FINAL_DISPLAY_BUDGET", "10"))
	values = sorted(set(pd.to_numeric(d.get("budget", pd.Series(dtype=float)), errors="coerce").dropna()))
	values = [int(v) for v in values if v > 0]
	# Choose using configured assay budgets only, never observed performance.
	return min(values, key=lambda v: (abs(v-requested), v)) if values else requested


def clean_name(name: str) -> str:
	return re.sub(r"[^\w.-]+", "_", name)[:160]


@dataclass
class Panel:
	name: str
	title: str
	png: Path
	csv: Path
	caption: str
	source_roles: list[str]


class Report:
	def __init__(self, arg):
		import matplotlib

		matplotlib.use("Agg")
		import matplotlib.pyplot as plt

		plt.rcParams["svg.fonttype"] = "none"
		plt.rcParams["pdf.fonttype"] = 42
		self.plt = plt
		self.arg = arg
		self.root = arg.analysis_root.resolve()
		self.out = (arg.out or self.root / "final").resolve()
		if self.out != (self.root / "final").resolve():
			raise ValueError(
				"Final output must be <analysis-root>/final; use --analysis-root to choose a project."
			)
		self.out.mkdir(parents=True, exist_ok=True)
		# Replace generated report views, including stale panels from prior arguments.
		# Fitted prediction data in outcome subdirectories remain the analysis inputs.
		for name in ("panels", "tables", "_previous_figures"):
			directory = self.out / name
			if directory.is_dir():
				shutil.rmtree(directory)
		for old in self.out.glob("Fig*.*"):
			if old.is_file():
				old.unlink()

		self.table_dir = self.out
		self.table_dir.mkdir(exist_ok=True)
		self.panel_dir = Path(tempfile.mkdtemp(prefix="le8-panels-", dir="/tmp"))
		self.panel_dir.mkdir(exist_ok=True)
		self.registry = []
		self.findings = []
		self.tables = {}
		self.panels = {}
		self.figures = []
		self.source_files = {}
		self.mode = (
			"verified aggregate excerpts"
			if (self.root.parent / "SNAPSHOT_NOTICE.md").exists()
			else "local aggregate results"
		)

	def flag(self, trait, layer, code, severity, detail, roles=""):
		self.findings.append(
			dict(
				trait=trait,
				layer=layer,
				code=code,
				severity=severity,
				detail=detail,
				source_roles=roles,
			)
		)

	def read(
		self, trait: str, layer: str, role: str, paths: list[Path]
	) -> pd.DataFrame:
		key = f"{trait}.{layer}.{role}"
		available = [p for p in paths if p.is_file()]
		chosen = available[0] if available else paths[0]
		status = "unavailable"
		d = pd.DataFrame()
		source_hash = ""
		if available:
			try:
				# Explicit aggregate allowlist, never participant-level RDS/parquet.
				d = pd.read_csv(chosen, sep="\t" if ".tsv" in chosen.name else ",", low_memory=False)
				source_hash = digest(chosen)
				status = "available" if len(d) else "header_only"
			except pd.errors.EmptyDataError:
				status = "empty"
			except Exception as exc:
				status = "unreadable"
				self.flag(
					trait,
					layer,
					"TABLE_READ_ERROR",
					"blocking",
					f"{chosen}: {exc}",
					role,
				)
			if len(available) > 1 and len({digest(p) for p in available}) > 1:
				self.flag(
					trait,
					layer,
					"ALTERNATIVE_OUTPUT_CONFLICT",
					"warning",
					f"Preferred {chosen}; differing alternatives exist. No merging across runs.",
					role,
				)
		self.registry.append(
			dict(
				key=key,
				trait=trait,
				layer=layer,
				role=role,
				source=str(chosen),
				status=status,
				rows=len(d),
				sha256=source_hash,
				alternatives=";".join(map(str, available[1:])),
			)
		)
		self.tables[key] = d
		self.source_files[key] = str(chosen)
		return d

	def get(self, trait, layer, role):
		return self.tables.get(f"{trait}.{layer}.{role}", pd.DataFrame())

	def load(self):
		for trait in self.arg.Y.split(","):
			for layer in self.arg.biom.split(","):
				b = self.root / trait / layer
				for role, (module, filename) in ROLES.items():
					if layer == "met":
						filename = filename.replace("pwas_", "mwas_")
					paths = [result_file(self.root, trait, layer, module, filename)]
					if role == "coloc":
						paths.append(b / module / "c3.credible_set_by_locus.csv")
					self.read(trait, layer, role, paths)
				for role in [
					"components",
					"contrasts",
					"bootstrap",
					"source_status",
					"calibration",
				]:
					self.read(
						trait,
						layer,
						"pgs_" + role,
						[
							final_dir(self.root, trait)
							/ f"Fig1.PGS_integrated.{layer}_{role}.csv"
						],
					)
				for role, path in abm_table_paths(self.root, trait, layer).items():
					self.read(trait, layer, role, [path])
				cell_status = self.get(trait, layer, "cell_status")
				if has(cell_status, "analysis", "status"):
					for analysis, roles in [
						("cell_expression", ["cell_enrichment", "cell_coverage"]),
						("CIGMA", ["cigma"]),
					]:
						current = cell_status[cell_status.analysis.eq(analysis)]
						if current.empty or not current.status.eq("completed").all():
							for role in roles:
								if not self.get(trait, layer, role).empty:
									self.flag(
										trait,
										layer,
										"CELL_RESULT_NOT_CURRENT",
										"blocking",
										role
										+ ": suppressed because the current status is not completed.",
										role,
									)
									self.tables[f"{trait}.{layer}.{role}"] = (
										pd.DataFrame()
									)
				self.audit(trait, layer)
		# Explicit completion status: code and cohort inventories are not completed validation.
		for trait in self.arg.Y.split(","):
			joint = [final_dir(self.root, trait, kind="joint") / "metrics.csv"]
			self.read(trait, "joint", "joint_metrics", joint)
			for nm in [
				"paired_contrasts",
				"PRS_omics_strata",
				"biomarker_PGS_omics_strata",
				"panel_counts",
				"summary_provenance",
				"primary_status",
				"prs_provenance",
				"fold_C_index",
				"calibration",
			]:
				self.read(
					trait,
					"joint",
					nm,
					[
						final_dir(self.root, trait, kind="joint") / f"{nm}.csv",
					],
				)
			if self.get(trait, "joint", "joint_metrics").empty:
				self.flag(
					trait,
					"joint",
					"JOINT_VALIDATION_UNAVAILABLE",
					"gap",
					"No recognized completed common-cohort joint metric table. Do not infer PRS/omics complementarity from separate cohorts. Existing native fitting is opt-in.",
				)

	def audit(self, trait, layer):
		if layer == "met":
			pairs = self.get(trait, layer, "met_pgs_pairs")
			if len(pairs):
				rows = metabolite_concordance(pairs)
				self.flag(trait, layer, "METABOLITE_INHERITED_MEASURED_CONCORDANCE", "finding",
					"; ".join(f"{r.analysis}: N={r.N_pairs}, Spearman rho={r.spearman_rho:.4g}, same-direction fraction={r.same_direction_fraction:.4g}" for r in rows.itertuples()) + ". " + MET_BOUNDARY,
					"met_pgs_pairs")
		d = self.get(trait, layer, "performance")
		if has(d, "N", "events", "C_index"):
			if d.N.nunique() > 1 or d.events.nunique() > 1:
				self.flag(
					trait,
					layer,
					"UNMATCHED_LEGACY_COMPARISON",
					"blocking",
					"Legacy model rows do not share N/events; no paired improvement inferred.",
					"performance",
				)
			if has(d, "n_selected") and num(d, "n_selected").nunique() > 1:
				self.flag(
					trait,
					layer,
					"UNEQUAL_ASSAY_BUDGETS",
					"warning",
					"Legacy panel comparison uses unequal assay budgets. It is not the equal-budget test.",
					"performance",
				)
		co = self.get(trait, layer, "focus_coefficients")
		members = self.get(trait, layer, "focus_members")
		effective = []
		if has(co, "model", "variable", "beta") and has(members, "model", "feature"):
			for model, sub in members.dropna(subset=["feature"]).groupby("model"):
				fs = set(sub.feature.astype(str))
				cc = co[co.model == model].copy()
				if model.startswith("YSconcept"):
					# The deployed risk coefficients address concepts; assay loadings
					# belong to the first stage and are shown in their own table.
					cb = cc[cc.variable.astype(str).str.startswith("concept_")]
					dg = self.get(trait, layer, "focus_diagnostics")
					dg = dg[dg.model.eq(model)] if "model" in dg else pd.DataFrame()
					effective.append(dict(model=model, nominal=len(fs), effective=np.nan,
						effective_concepts=int((num(cb,"beta").abs()>self.arg.coefficient_epsilon).sum()),
						molecular_lp_sd_test=float(dg.molecular_lp_sd_test.iloc[0]) if len(dg) and "molecular_lp_sd_test" in dg else np.nan,
						interpretation="Concept risk coefficients and assay loadings are distinct; use saved LP contributions to assess molecular contribution"))
					continue
				# Includes dummy columns for categorical features; ordinary omics are numeric.
				use = cc.variable.astype(str).map(
					lambda x: any(x == f or x.startswith(f + "__") for f in fs)
				)
				bc = cc[use]
				beta = num(bc, "beta")
				effective_features = {
					f
					for f in fs
					if any(
						(
							bc.variable.astype(str).eq(f)
							| bc.variable.astype(str).str.startswith(f + "__")
						)
						& (beta.abs() > self.arg.coefficient_epsilon)
					)
				}
				row = dict(
					model=model,
					nominal=len(fs),
					effective=len(effective_features) if len(bc) else np.nan,
					max_abs_beta=beta.abs().max() if len(beta) else np.nan,
					coefficient_threshold=self.arg.coefficient_epsilon,
					panel_hash=hashlib.sha256(
						"\n".join(sorted(fs)).encode()
					).hexdigest(),
					coefficient_rows=len(bc),
				)
				effective.append(row)
				if len(fs) and not len(bc):
					self.flag(
						trait,
						layer,
						"MOLECULAR_COEFFICIENT_ROWS_MISSING",
						"gap",
						f"{model}: panel exists but molecular coefficient rows are unavailable; effective count is NA, not zero.",
						"focus_coefficients;focus_members",
					)
				if len(fs) and len(bc) and not effective_features:
					self.flag(
						trait,
						layer,
						"MOLECULAR_COEFFICIENT_COLLAPSE",
						"blocking",
						f"{model}: nominal {len(fs)} assays, effective 0 at |beta|>{self.arg.coefficient_epsilon:g}. This fitted model provides no appreciable omics contribution; do not claim successful compact molecular prediction.",
						"focus_coefficients;focus_members",
					)
			eff = pd.DataFrame(effective)
			self.tables[f"{trait}.{layer}.effective_assays"] = eff
			eff.to_csv(
				self.table_dir / f"{trait}.{layer}.effective_assays.csv", index=False
			)
		loc = self.get(trait, layer, "coloc_by_locus")
		if has(loc, "feature", "locus", "PP.H4_robust_min"):
			parsed = []
			for _, rr in loc.iterrows():
				match = re.fullmatch(r"chr([^:]+):(\d+)-(\d+)", str(rr.locus))
				if match:
					chrom, lo, hi = match.groups()
					parsed.append(
						(
							chrom,
							int(lo),
							int(hi),
							str(rr.feature),
							str(rr.get("policy_coloc_pass", "")).lower() in {"true", "1"},
						)
					)
			windows = []
			for chrom in sorted({v[0] for v in parsed}):
				group = sorted((v for v in parsed if v[0] == chrom), key=lambda v: v[1])
				for ch, lo, hi, feat, pp in group:
					if (
						windows
						and windows[-1]["chrom"] == ch
						and lo <= windows[-1]["end"]
					):
						windows[-1]["end"] = max(windows[-1]["end"], hi)
						windows[-1]["features"].add(feat)
						windows[-1]["entries"] += 1
						windows[-1]["robust_entries"] += int(pp)
					else:
						windows.append(
							dict(
								chrom=ch,
								start=lo,
								end=hi,
								features={feat},
								entries=1,
								robust_entries=int(pp),
							)
						)
			for win in windows:
				win["n_features"] = len(win["features"])
				win["features"] = ";".join(sorted(win["features"]))
				win["scope"] = (
					"Overlap-merged tested intervals, not LD-defined independent causal loci"
				)
			ww = pd.DataFrame(windows)
			if len(ww):
				self.tables[f"{trait}.{layer}.coloc_merged_windows"] = ww
				ww.to_csv(
					self.table_dir / f"{trait}.{layer}.coloc_merged_windows.csv",
					index=False,
				)
				if len(ww) < len(loc):
					self.flag(
						trait,
						layer,
						"COLOC_REGIONAL_REDUNDANCY",
						"warning",
						f"{len(loc)} feature-locus entries map to {len(ww)} overlap-merged tested intervals; do not count correlated entries as independent mechanisms.",
						"coloc_by_locus",
					)
		contrasts = self.get(trait, layer, "focus_contrasts")
		if has(contrasts, "delta_AUC", "delta_lo", "delta_hi"):
			zeros = (
				(num(contrasts, "delta_AUC").abs() < 1e-14)
				& (num(contrasts, "delta_lo").abs() < 1e-14)
				& (num(contrasts, "delta_hi").abs() < 1e-14)
			)
			if zeros.any():
				self.flag(
					trait,
					layer,
					"IDENTICAL_CONTRASTS",
					"warning",
					f"{int(zeros.sum())} contrasts have exactly zero delta and interval. Check identical panels or clinical-only shrinkage; not independent replication.",
					"focus_contrasts;effective_assays",
				)
		pillars = self.get(trait, layer, "focus_pillars")
		if has(pillars, "component", "n"):
			empty = pillars[num(pillars, "n") == 0]
			if len(empty):
				self.flag(
					trait,
					layer,
					"INCOMPLETE_LE8_COVERAGE",
					"warning",
					"No qualifying proxies for: "
					+ ", ".join(sorted(set(empty.component.astype(str))))
					+ ". Do not force assignment to all eight pillars.",
					"focus_pillars",
				)
		dan = self.get(trait, layer, "dandelion")
		if not dan.empty and (has(dan, "metric", "value") or "primary_eligible" in dan):
			v = (
				dict(zip(dan.metric.astype(str), dan.value.astype(str)))
				if "metric" in dan
				else dan.iloc[0].astype(str).to_dict()
			)
			if v.get("primary_eligible", "").upper() != "TRUE":
				self.flag(
					trait,
					layer,
					"DANDELION_EXPLORATORY_ONLY",
					"warning",
					v.get("analysis_class", "Primary eligibility unverified")
					+ "; exclude from independent causal confirmation counts.",
					"dandelion;dandelion_lolo;dandelion_native;state_projection;age_models",
				)
		components = self.get(trait, layer, "pgs_components")
		boot = self.get(trait, layer, "pgs_bootstrap")
		if has(
			components,
			"feature",
			"model",
			"term",
			"sample_hash",
			"covariates",
			"N",
			"events",
		):
			joined = components[components.model.eq("joint")]
			for feature, z in joined.groupby("feature"):
				for groupcols in [["scope", "landmark", "end"]]:
					actual = [c for c in groupcols if c in z]
					groups = z.groupby(actual, dropna=False) if actual else [(None, z)]
					for _, zz in groups:
						terms = set(zz.term.astype(str))
						if not {".G", ".R"} <= terms:
							self.flag(
								trait,
								layer,
								"PGS_MISSING_COMPONENT",
								"blocking",
								str(feature),
								"pgs_components",
							)
						if any(
							zz[c].nunique(dropna=False) > 1
							for c in ["sample_hash", "covariates", "N", "events"]
						):
							self.flag(
								trait,
								layer,
								"PGS_UNMATCHED_COMPONENTS",
								"blocking",
								str(feature),
								"pgs_components",
							)
			wanted = {"GDF15", "MMP12"} & set(joined.feature)
			available = set(boot.feature) if "feature" in boot else set()
			if wanted - available:
				self.flag(
					trait,
					layer,
					"ANCHOR_REFIT_BOOTSTRAP_MISSING",
					"gap",
					"Calibration-refit bootstrap unavailable for "
					+ ", ".join(sorted(wanted - available))
					+ "; analytic conditional intervals are not refitted-bootstrap confirmation.",
					"pgs_bootstrap;pgs_components",
				)
			self.flag(
				trait,
				layer,
				"PGS_DISCOVERY_INDEPENDENCE",
				"warning",
				"A cross-fitted calibration does not remove overlap in original GWAS weights. Genetic/residual differences are not automatically causal or pure lifestyle effects.",
				"pgs_components",
			)
		ss = self.get(trait, layer, "pgs_source_status")
		if (
			has(ss, "status")
			and ss.status.astype(str).str.contains("unavailable", case=False).any()
		):
			self.flag(
				trait,
				layer,
				"SOURCE_SCORE_UNAVAILABLE",
				"gap",
				"Source-specific score analysis is marked unavailable; no cis/trans/MHC-exclusion null finding can be claimed from these placeholders.",
				"pgs_source_status",
			)
		state = twin_summary(self.get(trait, layer, 'abm_metrics'))
		self.flag(trait, layer, 'ABM_TWIN_' + state['status'].upper(), 'finding',
			state['key_result'] + ' ' + ABM_BOUNDARY, 'abm_metrics;abm_paired;abm_support')

	def panel(
		self,
		name: str,
		title: str,
		d: pd.DataFrame,
		draw: Callable,
		caption: str,
		roles: list[str],
		size=(7.2, 4.7),
	):
		if d.empty:
			return None
		name = clean_name(name)
		fig, ax = self.plt.subplots(figsize=size)
		draw(ax, d.copy())
		ax.set_title(title, loc="left", fontsize=11, pad=12)
		ax.spines[["top", "right"]].set_visible(False)
		ax.tick_params(labelsize=8)
		fig.tight_layout(pad=1.3)
		paths = {
			"png": self.panel_dir / (name + ".png"),
			"csv": self.out / (name + ".csv"),
		}
		fig.savefig(paths["png"], dpi=170, bbox_inches="tight")
		self.plt.close(fig)
		d.to_csv(paths["csv"], index=False)
		p = Panel(name, title, paths["png"], paths["csv"], caption, roles)
		self.panels[name] = p
		return p

	def c4_prediction_panels(self, trait, layer):
		metrics = self.get(trait, layer, "focus_metrics")
		if not has(
			metrics,
			"model",
			"budget",
			"stratum",
			"landmark",
			"horizon",
			"AUC",
			"N",
			"status",
		):
			return [], pd.DataFrame()
		budget = display_budget(metrics)
		models = [
			"Clinical",
			"NS_10",
			"YS_Yin_10",
			"YS_YinYang_10",
			"YSplus_Yin_10",
			"YSplus_YinYang_10",
			"YSbalanced_Yin_10",
			"YSbalanced_YinYang_10",
		]
		models = [m.replace("_10", f"_{budget}") for m in models]
		z = metrics[
			metrics.stratum.eq("All")
			& num(metrics, "landmark").eq(0)
			& num(metrics, "horizon").eq(10)
			& num(metrics, "budget").isin([0, budget])
			& metrics.model.isin(models)
			& metrics.status.eq("ok")
		].copy()
		if z.empty:
			return [], z
		if z.model.duplicated().any() or num(z, "N").nunique() != 1:
			self.flag(
				trait,
				layer,
				"C4_DISPLAY_COHORT_MISMATCH",
				"blocking",
				"The fixed display slice has duplicate models or unequal validation N.",
				"focus_metrics",
			)
			return [], pd.DataFrame()
		z = (
			z.assign(_order=z.model.map({m: i for i, m in enumerate(models)}))
			.sort_values("_order")
			.drop(columns="_order")
		)
		z["trait"], z["layer"] = trait, layer

		def labels(d):
			return d.model.str.replace(rf"_{budget}$", "", regex=True).str.replace(
				"_", " ", regex=False
			)

		def intervals(ax, data, estimate, lower, upper):
			y = np.arange(len(data))
			point = num(data, estimate)
			ax.scatter(point, y, s=22, zorder=3)
			if {lower, upper} <= set(data):
				lo, hi = num(data, lower), num(data, upper)
				valid = np.isfinite(lo) & np.isfinite(hi) & lo.le(hi)
				ax.hlines(y[valid], lo[valid], hi[valid], linewidth=1.3)
			ax.set_yticks(y, labels(data))
			ax.invert_yaxis()

		def discrimination(ax, data):
			intervals(ax, data, "AUC", "AUC_lo", "AUC_hi")
			clinical = data[data.model.eq("Clinical")]
			if len(clinical):
				ax.axvline(
					float(clinical.AUC.iloc[0]),
					color="grey",
					linestyle="--",
					linewidth=0.8,
				)
			ax.set_xlabel("10-year IPCW AUC; 95% bootstrap interval")

		panels = [
			self.panel(
				f"Fig4_{trait}_{layer}",
				f"{trait} {layer}: {budget}-assay budget, N={int(num(z, 'N').iloc[0])}",
				z,
				discrimination,
				"All validation participants, baseline landmark. Clinical uses zero assays; other rows retain their reported actual assay counts.",
				[f"{trait}.{layer}.focus_metrics"],
			)
		]
		contrasts = self.get(trait, layer, "focus_contrasts")
		if has(
			contrasts,
			"model",
			"reference",
			"stratum",
			"landmark",
			"horizon",
			"delta_AUC",
		):
			d = contrasts[
				contrasts.stratum.eq("All")
				& num(contrasts, "landmark").eq(0)
				& num(contrasts, "horizon").eq(10)
				& contrasts.reference.eq(f"NS_{budget}")
				& contrasts.model.isin(z.model)
			].copy()
			if d.model.duplicated().any():
				self.flag(
					trait,
					layer,
					"C4_DISPLAY_CONTRAST_DUPLICATE",
					"blocking",
					"Duplicate paired contrasts in the fixed display slice.",
					"focus_contrasts",
				)
			else:
				d = (
					d.assign(_order=d.model.map({m: i for i, m in enumerate(models)}))
					.sort_values("_order")
					.drop(columns="_order")
				)

				def paired(ax, data):
					intervals(ax, data, "delta_AUC", "delta_lo", "delta_hi")
					ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
					ax.set_xlabel("Paired AUC difference versus NS; 95% interval")

				panels.append(
					self.panel(
						f"Fig4_{trait}_{layer}_paired",
						f"{trait} {layer}: supervision versus NS",
						d,
						paired,
						"Frozen-fit paired bootstrap; zero and adverse differences are retained.",
						[f"{trait}.{layer}.focus_contrasts"],
					)
				)
		return [p for p in panels if p is not None], z

	def compose(self, name: str, title: str, panels: list[Panel | None], caption: str):
		panels = [p for p in panels if p is not None]
		if not panels:
			for ext in ["pdf", "png"]:
				old = self.out / (name + "." + ext)
				if old.is_file():
					old.unlink()
			self.figures.append(
				dict(
					figure=name,
					status="unavailable",
					panels=0,
					path="",
					caption=caption,
				)
			)
			return
		# Compose the individual panels into a single PNG figure.
		columns = 2 if len(panels) > 1 else 1
		rows = math.ceil(len(panels) / columns)
		width = 1080
		cell_h = 305
		footer_lines = textwrap.wrap(caption, 165)
		height = 65 + rows * cell_h + max(55, len(footer_lines) * 12 + 15)
		from PIL import Image, ImageDraw, ImageFont
		from matplotlib import font_manager

		factor = 2
		image = Image.new("RGB", (width * factor, int(height * factor)), "white")
		draw = ImageDraw.Draw(image)
		regular = font_manager.findfont(
			font_manager.FontProperties(family="DejaVu Sans")
		)
		bold = font_manager.findfont(
			font_manager.FontProperties(family="DejaVu Sans", weight="bold")
		)
		draw.text(
			(60, 16),
			name + " | " + title,
			fill="black",
			font=ImageFont.truetype(bold, 29),
		)
		draw.text(
			(60, 65),
			"DRAFT - internal results; " + self.mode + "; no new model fitting",
			fill="black",
			font=ImageFont.truetype(regular, 16),
		)
		for i, panel in enumerate(panels):
			col = i % columns
			row = i // columns
			cell_w = (width - 60) / columns
			tile = Image.open(panel.png).convert("RGB")
			tile.thumbnail(
				(int((cell_w - 12) * factor), int((cell_h - 24) * factor)),
				Image.Resampling.LANCZOS,
			)
			x = int((30 + col * cell_w + 7) * factor)
			y = int((65 + row * cell_h + 18) * factor)
			image.paste(tile, (x, y))
			draw.text(
				(int((30 + col * cell_w) * factor), int((65 + row * cell_h) * factor)),
				chr(65 + i),
				fill="black",
				font=ImageFont.truetype(bold, 23),
			)
		for i, line in enumerate(footer_lines):
			draw.text(
				(60, int((height - 18 - (len(footer_lines) - i) * 12) * factor)),
				line,
				fill="black",
				font=ImageFont.truetype(regular, 15),
			)
		image.save(self.out / (name + ".png"), dpi=(144, 144))
		self.figures.append(
			dict(
				figure=name,
				status="draft_generated",
				panels=len(panels),
				path=str(self.out / (name + ".png")),
				caption=caption,
			)
		)

	def make_figures(self):
		traits = self.arg.Y.split(",")
		layers = self.arg.biom.split(",")
		first = traits[0]
		# Fig1: column-based association, without pretending selected excerpts are a full atlas.
		f1 = []
		counts = []
		for tr in traits:
			for la in layers:
				d = self.get(tr, la, "incident")
				if has(d, "beta", "p.value"):
					z = d.copy()
					z["logp"] = -np.log10(num(z, "p.value").clip(lower=1e-300))
					z["bonferroni"] = num(z, "p.value") < 0.05 / len(z)
					counts.append(
						dict(
							trait=tr,
							layer=la,
							tested=len(z),
							Bonferroni=int(z.bonferroni.sum()),
						)
					)

					def volcano(ax, z):
						for label, g in z.groupby("bonferroni"):
							ax.scatter(
								num(g, "beta"),
								g.logp,
								s=8,
								alpha=0.55,
								label="Bonferroni" if label else "Other",
							)
						ax.axhline(
							-np.log10(0.05 / len(z)), linestyle="--", linewidth=0.7
						)
						ax.set_xlabel("Log hazard ratio per SD")
						ax.set_ylabel("-log10(P)")
						ax.legend(fontsize=7)

					f1.append(
						self.panel(
							f"Fig1_{tr}_{la}",
							f"{tr}: {la} incident associations",
							z,
							volcano,
							"Within-scan Bonferroni; no significance-based selection of the test set.",
							[f"{tr}.{la}.incident"],
						)
					)
		if not f1:

			def forest(ax, z):
				z = z.reset_index(drop=True)
				ax.errorbar(
					z.beta,
					np.arange(len(z)),
					xerr=[z.beta - z.lo, z.hi - z.beta],
					fmt="o",
					capsize=3,
				)
				ax.set_yticks(np.arange(len(z)), z.feature)
				ax.axvline(0, linestyle="--", linewidth=0.7)
				ax.set_xlabel("Measured-biomarker log HR (95% CI)")

			z = self.get(first, "prot", "pgs_components")
			if has(z, "model", "feature", "beta", "lo", "hi"):
				z = z[
					z.model.eq("measured")
					& z.feature.isin(["PCSK9", "LPA", "GDF15", "MMP12", "CCL19"])
				]
				f1.append(
					self.panel(
						"Fig1_anchor_excerpt",
						"CAD: measured biomarker anchors",
						z,
						forest,
						"Selected anchors in the matched PGS cohort; NOT a complete PWAS atlas.",
						[f"{first}.prot.pgs_components"],
					)
				)
			z = self.get("ra", "prot", "paired")
			if has(z, "model", "component", "feature", "beta", "std.error"):
				z = z[
					z.model.str.startswith("separate") & z.component.eq("Measured")
				].copy()
				z["lo"] = z.beta - 1.96 * z["std.error"]
				z["hi"] = z.beta + 1.96 * z["std.error"]
				f1.append(
					self.panel(
						"Fig1_RA_anchor_excerpt",
						"RA: measured biomarker anchors",
						z,
						forest,
						"Matched people within this RA analysis; distinct cohort from legacy prediction.",
						["ra.prot.paired"],
					)
				)
		# The observed-to-MR and MR-to-observed denominators answer different
		# questions. Plot them separately, retaining unavailable tests as a category.
		for direction in ["MR to observational", "Observational to MR"]:
			parts = []
			for tr in traits:
				for la in layers:
					d = self.get(tr, la, "mr_support")
					if has(d, "endpoint", "direction", "status", "n", "denominator"):
						q = d[
							d.endpoint.eq("incident") & d.direction.eq(direction)
						].copy()
						q["cohort"] = tr + " / " + la
						parts.append(q)
			if parts:
				z = pd.concat(parts, ignore_index=True)

				def support_matrix(ax, z):
					grid = z.pivot(index="cohort", columns="status", values="n")
					cols = [
						c
						for c in [
							"Concordant",
							"Discordant",
							"No FDR support",
							"Unavailable",
						]
						if c in grid
					]
					grid[cols].plot.barh(stacked=True, ax=ax)
					ax.set_xlabel(
						"Reported aggregate entries (not independent mechanisms)"
					)
					ax.set_ylabel("")
					ax.legend(fontsize=7)

				f1.append(
					self.panel(
						"Fig1_" + direction,
						"Incident: " + direction,
						z,
						support_matrix,
						"Unavailable MR is not a null causal result; correlated biomarkers and shared loci are not independent discoveries.",
						[f"{tr}.{la}.mr_support" for tr in traits for la in layers],
					)
				)
		if counts:
			pd.DataFrame(counts).to_csv(
				self.table_dir / "association_counts.csv", index=False
			)
		self.compose(
			"Fig1",
			"Incident molecular associations",
			f1,
			"Conventional association results establish context, not causality. In excerpt mode, only verified anchor rows are shown; no full-scan yield is inferred. Prevalent OR and incident HR estimate different quantities.",
		)
		# Fig2: source-dependent interpretation on the reported matched scale.
		f2 = []
		for la, anchors in [
			("prot", ["PCSK9", "LPA", "GDF15", "MMP12", "CCL19"]),
			("met", ["L_VLDL_TG.pct", "L_VLDL_TG", "Total_TG", "ApoB"]),
		]:
			z = self.get(first, la, "pgs_components")
			if has(z, "model", "term", "feature", "beta", "lo", "hi"):
				z = z[z.model.eq("joint") & z.feature.isin(anchors)].copy()
				if "landmark" in z:
					z = z[num(z, "landmark").eq(0)]

				def components(ax, z):
					order = list(dict.fromkeys(z.feature))
					pos = {v: i for i, v in enumerate(order)}
					for j, (term, g) in enumerate(z.groupby("term")):
						y = g.feature.map(pos).to_numpy() + (j - 0.5) * 0.18
						ax.errorbar(
							g.beta,
							y,
							xerr=[g.beta - g.lo, g.hi - g.beta],
							fmt="o",
							capsize=2,
							label=(
								"Calibrated genetic component"
								if term == ".G"
								else "Remaining component"
							),
						)
					ax.set_yticks(range(len(order)), order)
					ax.axvline(0, linestyle="--", linewidth=0.7)
					ax.set_xlabel("Log HR on the saved component scale")
					ax.legend(fontsize=7, loc="best")

				f2.append(
					self.panel(
						f"Fig2_{la}",
						f"CAD {la}: genetic and remaining components",
						z,
						components,
						"Same participants/covariates within each contrast; calibration uncertainty and original GWAS overlap still matter.",
						[f"{first}.{la}.pgs_components"],
					)
				)
		zlist = []
		for la in layers:
			z = self.get(first, la, "pgs_contrasts")
			if has(z, "feature", "beta_difference", "lo", "hi"):
				keep = [
					"GDF15",
					"MMP12",
					"LPA",
					"PCSK9",
					"L_VLDL_TG.pct",
					"L_VLDL_TG",
					"Total_TG",
					"ApoB",
				]
				zz = z[z.feature.isin(keep)].copy()
				if "landmark" in zz:
					zz = zz[num(zz, "landmark").eq(0)]
				zz["layer"] = la
				zlist.append(zz)
		if zlist:
			z = pd.concat(zlist, ignore_index=True)

			def contrast(ax, z):
				z = z.reset_index(drop=True)
				ax.errorbar(
					z.beta_difference,
					range(len(z)),
					xerr=[z.beta_difference - z.lo, z.hi - z.beta_difference],
					fmt="o",
					capsize=2,
				)
				ax.set_yticks(range(len(z)), z.feature)
				ax.axvline(0, linestyle="--", linewidth=0.7)
				ax.set_xlabel("Genetic minus remaining log-HR coefficient (95% CI)")

			f2.append(
				self.panel(
					"Fig2_contrasts",
					"CAD: formal within-biomarker contrasts",
					z,
					contrast,
					"Saved conditional covariance-based intervals. An absent genetic association does not prove reverse causation.",
					[f"{first}.{la}.pgs_contrasts" for la in layers],
				)
			)
		z = self.get("ra", "met", "coloc_by_locus")
		if has(z, "feature", "locus", "PP.H4", "PP.H4_robust_min"):

			def coloc_check(ax, z):
				ax.scatter(z["PP.H4"], z["PP.H4_robust_min"], s=28)
				for _, r in z[z.feature.isin(["Lactate", "Acetate"])].iterrows():
					ax.annotate(
						r.feature,
						(r["PP.H4"], r["PP.H4_robust_min"]),
						xytext=(5, -12),
						textcoords="offset points",
						fontsize=8,
					)
				if "policy_h4" in z:
					for cutoff in pd.to_numeric(z.policy_h4, errors="coerce").dropna().unique():
						ax.axhline(cutoff, linestyle="--", linewidth=0.7)
				ax.set_xlim(-0.03, 1.04)
				ax.set_ylim(-0.03, 1.02)
				ax.set_xlabel("Default-prior PP(H4)")
				ax.set_ylabel("Minimum PP(H4) across prior sensitivity")

			f2.append(
				self.panel(
					"Fig2_RA_coloc",
					"RA metabolites: colocalization filters the MR narrative",
					z,
					coloc_check,
					"Each row is a feature-locus test, not an independent locus. Conditional SNP posterior concentration is not the overall probability of H4.",
					["ra.met.coloc_by_locus"],
				)
			)
		self.compose(
			"Fig2",
			"Inherited propensity and measured state",
			f2,
			"Genetic and remaining components are associations, not pure causal and lifestyle components. Ratio and absolute triglyceride traits are distinct. GDF15/MMP12 calibration-refit bootstrap and source-specific score availability must be checked in the audit.",
		)
		# Fig3: supervision should be judged on measurable domain coverage as well as prediction.
		f3 = []
		z = self.get(first, "prot", "focus_pillars")
		if has(z, "cohort", "component", "n"):

			def pillars(ax, z):
				grid = z.pivot(index="component", columns="cohort", values="n")
				grid.plot.barh(ax=ax)
				ax.set_xlabel("Eligible replicated proxy associations")
				ax.set_ylabel("LE8 component")
				ax.legend(fontsize=8)

			f3.append(
				self.panel(
					"Fig3_pillars",
					"CAD protein proxies are not evenly distributed",
					z,
					pillars,
					"Eligible associations are not disjoint protein counts unless the input enforces exclusive assignment. Yang enters proxy learning only.",
					[f"{first}.prot.focus_pillars"],
				)
			)
		z = self.get(first, "prot", "effective_assays")
		if has(z, "model", "nominal", "effective"):
			z = z[
				z.model.str.endswith(("_5", "_10", "_50")) & z.coefficient_rows.gt(0)
			].copy()

			def effective(ax, z):
				z = z.set_index("model")[["nominal", "effective"]]
				z.plot.barh(ax=ax)
				ax.set_xlabel("Nominal versus numerically active assays")
				ax.set_ylabel("")
				ax.legend(fontsize=7)

			f3.append(
				self.panel(
					"Fig3_effective",
					"Assay count does not imply omics contribution",
					z,
					effective,
					"Effective means |coefficient| above the configured numerical threshold, not clinical assay utility.",
					[f"{first}.prot.focus_coefficients", f"{first}.prot.focus_members"],
					size=(7.2, max(4.7, 0.24 * len(z))),
				)
			)
		z = self.get(first, "prot", "focus_proxy")
		if has(z, "component", "model", "delta_R2"):
			# A fixed display subset, not a ranking by test R2.
			names = [
				"NS_10",
				"YS_Yin_10",
				"YS_YinYang_10",
				"YSplus_YinYang_10",
				"YSbalanced_YinYang_10",
			]
			z = z[z.model.isin(names)].copy()

			def proxy(ax, z):
				grid = z.pivot_table(
					index="component",
					columns="model",
					values="delta_R2",
					aggfunc="first",
				)
				im = ax.imshow(np.ma.masked_invalid(grid.to_numpy()), aspect="auto")
				ax.set_xticks(
					range(len(grid.columns)), grid.columns, rotation=35, ha="right"
				)
				ax.set_yticks(range(len(grid.index)), grid.index)
				ax.figure.colorbar(im, ax=ax, label="Held-out incremental R2")

			f3.append(
				self.panel(
					"Fig3_proxy",
					"Do ten assays reconstruct LE8 domains?",
					z,
					proxy,
					"Post-hoc OLS panel reconstruction beyond basic covariates; deployed concept fidelity is evaluated separately.",
					[f"{first}.prot.focus_proxy"],
				)
			)
		z = self.get(first, "prot", "focus_proxy")
		if has(z, "model", "component", "N", "R2_basic_omics"):
			ns = z[z.model.eq("NS_10")][["component", "N", "R2_basic_omics"]]
			ys = z[z.model.eq("YS_YinYang_10")][["component", "N", "R2_basic_omics"]]
			zz = ys.merge(
				ns, on="component", suffixes=("_YS", "_NS"), validate="one_to_one"
			)
			zz = zz[zz.N_YS.eq(zz.N_NS)].copy()
			zz["delta_R2"] = zz.R2_basic_omics_YS - zz.R2_basic_omics_NS
			if len(zz):

				def explain(ax, z):
					ax.scatter(z.delta_R2, range(len(z)))
					ax.axvline(0, linestyle="--")
					ax.set_yticks(range(len(z)), z.component)
					ax.set_xlabel("Held-out R²: YS YinYang minus NS (10 assays)")

				f3.append(
					self.panel(
						"Fig3_explainability",
						"Equal-budget LE8 information beyond basic covariates",
						zz,
						explain,
						"All available LE8 components, including negative differences. Point contrasts only; not a claim of intervention responsiveness.",
						[f"{first}.prot.focus_proxy"],
					)
				)
		z = self.get(first, "prot", "deployed_concept_fidelity")
		if has(z, "component", "model", "metric", "estimate", "status"):
			z = z[z.metric.eq("R2") & z.status.eq("ok")].copy()
			if len(z):
				def fidelity(ax, z):
					grid = z.pivot_table(index="component", columns="model", values="estimate", aggfunc="first")
					im = ax.imshow(np.ma.masked_invalid(grid.to_numpy()), aspect="auto")
					ax.set_xticks(range(len(grid.columns)), grid.columns, rotation=35, ha="right")
					ax.set_yticks(range(len(grid.index)), grid.index)
					ax.figure.colorbar(im, ax=ax, label="Deployed concept R²")
				f3.append(self.panel("Fig3_deployed_fidelity", "Do the deployed concepts predict measured LE8 domains?",
					z, fidelity, "Exact frozen concepts used by the risk models; R² uses the Yin training mean. RMSE, calibration and family-bootstrap intervals are in the fidelity table.",
					[f"{first}.prot.deployed_concept_fidelity"]))
		self.compose(
			"Fig3",
			"LE8 supervision, coverage and effective panel size",
			f3,
			"YS is ranked by replicated LE8 proxy strength without disease labels; NS uses multivariate disease selection and YSplus keeps the same total assay budget. Unsupported domains remain unavailable. Near-zero protein coefficients indicate a clinical-only fit, not a successful compact molecular panel. No universal YinYang gain is assumed.",
		)
		# Fig4: completed C4 validation, followed by separately available prediction protocols.
		f4 = []
		comparisons = []
		focus_used = joint_used = False
		wanted = [
			"C4 NS",
			"C4 YS",
			"C4 YSplus",
			"Pradeep-style / glmnet",
			"All-metabolite / glmnet",
			"Yu-style / LightGBM",
			"MWAS-ranked / LightGBM",
			"Evidence-selected compact",
		]
		for tr in traits:
			for la in layers:
				panels, data = self.c4_prediction_panels(tr, la)
				if panels:
					f4.extend(panels)
					comparisons.append(data)
					focus_used = True
					continue
		for tr in traits:
			z = self.get(tr, "joint", "joint_metrics")
			if has(z, "model", "budget", "landmark", "arm", "AUC", "N", "status"):
				joint_used = True
				# Primary model budget/landmark, with zero-assay clinical comparators.
				z = z[
					(num(z, "budget").isin([0, 10]))
					& num(z, "landmark").eq(5)
					& z.arm.eq("all_assays")
					& z.status.eq("ok")
				].copy()
				selected = [
					"Clinical",
					"Clinical_PRS",
					"Clinical_Protein_sharedBudget_NS",
					"Clinical_Metabolite_sharedBudget_NS",
					"Clinical_ProtMet_NS",
					"Clinical_ProtMet_PRS_NS",
					"Clinical_ProtMet_YS_YinYang",
					"Clinical_ProtMet_PRS_YS_YinYang",
				]
				z = z[z.model.isin(selected)]
				if z.N.nunique() > 1 or z.model.duplicated().any():
					self.flag(
						tr,
						"joint",
						"JOINT_DISPLAY_COHORT_MISMATCH",
						"blocking",
						"Primary slice has unequal N or duplicate models; no pooled figure drawn.",
						"joint_metrics",
					)
					z = pd.DataFrame()

				def joint_perf(ax, z):
					ax.scatter(z.AUC, range(len(z)))
					ax.set_yticks(range(len(z)), z.model)
					ax.set_xlabel("IPCW AUC, 5-year landmark to year 10; common cohort")

				f4.append(
					self.panel(
						"Fig4_" + tr + "_joint",
						tr + ": measured omics plus disease PRS",
						z,
						joint_perf,
						"Only completed primary-budget/landmark/all-assay results; no replacement primary selected from testing.",
						[f"{tr}.joint.joint_metrics"],
					)
				)
		if comparisons:
			pd.concat(comparisons, ignore_index=True).to_csv(
				self.table_dir / "prediction_comparison.csv", index=False
			)
		notes = [
			"Each disease/omics panel uses its own held-out cohort; protein and metabolite results are not a matched modality comparison."
		]
		if focus_used:
			notes.append(
				"C4 display: configured assay budget shown in the panel, all validation participants, baseline to year 10. Dashed lines show clinical AUC or zero paired gain versus NS. Intervals condition on frozen fits; exploratory comparisons require external validation. Death is censored; these are not competing-risk cumulative incidences."
			)
		if joint_used:
			notes.append(
				"Joint PRS/omics panels use their separate completed common-cohort protocol and five-year landmark."
			)
		self.compose(
			"Fig4",
			"Prediction and the cost of interpretability",
			f4,
			" ".join(notes),
		)
		# Current full-cohort twin models and controls.
		f5 = []
		for tr in traits:
			for la in layers:
				z = self.get(tr, la, 'abm_metrics')
				if not has(z, 'model', 'subset', 'AUC_IPCW'):
					continue
				z = z[z['subset'].eq('all')].copy()
				def plot_twin(ax, data):
					ax.barh(data.model, data.AUC_IPCW)
					ax.set(xlim=(0,1), xlabel='IPCW AUC, all eligible test participants')
				f5.append(self.panel('Fig5_' + tr + '_' + la, tr + ' / ' + la + ': twin models',
					z, plot_twin, ABM_BOUNDARY, [f'{tr}.{la}.abm_metrics']))
		self.compose('Fig5', 'Twin confidence and model controls', f5, ABM_BOUNDARY)
		# Additional C4 panels use the current completed results.
		s2 = []
		for tr in traits:
			z = self.get(tr, "prot", "focus_contrasts")
			if has(
				z, "stratum", "landmark", "model", "delta_AUC", "delta_lo", "delta_hi"
			):
				z = z[z.stratum.eq("All") & num(z, "landmark").eq(0)].copy()
				if "contrast" in z:
					z = z[z.contrast.eq("Added Yang for proxy learning")]

				def deltac(ax, z):
					z = z.reset_index(drop=True)
					ax.errorbar(
						z.delta_AUC,
						range(len(z)),
						xerr=[z.delta_AUC - z.delta_lo, z.delta_hi - z.delta_AUC],
						fmt="o",
						capsize=2,
					)
					ax.set_yticks(range(len(z)), z.model)
					ax.axvline(0, linestyle="--", linewidth=0.8)
					ax.set_xlabel("YinYang minus Yin IPCW AUC")

				s2.append(
					self.panel(
						"FigS2_" + tr,
						tr + ": added Yang at the same budget",
						z,
						deltac,
						"Post hoc internal paired bootstrap; zero intervals are audited for identical panels or null molecular coefficients.",
						[f"{tr}.prot.focus_contrasts"],
					)
				)
		self.compose(
			"FigS2",
			"Does Yang add predictive information?",
			s2,
			"All available saved same-budget Yang-versus-Yin comparisons are shown; zero intervals do not independently establish equivalence.",
		)
		s3 = []
		for tr in traits:
			z = self.get(tr, "prot", "cell_enrichment")
			if has(z, "model", "cell_type", "FDR_all_panel_cell_tests", "hits"):
				z = z.copy()
				z["logq"] = -np.log10(
					num(z, "FDR_all_panel_cell_tests").clip(lower=1e-300)
				)
				cells = z.groupby("cell_type").hits.sum().nlargest(16).index
				z = z[z.cell_type.isin(cells)]

				def cell(ax, z):
					grid = z.pivot_table(
						index="cell_type",
						columns="model",
						values="logq",
						aggfunc="mean",
					)
					im = ax.imshow(np.ma.masked_invalid(grid), aspect="auto")
					ax.set_xticks(
						range(len(grid.columns)), grid.columns, rotation=40, ha="right"
					)
					ax.set_yticks(range(len(grid.index)), grid.index)
					ax.figure.colorbar(im, ax=ax, label="-log10(global adjusted P)")

				s3.append(
					self.panel(
						"FigS3_" + tr,
						tr + ": external cell-expression enrichment",
						z,
						cell,
						"External expression relevance, not secretion origin; pooled fold display is descriptive only.",
						[f"{tr}.prot.cell_enrichment"],
					)
				)
		self.compose(
			"FigS3",
			"C5 cellular context",
			s3,
			"Cell-expression labels and native CIGMA have separate interpretations. Plasma protein abundance alone does not identify the secreting cell or establish a cellular-aging clock.",
		)

		s4 = []
		parts = []
		for tr in traits:
			for la in layers:
				z = self.get(tr, la, "coloc_audit")
				if has(z, "metric", "value"):
					z = z[
						z.metric.isin(
							[
								"tested loci",
								"robust shared-signal loci",
								"one-SNP posterior concentration",
							]
						)
					].copy()
					z["cohort"] = tr + " / " + la
					parts.append(z)
		if parts:
			z = pd.concat(parts, ignore_index=True)

			def coloc_counts(ax, z):
				z.pivot(index="cohort", columns="metric", values="value").plot.barh(
					ax=ax
				)
				ax.set_xlabel("Feature-locus entries")
				ax.set_ylabel("")
				ax.legend(fontsize=7)

			s4.append(
				self.panel(
					"FigS4_coloc",
					"Colocalization coverage and posterior audit",
					z,
					coloc_counts,
					"Categories overlap: concentration is a diagnostic, not proof of invalidity or causal confirmation.",
					[f"{tr}.{la}.coloc_audit" for tr in traits for la in layers],
				)
			)
		self.compose(
			"FigS4",
			"Genetic evidence needs locus-level checks",
			s4,
			"Colocalization counts cover prespecified anchors, MR-selected candidates and explicitly planned disease loci; these are not unique independent genomic regions. Posterior concentration on one SNP under H4 does not imply H4 is likely.",
		)

		s5 = []
		for tr in traits:
			z = self.get(tr, "joint", "PRS_omics_strata")
			if has(z, "stratum", "N", "observed_risk_ipcw", "landmark", "arm"):
				z = z[num(z, "landmark").eq(5) & z.arm.eq("all_assays")].copy()
				if z.stratum.duplicated().any():
					self.flag(
						tr,
						"joint",
						"DUPLICATE_PRS_OMIC_STRATA",
						"blocking",
						"Duplicate primary-slice strata; the display is suppressed rather than selecting a run.",
						"PRS_omics_strata",
					)
					z = pd.DataFrame()

				def strata_plot(ax, z):
					ax.barh(range(len(z)), z.observed_risk_ipcw)
					ax.set_yticks(
						range(len(z)),
						[f"{st} [N={int(n)}]" for st, n in zip(z.stratum, z.N)],
					)
					ax.set_xlabel("IPCW observed risk, landmark year 5 to year 10")

				s5.append(
					self.panel(
						"FigS5_" + tr,
						tr + ": disease PRS x measured omic score",
						z,
						strata_plot,
						"High status is defined from training-fold 75th percentiles, not held-out optimized cutoffs. No unsupported confidence intervals are invented.",
						[f"{tr}.joint.PRS_omics_strata"],
					)
				)
		self.compose(
			"FigS5",
			"Complementary inherited and measured risk dimensions",
			s5,
			"Only completed common-cohort primary-slice strata are displayed. These are cause-specific net-risk summaries with death censored, not competing-risk cumulative incidence or validated clinical thresholds.",
		)

	def finish(self):
		registry = pd.DataFrame(self.registry)
		registry.to_csv(self.out / "source_registry.csv", index=False)
		findings = pd.DataFrame(
			self.findings,
			columns=["trait", "layer", "code", "severity", "detail", "source_roles"],
		)
		findings.to_csv(self.out / "audit_findings.csv", index=False)
		pd.DataFrame(self.figures).to_csv(self.out / "figure_manifest.csv", index=False)
		pmeta = []
		for p in self.panels.values():
			pmeta.append(
				dict(
					panel=p.name,
					title=p.title,
					plot_data=str(p.csv),
					caption=p.caption,
					source_roles=";".join(p.source_roles),
					source_paths=";".join(
						self.source_files.get(k, "derived audit table")
						for k in p.source_roles
					),
				)
			)
		pd.DataFrame(pmeta).to_csv(self.out / "panel_manifest.csv", index=False)
		# Index, but never move or delete, mature source figures. Excludes private paths
		# and this output directory. These remain candidates, not automatically main figures.
		ext = []
		for tr in self.arg.Y.split(","):
			for folder in [self.root / tr, final_dir(self.root, tr)]:
				if not folder.exists():
					continue
				for path in sorted(folder.rglob("*")):
					if (
						path.is_file()
						and path.suffix.lower() in {".pdf", ".png", ".svg"}
						and not any(
							x.startswith("_") or x.startswith(".")
							for x in path.relative_to(folder).parts
						)
					):
						if self.out == path:
							continue
						ext.append(
							dict(
								trait=tr,
								source=str(path),
								bytes=path.stat().st_size,
								status="retained_source_figure_not_reselected",
							)
						)
		pd.DataFrame(ext, columns=["trait", "source", "bytes", "status"]).to_csv(
			self.out / "supplement_source_inventory.csv", index=False
		)
		claims = [
			(
				"Conventional prediction",
				"requires native same-cohort paired evidence",
				"Inspect current same-cohort paired contrasts; improvement and reproduction of published papers cannot be assumed.",
			),
			(
				"LE8 supervision",
				"evaluate prediction and explainability separately",
				"Current equal-budget paired prediction and held-out proxy reconstruction are reported separately in the question overview; incomplete pillar coverage and shrinkage collapse must be disclosed.",
			),
			(
				"YinYang superiority",
				"requires paired evidence",
				"Report the equal-budget paired contrasts and the actual fitted molecular coefficients.",
			),
			(
				"Inherited versus remaining state",
				"exploratory mechanistic interpretation",
				"Within-biomarker contrasts can differ; residual is not pure environment and PGS overlap is not removed by calibration cross-fitting.",
			),
			(
				"Non-inflammatory CAD",
				"not established",
				"Low baseline inflammation is a measured stratum, not a mechanistically validated CAD subtype.",
			),
			(
				"Individual predictability",
				"not established by readiness alone",
				"Primary-reference controls, permutation and risk-adjusted comparative error are required; retain negative results.",
			),
			(
				"Cellular origin/aging",
				"not established by annotation",
				"Expression enrichment is not secretion tracing, deconvolution or aging-clock replication.",
			),
			(
				"PRS + omics complementarity",
				"requires completed common-cohort validation",
				"Never combine separate prot/met cohorts into a purported matched improvement.",
			),
		]
		claim_order = ["Individual predictability", "Inherited versus remaining state", "LE8 supervision"]
		claims.sort(key=lambda row: claim_order.index(row[0]) if row[0] in claim_order else len(claim_order))
		pd.DataFrame(claims, columns=["claim", "maturity", "boundary"]).to_csv(
			self.out / "claim_register.csv", index=False
		)
		captions = "\n\n".join(
			"## " + f["figure"] + " — " + f["status"] + "\n" + f["caption"]
			for f in self.figures
		)
		(self.out / "CAPTIONS.md").write_text(captions)
		nblock = int(findings.severity.eq("blocking").sum()) if len(findings) else 0
		(self.out / "README.md").unlink(missing_ok=True)
		body = "<h1>LE8 5C: aggregate wrap-up</h1><p>" + html.escape(self.mode) + "</p>"
		body += '<p><a href="index.html">研究问题总览：C1 Individual predictability → C2 Genetic vs acquired/state → C4 LE8-supervised explainability</a></p>'
		body += (
			"<p><b>Internal draft. No new model fitting.</b> "
			+ str(nblock)
			+ " blocking interpretation/audit flags.</p>"
		)
		body += "<h2>Claim register</h2>" + pd.DataFrame(
			claims, columns=["Claim", "Maturity", "Boundary"]
		).to_html(index=False, escape=True)
		body += "<h2>Audit findings</h2>" + findings.to_html(index=False, escape=True)
		for figure in self.figures:
			if figure["status"] != "draft_generated":
				continue
			body += (
				"<h2>"
				+ html.escape(figure["figure"])
				+ '</h2><img src="'
				+ html.escape(Path(figure["path"]).name)
				+ '"><p>'
				+ html.escape(figure["caption"])
				+ "</p>"
			)
		shutil.rmtree(self.panel_dir)
		css = "body{font:15px Arial,sans-serif;max-width:1200px;margin:30px auto;line-height:1.5}table{border-collapse:collapse;width:100%;font-size:12px}td,th{padding:8px;border:1px solid #ccc;text-align:left}img{max-width:850px;width:100%}"
		(self.out / "report.html").write_text(
			'<!doctype html><html><meta charset="utf-8"><style>'
			+ css
			+ "</style><body>"
			+ body
			+ "</body></html>"
		)
		print(f"Final report: {self.out}; blocking interpretation flags={nblock}")
		if self.arg.strict and nblock:
			raise ValueError(
				"Strict audit failed; outputs retained for inspection. See audit_findings.csv."
			)


# Systematic evidence tables and cell-annotation inputs.
def bh(p, family_size=None):
	p = np.asarray(p, dtype=float)
	ok = np.isfinite(p)
	out = np.full(len(p), np.nan)
	if not ok.any():
		return out
	n = max(len(p), int(family_size or 0))
	ids = np.flatnonzero(ok)
	order = np.argsort(p[ok], kind="stable")
	adjusted = p[ids[order]] * n / np.arange(1, len(ids) + 1)
	out[ids[order]] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
	return out


def read(path):
	if not path.exists():
		return pd.DataFrame()
	try:
		return pd.read_csv(path)
	except pd.errors.EmptyDataError:
		return pd.DataFrame()


def aggregate(root):
	inventory, summary, concordance, windows, trajectories, pillars, contrasts = (
		[],
		[],
		[],
		[],
		[],
		[],
		[],
	)
	for layer, prefix in (("prot", "pwas"), ("met", "mwas")):
		for module in (
			"c1_correlate",
			"c2_cause",
			"c3_coloc",
			"c4_connect",
			"c1_prediction",
		):
			directory = (
				root / layer / "c1_correlate" / "prediction"
				if module == "c1_prediction"
				else root / layer / module
			)
			files = sorted(directory.glob("*.csv"))
			if not files:
				inventory.append(
					dict(
						layer=layer,
						module=module,
						file="",
						rows=0,
						status="unavailable",
					)
				)
			for f in files:
				d = read(f)
				inventory.append(
					dict(
						layer=layer,
						module=module,
						file=str(f.relative_to(root.parent)),
						rows=len(d),
						status="available" if len(d) else "empty",
						sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
					)
				)
		c1 = root / layer / "c1_correlate"
		incident = read(c1 / f"{prefix}_incident_adj2.csv")
		pgs = read(c1 / f"{prefix}_pgs_incident_full_genetic.csv")
		for kind, d in (("measured", incident), ("biomarker_PGS", pgs)):
			if {"term", "p.value", "FDR"} <= set(d):
				summary.append(
					dict(
						layer=layer,
						analysis=kind,
						tested=len(d),
						finite_p=int(d["p.value"].notna().sum()),
						FDR05=int((d.FDR < 0.05).sum()),
						N_min=d.N_total.min(),
						N_max=d.N_total.max(),
						events_min=d.N_event.min(),
						events_max=d.N_event.max(),
					)
				)
		if len(incident) and len(pgs):
			cols = ["term", "beta", "std.error", "FDR", "N_total", "N_event"]
			z = incident[cols].merge(
				pgs[cols], on="term", how="outer", suffixes=("_measured", "_PGS")
			)
			z["layer"] = layer
			z["both_supported"] = (z.FDR_measured < 0.05) & (z.FDR_PGS < 0.05)
			z["same_direction"] = np.sign(z.beta_measured) == np.sign(z.beta_PGS)
			z["interpretation"] = (
				"Parallel associations; different SD units/covariates/cohorts; not MR or a causal test"
			)
			concordance.append(z)
		lm = read(c1 / f"{prefix}_incident_landmark_adj2.csv")
		if len(lm) and len(incident):
			# Legacy scans selected about 500 proteins using the same outcomes.
			# Retain original FDR and add the full assay x landmark family.
			lm["FDR_assay_landmark_family"] = bh(
				lm["p.value"], len(incident) * lm.landmark_years.nunique()
			)
			lm["subset_selected_before_landmark"] = lm.term.nunique() < len(incident)
			lm["layer"] = layer
			trajectories.append(lm)
		rw = read(c1 / f"{prefix}_diagnosis_window_riskset_adj2.csv")
		if len(rw):
			z = (
				rw.groupby(["side", "window_lo", "window_hi"], dropna=False)
				.agg(
					events_min=("N_event", "min"),
					events_max=("N_event", "max"),
					N_min=("N_total", "min"),
					N_max=("N_total", "max"),
				)
				.reset_index()
			)
			z["layer"] = layer
			z["warning"] = np.where(
				(z.side == "Post-baseline incident") & (z.events_max == 0),
				"No observed events; verify registry coverage before interpreting lead time",
				"",
			)
			windows.append(z)
		p = read(root / layer / "c4_connect/c4.focus.pillar_counts.csv")
		if len(p):
			p["layer"] = layer
			pillars.append(p)
		c = read(root / layer / "c4_connect/c4.focus.contrasts.csv")
		if len(c):
			c["layer"] = layer
			c["nominal_interval_excludes_zero"] = (c.delta_lo > 0) | (c.delta_hi < 0)
			c["interpretation"] = (
				"Exploratory; all tested contrasts retained; no subgroup-heterogeneity claim from separate CIs"
			)
			contrasts.append(c)
	combine = lambda xs: pd.concat(xs, ignore_index=True) if xs else pd.DataFrame()
	return dict(
		inventory=pd.DataFrame(inventory),
		c1_summary=pd.DataFrame(summary),
		measured_PGS=combine(concordance),
		event_windows=combine(windows),
		landmark_family=combine(trajectories),
		pillar_support=combine(pillars),
		all_c4_contrasts=combine(contrasts),
	)


def cell_annotation_inputs(root, output):
	"""Rebuild cell-annotation inputs from native saved associations/panels."""
	c1 = root / "prot/c1_correlate"
	d = read(c1 / "pwas_incident_adj2.csv")
	if d.empty or not {"term", "FDR"} <= set(d):
		return None
	universe = pd.DataFrame(
		{"assay": d.term, "gene": d.term.str.upper().replace({"NTPROBNP": "NPPB"})}
	)
	panels = [
		pd.DataFrame(
			{"feature": d.loc[d.FDR < 0.05, "term"], "model": "C1_incident_FDR05"}
		)
	]
	lm = read(c1 / "pwas_incident_landmark_adj2.csv")
	if len(lm):
		lm["q_full_family"] = bh(lm["p.value"], len(d) * lm.landmark_years.nunique())
		z = lm[(lm.landmark_years == 5) & (lm.q_full_family < 0.05)]
		panels.append(
			pd.DataFrame({"feature": z.term, "model": "C1_landmark5_full_family_FDR05"})
		)
	f = read(root / "prot/c4_connect/c4.focus.panel_members.csv")
	if {"feature", "model"} <= set(f):
		panels.append(f[["feature", "model"]].dropna())
	output.mkdir(parents=True, exist_ok=True)
	universe.to_csv(output / "c5.systematic.cell_universe.csv", index=False)
	pd.concat(panels, ignore_index=True).to_csv(
		output / "c5.systematic.cell_panels.csv", index=False
	)
	return (
		output / "c5.systematic.cell_universe.csv",
		output / "c5.systematic.cell_panels.csv",
	)


def cell_annotation(root, code_dir, outdir=None):
	"""Use the full assayed universe, not a significant-only background."""
	output = outdir or root
	paths = cell_annotation_inputs(root, output)
	if paths is None:
		return "C1 protein associations unavailable"
	cell = load_module("c5.cellulation.py")
	annotate, default_cell_atlas = cell.annotate, cell.default_cell_atlas

	annotate(
		paths[0],
		default_cell_atlas(),
		paths[1],
		output,
		"c5.systematic.cell",
	)
	return "External CellAge cell labels; full-assay background; no inferred tissue of release"


def plots(t, root, trait=None):
	trait = trait or root.name
	import matplotlib

	matplotlib.use("Agg")
	import matplotlib.pyplot as plt

	plt.rcParams.update(
		{
			"font.size": 9,
			"axes.spines.top": False,
			"axes.spines.right": False,
			"pdf.fonttype": 42,
			"svg.fonttype": "none",
		}
	)
	fig, axs = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
	for ax, title in zip(
		axs.flat,
		[
			"A  Tested associations",
			"B  Measured vs biomarker PGS",
			"C  Event ascertainment windows",
			"D  Supported LE8 domains",
			"E  Overall / subgroup comparisons",
			"F  Distal association support",
		],
	):
		ax.set_title(title, loc="left", fontweight="bold")
	d = t["c1_summary"]
	if len(d):
		labels = d.layer + " / " + d.analysis
		axs[0, 0].barh(labels, d.tested, color="#dce3eb", label="Tested")
		axs[0, 0].barh(labels, d.FDR05, color="#287a78", label="FDR < 0.05")
		axs[0, 0].legend(frameon=False)
	d = t["measured_PGS"]
	if len(d):
		for layer, color in (("prot", "#8d64ab"), ("met", "#df8d35")):
			z = d[d.layer == layer]
			axs[0, 1].scatter(
				z.beta_measured, z.beta_PGS, s=9, alpha=0.4, color=color, label=layer
			)
		z = d[d.both_supported].sort_values("FDR_PGS").head(3)
		if trait == "amr":
			z = pd.concat([z, d[d.term == "GlycA"]]).drop_duplicates("term")
		for _, row in z.iterrows():
			axs[0, 1].annotate(
				row.term,
				(row.beta_measured, row.beta_PGS),
				xytext=(4, 6),
				textcoords="offset points",
				fontsize=7,
			)
		axs[0, 1].axhline(0, color="grey", lw=0.7)
		axs[0, 1].axvline(0, color="grey", lw=0.7)
		axs[0, 1].set(
			xlabel="Measured log HR / measured SD", ylabel="PGS log HR / PGS SD"
		)
		axs[0, 1].legend(frameon=False)
	d = t["event_windows"]
	if len(d):
		for layer in d.layer.unique():
			z = d[(d.layer == layer) & (d.side == "Post-baseline incident")]
			axs[0, 2].plot(z.window_hi, z.events_max, "o-", label=layer)
		axs[0, 2].set(
			xlabel="Window upper bound (years)", ylabel="Events in window (feature max)"
		)
		axs[0, 2].legend(frameon=False)
	d = t["pillar_support"]
	if len(d):
		z = d.pivot_table(
			index="component", columns=["layer", "cohort"], values="n", aggfunc="first"
		).fillna(0)
		im = axs[1, 0].imshow(np.log1p(z), aspect="auto", cmap="Blues")
		axs[1, 0].set_yticks(range(len(z)), z.index)
		axs[1, 0].set_xticks(
			range(len(z.columns)),
			[" / ".join(x) for x in z.columns],
			rotation=25,
			ha="right",
		)
		for i in range(len(z)):
			for j in range(len(z.columns)):
				axs[1, 0].text(
					j,
					i,
					str(int(z.iloc[i, j])),
					ha="center",
					va="center",
					fontsize=8,
					color="white" if np.log1p(z.iloc[i, j]) > 4 else "black",
				)
	d = t["all_c4_contrasts"]
	if len(d):
		z = d[
			(d.layer == "prot")
			& (d.model == "YSplus_YinYang_50")
			& (d.reference == "NS_50")
		]
		labels = (
			z.stratum.str.replace(" baseline inflammation", "")
			+ " / L"
			+ z.landmark.astype(str)
		)
		axs[1, 1].errorbar(
			z.delta_AUC,
			range(len(z)),
			xerr=[
				np.maximum(0, z.delta_AUC - z.delta_lo),
				np.maximum(0, z.delta_hi - z.delta_AUC),
			],
			fmt="o",
			color="#287a78",
		)
		axs[1, 1].set_yticks(range(len(z)), labels, fontsize=7)
		axs[1, 1].axvline(0, color="grey", lw=0.7)
		axs[1, 1].set_xlabel("YSplus YY vs NS; paired delta AUC (50 assays)")
	d = t["landmark_family"]
	if len(d):
		for layer in d.layer.unique():
			z = (
				d[d.layer == layer]
				.groupby("landmark_years")["FDR_assay_landmark_family"]
				.apply(lambda x: (x < 0.05).sum())
			)
			axs[1, 2].plot(z.index, z.values, "o-", label=layer)
		axs[1, 2].set(
			xlabel="Landmark (years)", ylabel="FDR < .05, full assay × landmark family"
		)
		axs[1, 2].legend(frameon=False)
	# A C1-only trait gets meaningful C1 panels in the available space.
	if t["pillar_support"].empty and len(t["measured_PGS"]):
		for ax, layer in ((axs[1, 0], "prot"), (axs[1, 1], "met")):
			z = t["measured_PGS"].query("layer == @layer").nsmallest(8, "FDR_measured")
			pos = np.arange(len(z))
			ax.errorbar(
				z.beta_measured,
				pos,
				xerr=1.96 * z["std.error_measured"],
				fmt="o",
				color="#287a78",
				label="Measured",
			)
			ax.errorbar(
				z.beta_PGS,
				pos + 0.18,
				xerr=1.96 * z["std.error_PGS"],
				fmt="s",
				color="#a276b4",
				label="Matched PGS",
			)
			ax.set_yticks(pos, z.term, fontsize=8)
			ax.axvline(0, color="grey", lw=0.7)
			ax.set_title(
				("D  " if layer == "prot" else "E  ")
				+ layer
				+ " leading measured signals",
				loc="left",
				fontweight="bold",
			)
			ax.set_xlabel("Log HR per respective SD; different cohorts/covariates")
			ax.legend(frameon=False, fontsize=8)
	for ax in axs.flat:
		if not ax.has_data():
			ax.text(
				0.5,
				0.5,
				"Not available for this trait",
				transform=ax.transAxes,
				ha="center",
			)
	fig.suptitle(
		f"{trait}: systematic 5C evidence audit — descriptive, existing results",
		fontsize=14,
	)
	fig.savefig(root / "Fig1.evidence_atlas.png", dpi=220)
	plt.close(fig)


# Question-led synthesis and dashboard figures.
PILLARS = {
	"diet.pts": "饮食",
	"pa.pts": "身体活动",
	"smoke.pts": "烟草暴露",
	"sleep.pts": "睡眠",
	"bmi.pts": "BMI",
	"nonhdl.pts": "血脂",
	"hba1c.pts": "血糖",
	"bp.pts": "血压",
}


def truth(s):
	return s.astype(str).str.lower().isin(["true", "1", "t"])


def write_question_html(out, overview):
	import html

	page = overview.to_html(index=False, escape=True, border=0)
	content = (
		'<!doctype html><html lang="zh"><meta charset="utf-8"><title>LE8 研究问题</title>'
		"<style>body{font:16px system-ui;max-width:1500px;margin:40px auto;color:#193142}td,th{padding:12px;text-align:left;border-bottom:1px solid #d8e2e8}th{background:#eef5f6}</style>"
		"<h1>C1 Individual predictability → C2 Genetic vs acquired/state → C4 LE8-supervised explainability</h1>"
		"<p>C3 locus support and C5 cellular context are supporting evidence. All summaries are recalculated from current aggregate results; unavailable evidence is not a null finding.</p>"
	)
	modules = out / "final.questions.module_summary.csv"
	if modules.exists():
		content += pd.read_csv(modules).to_html(index=False, escape=True, border=0)
	content += "<h2>Result-driven claims</h2>" + page
	for name, title in [("abm_summary", "C1 · Twin full-cohort comparison and confidence"),
		("deviation_status", "C1 · Molecular deviation availability"),
		("met_pgs_concordance", "C1 · Metabolite inherited–measured concordance"),
		("met_temporal_counts", "C1 · Metabolite temporal classes"),
		("met_anchor", "C1 · L_VLDL_TG.pct prespecified anchor and burden/size comparators"),
		("met_pgs_conditional", "C1 · Conditional PGS sensitivity"),
		("met_measured_conditional", "C1 · Measured sensitivity"),
		("met_trajectory", "C1 · Risk-set trajectory"),
		("decomposition", "C2 · Genetic and remaining components"),
		("c4_endpoints", "C4 · Prediction and explainability are separate endpoints")]:
		path = out / f"final.questions.{name}.csv"
		if path.exists():
			d = pd.read_csv(path)
			content += f'<h2>{html.escape(title)}</h2><p>Complete aggregate tables are retained in the Final workbooks and Shiny downloads.</p>'
			if name.startswith("met_"):
				content += "<p>" + html.escape(MET_BOUNDARY + " " + ANCHOR_BOUNDARY) + "</p>"
			if name == "decomposition":
				content += f"<p>{len(d)} rows in the complete table; the first 50 rows are previewed below in source order.</p>"
				cols = [c for c in ("Y", "layer", "feature", "status", "incident_beta_observed", "incident_beta_genetic", "incident_beta_residual", "incident_beta_difference", "incident_difference_lo", "incident_difference_hi", "FDR_incident_p_difference", "prevalent_beta_difference", "FDR_prevalent_p_difference", "causal_fraction_identifiable") if c in d]
				d = d.loc[:, cols].head(50)
			content += d.drop(columns=[c for c in d if c.startswith("source_")], errors="ignore").to_html(index=False, escape=True, border=0)
	manifest = out / "final.questions.figure_manifest.csv"
	if manifest.exists():
		for _, r in pd.read_csv(manifest).iterrows():
			hashes = json.loads(r.get("source_hashes", "{}"))
			current = bool(hashes) and all(
				re.fullmatch(r"[a-z_]+", name)
				and (out / f"final.questions.{name}.csv").exists()
				and digest(out / f"final.questions.{name}.csv") == sha
				for name, sha in hashes.items()
			)
			if current and re.fullmatch(r"Fig\d+[.]question_[A-Za-z_]+", r.figure):
				content += (
					"<h2>"
					+ html.escape(r.figure)
					+ '</h2><img style="width:100%" src="'
					+ r.figure
					+ '.png"><p>'
					+ html.escape(r.caption)
					+ "</p>"
				)
			else:
				content += (
					"<p>"
					+ html.escape(str(r.figure))
					+ "：源表已更新或图版本未核实；运行 final 重新生成图。</p>"
				)
	(out / "index.html").write_text(content + "</html>", encoding="utf-8")


def interval_status(estimate, lo, hi):
	if not all(np.isfinite(x) for x in (estimate, lo, hi)):
		return "区间缺失"
	if lo > hi or estimate < lo - 1e-8 or estimate > hi + 1e-8:
		return "区间需核查"
	if estimate == lo == hi == 0:
		return "相同预测 / 零差异"
	if lo > 0:
		return "名义区间支持增加"
	if hi < 0:
		return "名义区间支持降低"
	return "区间包含零"


def panel_diagnostics(members, coefficients):
	"""Count measured assays and numerical coefficient use separately."""
	if not {"model", "feature"} <= set(members):
		return pd.DataFrame(
			columns=["model", "panel_assays", "effective_assays", "panel_hash"]
		)
	rows = []
	import hashlib

	for model, z in members.groupby("model", sort=False):
		fs = set(z.feature.dropna().astype(str))
		effective = np.nan
		if {"model", "variable", "beta"} <= set(coefficients) and fs:
			c = coefficients[coefficients.model.eq(model)]
			cc = c[c.variable.isin(fs)]
			if set(cc.variable) == fs:
				effective = int((cc.beta.abs() > 1e-8).sum())
		if not fs:
			effective = 0
		rows.append(
			dict(
				model=model,
				panel_assays=len(fs),
				effective_assays=effective,
				panel_hash=hashlib.sha256("\n".join(sorted(fs)).encode()).hexdigest(),
			)
		)
	return pd.DataFrame(rows)


def attach_panel_contrasts(contrasts, metrics, diagnostics):
	if contrasts.empty or not {
		"model",
		"reference",
		"stratum",
		"landmark",
		"horizon",
		"delta_AUC",
		"delta_lo",
		"delta_hi",
	} <= set(contrasts):
		return pd.DataFrame()
	d = contrasts.copy()
	keys = ["model", "stratum", "landmark", "horizon"]
	fields = [
		c
		for c in ["budget", "actual_assays", "N", "events_by_horizon", "status"]
		if c in metrics
	]
	if not set(keys) <= set(metrics) or metrics.duplicated(keys).any():
		d["comparison_valid"] = False
		d["interpretation"] = "Missing/duplicate model metrics; comparison unverified"
		return d
	for ref in (False, True):
		suffix = "_reference" if ref else ""
		m = metrics[keys + fields].rename(columns={c: c + suffix for c in fields})
		m = m.rename(columns={"model": "reference"}) if ref else m
		d = d.merge(
			m,
			on=["reference" if ref else "model", *keys[1:]],
			how="left",
			validate="many_to_one",
		)
		dg = diagnostics.rename(
			columns={
				c: ("reference" if ref else "model") if c == "model" else c + suffix
				for c in diagnostics
			}
		)
		d = d.merge(
			dg, on="reference" if ref else "model", how="left", validate="many_to_one"
		)
	d["comparison_valid"] = True
	for c in ["budget", "actual_assays", "N", "events_by_horizon"]:
		if c not in d or c + "_reference" not in d:
			d["comparison_valid"] = False
		else:
			d["comparison_valid"] &= d[c].notna() & d[c].eq(d[c + "_reference"])
	if "status" in d and "status_reference" in d:
		d["comparison_valid"] &= d.status.eq("ok") & d.status_reference.eq("ok")
	else:
		d["comparison_valid"] = False
	d["identical_panel"] = d.panel_hash.notna() & d.panel_hash.eq(
		d.panel_hash_reference
	)
	d["clinical_only"] = d.effective_assays.eq(0) | d.effective_assays_reference.eq(0)
	d["interval_status"] = [
		interval_status(*r) for r in d[["delta_AUC", "delta_lo", "delta_hi"]].to_numpy()
	]
	d.loc[~d.comparison_valid, "interval_status"] = "比较条件未核实"
	d["interpretation"] = (
		"Same native validation protocol; paired conditional bootstrap; nominal intervals across all budgets/strata, not multiplicity-adjusted confirmation"
	)
	return d




class Questions:
	def __init__(self, root, traits, layers, max_bytes=128 * 1024**2):
		self.root, self.traits, self.layers = Path(root), traits, layers
		self.max_bytes = max_bytes
		self.sources, self.tables, self.overview = [], {}, []

	def read(self, path, Y, layer, role):
		p = Path(path)
		root = self.root if p.is_relative_to(self.root) else abm_source_root(self.root)
		rel = str(p.relative_to(root))
		record = dict(
			Y=Y,
			layer=layer,
			role=role,
			source_file=rel,
			sha256="",
			rows=0,
			status="missing",
		)
		if p.exists():
			if not safe_path(p, root) or p.stat().st_size > self.max_bytes:
				raise ValueError(f"Unsafe or oversized question source: {rel}")
			record["sha256"] = digest(p)
			try:
				d = pd.read_csv(p, sep="\t" if ".tsv" in p.name else ",")
			except pd.errors.EmptyDataError:
				d = pd.DataFrame()
			if _index.has_private_columns(d.columns):
				raise ValueError(f"Participant data rejected: {rel}")
			record.update(rows=len(d), status="available" if len(d) else "empty")
		else:
			d = pd.DataFrame()
		self.sources.append(record)
		if len(d):
			d = d.copy()
			d["Y"], d["layer"] = Y, layer
			d["source_file"], d["source_sha256"] = rel, record["sha256"]
			d["source_row"] = np.arange(1, len(d) + 1)
		return d

	def add(self, name, d):
		if len(d):
			self.tables.setdefault(name, []).append(d)

	def inherited_sources(self, d, Y, layer, role):
		"""Retain provenance for frames already normalized by the evidence index."""
		if d.empty:
			return d
		for col in ["source_file", "source_row", "source_id"]:
			if col not in d and "_" + col in d:
				d[col] = d["_" + col]
		if "source_file" in d:
			hashes = {}
			for source in d.source_file.dropna().unique():
				self.read(self.root / source, Y, layer, role)
				hashes[source] = self.sources[-1]["sha256"]
			d["source_sha256"] = d.source_file.map(hashes)
		return d

	def claim(self, Y, layer, question, status, result, limitation, evidence):
		module = ("C1" if question in {"ABM 分层", "Twin confidence",
					"Molecular deviation", "Metabolite inherited–measured concordance", "L_VLDL_TG.pct illustrative anchor"}
			else "C2" if question == "遗传与实测差异" else "C5" if question == "Cellulation" else "C4")
		self.overview.append(
			dict(
				Y=Y,
				layer=layer,
				question=question,
				module=module,
				status=status,
				result=result,
				limitation=limitation,
				evidence=evidence,
			)
		)

	def layer(self, Y, layer, candidates, proxies):
		root = self.root / Y / layer
		read = lambda module, name, role: self.read(
			root / module / name, Y, layer, role
		)
		c4 = lambda name: read("c4_connect", "c4.focus." + name + ".csv", name)
		metrics, members, coef = (
			c4("metrics"),
			c4("panel_members"),
			c4("model_coefficients"),
		)
		diagnostics = panel_diagnostics(members, coef)
		ct = attach_panel_contrasts(c4("contrasts"), metrics, diagnostics)
		if len(metrics) and "model" in metrics:
			metrics = metrics.merge(
				diagnostics, on="model", how="left", validate="many_to_one"
			)
		heterogeneity = c4("heterogeneity")
		for name, d in [
			("prediction", metrics),
			("contrasts", ct),
			("members", members),
			("pillars", c4("panel_coverage")),
			("heterogeneity", heterogeneity),
			("inflammation_definition", c4("inflammation_definition")),
			("design", c4("design")),
			("fit", c4("fit_diagnostics")),
			("concept_coefficients", c4("concept_coefficients")),
			("concept_status", c4("concept_status")),
			("concept_fold_panels", c4("concept_fold_panels")),
		]:
			self.add(name, d)
		if len(heterogeneity) and {"landmark", "p_heterogeneity", "FDR"} <= set(
			heterogeneity
		):
			h = heterogeneity[heterogeneity.landmark.eq(0)]
			self.claim(
				Y,
				layer,
				"炎症分层异质性",
				"探索性检验",
				f"基线至10年：{int(h.p_heterogeneity.notna().sum())} 个可估计组间增益差异，{int((h.FDR < 0.05).sum())} 个 FDR<0.05。",
				"不同组各自显著与否不能替代直接异质性检验；低基线炎症不是已确立的疾病亚型。",
				"heterogeneity",
			)
		else:
			self.claim(
				Y,
				layer,
				"炎症分层异质性",
				"未提供分层结果",
				"此组未提供可用的低/高炎症组间比较。",
				"不能由未检验推断无异质性。",
				"heterogeneity",
			)
		primary = (
			ct[
				ct.stratum.eq("All")
				& ct.landmark.eq(0)
				& ct.reference.str.match(r"^NS_\d+$")
			]
			if len(ct)
			else pd.DataFrame()
		)
		if len(primary) and "comparison_valid" in primary:
			primary = primary[primary.comparison_valid]
		if len(primary):
			pos, neg = (
				int((primary.delta_lo > 0).sum()),
				int((primary.delta_hi < 0).sum()),
			)
			self.claim(
				Y,
				layer,
				"LE8 预测增益",
				"supported" if pos == len(primary) else "no universal superiority",
				("No universal prediction superiority was observed. " if pos < len(primary) else "") + f"基线至10年、全验证人群：{len(primary)} 个同预算配对 ΔAUC 比较，{pos} 个名义区间完全大于0，{neg} 个完全小于0。",
				f"{int(primary.clinical_only.sum())} 个比较至少一方分子系数趋近零；全部预算/方法均展示，未据测试表现选最佳。",
				"contrasts",
			)
		else:
			self.claim(
				Y,
				layer,
				"LE8 预测增益",
				"结果不足",
				"缺少可核实的同预算配对比较。",
				"需要相同测试人群、终点、预算及配对区间。",
				"contrasts",
			)
		proxy = (
			proxies[(proxies.Y == Y) & (proxies.layer == layer)].copy()
			if len(proxies)
			else pd.DataFrame()
		)
		proxy = self.inherited_sources(proxy, Y, layer, "proxy_point_source")
		if len(proxy) and "lo" in proxy and proxy.lo.notna().any():
			self.read(
				root / "c4_connect/c4.explain.contrasts.csv",
				Y,
				layer,
				"proxy_interval_source",
			)
			proxy["interval_source_file"] = self.sources[-1]["source_file"]
			proxy["interval_source_sha256"] = self.sources[-1]["sha256"]
		if len(proxy):
			proxy["pillar"] = proxy.component.map(PILLARS).fillna(proxy.component)
			pos = int((proxy.delta_R2_vs_NS > 0).sum())
			confirmed = int(
				(
					(
						proxy.get(
							"FDR_all_proxy_contrasts",
							pd.Series(np.nan, index=proxy.index),
						)
						< 0.05
					)
					& (proxy.get("lo", pd.Series(np.nan, index=proxy.index)) > 0)
				).sum()
			)
			self.claim(
				Y,
				layer,
				"LE8 可解释性",
				"已有重建结果",
				f"{pos}/{len(proxy)} 个 panel×component 比较的 ΔR²>0；{confirmed} 个同时满足配对区间>0和全比较 FDR<0.05。",
				"预测LE8成分是代理保真度证据；不等于生活方式干预响应。缺失区间不判为显著。",
				"proxy",
			)
		self.add("proxy", proxy)
		yy = (
			ct[
				ct.stratum.eq("All")
				& ct.landmark.eq(0)
				& ct.reference.str.contains("_Yin_")
			]
			if len(ct)
			else pd.DataFrame()
		)
		self.claim(
			Y,
			layer,
			"Yin–Yang 增益",
			"已有配对比较" if len(yy) else "结果不足",
			f"{len(yy)} 个 YinYang–Yin 比较；{int((yy.delta_lo > 0).sum()) if len(yy) else 0} 个名义区间支持增加。",
			"Prevalent donors 仅进入 proxy discovery；incident 风险验证保持独立。增益本身不确定因果方向。",
			"contrasts",
		)

		for name, path in abm_table_paths(self.root, Y, layer).items():
			d = self.read(path, Y, layer, name)
			if len(d):
				d['backend'] = 'twin' if name.startswith('abm_') else 'molecular_deviation'
				self.add(name, d)
		genetic = (
			candidates[(candidates.Y == Y) & (candidates.layer == layer)].copy()
			if len(candidates)
			else pd.DataFrame()
		)
		genetic = self.inherited_sources(genetic, Y, layer, "genetic_catalogue_source")
		if len(genetic):
			matched = genetic[truth(genetic.same_people)]
			nullsig = (matched.measured_p >= 0.05) & (matched.pgs_FDR < 0.05)
			diffs = matched.GR_difference_FDR < 0.05
			self.claim(
				Y,
				layer,
				"遗传与实测差异",
				"存在不同关联信息"
				if nullsig.any() or diffs.any()
				else "尚无充分差异证据",
				f"同人群 {len(matched)} 个 assay/scope 记录：{int(nullsig.sum())} 个实测 P≥0.05、PGS FDR<0.05；{matched.GR_difference_FDR.notna().sum()} 个有匹配 G−R 检验，其中 {int(diffs.sum())} 个 FDR<0.05。C2 完整分解见独立结果表。",
				"PGS 是遗传倾向，不是出生时浓度；R包含未捕获遗传、环境、疾病、治疗和误差。差异不是因果证明。",
				"genetic",
			)
		self.add("genetic", genetic)
		for name, module, filename in [
			*[(name, "c1_correlate", filename) for name, filename in MET_TABLES.items() if layer == "met"],
			("temporal", "c1_correlate", "c1.directionality_triage.csv"),
			("cohort", "c1_correlate", "c1.cohort.csv"),
			("mediation", "c4_connect", "c4.mediation_all.csv"),
			("modules", "c4_connect", "c4.supervised_module_membership.csv"),
			("deployed_concept_fidelity", "c4_connect", "c4.focus.deployed_concept_fidelity.csv"),
			("bidirectional_mr", "c2_cause", "c2.bidirectional_mr.csv"),
			("nonlinear", "c4_connect", "c4.nonlin_tests.csv"),
			("nonlinear_curves", "c4_connect", "c4.nonlin_curves.csv"),
			("same_locus", "c3_coloc", "c3.same_locus_evidence.csv"),
			("coloc_robust", "c3_coloc", "c3.coloc_summary.csv"),
			("coloc_prior", "c3_coloc", "c3.credible_set_by_locus.csv"),
			("susie_pairs", "c3_coloc", "c3.susie_signal_pairs.csv"),
			("mr_signal_evidence", "c3_coloc", "c3.same_locus_signal_evidence.csv"),
			("susie_diagnostics", "c3_coloc", "c3.susie_diagnostics.csv"),
			("cell", "c5_cellulation", "c5.cell.enrichment.csv"),
			("cell_contrasts", "c5_cellulation", "c5.cell.panel_contrasts.csv"),
			("cell_evidence", "c5_cellulation", "c5.evidence_long.csv"),
			("cell_status", "c5_cellulation", "c5.cellulation_status.csv"),
			("cigma", "c5_cellulation", "c5.CIGMA_annotation.csv"),
			("mr_scope", "c2_cause", "c2.method_scope_audit.csv"),
			("dandelion", "c2_cause", "c2.dandelion_input_audit.csv"),
			("dandelion_lolo", "c2_cause", "c2.dandelion_leave_locus_out.csv"),
			("dandelion_native", "c2_cause", "c2.dandelion_native_diagnostics.csv"),
			("state_projection", "c2_cause", "c2.state_projection.csv"),
			("decomposition", "c2_cause", "c2.individual_genetic_decomposition.csv"),
			("age_models", "c4_connect", "c4.age_models.csv"),
		]:
			d = read(module, filename, name)
			self.add(name, d)
			if name == "cell_status":
				details = (
					"; ".join(d.analysis.astype(str) + ": " + d.status.astype(str))
					if len(d)
					else "Cellulation 结果缺失"
				)
				self.claim(
					Y,
					layer,
					"Cellulation",
					"以实际完成状态为准",
					details,
					"表达富集不证明分泌来源；代谢物不直接映射至单一编码基因；CIGMA需要相应数据或完成的原生输出。",
					"cell_status",
				)
		self.scientific_summaries(Y, layer, ct, proxy)

	def current(self, name, Y, layer):
		frames = self.tables.get(name, [])
		d = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
		return d.loc[d.Y.eq(Y) & d.layer.eq(layer)].copy() if {"Y", "layer"} <= set(d) else d

	def scientific_summaries(self, Y, layer, contrasts, proxy):
		get = lambda name: self.current(name, Y, layer)
		s = twin_summary(get('abm_metrics'))
		self.add('abm_summary', pd.DataFrame([dict(Y=Y, layer=layer, **s)]))
		self.claim(Y, layer, 'Twin confidence', s['status'], s['key_result'], ABM_BOUNDARY,
			'abm_metrics;abm_paired;abm_support;abm_confidence')
		status = get('deviation_status')
		text = '; '.join(str(v) for v in status.get('validation_status', pd.Series(dtype=str)).dropna().unique())
		self.claim(Y, layer, 'Molecular deviation', text or 'unavailable',
			text or 'Current molecular deviation results unavailable.', ABM_BOUNDARY, 'deviation_metrics;deviation_status')
		if layer == "met":
			concordance = metabolite_concordance(get("met_pgs_pairs")).assign(Y=Y, layer=layer)
			self.add("met_pgs_concordance", concordance)
			fmt = lambda v: f"{v:.3f}" if pd.notna(v) else "unavailable"
			concordance_text = "; ".join(f"{r.analysis}: N={r.N_pairs}, Spearman rho={fmt(r.spearman_rho)}, same-direction fraction={fmt(r.same_direction_fraction)}" for r in concordance.itertuples())
			self.claim(Y, layer, "Metabolite inherited–measured concordance",
				"exploratory" if concordance.N_pairs.gt(0).any() else "unavailable", concordance_text, MET_BOUNDARY, "met_pgs_pairs;met_pgs_concordance")
			anchor = get("met_anchor")
			anchor = anchor[anchor.feature.eq("L_VLDL_TG.pct")] if "feature" in anchor else pd.DataFrame()
			anchor_text = "; ".join(f"{r.analysis} / {r.model}: beta={r.beta:.4g}" for r in anchor.itertuples() if np.isfinite(r.beta)) if {"analysis", "model", "beta"} <= set(anchor) else ""
			self.claim(Y, layer, "L_VLDL_TG.pct illustrative anchor", "exploratory" if anchor_text else "unavailable",
				anchor_text or "Anchor effects unavailable.", ANCHOR_BOUNDARY + " " + MET_BOUNDARY,
				"met_anchor;met_pgs_conditional;met_measured_conditional;met_trajectory")
			temporal = get("temporal")
			if "direction_class" in temporal:
				counts = temporal.groupby("direction_class", dropna=False).size().rename("N_features").reset_index()
				self.add("met_temporal_counts", counts.assign(Y=Y, layer=layer))
			if concordance.N_pairs.gt(0).any():
				s["key_result"] += " Metabolite inherited–measured concordance: " + concordance_text + "."
				if s["status"] == "unavailable": s["status"] = "exploratory"
		c4 = c4_endpoint_summary(contrasts, proxy)
		self.add("c4_endpoints", pd.DataFrame([dict(Y=Y, layer=layer, **c4)]))
		gr, genetic = get("decomposition"), get("genetic")
		gr_tests = finite_column(gr, "FDR_incident_p_difference")
		if gr_tests.notna().sum() == 0:
			gr_tests = finite_column(genetic, "GR_difference_FDR")
		genetic_claims = [c["result"] for c in self.overview if c["Y"] == Y and c["layer"] == layer and c["module"] == "C2"]
		cell = get("cell_status")
		cell_text = "; ".join(cell.analysis.astype(str) + ": " + cell.status.astype(str)) if {"analysis", "status"} <= set(cell) else "Cellulation status unavailable"
		coloc_tables = ["same_locus", "coloc_robust", "coloc_prior", "susie_pairs", "susie_diagnostics", "mr_signal_evidence"]
		coloc_text = "; ".join(f"{name}: {len(get(name))} rows" if len(get(name)) else f"{name}: unavailable" for name in coloc_tables)
		for module, status, result, boundary, sources in [
			("C1", s["status"], s["key_result"], ABM_BOUNDARY + (" " + MET_BOUNDARY + " " + ANCHOR_BOUNDARY if layer == "met" else ""),
				["abm_metrics", "abm_paired", "abm_support"] + (list(MET_TABLES) + ["temporal"] if layer == "met" else [])),
			("C2", "exploratory" if len(gr) or len(genetic) else "unavailable",
				" ".join(genetic_claims) + f" Decomposition: {len(gr)} rows; direct G−R tests: {gr_tests.notna().sum()}, FDR<0.05: {int((gr_tests < .05).sum())}. Discordant measured/PGS/MR evidence remains visible in the source tables.",
				"G and R share the training-fold measured-omic scale. R includes uncaptured genetics, environment, treatment, disease processes and assay error. Predictive biomarkers need not be causal antecedents; no causal/consequence percentage is identifiable. MR / reverse MR are separate evidence.",
				["genetic", "decomposition", "temporal", "bidirectional_mr"]),
			("C4", c4["prediction_status"], c4["key_result"],
				"Prediction: paired ΔAUC versus same-budget NS. Explainability: held-out ΔR² for LE8 proxy reconstruction. LE8 supervision can impose modifiable/interpretable biological structure while preserving broadly similar prediction; proxy fidelity does not prove intervention response.", ["contrasts", "proxy"]),
			("C3", "exploratory" if any(len(get(n)) for n in coloc_tables) else "unavailable", coloc_text,
				"Colocalization supports a shared regional signal; it does not prove mediation, a unique causal variant or a causal protein.", coloc_tables),
			("C5", "exploratory" if len(get("cell")) or len(get("cigma")) else "unavailable", cell_text + (f"; CIGMA: {len(get('cigma'))} rows" if len(get("cigma")) else "; CIGMA unavailable"),
				"Cell-expression enrichment is cellular context, not automatically secretion origin, causal cell type or a cellular aging clock.", ["cell_status", "cell", "cell_contrasts", "cigma"]),
		]:
			files = sorted({v for name in sources for v in get(name).get("source_file", pd.Series(dtype=str)).dropna().astype(str)})
			self.add("module_summary", pd.DataFrame([dict(Y=Y, layer=layer, module=module,
				priority="PRIMARY" if module in MODULE_ORDER[:3] else "SUPPORT", title=MODULE_TITLES[module], status=status,
				key_result=result, interpretation_boundary=boundary, evidence=";".join(sources), source=";".join(files))]))

	def finish(self, out, index_out):
		order = {
			name: i
			for i, name in enumerate(
				[
					"Twin confidence",
					"Molecular deviation",
					"ABM 分层",
					"Metabolite inherited–measured concordance",
					"L_VLDL_TG.pct illustrative anchor",
					"遗传与实测差异",
					"LE8 可解释性",
					"LE8 预测增益",
					"Yin–Yang 增益",
					"炎症分层异质性",
					"Cellulation",
				]
			)
		}
		overview = pd.DataFrame(self.overview)
		if len(overview):
			overview = (
				overview.assign(_order=overview.question.map(order))
				.sort_values(["Y", "layer", "_order"])
				.drop(columns="_order")
			)
		self.tables["overview"] = [overview]
		self.tables["sources"] = [pd.DataFrame(self.sources)]
		names = [
			"module_summary", "abm_summary", "c4_endpoints",
			*MET_TABLES, "met_pgs_concordance", "met_temporal_counts", "coloc_robust", "coloc_prior", "cigma",
			"prediction",
			"contrasts",
			"members",
			"pillars",
			"heterogeneity",
			"inflammation_definition",
			"design",
			"fit",
			"proxy",
			"abm_metrics",
			"abm_coverage",
			"abm_paired",
			"abm_support", "abm_training", "abm_confidence", "deviation_metrics", "deviation_status", "deviation_support",
			"genetic",
			"temporal",
			"cohort",
			"mediation",
			"modules", "deployed_concept_fidelity", "bidirectional_mr",
			"nonlinear",
			"nonlinear_curves",
			"same_locus",
			"susie_pairs", "susie_diagnostics", "mr_signal_evidence",
			"cell",
			"cell_status",
			"cell_contrasts", "cell_evidence", "concept_coefficients", "concept_status", "concept_fold_panels",
			"mr_scope",
			"dandelion", "dandelion_lolo", "dandelion_native", "state_projection", "decomposition", "age_models",
			"overview",
			"sources",
		]
		for name in names:
			frames = self.tables.get(name, [])
			d = (
				pd.concat(frames, ignore_index=True)
				if frames
				else pd.DataFrame(columns=["Y", "layer", "status"])
			)
			atomic_csv(d, out / ("final.questions." + name + ".csv"))
			atomic_csv(d, index_out / ("question_" + name + ".csv"))
		write_question_html(out, overview)
		return names


def build_questions(
	root, traits, layers, index_out, candidates, proxies, max_bytes=128 * 1024**2
):
	q = Questions(root, traits, layers, max_bytes)
	for Y in traits:
		for layer in layers:
			q.layer(Y, layer, candidates, proxies)
	out = Path(root) / "final"
	out.mkdir(parents=True, exist_ok=True)
	return q.finish(out, Path(index_out))


def question_figures(out, traits, layers):
	"""Fixed display specifications, never a best-performing budget search."""
	import matplotlib

	matplotlib.use("Agg")
	import matplotlib.pyplot as plt

	out = Path(out)
	read = lambda name: pd.read_csv(out / ("final.questions." + name + ".csv"))
	pairs = [(y, l) for y in traits for l in layers]

	def select(d, y, l):
		return d[d.Y.eq(y) & d.layer.eq(l)].copy() if len(d) else d

	def absent(ax, message):
		ax.text(
			0.5,
			0.5,
			message,
			ha="center",
			va="center",
			transform=ax.transAxes,
			wrap=True,
		)
		ax.set_axis_off()

	manifest = []

	def save(fig, name, sources, caption):
		fig.savefig(out / (name + ".png"), dpi=170, bbox_inches="tight")
		plt.close(fig)
		manifest.append(
			dict(
				figure=name,
				source_tables=";".join(sources),
				caption=caption,
				source_hashes=json.dumps(
					{s: digest(out / f"final.questions.{s}.csv") for s in sources},
					sort_keys=True,
				),
			)
		)

	plt.rcParams.update(
		{
			"font.size": 9,
			"axes.spines.top": False,
			"axes.spines.right": False,
			"pdf.fonttype": 42,
		}
	)
	ct, px = read("contrasts"), read("proxy")
	fig, axs = plt.subplots(
		2,
		len(pairs),
		figsize=(5 * len(pairs), 9),
		squeeze=False,
		constrained_layout=True,
	)
	for col, (y, l) in enumerate(pairs):
		d = select(ct, y, l)
		budget = display_budget(d)
		if len(d):
			d = d[
				d.stratum.eq("All")
				& d.landmark.eq(0)
				& d.reference.eq(f"NS_{budget}")
				& truth(d.comparison_valid)
			]
		ax = axs[0, col]
		if len(d):
			d = d.sort_values("model")
			for i, (_, r) in enumerate(d.iterrows()):
				ax.plot([r.delta_lo, r.delta_hi], [i, i], color="#247b83")
				ax.plot(r.delta_AUC, i, "o", color="#247b83")
			ax.set_yticks(range(len(d)), d.model.str.replace(f"_{budget}", "", regex=False))
			ax.axvline(0, ls="--", color="grey", lw=0.8)
			ax.set_xlabel("Paired delta AUC vs NS (nominal 95% CI)")
		else:
			absent(ax, "Equal-budget paired prediction unavailable")
		ax.set_title(f"{y} | {l} | {budget} assays, all validation participants")
		p = select(px, y, l)
		if len(p):
			p = p[p.model.eq(f"YS_YinYang_{budget}")]
		ax = axs[1, col]
		if len(p):
			p = p.set_index("component").reindex(list(PILLARS)).reset_index()
			for i, r in p.iterrows():
				if "lo" in r and np.isfinite(r.lo) and np.isfinite(r.hi):
					ax.plot([r.lo, r.hi], [i, i], color="#247b83")
				ax.plot(r.delta_R2_vs_NS, i, "o", color="#247b83")
			ax.set_yticks(
				range(8),
				[
					"Diet",
					"Activity",
					"Smoking",
					"Sleep",
					"BMI",
					"Lipids",
					"Glucose",
					"Blood pressure",
				],
			)
			ax.axvline(0, ls="--", color="grey", lw=0.8)
			ax.set_xlabel("Held-out delta R² vs NS (paired 95% CI if available)")
		else:
			absent(ax, "LE8 reconstruction unavailable")
		ax.set_title(f"YS YinYang, {budget} assays: all eight LE8 components")
	save(
		fig,
		"Fig6.question_LE8",
		["contrasts", "proxy"],
		"Display uses the configured budget (default 10), or the closest available assay budget; never chosen by model performance. All budgets retained in workbooks/Shiny. Prediction and reconstruction are distinct endpoints. Intervals condition on fitted panels; clinical-only shrinkage must be checked in the contrast table. No claim of intervention responsiveness.",
	)

	metrics, coverage = read("abm_metrics"), read("abm_coverage")
	fig, axs = plt.subplots(
		2,
		len(pairs),
		figsize=(4.5 * len(pairs), 8),
		squeeze=False,
		constrained_layout=True,
	)
	for col, (y, l) in enumerate(pairs):
		d = select(metrics, y, l)
		ax = axs[0, col]
		if has(d, 'model', 'subset', 'AUC_IPCW'):
			for model in sorted(d.model.unique()):
				z = d[d.model.eq(model)].set_index('subset').reindex(['all','selected','rejected'])
				ax.plot(range(3), z.AUC_IPCW, 'o-', label=model)
			ax.set_xticks(range(3), ['All','Selected','Rejected'])
			ax.set_ylabel('IPCW AUC')
			ax.legend(fontsize=6)
		else: absent(ax, 'Twin metrics unavailable')
		ax.set_title(f'{y} | {l}')
		ax = axs[1, col]
		contrasts = select(read('abm_paired'), y, l)
		if has(contrasts, 'contrast', 'metric', 'delta', 'lower', 'upper'):
			z = contrasts[contrasts.metric.eq('AUC_IPCW')]
			for j, (_, row) in enumerate(z.iterrows()):
				ax.plot([row.lower,row.upper], [j,j], color='#247b83')
				ax.plot(row.delta,j,'o',color='#247b83')
			ax.set_yticks(range(len(z)),z.contrast)
			ax.axvline(0,color='grey',ls='--')
			ax.set_xlabel('Paired AUC contrast; predeclared comparison masks')
		else: absent(ax, 'Paired twin contrasts unavailable')
	save(fig, 'Fig7.question_ABM', ['abm_metrics','abm_paired','abm_support'], ABM_BOUNDARY)

	genetic = read("genetic")
	fig, axs = plt.subplots(
		1,
		len(pairs),
		figsize=(4.7 * len(pairs), 4.5),
		squeeze=False,
		constrained_layout=True,
	)
	for col, (y, l) in enumerate(pairs):
		d = select(genetic, y, l)
		if len(d):
			d = d[
				truth(d.same_people)
				& d.scope.eq("existing_full")
				& d.adjustment.eq("basic")
			]
		ax = axs[0, col]
		if len(d):
			flag = (d.measured_p >= 0.05) & (d.pgs_FDR < 0.05)
			ax.scatter(
				d.measured_beta,
				d.pgs_beta,
				c=np.where(flag, "#c26427", "#247b83"),
				alpha=0.4,
				s=10,
			)
			for _, r in d[
				d.feature.isin(["GDF15", "MMP12", "ApoB", "L_VLDL_TG.pct"])
			].iterrows():
				offset = {
					"GDF15": (-5, 12),
					"MMP12": (-5, -16),
					"ApoB": (4, 7),
					"L_VLDL_TG.pct": (4, -13),
				}[r.feature]
				ax.annotate(
					r.feature,
					(r.measured_beta, r.pgs_beta),
					xytext=offset,
					textcoords="offset points",
					fontsize=8,
					ha="right" if offset[0] < 0 else "left",
					arrowprops=dict(arrowstyle="-", lw=0.5, color="grey"),
				)
			ax.axhline(0, color="grey", lw=0.6)
			ax.axvline(0, color="grey", lw=0.6)
			ax.set_xlabel("Measured beta (per own SD)")
			ax.set_ylabel("PGS beta (per own SD)")
			if d.pgs_beta.abs().max() > 0.5:
				# Retain extreme estimates without compressing the dense central
				# cloud into a line; the transformed axis is explicitly labelled.
				ax.set_yscale("symlog", linthresh=0.05)
				ax.set_ylabel("PGS beta (per own SD; symlog axis)")
		else:
			absent(ax, "Matched full-scope basic-adjustment data unavailable")
		ax.set_title(f"{y} | {l} | same people")
	save(
		fig,
		"Fig8.question_genetic",
		["genetic", "temporal", "same_locus"],
		"All matched assays in existing_full/basic; orange = measured P>=0.05 and PGS full-family FDR<0.05. Distinct own-SD scales, not an effect-difference test. Named markers are prespecified examples. PGS is inherited propensity, not birth concentration; G/R differences, MR and exact-locus coloc require separate evaluation.",
	)

	met_pairs, met_summary = read("met_pgs_pairs"), read("met_pgs_concordance")
	if len(met_pairs):
		fig, axs = plt.subplots(len(traits), 3, figsize=(15, 4.5 * len(traits)), squeeze=False, constrained_layout=True)
		for row, y in enumerate(traits):
			for col, analysis in enumerate(MET_ANALYSES):
				ax = axs[row, col]
				d = select(met_pairs, y, "met")
				d = d[d.analysis.eq(analysis)] if "analysis" in d else pd.DataFrame()
				x, measured = finite_column(d, "pgs_full_beta"), finite_column(d, "observed_beta")
				valid = x.notna() & measured.notna()
				if valid.any():
					ax.scatter(x[valid], measured[valid], s=13, alpha=.6, color="#247b83")
					ax.axhline(0, color="grey", lw=.5); ax.axvline(0, color="grey", lw=.5)
					stats = select(met_summary, y, "met")
					stats = stats[stats.analysis.eq(analysis)].iloc[0]
					ax.text(.03, .97, f"N={stats.N_pairs}; Spearman rho={stats.spearman_rho:.3f}", transform=ax.transAxes, va="top")
					ax.set(xlabel="Inherited omic PGS effect", ylabel="Measured adult omic effect")
				else: absent(ax, "Concordance unavailable")
				ax.set_title(f"{y} | {analysis}")
		save(fig, "Fig10.question_metabolite", ["met_pgs_pairs", "met_pgs_concordance"], MET_BOUNDARY)
	figure_order = {name: i for i, name in enumerate(["Fig7.question_ABM", "Fig10.question_metabolite", "Fig8.question_genetic", "Fig6.question_LE8"])}
	manifest.sort(key=lambda row: figure_order.get(row["figure"], len(figure_order)))
	pd.DataFrame(manifest).to_csv(
		out / "final.questions.figure_manifest.csv", index=False
	)
	write_question_html(out, read("overview"))


# Public command-line entry points.
def report_main(argv=None):
	ap = argparse.ArgumentParser(
		description="Assemble publication figures from existing aggregate results",
		prog="final.py report",
	)
	ap.add_argument("--Y", default="cvd_cad,ra")
	ap.add_argument("--biom", default="prot,met")
	ap.add_argument("--analysis-root", type=Path, required=True)
	ap.add_argument("--out", type=Path)
	ap.add_argument("--coefficient-epsilon", type=float, default=1e-8)
	ap.add_argument("--strict", action="store_true")
	arg = ap.parse_args(argv)
	if not math.isfinite(arg.coefficient_epsilon) or arg.coefficient_epsilon <= 0:
		ap.error("Invalid coefficient epsilon")
	for x in arg.Y.split(",") + arg.biom.split(","):
		if not x or "/" in x or ".." in x:
			ap.error("Unsafe trait/layer name")
	if not set(arg.biom.split(",")) <= {"prot", "met"}:
		ap.error("Invalid molecular layer")
	report = Report(arg)
	report.load()
	report.make_figures()
	report.finish()


def tables_main(argv=None):
	ap = argparse.ArgumentParser(
		description="Build a systematic evidence atlas for one trait",
		prog="final.py tables",
	)
	ap.add_argument("--root", type=Path, required=True)
	ap.add_argument(
		"--outdir",
		type=Path,
		help="Separate report location; source tables never modified",
	)
	ap.add_argument("--no-plots", action="store_true")
	ap.add_argument("--cell-annotation", action="store_true")
	args = ap.parse_args(argv)
	root = args.root.resolve()
	out = (
		args.outdir or final_dir(root.parent, root.name, kind="systematic")
	).resolve()
	if not out.is_relative_to(root.parent / "final"):
		raise ValueError("Final tables must be written inside <analysis-root>/final")
	out.mkdir(parents=True, exist_ok=True)
	tables = aggregate(root)
	for name, d in tables.items():
		d.to_csv(out / f"{name}.csv", index=False)
	if not args.no_plots:
		plots(tables, out, root.name)
	status = dict(
		source_root=str(root),
		interpretation="Descriptive consolidation; no independent validation or causal voting",
	)
	if args.cell_annotation:
		status["cell_annotation"] = cell_annotation(
			root, Path(__file__).resolve().parent, out
		)
	(out / "status.json").write_text(json.dumps(status, indent=2) + "\n")


def main(argv=None):
	argv = list(sys.argv[1:] if argv is None else argv)
	commands = {"report": report_main, "tables": tables_main}
	if argv and argv[0] in commands:
		return commands[argv[0]](argv[1:])
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument(
		"command", choices=commands, help="Run COMMAND --help for its options"
	)
	parser.parse_args(argv)


if __name__ == "__main__":
	try:
		main()
	except Exception as exc:
		print("ERROR: " + str(exc), file=sys.stderr)
		raise
