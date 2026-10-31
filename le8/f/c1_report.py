#!/usr/bin/env python3
"""C1 saved-model comparisons and publication figures; no risk-model training."""
from __future__ import annotations


# Conventional and twin-model comparison
import gzip, json, re, hashlib
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

if "le8_common" not in sys.modules:
	_io_spec = spec_from_file_location("le8_common", Path(__file__).with_name("0.common.py"))
	_io_module = module_from_spec(_io_spec)
	sys.modules[_io_spec.name] = _io_module
	try:
		_io_spec.loader.exec_module(_io_module)
	except BaseException:
		sys.modules.pop(_io_spec.name, None)
		raise

import numpy as np
import pandas as pd
from c1_evaluation import infer_selector
from c1_protocol import RELEASE, eid_hash, grouping_audit
# C1_LOCKPAPER::2026-10-17.c1-lockpaper-v1
from sklearn.metrics import roc_auc_score
from abm import uno_c, harrell_c, ipcw_information


def _metrics(time, event, score, km, probability=False):
	time = np.asarray(time, float)
	event = np.asarray(event, int)
	score = np.asarray(score, float)
	y, w = km.labels_weights(time, event)
	ok = (w > 0) & np.isfinite(score)
	finite = np.isfinite(score)
	row = dict(
	 N=int(finite.sum()),
	 known_N=int(ok.sum()),
	 cases=int(y[ok].sum()),
	 controls=int(np.sum(ok & (y == 0))),
	 Brier_normalization="Hajek: sum(IPCW*loss)/sum(IPCW)",
	 **ipcw_information(y[finite], w[finite]),
	)
	row["AUC_IPCW"] = (
	 float(roc_auc_score(y[ok], score[ok], sample_weight=w[ok])) if len(np.unique(y[ok])) == 2 else np.nan
	)
	row["Harrell_C"], row["Harrell_comparable_pairs"] = harrell_c(time[finite], event[finite], score[finite])
	# Early-censored people can still contribute valid comparable pairs to Uno C.
	row["Uno_C_horizon"], row["comparable_pairs"] = uno_c(time[finite], event[finite], score[finite], km)
	if probability and ok.any():
		p = np.clip(score[ok], 1e-7, 1 - 1e-7)
		row["Brier_IPCW"] = float(np.average((y[ok] - p) ** 2, weights=w[ok]))
		row["LogLoss_IPCW"] = float(
		 np.average(-(y[ok] * np.log(p) + (1 - y[ok]) * np.log(1 - p)), weights=w[ok])
		)
	else:
		row["Brier_IPCW"] = np.nan
		row["LogLoss_IPCW"] = np.nan
	return row


def benchmark_rows(benchmark, trait, layer, horizon, outer_roster_sha256=None):
	"""Validate saved C1 metadata and unique model/EID keys without fitting."""
	if benchmark.get("format") != "le8-c1-benchmark-v1":
		raise ValueError("Unknown C1 benchmark format")
	if benchmark.get("trait") != trait or benchmark.get("layer") != layer:
		raise ValueError("trait_or_layer_mismatch")
	if not np.isclose(float(benchmark["horizon"]), horizon, rtol=0, atol=1e-10):
		raise ValueError("prediction_horizon_mismatch")
	if outer_roster_sha256 is not None and benchmark.get("lockpaper_protocol", {}).get("outer_roster_sha256") != outer_roster_sha256:
		raise ValueError("C1 outer_roster_mismatch")
	c1 = pd.DataFrame(benchmark["predictions"])
	if c1.empty:
		return c1
	if not {"eid", "biom_set", "model", "time", "event", "score"} <= set(c1):
		raise ValueError("Incomplete C1 benchmark predictions")
	c1["eid"] = c1.eid.astype(str)
	if c1.duplicated(["eid", "biom_set", "model"]).any():
		raise ValueError("Duplicate C1 model/participant predictions")
	return c1


def conventional_scores(c1, ids, time, event):
	"""Map all available conventional arms to one roster; absent rows stay NaN."""
	out = {}
	for (panel, model), rows in c1.groupby(["biom_set", "model"], sort=False):
		aligned = rows.set_index("eid").reindex(np.asarray(ids, str))
		present = pd.Index(rows.eid).get_indexer(np.asarray(ids, str)) >= 0
		if not np.allclose(aligned.time.to_numpy()[present], np.asarray(time)[present], rtol=0, atol=1e-8) or not np.array_equal(aligned.event.to_numpy()[present], np.asarray(event)[present]):
			raise ValueError("C1 outcome_definition_mismatch")
		out[f"C1_{model}::{panel}"] = aligned.score.to_numpy(float)
	return out


def _parse_mask(name):
	m = re.match(r"^(.+):([0-9.]+)$", str(name))
	if not m:
		return str(name), np.nan
	return m.group(1), float(m.group(2))


def _aurc(d, metric):
	g = d[np.isfinite(d[metric]) & np.isfinite(d.coverage)].sort_values("coverage")
	if len(g) < 2:
		return np.nan
	lo, hi = g.coverage.min(), g.coverage.max()
	if hi <= lo:
		return np.nan
	return float(np.trapezoid(g[metric], g.coverage) / (hi - lo))


def _stable_hash(value, salt):
	return int.from_bytes(hashlib.sha256(f"{salt}:{value}".encode()).digest()[:8], "big")


def _risk_matched_control(selected, qc, risk, ids, salt, bins=10):
	"""Outcome-blind exact-count matching on frozen predicted-risk strata."""
	risk = np.asarray(risk, float)
	ids = np.asarray(ids, str)
	ok = np.asarray(qc, bool) & np.isfinite(risk)
	selected = np.asarray(selected, bool) & ok
	q = np.unique(np.quantile(risk[ok], np.linspace(0, 1, bins + 1))) if ok.any() else np.array([])
	if len(q) < 3:
		labels = np.zeros(len(ids), int)
	else:
		q[0], q[-1] = -np.inf, np.inf
		labels = np.searchsorted(q[1:-1], risk, side="right")
	out = np.zeros(len(ids), bool)
	for b in np.unique(labels[ok]):
		pool = np.flatnonzero(ok & (labels == b))
		count = int(np.sum(selected & (labels == b)))
		order = sorted(pool, key=lambda j: _stable_hash(ids[j], salt))
		out[np.asarray(order[:count], int)] = True
	return out


def compare(benchmark, prediction, truth, km, config):
	def unavailable(reason):
		return {
		 "status": pd.DataFrame([dict(status=reason)]),
		 "metrics": pd.DataFrame(),
		 "coverage_curve": pd.DataFrame(),
		 "selector_summary": pd.DataFrame(),
		 "same_subset_comparison": pd.DataFrame(),
		}

	if benchmark.get("format") != "le8-c1-benchmark-v1":
		raise ValueError("Unknown C1 benchmark format")
	if benchmark.get("trait") != config.trait or benchmark.get("layer") != config.layer:
		return unavailable("trait_or_layer_mismatch")
	if not np.isclose(float(benchmark["horizon"]), config.horizon, rtol=0, atol=1e-10):
		return unavailable("prediction_horizon_mismatch")
	ids = np.asarray(prediction["eid"], str)
	truth_ids = np.asarray(truth["eid"], str)
	if len(set(ids)) != len(ids) or not np.array_equal(ids, truth_ids):
		raise ValueError("Invalid ABM prediction/truth roster")
	time, event = np.asarray(truth["time"], float), np.asarray(truth["event"], int)
	c1 = benchmark_rows(benchmark, config.trait, config.layer, config.horizon)
	if c1.empty:
		return unavailable("C1_predictions_unavailable")
	joined = c1.merge(pd.DataFrame({"eid": ids, "abm_row": np.arange(len(ids))}), on="eid")
	if joined.empty:
		return unavailable("no_common_heldout_participants")
	jj = joined.abm_row.to_numpy()
	if not np.allclose(joined.time, time[jj], rtol=0, atol=1e-8) or not np.array_equal(
	 joined.event.to_numpy(), event[jj]
	):
		return unavailable("outcome_definition_mismatch")

	qc = np.asarray(prediction["diagnostics"]["technical_qc"], bool)
	masks = {
	 "all": qc,
	 "selected": np.asarray(prediction["selected"], bool) & qc,
	 "rejected": qc & ~np.asarray(prediction["selected"], bool),
	}
	masks.update({k: np.asarray(v, bool) & qc for k, v in prediction.get("masks", {}).items()})
	abm_models = [
	 n
	 for n in dict.fromkeys(
	  [
	   getattr(config, "primary", "elasticnet_clean"),
	   "twin_clean",
	   "twin_full",
	   "twin_random",
	   "twin_integrated_clean",
	   "twin_integrated_full",
	   "twin_integrated_random",
	   "elasticnet_full",
	   "elasticnet_clean",
	   "elasticnet_random",
	   "lightgbm_full",
	   "lightgbm_clean",
	   "knn_unweighted_full",
	   "knn_unweighted_clean",
	   "mixture_experts",
	   "fallback_policy",
	  ]
	 )
	 if n in prediction["probability"]
	]

	# 1) Fair pairwise comparison: exact same people for C1 and ABM.
	rows = []
	baseline_cache = {}
	abm_cache = {}
	for (panel, base_name), d in joined.groupby(["biom_set", "model"], sort=False):
		d = d.sort_values("abm_row")
		ix0 = d.abm_row.to_numpy()
		base = d.score.to_numpy(float)
		for subset, mask in masks.items():
			route, cov = _parse_mask(subset)
			for am in abm_models:
				ap = np.asarray(prediction["probability"][am], float)[ix0]
				use = mask[ix0] & np.isfinite(base) & np.isfinite(ap)
				ix = ix0[use]
				# Identical cohorts recur for every arm/panel. Cache the exact
				# metric calculation; never change its rows or probability values.
				roster_key = hashlib.sha256(ix.tobytes()).digest()
				bk = (panel, base_name, roster_key)
				ak = (am, roster_key)
				if bk not in baseline_cache:
					baseline_cache[bk] = _metrics(time[ix], event[ix], base[use], km, False)
				if ak not in abm_cache:
					abm_cache[ak] = _metrics(time[ix], event[ix], ap[use], km, True)
				bm = baseline_cache[bk]
				aa = abm_cache[ak]
				rows.append(
				 dict(
				  baseline_panel=panel,
				  baseline_model=base_name,
				  ABM_model=am,
				  model=am,
				  subset=subset,
				  selector_route=route,
				  coverage_target=cov,
				  N=aa["N"],
				  known_N=aa["known_N"],
				  cases=aa["cases"],
				  controls=aa["controls"],
				  common_holdout_N=len(d),
				  ABM_test_N=len(ids),
				  comparable_pairs=aa["comparable_pairs"],
				  status=(
				   "evaluated"
				   if np.isfinite(aa["AUC_IPCW"])
				   else "insufficient_case_control_information"
				  ),
				  coverage_of_ABM_test=float(use.sum() / len(ids)),
				  same_people=True,
				  comparison_N=len(ix),
				  events=int(event[ix].sum()),
				  families=len(np.unique(np.asarray(prediction["group"])[ix])) if "group" in prediction else np.nan,
				  comparison_eid_hash=eid_hash(ids[ix]),
				  outer_roster_sha256=benchmark.get("lockpaper_protocol", {}).get("outer_roster_sha256", ""),
				  C1_Harrell_C=bm["Harrell_C"],
				  ABM_Harrell_C=aa["Harrell_C"],
				  delta_Harrell_C=aa["Harrell_C"] - bm["Harrell_C"],
				  delta_AUC_IPCW=aa["AUC_IPCW"] - bm["AUC_IPCW"],
				  C1_AUC_IPCW=bm["AUC_IPCW"],
				  ABM_AUC_IPCW=aa["AUC_IPCW"],
				  delta_AUC=aa["AUC_IPCW"] - bm["AUC_IPCW"],
				  C1_Uno_C_horizon=bm["Uno_C_horizon"],
				  ABM_Uno_C_horizon=aa["Uno_C_horizon"],
				  delta_Uno_C=aa["Uno_C_horizon"] - bm["Uno_C_horizon"],
				  ABM_Brier_IPCW=aa["Brier_IPCW"],
				  ABM_LogLoss_IPCW=aa["LogLoss_IPCW"],
				 )
				)
	pairwise = pd.DataFrame(rows)

	# 2) Selector quality independent of which new risk algorithm is used:
	# apply each selector mask to the SAME frozen elasticnet_full probability.
	curve = []
	ref = np.asarray(prediction["probability"]["elasticnet_full"], float)
	y_all, w_all = km.labels_weights(time, event)
	total_cases = max(1, int(y_all[qc].sum()))
	for subset, mask in masks.items():
		route, cov = _parse_mask(subset)
		if not np.isfinite(cov):
			continue
		use = mask & np.isfinite(ref)
		mm = _metrics(time[use], event[use], ref[use], km, True)
		rm = _risk_matched_control(mask, qc, ref, ids, f"{route}:{cov}:riskmatch")
		ruse = rm & np.isfinite(ref)
		rr = _metrics(time[ruse], event[ruse], ref[ruse], km, True)
		curve.append(
		 dict(
		  selector_route=route,
		  coverage=float(use.sum() / len(ids)),
		  coverage_target=cov,
		  N=mm["N"],
		  known_N=mm["known_N"],
		  cases=mm["cases"],
		  case_coverage=float(y_all[use].sum() / total_cases),
		  Harrell_C=mm["Harrell_C"],
		  AUC_IPCW=mm["AUC_IPCW"],
		  Uno_C_horizon=mm["Uno_C_horizon"],
		  Brier_IPCW=mm["Brier_IPCW"],
		  LogLoss_IPCW=mm["LogLoss_IPCW"],
		  riskmatched_N=rr["N"],
		  riskmatched_cases=rr["cases"],
		  riskmatched_Harrell_C=rr["Harrell_C"],
		  riskmatched_AUC_IPCW=rr["AUC_IPCW"],
		  riskmatched_Brier_IPCW=rr["Brier_IPCW"],
		  riskmatched_LogLoss_IPCW=rr["LogLoss_IPCW"],
		  delta_AUC_vs_riskmatched=mm["AUC_IPCW"] - rr["AUC_IPCW"],
		  delta_Brier_vs_riskmatched=mm["Brier_IPCW"] - rr["Brier_IPCW"],
		  delta_LogLoss_vs_riskmatched=mm["LogLoss_IPCW"] - rr["LogLoss_IPCW"],
		  fixed_risk_model="elasticnet_full",
		 )
		)
	curve = pd.DataFrame(curve)
	summ = []
	if not curve.empty:
		for route, g in curve.groupby("selector_route", sort=False):
			summ.append(
			 dict(
			  selector_route=route,
			  AURC_Brier=_aurc(g, "Brier_IPCW"),
			  AURC_LogLoss=_aurc(g, "LogLoss_IPCW"),
			  AURC_delta_Brier_vs_riskmatched=_aurc(g, "delta_Brier_vs_riskmatched"),
			  AURC_delta_LogLoss_vs_riskmatched=_aurc(g, "delta_LogLoss_vs_riskmatched"),
			  mean_delta_AUC_vs_riskmatched=float(np.nanmean(g.delta_AUC_vs_riskmatched)),
			  mean_AUC=float(np.nanmean(g.AUC_IPCW)),
			  mean_UnoC=float(np.nanmean(g.Uno_C_horizon)),
			  min_coverage=float(g.coverage.min()),
			  max_coverage=float(g.coverage.max()),
			  interpretation="Primary: negative delta loss versus risk-matched outcome-blind control; raw AURC/discrimination are secondary",
			 )
			)
	summary = (
	 pd.DataFrame(summ).sort_values("AURC_delta_Brier_vs_riskmatched", na_position="last")
	 if summ
	 else pd.DataFrame()
	)

	status = pd.DataFrame(
	 [
	  dict(
	   status="evaluated",
	   scope="intersection of C1 outer test and ABM test; every delta uses identical people",
	   common_holdout_N=int(joined.eid.nunique()),
	   ABM_test_N=len(ids),
	   C1_test_N=int(c1.eid.nunique()),
	   C1_validation_eid_hash=eid_hash(c1.eid),
	   ABM_test_eid_hash=eid_hash(ids),
	   C1_common_eid_hash=eid_hash(c1.loc[c1.eid.isin(ids), "eid"]),
	   ABM_common_eid_hash=eid_hash(ids[pd.Index(c1.eid.unique()).get_indexer(ids) >= 0]),
	   common_holdout_eid_hash=eid_hash(joined.eid),
	   outer_roster=benchmark.get("lockpaper_protocol", {}).get("outer_roster", ""),
	   outer_roster_sha256=benchmark.get("lockpaper_protocol", {}).get("outer_roster_sha256", ""),
	   horizon=config.horizon,
	   primary_selective_metric="Delta IPCW Brier/log-loss versus predicted-risk-matched outcome-blind control across coverage; raw AURC secondary",
	   co_primary_metrics="Harrell C-index over follow-up; 10-year IPCW AUC",
	   secondary_metrics="horizon-truncated Uno C; 10-year IPCW Brier; 10-year IPCW log loss",
	   caution="selected-vs-all C-index is descriptive selection composition, not algorithm superiority",
	   primary_selector=str(
	    prediction.get(
	     "primary_predictability_selector", getattr(config, "selector", "unspecified")
	    )
	   ),
	   selection_semantics=str(prediction.get("selection_semantics", "legacy core selector")),
	   test_outcomes_used_for_selection=False,
	   grouping_unit=getattr(config, "grouping_unit", "individual_unverified"),
	   independent_family_support_verified=getattr(config, "grouping_unit", "individual_unverified") == "family",
	  )
	 ]
	)
	if status.C1_common_eid_hash.iloc[0] != status.ABM_common_eid_hash.iloc[0]:
		raise ValueError("C1 and ABM common held-out eid hashes differ")
 # Keep the existing cumulative-renderer contract compact and discrimination-only.
 # Detailed route/loss results remain in their dedicated workbook sheets.
	legacy = pairwise[pairwise["subset"].isin(["all", "selected", "rejected"])].drop(
	 columns=["ABM_Brier_IPCW", "ABM_LogLoss_IPCW"]
	)
	probabilities = {"CPU_" + name: np.asarray(prediction["probability"][name], float) for name in abm_models}
	probabilities.update(conventional_scores(c1, ids, time, event))
	finite = np.logical_and.reduce([np.isfinite(v) for v in probabilities.values()])
	data = dict(eid=ids, time=time, event=event, y=y_all, w=w_all,
	 groups=np.asarray(prediction.get("group", ids)), family_available="group" in prediction,
	 grouping_unit=getattr(config, "grouping_unit", "individual_unverified"),
	 requested_masks={k: v.copy() for k, v in masks.items() if k in ("all", "selected", "rejected")},
	 probability=probabilities, masks={k: v & finite for k, v in masks.items() if k in ("all", "selected", "rejected")},
	 trait=config.trait, layer=config.layer,
	 outer_roster_sha256=benchmark.get("lockpaper_protocol", {}).get("outer_roster_sha256", ""))
	perf, delta = performance(data, km, getattr(config, "bootstrap", 0) if "group" in prediction else 0, getattr(config, "seed", 2026))
	primary = primary_tables(perf, delta)
	return {
	 **primary,
	 "status": status,
	 "metrics": legacy,
	 "same_subset_comparison": pairwise,
	 "coverage_curve": curve,
	 "selector_summary": summary,
	}


def comparison_from_paths(benchmark, twin, *, selector_inference=False):
	from abm import VERSION, model_hash, load_model

	if not benchmark.is_file():
		return {"status": pd.DataFrame([dict(status="C1_benchmark_unavailable")])}
	required = [
	 "fit/model_LOCAL_ONLY.joblib",
	 "test_predictions_LOCAL_ONLY.joblib",
	 "test_truth_LOCAL_ONLY.npz",
	]
	if not all((twin / x).is_file() for x in required):
		return {"status": pd.DataFrame([dict(status="ABM_frozen_predictions_unavailable")])}
	import joblib

	with gzip.open(benchmark, "rt") as f:
		base = json.load(f)
	model = load_model(twin / required[0])
	if getattr(model, "molecular_qc_version", None) != VERSION or base.get("lockpaper_protocol", {}).get("release") != RELEASE:
		return {"status": pd.DataFrame([dict(status="C1_LOCKPAPER_refit_required", detail="Refit both C1 modules with --replace TRUE")])}
	contract = json.loads((twin / "host_contract.json").read_text())
	if contract.get("outer_roster_sha256") != base["lockpaper_protocol"].get("outer_roster_sha256"):
		return {"status": pd.DataFrame([dict(status="outer_roster_mismatch", detail="C1 and ABM used different outer rosters")])}
	pred = abm.load_model(twin / required[1])
	if pred["frozen_model_hash"] != model_hash(model):
		raise ValueError("ABM comparison predictions do not match frozen model")
	with np.load(twin / required[2], allow_pickle=False) as truth:
		if not np.array_equal(np.asarray(pred["eid"], str), np.asarray(truth["eid"], str)):
			raise ValueError("Invalid ABM prediction/truth roster")
		ans=compare(base,pred,truth,model.km,model.config)
		audit = grouping_audit(contract)
		if audit:
			if audit.get("grouping_unit") != getattr(model.config, "grouping_unit", "individual_unverified"):
				raise ValueError("Frozen model and grouping audit disagree")
			ans["grouping_audit"] = pd.DataFrame([audit])
		if selector_inference:
			curve,inference=infer_selector(pred,truth["time"],truth["event"],pred["group"],model.km,
			 coverages=getattr(model.config,"coverages",(0.2,0.4,0.5,0.6,0.8,1.0)),
			 B=getattr(model.config,"bootstrap",1000),seed=getattr(model.config,"seed",2026),
			 grouping_unit=getattr(model.config,"grouping_unit","individual_unverified"))
			ans["lockpaper_selector_curve"]=curve
			ans["lockpaper_selector_inference"]=inference
			for key, value in audit.items():
				inference[key] = value
		return ans


def run_comparison(args, scope, work, publish, completed, finish):
	from abm import json_write
	from le8_common import write_tables, sha256

	benchmark = (
	 args.analysis_root / args.Y / args.biom / "c1_correlate/prediction/_LOCAL_ONLY/benchmark.json.gz"
	)
	twin = scope / "abm_twin"
	target = scope / "baseline_comparison"
	stage = work / "baseline_comparison"
	inputs = dict(
	 C1=benchmark,
	 model=twin / "fit/MODEL_FROZEN.json",
	 predictions=twin / "test_predictions_LOCAL_ONLY.joblib",
	 truth=twin / "test_truth_LOCAL_ONLY.npz",
	 host_contract=twin / "host_contract.json",
	)
	signature = {k: sha256(p) if p.is_file() else "unavailable" for k, p in inputs.items()}
	signature["lockpaper_release"] = RELEASE
	if completed(target) and (target / "sources.json").is_file():
		if {k:v for k,v in json.loads((target / "sources.json").read_text()).items() if k != "comparison_code"} == signature:
			return
	tables = comparison_from_paths(benchmark, twin, selector_inference=True)
	write_tables(stage, tables, "ABM twin", args.biom, args.Y, workbook="C1_ABM_comparison.xlsx")
	json_write(stage / "sources.json", signature)
	finish(stage, dict(source="same-roster C1 versus ABM selective prediction", source_hashes=signature))
	publish(stage, target)
	print("[LE8 twin] DONE fair C1/ABM comparison", flush=True)


# Saved CPU and GPU comparison
from pathlib import Path
import argparse
import gzip
import hashlib
import html
import json
import math
import os
import shutil
import tempfile
import uuid
import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import abm
from c1_protocol import eid_hash
from le8_common import write_tables, sha256

VERSION = 'c1-final-saved-comparison-v3'


def read_json(path):
    return json.loads(Path(path).read_text())


def checked_completion(root):
    path = root / 'COMPLETE.json'
    if not path.is_file(): raise FileNotFoundError(f'Completed results unavailable: {root}')
    marker = read_json(path)
    if marker.get('status') != 'completed' or not marker.get('files'):
        raise ValueError(f'Invalid completion marker: {path}')
    for name, digest in marker['files'].items():
        p = root / name
        if not p.resolve().is_relative_to(root.resolve()): raise ValueError('Invalid receipt path')
        if not p.is_file() or sha256(p) != digest: raise ValueError(f'Completed file missing or changed: {p}')
    return [path, *(root / name for name in marker['files'])]


def ids(values, label):
    x = np.asarray(values, str)
    if len(np.unique(x)) != len(x) or np.any(x == ''): raise ValueError(f'Duplicate/empty IDs: {label}')
    return x


def boolean(values):
    x = pd.Series(values)
    if not x.astype(str).str.lower().isin(['true','false','1','0']).all(): raise ValueError('Invalid saved selection mask')
    return x.astype(str).str.lower().isin(['true','1']).to_numpy()


def paired_cohort(cpu, truth, gpu, km, conventional=None):
    ci = ids(cpu['eid'], 'CPU'); ti = ids(truth['eid'], 'truth'); gi = ids(gpu.eid, 'GPU')
    if not np.array_equal(ci, ti): raise ValueError('CPU predictions/truth ID order mismatch')
    index = pd.Index(ci).get_indexer(gi); keep = index >= 0
    ix = index[keep]; g = gpu.loc[keep].reset_index(drop=True)
    if not len(ix): raise ValueError('No common test participants')
    if not np.array_equal(np.asarray(cpu['group'], str)[ix], g.family.astype(str)):
        raise ValueError('CPU/GPU family assignments differ')
    t, e = np.asarray(truth['time'])[ix], np.asarray(truth['event'])[ix]
    if not np.allclose(t, g.time, rtol=0, atol=1e-8, equal_nan=False) or not np.array_equal(e,g.event):
        raise ValueError('CPU/GPU endpoint or follow-up differs')
    y,w = km.labels_weights(t,e)
    if 'y' in g and not np.array_equal(y,g.y): raise ValueError('CPU/GPU horizon labels differ')
    probabilities = {f'CPU_{name}':np.asarray(p,float)[ix] for name,p in cpu['probability'].items()}
    probabilities['GPU_Transformer'] = g.probability.to_numpy(float)
    if conventional is not None and not conventional.empty:
        probabilities.update(conventional_scores(conventional, g.eid.to_numpy(), t, e))
    finite = np.logical_and.reduce([np.isfinite(v) for v in probabilities.values()])
    if not finite.any(): raise ValueError('No shared finite predictions across models')
    cpu_qc = np.asarray(cpu.get('diagnostics', {}).get('technical_qc'))
    if cpu_qc.shape != ci.shape or not np.isin(cpu_qc, [True, False, 0, 1]).all():
        raise ValueError('Frozen CPU technical-QC mask is required for the shared comparison')
    qc = cpu_qc.astype(bool)[ix]
    # Legacy GPU exports lack a QC column. Their shared estimand uses the
    # documented CPU reference QC; do not infer QC from zero IPCW weights,
    # since an early-censored participant also has zero horizon weight.
    gpu_qc_status = 'not_saved_using_frozen_CPU_reference_QC'
    if 'technical_qc' in g:
        qc &= boolean(g.technical_qc)
        gpu_qc_status = 'saved_GPU_QC_intersected_with_CPU_QC'
    evaluated = finite & qc
    if not evaluated.any(): raise ValueError('No shared QC-passing finite predictions')
    for name,v in probabilities.items():
        if not name.startswith('C1_') and np.any((v[finite]<0)|(v[finite]>1)): raise ValueError('Invalid risk probability: '+name)
    requested_masks = dict(all_common=np.ones(len(ix), bool),
                           CPU_selected=np.asarray(cpu['selected'],bool)[ix],
                           GPU_selected=boolean(g.selected))
    requested_masks['both_selected'] = requested_masks['CPU_selected'] & requested_masks['GPU_selected']
    masks = {name: mask & evaluated for name, mask in requested_masks.items()}
    audit = dict(CPU_test_N=len(ci),GPU_test_N=len(gi),common_test_N=len(ix),
        CPU_only_N=len(ci)-len(ix),GPU_only_N=len(gi)-len(ix),excluded_nonfinite_N=int((~finite).sum()),
        excluded_QC_N=int((~qc).sum()),excluded_QC_or_nonfinite_N=int((~evaluated).sum()),
        comparison_N=int(evaluated.sum()),family_N=len(np.unique(g.family.to_numpy()[evaluated])),
        common_test_roster_sha256=eid_hash(g.eid),comparison_roster_sha256=eid_hash(g.eid[evaluated]),
        QC_rule='frozen CPU reference QC intersected with saved GPU QC when available',GPU_QC_status=gpu_qc_status,
        evaluation_weights='same frozen CPU training-derived censoring KM for every model',
        selection_scope='Every row compares models on exactly the same saved subset; no selection refitting',
        family_match=True,endpoint_match=True)
    if 'weight' in g:
        audit['saved_GPU_vs_shared_weight_max_difference'] = float(np.max(np.abs(g.weight.to_numpy()-w)))
    return dict(eid=g.eid.to_numpy(),time=t,event=e,y=y,w=w,groups=g.family.to_numpy(),
                probability=probabilities,masks=masks,requested_masks=requested_masks,qc_mask=qc,audit=audit)


def weighted_metrics(y, w, p, probability=True):
    use = (w > 0) & np.isfinite(p)
    result = dict(AUC_IPCW=np.nan, Brier_IPCW=np.nan, LogLoss_IPCW=np.nan)
    if not use.any(): return result
    y, w, p = y[use], w[use], p[use]
    # Conventional Cox scores are unrestricted risks, not probabilities.
    if len(np.unique(y)) == 2:
        result['AUC_IPCW'] = float(roc_auc_score(y, p, sample_weight=w))
    if probability:
        p = np.clip(p, 1e-7, 1-1e-7)
        result.update(Brier_IPCW=float(np.average((y-p)**2, weights=w)),
            LogLoss_IPCW=float(np.average(-(y*np.log(p)+(1-y)*np.log(1-p)), weights=w)))
    return result


def performance(data, km, bootstrap, seed):
    if int(bootstrap) != bootstrap or bootstrap < 0:
        raise ValueError('Bootstrap replicates must be a nonnegative integer')
    bootstrap = int(bootstrap)
    minimum_valid = max(20, math.ceil(.8 * bootstrap))
    probabilities = {name: np.asarray(p, float) for name, p in data['probability'].items()}
    if not probabilities or any(p.shape != (len(data['eid']),) for p in probabilities.values()):
        raise ValueError('Every model must provide one score per comparison participant')
    common_finite = np.logical_and.reduce([np.isfinite(p) for p in probabilities.values()])
    grouping_unit = data.get('grouping_unit', 'individual_unverified')
    if grouping_unit not in ('family', 'individual_unverified'):
        raise ValueError('Invalid comparison grouping unit')
    rows, deltas = [], []
    for subset, mask in data['masks'].items():
        # Enforce one mask even for callers supplying unfiltered subset masks.
        requested = np.asarray(data.get('requested_masks', {}).get(subset, mask), bool)
        mask = np.asarray(mask, bool)
        if requested.shape != common_finite.shape or mask.shape != common_finite.shape or np.any(mask & ~requested):
            raise ValueError('Invalid requested/evaluated comparison mask: ' + subset)
        mask = mask & common_finite
        y, w = data['y'][mask], data['w'][mask]
        t, e = data['time'][mask], data['event'][mask]
        p = {k:v[mask] for k,v in probabilities.items()}
        _, codes = np.unique(data['groups'][mask], return_inverse=True)
        families = len(np.unique(codes))
        family_available = data.get('family_available', True)
        common = dict(trait=data.get('trait', ''), layer=data.get('layer', ''), subset=subset,
            comparison_N=int(mask.sum()), requested_N=int(requested.sum()), evaluated_N=int(mask.sum()),
            excluded_N=int(requested.sum() - mask.sum()),
            excluded_nonfinite_N=int(np.sum(requested & ~common_finite)),
            excluded_QC_N=int(np.sum(requested & ~np.asarray(data.get('qc_mask', np.ones(len(requested), bool)), bool))),
            events=int(e.sum()), bootstrap_cluster_N=families if family_available else np.nan,
            families=families if family_available and grouping_unit == 'family' else np.nan,
            requested_eid_hash=eid_hash(np.asarray(data['eid'])[requested]),
            comparison_eid_hash=eid_hash(np.asarray(data['eid'])[mask]),
            outer_roster_sha256=data.get('outer_roster_sha256', ''), same_people=True,
            grouping_unit=grouping_unit, Harrell_time_scope='full_follow_up_strict_time',
            IPCW_horizon=km.horizon)
        measures = {}
        for name, risk in p.items():
            measure = _metrics(t, e, risk, km, not name.startswith('C1_'))
            measures[name] = measure
            rows.append(dict(**common, model=name, model_family=name.split('_')[0], **measure,
                status='evaluated' if np.isfinite(measure['AUC_IPCW']) and np.isfinite(measure['Harrell_C']) else 'insufficient_information'))
        pairs = [(a,b) for a in p for b in p
                 if (a.startswith('CPU_') and b.startswith('C1_')) or
                    (a == 'GPU_Transformer' and b != a)]
        # Planned CPU contrasts share the exact same subset as every external
        # benchmark; cleaning effects need both co-primary metric intervals too.
        cpu_pairs = [('elasticnet_clean', 'elasticnet_full'), ('elasticnet_clean', 'elasticnet_random'),
                     ('twin_clean', 'twin_full'), ('twin_clean', 'twin_random'),
                     ('twin_integrated_clean', 'twin_integrated_full'), ('twin_integrated_clean', 'twin_integrated_random'),
                     ('twin_clean', 'twin_integrated_clean'), ('fallback_policy', 'elasticnet_full'),
                     ('lightgbm_clean', 'lightgbm_full'), ('knn_unweighted_clean', 'knn_unweighted_full'),
                     ('mixture_experts', 'elasticnet_full')]
        pairs.extend(('CPU_' + a, 'CPU_' + b) for a, b in cpu_pairs if 'CPU_' + a in p and 'CPU_' + b in p)
        metric_names = ('Harrell_C','AUC_IPCW','Brier_IPCW','LogLoss_IPCW')
        samples = {pair:{m:[] for m in metric_names} for pair in pairs}
        # Fixed scores/time never change during bootstrap: cache rank sweeps,
        # but apply each fresh family draw identically to every model.
        harrell_plans = {name:abm._harrell_bootstrap_plan(t, e, risk) for name,risk in p.items()} if bootstrap and families and family_available else {}
        rng = np.random.default_rng(seed)
        for _ in range(bootstrap if families and family_available else 0):
            multiplicity = np.bincount(rng.integers(families, size=families), minlength=families)[codes]
            estimates = {}
            for name, risk in p.items():
                estimates[name] = weighted_metrics(y, w*multiplicity, risk, not name.startswith('C1_'))
                estimates[name]['Harrell_C'] = harrell_plans[name](multiplicity)[0]
            for (a,b), metrics in samples.items():
                for metric in metrics:
                    difference = estimates[a][metric] - estimates[b][metric]
                    if np.isfinite(difference): metrics[metric].append(difference)
        for (a,b), metrics in samples.items():
            for metric, values in metrics.items():
                enough = len(values) >= minimum_valid
                ci = np.quantile(values,[.025,.975]) if enough else [np.nan,np.nan]
                difference = measures[a][metric] - measures[b][metric]
                inference_status = ('insufficient_information' if not np.isfinite(difference) else
                    'grouping_unavailable' if not family_available else
                    'bootstrap_not_requested' if not bootstrap else
                    'evaluated' if enough else 'insufficient_bootstrap')
                deltas.append(dict(**common, model_A=a, model_B=b, CPU_model=b,
                    metric=metric, N=int(mask.sum()), delta=difference,
                    GPU_minus_CPU=difference if a == 'GPU_Transformer' and b.startswith('CPU_') else np.nan,
                    lower95=ci[0], upper95=ci[1], bootstrap_requested=bootstrap, bootstrap_valid=len(values),
                    bootstrap_minimum_valid=minimum_valid,
                    bootstrap_unit=grouping_unit if family_available else 'unavailable',
                    better_GPU_direction='positive' if metric in ('AUC_IPCW','Harrell_C') else 'negative',
                    status=inference_status, inference_status=inference_status,
                    point_status='evaluated' if np.isfinite(difference) else 'insufficient_information',
                    inference='exploratory paired fixed-prediction group bootstrap; no refitting'))
    return pd.DataFrame(rows), pd.DataFrame(deltas)


def primary_tables(perf, deltas):
    performance_table = perf.rename(columns={'AUC_IPCW':'AUC_IPCW_10y',
        'Brier_IPCW':'Brier_IPCW_10y','LogLoss_IPCW':'LogLoss_IPCW_10y'})
    rows = []
    if not deltas.empty:
        for _, d in deltas.groupby(['subset','model_A','model_B'], sort=False):
            first = d.iloc[0]
            row = {k:first[k] for k in ('trait','layer','subset','model_A','model_B',
                'comparison_N','requested_N','evaluated_N','excluded_N','excluded_nonfinite_N',
                'events','families','bootstrap_cluster_N','requested_eid_hash','comparison_eid_hash',
                'outer_roster_sha256','same_people','bootstrap_unit','Harrell_time_scope','IPCW_horizon')}
            row['status'] = 'evaluated'
            for metric, label, ci_label in [('Harrell_C','Harrell_C','Harrell_C'),('AUC_IPCW','AUC_IPCW_10y','AUC')]:
                m = d.loc[d.metric == metric].iloc[0]
                row['delta_'+label] = m['delta']
                row[ci_label+'_lower95'], row[ci_label+'_upper95'] = m.lower95, m.upper95
                row[ci_label+'_bootstrap_valid'] = m.bootstrap_valid
                row[ci_label+'_inference_status'] = m.inference_status
                row[ci_label+'_point_status'] = m.point_status
                if m.status != 'evaluated': row['status'] = m.status
            rows.append(row)
    return dict(same_subset_primary_performance=performance_table,
                same_subset_primary_deltas=pd.DataFrame(rows))


def numerical_parity(cpu, saved, network):
    ci=ids(cpu['eid'],'CPU'); gi=ids(saved['eid'],'parity')
    ix=pd.Index(gi).get_indexer(ci)
    if len(ci)!=len(gi) or np.any(ix<0): raise ValueError('Parity test roster differs from CPU')
    if not np.array_equal(np.asarray(saved['group'],str)[ix],np.asarray(cpu['group'],str)):
        raise ValueError('Parity family mismatch')
    rows=[]
    exact={'eligible_match_count','matched_count','family_count','supported','technical_qc'}
    for kind in ('probability','diagnostics','masks'):
        actual=saved.get(kind,{})
        if not actual: raise ValueError('Missing saved parity fields: '+kind)
        if kind=='probability' and not {'twin_clean','twin_full'} <= set(actual):
            raise ValueError('Saved parity lacks core twin model predictions')
        for field,a in actual.items():
            if field not in cpu[kind]: raise ValueError('Unknown parity field: '+field)
            a=np.asarray(a)[ix];b=np.asarray(cpu[kind][field]);delta=np.abs(a.astype(float)-b.astype(float))
            tolerance=0. if kind=='masks' or field in exact else float(network['parity_probability_tolerance' if kind=='probability' else 'parity_continuous_tolerance'])
            finite=np.isfinite(delta).all();maximum=float(delta.max()) if finite else float('inf')
            rows.append(dict(kind=kind,field=field,N=len(ci),max_abs_difference=maximum,
                tolerance=tolerance,mismatches=int(np.count_nonzero(a!=b)),passed=bool(finite and maximum<=tolerance)))
    a=np.asarray(saved['selected'])[ix];b=np.asarray(cpu['selected'])
    rows.append(dict(kind='selector',field='selected',N=len(ci),max_abs_difference=float(np.any(a!=b)),
                     tolerance=0.,mismatches=int(np.count_nonzero(a!=b)),passed=bool(np.array_equal(a,b))))
    return pd.DataFrame(rows)


def timing_tables(gpu_root, transformer):
    rows=[]
    parity=gpu_root/'parity/timing.tsv'
    if parity.is_file():
        for row in pd.read_csv(parity,sep='\t').to_dict('records'):
            rows.append(dict(category='same_algorithm_matching',status='recorded',**row,
                interpretation='Recorded warmed-up median on identical queries; matching speedup only, not training speedup'))
    steps=transformer/'performance/step_timings.tsv'
    if steps.is_file():
        for row in pd.read_csv(steps,sep='\t').to_dict('records'):
            rows.append(dict(category='different_model_step',status='recorded',**row,
                interpretation='Recorded distinct-model core inference/training; not algorithm-equivalent speedup'))
    history=transformer/'training_history.tsv'
    if history.is_file():
        d=pd.read_csv(history,sep='\t')
        if 'seconds' in d:
            rows.append(dict(category='GPU_recorded_epochs',status='recorded',seconds=float(d.seconds.sum()),
                interpretation='Sum of recorded epoch times; excludes preparation, reporting and unrecorded work'))
    report=gpu_root/'GPU_REPORT.json'
    if report.is_file():
        for key,value in read_json(report).get('timing',{}).items():
            if key.endswith('_seconds') and value is not None:
                rows.append(dict(category='GPU_workflow',step=key,seconds=value,status='recorded',
                    interpretation='Recorded workflow wall time; not comparable to CPU training alone'))
    rows.append(dict(category='CPU_training',status='not_recorded',seconds=np.nan,
                     interpretation='CPU artifacts have no reliable total training timer; no training speedup inferred'))
    if not parity.is_file(): rows.append(dict(category='same_algorithm_matching',status='unavailable',interpretation='No saved CPU/GPU matching benchmark'))
    if not steps.is_file(): rows.append(dict(category='different_model_step',status='unavailable',interpretation='No saved paired inference timings'))
    return pd.DataFrame(rows)


def load_inputs(scope):
    cpu_root=scope/'c1_abm/abm_twin';gpu_root=scope/'c1_abm_GPU'
    files=checked_completion(cpu_root)+checked_completion(gpu_root)
    transformer=gpu_root/'transformer' if (gpu_root/'transformer').is_dir() else gpu_root
    paths=[cpu_root/'fit/model_LOCAL_ONLY.joblib',cpu_root/'test_predictions_LOCAL_ONLY.joblib',
           cpu_root/'test_truth_LOCAL_ONLY.npz',cpu_root/'host_contract.json',gpu_root/'host_contract.json',
           transformer/'test_predictions_LOCAL_ONLY.tsv.gz',transformer/'MODEL_FROZEN.json']
    for p in paths:
        if not p.is_file(): raise FileNotFoundError('Saved comparison input unavailable: '+str(p))
    files+=paths
    ca,ga=read_json(cpu_root/'host_contract.json'),read_json(gpu_root/'host_contract.json')
    for key in ('outer_roster_sha256','group_map_sha256'):
        if not ca.get(key) or ca.get(key)!=ga.get(key): raise ValueError('CPU/GPU protocol mismatch: '+key)
    if ca.get('grouping_unit')!='family' or ga.get('grouping_unit')!='family': raise ValueError('Verified family grouping required')
    # New GPU runs retain exact CPU artifact hashes. Historical independent
    # Transformer runs are instead verified by endpoint/roster/family alignment.
    for name,digest in ga.get('cpu_source',{}).get('files',{}).items():
        path=cpu_root/name
        if not path.is_file() or sha256(path)!=digest: raise ValueError('GPU used different CPU artifacts: '+name)
        files.append(path)
    model=abm.load_model(paths[0]);cpu=abm.load_model(paths[1])
    if cpu['frozen_model_hash']!=abm.model_hash(model): raise ValueError('CPU frozen model/prediction hash mismatch')
    with np.load(paths[2],allow_pickle=False) as t: truth={k:t[k] for k in t.files}
    gpu=pd.read_csv(paths[-2],sep='\t',dtype={'eid':str,'family':str})
    marker=read_json(paths[-1])
    if not np.isclose(marker['horizon'],model.config.horizon,rtol=0,atol=1e-10): raise ValueError('CPU/GPU prediction horizons differ')
    optional=[gpu_root/'parity/parity_predictions_LOCAL_ONLY.joblib',gpu_root/'parity/timing.tsv',
        gpu_root/'parity/timing_summary.json',gpu_root/'GPU_REPORT.json',
        transformer/'performance/step_timings.tsv',transformer/'performance/timing_scope.json',transformer/'training_history.tsv']
    files += [p for p in optional if p.is_file()]
    benchmark_path = scope/'c1_correlate/prediction/_LOCAL_ONLY/benchmark.json.gz'
    conventional = None
    if benchmark_path.is_file():
        with gzip.open(benchmark_path, 'rt') as stream:
            benchmark = json.load(stream)
        conventional = benchmark_rows(benchmark, scope.parent.name, scope.name,
                                     model.config.horizon, ca['outer_roster_sha256'])
        files.append(benchmark_path)
    ga = dict(ga, conventional_status='available' if conventional is not None and not conventional.empty else 'unavailable')
    fingerprint={str(p):sha256(p) for p in sorted(set(files))}
    return model,cpu,truth,gpu,ga,transformer,fingerprint,conventional


def publish(target, tables, metadata):
    target.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='_comparison_',dir=target.parent))
    try:
        write_tables(stage,tables,'CPU/GPU saved comparison',target.parent.name,metadata['trait'],workbook='CPU_GPU_comparison.xlsx')
        (stage/'comparison.json').write_text(json.dumps(metadata,indent=2,allow_nan=False))
        note=('All performance rows use the same participants within each subset and the frozen CPU censoring model. '
              'GPU Transformer is a different model; numerical equivalence applies only to the saved frozen-model parity branch. '
              'Timing entries retain their original scope. Missing timing is not estimated. No model was trained or projected.')
        body=''.join('<h2>'+html.escape(name)+'</h2>'+d.to_html(index=False,na_rep='Unavailable',escape=True) for name,d in tables.items())
        (stage/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>CPU/GPU comparison</title>'
            '<style>body{font:15px system-ui;margin:32px}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border:1px solid #ddd}h2{margin-top:32px}</style>'
            '<h1>CPU ABM / GPU comparison</h1><p>'+note+'</p>'+body)
        for path,digest in metadata.get('sources',{}).items():
            if sha256(path)!=digest: raise ValueError('Source changed during comparison: '+path)
        manifest={p.name:sha256(p) for p in stage.iterdir() if p.is_file()}
        (stage/'COMPLETE.json').write_text(json.dumps(dict(status='completed',files=manifest),indent=2))
        backup=None
        if target.exists():
            backup=target.parent/'_history'/(target.name+'_'+uuid.uuid4().hex[:10]);backup.parent.mkdir(exist_ok=True);target.rename(backup)
        try: stage.rename(target)
        except BaseException:
            if backup: backup.rename(target)
            raise
    finally:
        if stage.exists(): shutil.rmtree(stage)


def compare_scope(root, trait, layer, bootstrap=200, seed=2026, check=False):
    scope=Path(root)/trait/layer;target=scope/'c1_abm_comparison'
    model,cpu,truth,gpu,contract,transformer,sources,conventional=load_inputs(scope)
    data=paired_cohort(cpu,truth,gpu,model.km,conventional)
    data.update(trait=trait, layer=layer, outer_roster_sha256=contract['outer_roster_sha256'],
                grouping_unit=getattr(model.config, 'grouping_unit', 'individual_unverified'))
    metadata=dict(version=VERSION,trait=trait,layer=layer,bootstrap=bootstrap,seed=seed,
        sources=sources,
        conventional_status=contract['conventional_status'],
        co_primary_metrics='Harrell C-index over follow-up; 10-year IPCW AUC',
        secondary_metrics='horizon-truncated Uno C; 10-year IPCW Brier; 10-year IPCW log loss',
        caution='selected-vs-all performance is descriptive selection composition, not algorithm superiority',
        training_performed=False,inference_performed=False)
    if check: return dict(status='inputs_verified',**data['audit'])
    if (target/'comparison.json').is_file() and read_json(target/'comparison.json')==metadata:
        try:
            checked_completion(target)
            print(f'[LE8 comparison] SKIP {trait}/{layer}: verified unchanged report',flush=True)
            return target
        except (ValueError,FileNotFoundError): pass
    metrics,deltas=performance(data,model.km,bootstrap,seed)
    parity_path=scope/'c1_abm_GPU/parity/parity_predictions_LOCAL_ONLY.joblib'
    if parity_path.is_file():
        if not contract.get('cpu_source'): raise ValueError('Parity lacks frozen CPU source receipt')
        parity=numerical_parity(cpu,abm.load_model(parity_path),contract['gpu_config'])
        parity_status='passed' if parity.passed.all() else 'failed'
    else:
        parity=pd.DataFrame([dict(status='unavailable',reason='No saved frozen-model CPU/GPU parity predictions; Transformer predictions are not an equivalence test')])
        parity_status='unavailable'
    tables=dict(**primary_tables(metrics,deltas),status=pd.DataFrame([dict(status='evaluated',trait=trait,layer=layer,parity_status=parity_status,
        exploratory=True,conventional_status=contract['conventional_status'],
        co_primary_metrics=metadata['co_primary_metrics'],secondary_metrics=metadata['secondary_metrics'],
        caution=metadata['caution'],training_performed=False,inference_performed=False)]),cohort=pd.DataFrame([data['audit']]),
        performance=metrics,paired_differences=deltas,numerical_parity=parity,
        timings=timing_tables(scope/'c1_abm_GPU',transformer))
    publish(target,tables,metadata)
    print(f'[LE8 comparison] DONE {trait}/{layer}: {target}/CPU_GPU_comparison.xlsx',flush=True)
    return target


def comparison_main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis-root',type=Path,default=Path('/mnt/d/analysis/le8'))
    p.add_argument('--Y',required=True);p.add_argument('--biom',default='prot,met')
    p.add_argument('--bootstrap',type=int,default=200);p.add_argument('--seed',type=int,default=2026)
    p.add_argument('--check',action='store_true');a=p.parse_args()
    if a.bootstrap<0: p.error('--bootstrap must be nonnegative')
    pending=[]
    from threadpoolctl import threadpool_limits
    with threadpool_limits(2):
        for trait in a.Y.split(','):
            for layer in a.biom.split(','):
                try: compare_scope(a.analysis_root,trait,layer,a.bootstrap,a.seed,a.check)
                except FileNotFoundError as exc:
                    pending.append(str(exc));print('[LE8 comparison] PENDING '+str(exc),flush=True)
    if pending: raise SystemExit(2)


# Publication figures and tables
# Shared PNG byte limit; PDF/SVG export settings are preserved.
import sys as _png_sys
from pathlib import Path as _PngPath
_png_sys.path.append(str(_PngPath(__file__).resolve().parents[2] / "0f"))
from png_output import save_png_figure

import argparse, gzip, hashlib, io, json, math, zipfile
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from c1_protocol import grouping_audit, GROUPING_SOURCE

RELEASE="2026-10-17.c1-lockpaper-v1"
REPORT_VERSION="2026-10-21.same-people-reporting-v1"
LAYERS=("prot", "met")
LAYER_COLORS={"prot":"#1f77b4", "met":"#ff7f0e"}
CLEANING_CONTRASTS=("twin_library_cleaning", "twin_integrated_library_cleaning",
                    "training_effect_same_subset", "pruning_specificity")

def read_table(path):
    path=Path(path)
    data = path.read_bytes() if path.is_file() else None
    if data is None:
        for book in sorted(path.parent.glob("*.xlsx")):
            with zipfile.ZipFile(book) as archive:
                if "le8/manifest.json" not in archive.namelist(): continue
                entry=json.loads(archive.read("le8/manifest.json")).get("files",{}).get(path.name)
                if not entry: continue
                data=archive.read(entry["part"])
                if hashlib.sha256(data).hexdigest()!=entry["sha256"]:
                    raise ValueError("Corrupt frozen table: " + str(path))
                break
    if data is None: return pd.DataFrame()
    if path.suffix==".gz": data=gzip.decompress(data)
    try:
        return pd.read_csv(io.BytesIO(data),sep="\t" if ".tsv" in path.name else ",")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

def savefig(fig,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    save_png_figure(fig, path.with_suffix(".pdf"),bbox_inches="tight")
    save_png_figure(fig, path.with_suffix(".png"),dpi=300,bbox_inches="tight")
    plt.close(fig)

def unavailable_figure(out, name, title, detail):
    """Replace a previous figure with an explicit unavailable state."""
    fig,ax=plt.subplots(figsize=(10,4));ax.axis("off")
    ax.set_title(title)
    ax.text(.5,.55,"Results not available",ha="center",fontsize=15)
    ax.text(.5,.35,detail,ha="center",va="center",wrap=True,fontsize=10)
    savefig(fig,Path(out)/name)

def fig1_conceptual(out):
    fig,ax=plt.subplots(figsize=(14,8));ax.set_xlim(0,14);ax.set_ylim(0,8);ax.axis("off")
    def box(x,y,w,h,title,body,rounding=.12):
        p=FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.03,rounding_size={rounding}",
                         fill=False,linewidth=1.5)
        ax.add_patch(p);ax.text(x+.18,y+h-.35,title,fontsize=13,fontweight="bold",va="top")
        ax.text(x+.18,y+h-.82,body,fontsize=10,va="top",wrap=True)
    def arr(a,b):
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=15,linewidth=1.3))
    ax.text(.2,7.65,"Figure 1. From population-average risk to individual molecular predictability",
            fontsize=16,fontweight="bold",va="top")
    box(.3,4.7,3.7,2.2,"A  Population-average prediction",
        "Conventional glmnet / LightGBM\nlearn one global mapping\nX → disease risk\n\nPopulation summary: AUC / C-index")
    box(5.0,4.7,3.7,2.2,"B  Molecular neighbourhood",
        "Each participant is embedded in an\nomics-defined neighbourhood.\nReference agents are evaluated for\nsupport, stability and reliability.")
    box(10.0,4.7,3.7,2.2,"C  Dual individual output",
        "Disease risk: 17%\nMolecular predictability: HIGH\n\nRisk answers: how high?\nPredictability answers: how reliable?")
    arr((4.05,5.8),(4.9,5.8));arr((8.75,5.8),(9.9,5.8))
    box(.3,1.0,4.0,2.2,"D  Locked validation design",
        "One outcome-blind outer test roster\nshared by conventional C1 and ABM.\nABM development is split internally into\nbuild / tune / gate / calibration roles.")
    box(5.0,1.0,4.0,2.2,"E  Primary predictability estimand",
        "At 20/40/50/60/80/100% coverage:\ncompare frozen-selector participants with\npredicted-risk–matched outcome-blind controls.\nPrimary: ΔBrier / ΔLogLoss AURC.")
    box(9.7,1.0,4.0,2.2,"F  Reference-agent specificity",
        "Twin distance: molecular-only (primary).\nClinical+omics geometry: sensitivity.\nTest whether cleaning helps local twin inference\nbut not global population regression.")
    arr((4.3,2.1),(4.9,2.1));arr((9.0,2.1),(9.6,2.1))
    ax.text(.3,.42,"Conceptual figure; values shown for illustration only. No test outcome enters individual selection.",
            fontsize=9)
    savefig(fig,Path(out)/"Fig1_conceptual_framework")

def _benchmark(perf,layer,subset="all_common",trait=None):
    """Extract conventional predictions evaluated on verified identical people.

    Native C1 holdout performance is deliberately not a fallback: sharing the
    outer split does not imply sharing the same complete-case participants.
    """
    if perf.empty:return pd.DataFrame()
    required={"subset","model","N","comparison_N","comparison_eid_hash",
              "same_people","status","Harrell_C","AUC_IPCW_10y"}
    if not required<=set(perf):
        raise ValueError("Same-person benchmark metadata missing: "+", ".join(sorted(required-set(perf))))
    d=perf[perf.subset.eq(subset)].copy()
    if "layer" in d:d=d[d.layer.eq(layer)]
    if trait is not None and "trait" in d:d=d[d.trait.eq(trait)]
    preferred="Pradeep-style / glmnet" if layer=="prot" else "All-metabolite / glmnet"
    models=["Model 0","Biomarkers","Combined"]
    names={"C1_"+m+"::"+preferred:m for m in models}
    d["model"]=d.model.astype(str).str.strip()
    q=d[d.model.isin(names)].copy()
    if len(q)!=3 or q.model.nunique()!=3:
        raise ValueError(f"Expected three unique benchmarks for {layer}/{subset}/{preferred}")
    if not q.same_people.map(lambda x:str(x).strip().lower() in ("true","1")).all():
        raise ValueError("Benchmark does not certify identical participants")
    if not q.status.isin(["evaluated","completed"]).all():
        raise ValueError("Benchmark evaluation is incomplete")
    n=pd.to_numeric(q.N,errors="coerce");comparison_n=pd.to_numeric(q.comparison_N,errors="coerce")
    hashes=q.comparison_eid_hash.fillna("").astype(str).str.strip()
    if (not np.isfinite(n).all() or (n<=0).any() or (n%1!=0).any()
            or n.nunique()!=1 or not n.eq(comparison_n).all()
            or hashes.nunique()!=1 or not hashes.str.fullmatch(r"[0-9a-fA-F]{64}").all()):
        raise ValueError("Benchmark participant counts or participant hashes differ")
    if "outer_roster_sha256" in q and q.outer_roster_sha256.nunique(dropna=False)!=1:
        raise ValueError("Benchmark outer rosters differ")
    rows=[]
    for _,r in q.iterrows():
        for metric,column in [("IPCW AUC","AUC_IPCW_10y"),("Harrell C","Harrell_C")]:
            value=float(r[column])
            if not np.isfinite(value) or not 0<=value<=1:
                raise ValueError(f"Invalid {column} for {layer}/{r.model}")
            row=r.to_dict()
            row.update(model=names[r.model],source_model=r.model,panel=preferred,
                       metric=metric,value=value,layer=layer,
                       stage="c1_correlate",family="C1 same-person holdout")
            rows.append(row)
    return pd.DataFrame(rows)

def _load_benchmark(root,layer,trait=None):
    candidates=[(Path(root)/layer/"c1_abm_comparison/same_subset_primary_performance.tsv.gz","all_common"),
                (Path(root)/layer/"c1_abm/baseline_comparison/same_subset_primary_performance.tsv.gz","all")]
    for path,subset in candidates:
        try:
            perf=read_table(path)
            if perf.empty:continue
            result=_benchmark(perf,layer,subset,trait)
            result["source_table"]=str(path)
            return result,"Verified identical participants within this omic layer"
        except (ValueError,KeyError) as exc:
            # A contradictory preferred table must not be hidden by a fallback.
            return pd.DataFrame(),f"Invalid same-person table {path.name}: {exc}"
    return pd.DataFrame(),"Completed same-person holdout results unavailable; training screening AUC is not validation"

def fig2_population(root,out,trait=None):
    rows=[];availability={}
    for layer in LAYERS:
        b,detail=_load_benchmark(root,layer,trait)
        availability[layer]=detail
        if len(b):rows.append(b)
    d=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
    d.attrs["availability"]=availability
    if d.empty:
        unavailable_figure(out,"Fig2_population_average_prediction","Figure 2. Population prediction",
                           "Completed validation on identical C1/ABM participants is unavailable.\nTraining screening AUCs are not validation estimates.")
        return d
    fig,axs=plt.subplots(1,2,figsize=(12,5.5))
    models=["Model 0","Biomarkers","Combined"];x=np.arange(len(LAYERS));w=.23
    labels=[]
    for layer,label in zip(LAYERS,["Proteomics","Metabolomics"]):
        q=d[d.layer.eq(layer)]
        labels.append(label+(f"\nN = {int(q.N.iloc[0]):,}" if len(q) else "\nNot available"))
    for ax,metric,ylabel in zip(axs,["Harrell C","IPCW AUC"],["Harrell C-index","10-year IPCW AUC"]):
        values=d[d.metric.eq(metric)]
        for j,m in enumerate(models):
            vals=[]
            for layer in LAYERS:
                q=values[(values.layer==layer)&(values.model==m)]
                vals.append(float(q.value.iloc[0]) if len(q) else np.nan)
            ax.bar(x+(j-1)*w,vals,width=w,label={"Model 0":"Clinical","Biomarkers":"Omics","Combined":"Clinical + omics"}[m])
        ax.set_xticks(x,labels);ax.set_ylabel(ylabel)
        finite=pd.to_numeric(values.value,errors="coerce").dropna()
        if len(finite):ax.set_ylim(min(.5,finite.min()-.03),min(1,finite.max()+.03))
    axs[0].legend(frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(1.08,-.28))
    fig.suptitle("Figure 2. Population prediction on common C1/ABM participants")
    fig.text(.5,.015,"Identical participants across models within each omic layer; omic layers have different cohorts.",ha="center",fontsize=9)
    fig.subplots_adjust(bottom=.24,wspace=.3)
    savefig(fig,Path(out)/"Fig2_population_average_prediction")
    return d

def fig3_predictability(root,out):
    allcurve=[];allsum=[];strata=[]
    for layer in ["prot","met"]:
        base=Path(root)/layer/"c1_abm"/"baseline_comparison"
        curve=read_table(base/"lockpaper_selector_curve.tsv.gz")
        summ=read_table(base/"lockpaper_selector_inference.tsv.gz")
        if len(curve):
            z=curve[curve.route=="risk_resid_error"].copy();z["layer"]=layer;allcurve.append(z)
        if len(summ):
            z=summ[summ.route=="risk_resid_error"].copy();z["layer"]=layer;allsum.append(z)
        met=read_table(Path(root)/layer/"c1_abm"/"abm_twin"/"test_report"/"metrics.tsv.gz")
        if len(met):
            z=met[(met.model=="elasticnet_full")&met.subset.isin(["selected","all","rejected"])].copy()
            z=z.drop_duplicates(["model","subset"]);z["layer"]=layer;strata.append(z)
    c=pd.concat(allcurve,ignore_index=True) if allcurve else pd.DataFrame()
    s=pd.concat(allsum,ignore_index=True) if allsum else pd.DataFrame()
    st=pd.concat(strata,ignore_index=True) if strata else pd.DataFrame()
    if c.empty:
        unavailable_figure(out,"Fig3_predictability_coverage","Figure 3. Frozen predictability",
                           "Frozen-selector holdout coverage and error estimates are unavailable.")
        return c,s
    fig,axs=plt.subplots(1,3,figsize=(16,6.7))
    if len(st):
        order=["selected","all","rejected"];x=np.arange(3);w=.34
        for j,layer in enumerate(["prot","met"]):
            vals=[]
            for subset in order:
                q=st[(st.layer==layer)&(st.subset==subset)]
                vals.append(float(q.AUC_IPCW.iloc[0]) if len(q) else np.nan)
            axs[0].bar(x+(j-.5)*w,vals,width=w,color=LAYER_COLORS[layer],label="Proteomics" if layer=="prot" else "Metabolomics")
        axs[0].set_xticks(x,["Selected","All","Rejected"])
        axs[0].set_ylabel("10-year IPCW AUC");axs[0].set_title("A  Frozen selector strata")
        axs[0].legend(frameon=False)
    for layer,g in c.groupby("layer"):
        g=g.sort_values("coverage_target");label="Proteomics" if layer=="prot" else "Metabolomics"
        axs[1].plot(100*g.coverage_target,g.delta_Brier,marker="o",color=LAYER_COLORS[layer],label=label)
        axs[2].plot(100*g.coverage_target,g.delta_LogLoss,marker="o",color=LAYER_COLORS[layer],label=label)
    for ax,title,ylabel in [(axs[1],"B  Error beyond risk composition","Δ IPCW Brier vs risk-matched"),
                            (axs[2],"C  Log-loss beyond risk composition","Δ IPCW log-loss vs risk-matched")]:
        ax.axhline(0,ls="--",linewidth=1);ax.set_xlabel("Target coverage (%)");ax.set_ylabel(ylabel);ax.set_title(title);ax.legend(frameon=False)
    fig.suptitle("Figure 3. Frozen individual predictability: discrimination and error",fontweight="bold")
    targets=sorted(pd.to_numeric(c.coverage_target,errors="coerce").dropna().unique())
    coverage_rows=[];coverage_labels=[]
    for layer,label in zip(LAYERS,["Proteomics","Metabolomics"]):
        g=c[c.layer.eq(layer)]
        if g.empty:continue
        row=[]
        for target in targets:
            q=g[np.isclose(pd.to_numeric(g.coverage_target,errors="coerce"),target)]
            if len(q)==1 and {"coverage","selected_N"}<=set(q):
                actual=float(q.coverage.iloc[0]);selected_n=float(q.selected_N.iloc[0])
                row.append(f"{100*actual:.1f}%\nN={int(selected_n):,}" if np.isfinite(actual) and np.isfinite(selected_n) else "Not available")
            else:row.append("Not available")
        coverage_rows.append(row);coverage_labels.append(label)
    ax_table=fig.add_axes([.40,.035,.57,.15]);ax_table.axis("off")
    if coverage_rows:
        table=ax_table.table(cellText=coverage_rows,rowLabels=coverage_labels,
                             colLabels=[f"Target {100*x:g}%" for x in targets],loc="center",cellLoc="center")
        table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.65)
    fig.text(.40,.225,"Actual population coverage and selected N at each target",fontsize=9)
    fig.subplots_adjust(bottom=.35,top=.86,wspace=.35)
    savefig(fig,Path(out)/"Fig3_predictability_coverage")
    return c,s


def _contrast(root,layer):
    d=read_table(Path(root)/layer/"c1_abm"/"abm_twin"/"test_report"/"contrasts.tsv.gz")
    if d.empty:return d
    required={"metric","contrast","subset","reference_subset","model","reference","delta","lower","upper"}
    if not required<=set(d):return pd.DataFrame()
    # Both local and global changes are conditional on the same selected people.
    z=d[(d.metric=="AUC_IPCW")&d.contrast.isin(CLEANING_CONTRASTS)
        &d.subset.eq("selected")&d.reference_subset.eq("selected")].copy()
    metrics=read_table(Path(root)/layer/"c1_abm/abm_twin/test_report/metrics.tsv.gz")
    if metrics.empty or not {"model","subset","N","AUC_IPCW"}<=set(metrics):return pd.DataFrame()
    rows=[]
    for _,r in z.iterrows():
        a=metrics[(metrics.model==r.model)&metrics.subset.eq("selected")]
        b=metrics[(metrics.model==r.reference)&metrics.subset.eq("selected")]
        if len(a)!=1 or len(b)!=1:continue
        a,b=a.iloc[0],b.iloc[0]
        if not np.isfinite(float(a.N)) or float(a.N)<=0 or a.N!=b.N:continue
        if any(col in metrics and a[col]!=b[col] for col in ("known_N","cases","controls","requested_N")):continue
        if not np.isclose(float(r.delta),float(a.AUC_IPCW)-float(b.AUC_IPCW),rtol=1e-7,atol=1e-10):continue
        row=r.to_dict();row.update(layer=layer,comparison_N=int(a.N),comparison_scope="selected ABM test participants")
        if "Harrell_C" in metrics:
            row["delta_Harrell_C"]=float(a.Harrell_C)-float(b.Harrell_C)
        # No Harrell interval is inferred from a stored Uno C interval.
        rows.append(row)
    z=pd.DataFrame(rows)
    if len(z) and (z.contrast.duplicated().any() or z.comparison_N.nunique()!=1):
        return pd.DataFrame()
    return z

def fig4_cleaning(root,out):
    z=pd.concat([_contrast(root,x) for x in ["prot","met"]],ignore_index=True)
    if z.empty:
        unavailable_figure(out,"Fig4_reference_agent_cleaning","Figure 4. Reference-agent cleaning",
                           "Paired local and global contrasts on the same selected participants are unavailable.")
        return z
    order=CLEANING_CONTRASTS
    labels={"twin_library_cleaning":"Molecular-only twin: clean − full",
            "twin_integrated_library_cleaning":"Clinical+omics twin: clean − full",
            "training_effect_same_subset":"Global elastic-net: clean − full",
            "pruning_specificity":"Global elastic-net: profile-pruned − random-pruned"}
    fig,ax=plt.subplots(figsize=(11,6))
    y=0;yt=[];yl=[]
    for layer in ["prot","met"]:
        for con in order:
            q=z[(z.layer==layer)&(z.contrast==con)]
            if not len(q):continue
            r=q.iloc[0]
            lo,hi,pt=float(r.lower),float(r.upper),float(r.delta)
            ax.hlines(y,lo,hi,linewidth=1.4,color=LAYER_COLORS[layer]);ax.plot(pt,y,"o",color=LAYER_COLORS[layer])
            yt.append(y);yl.append(("Prot: " if layer=="prot" else "Met: ")+labels[con]+f" (N={int(r.comparison_N):,})");y+=1
        y+=.6
    ax.axvline(0,ls="--",linewidth=1);ax.set_yticks(yt,yl);ax.set_xlabel("Paired ΔIPCW AUC (95% group-bootstrap CI)")
    ax.set_title("Figure 4. Reference-agent cleaning on the same selected participants")
    savefig(fig,Path(out)/"Fig4_reference_agent_cleaning")
    return z

def supplement_selectors(root,out):
    rows=[]
    for layer in ["prot","met"]:
        d=read_table(Path(root)/layer/"c1_abm"/"baseline_comparison"/"lockpaper_selector_inference.tsv.gz")
        if len(d): d=d.copy();d["layer"]=layer;rows.append(d)
    z=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
    if z.empty:
        unavailable_figure(out,"FigS_selector_multiroute","Supplementary Figure. Selector comparison",
                           "Held-out selector inference is unavailable.")
        return z
    piv=z.pivot(index="route",columns="layer",values="AURC_delta_Brier")
    fig,ax=plt.subplots(figsize=(8,7));im=ax.imshow(piv.fillna(0).to_numpy(),aspect="auto")
    ax.set_yticks(range(len(piv)),piv.index);ax.set_xticks(range(len(piv.columns)),["Metabolomics" if x=="met" else "Proteomics" for x in piv.columns])
    for i in range(len(piv)):
        for j in range(len(piv.columns)):
            v=piv.iloc[i,j];ax.text(j,i,"NA" if not np.isfinite(v) else f"{v:.4f}",ha="center",va="center",fontsize=8)
    ax.set_title("Supplementary Figure. Selector-specific error reduction vs risk-matched controls")
    fig.colorbar(im,ax=ax,label="ΔBrier AURC (negative is better)")
    savefig(fig,Path(out)/"FigS_selector_multiroute")
    return z

def manuscript_grouping_audit(root):
    rows = []
    for layer in ("prot", "met"):
        path = Path(root) / layer / "c1_abm/abm_twin/host_contract.json"
        if path.is_file():
            # Describe the actual frozen fit, never a newly created global roster.
            audit = grouping_audit(json.loads(path.read_text()))
            rows.append(dict(layer=layer, **audit))
    return pd.DataFrame(rows)


def reporting_status(root,bench,curve,summ,clean,selectors):
    """Describe evidence availability independently of intermediate filenames."""
    rows=[]
    modules=[("population_prediction",bench,6),("predictability_curve",curve,6),
             ("predictability_inference",summ,1),("reference_cleaning",clean,4),
             ("selector_sensitivity",selectors,1)]
    for layer in LAYERS:
        for module,frame,minimum in modules:
            q=frame[frame.layer.eq(layer)] if "layer" in frame else pd.DataFrame()
            status="complete" if len(q)>=minimum else ("incomplete" if len(q) else "not_available")
            detail="Frozen validation estimates available" if status=="complete" else "Completed validation estimates unavailable"
            if module=="population_prediction":
                detail=bench.attrs.get("availability",{}).get(layer,detail)
                if status!="complete" and (Path(root)/layer/"c1_correlate/prediction/prediction.results.xlsx").is_file():
                    detail+="; an intermediate prediction workbook does not establish completed validation"
            if module=="predictability_curve" and len(q):
                if not {"coverage","selected_N"}<=set(q) or q[["coverage","selected_N"]].isna().any().any():
                    status="incomplete";detail="Actual coverage or selected participant counts unavailable"
            rows.append(dict(layer=layer,module=module,status=status,rows=len(q),detail=detail))
    return pd.DataFrame(rows)

def tables(out,bench,curve,summ,clean,selectors,grouping=None,reporting=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    with pd.ExcelWriter(out/"C1_LOCKPAPER_main_tables.xlsx",engine="openpyxl") as w:
        if reporting is not None:
            reporting.to_excel(w,sheet_name="Table0_reporting_status",index=False)
        bench.to_excel(w,sheet_name="Table1_population_prediction",index=False)
        summ.to_excel(w,sheet_name="Table2_predictability",index=False)
        clean.to_excel(w,sheet_name="Table3_reference_cleaning",index=False)
        curve.to_excel(w,sheet_name="Fig3_source",index=False)
        selectors.to_excel(w,sheet_name="Supplement_selectors",index=False)
        if grouping is not None:
            grouping.to_excel(w,sheet_name="Supplement_kinship_audit",index=False)
    bench.to_csv(out/"Table1_population_prediction.tsv",sep="\t",index=False)
    summ.to_csv(out/"Table2_predictability.tsv",sep="\t",index=False)
    clean.to_csv(out/"Table3_reference_cleaning.tsv",sep="\t",index=False)
    if reporting is not None:
        reporting.to_csv(out/"Table0_reporting_status.tsv",sep="\t",index=False)

RELIABILITY_TABLES = ("reliability_status", "reliability_population_performance", "risk_coverage", "aurc_contrasts",
    "reliability_increment_50", "risk_standardized_selection", "support_trend", "reliability_prevalent_support_audit")


def reliability_tables(scope):
    """Only completed, checksum-verified aggregate files enter publication reports."""
    directory=Path(scope)/"reliability"/"FINAL"
    marker=directory/"COMPLETE.json"
    if not marker.is_file():return {},{}
    try:
        receipt=json.loads(marker.read_text());files=receipt.get("files",{})
        contract=directory/"reliability_contract.json"
        requested=[contract,*[directory/(name+".tsv") for name in RELIABILITY_TABLES if (directory/(name+".tsv")).is_file()]]
        for file in requested:
            if files.get(file.name)!=hashlib.sha256(file.read_bytes()).hexdigest():
                raise ValueError(f"Reliability aggregate integrity failed: {file}")
        meta=json.loads(contract.read_text())
        if meta.get("version")!="LE8-C1-reliability-2026-10-26-v2":return {},{}
        tables={name:pd.read_csv(directory/(name+".tsv"),sep="\t") for name in RELIABILITY_TABLES if (directory/(name+".tsv")).is_file()}
        return tables,meta
    except (OSError,json.JSONDecodeError,pd.errors.EmptyDataError) as error:
        raise ValueError(f"Cannot read reliability results: {directory}") from error


def reliability_summary_row(trait, tables, contract):
    row=dict(trait=trait,discovery_or_replication=contract.get("discovery_or_replication","unavailable"))
    status=tables.get("reliability_status",pd.DataFrame())
    for prefix,name in (("P0","P0_risk_only"),("P1","P1_incident_support"),("P2","P2_incident_plus_prevalent")):
        values=status[status["head"].eq(name)] if "head" in status else pd.DataFrame()
        row[prefix+"_status"]=values.iloc[0]["status"] if len(values) else "unavailable"
    primary=[]
    def extract(table,model,reference,metric,prefix):
        frame=tables.get(table,pd.DataFrame());result=None
        if {"model","reference","metric"}<=set(frame):
            z=frame[frame.model.eq(model)&frame.reference.eq(reference)&frame.metric.eq(metric)]
            if len(z):result=z.iloc[0]
        row[prefix]=float(result["difference"]) if result is not None else np.nan
        row[prefix+"_CI"]=f"[{result['lower']:.6g}, {result['upper']:.6g}]" if result is not None and np.isfinite(result.get("lower",np.nan)) and np.isfinite(result.get("upper",np.nan)) else "unavailable"
        return result
    for metric in ("Brier","LogLoss"):
        primary.append(extract("aurc_contrasts","P1_incident_support","P0_risk_only","AURC_"+metric,"P1_minus_P0_AURC_"+metric))
        primary.append(extract("reliability_increment_50","P1_incident_support","P0_risk_only",metric,"P1_minus_P0_50pct_"+metric))
    extract("aurc_contrasts","P2_incident_plus_prevalent","P1_incident_support","AURC_Brier","P2_minus_P1_AURC_Brier")
    trend=tables.get("support_trend",pd.DataFrame())
    for support,prefix in (("eligible_match_count","twin_count_Brier_slope"),("family_ess","family_ESS_Brier_slope")):
        z=trend[trend.support_metric.eq(support)&trend.error_metric.eq("Brier")] if {"support_metric","error_metric"}<=set(trend) else pd.DataFrame()
        row[prefix]=float(z.iloc[0].slope) if len(z) else np.nan
        row[prefix+"_CI"]=f"[{z.iloc[0].lower:.6g}, {z.iloc[0].upper:.6g}]" if len(z) and np.isfinite(z.iloc[0].lower) and np.isfinite(z.iloc[0].upper) else "unavailable"
    if contract.get("status")!="completed_primary_P1_vs_P0":row["primary_conclusion"]="primary_unavailable"
    elif any(r is None or r.get("inference_status")!="evaluated" for r in primary):row["primary_conclusion"]="insufficient_primary_inference"
    elif all(r["upper"]<0 for r in primary):row["primary_conclusion"]="lower_error_supported_across_primary_contrasts"
    else:row["primary_conclusion"]="primary_support_increment_not_consistently_established"
    return row


def write_pgs_cross_trait_summary(root, layer):
    """Summarize saved C1 PGS scans, including traits outside the current batch.

    The run/skip switch controls fitting only. Missing/empty/all-NA scans are
    availability records, never tested features or negative findings.
    """
    import os, tempfile, fcntl
    from final import PGS_SCANS, PGS_SUMMARY_COLUMNS, saved_pgs_tables
    if layer not in {"prot", "met"}:
        raise ValueError("PGS summary layer must be prot or met")
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    from analysis_layout import report_directory, summary_lock_directory
    output = report_directory(root)
    output.mkdir(parents=True, exist_ok=True)
    lockdir = summary_lock_directory(root)
    lockdir.mkdir(parents=True, exist_ok=True)
    def read(path):
        if not path.is_file() or not path.stat().st_size:
            return pd.DataFrame()
        try:
            return pd.read_csv(path)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()
    def publish(path, write):
        with tempfile.NamedTemporaryFile(dir=output, prefix=".pgs-summary-", delete=False) as handle:
            temp = Path(handle.name)
        try:
            write(temp)
            os.replace(temp, path)
        finally:
            temp.unlink(missing_ok=True)
    with (lockdir / ("pgs_summary_" + layer + ".lock")).open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        collected = {name: [] for name in ("summary", "status", "associations")}
        prefix = "pwas" if layer == "prot" else "mwas"
        for scope in sorted(root.glob("*/" + layer + "/c1_correlate")):
            scans = {}
            for role in PGS_SCANS:
                path = scope / (prefix + "_" + role + ".csv")
                d = read(path)
                if len(d):
                    d = d.assign(source_file=str(path.relative_to(root)))
                scans[role] = d
            saved = saved_pgs_tables(scans, read(scope / "c1.pgs_status.csv"))
            for name, d in saved.items():
                collected[name].append(d.assign(Y=scope.parent.parent.name, layer=layer))
        empty_columns = dict(summary=PGS_SUMMARY_COLUMNS,
            status=["status", "declared_status", "analyses_summarized", "detail"],
            associations=["term", "analysis", "cohort_scope", "p.value", "FDR", "source_file"])
        result = {}
        for name, frames in collected.items():
            nonempty = [d for d in frames if len(d)]
            d = pd.concat(nonempty, ignore_index=True) if nonempty else pd.DataFrame(columns=[*empty_columns[name], "Y", "layer"])
            result[name] = d
            target = output / (f"c1_pgs_cross_trait_{name}.{layer}.tsv")
            publish(target, lambda path, frame=d: frame.to_csv(path, sep="\t", index=False))
        body = '<!doctype html><meta charset="utf-8"><title>C1 saved PGS summary</title><h1>C1 saved PGS analyses</h1>'
        body += '<p>Saved estimates and original within-scan FDR; no refitting. Missing PGS is omitted, not a negative result. Original cohort and model definitions apply.</p>'
        body += result["status"].to_html(index=False, escape=True)
        if len(result["summary"]):
            body += '<h2>Completed analyses</h2>' + result["summary"].to_html(index=False, escape=True)
        target = output / f"c1_pgs_cross_trait_summary.{layer}.html"
        publish(target, lambda path: path.write_text(body, encoding="utf-8"))
        print(f"[LE8 C1] PGS summary: {len(result['summary'])} saved analyses retained; missing PGS skipped; {target}", flush=True)
        return result


def write_reliability_cross_trait_summary(root, layer):
    import os,tempfile,fcntl
    from analysis_layout import report_directory, summary_lock_directory
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    output=report_directory(root);output.mkdir(parents=True,exist_ok=True)
    lockdir=summary_lock_directory(root);lockdir.mkdir(parents=True,exist_ok=True)
    with (lockdir/("reliability_summary_"+layer+".lock")).open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        rows=[]
        for scope in sorted(root.glob("*/"+layer+"/c1_correlate")):
            tables,contract=reliability_tables(scope)
            if tables:rows.append(reliability_summary_row(scope.parent.parent.name,tables,contract))
        frame=pd.DataFrame(rows)
        if not rows:return frame
        target=output/("c1_reliability_cross_trait_summary."+layer+".tsv")
        with tempfile.NamedTemporaryFile(dir=output,prefix=".reliability-summary-",delete=False) as handle:temp=Path(handle.name)
        try:frame.to_csv(temp,sep="\t",index=False);os.replace(temp,target)
        finally:temp.unlink(missing_ok=True)
        return frame


def write_reliability_report(root, traits, layers, out):
    """Primary report for FINAL reliability; historical summaries remain supplementary."""
    import html
    out=Path(out);root=Path(root);collected={};summaries=[]
    for trait in traits:
        for layer in layers:
            tables,contract=reliability_tables(root/trait/layer/"c1_correlate")
            if not tables:continue
            summaries.append(dict(reliability_summary_row(trait,tables,contract),layer=layer))
            for name,frame in tables.items():collected.setdefault(name,[]).append(frame.assign(trait=trait,layer=layer))
    if not summaries:return False
    out.mkdir(parents=True,exist_ok=True)
    merged={k:pd.concat(v,ignore_index=True) for k,v in collected.items()}
    merged["cross_trait_summary"]=pd.DataFrame(summaries)
    with pd.ExcelWriter(out/"C1_FINAL_reliability_tables.xlsx",engine="openpyxl") as writer:
        for name,frame in merged.items():
            frame.to_excel(writer,sheet_name=name[:31],index=False)
            frame.to_csv(out/("c1."+name+".tsv"),sep="\t",index=False)
    body=['<!doctype html><html lang="en"><meta charset="utf-8"><title>C1 FINAL reliability</title>',
        '<style>body{font:15px system-ui;margin:32px;line-height:1.5}table{border-collapse:collapse;font-size:12px}td,th{padding:6px;border:1px solid #ddd}.table{overflow:auto}</style>',
        '<h1>C1 FINAL: reliability of frozen glmnet risk</h1>',
        '<p>Primary evidence: direct P1−P0 selected-set losses at 50% coverage and AURC on a common target grid. Negative differences favor molecular support. Intervals use paired family bootstrap conditional on frozen models and heads. Selected-vs-nonselected tables are supplementary.</p>']
    def section(title,name,head=None):
        frame=merged.get(name,pd.DataFrame())
        if head is not None and "head" in frame:frame=frame[frame["head"].eq(head)]
        if head is not None and "model" in frame:frame=frame[frame.model.eq(head)]
        body.extend(["<h2>"+html.escape(title)+"</h2>",'<div class="table">'+(frame.to_html(index=False,escape=True) if len(frame) else '<p>Unavailable; see head status.</p>')+'</div>'])
    section("1. All-person frozen glmnet performance","reliability_population_performance")
    section("2. P0 risk-only reliability","risk_coverage","P0_risk_only")
    section("3. P1 risk and incident molecular support","risk_coverage","P1_incident_support")
    section("4. Primary P1−P0: direct 50% contrast","reliability_increment_50","P1_incident_support")
    section("Primary and secondary target-grid AURC","aurc_contrasts")
    section("5. P2 risk, incident and prevalent support","risk_coverage","P2_incident_plus_prevalent")
    section("6. Secondary P2−P1: direct 50% contrast","reliability_increment_50","P2_incident_plus_prevalent")
    section("7. Continuous support-error trend, adjusted for frozen risk strata","support_trend")
    section("8. Discovery / replication and primary conclusions","cross_trait_summary")
    section("Availability of every head","reliability_status")
    body.append('<p>HTN is discovery. COPD is locked replication only with documented prior non-inspection; otherwise it remains an unverified replication candidate. No K/PCA/radius/architecture/coverage tuning using test outcomes.</p></html>')
    (out/"C1_FINAL_reliability.html").write_text("\n".join(body))
    (out/"C1_FINAL_REPORT.json").write_text(json.dumps(dict(primary_report="C1_FINAL_reliability.html",source="verified FINAL aggregate tables",risk_model_refit=False,scopes=summaries),indent=2,default=str)+"\n")
    return True


def article_main(argv=None):
    p=argparse.ArgumentParser();p.add_argument("--analysis-root",required=True);p.add_argument("--Y",default="cvd_cad");p.add_argument("--out")
    a=p.parse_args(argv);root=Path(a.analysis_root)/a.Y;out=Path(a.out) if a.out else root/"paper_c1"
    out.mkdir(parents=True,exist_ok=True)
    if write_reliability_report(a.analysis_root,[a.Y],LAYERS,out):
        print(out/"C1_FINAL_reliability.html")
        return 0
    fig1_conceptual(out)
    bench=fig2_population(root,out,trait=a.Y);curve,summ=fig3_predictability(root,out);clean=fig4_cleaning(root,out);selectors=supplement_selectors(root,out)
    grouping = manuscript_grouping_audit(root)
    reporting=reporting_status(root,bench,curve,summ,clean,selectors)
    tables(out,bench,curve,summ,clean,selectors,grouping,reporting)
    grouped = (len(grouping) == 2 and {"layer", "grouping_unit", "grouping_source", "family_disjoint_outer_split", "family_disjoint_ABM_roles"} <= set(grouping)
               and set(grouping.layer) == {"prot", "met"}
               and grouping.grouping_unit.eq("family").all()
               and grouping.grouping_source.eq(GROUPING_SOURCE).all()
               and grouping.family_disjoint_outer_split.eq(True).all()
               and grouping.family_disjoint_ABM_roles.eq(True).all())
    methods = (
        "Genetically related participants were grouped using the UK Biobank KING relatedness data. "
        "All participants belonging to the same kinship-connected component were assigned to the same "
        "outer data split, preventing close relatives from crossing the development-test boundary. "
        "Components were constructed from all released related pairs before intersection with the analytical cohort. "
        "Participants absent from the relatedness graph received singleton groups. No participants were removed "
        "because of relatedness. The same groups were used for ABM development-role assignment and "
        "kinship-component block bootstrap inference. These groups represent genetic connectivity, not known pedigrees."
        if grouped else "Kinship-aware Methods language is pending verified kinship-aware fits for both omic layers."
    )
    (out/"Methods_kinship.txt").write_text(methods + "\n")
    required=[f"{name}.{ext}" for name in ("Fig1_conceptual_framework","Fig2_population_average_prediction",
        "Fig3_predictability_coverage","Fig4_reference_agent_cleaning","FigS_selector_multiroute") for ext in ("pdf","png")]
    missing=[name for name in required if not (out/name).is_file()]
    complete=grouped and not missing and reporting.status.eq("complete").all()
    (out/"LOCKPAPER_REPORT.json").write_text(json.dumps(dict(release=RELEASE,reporting_version=REPORT_VERSION,
        status="complete" if complete else "incomplete",components=reporting.to_dict("records"),
        kinship_validation="complete" if grouped else "pending",
        kinship_grouping_source=GROUPING_SOURCE if grouped else "unverified",
        missing_outputs=missing,source="frozen aggregate outputs only; no fitting or method selection"),indent=2)+"\n")
    (out/"LOCKPAPER_RELEASE.txt").write_text(RELEASE+"\n")
    print(out)


def main():
    import sys
    actions = {"article": article_main, "compare": comparison_main}
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        print("C1 internal reporting: c1_report.py {article|compare} [options]")
        return 0
    action = sys.argv.pop(1)
    if action not in actions:
        raise SystemExit("Unknown C1 report action: " + action)
    return actions[action]()


if __name__ == "__main__":
    import sys
    sys.exit(main())
