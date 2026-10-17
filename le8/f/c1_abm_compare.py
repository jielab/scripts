"""FINAL fair comparison of conventional C1 prediction and TwinABM.

Every model delta is computed on identical held-out participants. Selective
prediction is summarized route-by-route across frozen coverage levels rather than
by comparing a small ABM subset with a full-population conventional C-index.
"""

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
from c1_selector_inference import infer_selector
from c1_lock_roster import RELEASE, eid_hash, grouping_audit
# C1_LOCKPAPER::2026-10-17.c1-lockpaper-v1
from sklearn.metrics import roc_auc_score
from abm import uno_c


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
	)
	row["AUC_IPCW"] = (
		float(roc_auc_score(y[ok], score[ok], sample_weight=w[ok])) if len(np.unique(y[ok])) == 2 else np.nan
	)
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
	c1 = pd.DataFrame(benchmark["predictions"])
	req = {"eid", "biom_set", "model", "time", "event", "score"}
	if c1.empty:
		return unavailable("C1_predictions_unavailable")
	if not req <= set(c1):
		raise ValueError("Incomplete C1 benchmark predictions")
	c1["eid"] = c1.eid.astype(str)
	if c1.duplicated(["eid", "biom_set", "model"]).any():
		raise ValueError("Duplicate C1 model/participant predictions")
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
				"twin_integrated_clean",
				"twin_integrated_full",
				"elasticnet_full",
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
				AUC_IPCW=mm["AUC_IPCW"],
				Uno_C_horizon=mm["Uno_C_horizon"],
				Brier_IPCW=mm["Brier_IPCW"],
				LogLoss_IPCW=mm["LogLoss_IPCW"],
				riskmatched_N=rr["N"],
				riskmatched_cases=rr["cases"],
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
				ABM_common_eid_hash=eid_hash(ids[np.isin(ids, c1.eid)]),
				common_holdout_eid_hash=eid_hash(joined.eid),
				outer_roster=benchmark.get("lockpaper_protocol", {}).get("outer_roster", ""),
				outer_roster_sha256=benchmark.get("lockpaper_protocol", {}).get("outer_roster_sha256", ""),
				horizon=config.horizon,
				primary_selective_metric="Delta IPCW Brier/log-loss versus predicted-risk-matched outcome-blind control across coverage; raw AURC secondary",
				secondary_metrics="same-subset IPCW AUC and horizon-truncated Uno C",
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
	return {
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
	pred = joblib.load(twin / required[1])
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
		if json.loads((target / "sources.json").read_text()) == signature:
			return
	tables = comparison_from_paths(benchmark, twin, selector_inference=True)
	write_tables(stage, tables, "ABM twin", args.biom, args.Y, workbook="C1_ABM_comparison.xlsx")
	json_write(stage / "sources.json", signature)
	finish(stage, dict(source="same-roster C1 versus ABM selective prediction", source_hashes=signature))
	publish(stage, target)
	print("[LE8 twin] DONE fair C1/ABM comparison", flush=True)
