"""GRID: frozen, individual-reference prediction with training-only reliability gates.

Public API
----------
fit_predict(data, feature_groups, config=None) -> result
predict_new(bundle, data, retain_matches=True) -> result
materialize_matches(result, arm=None, query_indices=None) -> pandas.DataFrame

The caller owns phenotype construction, the common 50:50 outer roster, SNP
annotation, formal file publication and final test evaluation.  This module never
reads a test outcome.  It adapts LE8's reference-borrowing and selective-prediction
principles. The default selective-attention backend adapts LE8's module-token
Transformer/QK donor attention to genetic features and OOF residual values.

All learning uses the development half.  Its default roles are build/tune_model/
tune_gate/calibration = .60/.15/.15/.10.  Calibration is split once into threshold
fit and independent audit.  Only build participants form the donor bank (about
30% of the complete cohort).  Donor residuals are genuine family-held-out
predictions from baseline cross-fitting, not in-sample residuals.
"""

from __future__ import annotations

import hashlib
import math
import importlib.util
from pathlib import Path
import sys
from collections import Counter
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.neighbors import KDTree


# 🚩 Configuration and data contracts


DEFAULTS = {
	"abm_backend": "selective_attention",
	"device": "cuda",
	"attention_epochs": 20,
	"attention_patience": 4,
	"attention_batch": 256,
	"attention_width": 64,
	"attention_heads": 4,
	"attention_layers": 2,
	"attention_dropout": .1,
	"attention_lr": .0005,
	"attention_reconstruction": .1,
	"donor_block": 8192,
	"trait_type": "continuous",
	"seed": 20261007,
	"train_role_fractions": (0.60, 0.15, 0.15, 0.10),
	"folds": 5,
	"min_role_n": 20,
	"coverage": 0.50,
	"coverages": (0.20, 0.40, 0.50, 0.60, 0.80, 1.0),
	"primary_arm": "auto",
	"max_dims": 12,
	"ancestry_dims": 2,
	"csx_dims": 4,
	"frequency_dims": 2,
	"evolution_dims": 4,
	"ancestry_weight": 0.25,
	"metric_ridge": 10.0,
	"metric_weight_floor": 0.25,
	"k_grid": (16, 32, 64),
	"alpha_grid": (0.0, 0.25, 0.50, 0.75, 1.0),
	"radius_grid": (1.0, 2.0, 4.0),
	"prior_strength": 20.0,
	"radius_ref_n": 2048,
	"radius_neighbors": 5,
	"min_matches": 3,
	"min_ess": 2.0,
	"max_missing_fraction": 0.20,
	"csx_alpha_grid": (0.0, 1.0, 10.0, 100.0),
	"ancestry_min_n": 100,
	"ancestry_min_groups": 20,
	"ancestry_min_cases": 10,
	"ridge_alpha_grid": (1.0, 10.0, 100.0, 1000.0),
	"hgb_leaves_grid": (7, 15, 31),
	"hgb_iterations": 150,
	"hgb_l2": 10.0,
	"hgb_min_leaf": 30,
	"gate_iterations": 100,
	"gate_min_leaf": 20,
	"retrieval": "cuda",
	"query_batch": 1024,
	"ann_candidate_multiplier": 4,
	"ann_ef": 256,
	"ann_m": 24,
	"ann_audit_n": 64,
	"ann_recall_min": 0.95,
	"audit_bootstrap": 500,
	"audit_min_n": 30,
	"audit_min_groups": 20,
	"audit_min_cases": 10,
	"retain_match_arms": "primary",
	"enforce_half_split": True,
	"half_split_tolerance": 0.02,
	"strict_new_ids": True,
	"full_training_global_controls": True,
	"verbose": True,
}


def _attention():
	key = "grid_grid_attention"
	if key not in sys.modules:
		spec = importlib.util.spec_from_file_location(key, Path(__file__).with_name("attention.py"))
		module = importlib.util.module_from_spec(spec)
		sys.modules[key] = module
		try:
			spec.loader.exec_module(module)
		except Exception:
			sys.modules.pop(key, None)
			raise
	return sys.modules[key]


def preflight(config):
	c = _configuration(config)
	if c["abm_backend"] == "selective_attention" or c["retrieval"] == "cuda":
		return _attention().preflight(c["device"])
	return {"abm_backend": "reference", "device": "cpu", "neural_training": False}


def _log(config, status, stage, detail=""):
	if config.get("verbose", True):
		print(f"[GRID] {status} {stage}" + (f" | {detail}" if detail else ""), flush=True)


def _seed(value, seed):
	return int.from_bytes(hashlib.blake2b(
		(f"{seed}:{value}").encode(), digest_size=8
	).digest(), "little")


def _unique(values):
	return list(dict.fromkeys(values))


def _configuration(config):
	c = {**DEFAULTS, **(config or {})}
	if c["abm_backend"] not in {"reference", "selective_attention"}:
		raise ValueError("abm_backend must be reference or selective_attention")
	if c["device"] != "cpu" and not (c["device"] == "cuda" or c["device"].startswith("cuda:")):
		raise ValueError("device must be cpu or cuda[:index]")
	for key in ("attention_epochs", "attention_patience", "attention_batch", "attention_width", "attention_heads", "attention_layers", "donor_block"):
		if int(c[key]) != c[key] or c[key] < 1:
			raise ValueError(f"{key} must be a positive integer")
	if c["attention_width"] % c["attention_heads"] or not 0 <= c["attention_dropout"] < 1:
		raise ValueError("Attention width must be divisible by heads; dropout must be in [0,1)")
	if not np.isfinite(c["attention_lr"]) or c["attention_lr"] <= 0 or not np.isfinite(c["attention_reconstruction"]) or c["attention_reconstruction"] < 0:
		raise ValueError("Invalid attention learning rate or reconstruction weight")
	if c["retrieval"] == "cuda" and c["device"] == "cpu":
		raise ValueError("Explicit CPU execution requires --retrieval kd_tree or hnsw")
	if c["trait_type"] not in {"continuous", "binary"}:
		raise ValueError("trait_type must be continuous or binary")
	c["train_role_fractions"] = tuple(map(float, c["train_role_fractions"]))
	if len(c["train_role_fractions"]) != 4 or min(c["train_role_fractions"]) <= 0 or not np.isclose(sum(c["train_role_fractions"]), 1):
		raise ValueError("Four positive train_role_fractions summing to one are required")
	for name in ("k_grid", "csx_alpha_grid", "ridge_alpha_grid", "hgb_leaves_grid"):
		c[name] = tuple(sorted(set(c[name])))
	for name in ("alpha_grid", "radius_grid", "coverages"):
		c[name] = tuple(sorted(set(map(float, c[name]))))
	if not 0 < c["coverage"] <= 1 or not c["coverages"] or not all(0 < q <= 1 for q in c["coverages"]):
		raise ValueError("Coverage must be in (0,1]")
	c["coverages"] = tuple(sorted(set(c["coverages"]) | {float(c["coverage"])}))
	if not c["k_grid"] or min(c["k_grid"]) < 1 or any(int(k) != k for k in c["k_grid"]):
		raise ValueError("k_grid must contain positive integers")
	if not c["alpha_grid"] or not all(0 <= x <= 1 for x in c["alpha_grid"]):
		raise ValueError("alpha_grid must lie in [0,1]")
	if 0.0 not in c["alpha_grid"]:
		c["alpha_grid"] = (0.0, *c["alpha_grid"])
	if not c["radius_grid"] or min(c["radius_grid"]) <= 0:
		raise ValueError("radius_grid must be positive")
	if not 1 <= c["max_dims"] <= 12:
		raise ValueError("This implementation supports an explicit retrieval dimension cap of 1..12")
	budgets = [c[name] for name in ("ancestry_dims", "csx_dims", "frequency_dims", "evolution_dims")]
	if any(int(v) != v or v < 0 for v in budgets) or sum(budgets) > c["max_dims"]:
		raise ValueError("Fixed ancestry/csx/frequency/evolution dimension budgets must be nonnegative integers summing to <=max_dims; default 2+4+2+4=12")
	if c["folds"] < 2 or c["min_role_n"] < 2:
		raise ValueError("At least two cross-fitting folds and two participants per role are required")
	if c["retrieval"] not in {"kd_tree", "hnsw", "cuda"}:
		raise ValueError("retrieval must be kd_tree, hnsw or cuda")
	if c["prior_strength"] < 0 or c["query_batch"] < 1 or c["min_ess"] < 1:
		raise ValueError("Invalid borrowing or batching configuration")
	if not 0 <= c["max_missing_fraction"] < 1:
		raise ValueError("max_missing_fraction must lie in [0,1)")
	if not 0 < c["ann_recall_min"] <= 1 or c["ann_candidate_multiplier"] < 1:
		raise ValueError("Invalid ANN recall/candidate setting")
	if c["audit_bootstrap"] < 20:
		raise ValueError("At least 20 audit bootstrap draws are required")
	return c


def _feature_contract(feature_groups):
	allowed = {"covariates", "csx", "ancestry", "evolution", "frequency", "permuted_evolution", "disco"}
	unknown = set(feature_groups) - allowed
	if unknown:
		raise ValueError(f"Unknown feature groups: {sorted(unknown)}")
	groups = {}
	for key in allowed:
		value = feature_groups.get(key, [])
		if isinstance(value, str) or not isinstance(value, (list, tuple)):
			raise TypeError(f"feature_groups[{key}] must be a list of column names")
		groups[key] = _unique(map(str, value))
	if not groups["csx"]:
		raise ValueError("At least one CSx source score is required")
	forbidden = {"y", "Y", "outcome", "event", "time", "split", "role", "eid", "family_id", "ancestry"}
	for key, columns in groups.items():
		if forbidden & set(columns):
			raise ValueError(f"Non-feature columns in {key}: {sorted(forbidden & set(columns))}")
	return groups


def _ids(frame, column):
	if column not in frame or frame[column].isna().any():
		raise ValueError(f"Complete {column} is required")
	values = frame[column].astype(str).to_numpy()
	if np.any(np.char.str_len(np.char.strip(values.astype(str))) == 0):
		raise ValueError(f"Blank {column}")
	return values


def _validate_features(frame, groups):
	columns = _unique(x for names in groups.values() for x in names)
	missing = set(columns) - set(frame.columns)
	if missing:
		raise ValueError(f"Missing predictor columns: {sorted(missing)}")
	for name in columns:
		if not pd.api.types.is_numeric_dtype(frame[name]):
			raise ValueError(f"Predictor {name} must be numeric; encode categories before calling GRID")
	return columns


def _roles(frame, config):
	groups = _ids(frame, "family_id")
	counts = Counter(groups)
	levels = sorted(counts, key=lambda g: _seed(g, config["seed"] + 11))
	boundaries = np.cumsum(config["train_role_fractions"])
	role_names = np.array(["build", "tune_model", "tune_gate", "calibration"])
	lookup = {}
	position = 0
	for g in levels:
		midpoint = (position + counts[g] / 2) / len(frame)
		lookup[g] = role_names[min(int(np.searchsorted(boundaries, midpoint)), 3)]
		position += counts[g]
	result = np.array([lookup[g] for g in groups])
	for name in role_names:
		if np.sum(result == name) < config["min_role_n"]:
			raise ValueError(f"Insufficient independent data for internal role {name}")
	return result


def _folds(groups, count, seed):
	levels = sorted(set(groups), key=lambda g: _seed(g, seed))
	if len(levels) < count:
		raise ValueError("Too few independent families for cross-fitting")
	lookup = {g: j % count for j, g in enumerate(levels)}
	return np.array([lookup[g] for g in groups], dtype=int)


def _labels(frame, kind, context, require_two_classes=True):
	# Called only on an explicitly selected development role; never on outer test.
	y = pd.to_numeric(frame["y"], errors="coerce").to_numpy(float)
	if not np.isfinite(y).all():
		raise ValueError(f"{context}: missing/nonfinite development labels; filter eligibility before splitting")
	if kind == "binary" and not np.isin(y, [0, 1]).all():
		raise ValueError(f"{context}: binary labels must be 0/1")
	if kind == "binary" and require_two_classes and len(np.unique(y)) < 2:
		raise ValueError(f"{context}: both outcome classes are required")
	return y


# 🚩 Training-only baseline models and out-of-fold donor residuals


@dataclass
class NumericDesign:
	columns: list[str]

	def fit(self, frame):
		x = self.raw(frame)
		with np.errstate(invalid="ignore"):
			self.median = np.array([
				np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else 0.
				for v in x.T
			])
		filled = np.where(np.isfinite(x), x, self.median)
		self.center = filled.mean(0)
		self.scale = np.maximum(filled.std(0), 1e-8)
		return self

	def raw(self, frame):
		if not self.columns:
			return np.zeros((len(frame), 1), dtype=float)
		return frame[self.columns].to_numpy(float)

	def transform(self, frame):
		x = self.raw(frame)
		return (np.where(np.isfinite(x), x, self.median) - self.center) / self.scale


class FrozenPredictor:
	def __init__(self, columns, kind, model_kind="ridge", parameter=1., config=None):
		self.columns = list(columns)
		self.kind, self.model_kind, self.parameter = kind, model_kind, parameter
		self.config = dict(config or DEFAULTS)

	def fit(self, frame, y):
		self.design = NumericDesign(self.columns).fit(frame)
		x = self.design.transform(frame)
		if self.model_kind == "hgb":
			common = dict(max_leaf_nodes=int(self.parameter), max_iter=self.config["hgb_iterations"],
				min_samples_leaf=min(self.config["hgb_min_leaf"], max(2, len(frame) // 10)),
				l2_regularization=self.config["hgb_l2"], learning_rate=0.05,
				early_stopping=False, random_state=self.config["seed"])
			self.model = (HistGradientBoostingClassifier(**common) if self.kind == "binary"
				else HistGradientBoostingRegressor(**common))
		elif self.kind == "binary":
			self.model = LogisticRegression(C=1. / max(float(self.parameter), 1e-6),
				solver="lbfgs", max_iter=3000, random_state=self.config["seed"])
		else:
			self.model = Ridge(alpha=float(self.parameter), solver="lsqr", tol=1e-8)
		self.model.fit(x, y)
		return self

	def predict(self, frame):
		x = self.design.transform(frame)
		return (self.model.predict_proba(x)[:, 1] if self.kind == "binary"
			else self.model.predict(x))

	def decompose(self, frame):
		"""Exact additive baseline terms, on the logit scale for binary outcomes."""
		if self.model_kind != "ridge":
			raise ValueError("An additive decomposition is defined only for the linear/logistic baseline")
		x = self.design.transform(frame)
		coef = np.asarray(self.model.coef_, float).reshape(-1)
		intercept = float(np.asarray(self.model.intercept_).reshape(-1)[0])
		terms = x * coef
		linear = intercept + terms.sum(axis=1)
		prediction = 1. / (1. + np.exp(-np.clip(linear, -700., 700.))) if self.kind == "binary" else linear
		if not np.allclose(prediction, self.predict(frame), rtol=1e-8, atol=1e-10):
			raise AssertionError("Baseline contribution reconstruction failed")
		out = pd.DataFrame({"baseline_intercept": np.full(len(frame), intercept),
			"baseline_linear_predictor": linear,
			"baseline_link": "logit" if self.kind == "binary" else "identity"})
		for j, column in enumerate(self.columns):
			out["baseline_contribution." + column] = terms[:, j]
		return out


class AncestryPredictor:
	"""A strong CSx comparator: one score-combination model per target ancestry.

	Every group uses the same validation-selected regularization value. A pooled
	model is available for small groups or unseen ancestry labels. Group labels
	are used here for the established comparator, not to define GRID neighbors.
	"""
	def __init__(self, columns, kind, parameter=1., config=None):
		self.columns, self.kind = list(columns), kind
		self.model_kind, self.parameter = "ancestry", parameter
		self.config = dict(config or DEFAULTS)

	def fit(self, frame, y):
		self.pooled = FrozenPredictor(self.columns, self.kind, "ridge", self.parameter, self.config).fit(frame, y)
		self.models, self.records = {}, []
		ancestry = _ids(frame, "ancestry")
		family = _ids(frame, "family_id")
		for name in sorted(set(ancestry)):
			ix = np.flatnonzero(ancestry == name)
			n, ng = len(ix), len(set(family[ix]))
			reason = ""
			if n < self.config["ancestry_min_n"]:
				reason = "insufficient_build_n"
			elif ng < self.config["ancestry_min_groups"]:
				reason = "insufficient_independent_build_groups"
			elif self.kind == "binary" and min(np.sum(y[ix] == 1), np.sum(y[ix] == 0)) < self.config["ancestry_min_cases"]:
				reason = "insufficient_build_cases_or_controls"
			if not reason:
				self.models[name] = FrozenPredictor(self.columns, self.kind, "ridge", self.parameter, self.config).fit(frame.iloc[ix], y[ix])
			self.records.append(dict(ancestry=name, build_n=n, build_groups=ng,
				build_cases=int(np.sum(y[ix] == 1)) if self.kind == "binary" else np.nan,
				build_controls=int(np.sum(y[ix] == 0)) if self.kind == "binary" else np.nan,
				status="pooled_fallback" if reason else "ancestry_specific",
				fallback_reason=reason or "none", shared_alpha=float(self.parameter)))
		return self

	def predict(self, frame):
		prediction = self.pooled.predict(frame)
		ancestry = _ids(frame, "ancestry")
		for name, model in self.models.items():
			ix = np.flatnonzero(ancestry == name)
			if len(ix):
				prediction[ix] = model.predict(frame.iloc[ix])
		return prediction

	def group_audit(self, tune, y):
		records = {row["ancestry"]: dict(row) for row in self.records}
		ancestry = _ids(tune, "ancestry")
		for name in sorted(set(ancestry) | set(records)):
			row = records.setdefault(name, dict(ancestry=name, build_n=0, build_groups=0,
				build_cases=0 if self.kind == "binary" else np.nan,
				build_controls=0 if self.kind == "binary" else np.nan,
				status="pooled_fallback", fallback_reason="no_build_participants",
				shared_alpha=float(self.parameter)))
			ix = np.flatnonzero(ancestry == name)
			row.update(tune_n=len(ix),
				tune_cases=int(np.sum(y[ix] == 1)) if self.kind == "binary" else np.nan,
				tune_controls=int(np.sum(y[ix] == 0)) if self.kind == "binary" else np.nan)
		return pd.DataFrame(records.values())

	def decompose(self, frame):
		out = self.pooled.decompose(frame)
		ancestry = _ids(frame, "ancestry")
		for name, model in self.models.items():
			ix = np.flatnonzero(ancestry == name)
			if len(ix):
				out.iloc[ix] = model.decompose(frame.iloc[ix]).to_numpy()
		return out


def _make_predictor(columns, kind, model_kind, parameter, config):
	if model_kind == "ancestry":
		return AncestryPredictor(columns, kind, parameter, config)
	return FrozenPredictor(columns, kind, model_kind, parameter, config)


def _mse(y, prediction):
	return float(np.mean((np.asarray(y) - np.asarray(prediction)) ** 2))


def _tune_predictor(name, columns, model_kind, parameters, build, tune, yb, yt, config):
	best = None
	rows = []
	for value in parameters:
		model = _make_predictor(columns, config["trait_type"], model_kind, value, config).fit(build, yb)
		score = _mse(yt, model.predict(tune))
		rows.append(dict(stage="global_model", model=name, parameter=float(value),
			validation_loss=score, loss="Brier" if config["trait_type"] == "binary" else "MSE"))
		if best is None or score < best[0] - 1e-12:
			best = (score, model, value)
	for row in rows:
		row["selected"] = row["parameter"] == float(best[2])
	return best[1], rows


def _oof_residuals(build, y, baseline, config):
	groups = _ids(build, "family_id")
	fold = _folds(groups, config["folds"], config["seed"] + 71)
	prediction = np.full(len(build), np.nan)
	audit = []
	for j in range(config["folds"]):
		fit = fold != j
		held = fold == j
		if config["trait_type"] == "binary" and len(np.unique(y[fit])) != 2:
			raise ValueError("An OOF training fold has only one class; fewer folds or more events are required")
		fit_indices = np.flatnonzero(fit)
		inner_count = min(5, len(set(groups[fit])))
		if inner_count < 2:
			raise ValueError("Too few independent families for nested donor-baseline tuning")
		inner = _folds(groups[fit], inner_count, config["seed"] + 72 + j)
		inner_fit = fit_indices[inner != 0]
		inner_tune = fit_indices[inner == 0]
		if config["trait_type"] == "binary" and len(np.unique(y[inner_fit])) != 2:
			raise ValueError("A nested donor-baseline training fold has only one class")
		candidate, _ = _tune_predictor("donor_nested_CSx", baseline.columns,
			baseline.model_kind, config["csx_alpha_grid"], build.iloc[inner_fit],
			build.iloc[inner_tune], y[inner_fit], y[inner_tune], config)
		model = _make_predictor(baseline.columns, baseline.kind, baseline.model_kind,
			candidate.parameter, config).fit(build.iloc[fit_indices], y[fit])
		prediction[held] = model.predict(build.iloc[np.flatnonzero(held)])
		if set(groups[fit]) & set(groups[held]):
			raise AssertionError("A family crossed a donor OOF fold")
		audit.append(dict(fold=j, fit_n=int(fit.sum()), donor_n=int(held.sum()),
			fit_groups=len(set(groups[fit])), donor_groups=len(set(groups[held])),
			oof_loss=_mse(y[held], prediction[held]), family_overlap=False,
			selected_alpha=float(candidate.parameter), inner_tune_n=len(inner_tune),
			heldout_outcome_used_for_tuning=False))
	if not np.isfinite(prediction).all():
		raise AssertionError("Incomplete out-of-fold donor baseline predictions")
	return y - prediction, prediction, fold, audit


# 🚩 Compact supervised matching coordinates with explicit feature blocks


class MatchGeometry:
	def fit(self, build, residual, blocks, config):
		self.blocks = {}
		self.config = dict(config)
		self.dimensions = int(config["max_dims"])
		names = [name for name, cols in blocks.items() if cols]
		# Fixed slots ensure that adding a frequency block never squeezes away
		# evolution dimensions. Missing blocks leave their own slots at zero.
		allocation, starts = {}, {}
		start = 0
		for name in ("ancestry", "csx", "frequency", "evolution"):
			allocation[name] = int(config[name + "_dims"])
			starts[name] = start
			start += allocation[name]
		allocation["permuted_evolution"] = allocation["evolution"]
		starts["permuted_evolution"] = starts["evolution"]
		effective = 0
		audit = []
		for name in names:
			start = starts[name]
			if allocation[name] < 1:
				raise ValueError(f"A positive dimension budget is required for supplied block {name}")
			columns = list(blocks[name])
			design = NumericDesign(columns).fit(build)
			x = design.transform(build)
			keep = np.std(x, axis=0) > 1e-8
			if not keep.any():
				self.blocks[name] = dict(design=design, keep=keep, empty=True,
					start=start, stop=start, columns=columns)
				for column in columns:
					audit.append(dict(block=name, feature=column, metric_weight=0.,
						status="constant_in_build"))
				continue
			x = x[:, keep]
			# A small, explicit ridge metric learns useful axes from donor OOF
			# residuals. Query outcomes are not accepted by this transform.
			teacher = Ridge(alpha=config["metric_ridge"], solver="lsqr", tol=1e-8).fit(x, residual)
			strength = np.asarray(teacher.coef_, float) ** 2
			if strength.mean() > 1e-16:
				strength /= strength.mean()
			else:
				strength.fill(0.)
			weights = np.sqrt(np.clip(config["metric_weight_floor"] + strength, 0.05, 5.0))
			x = x * weights
			n = min(allocation[name], x.shape[1], max(1, len(x) - 1))
			pca = PCA(n_components=n, svd_solver="full" if x.shape[1] < 100 else "randomized",
				random_state=config["seed"] + 83).fit(x)
			z = pca.transform(x)
			scale = max(float(np.sqrt(np.mean(np.sum(z * z, axis=1)))), 1e-8)
			block_weight = config["ancestry_weight"] if name == "ancestry" else 1.
			self.blocks[name] = dict(design=design, keep=keep, empty=False,
				weights=weights, pca=pca, scale=scale, block_weight=block_weight,
				start=start, stop=start+n, columns=columns)
			effective += n
			weight_map = dict(zip(np.array(columns)[keep], weights))
			for column in columns:
				audit.append(dict(block=name, feature=column, metric_weight=float(weight_map.get(column, 0.)),
					status="retained" if column in weight_map else "constant_in_build"))
		self.effective_dimensions = effective
		self.audit = pd.DataFrame(audit)
		if not effective:
			raise ValueError("No variable matching coordinates in the build set")
		return self

	def transform(self, frame):
		z = np.zeros((len(frame), self.dimensions), dtype=np.float64)
		for name, block in self.blocks.items():
			if block["empty"]:
				continue
			x = block["design"].transform(frame)[:, block["keep"]] * block["weights"]
			z[:, block["start"]:block["stop"]] = (
				block["pca"].transform(x) / block["scale"] * math.sqrt(block["block_weight"])
			)
		if not np.isfinite(z).all():
			raise ValueError("Nonfinite matching coordinates")
		return z

	def component_distances(self, query_z, bank_z, indices):
		out = {}
		safe = np.maximum(indices, 0)
		for name, b in self.blocks.items():
			if b["empty"]:
				out[name] = np.zeros(indices.shape, dtype=float)
			else:
				lo, hi = b["start"], b["stop"]
				delta = query_z[:, None, lo:hi] - bank_z[safe, lo:hi]
				out[name] = np.sqrt(np.sum(delta * delta, axis=2))
				out[name][indices < 0] = np.nan
		return out


# 🚩 Indexed retrieval, absolute training radius and donor-wise decomposition


class ReferenceBank:
	def fit(self, build, residual, baseline_oof, folds, blocks, config):
		self.config = dict(config)
		self.attention = None
		self.ids = _ids(build, "eid")
		self.groups = _ids(build, "family_id")
		self.ancestry = build["ancestry"].astype(str).to_numpy()
		self.residual = np.asarray(residual, float)
		self.baseline_oof = np.asarray(baseline_oof, float)
		self.observed_y = build["y"].to_numpy(float).copy()
		self.donor_folds = np.asarray(folds, int)
		self.geometry = MatchGeometry().fit(build, residual, blocks, config)
		self.z = self.geometry.transform(build)
		self.tie = np.array([_seed(v, config["seed"] + 89) for v in self.ids], dtype=np.uint64)
		self.group_counts = Counter(self.groups)
		self._tree = None
		self._ann = None
		self._count_tree = None
		self._count_group_order = None
		self._make_index()
		selected = np.argsort(self.tie, kind="stable")[:min(config["radius_ref_n"], len(build))]
		neighbors = min(config["radius_neighbors"], max(1, len(build) - max(self.group_counts.values())))
		jj, dd = self.query(self.z[selected], self.ids[selected], self.groups[selected], neighbors)
		kth = dd[:, -1]
		usable = kth[np.isfinite(kth) & (kth > 1e-12)]
		# Zero distances alone do not demonstrate a supported genetic geometry.
		self.radius = float(np.median(usable)) if len(usable) else 0.
		self.radius_status = "available" if self.radius > 0 else "constant_or_unavailable"
		self.retrieval_audit = {
			"backend": config["retrieval"], "donor_n": len(self.ids),
			"donor_groups": len(self.group_counts), "dimension_budget": config["max_dims"],
			"effective_dimensions": self.geometry.effective_dimensions,
			"radius_reference_n": len(selected), "radius": self.radius,
			"radius_status": self.radius_status, "self_family_excluded": True,
			"eligible_count_method": "exact_cuda_blocked_radius_count" if config["retrieval"] == "cuda" else "exact_KDTree_radius_count_minus_same_ID_or_family",
			"matched_count_scope": "retrieved_donors_used_up_to_k",
		}
		if config["retrieval"] == "hnsw":
			self._audit_ann(selected[:config["ann_audit_n"]])
		return self

	def __getstate__(self):
		return {k: v for k, v in vars(self).items()
			if k not in {"_tree", "_ann"} and not k.startswith("_count_")}

	def __setstate__(self, state):
		self.__dict__.update(state)
		self._tree = None
		self._ann = None
		self._count_tree = None
		self._count_group_order = None

	def _make_index(self):
		if self.config["retrieval"] == "cuda":
			return
		if self.config["retrieval"] == "kd_tree":
			if self._tree is None:
				self._tree = KDTree(self.z, leaf_size=40, metric="euclidean")
		elif self._ann is None:
			try:
				import hnswlib
			except ImportError as exc:
				raise ImportError("Optional retrieval=hnsw requires hnswlib; use kd_tree or install it explicitly") from exc
			self._ann = hnswlib.Index(space="l2", dim=self.z.shape[1])
			self._ann.init_index(max_elements=len(self.ids), ef_construction=max(200, self.config["ann_ef"]),
				M=self.config["ann_m"], random_seed=self.config["seed"])
			self._ann.add_items(self.z.astype(np.float32), np.arange(len(self.ids)), num_threads=1)
			self._ann.set_ef(max(self.config["ann_ef"], max(self.config["k_grid"]) * self.config["ann_candidate_multiplier"]))
			self._ann.set_num_threads(1)

	def _make_count_index(self):
		# Exact support counts remain exact even when nearest-neighbor retrieval
		# uses HNSW. Reuse the existing KDTree when that is the search backend.
		if self._count_tree is None:
			if self.config["retrieval"] == "kd_tree":
				self._make_index()
				self._count_tree = self._tree
			else:
				self._count_tree = KDTree(self.z, leaf_size=40, metric="euclidean")
		if self._count_group_order is None:
			self._count_group_order = np.argsort(self.groups, kind="stable")
			ordered = self.groups[self._count_group_order]
			self._count_groups, self._count_starts, self._count_sizes = np.unique(
				ordered, return_index=True, return_counts=True)
			self._count_id_order = np.argsort(self.ids, kind="stable")
			self._count_ids = self.ids[self._count_id_order]

	def count_eligible(self, z, ids, groups, caliper, query_qc=None):
		"""Count every admissible donor, not just the first k retrieved.

		Only scalar radius counts are retained. Known same-ID/family members
		are subtracted using their indexed positions; no query-by-bank distance
		matrix or radius-neighbor list is materialized.
		"""
		ids, groups = np.asarray(ids, str), np.asarray(groups, str)
		qc = np.ones(len(z), dtype=bool) if query_qc is None else np.asarray(query_qc, bool)
		if ids.shape != (len(z),) or groups.shape != (len(z),) or qc.shape != (len(z),):
			raise ValueError("Eligible-count IDs, families and QC must match query rows")
		counts = np.zeros(len(z), dtype=np.int64)
		if caliper <= 0 or not qc.any():
			return counts
		if np.isnan(caliper):
			raise ValueError("Eligible-count caliper cannot be NaN")
		if self.config["retrieval"] == "cuda":
			rows = np.flatnonzero(qc)
			counts[rows] = _attention().neighbors(self.z, z[rows], self.ids, ids[rows],
				self.groups, groups[rows], self.tie, 0, self.config["device"],
				min(self.config["query_batch"], self.config["attention_batch"]), self.config["donor_block"], radius=caliper)[2]
			return counts
		self._make_count_index()
		for begin in range(0, len(z), self.config["query_batch"]):
			stop = min(len(z), begin + self.config["query_batch"])
			rows = np.flatnonzero(qc[begin:stop]) + begin
			if not len(rows):
				continue
			counts[rows] = self._count_tree.query_radius(z[rows], r=caliper, count_only=True)
			family_positions = np.searchsorted(self._count_groups, groups[rows])
			id_positions = np.searchsorted(self._count_ids, ids[rows])
			family_known = family_positions < len(self._count_groups)
			family_known[family_known] &= (
				self._count_groups[family_positions[family_known]] == groups[rows[family_known]])
			id_known = id_positions < len(self._count_ids)
			id_known[id_known] &= self._count_ids[id_positions[id_known]] == ids[rows[id_known]]
			# Ordinary new families/IDs never enter this loop. For overlapping
			# queries, compute distances only to their small excluded subset.
			for local in np.flatnonzero(family_known | id_known):
				row = rows[local]
				excluded = np.empty(0, dtype=np.int64)
				if family_known[local]:
					position = family_positions[local]
					start, size = self._count_starts[position], self._count_sizes[position]
					excluded = self._count_group_order[start:start + size]
				if id_known[local]:
					individual = self._count_id_order[id_positions[local]]
					# Self may already be among family exclusions; subtract once.
					if not family_known[local] or self.groups[individual] != groups[row]:
						excluded = np.append(excluded, individual)
				distance = np.sqrt(np.sum((self.z[excluded] - z[row]) ** 2, axis=1))
				counts[row] -= int(np.count_nonzero(distance <= caliper))
		if (counts < 0).any():
			raise AssertionError("Exclusion count exceeded exact radius support")
		return counts

	def _retrieve(self, z, count):
		self._make_index()
		if self.config["retrieval"] == "kd_tree":
			d, j = self._tree.query(z, k=count, return_distance=True, sort_results=True)
		else:
			j, _ = self._ann.knn_query(z.astype(np.float32), k=count, num_threads=1)
			# Recompute exact distances inside the ANN candidate pool.
			d = np.sqrt(np.sum((z[:, None, :] - self.z[j]) ** 2, axis=2))
		return j, d

	def query(self, z, ids, groups, k):
		k = min(int(k), len(self.ids))
		if k < 1:
			raise ValueError("Empty reference bank")
		if self.config["retrieval"] == "cuda":
			j, d, _ = _attention().neighbors(self.z, z, self.ids, ids, self.groups, groups,
				self.tie, k, self.config["device"], min(self.config["query_batch"], self.config["attention_batch"]), self.config["donor_block"])
			return j, d
		indices = np.full((len(z), k), -1, dtype=np.int32)
		distances = np.full((len(z), k), np.inf, dtype=np.float64)
		for begin in range(0, len(z), self.config["query_batch"]):
			stop = min(len(z), begin + self.config["query_batch"])
			block_groups = groups[begin:stop]
			excluded = max([self.group_counts.get(str(g), 0) for g in block_groups] + [0]) + 1
			count = min(len(self.ids), k + excluded)
			if self.config["retrieval"] == "hnsw":
				count = min(len(self.ids), max(count, k * self.config["ann_candidate_multiplier"]))
			j, d = self._retrieve(z[begin:stop], count)
			allowed = (self.ids[j] != ids[begin:stop, None]) & (self.groups[j] != block_groups[:, None])
			d[~allowed] = np.inf
			for row in range(stop - begin):
				order = np.lexsort((self.tie[j[row]], d[row]))[:k]
				take = np.isfinite(d[row, order])
				n = int(take.sum())
				indices[begin+row, :n] = j[row, order[take]]
				distances[begin+row, :n] = d[row, order[take]]
		return indices, distances

	def _audit_ann(self, sample):
		k = min(16, len(self.ids) - max(self.group_counts.values()))
		if k < 1 or not len(sample):
			raise ValueError("Insufficient references for ANN recall audit")
		q = self.z[sample]
		ann_j, ann_d = self.query(q, self.ids[sample], self.groups[sample], k)
		exact = KDTree(self.z, leaf_size=40)
		count = min(len(self.ids), k + max(self.group_counts.values()) + 1)
		dist, ix = exact.query(q, k=count)
		recall = []
		for row, s in enumerate(sample):
			good = (self.ids[ix[row]] != self.ids[s]) & (self.groups[ix[row]] != self.groups[s])
			expected = ix[row, good][:k]
			found = ann_j[row, np.isfinite(ann_d[row])]
			recall.append(len(set(expected) & set(found)) / max(len(expected), 1))
		value = float(np.mean(recall))
		self.retrieval_audit.update(ann_recall_at_k=value, ann_audit_k=k, ann_audit_n=len(sample),
			ann_audit_source="build_outcome_blind")
		if value < self.config["ann_recall_min"]:
			raise ValueError(f"ANN mean recall {value:.4f} is below the required {self.config['ann_recall_min']}; increase ann_ef or use kd_tree")

	def borrow(self, frame, baseline, *, k, alpha, radius_multiplier, retain=False, retrieved=None,
			eligible_counts=None):
		z = self.geometry.transform(frame)
		matching_columns = _unique(col for block in self.geometry.blocks.values() for col in block["columns"])
		missing_fraction = 1. - np.isfinite(frame[matching_columns].to_numpy(float)).mean(1)
		query_qc = missing_fraction <= self.config["max_missing_fraction"]
		if retrieved is None:
			j, d = self.query(z, _ids(frame, "eid"), _ids(frame, "family_id"), k)
		else:
			j, d = retrieved[0][:, :k], retrieved[1][:, :k]
		valid = (j >= 0) & np.isfinite(d)
		valid &= query_qc[:, None]
		caliper = self.radius * radius_multiplier
		valid &= (d <= caliper) if caliper > 0 else False
		if eligible_counts is None:
			eligible_count = self.count_eligible(z, _ids(frame, "eid"), _ids(frame, "family_id"), caliper, query_qc)
		else:
			eligible_count = np.asarray(eligible_counts)
			if (eligible_count.shape != (len(frame),) or not np.isfinite(eligible_count).all()
					or (eligible_count < 0).any() or (eligible_count > len(self.ids)).any()
					or (eligible_count != np.floor(eligible_count)).any()):
				raise ValueError("Cached eligible counts must be integer counts aligned to the query rows")
			eligible_count = np.where(query_qc, eligible_count, 0).astype(np.int64)
		safe = np.maximum(j, 0)
		logits = np.where(valid, -0.5 * (d / max(self.radius, 1e-12)) ** 2, -np.inf)
		maximum = np.max(logits, axis=1, keepdims=True)
		maximum[~np.isfinite(maximum)] = 0.
		weights = np.exp(logits - maximum)
		denom = weights.sum(1, keepdims=True)
		weights = np.divide(weights, denom, out=np.zeros_like(weights), where=denom > 0)
		if getattr(self, "attention", None) is not None:
			weights = self.attention.weights(frame, j, d, valid, self.radius)
		ss = (weights * weights).sum(1)
		ess = np.divide(1., ss, out=np.zeros(len(frame)), where=ss > 0)
		family_ess = ess.copy()
		family_count = (weights > 0).sum(1)
		# Most query banks have no repeated donor family. Only recompute those
		# rows that need group aggregation, rather than looping over every person.
		_, group_codes = np.unique(self.groups, return_inverse=True)
		codes = np.where(weights > 0, group_codes[safe], -1)
		ordered = np.sort(codes, axis=1)
		repeated = np.any((ordered[:, 1:] == ordered[:, :-1]) & (ordered[:, 1:] >= 0), axis=1)
		for row in np.flatnonzero(repeated):
			good = weights[row] > 0
			if not good.any():
				continue
			codes, inv = np.unique(self.groups[safe[row, good]], return_inverse=True)
			mass = np.bincount(inv, weights=weights[row, good])
			family_count[row] = len(codes)
			family_ess[row] = 1. / np.sum(mass * mass)
		nearest = np.min(np.where(valid, d, np.inf), axis=1)
		farthest = np.max(np.where(valid, d, 0.), axis=1)
		count = valid.sum(1)
		if (eligible_count < count).any():
			raise AssertionError("Used donor count exceeded exact eligible support")
		support = np.exp(-0.5 * (nearest / max(self.radius, 1e-12)) ** 2)
		enough = (count >= self.config["min_matches"]) & (family_ess >= self.config["min_ess"])
		fraction = alpha * support * family_ess / (family_ess + self.config["prior_strength"])
		fraction[~enough] = 0.
		contributions = fraction[:, None] * weights * self.residual[safe]
		raw = np.asarray(baseline) + contributions.sum(1)
		prediction = np.clip(raw, 0., 1.) if self.config["trait_type"] == "binary" else raw
		clipping = prediction - raw
		if not np.allclose(prediction, np.asarray(baseline) + contributions.sum(1) + clipping, rtol=1e-10, atol=1e-10):
			raise AssertionError("Donor contributions do not reconstruct the delivered prediction")
		nearest[~np.isfinite(nearest)] = np.nan
		farthest[count == 0] = np.nan
		reason = np.select(
			[~query_qc, count == 0, count < self.config["min_matches"], family_ess < self.config["min_ess"],
				np.full(len(frame), alpha == 0)],
			["excess_missing_predictors", "no_donor_within_caliper", "too_few_matched_donors", "low_effective_family_n",
				"tuning_selected_baseline"],
			default="local_correction",
		)
		diagnostics = pd.DataFrame({
			"matched_count": count, "matched_family_count": family_count,
			"eligible_match_count": eligible_count, "matches_truncated": eligible_count > count,
			"retrieval_k": np.full(len(frame), j.shape[1], dtype=np.int64),
			"donor_ess": ess, "family_ess": family_ess,
			"max_donor_weight": weights.max(1), "nearest_distance": nearest,
			"farthest_distance": farthest, "nearest_radius_ratio": nearest / max(self.radius, 1e-12),
			"support": support, "borrow_fraction": fraction,
			"weighted_residual": (weights * self.residual[safe]).sum(1),
			"local_correction": contributions.sum(1), "clipping_correction": clipping,
			"matching_supported": enough, "fallback_reason": reason,
			"matching_missing_fraction": missing_fraction, "matching_qc_pass": query_qc,
			"caliper": np.full(len(frame), caliper),
		})
		packed = None
		if retain:
			actual = np.where(valid, j, -1).astype(np.int32)
			distances = np.where(valid, d, np.nan)
			packed = {
				"query_eids": _ids(frame, "eid"), "donor_eids": self.ids.copy(),
				"donor_family_ids": self.groups.copy(), "donor_ancestry": self.ancestry.copy(),
				"donor_oof_baseline": self.baseline_oof.copy(),
				"donor_observed_y": self.observed_y.copy(),
				"donor_oof_residual": self.residual.copy(), "donor_fold": self.donor_folds.copy(),
				"neighbor_indices": actual, "distance": distances, "weights": weights,
				"weighted_residual_contribution": contributions,
				"baseline": np.asarray(baseline).copy(), "prediction": prediction.copy(),
				"clipping_correction": clipping, "borrow_fraction": fraction,
				"eligible_match_count": eligible_count, "matches_truncated": eligible_count > count,
				"retrieval_k": int(j.shape[1]),
				"distance_components": self.geometry.component_distances(z, self.z, actual),
			}
		return prediction, diagnostics, packed


def _tune_bank(name, bank, tune, baseline, y, config):
	z = bank.geometry.transform(tune)
	maximum = min(max(config["k_grid"]), len(bank.ids))
	ids, families = _ids(tune, "eid"), _ids(tune, "family_id")
	j, d = bank.query(z, ids, families, maximum)
	matching_columns = _unique(col for block in bank.geometry.blocks.values() for col in block["columns"])
	query_qc = (1. - np.isfinite(tune[matching_columns].to_numpy(float)).mean(1)
		<= config["max_missing_fraction"])
	# Support depends on radius, never k or alpha. Reuse one exact count per
	# radius across the entire hyperparameter grid, retaining only N scalars.
	eligible_by_radius = {
		radius: bank.count_eligible(z, ids, families, bank.radius * radius, query_qc)
		for radius in config["radius_grid"]
	}
	rows = []
	best = None
	for k in config["k_grid"]:
		if k > maximum:
			continue
		for radius in config["radius_grid"]:
			# Distances and donor weights do not depend on alpha. Reuse one local
			# decomposition for this k/radius instead of repeating every retrieval.
			_, unit, _ = bank.borrow(tune, baseline, k=int(k), alpha=1.,
				radius_multiplier=radius, retrieved=(j, d), eligible_counts=eligible_by_radius[radius])
			correction = unit.local_correction.to_numpy()
			for alpha in config["alpha_grid"]:
				pred = baseline + alpha * correction
				if config["trait_type"] == "binary":
					pred = np.clip(pred, 0., 1.)
				score = _mse(y, pred)
				rows.append(dict(stage="matching", model=name, k=int(k), alpha=float(alpha),
					radius_multiplier=float(radius), validation_loss=score,
					loss="Brier" if config["trait_type"] == "binary" else "MSE",
					matched_fraction=float((unit.matched_count > 0).mean())))
				# Grids are sorted: exact ties select the smaller, more conservative setting.
				if best is None or score < best[0] - 1e-12:
					best = (score, dict(k=int(k), alpha=float(alpha), radius_multiplier=float(radius)))
	if best is None:
		raise ValueError("No neighbor setting fits the build bank; lower k_grid")
	for row in rows:
		row["selected"] = all(row[k] == v for k, v in best[1].items())
	return best[1], rows


# 🚩 Frozen gain/absolute-error gate and independent internal policy audit


class ReliabilityGate:
	def fit(self, table, y, primary, config):
		self.columns = list(table.columns)
		self.design = NumericDesign(self.columns).fit(table)
		x = self.design.transform(table)
		baseline = table["CSx"].to_numpy(float)
		self.models = {}
		targets = {
			"gain": (y - baseline) ** 2 - (y - primary) ** 2,
			"absolute_error": (y - primary) ** 2,
		}
		for name, target in targets.items():
			model = HistGradientBoostingRegressor(max_leaf_nodes=5,
				learning_rate=0.04, max_iter=config["gate_iterations"],
				min_samples_leaf=min(config["gate_min_leaf"], max(2, len(y) // 10)),
				l2_regularization=10., early_stopping=False,
				random_state=config["seed"] + 233)
			model.fit(x, target)
			self.models[name] = model
		return self

	def predict(self, table):
		x = self.design.transform(table)
		return {
			"gain": self.models["gain"].predict(x),
			"absolute_error": np.maximum(self.models["absolute_error"].predict(x), 0.),
		}


class CoverageRule:
	def fit(self, score, ids, coverage, seed):
		score = np.asarray(score, float)
		self.coverage, self.seed = float(coverage), int(seed)
		finite = np.flatnonzero(np.isfinite(score))
		if not len(finite):
			raise ValueError("No finite training-derived reliability scores")
		ties = np.array([_seed(v, seed) for v in ids], dtype=np.uint64)
		order = finite[np.lexsort((ties[finite], -score[finite]))]
		j = order[max(0, math.ceil(len(order) * coverage) - 1)]
		self.threshold, self.tie = float(score[j]), int(ties[j])
		return self

	def apply(self, score, ids):
		score = np.asarray(score, float)
		if self.coverage == 1:
			return np.isfinite(score)
		ties = np.array([_seed(v, self.seed) for v in ids], dtype=np.uint64)
		return np.isfinite(score) & (
			(score > self.threshold) | ((score == self.threshold) & (ties <= self.tie))
		)


def _gate_table(table, primary):
	suffix = "evolution" if "Ridge_evolution" in table else "no_evolution"
	ridge, hgb = "Ridge_" + suffix, "HGB_" + suffix
	columns = ["CSx", ridge, hgb, primary]
	z = table[columns].copy()
	z["primary_minus_csx"] = table[primary] - table["CSx"]
	z["ridge_minus_hgb"] = table[ridge] - table[hgb]
	for name in ["matched_count", "matched_family_count", "eligible_match_count", "matches_truncated", "retrieval_k", "donor_ess", "family_ess",
			"max_donor_weight", "nearest_radius_ratio", "support", "borrow_fraction",
			"local_correction", "missing_fraction", "matching_supported"]:
		z[name] = table[name].astype(float)
	for name in ["matched_count", "matched_family_count", "eligible_match_count", "donor_ess", "family_ess"]:
		z[name] = np.log1p(z[name])
	return z


def _raw_predictions(bundle, frame, retain_matches=False):
	c = bundle["config"]
	ids = _ids(frame, "eid")
	result = {
		"eid": ids, "family_id": _ids(frame, "family_id"),
		"ancestry": frame["ancestry"].astype(str).to_numpy(),
	}
	for name, model in bundle["models"].items():
		result[name] = model.predict(frame)
	for name, values in bundle["models"]["CSx"].decompose(frame).items():
		result[name] = values.to_numpy()
	result["csx_baseline_source"] = np.where(
		frame["ancestry"].astype(str).isin(bundle["models"]["CSx"].models),
		"ancestry_specific", "pooled_fallback",
	)
	packed = {}
	primary = bundle["primary_arm"]
	for arm, bank in bundle["banks"].items():
		retain = retain_matches and (
			c["retain_match_arms"] == "all" or arm == primary or
			isinstance(c["retain_match_arms"], (list, tuple)) and arm in c["retain_match_arms"]
		)
		pred, diag, evidence = bank.borrow(frame, np.asarray(result["CSx"]),
			**bundle["bank_parameters"][arm], retain=retain)
		result[arm] = pred
		for column in diag:
			result[f"{arm}.{column}"] = diag[column].to_numpy()
			if arm == primary:
				result[column] = diag[column].to_numpy()
		if evidence is not None:
			packed[arm] = evidence
	columns = _unique(bundle["feature_groups"]["covariates"] + [
		col for block in bundle["banks"][primary].geometry.blocks.values() for col in block["columns"]
	])
	x = frame[columns].to_numpy(float)
	result["missing_fraction"] = 1. - np.isfinite(x).mean(1) if columns else 0.
	result["technical_qc_pass"] = result["missing_fraction"] <= c["max_missing_fraction"]
	return pd.DataFrame(result), packed


def _apply_policy(bundle, raw):
	table = raw.copy()
	gate = bundle["gate"].predict(_gate_table(table, bundle["primary_arm"]))
	ids = table["eid"].to_numpy(str)
	c = bundle["config"]
	table["expected_gain_over_csx"] = gate["gain"]
	table["expected_squared_error"] = gate["absolute_error"]
	table["expected_rmse"] = np.sqrt(gate["absolute_error"])
	ref = np.asarray(bundle.get("confidence_reference", []), dtype=float)
	if len(ref):
		pos = np.searchsorted(ref, gate["absolute_error"], side="right")
		table["confidence_percentile"] = 1.0 - pos / float(len(ref))
	else:
		table["confidence_percentile"] = np.nan
	audit_ok = all(bundle.get(name, {}).get("status") == "supported_lower_group_error"
		for name in ("low_error_audit", "confidence_audit"))
	table["confidence_validated"] = bool(audit_ok)
	table["confidence_candidate"] = (
		table["technical_qc_pass"].to_numpy(bool)
		& table["matching_supported"].to_numpy(bool)
		& (table["confidence_percentile"].to_numpy(float) >= 0.80)
	)
	table["high_confidence"] = audit_ok & table["confidence_candidate"]
	table["confidence_label"] = np.select(
		[~np.full(len(table), audit_ok, dtype=bool),
		 table["high_confidence"].to_numpy(bool),
		 table["confidence_percentile"].to_numpy(float) >= 0.50],
		["unvalidated", "high", "moderate"], default="low"
	)

	table["low_error_screen_audit_status"] = bundle.get("low_error_audit", {}).get("status", "not_evaluated")
	table["gain_audit_scope"] = "pooled_internal_calibration; not a per-person guarantee"
	ancestry_audits = bundle.get("ancestry_gain_audit", {})
	table["gain_audit_ancestry_status"] = table.ancestry.astype(str).map(
		{name: value["status"] for name, value in ancestry_audits.items()}).fillna("not_evaluated")
	table["research_selected"] = bundle["coverage_rules"][c["coverage"]].apply(gate["gain"], ids)
	table["selected"] = table["research_selected"] & table["technical_qc_pass"]
	table["selected_absolute_error"] = bundle["absolute_error_rule"].apply(-gate["absolute_error"], ids) & table["technical_qc_pass"]
	table["candidate"] = table["selected"] & (gate["gain"] > 0) & table["matching_supported"].to_numpy(bool)
	table["released"] = table["candidate"] & (bundle["audit"]["status"] == "supported_gain")
	table["GRID_policy"] = np.where(table["released"], table[bundle["primary_arm"]], table["CSx"])
	table["policy_source"] = np.where(table["released"], bundle["primary_arm"], "CSx")
	table["selection_reason"] = np.select(
		[~table["technical_qc_pass"], ~table["selected"], gate["gain"] <= 0, ~table["matching_supported"].to_numpy(bool), ~table["released"]],
		["technical_predictor_QC", "outside_fixed_coverage", "nonpositive_expected_gain", "insufficient_matching_support", "internal_audit_not_supportive"],
		default="released",
	)
	for coverage, rule in bundle["coverage_rules"].items():
		table[f"selected.q{coverage:g}"] = rule.apply(gate["gain"], ids) & table["technical_qc_pass"]
	for label, rules in bundle["control_rules"].items():
		if label == "absolute_error":
			score = -gate["absolute_error"]
		elif label == "support_only":
			score = table["support"].to_numpy() * table["family_ess"].to_numpy()
		elif label == "clinical_lowrisk":
			score = -table["covariate_baseline"].to_numpy()
		else:
			score = np.array([_seed(v, c["seed"] + 313) / 2**64 for v in ids])
		for coverage, rule in rules.items():
			table[f"selected.{label}.q{coverage:g}"] = rule.apply(score, ids) & table["technical_qc_pass"]
	return table


def _internal_audit(bundle, frame, y):
	raw, _ = _raw_predictions(bundle, frame)
	table = _apply_policy(bundle, raw)
	mask = table["candidate"].to_numpy(bool)
	groups = table["family_id"].to_numpy(str)[mask]
	c = bundle["config"]
	delta = ((y - table[bundle["primary_arm"]].to_numpy()) ** 2
		- (y - table["CSx"].to_numpy()) ** 2)[mask]
	levels, inverse = np.unique(groups, return_inverse=True)
	result = {
		"status": "insufficient_audit_information", "source": "independent_calibration_audit",
		"primary_arm": bundle["primary_arm"], "reference": "CSx",
		"candidate_n": int(mask.sum()), "candidate_groups": len(levels),
		"delta_loss": float(np.mean(delta)) if len(delta) else np.nan,
		"lower": np.nan, "upper": np.nan, "bootstrap": c["audit_bootstrap"],
		"metric": "Brier" if c["trait_type"] == "binary" else "MSE",
		"test_outcomes_used": False,
	}
	enough = len(delta) >= c["audit_min_n"] and len(levels) >= c["audit_min_groups"]
	if c["trait_type"] == "binary":
		cases = int(y[mask].sum())
		result.update(cases=cases, controls=int(mask.sum()) - cases)
		enough &= min(cases, int(mask.sum()) - cases) >= c["audit_min_cases"]
	if enough:
		group_sum = np.bincount(inverse, weights=delta)
		group_n = np.bincount(inverse)
		rng = np.random.default_rng(c["seed"] + 337)
		draws = []
		for _ in range(c["audit_bootstrap"]):
			sample = rng.integers(0, len(levels), size=len(levels))
			draws.append(group_sum[sample].sum() / group_n[sample].sum())
		result["lower"], result["upper"] = map(float, np.quantile(draws, [.025, .975]))
		result["status"] = "supported_gain" if result["upper"] < 0 else "gain_not_established"
	return result


def _internal_low_error_audit(bundle, frame, y, selector="selected_absolute_error"):
	"""Evaluate the already frozen low-error selector against all audit people.

	This is a second, descriptive calibration audit. It never changes the gate,
	coverage thresholds, model or gain policy. It supports a group-average
	error statement, not a guarantee for any individual.
	"""
	raw, _ = _raw_predictions(bundle, frame)
	table = _apply_policy(bundle, raw)
	mask = table[selector].to_numpy(bool)
	loss = (np.asarray(y) - table[bundle["primary_arm"]].to_numpy(float)) ** 2
	levels, inverse = np.unique(table.family_id.astype(str), return_inverse=True)
	n, ns = len(mask), int(mask.sum())
	c = bundle["config"]
	result = {"status": "insufficient_audit_information", "source": "independent_calibration_audit",
		"selector": "frozen_confidence_top20_qc_supported" if selector == "confidence_candidate" else "frozen_expected_squared_error", "reference": "all_audit_individuals_same_model",
		"n": n, "selected_n": ns, "coverage": ns / n if n else np.nan,
		"selected_loss": float(loss[mask].mean()) if ns else np.nan,
		"all_loss": float(loss.mean()) if n else np.nan,
		"rejected_loss": float(loss[~mask].mean()) if n > ns else np.nan,
		"delta_selected_minus_all": float(loss[mask].mean() - loss.mean()) if ns else np.nan,
		"lower": np.nan, "upper": np.nan, "bootstrap": c["audit_bootstrap"],
		"metric": "Brier" if c["trait_type"] == "binary" else "MSE",
		"test_outcomes_used": False, "changes_prediction_policy": False,
		"interpretation": "selected-group average error; not individual accuracy or pure genetic accuracy"}
	ngs = len(np.unique(inverse[mask]))
	ngr = len(np.unique(inverse[~mask]))
	result.update(selected_families=ngs, rejected_families=ngr)
	enough = min(ns, n - ns) >= c["audit_min_n"] and min(ngs, ngr) >= c["audit_min_groups"]
	if c["trait_type"] == "binary":
		cases = int(np.asarray(y)[mask].sum())
		result.update(selected_cases=cases, selected_controls=ns-cases,
			case_coverage=cases / float(np.asarray(y).sum()) if np.asarray(y).sum() else np.nan)
		enough &= min(cases, ns-cases) >= c["audit_min_cases"]
	if enough:
		counts = np.bincount(inverse, minlength=len(levels))
		sums = np.bincount(inverse, weights=loss, minlength=len(levels))
		selected_counts = np.bincount(inverse, weights=mask.astype(float), minlength=len(levels))
		selected_sums = np.bincount(inverse, weights=loss * mask, minlength=len(levels))
		rng = np.random.default_rng(c["seed"] + 347)
		draws = []
		for _ in range(c["audit_bootstrap"]):
			ix = rng.integers(0, len(levels), size=len(levels))
			denom = selected_counts[ix].sum()
			if denom:
				draws.append(selected_sums[ix].sum()/denom - sums[ix].sum()/counts[ix].sum())
		if len(draws) >= 20:
			result["lower"], result["upper"] = map(float, np.quantile(draws, [.025, .975]))
			result["status"] = "supported_lower_group_error" if result["upper"] < 0 else "lower_error_not_established"
	return result


# 🚩 Public fit/predict interface: the outer test label is never inspected


def fit_predict(data: pd.DataFrame, feature_groups: dict, config: dict | None = None) -> dict:
	c = _configuration(config)
	preflight(c)
	groups = _feature_contract(feature_groups)
	if c["primary_arm"] == "auto":
		c["primary_arm"] = "GRID_evolution" if groups["evolution"] else "GRID_frequency_only" if groups["frequency"] else "GRID_no_evolution"
	if not isinstance(data, pd.DataFrame) or "split" not in data or "y" not in data:
		raise ValueError("data must include split and y, plus eid/family_id/ancestry")
	ids = _ids(data, "eid")
	family = _ids(data, "family_id")
	if len(set(ids)) != len(ids):
		raise ValueError("Duplicate participant IDs")
	_ids(data, "ancestry")
	all_features = _validate_features(data, groups)
	split = data["split"].astype(str).replace({"training": "train", "testing": "test"}).to_numpy()
	if set(split) != {"train", "test"}:
		raise ValueError("A fixed common train/test outer split is required")
	if set(family[split == "train"]) & set(family[split == "test"]):
		raise ValueError("A family crosses the outer train/test split")
	train_fraction = float(np.mean(split == "train"))
	if c["enforce_half_split"] and abs(train_fraction - .5) > c["half_split_tolerance"]:
		raise ValueError(f"Expected a 50:50 outer roster; observed training fraction {train_fraction:.4f}")
	development = data.loc[split == "train"].sort_values("eid", kind="stable").reset_index(drop=True).copy()
	# Test data are restricted to predictor/identity columns immediately.  No
	# access to test y occurs anywhere below, including split or gate construction.
	test_columns = _unique(["eid", "family_id", "ancestry", *all_features])
	test = data.loc[split == "test", test_columns].reset_index(drop=True).copy()
	y_development = _labels(development, c["trait_type"], "development")
	for column in all_features:
		v = development[column].to_numpy(float)
		if np.isfinite(v).all() and np.array_equal(v, y_development):
			raise ValueError(f"Predictor {column} is identical to development y; possible outcome leakage")
	role = _roles(development, c)
	part = {name: np.flatnonzero(role == name) for name in ("build", "tune_model", "tune_gate", "calibration")}
	frame = {name: development.iloc[ix].reset_index(drop=True) for name, ix in part.items()}
	y = {name: _labels(frame[name], c["trait_type"], name) for name in ("build", "tune_model", "tune_gate")}
	build, tune = frame["build"], frame["tune_model"]
	_log(c, "START", "global_models", f"build={len(build)} tune_model={len(tune)} outer_test={len(test)}")
	baseline_columns = _unique(groups["covariates"] + groups["csx"])
	full_columns = _unique(groups["covariates"] + groups["csx"] + groups["ancestry"]
		+ groups["frequency"] + groups["evolution"])
	no_evolution_columns = _unique(groups["covariates"] + groups["csx"] + groups["ancestry"] + groups["frequency"])
	models = {}
	tuning = []
	models["covariate_baseline"] = FrozenPredictor(groups["covariates"], c["trait_type"],
		"ridge", 0., c).fit(build, y["build"])
	model_specs = [
		("CSx_pooled", baseline_columns, "ridge", c["csx_alpha_grid"]),
		("CSx", baseline_columns, "ancestry", c["csx_alpha_grid"]),
		("Ridge_no_evolution", no_evolution_columns, "ridge", c["ridge_alpha_grid"]),
		("HGB_no_evolution", no_evolution_columns, "hgb", c["hgb_leaves_grid"]),
	]
	if groups["evolution"]:
		model_specs += [("Ridge_evolution", full_columns, "ridge", c["ridge_alpha_grid"]),
			("HGB_evolution", full_columns, "hgb", c["hgb_leaves_grid"])]
	if groups["disco"]:
		model_specs.append(("DiscoDivas_calibrated_reference", _unique(groups["covariates"] + groups["disco"]),
			"ridge", c["csx_alpha_grid"]))
	for name, columns, model_kind, grid in model_specs:
		model, rows = _tune_predictor(name, columns, model_kind, grid,
			build, tune, y["build"], y["tune_model"], c)
		models[name] = model
		tuning.extend(rows)
	_log(c, "DONE", "global_models")
	_log(c, "START", "donor_crossfit", f"family_folds={c['folds']}")
	residual, oof, folds, crossfit_audit = _oof_residuals(build, y["build"], models["CSx"], c)
	_log(c, "DONE", "donor_crossfit", f"OOF_MSE={_mse(y['build'], oof):.6g}")
	common = {"ancestry": groups["ancestry"], "csx": groups["csx"]}
	arms = {"GRID_no_evolution": common.copy()}
	if groups["frequency"]:
		arms["GRID_frequency_only"] = {**common, "frequency": groups["frequency"]}
	if groups["evolution"]:
		arms["GRID_evolution"] = {**common, "frequency": groups["frequency"], "evolution": groups["evolution"]}
	if groups["permuted_evolution"]:
		arms["GRID_permuted_evolution"] = {**common, "frequency": groups["frequency"], "permuted_evolution": groups["permuted_evolution"]}
	if c["primary_arm"] not in arms:
		raise ValueError(f"Primary arm {c['primary_arm']} has no supplied features; explicit evolution features are required for GRID_evolution")
	banks, parameters, retrieval_rows, metric_rows = {}, {}, [], []
	tune_baseline = models["CSx"].predict(tune)
	for name, blocks in arms.items():
		_log(c, "START", "matching", name)
		bank = ReferenceBank().fit(build, residual, oof, folds, blocks, c)
		if c["abm_backend"] == "selective_attention":
			_log(c, "START", "selective_attention", f"{name} device={c['device']}")
			bank.attention = _attention().SelectiveAttention().fit(bank, build, tune, tune_baseline, y["tune_model"], blocks, c)
			bank.retrieval_audit.update(bank.attention.fit_status)
			tuning.extend(dict(stage="attention", model=name, **row,
				selected=row["epoch"] == bank.attention.fit_status["selected_epoch"]) for row in bank.attention.history)
		selected, rows = _tune_bank(name, bank, tune, tune_baseline, y["tune_model"], c)
		banks[name], parameters[name] = bank, selected
		tuning.extend(rows)
		retrieval_rows.append({"arm": name, **bank.retrieval_audit})
		metric_rows.append(bank.geometry.audit.assign(arm=name))
		_log(c, "DONE", "matching", f"{name} k={selected['k']} alpha={selected['alpha']} radius={selected['radius_multiplier']}")
	bundle = {
		"format": "GRID-ABM-1", "config": c, "feature_groups": groups, "all_features": all_features,
		"models": models, "banks": banks, "bank_parameters": parameters, "primary_arm": c["primary_arm"],
		"development_eids": _ids(development, "eid"), "development_family_ids": _ids(development, "family_id"),
		"donor_eids": _ids(build, "eid"), "audit": {"status": "not_evaluated"},
		"trait_type": c["trait_type"],
		"baseline_scope": "target_ancestry_specific_with_pooled_fallback",
		"disco_scope": "calibration of supplied reference-interpolation scalar; ancestry anchors are not phenotype-refitted here",
		"csx_scope": "target-cohort calibration of supplied posterior scores; upstream phi is not retuned here",
		"interpretation": "individual-reference OOF residual correction; explanations are associative, not causal",
	}
	_log(c, "START", "reliability_gate", f"tune_gate={len(frame['tune_gate'])}")
	gate_raw, _ = _raw_predictions(bundle, frame["tune_gate"])
	bundle["gate"] = ReliabilityGate().fit(_gate_table(gate_raw, c["primary_arm"]),
		y["tune_gate"], gate_raw[c["primary_arm"]].to_numpy(), c)
	calibration = frame["calibration"]
	cal_family = _ids(calibration, "family_id")
	cal_fold = _folds(cal_family, 2, c["seed"] + 251)
	threshold = calibration.iloc[np.flatnonzero(cal_fold == 0)].reset_index(drop=True)
	audit_frame = calibration.iloc[np.flatnonzero(cal_fold == 1)].reset_index(drop=True)
	if min(len(threshold), len(audit_frame)) < max(2, c["min_role_n"] // 2):
		raise ValueError("Insufficient independent calibration fit/audit participants")
	threshold_raw, _ = _raw_predictions(bundle, threshold)
	values = bundle["gate"].predict(_gate_table(threshold_raw, c["primary_arm"]))
	threshold_ids = _ids(threshold, "eid")
	bundle["coverage_rules"] = {
		q: CoverageRule().fit(values["gain"], threshold_ids, q, c["seed"] + 271)
		for q in c["coverages"]
	}
	bundle["absolute_error_rule"] = CoverageRule().fit(-values["absolute_error"], threshold_ids,
		c["coverage"], c["seed"] + 277)
	bundle["confidence_reference"] = np.sort(np.asarray(values["absolute_error"], dtype=float))
	control_scores = {
		"absolute_error": -values["absolute_error"],
		"support_only": threshold_raw["support"].to_numpy() * threshold_raw["family_ess"].to_numpy(),
		"random": np.array([_seed(v, c["seed"] + 313) / 2**64 for v in threshold_ids]),
	}
	if c["trait_type"] == "binary":
		control_scores["clinical_lowrisk"] = -threshold_raw["covariate_baseline"].to_numpy()
	bundle["control_rules"] = {
		name: {q: CoverageRule().fit(score, threshold_ids, q, c["seed"] + 281)
			for q in c["coverages"]} for name, score in control_scores.items()
	}
	audit_y = _labels(audit_frame, c["trait_type"], "calibration_audit", require_two_classes=False)
	bundle["audit"] = _internal_audit(bundle, audit_frame, audit_y)
	bundle["low_error_audit"] = _internal_low_error_audit(bundle, audit_frame, audit_y)
	bundle["confidence_audit"] = _internal_low_error_audit(bundle, audit_frame, audit_y, "confidence_candidate")
	bundle["ancestry_gain_audit"], bundle["ancestry_low_error_audit"] = {}, {}
	for ancestry in sorted(audit_frame.ancestry.astype(str).unique()):
		mask = audit_frame.ancestry.astype(str).eq(ancestry).to_numpy()
		subset = audit_frame.loc[mask].reset_index(drop=True)
		bundle["ancestry_gain_audit"][ancestry] = _internal_audit(bundle, subset, audit_y[mask])
		bundle["ancestry_low_error_audit"][ancestry] = _internal_low_error_audit(bundle, subset, audit_y[mask])
	_log(c, "DONE", "reliability_gate", f"independent_audit={bundle['audit']['status']}")
	# These stronger competitors use every label in the outer training half.
	# They are fitted AFTER the matching policy audit and never replace its
	# baseline, gate inputs, thresholds, or audit target. Their hyperparameters
	# are already selected from build/tune_model, so no test information is used.
	full_sources = ["CSx"] + (["DiscoDivas_calibrated_reference"] if "DiscoDivas_calibrated_reference" in models else [])
	if c["full_training_global_controls"]:
		full_sources += [name for name in ("Ridge_no_evolution", "HGB_no_evolution", "Ridge_evolution", "HGB_evolution") if name in models]
	_log(c, "START", "full_training_comparators", f"development={len(development)}")
	for source in full_sources:
		original = models[source]
		models[source + "_full_training"] = _make_predictor(original.columns, original.kind,
			original.model_kind, original.parameter, c).fit(development, y_development)
	_log(c, "DONE", "full_training_comparators")
	_log(c, "START", "frozen_test_prediction", f"n={len(test)}")
	raw, matches = _raw_predictions(bundle, test, retain_matches=True)
	predictions = _apply_policy(bundle, raw)
	predictions["split"] = "test"
	_log(c, "DONE", "frozen_test_prediction", f"selected={int(predictions.selected.sum())} released={int(predictions.released.sum())}")
	role_detail = development[["eid", "family_id", "ancestry"]].copy()
	role_detail["role"] = role
	cal_ix = part["calibration"]
	role_detail.loc[cal_ix[cal_fold == 0], "role"] = "calibration_fit"
	role_detail.loc[cal_ix[cal_fold == 1], "role"] = "calibration_audit"
	role_counts = role_detail.groupby("role", sort=False).agg(
		n=("eid", "size"), families=("family_id", "nunique")
	).reset_index()
	role_counts["fraction_of_development"] = role_counts["n"] / len(development)
	role_counts["fraction_of_outer_cohort"] = role_counts["n"] / len(data)
	bundle["training_diagnostics"] = {
		"abm_backend": c["abm_backend"], "attention_device": c["device"] if c["abm_backend"] == "selective_attention" else None,
		"retrieval_backend": c["retrieval"],
		"outer_train_n": len(development), "outer_test_n": len(test), "outer_train_fraction": train_fraction,
		"donor_n": len(build), "donor_fraction_of_outer_cohort": len(build) / len(data),
		"test_y_accessed": False, "primary_coverage": c["coverage"],
		"family_disjoint": True, "donor_residual_source": "family_out_of_fold_CSx_baseline",
		"baseline_scope": bundle["baseline_scope"],
		"gate_reference": "CSx fitted only in build, with alpha selected in tune_model",
		"full_training_comparator_changes_matching_policy": False,
	}
	budgets = []
	for name in models:
		full = name.endswith("_full_training")
		tuned = name != "covariate_baseline"
		budgets.append(dict(method=name, fit_n=len(development) if full else len(build),
			tune_n=len(tune) if tuned else 0, donor_n=0, gate_n=0, independent_audit_n=0,
			threshold_x_only_n=0, unique_development_labels_used=len(development) if full else len(build) + (len(tune) if tuned else 0),
			outer_development_n=len(development), scope="full_training_refit_fixed_hyperparameters" if full else "build_fit"))
	for name in banks:
		budgets.append(dict(method=name, fit_n=len(build), tune_n=len(tune), donor_n=len(build),
			gate_n=0, independent_audit_n=0, threshold_x_only_n=0,
			unique_development_labels_used=len(build) + len(tune), outer_development_n=len(development),
			scope="build_metric_and_nested_oof_donor_residuals"))
	budgets.append(dict(method="GRID_policy", fit_n=len(build), tune_n=len(tune), donor_n=len(build),
		gate_n=len(frame["tune_gate"]), independent_audit_n=len(audit_frame), threshold_x_only_n=len(threshold),
		unique_development_labels_used=len(build) + len(tune) + len(frame["tune_gate"]) + len(audit_frame),
		outer_development_n=len(development), scope="frozen_matching_policy_with_independent_internal_audit"))
	budget_table = pd.DataFrame(budgets)
	bundle["method_label_budget"] = budget_table
	baseline_groups = pd.concat([
		models["CSx"].group_audit(tune, y["tune_model"]).assign(model="CSx", fit_role="build"),
		models["CSx_full_training"].group_audit(tune, y["tune_model"]).assign(model="CSx_full_training", fit_role="entire_development"),
	], ignore_index=True).rename(columns={"build_n": "fit_n", "build_groups": "fit_groups",
		"build_cases": "fit_cases", "build_controls": "fit_controls"})
	return {
		"bundle": bundle, "predictions": predictions, "matches": matches,
		"roles": role_detail, "role_counts": role_counts, "tuning": pd.DataFrame(tuning),
		"crossfit_audit": pd.DataFrame(crossfit_audit),
		"donor_oof": pd.DataFrame({"eid": _ids(build, "eid"), "family_id": _ids(build, "family_id"),
			"fold": folds, "baseline_oof": oof, "residual_oof": residual}),
		"retrieval_audit": pd.DataFrame(retrieval_rows),
		"baseline_groups": baseline_groups, "method_label_budget": budget_table,
		"metric_features": pd.concat(metric_rows, ignore_index=True),
		"calibration_audit": pd.DataFrame([{**bundle["audit"], "ancestry":"ALL", "used_for_policy":True},
			*[{**row, "ancestry":name, "used_for_policy":False} for name,row in bundle["ancestry_gain_audit"].items()]]),
		"confidence_audit": pd.DataFrame([bundle["confidence_audit"]]),
		"low_error_audit": pd.DataFrame([{**bundle["low_error_audit"], "ancestry":"ALL"},
			*[{**row, "ancestry":name} for name,row in bundle["ancestry_low_error_audit"].items()]]),
		"training_diagnostics": bundle["training_diagnostics"],
	}


def predict_new(bundle: dict, data: pd.DataFrame, retain_matches: bool = True) -> dict:
	"""Frozen inference on identity/predictor columns; no y column is required."""
	if bundle.get("format") != "GRID-ABM-1":
		raise ValueError("Unrecognized GRID bundle format")
	ids = _ids(data, "eid")
	family = _ids(data, "family_id")
	_ids(data, "ancestry")
	if len(set(ids)) != len(ids):
		raise ValueError("Duplicate projection IDs")
	_validate_features(data, bundle["feature_groups"])
	if bundle["config"]["strict_new_ids"] and set(ids) & set(bundle["development_eids"]):
		raise ValueError("Projection includes development individuals; genuine new-person inference requires unseen IDs")
	columns = _unique(["eid", "family_id", "ancestry", *bundle["all_features"]])
	frame = data[columns].reset_index(drop=True)
	raw, matches = _raw_predictions(bundle, frame, retain_matches=retain_matches)
	predictions = _apply_policy(bundle, raw)
	predictions["development_family_overlap"] = np.isin(family, bundle["development_family_ids"])
	predictions["split"] = "projection"
	return {"bundle": bundle, "predictions": predictions, "matches": matches}


def materialize_matches(result: dict, arm: str | None = None, query_indices=None) -> pd.DataFrame:
	"""Expand actual matches for requested query rows; avoid whole-cohort expansion
	when only an individual's explanation is needed. Packed arrays retain all
	donors for every query in each retained arm.
	"""
	if arm is None:
		arm = result["bundle"]["primary_arm"]
	if arm not in result.get("matches", {}):
		raise ValueError(f"Packed matches were not retained for {arm}")
	m = result["matches"][arm]
	query = np.arange(len(m["query_eids"])) if query_indices is None else np.asarray(query_indices, int)
	ix = m["neighbor_indices"][query]
	rr, cc = np.nonzero(ix >= 0)
	qi, di = query[rr], ix[rr, cc]
	frame = pd.DataFrame({
		"eid": m["query_eids"][qi], "arm": arm, "rank": cc + 1,
		"donor_eid": m["donor_eids"][di], "donor_family_id": m["donor_family_ids"][di],
		"donor_ancestry": m["donor_ancestry"][di], "distance": m["distance"][qi, cc],
		"weight": m["weights"][qi, cc],
		"donor_oof_baseline": m["donor_oof_baseline"][di],
		"donor_observed_y": m["donor_observed_y"][di],
		"donor_oof_residual": m["donor_oof_residual"][di],
		"donor_crossfit_fold": m["donor_fold"][di],
		"borrow_fraction": m["borrow_fraction"][qi],
		"weighted_residual_contribution": m["weighted_residual_contribution"][qi, cc],
	})
	for name, distances in m["distance_components"].items():
		frame[f"distance_{name}"] = distances[qi, cc]
	return frame
