#!/usr/bin/env python3
"""abm utilities and command entry points."""
from __future__ import annotations
import sys

sys.dont_write_bytecode = True
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

# 🚩 Persisted model types
# These schema identifiers keep frozen content hashes stable across file moves.
# No saved model, prediction, or receipt is rewritten during loading.
MODEL_TYPES = {
	"abm_twin_core": [
		"Config",
		"CensoringKM",
		"NumericPrep",
		"ModelProfileGeometry",
		"FullProfileGeometry",
		"Calibrator",
		"TwinBank",
		"FrozenRule",
		"TwinABM",
	],
	"abm_predictability_final": [
		"UnweightedKNNBank",
		"ClusterMixtureExperts",
		"AlternativeRiskRoutes",
		"PredictabilityExtension",
	],
	"molecular_deviation": ["DeviationConfig", "CovPrep", "MolecularNormalizer", "DeviationStudy"],
}


# Persisted C1 schema names survive source consolidation. These are class-name
# mappings, not source checksums; loading never rewrites a frozen artifact.
MODEL_TYPES.update({
    "c1_unified_models": ["IndividualAssessment"],
    "c1_prevalent": ["PrevalentConfig", "PrevalentPrototypeTransformer"],
    "c1_abm_gpu": ["MolecularTransformer", "GpuFrozenTwinBank"],
    "c1_reliability_eval": ["ReliabilityHead"],
    "c1_registry": ["MethodSpec"],
})
MODEL_MODULES = {module: ("deviation" if module == "molecular_deviation" else "abm")
                 for module in MODEL_TYPES}
MODEL_MODULES.update({"c1_unified_models": "c1_models", "c1_prevalent": "c1_models",
                      "c1_abm_gpu": "c1_deep", "c1_reliability_eval": "c1_evaluation",
                      "c1_registry": "c1_protocol"})


def frozen_type_name(value):
	cls = type(value)
	for module, names in MODEL_TYPES.items():
		if cls.__name__ in names and cls.__module__ in {module, MODEL_MODULES[module]}:
			return module + "." + cls.__qualname__
	return cls.__module__ + "." + cls.__qualname__


def load_model(path):
	"""Read trusted joblib results, resolving the original model schema names."""
	from joblib.numpy_pickle import NumpyUnpickler, _validate_fileobject_and_memmap
	import importlib

	class ModelUnpickler(NumpyUnpickler):
		def find_class(self, module, name):
			if name in MODEL_TYPES.get(module, ()):
				destination = MODEL_MODULES[module]
				return getattr(importlib.import_module(destination), name)
			return super().find_class(module, name)

	with open(path, "rb") as stream:
		with _validate_fileobject_and_memmap(stream, str(path), None) as (handle, _):
			if isinstance(handle, str):
				raise ValueError("Unsupported joblib format; expected a current LE8 model")
			return ModelUnpickler(str(path), handle, ensure_native_byte_order=True).load()


from le8_common import save_figure

# 🚩 abm twin core
# Research implementation of training-donor pruning + outcome-blind twin selection.
#
# This is a new, independently runnable numeric core, NOT a replacement for all of
# 0.host.py. The host keeps phenotype definitions, assay mapping, family rosters,
# and its existing output exporter. No test labels are accepted by predict().
#
# Estimand: fixed-horizon net risk under independent censoring. Death is treated as
# censoring when represented as event=0. This is NOT a competing-risk estimator.
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from typing import Any
import warnings

import joblib
import numpy as np
from scipy.special import expit, logit
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.model_selection import StratifiedGroupKFold

VERSION = "ABM-TWIN-2026-10-17-c1-lockpaper"
ROLES = ("build", "tune_model", "tune_gate", "calibration_fit", "calibration_audit")


@dataclass(frozen=True)
class Config:
	trait: str = "cvd_cad"
	layer: str = "prot"
	seed: int = 12345
	horizon: float = 10.0
	folds: int = 5
	repeats: int = 2
	keep_fraction: float = 0.95
	pruning: str = "profile_concordance"  # literal model--outcome agreement; class/risk control retained
	# Opt-in reference experiment. Legacy defaults preserve the original protocol.
	reference_radius_policy: str = "arm_calibrated"
	twin_random_control: bool = False
	min_retained_cases: int = 0
	min_case_retention: float = 0.0
	reference_experiment: str = "legacy"
	# Legacy defaults are deliberate: old frozen models retain their original behavior.
	# New public method presets select case-preserving policies explicitly.
	geometry: str = "coefficient_panel"
	panel_selection: str = "pilot"
	technical_qc_scope: str = "legacy_panel"
	panel_size: int = 10
	always_include_indices: tuple[int, ...] = ()
	profile_kind: str = "molecular_predictors"
	risk_bins: int = 5
	soft_floor: float = 0.1
	primary: str = "elasticnet_clean"
	selector: str = "support_quality"
	primary_coverage: float = 0.5
	coverages: tuple[float, ...] = (0.2, 0.4, 0.5, 0.6, 0.8, 1.0)
	c_grid: tuple[float, ...] = (0.01, 0.1, 1.0)
	pilot_c: float = 0.1
	max_iter: int = 2000
	max_iter_retries: int = 2
	lightgbm: bool = False
	tree_iterations: int = 150
	cores: int = 4
	parallel_workers: int = 4
	feature_missing: float = 0.3
	sample_missing: float = 0.2
	dimensions: int = 24
	k: int = 64
	# Zero reproduces historic k-dependent radii. New methods freeze K_ref=64
	# independently of the maximum K used to estimate a probability.
	radius_k_ref: int = 0
	min_twins: int = 20
	min_family_ess: float = 15.0
	radius_quantile: float = 0.9
	prior_strength: float = 20.0
	gate_iterations: int = 120
	min_censor_survival: float = 0.05
	min_events: int = 20
	audit_min_cases: int = 30
	audit_min_controls: int = 30
	bootstrap: int = 1000
	grouping_unit: str = "family"
	final_multiroute_release: str = ""

	def validate(self) -> None:
		from c1_protocol import RELEASE, COVERAGES, validate_seed
		if self.final_multiroute_release == RELEASE:
			validate_seed(self.seed)
			if not np.isclose(self.primary_coverage, 0.5, rtol=0, atol=1e-12):
				raise ValueError("C1 LOCKPAPER fixes primary_coverage at 0.50")
		if self.final_multiroute_release == RELEASE and (tuple(self.coverages) != COVERAGES or self.bootstrap != 1000):
			raise ValueError("C1 LOCKPAPER fixes six coverage points and 1000 family/group bootstrap replicates")
		if os.getenv("ABM_FINAL_SELECTOR", "risk_resid_error").strip() != "risk_resid_error":
			raise ValueError("C1 LOCKPAPER fixes primary selector risk_resid_error")
		if self.grouping_unit not in {"family", "individual_unverified"}:
			raise ValueError("Invalid grouping unit")
		if not (0 < self.keep_fraction <= 1 and self.horizon > 0):
			raise ValueError("Invalid retention fraction or horizon")
		if self.pruning not in {"profile_concordance", "raw_loss", "class_risk_stratified", "case_preserving", "none"}:
			raise ValueError("Unknown pruning mode")
		if self.geometry not in {"coefficient_panel", "all_molecular_pca", "all_molecular"}:
			raise ValueError("Unknown molecular geometry")
		if self.panel_selection not in {"pilot", "all"}:
			raise ValueError("Unknown panel selection policy")
		if self.technical_qc_scope not in {"legacy_panel", "assay"}:
			raise ValueError("Unknown technical QC scope")
		if self.reference_radius_policy not in {"arm_calibrated", "full_fixed"}:
			raise ValueError("Unknown reference radius policy")
		if not isinstance(self.twin_random_control, bool):
			raise ValueError("twin_random_control must be boolean")
		if self.min_retained_cases < 0 or not 0 <= self.min_case_retention <= 1:
			raise ValueError("Invalid predeclared case-retention requirement")
		if self.selector not in {"support_quality", "absolute_error", "rank_error", "support_only"}:
			raise ValueError("Unknown outcome-blind selector")
		if self.primary_coverage not in self.coverages or any(not 0 < q <= 1 for q in self.coverages):
			raise ValueError("Primary coverage must be predeclared in coverages")
		if self.panel_size < 2 or any(j < 0 for j in self.always_include_indices):
			raise ValueError("Panel size >=2 and valid fixed covariate indices required")
		if self.folds < 2 or self.repeats < 1 or self.k < self.min_twins:
			raise ValueError("Invalid cross-fitting/matching parameters")
		if type(self.radius_k_ref) is not int or self.radius_k_ref < 0:
			raise ValueError("radius_k_ref must be a nonnegative integer (0 means legacy k)")
		if not (0 < self.radius_quantile < 1 and 0 < self.soft_floor <= 1):
			raise ValueError("Invalid radius/soft-weight setting")
		if (
			self.risk_bins < 1
			or self.dimensions < 1
			or type(self.parallel_workers) is not int
			or not 1 <= self.parallel_workers <= 4
			or self.cores < 1
			or self.gate_iterations < 1
			or self.min_family_ess <= 0
			or self.prior_strength < 0
			or self.min_events < 2
			or not 0 < self.feature_missing <= 1
			or not 0 <= self.sample_missing < 1
			or not 0 < self.min_censor_survival <= 1
			or any(c <= 0 for c in self.c_grid)
		):
			raise ValueError("Invalid fitting/quality thresholds")
		if self.audit_min_cases < 2 or self.audit_min_controls < 2 or self.bootstrap < 1:
			raise ValueError("Invalid audit information requirements")
		if (
			isinstance(self.max_iter, bool)
			or not isinstance(self.max_iter, int)
			or self.max_iter < 1
			or isinstance(self.max_iter_retries, bool)
			or not isinstance(self.max_iter_retries, int)
			or not 0 <= self.max_iter_retries <= 5
		):
			raise ValueError(
				"max_iter must be a positive integer; max_iter_retries must be an integer from 0 to 5"
			)


def tie_hash(value: str, seed: int) -> int:
	return int.from_bytes(hashlib.sha256(f"{seed}:{value}".encode()).digest()[:8], "big")


def fingerprint(ids: np.ndarray) -> str:
	return hashlib.sha256("\n".join(map(str, ids)).encode()).hexdigest()


def model_hash(value: Any) -> str:
	"""Canonical content hash, invariant to pickle's pointer/alias bookkeeping."""

	def digest(x):
		h = hashlib.sha256()

		def put(v):
			h.update(v if isinstance(v, bytes) else str(v).encode())

		if x is None:
			put("none")
		elif isinstance(x, np.generic):
			return digest(x.item())
		elif isinstance(x, np.ndarray):
			put("array:")
			put(x.dtype.str)
			put(repr(x.shape))
			if x.dtype.hasobject:
				for item in x.ravel(order="C"):
					put(digest(item))
			else:
				put(np.ascontiguousarray(x).tobytes())
		elif isinstance(x, (str, bytes, int, float, bool)):
			put(type(x).__name__)
			put(repr(x))
		elif isinstance(x, dict):
			put("dict")
			for k in sorted(x, key=lambda v: (type(v).__name__, repr(v))):
				put(digest(k))
				put(digest(x[k]))
		elif isinstance(x, (list, tuple)):
			put(type(x).__name__)
			for item in x:
				put(digest(item))
		elif isinstance(x, (set, frozenset)):
			put(type(x).__name__)
			for item in sorted(digest(v) for v in x):
				put(item)
		elif type(x).__module__.startswith("lightgbm") and hasattr(x, "model_to_string"):
			put("lightgbm_booster")
			put(x.model_to_string())
		elif hasattr(x, "__dict__"):
			put(frozen_type_name(x))
			put(digest(vars(x)))
		elif hasattr(x, "__getstate__"):
			put(frozen_type_name(x))
			put(digest(x.__getstate__()))
		else:
			raise TypeError("Unsupported fitted-state hash type: " + repr(type(x)))
		return h.digest()

	return digest(value).hex()


def json_write(path: Path, value: Any) -> None:
	def clean(x):
		if isinstance(x, Path):
			return str(x)
		if isinstance(x, dict):
			return {str(k): clean(v) for k, v in x.items()}
		if isinstance(x, (tuple, list)):
			return [clean(v) for v in x]
		if isinstance(x, np.ndarray):
			return clean(x.tolist())
		if isinstance(x, np.generic):
			return clean(x.item())
		if isinstance(x, float) and not np.isfinite(x):
			return None
		return x

	path.write_text(json.dumps(clean(value), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


class CensoringKM:
	"""Training-only reverse KM; fail explicitly when censoring positivity fails."""

	def fit(self, time, event, horizon: float, minimum: float):
		t, e = np.asarray(time, float), np.asarray(event, int)
		self.horizon, self.minimum = horizon, minimum
		self.times, inverse, totals = np.unique(t, return_inverse=True, return_counts=True)
		censored = np.bincount(inverse, weights=(e == 0), minlength=len(totals))
		risk = np.cumsum(totals[::-1])[::-1]
		self.survival = np.cumprod(1 - censored / risk)
		return self

	def at(self, time, left=False):
		pos = np.searchsorted(self.times, time, side="left" if left else "right") - 1
		return np.where(pos < 0, 1.0, self.survival[np.maximum(pos, 0)])

	def labels_weights(self, time, event):
		t, e = np.asarray(time, float), np.asarray(event, int)
		case = (e == 1) & (t <= self.horizon)
		control = t > self.horizon
		# A censoring exactly at horizon is conservatively not treated as a control.
		y = case.astype(int)
		w = np.zeros(len(t))
		gc = self.at(t[case], left=True)
		gh = float(self.at(np.array([self.horizon]))[0])
		if (len(gc) and np.min(gc) < self.minimum) or (control.any() and gh < self.minimum):
			raise ValueError("Insufficient training censoring support at the chosen horizon")
		w[case] = 1 / gc
		w[control] = 1 / gh
		return y, w


class NumericPrep:
	"""Medians, retained columns and scale learned exclusively in this fit's data."""

	def fit(self, X, config: Config):
		x = np.asarray(X, float)
		self.input_columns = x.shape[1]
		self.clinical_columns = tuple(config.always_include_indices)
		self.sample_missing = config.sample_missing
		self.technical_qc_scope = config.technical_qc_scope
		if not set(self.clinical_columns) <= set(range(self.input_columns)):
			raise ValueError("Mandatory clinical column index outside input schema")
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", RuntimeWarning)
			variance = np.nanvar(x, axis=0)
		self.keep = np.flatnonzero(
			(np.mean(np.isnan(x), axis=0) < config.feature_missing) & (variance > 1e-12)
		)
		if len(self.keep) < 2:
			raise ValueError("Fewer than two usable numeric features")
		# Explicit schema identities, never name heuristics or clinical dummy counts.
		self.molecular_qc_columns = self.molecular_columns.copy()
		good = self.technical_qc(x)
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", RuntimeWarning)
			variance = np.nanvar(x[good][:, self.keep], axis=0)
		self.keep = self.keep[np.isfinite(variance) & (variance > 1e-12)]
		if len(self.keep) < 2 or not len(self.molecular_columns):
			raise ValueError("No usable molecular panel after sample QC")
		self.fit_statistics(x[self.technical_qc(x)])
		return self

	@property
	def molecular_columns(self):
		return np.asarray([j for j in self.keep if j not in self.clinical_columns], int)

	def fit_statistics(self, X):
		"""Freeze imputation/scale using only technically eligible training rows."""
		z = np.asarray(X, float)[:, self.keep]
		if not len(z):
			raise ValueError("No molecular-QC eligible training samples")
		self.median = np.nanmedian(z, axis=0)
		if not np.isfinite(self.median).all():
			raise ValueError("Unmeasured predictor after molecular QC")
		z = np.where(np.isfinite(z), z, self.median)
		molecular = np.isin(self.keep, self.molecular_columns)
		if not np.any(np.var(z[:, molecular], axis=0) > 1e-12):
			raise ValueError("No variable molecular predictors after panel sample QC")
		self.mean, self.scale = z.mean(0), np.maximum(z.std(0), 1e-8)

	def transform(self, X):
		x = np.asarray(X, float)
		if x.ndim != 2 or x.shape[1] != self.input_columns or np.isinf(x).any():
			raise ValueError("Feature schema mismatch or infinite input")
		z = x[:, self.keep]
		return (np.where(np.isfinite(z), z, self.median) - self.mean) / self.scale

	def missing(self, X, *, panel=True):
		columns = self.molecular_columns if panel else self.molecular_qc_columns
		if not len(columns):
			raise ValueError("No usable molecular columns for QC")
		return np.mean(~np.isfinite(np.asarray(X)[:, columns]), axis=1)

	def technical_qc(self, X):
		# New method comparisons use one broad assay rule independent of panel size.
		if getattr(self, "technical_qc_scope", "legacy_panel") == "assay":
			return self.missing(X, panel=False) <= self.sample_missing
		# Legacy fits additionally checked panel-specific missingness.
		return (self.missing(X, panel=False) <= self.sample_missing) & (
			self.missing(X) <= self.sample_missing
		)


def check_classes(y, w, minimum=2):
	if min(np.sum((w > 0) & (y == k)) for k in (0, 1)) < minimum:
		raise ValueError("Insufficient known cases/controls; do not treat censored rows as controls")


def fit_classifier(X, y, w, config: Config, kind="elasticnet", parameter=0.1):
	ok = np.asarray(w) > 0
	check_classes(y, w, config.min_events)
	if kind == "elasticnet":
		model = LogisticRegression(
			C=float(parameter),
			solver="saga",
			penalty="elasticnet",
			l1_ratio=0.5,
			max_iter=config.max_iter,
			random_state=config.seed,
		)
	elif kind == "lightgbm":
		from lightgbm import LGBMClassifier

		model = LGBMClassifier(
			n_estimators=config.tree_iterations,
			num_leaves=int(parameter),
			learning_rate=0.03,
			min_child_samples=40,
			reg_lambda=10,
			n_jobs=config.cores,
			random_state=config.seed,
			verbosity=-1,
			deterministic=True,
			force_col_wise=True,
		)
	else:
		raise ValueError(kind)
	x_fit, y_fit = np.asarray(X)[ok], np.asarray(y)[ok]
	weight = np.asarray(w)[ok] / np.mean(np.asarray(w)[ok])
	if kind == "elasticnet":
		# Let fit finish so coefficients are available to warm-start a retry.
		# Escalating the warning to an exception inside sklearn loses that state.
		attempts = []
		for attempt in range(config.max_iter_retries + 1):
			budget = config.max_iter * 2**attempt
			model.set_params(max_iter=budget, warm_start=attempt > 0)
			with warnings.catch_warnings(record=True) as caught:
				warnings.simplefilter("always", ConvergenceWarning)
				model.fit(x_fit, y_fit, sample_weight=weight)
			convergence = [v for v in caught if issubclass(v.category, ConvergenceWarning)]
			for warning in caught:
				if not issubclass(warning.category, ConvergenceWarning):
					warnings.warn_explicit(
						warning.message, warning.category, warning.filename, warning.lineno
					)
			attempts.append(
				dict(max_iter=budget, n_iter=int(np.max(model.n_iter_)), converged=not convergence)
			)
			total = sum(v["n_iter"] for v in attempts)
			if not convergence:
				model.convergence_audit_ = attempts
				if attempt:
					print(
						f"[LE8] DONE ABM elastic-net converged: C={parameter:g}; "
						f"attempts={len(attempts)}; total_iterations={total}",
						flush=True,
					)
				return model
			if attempt < config.max_iter_retries:
				print(
					f"[LE8] PROGRESS ABM elastic-net convergence retry {attempt+1}/{config.max_iter_retries}: "
					f"C={parameter:g}; rows={len(x_fit)}; features={x_fit.shape[1]}; "
					f"total_iterations={total}; next_max_iter={budget*2}; solver=saga; tol={model.tol:g}",
					flush=True,
				)
		raise RuntimeError(
			f"ABM elastic-net did not converge: solver=saga, C={parameter:g}, tol={model.tol:g}, "
			f"iterations={total}, attempts={len(attempts)}. "
			"Review training scales or increase max_iter/max_iter_retries in --abm-config; "
			"no unconverged fit was accepted."
		)
	with warnings.catch_warnings():
		warnings.simplefilter("error", ConvergenceWarning)
		model.fit(x_fit, y_fit, sample_weight=weight)
	return model


def logloss_rows(y, p):
	p = np.clip(p, 1e-7, 1 - 1e-7)
	return -(y * np.log(p) + (1 - y) * np.log1p(-p))


def choose_model_panel(prep, X, y, w, config):
	"""Selection uses only this training role. Never import full-cohort ProtWAS hits.

    Pilot elastic-net coefficients rank omics within an a-priori assay budget.
    A constant fitted pilot is a failure, not permission to select arbitrary hits.
    Clinical columns named by explicit indices are retained outside the omic budget.
    """
	if not set(config.always_include_indices) <= set(range(np.asarray(X).shape[1])):
		raise ValueError("Mandatory clinical column index outside input schema")
	if config.panel_selection == "all":
		# No global coefficient or univariate-significance screen precedes local learning.
		# Assay missingness/variance QC is still frozen from this training role.
		prep.panel_source = "all QC-passed predictors in this fit training role"
		prep.panel_pilot_coefficients = np.full(len(prep.keep), np.nan)
		prep.panel_pilot_zero_count = 0
		prep.panel_pilot_convergence = []
		prep.fit_statistics(np.asarray(X)[prep.technical_qc(X)])
		return prep.keep.copy()
	fit = fit_classifier(prep.transform(X), y, w * prep.technical_qc(X), config, parameter=config.pilot_c)
	importance = np.abs(fit.coef_[0])
	mandatory = np.array([j for j, k in enumerate(prep.keep) if k in config.always_include_indices], int)
	candidate = np.array([j for j, k in enumerate(prep.keep) if k not in config.always_include_indices], int)
	if not len(candidate):
		raise ValueError("No usable molecular predictors for panel selection")
	if len(candidate) and not np.any(importance[candidate] > 1e-10):
		raise ValueError(
			"Panel pilot has no nonzero molecular coefficients; change the development-only specification"
		)
	order = candidate[np.lexsort((prep.keep[candidate], -importance[candidate]))]
	chosen = np.sort(np.r_[mandatory, order[: config.panel_size]])
	if len(chosen) < 2:
		raise ValueError("Too few selected predictors")
	original = prep.keep[chosen].copy()
	for name in ("keep", "median", "mean", "scale"):
		setattr(prep, name, getattr(prep, name)[chosen])
	prep.panel_source = "this fit's training data only"
	prep.panel_pilot_coefficients = fit.coef_[0][chosen].copy()
	prep.panel_pilot_zero_count = int(np.sum(np.abs(prep.panel_pilot_coefficients) < 1e-10))
	prep.panel_pilot_convergence = fit.convergence_audit_
	prep.fit_statistics(np.asarray(X)[prep.technical_qc(X)])
	return original


class ModelProfileGeometry:
	"""Coefficient-weighted profile geometry with an explicit predictor subset."""

	def fit(self, X, coefficients, columns=None):
		X = np.asarray(X)
		beta = np.abs(np.asarray(coefficients, float))
		self.columns = np.arange(X.shape[1], dtype=int) if columns is None else np.asarray(columns, int)
		if len(self.columns) < 2:
			raise ValueError("Twin geometry requires at least two predictors")
		if np.any(self.columns < 0) or np.any(self.columns >= X.shape[1]):
			raise ValueError("Invalid geometry columns")
		b = beta[self.columns]
		self.weights = np.sqrt(0.25 + 0.75 * b / max(b.mean(), 1e-8))
		self.n_components_ = len(self.columns)
		return self

	def transform(self, X):
		return np.asarray(X)[:, self.columns] * self.weights



class FullProfileGeometry:
	"""Outcome-blind complete-profile geometry; learned on eligible build X only.

	PCA is an optional dimension reduction, not discovery of causal disease paths.
	The full-space option retains all QC-passed features with equal initial scale.
	"""

	def fit(self, X, columns, config: Config):
		x = np.asarray(X, float)
		self.columns = np.asarray(columns, int).copy()
		if len(self.columns) < 2 or len(np.unique(self.columns)) != len(self.columns):
			raise ValueError("Complete profile geometry requires at least two distinct columns")
		if np.any(self.columns < 0) or np.any(self.columns >= x.shape[1]):
			raise ValueError("Invalid complete profile geometry columns")
		self.input_columns = x.shape[1]
		self.kind = config.geometry
		self.pca = None
		if self.kind == "all_molecular_pca":
			self.pca = PCA(n_components=min(config.dimensions, len(self.columns), len(x) - 1),
				svd_solver="randomized", random_state=config.seed).fit(x[:, self.columns])
			self.n_components_ = self.pca.n_components_
		elif self.kind == "all_molecular":
			self.n_components_ = len(self.columns)
		else:
			raise ValueError("FullProfileGeometry requires a complete-profile configuration")
		return self

	def transform(self, X):
		x = np.asarray(X, float)
		if x.ndim != 2 or x.shape[1] != self.input_columns:
			raise ValueError("Complete profile geometry input schema changed")
		z = x[:, self.columns]
		return self.pca.transform(z) if self.pca is not None else z

def _parallel_fits(function, tasks, config, matrix_bytes, label):
	"""Ordered process results, read-only shared arrays, and no nested thread pools.

	The phase cgroup is the hard process-tree cap. This estimate limits copies
	before dispatch; no independent worker receives the full phase budget.
	"""
	from joblib import Parallel, delayed, parallel_config
	from threadpoolctl import threadpool_limits
	import tempfile
	limit = min(32.0, float(os.getenv("LE8_PHASE_MEMORY_GIB", "32")))
	# Parent working arrays, Python/scientific imports, and worker slicing/QC copies.
	parent = 2.0 + matrix_bytes * 3 / 2**30
	per_worker = 1.0 + matrix_bytes * 8 / 2**30
	workers = min(config.parallel_workers, config.cores, len(tasks),
	              max(1, int((limit - parent - 2.0) / per_worker)))
	print(f"[LE8] ABM {label}: {len(tasks)} tasks, {workers} processes; "
	      f"one thread/worker; phase budget={limit:g} GiB", flush=True)
	with tempfile.TemporaryDirectory(prefix="le8-abm-workers-", dir="/tmp") as tmp:
		with threadpool_limits(limits=1), parallel_config(backend="loky", inner_max_num_threads=1):
			return Parallel(n_jobs=workers, max_nbytes="1M", mmap_mode="r",
			                temp_folder=tmp, pre_dispatch=workers)(
				delayed(function)(*args) for args in tasks)


def _quality_fold(X, time, event, ids, groups, config, repeat, fold, tr, va):
	import time as _clock
	started = _clock.monotonic()
	print(f"[LE8] ABM fold START repeat={repeat+1} fold={fold+1} pid={os.getpid()}", flush=True)
	logs = []
	if set(groups[tr]) & set(groups[va]):
		raise AssertionError("Family leakage in OOF")
	prep = NumericPrep().fit(X[tr], config)
	train_qc = prep.technical_qc(X[tr])
	km = CensoringKM().fit(
		time[tr][train_qc], event[tr][train_qc], config.horizon, config.min_censor_survival
	)
	y, w = km.labels_weights(time[tr], event[tr])
	panel = choose_model_panel(prep, X[tr], y, w, config)
	train_qc = prep.technical_qc(X[tr])
	km = CensoringKM().fit(
		time[tr][train_qc], event[tr][train_qc], config.horizon, config.min_censor_survival
	)
	y, w = km.labels_weights(time[tr], event[tr])
	model = fit_classifier(prep.transform(X[tr]), y, w * train_qc, config, parameter=config.pilot_c)
	p = model.predict_proba(prep.transform(X[va]))[:, 1]
	yy, ww = km.labels_weights(time[va], event[va])
	technical = prep.technical_qc(X[va])
	prob = p
	losses = np.where((ww > 0) & technical, logloss_rows(yy, p), np.nan)
	logs.append(
		dict(
			repeat=repeat,
			fold=fold,
			train_n=len(tr),
			recipient_n=len(va),
			train_technical_QC_excluded=int((~train_qc).sum()),
			recipient_technical_QC_excluded=int((~technical).sum()),
			molecular_qc_input_indices=prep.molecular_qc_columns.tolist(),
			molecular_panel_input_indices=prep.molecular_columns.tolist(),
			recipient_hash=fingerprint(ids[va]),
			family_disjoint=True,
			panel_input_indices=panel.tolist(),
			pilot_iterations=sum(v["n_iter"] for v in prep.panel_pilot_convergence),
			pilot_fit_attempts=len(prep.panel_pilot_convergence),
			model_iterations=sum(v["n_iter"] for v in model.convergence_audit_),
			model_fit_attempts=len(model.convergence_audit_),
			selection_in_teacher_only=True,
		)
	)
	print(f"[LE8] ABM fold DONE repeat={repeat+1} fold={fold+1} seconds={_clock.monotonic()-started:.1f}", flush=True)
	return repeat, va, prob, losses, technical, logs[0]


class InsufficientReferenceSupport(ValueError):
	"""Donor pruning left a model arm without its required training classes."""

	def __init__(self, audit, composition):
		self.audit = audit
		self.composition = composition
		failures = [f"{row['arm']}: cases={row['known_cases']}, controls={row['known_controls']} "
			f"(minimum={row['minimum_per_class']})" for row in audit if not row["sufficient"]]
		super().__init__("Insufficient known cases/controls after reference pruning: " + "; ".join(failures) +
			". Censored rows remain unknown; no replacement model was fitted.")


def check_reference_support(y, w, multipliers, config, composition):
	"""Check all arms before dispatching any expensive final fits."""
	audit = []
	for arm, mult in multipliers.items():
		known = (w * mult) > 0
		cases = int(np.sum(known & (y == 1)))
		controls = int(np.sum(known & (y == 0)))
		audit.append(dict(arm=arm, known_cases=cases, known_controls=controls,
			minimum_per_class=config.min_events, sufficient=min(cases, controls) >= config.min_events))
	if any(not row["sufficient"] for row in audit):
		raise InsufficientReferenceSupport(audit, composition)
	return audit


def _fit_arm(z, y, w, b, t, a, kind, arm, mult):
	tuning = []
	best = None
	for par in a.c_grid if kind == "elasticnet" else (7, 15, 31):
		model = fit_classifier(z[b], y[b], w[b] * mult, a, kind, par)
		p = model.predict_proba(z[t])[:, 1]
		score = float(np.average(logloss_rows(y[t], p), weights=w[t]))
		fit_audit = getattr(model, "convergence_audit_", [])
		tuning.append(
			dict(
				model=kind + "_" + arm,
				parameter=par,
				logloss=score,
				iterations=sum(v["n_iter"] for v in fit_audit) if fit_audit else None,
				fit_attempts=len(fit_audit) if fit_audit else None,
			)
		)
		if best is None or score < best[0]:
			best = score, model
	return kind + "_" + arm, best[1], tuning


def crossfit_quality(X, time, event, ids, groups, config: Config):
	"""Own outcomes are NEVER in the teacher fit. They ARE used to score donors.

    This use of build labels is deliberate supervised data selection. It is not
    evidence that discarded labels are erroneous. Fold preprocessing and KM are
    refitted; hyperparameters are fixed pilot settings, not tuned on recipient labels.
    """
	n = len(ids)
	prob = np.full((config.repeats, n), np.nan)
	losses = np.full_like(prob, np.nan)
	technical = np.zeros((config.repeats, n), bool)
	logs = []
	strat = ((np.asarray(event) == 1) & (np.asarray(time) <= config.horizon)).astype(int)
	tasks = []
	for repeat in range(config.repeats):
		splitter = StratifiedGroupKFold(config.folds, shuffle=True, random_state=config.seed + repeat)
		for fold, (tr, va) in enumerate(splitter.split(X, strat, groups)):
			tasks.append((X, time, event, ids, groups, config, repeat, fold, tr, va))
	for repeat, va, p, loss, qc, audit in _parallel_fits(_quality_fold, tasks, config, X.nbytes, "crossfit"):
		prob[repeat, va], losses[repeat, va], technical[repeat, va] = p, loss, qc
		logs.append(audit)

	if not np.isfinite(prob).all():
		raise AssertionError("Incomplete OOF predictions")
	with warnings.catch_warnings():
		warnings.simplefilter("ignore", RuntimeWarning)
		quality = np.nanmedian(losses, axis=0)
	quality[~technical.all(axis=0)] = np.nan
	return dict(
		oof_probability=prob.mean(0),
		oof_sd=prob.std(0),
		oof_loss=quality,
		model_concordance=np.exp(-quality),
		repeat_losses=losses,
		technical_qc=technical.all(axis=0),
		audit=logs,
	)


def prune_donors(y, w, quality, ids, config: Config, incident_event=None):
	"""Return training-reference masks/weights; never filter evaluation people.

	``case_preserving`` protects *all observed incident events*, including events
	after the target horizon. Their horizon labels remain correct (late events are
	horizon controls); protection never relabels them as horizon cases. Only known
	event-free controls can be removed by loss. Technical assay QC and unknown
	fixed-horizon labels remain distinct from supervised reference selection.
	"""
	y, w = np.asarray(y), np.asarray(w)
	protective = config.pruning in {"case_preserving", "none"}
	if incident_event is None:
		if protective:
			raise ValueError("Case-preserving pruning requires all-follow-up incident_event labels")
		event = y.copy()
	else:
		event = np.asarray(incident_event)
		if event.shape != y.shape or not np.isin(event, [0, 1]).all():
			raise ValueError("Invalid incident_event labels for reference selection")
	final_qc = np.asarray(quality.get("reference_technical_qc", quality["technical_qc"]), bool)
	finite_quality = np.isfinite(quality["oof_loss"])
	valid = ((w > 0) & final_qc if protective else
		(w > 0) & quality["technical_qc"] & finite_quality)
	if not valid.any():
		raise ValueError("No QC-eligible donors with known outcomes")
	p, err = np.asarray(quality["oof_probability"]), np.asarray(quality["oof_loss"])
	if not np.isfinite(p[valid]).all():
		raise ValueError("Reference probabilities must be finite; cannot silently discard unscored cases")
	edges = np.unique(np.quantile(p[valid], np.linspace(0, 1, config.risk_bins + 1)))[1:-1]
	bins = np.searchsorted(edges, p, side="right")
	strata = (y * (config.risk_bins + 1) + bins
		if config.pruning == "class_risk_stratified" else np.zeros(len(y), int))
	keep, soft, control = np.zeros(len(y), bool), np.zeros(len(y)), np.zeros(len(y), bool)
	ties = np.asarray([tie_hash(s, config.seed) for s in ids], np.uint64)
	audit = []
	protected = valid & (event == 1) if protective else np.zeros(len(y), bool)
	# Missing OOF quality is uncertainty about the teacher, never proof of a bad row.
	unscored = valid & ~finite_quality if protective else np.zeros(len(y), bool)
	if config.pruning == "none":
		keep[valid], soft[valid] = True, 1.0
		audit.append(dict(stratum=0, total=int(valid.sum()), retained=int(valid.sum()), policy="no_loss_filter"))
	else:
		fixed = protected | unscored
		keep[fixed], soft[fixed] = True, 1.0
		for label in np.unique(strata[valid]):
			idx = np.flatnonzero(valid & ~fixed & (strata == label))
			if not len(idx):
				continue
			order = idx[np.lexsort((ties[idx], err[idx]))]
			count = max(1, math.ceil(config.keep_fraction * len(idx)))
			keep[order[:count]] = True
			rank = (np.arange(len(order)) + 0.5) / len(order)
			soft[order] = config.soft_floor + (1 - config.soft_floor) * (1 - rank)
			audit.append(dict(stratum=int(label), total=len(idx), retained=count))
		if protective:
			audit.append(dict(stratum=-1, total=int(fixed.sum()), retained=int(fixed.sum()),
				policy="protect_all_incident_events_and_unscored_rows"))
	# Same-count, same case x risk-bin random reference control. Protected late
	# events have horizon y=0 and must remain protected in this arm as well.
	control[protected | unscored] = True
	for label in np.unique((y * (config.risk_bins + 1) + bins)[valid]):
		idx = np.flatnonzero(valid & ((y * (config.risk_bins + 1) + bins) == label))
		pool = idx[~control[idx]]
		r = sorted(pool, key=lambda j: tie_hash(ids[j], config.seed + 701))
		n = int(keep[idx].sum() - control[idx].sum())
		if n < 0:
			raise AssertionError("Random reference violates protected donor count")
		control[np.asarray(r[:n], int)] = True
	case_n = int(np.sum(valid & (y == 1)))
	retained_case_n = int(np.sum(keep & (y == 1)))
	case_retention = retained_case_n / case_n if case_n else np.nan
	composition = dict(eligible_N=int(valid.sum()), retained_N=int(keep.sum()),
		eligible_cases=case_n, retained_cases=retained_case_n,
		case_retention=case_retention, removed_cases=case_n - retained_case_n,
		removed_controls=int(np.sum(valid & ~keep & (y == 0))),
		eligible_all_followup_events=int(np.sum(valid & (event == 1))),
		retained_all_followup_events=int(np.sum(keep & (event == 1))),
		protected_all_followup_events=int(protected.sum()),
		retained_unscored_N=int(np.sum(keep & ~finite_quality)),
		pruning_mode=config.pruning,
		retention_denominator="known-horizon technically eligible build donors",
		interpretation="Supervised reference selection; high loss does not establish noise or an erroneous outcome")
	if protective and (np.any(protected & ~keep) or np.any(protected & ~control)
		or np.any(soft[protected] != 1)):
		raise AssertionError("Case-preserving policy erased or downweighted an incident event")
	if retained_case_n < config.min_retained_cases or case_retention < config.min_case_retention:
		raise ValueError("Reference pruning violates the predeclared case-retention requirement: "
			f"{retained_case_n}/{case_n} cases retained. Use a separately named case-preserving "
			"reference experiment or lower pruning strength; do not relabel removed cases as noise.")
	if case_n and case_retention < 0.5:
		warnings.warn(f"Reference pruning retains only {retained_case_n}/{case_n} cases; "
			"report this composition change alongside local prediction performance.", RuntimeWarning)
	# Nonfinite teacher quality in a protected donor is represented conservatively
	# for diagnostics only. It does not change outcomes or probability weights.
	finite_values = err[valid & finite_quality]
	unscored_loss = max(10.0, float(np.max(finite_values))) if len(finite_values) else 10.0
	donor_loss = np.where(finite_quality, err, unscored_loss)
	return dict(keep=keep, soft=soft, random_keep=control, eligible=valid, strata=strata, audit=audit,
		composition=composition, donor_loss=donor_loss, protected=protected)


class Calibrator:
	def fit(self, p, y, w):
		check_classes(y, w)
		ok = w > 0
		z = logit(np.clip(p, 1e-7, 1 - 1e-7)).reshape(-1, 1)
		self.model = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
		self.model.fit(z[ok], y[ok], sample_weight=w[ok] / w[ok].mean())
		return self

	def predict(self, p):
		return self.model.predict_proba(logit(np.clip(p, 1e-7, 1 - 1e-7)).reshape(-1, 1))[:, 1]


class TwinBank:
	"""All-radius candidates counted exactly; only <=k are used for borrowing.

    The supplied frozen model-panel geometry defines proximity, not query outcomes. Radius is learned
    from independent tune_model queries; family and self matches are excluded.
    """

	def fit(self, X, y, w, ids, groups, losses, config: Config, geometry=None):
		self.config = config
		self.pca = geometry or PCA(
			n_components=min(config.dimensions, X.shape[1], len(X) - 1),
			svd_solver="randomized",
			random_state=config.seed,
		).fit(X)
		self.Z = self.pca.transform(X) / np.sqrt(self.pca.n_components_)
		self.y, self.w = np.asarray(y), np.asarray(w)
		self.ids, self.groups, self.losses = ids.copy(), groups.copy(), losses.copy()
		from sklearn.neighbors import KDTree

		self.index = KDTree(self.Z, leaf_size=40) if self.Z.shape[1] <= 24 else None
		self.radius = None
		return self

	def calibrate_radius(self, X, ids, groups):
		z = self.pca.transform(X) / np.sqrt(self.pca.n_components_)
		k_ref = getattr(self.config, "radius_k_ref", 0) or self.config.k
		kth = []
		for start in range(0, len(z), 128):
			query = z[start : start + 128]
			if self.index is not None:
				distance, near = self.index.query(query, k=min(k_ref + 64, len(self.Z)))
				for j, (ds, nn) in enumerate(zip(distance, near)):
					ok = (self.groups[nn] != groups[start + j]) & (self.ids[nn] != ids[start + j])
					good = ds[ok]
					if len(good) >= k_ref:
						kth.append(good[k_ref - 1])
						continue
					# Unusually large families: exact full-radius fallback for this query.
					dd = euclidean_distances(query[j : j + 1], self.Z)[0]
					dd[(self.groups == groups[start + j]) | (self.ids == ids[start + j])] = np.inf
					kth.append(
						np.partition(dd, min(k_ref, len(dd)) - 1)[min(k_ref, len(dd)) - 1]
					)
			else:
				d = euclidean_distances(query, self.Z)
				valid = (groups[start : start + 128, None] != self.groups[None, :]) & (
					ids[start : start + 128, None] != self.ids[None, :]
				)
				d[~valid] = np.inf
				kk = min(k_ref, d.shape[1]) - 1
				kth.extend(np.partition(d, kk, axis=1)[:, kk])
		finite = np.asarray(kth)[np.isfinite(kth)]
		if len(finite) < 10:
			raise ValueError("Not enough independent queries to calibrate matching radius")
		self.radius = float(np.quantile(finite, self.config.radius_quantile))
		if self.radius <= 1e-12:
			raise ValueError("Collapsed matching geometry")
		return self

	def predict(self, X, ids, groups, fallback):
		if self.radius is None:
			raise RuntimeError("Radius is not frozen")
		z = self.pca.transform(X) / np.sqrt(self.pca.n_components_)
		n = len(z)
		out = {
			k: np.zeros(n)
			for k in (
				"eligible_match_count",
				"matched_count",
				"eligible_reference_count",
				"twin_fraction",
				"matching_radius",
				"family_count",
				"family_ess",
				"donor_ess",
				"donor_loss",
				"radius_ratio",
				"nearest_radius_ratio",
				"label_disagreement",
				"borrow_fraction",
			)
		}
		out["donor_loss"].fill(10.0)
		out["radius_ratio"].fill(100.0)
		out["nearest_radius_ratio"].fill(100.0)
		out["matching_radius"].fill(self.radius)
		# Denominator excludes self and the query's family, exactly as matching does.
		# Compute lookup maps once, including when reading a legacy frozen bank.
		valid_reference = np.asarray(self.w) > 0
		ref_groups, ref_n = np.unique(np.asarray(self.groups)[valid_reference], return_counts=True)
		group_size = dict(zip(ref_groups.tolist(), ref_n.tolist()))
		id_groups = dict(zip(np.asarray(self.ids)[valid_reference].tolist(), np.asarray(self.groups)[valid_reference].tolist()))
		reference_n = int(valid_reference.sum())
		for i, (query_id, query_group) in enumerate(zip(ids, groups)):
			denom = reference_n - group_size.get(query_group, 0)
			if query_id in id_groups and id_groups[query_id] != query_group:
				denom -= 1
			out["eligible_reference_count"][i] = denom
		prob = np.asarray(fallback, float).copy()
		for start in range(0, n, 128):
			query = z[start : start + 128]
			if self.index is not None:
				near, distances = self.index.query_radius(query, r=self.radius, return_distance=True)
			else:
				dense = euclidean_distances(query, self.Z)
				near = [np.flatnonzero(dd <= self.radius) for dd in dense]
				distances = [dd[ii] for dd, ii in zip(dense, near)]
			for row, (nn, ds) in enumerate(zip(near, distances)):
				i = start + row
				allowed = (self.groups[nn] != groups[i]) & (self.ids[nn] != ids[i]) & (self.w[nn] > 0)
				idx = nn[allowed]
				ds = ds[allowed]
				out["eligible_match_count"][i] = len(idx)
				denom = out["eligible_reference_count"][i]
				out["twin_fraction"][i] = len(idx) / denom if denom > 0 else np.nan
				if not len(idx):
					continue
				# Stable tie order is independent of donor outcomes.
				order = np.lexsort((self.ids[idx], ds))[: self.config.k]
				idx = idx[order]
				ds = ds[order]
				weight = np.exp(-0.5 * (ds / self.radius) ** 2) * self.w[idx]
				weight /= weight.sum()
				_, inv = np.unique(self.groups[idx], return_inverse=True)
				fw = np.bincount(inv, weights=weight)
				ess, fess = 1 / np.sum(weight**2), 1 / np.sum(fw**2)
				local = float(weight @ self.y[idx])
				borrowing = fess / (fess + self.config.prior_strength)
				prob[i] = borrowing * local + (1 - borrowing) * fallback[i]
				for name, value in dict(
					matched_count=len(idx),
					family_count=len(fw),
					family_ess=fess,
					donor_ess=ess,
					donor_loss=weight @ self.losses[idx],
					radius_ratio=ds[-1] / self.radius,
					nearest_radius_ratio=ds[0] / self.radius,
					label_disagreement=local * (1 - local),
					borrow_fraction=borrowing,
				).items():
					out[name][i] = value
		out["supported"] = (out["matched_count"] >= self.config.min_twins) & (
			out["family_ess"] >= self.config.min_family_ess
		)
		return np.clip(prob, 1e-7, 1 - 1e-7), out


def pairwise_error_target(y, w, probability):
	"""Per-person, weighted pair-discordance on independent gate-training data.

    This is an AUC-oriented surrogate, not a guarantee about a jointly selected
    subpopulation's AUC. O(n log n), no n-by-n pair matrix.
    """
	p, y, w = np.asarray(probability), np.asarray(y), np.asarray(w)
	target = np.full(len(p), 0.5)
	check_classes(y, w)
	for label in (0, 1):
		other = (y != label) & (w > 0)
		order = np.argsort(p[other], kind="stable")
		po, wo = p[other][order], w[other][order]
		cumulative = np.r_[0.0, np.cumsum(wo)]
		i = np.flatnonzero((y == label) & (w > 0))
		left = cumulative[np.searchsorted(po, p[i], side="left")]
		right = cumulative[np.searchsorted(po, p[i], side="right")]
		total = cumulative[-1]
		target[i] = (
			((total - right) + 0.5 * (right - left)) / total
			if label == 1
			else (left + 0.5 * (right - left)) / total
		)
	return target


class FrozenRule:
	def fit(self, score, eligible, ids, coverage, seed):
		self.seed, self.coverage = seed, coverage
		idx = np.flatnonzero(eligible & np.isfinite(score))
		self.support_limited = len(idx) < math.ceil(coverage * len(ids))
		# Coverage is relative to ALL calibration queries, never only supported ones.
		count = min(len(idx), math.ceil(coverage * len(ids)))
		self.threshold, self.tie = np.inf, 0
		if count:
			ties = np.array([tie_hash(v, seed) for v in ids], np.uint64)
			order = idx[np.lexsort((ties[idx], -np.asarray(score)[idx]))]
			boundary = order[count - 1]
			self.threshold, self.tie = float(score[boundary]), int(ties[boundary])
			if count == len(idx):
				self.threshold = -np.inf
		return self

	def apply(self, score, eligible, ids):
		ties = np.array([tie_hash(v, self.seed) for v in ids], np.uint64)
		score = np.asarray(score)
		return (
			eligible
			& np.isfinite(score)
			& ((score > self.threshold) | ((score == self.threshold) & (ties <= self.tie)))
		)


class TwinABM:
	def __init__(self, config: Config = Config()):
		config.validate()
		self.config = config

	@staticmethod
	def validate_features(X, ids, groups):
		x = np.asarray(X, float)
		ids, groups = np.asarray(ids).astype(str), np.asarray(groups).astype(str)
		if x.ndim != 2 or len(x) != len(ids) or len(ids) != len(groups) or np.isinf(x).any():
			raise ValueError("Invalid numeric feature input")
		invalid = {"", "nan", "none", "na", "<na>"}
		if (
			len(set(ids)) != len(ids)
			or any(v.strip().lower() in invalid for v in ids)
			or any(v.strip().lower() in invalid for v in groups)
		):
			raise ValueError("IDs must be unique and all family groups explicitly defined")
		return x, ids, groups

	def fit(self, X, time, event, ids, groups, roles, feature_names):
		a = self.config
		X, ids, groups = self.validate_features(X, ids, groups)
		time, event, roles = np.asarray(time, float), np.asarray(event), np.asarray(roles).astype(str)
		if (
			len(time) != len(X)
			or len(event) != len(X)
			or len(roles) != len(X)
			or not np.isfinite(time).all()
			or (time <= 0).any()
			or not np.isin(event, [0, 1]).all()
		):
			raise ValueError("Invalid survival outcomes")
		if set(roles) != set(ROLES):
			raise ValueError("Development must contain exactly the five declared roles, never test")
		group_roles = {}
		for g, r in zip(groups, roles):
			if g in group_roles and group_roles[g] != r:
				raise ValueError("Family crosses development roles")
			group_roles[g] = r
		self.feature_names = list(map(str, feature_names))
		if len(self.feature_names) != X.shape[1] or len(set(self.feature_names)) != len(self.feature_names):
			raise ValueError("Invalid named numeric schema")
		self.development_ids, self.development_groups = ids.copy(), groups.copy()
		self.audit_ids = ids[roles == "calibration_audit"].copy()
		self.audit_groups = groups[roles == "calibration_audit"].copy()
		parts = {r: np.flatnonzero(roles == r) for r in ROLES}
		b, t, gate, cal = (parts[r] for r in ROLES[:4])
		self.prep = NumericPrep().fit(X[b], a)
		if a.geometry != "coefficient_panel":
			# Keep full-X distance independent of subsequent global-model feature ranking.
			import copy
			self.match_prep = copy.deepcopy(self.prep)
		initial_qc = self.prep.technical_qc(X[b])
		self.km = CensoringKM().fit(
			time[b][initial_qc], event[b][initial_qc], a.horizon, a.min_censor_survival
		)
		# Audit outcomes are not even passed to label/weight computation during fit.
		non_audit = roles != "calibration_audit"
		y, w = np.zeros(len(X), int), np.zeros(len(X))
		y[non_audit], w[non_audit] = self.km.labels_weights(time[non_audit], event[non_audit])
		self.panel_indices = choose_model_panel(self.prep, X[b], y[b], w[b], a)
		self.panel_names = [self.feature_names[j] for j in self.panel_indices]
		technical = self.prep.technical_qc(X)
		self.molecular_qc_version = VERSION
		self.km = CensoringKM().fit(
			time[b][technical[b]], event[b][technical[b]], a.horizon, a.min_censor_survival
		)
		y[non_audit], w[non_audit] = self.km.labels_weights(time[non_audit], event[non_audit])
		horizon_known = w > 0
		w *= technical
		for r in ROLES[:4]:
			check_classes(y[parts[r]], w[parts[r]], a.min_events)
		z = self.prep.transform(X)
		quality = crossfit_quality(X[b], time[b], event[b], ids[b], groups[b], a)
		oof_qc = quality["technical_qc"].copy()
		quality["technical_qc"] &= technical[b]
		quality["reference_technical_qc"] = technical[b].copy()
		pruning = prune_donors(y[b], w[b], quality, ids[b], a, incident_event=event[b])
		reference_qc = technical[b] if a.pruning in {"case_preserving", "none"} else quality["technical_qc"]
		multipliers = {
			"full": np.ones(len(b)),
			"clean": pruning["keep"].astype(float),
			"soft": pruning["soft"],
			"random": pruning["random_keep"].astype(float),
		}
		check_reference_support(y[b], w[b], multipliers, a, pruning["composition"])
		self.models, self.model_tuning = {}, []
		tasks = [(z, y, w, b, t, replace(a, cores=1), kind, arm, mult)
		         for kind in (["elasticnet", "lightgbm"] if a.lightgbm else ["elasticnet"])
		         for arm, mult in multipliers.items()]
		for name, model, tuning in _parallel_fits(_fit_arm, tasks, a, z.nbytes, "model arms"):
			self.models[name] = model
			self.model_tuning.extend(tuning)

		# C1_LOCKPAPER::2026-10-17.c1-lockpaper-v1
		# Risk models keep clinical+omics. Primary twin distance is molecular-only.
		eligible = pruning["eligible"]
		self.banks, self.banks_integrated = {}, {}
		if a.geometry == "coefficient_panel":
			match_z = z
			mol_positions = np.flatnonzero(~np.isin(self.prep.keep, np.asarray(self.prep.clinical_columns, int)))
			geometry_names = self.panel_names
			coef = self.models["elasticnet_full"].coef_[0]
			mol_geometry = ModelProfileGeometry().fit(match_z[b], coef, columns=mol_positions)
			integrated_geometry = ModelProfileGeometry().fit(match_z[b], coef, columns=np.arange(match_z.shape[1]))
		else:
			match_z = self.match_prep.transform(X)
			mol_positions = np.flatnonzero(~np.isin(self.match_prep.keep, np.asarray(self.match_prep.clinical_columns, int)))
			geometry_names = [self.feature_names[j] for j in self.match_prep.keep]
			# Imputation, scaling, feature QC and PCA use the build role only.
			fit_rows = b[self.match_prep.technical_qc(X[b])]
			mol_geometry = FullProfileGeometry().fit(match_z[fit_rows], mol_positions, a)
			integrated_geometry = FullProfileGeometry().fit(match_z[fit_rows], np.arange(match_z.shape[1]), a)
		if len(mol_positions) < 2:
			raise ValueError("Fewer than two molecular predictors in twin geometry")
		self.geometry_molecular_model_positions = mol_positions.copy()
		self.geometry_molecular_features = [geometry_names[j] for j in mol_positions]
		reference_arms = [("full", eligible), ("clean", pruning["keep"])]
		if a.twin_random_control:
			reference_arms.append(("random", pruning["random_keep"]))
		for arm, select in reference_arms:
			bb, tt = b[select], t[technical[t]]
			for banks, geometry in ((self.banks, mol_geometry), (self.banks_integrated, integrated_geometry)):
				bank = TwinBank().fit(match_z[bb], y[bb], w[bb], ids[bb], groups[bb], pruning["donor_loss"][select], a, geometry)
				if arm != "full" and a.reference_radius_policy == "full_fixed":
					# Freeze the full-bank radius for a reference-composition ablation.
					# Re-estimating it per arm would change both reference set and neighbourhood.
					bank.radius = banks["full"].radius
				else:
					bank.calibrate_radius(match_z[tt], ids[tt], groups[tt])
				banks[arm] = bank
		self.calibrators = {}
		pred, _ = self._raw(X[cal], ids[cal], groups[cal])
		self.calibrators = {name: Calibrator().fit(p, y[cal], w[cal]) for name, p in pred.items()}
		if a.primary not in self.calibrators:
			raise ValueError("Requested primary model was not fitted")
		# Gate is trained on independent gate-role outcomes, never build residual labels.
		gp, gd = self._calibrated(X[gate], ids[gate], groups[gate])
		gx = self._gate_features(gp, gd)
		self.gate_models = {}
		targets = {
			"absolute_error": (y[gate] - gp[a.primary]) ** 2,
			"rank_error": pairwise_error_target(y[gate], w[gate], gp[a.primary]),
		}
		known = w[gate] > 0
		for name, target in targets.items():
			model = HistGradientBoostingRegressor(
				max_iter=a.gate_iterations,
				max_leaf_nodes=5,
				learning_rate=0.04,
				min_samples_leaf=max(20, min(100, int(known.sum() / 10))),
				l2_regularization=10,
				early_stopping=False,
				random_state=a.seed + 903,
			)
			model.fit(gx[known], target[known], sample_weight=w[gate][known] / w[gate][known].mean())
			self.gate_models[name] = model
		cp, cd = self._calibrated(X[cal], ids[cal], groups[cal])
		score, _ = self._scores(cp, cd, ids[cal])
		self.rules = {}
		calibration_qc = cd["technical_qc"]
		for name, values in score.items():
			allowed = (
				calibration_qc
				if name in ("confidence", "lowrisk", "random")
				else calibration_qc & cd["supported"]
			)
			for q in a.coverages:
				self.rules[(name, q)] = FrozenRule().fit(values, allowed, ids[cal], q, a.seed + 91)
		self.risk_edges = np.unique(
			np.quantile(cp["elasticnet_full"][calibration_qc], np.linspace(0, 1, 11))
		)[1:-1]
		self.quality_edges = np.unique(np.quantile(score[a.selector][calibration_qc], np.linspace(0, 1, 6)))[
			1:-1
		]
		self.training_audit = dict(
			oof=quality["audit"],
			pruning=pruning["audit"],
			pruning_composition=pruning["composition"],
			reference_experiment=a.reference_experiment,
			reference_radius_policy=a.reference_radius_policy,
			twin_random_control=a.twin_random_control,
			reference_bank_radii={"molecular": {k: v.radius for k, v in self.banks.items()},
				"integrated": {k: v.radius for k, v in self.banks_integrated.items()}},
			known_donors=int(eligible.sum()),
			retained_donors=int(pruning["keep"].sum()),
			build_n=len(b),
			molecular_qc_version=VERSION,
			technical_QC_excluded=int((~reference_qc).sum()),
			oof_teacher_QC_excluded=int((~quality["technical_qc"]).sum()),
			profile_discordance_excluded=int((eligible & ~pruning["keep"]).sum()),
			outcome_unavailable_N=int((~horizon_known[b]).sum()),
			molecular_qc_features=[self.feature_names[j] for j in self.prep.molecular_qc_columns],
			molecular_panel_features=[self.feature_names[j] for j in self.prep.molecular_columns],
			sample_missing_threshold=a.sample_missing,
			selected_panel=self.panel_names,
			panel_budget=a.panel_size,
			selected_zero_pilot_coefficients=self.prep.panel_pilot_zero_count,
			panel_pilot_iterations=sum(v["n_iter"] for v in self.prep.panel_pilot_convergence),
			panel_pilot_fit_attempts=len(self.prep.panel_pilot_convergence),
			final_nonzero_coefficients={
				n: int(np.sum(np.abs(m.coef_[0]) > 1e-10))
				for n, m in self.models.items()
				if hasattr(m, "coef_")
			},
			matching_space=("selected molecular predictors; coefficient-weighted distances"
				if a.geometry == "coefficient_panel" else a.geometry + "; build-only complete molecular profile"),
			geometry=a.geometry, panel_selection=a.panel_selection,
			technical_qc_scope=a.technical_qc_scope,
			matching_space_sensitivity=("selected molecular + mandatory baseline clinical predictors"
				if a.geometry == "coefficient_panel" else "complete molecular + baseline clinical profile"),
			geometry_molecular_features=self.geometry_molecular_features,
			clinical_predictors_excluded_from_distance=[self.feature_names[j] for j in self.prep.clinical_columns],
			geometry_integrated_features=geometry_names,
			primary_predictability_selector="risk_resid_error",
			release=a.final_multiroute_release,
			seed=a.seed,
			primary_coverage=a.primary_coverage,
			grouping_unit=a.grouping_unit,
			independent_family_support_verified=a.grouping_unit == "family",
			coverage_grid=list(a.coverages),
			selector_bootstrap_requested=a.bootstrap,
			confidence_definition="100*geometric_mean(proximity,independent_support,donor_concordance,profile_radius); zero without support",
			retained_cases=int(np.sum(y[b][pruning["keep"]])),
			known_cases=int(np.sum(y[b])),
			incident_events_build=int(np.sum(event[b] == 1)),
			incident_events_after_technical_qc=int(np.sum((event[b] == 1) & technical[b])),
			retained_incident_events=int(np.sum((event[b] == 1) & pruning["keep"])),
			case_preservation_policy=a.pruning in {"case_preserving", "none"},
			literal_pruning_is_supervised_selection_not_proof_of_noise=True,
			roles={
				r: dict(
					N=len(parts[r]),
					cases=int(y[parts[r]].sum()) if r != "calibration_audit" else None,
					technical_QC_excluded=int((~technical[parts[r]]).sum()),
					id_hash=fingerprint(ids[parts[r]]),
				)
				for r in ROLES
			},
		)
		self.donor_table = dict(
			eid=ids[b],
			group=groups[b],
			oof_probability=quality["oof_probability"],
			oof_sd=quality["oof_sd"],
			oof_loss=quality["oof_loss"],
			model_concordance=quality["model_concordance"],
			retained=pruning["keep"],
			soft_weight=pruning["soft"],
			random_retained=pruning["random_keep"],
			incident_event=event[b],
			protected_incident_event=pruning["protected"],
			horizon_known=horizon_known[b],
			technical_qc=reference_qc,
			oof_technical_qc=oof_qc,
			panel_technical_qc=technical[b],
			molecular_missing_fraction=self.prep.missing(X[b]),
			donor_eligible=eligible,
			technical_QC_excluded=~reference_qc,
			profile_discordance_excluded=eligible & ~pruning["keep"],
		)
		self.frozen = True
		# calibration_audit labels have not been used for fitting or threshold choice.
		return self

	def _raw(self, X, ids, groups):
		if getattr(self, "molecular_qc_version", None) != VERSION:
			raise RuntimeError(
				"ABM molecular QC requires an explicit refit in a new run; preserve the old snapshot"
			)
		qc = self.prep.technical_qc(X)
		z = self.prep.transform(X)
		match_z = self.match_prep.transform(X) if hasattr(self, "match_prep") else z
		pred = {name: model.predict_proba(z)[:, 1] for name, model in self.models.items()}
		diag = None
		for arm, bank in self.banks.items():
			prob, d = bank.predict(match_z[qc], ids[qc], groups[qc], pred["elasticnet_full"][qc])
			pred["twin_" + arm] = pred["elasticnet_full"].copy()
			pred["twin_" + arm][qc] = prob
			pred["twin_molecular_" + arm] = pred["twin_" + arm].copy()
			if arm == "clean":
				diag = {name: np.zeros(len(X), dtype=value.dtype) for name, value in d.items()}
				for name, value in d.items():
					diag[name][qc] = value
				for name in ("radius_ratio", "nearest_radius_ratio"):
					diag[name][~qc] = 100.0
				diag["donor_loss"][~qc] = 10.0
		for arm, bank in self.banks_integrated.items():
			prob, _ = bank.predict(match_z[qc], ids[qc], groups[qc], pred["elasticnet_full"][qc])
			pred["twin_integrated_" + arm] = pred["elasticnet_full"].copy()
			pred["twin_integrated_" + arm][qc] = prob
		diag["technical_qc"] = qc
		diag["molecular_missing_fraction"] = self.prep.missing(X)
		diag["frozen_molecular_missing_fraction"] = self.prep.missing(X, panel=False)
		return pred, diag

	def _calibrated(self, X, ids, groups):
		pred, diag = self._raw(X, ids, groups)
		return {name: self.calibrators[name].predict(p) for name, p in pred.items()}, diag

	def _gate_features(self, p, d):
		return np.c_[
			p[self.config.primary],
			p["elasticnet_full"],
			p["twin_clean"],
			np.log1p(d["eligible_match_count"]),
			np.log1p(d["family_ess"]),
			d["radius_ratio"],
			d["donor_loss"],
			d["label_disagreement"],
			d["borrow_fraction"],
		]

	def _scores(self, p, d, ids):
		gx = self._gate_features(p, d)
		est = {n: m.predict(gx) for n, m in self.gate_models.items()}
		proximity = np.exp(-np.minimum(d["nearest_radius_ratio"], 100) ** 2)
		support = np.clip(d["family_ess"] / self.config.k, 0, 1)
		concordance = np.exp(-np.maximum(d["donor_loss"], 0))
		profile_radius = np.exp(-np.minimum(d["radius_ratio"], 100) ** 2)
		confidence = np.where(
			d["supported"] & d["technical_qc"],
			(proximity * support * concordance * profile_radius) ** 0.25,
			0.0,
		)
		d.update(
			proximity_score=proximity,
			independent_support_score=support,
			donor_concordance_score=concordance,
			profile_radius_score=profile_radius,
			pre_prediction_confidence=100 * confidence,
		)
		score = {
			"absolute_error": -est["absolute_error"],
			"rank_error": -est["rank_error"],
			"support_quality": confidence,
			"support_only": -d["radius_ratio"] + 0.01 * np.log1p(d["family_ess"]),
			"confidence": np.abs(p["elasticnet_full"] - 0.5),
			"lowrisk": -p["elasticnet_full"],
			"random": np.array([tie_hash(v, self.config.seed + 99) / 2**64 for v in ids]),
		}
		return score, est

	def predict(self, X, ids, groups, feature_names, *, allow_development_audit=False):
		"""Outcome-blind API. There is intentionally no y/time/event argument."""
		X, ids, groups = self.validate_features(X, ids, groups)
		if list(map(str, feature_names)) != self.feature_names:
			raise ValueError("Feature names/order mismatch")
		if allow_development_audit and (
			not set(ids) <= set(self.audit_ids) or not set(groups) <= set(self.audit_groups)
		):
			raise ValueError("Audit override is restricted to the predeclared audit roster")
		if not allow_development_audit and (
			set(ids) & set(self.development_ids) or set(groups) & set(self.development_groups)
		):
			raise ValueError("Projection overlaps development participants/families")
		p, d = self._calibrated(X, ids, groups)
		scores, est = self._scores(p, d, ids)
		masks = {}
		for (name, q), rule in self.rules.items():
			eligible = (
				d["technical_qc"]
				if name in ("confidence", "lowrisk", "random")
				else d["technical_qc"] & d["supported"]
			)
			masks[f"{name}:{q:g}"] = rule.apply(scores[name], eligible, ids)
		selected = masks[f"{self.config.selector}:{self.config.primary_coverage:g}"]
		p["fallback_policy"] = np.where(
			d["technical_qc"], np.where(selected, p[self.config.primary], p["elasticnet_full"]), np.nan
		)
		qc_reason = np.where(
			d["molecular_missing_fraction"] == 1,
			"all_panel_molecules_missing",
			"molecular_missingness_exceeds_threshold",
		)
		reason = np.where(
			~d["technical_qc"],
			qc_reason,
			np.where(
				~d["supported"],
				"no_adequate_twins",
				np.where(selected, "selected", "below_frozen_quality_cutoff"),
			),
		)
		return dict(
			eid=ids,
			group=groups,
			probability=p,
			diagnostics=d,
			estimated_error=est,
			masks=masks,
			selected=selected,
			reason=reason,
			selector_scores=scores,
			frozen_model_hash=model_hash(self),
			version=VERSION,
			confidence_semantics="reference support, not probability of correct outcome",
			evaluation_scope="calibration_audit" if allow_development_audit else "heldout_projection",
			profile_kind=self.config.profile_kind,
			panel_features=self.panel_names,
		)


import os
from sklearn.neighbors import KDTree
from sklearn.cluster import KMeans

EPS = 1e-12
ROUTE_NAMES = (
	"risk_resid_error",
	"logloss_resid",
	"local_stability",
	"model_consensus",
	"hybrid_predictability",
)


def _weighted_bin_residual(target, risk, weight, bins=10):
	"""Residualize row loss within frozen risk strata.

    This is the central safeguard against a gate that merely learns
    'low disease risk = easy/predictable'.
    """
	target = np.asarray(target, float)
	risk = np.asarray(risk, float)
	weight = np.asarray(weight, float)
	ok = np.isfinite(target) & np.isfinite(risk) & np.isfinite(weight) & (weight > 0)
	if ok.sum() < 50:
		return np.full(len(target), np.nan), np.array([-np.inf, np.inf]), np.array([np.nan])
	q = np.unique(np.quantile(risk[ok], np.linspace(0, 1, int(bins) + 1)))
	if len(q) < 3:
		mu = np.average(target[ok], weights=weight[ok])
		out = target - mu
		return out, np.array([-np.inf, np.inf]), np.array([mu])
	q[0], q[-1] = -np.inf, np.inf
	label = np.searchsorted(q[1:-1], risk, side="right")
	means = []
	global_mu = np.average(target[ok], weights=weight[ok])
	for j in range(len(q) - 1):
		use = ok & (label == j)
		means.append(np.average(target[use], weights=weight[use]) if use.any() else global_mu)
	means = np.asarray(means, float)
	return target - means[label], q, means


def _pct(values, reference):
	"""Frozen empirical percentile, higher is better."""
	values = np.asarray(values, float)
	reference = np.asarray(reference, float)
	reference = np.sort(reference[np.isfinite(reference)])
	ans = np.full(len(values), np.nan)
	if not len(reference):
		return ans
	ok = np.isfinite(values)
	ans[ok] = np.searchsorted(reference, values[ok], side="right") / len(reference)
	return ans


def _consensus(probability):
	names = [
		n
		for n in (
			"elasticnet_full",
			"elasticnet_clean",
			"elasticnet_soft",
			"twin_clean",
			"twin_full",
			"lightgbm_full",
			"lightgbm_clean",
			"knn_unweighted_full",
			"knn_unweighted_clean",
			"mixture_experts",
		)
		if n in probability
	]
	if len(names) < 2:
		n = len(next(iter(probability.values())))
		return np.zeros(n), names
	x = np.vstack([np.asarray(probability[n], float) for n in names]).T
	return -np.nanstd(x, axis=1), names


def _local_stability(d):
	n = len(d["technical_qc"])
	donor_loss = np.asarray(d.get("donor_loss", np.full(n, 10.0)), float)
	disagree = np.asarray(d.get("label_disagreement", np.ones(n)), float)
	radius = np.asarray(d.get("radius_ratio", np.full(n, 100.0)), float)
	nearest = np.asarray(d.get("nearest_radius_ratio", radius), float)
	fam = np.asarray(d.get("family_ess", np.zeros(n)), float)
	match = np.asarray(d.get("matched_count", np.zeros(n)), float)
	# No outcome from the query enters this formula. Donor loss was derived from
	# donor OOF performance in development and is frozen before query projection.
	return (
		-(donor_loss + disagree + np.log1p(np.maximum(radius, 0)) + 0.5 * np.log1p(np.maximum(nearest, 0)))
		+ 0.15 * np.log1p(np.maximum(fam, 0))
		+ 0.05 * np.log1p(np.maximum(match, 0))
	)


class UnweightedKNNBank:
	"""Alternative molecular geometry independent of TwinBank coefficient weighting."""

	def fit(self, z, y, w, keep, dimensions, seed):
		keep = np.asarray(keep, bool) & (np.asarray(w) > 0)
		if keep.sum() < 100:
			raise ValueError("Too few KNN donors")
		d = min(int(dimensions), z.shape[1], keep.sum() - 1)
		self.pca = PCA(n_components=max(1, d), svd_solver="randomized", random_state=seed).fit(z[keep])
		self.Z = self.pca.transform(z[keep]) / np.sqrt(self.pca.n_components_)
		self.y = np.asarray(y)[keep].astype(float)
		self.w = np.asarray(w)[keep].astype(float)
		self.tree = KDTree(self.Z, leaf_size=40)
		return self

	def predict(self, z, k):
		q = self.pca.transform(z) / np.sqrt(self.pca.n_components_)
		k = max(3, min(int(k), len(self.Z)))
		dist, ix = self.tree.query(q, k=k)
		# Each query uses its own kth-neighbour radius. Other query rows must
		# never change a frozen participant prediction or its selector mask.
		ww = np.exp(-0.5 * np.square(dist / (dist[:, -1:] + 1e-6))) * self.w[ix]
		den = ww.sum(1)
		den = np.where(den > 0, den, 1.0)
		return np.sum(ww * self.y[ix], axis=1) / den


class ClusterMixtureExperts:
	"""Unsupervised profile clusters with cluster-specific disease models.

    This directly tests the 'different roads to the same disease' hypothesis.
    Cluster definition is X-only; outcome labels only fit the experts in build.
    """

	def fit(self, z, y, w, n_clusters, seed, min_events=10):
		ok = np.asarray(w) > 0
		self.kmeans = KMeans(n_clusters=int(n_clusters), n_init=10, random_state=seed).fit(z[ok])
		self.global_model = LogisticRegression(C=0.1, solver="lbfgs", max_iter=2000).fit(
			z[ok], y[ok], sample_weight=w[ok] / w[ok].mean()
		)
		lab = self.kmeans.predict(z)
		self.models = {}
		self.fallback_clusters = []
		for c in range(int(n_clusters)):
			use = ok & (lab == c)
			yy = np.asarray(y)[use]
			if use.sum() < 80 or yy.sum() < min_events or (len(yy) - yy.sum()) < min_events:
				self.models[c] = None
				self.fallback_clusters.append(c)
				continue
			fit = LogisticRegression(C=0.1, solver="lbfgs", max_iter=2000)
			fit.fit(z[use], yy, sample_weight=np.asarray(w)[use] / np.mean(np.asarray(w)[use]))
			self.models[c] = fit
		return self

	def predict(self, z, return_cluster=False):
		lab = self.kmeans.predict(z)
		p = np.empty(len(z), float)
		for c in np.unique(lab):
			use = lab == c
			m = self.models.get(int(c))
			m = self.global_model if m is None else m
			p[use] = m.predict_proba(z[use])[:, 1]
		return (p, lab) if return_cluster else p


@dataclass
class AlternativeRiskRoutes:
	knn_full: object
	knn_clean: object
	knn_full_k: int
	knn_clean_k: int
	moe: object
	moe_k: int
	calibrators: dict
	audit: dict

	def predict(self, model, X):
		z = model.prep.transform(X)
		raw = {
			"knn_unweighted_full": self.knn_full.predict(z, self.knn_full_k),
			"knn_unweighted_clean": self.knn_clean.predict(z, self.knn_clean_k),
		}
		raw["mixture_experts"], cluster = self.moe.predict(z, return_cluster=True)
		out = {k: self.calibrators[k].predict(v) for k, v in raw.items()}
		return out, cluster


def _fit_alternative_risk_routes(model, data):
	role = np.asarray(data["role"]).astype(str)
	X = np.asarray(data["X"], float)
	time = np.asarray(data["time"], float)
	event = np.asarray(data["event"], int)
	b = np.flatnonzero(role == "build")
	t = np.flatnonzero(role == "tune_model")
	cal = np.flatnonzero(role == "calibration_fit")
	z = model.prep.transform(X)
	yb, wb = model.km.labels_weights(time[b], event[b])
	yt, wt = model.km.labels_weights(time[t], event[t])
	yc, wc = model.km.labels_weights(time[cal], event[cal])
	eligible = np.asarray(model.donor_table.get("donor_eligible", np.ones(len(b), bool)), bool)
	clean = np.asarray(model.donor_table.get("retained", eligible), bool) & eligible
	wb = wb * eligible
	wt = wt * model.prep.technical_qc(X[t])
	wc = wc * model.prep.technical_qc(X[cal])
	fullbank = UnweightedKNNBank().fit(
		z[b], yb, wb, eligible, model.config.dimensions, model.config.seed + 2001
	)
	cleanbank = UnweightedKNNBank().fit(
		z[b], yb, wb, clean, model.config.dimensions, model.config.seed + 2002
	)
	kgrid = [k for k in (16, 32, 64, 128) if k <= eligible.sum()]
	if not kgrid:
		kgrid = [min(16, int(eligible.sum()))]

	def tune_bank(bank):
		best = None
		for k in kgrid:
			p = bank.predict(z[t], k)
			ok = wt > 0
			score = float(np.average(logloss_rows(yt[ok], p[ok]), weights=wt[ok]))
			if best is None or score < best[0]:
				best = (score, k)
		return best

	sf, kf = tune_bank(fullbank)
	sc, kc = tune_bank(cleanbank)
	# Mixture-of-experts cluster count is chosen on tune_model only.
	cand = []
	for q in os.getenv("ABM_MOE_CLUSTERS", "4,8,12").split(","):
		try:
			q = int(q.strip())
			if q >= 2 and q < max(3, eligible.sum() // 100):
				cand.append(q)
		except Exception:
			pass
	if not cand:
		cand = [4]
	mbest = None
	for q in cand:
		mm = ClusterMixtureExperts().fit(
			z[b], yb, wb, q, model.config.seed + 2100 + q, max(5, model.config.min_events // 2)
		)
		pp = mm.predict(z[t])
		ok = wt > 0
		score = float(np.average(logloss_rows(yt[ok], pp[ok]), weights=wt[ok]))
		if mbest is None or score < mbest[0]:
			mbest = (score, q, mm)
	_, mq, moe = mbest
	rawcal = {
		"knn_unweighted_full": fullbank.predict(z[cal], kf),
		"knn_unweighted_clean": cleanbank.predict(z[cal], kc),
		"mixture_experts": moe.predict(z[cal]),
	}
	calibrators = {n: Calibrator().fit(p, yc, wc) for n, p in rawcal.items()}
	audit = {
		"knn_full_k": int(kf),
		"knn_clean_k": int(kc),
		"knn_full_tune_logloss": float(sf),
		"knn_clean_tune_logloss": float(sc),
		"moe_clusters": int(mq),
		"moe_tune_logloss": float(mbest[0]),
		"moe_fallback_clusters": list(map(int, moe.fallback_clusters)),
		"test_outcomes_used": False,
		"geometry": "unweighted PCA KNN and X-only KMeans mixture-of-experts",
	}
	model.final_alt_risk = AlternativeRiskRoutes(
		fullbank, cleanbank, int(kf), int(kc), moe, int(mq), calibrators, audit
	)
	model.training_audit["final_alternative_risk_routes"] = audit
	return model


def _augment_alt_probability(model, X, probability):
	alt = getattr(model, "final_alt_risk", None)
	if alt is None:
		return dict(probability), None
	p = dict(probability)
	add, cluster = alt.predict(model, X)
	p.update(add)
	return p, cluster


@dataclass
class PredictabilityExtension:
	primary: str
	risk_bins: int
	models: dict
	rules: dict
	percentile_reference: dict
	risk_edges: dict
	risk_means: dict
	audit: dict
	consensus_models: tuple

	def raw_scores(self, model, probability, diagnostics):
		gx = model._gate_features(probability, diagnostics)
		out = {
			"risk_resid_error": -self.models["risk_resid_error"].predict(gx),
			"logloss_resid": -self.models["logloss_resid"].predict(gx),
			"local_stability": _local_stability(diagnostics),
		}
		out["model_consensus"], _ = _consensus(probability)
		return out

	def scores(self, model, probability, diagnostics):
		raw = self.raw_scores(model, probability, diagnostics)
		pct = {k: _pct(v, self.percentile_reference[k]) for k, v in raw.items()}
		hybrid = np.nanmean(
			np.vstack(
				[pct[k] for k in ("risk_resid_error", "logloss_resid", "local_stability", "model_consensus")]
			),
			axis=0,
		)
		return {**raw, "hybrid_predictability": hybrid}


def fit_predictability_extension(model, data):
	"""Attach frozen extension to an already fitted TwinABM.

    tune_gate outcomes fit the error meta-models; calibration_fit fixes score
    distributions and coverage thresholds. calibration_audit and test remain
    untouched by this function.
    """
	model = _fit_alternative_risk_routes(model, data)
	role = np.asarray(data["role"]).astype(str)
	X = np.asarray(data["X"], float)
	time = np.asarray(data["time"], float)
	event = np.asarray(data["event"], int)
	ids = np.asarray(data["eid"]).astype(str)
	groups = np.asarray(data["group"]).astype(str)
	a = model.config
	gate = role == "tune_gate"
	cal = role == "calibration_fit"
	if gate.sum() < 100 or cal.sum() < 100:
		raise ValueError("Predictability extension requires tune_gate and calibration_fit roles")

	pg, dg = model._calibrated(X[gate], ids[gate], groups[gate])
	pg, _ = _augment_alt_probability(model, X[gate], pg)
	y, w = model.km.labels_weights(time[gate], event[gate])
	known = (w > 0) & np.asarray(dg["technical_qc"], bool) & np.asarray(dg["supported"], bool)
	if known.sum() < max(100, 2 * a.min_events) or y[known].sum() < a.min_events:
		raise ValueError("Insufficient supported tune_gate outcomes for predictability extension")
	predictability_target_model = "elasticnet_full"
	base_risk = np.asarray(pg[predictability_target_model], float)
	sq = np.square(y - base_risk)
	ll = logloss_rows(y, base_risk)
	sqr, sq_edges, sq_means = _weighted_bin_residual(sq, base_risk, w * known, bins=max(5, 2 * a.risk_bins))
	llr, ll_edges, ll_means = _weighted_bin_residual(ll, base_risk, w * known, bins=max(5, 2 * a.risk_bins))
	gx = model._gate_features(pg, dg)

	models = {}
	for name, target in (("risk_resid_error", sqr), ("logloss_resid", llr)):
		fit = HistGradientBoostingRegressor(
			max_iter=max(150, a.gate_iterations),
			max_leaf_nodes=7,
			learning_rate=0.035,
			min_samples_leaf=max(20, min(150, int(known.sum() / 12))),
			l2_regularization=15,
			early_stopping=False,
			random_state=a.seed + (1701 if name == "risk_resid_error" else 1702),
		)
		fit.fit(gx[known], target[known], sample_weight=w[known] / np.mean(w[known]))
		models[name] = fit

	# Build a temporary extension to score calibration_fit, then freeze empirical
	# percentile transforms and all coverage rules without using audit/test labels.
	pc, dc = model._calibrated(X[cal], ids[cal], groups[cal])
	pc, _ = _augment_alt_probability(model, X[cal], pc)
	consensus_raw, consensus_names = _consensus(pc)
	tmp = PredictabilityExtension(
		primary="risk_resid_error",
		risk_bins=max(5, 2 * a.risk_bins),
		models=models,
		rules={},
		percentile_reference={},
		risk_edges={"squared": sq_edges, "logloss": ll_edges},
		risk_means={"squared": sq_means, "logloss": ll_means},
		audit={},
		consensus_models=tuple(consensus_names),
	)
	raw_cal = {
		"risk_resid_error": -models["risk_resid_error"].predict(model._gate_features(pc, dc)),
		"logloss_resid": -models["logloss_resid"].predict(model._gate_features(pc, dc)),
		"local_stability": _local_stability(dc),
		"model_consensus": consensus_raw,
	}
	eligible = np.asarray(dc["technical_qc"], bool) & np.asarray(dc["supported"], bool)
	refs = {k: np.sort(np.asarray(v)[eligible & np.isfinite(v)]) for k, v in raw_cal.items()}
	tmp.percentile_reference = refs
	scores_cal = tmp.scores(model, pc, dc)

	eligible = np.asarray(dc["technical_qc"], bool) & np.asarray(dc["supported"], bool)
	rules = {}
	for name, values in scores_cal.items():
		for q in a.coverages:
			rules[(name, float(q))] = FrozenRule().fit(values, eligible, ids[cal], q, a.seed + 1811)

	primary = os.getenv("ABM_FINAL_SELECTOR", "risk_resid_error").strip()
	if primary != "risk_resid_error":
		raise ValueError("C1 LOCKPAPER fixes primary selector risk_resid_error; alternative selectors remain secondary")
	ext = PredictabilityExtension(
		primary=primary,
		risk_bins=max(5, 2 * a.risk_bins),
		models=models,
		rules=rules,
		percentile_reference=refs,
		risk_edges={"squared": sq_edges, "logloss": ll_edges},
		risk_means={"squared": sq_means, "logloss": ll_means},
		audit={
			"primary": primary,
			"primary_predictability_selector": primary,
			"predictability_target_model": predictability_target_model,
			"predictability_target_definition": "OOF/frozen-model prediction error residualized within frozen predicted-risk strata",
			"routes": list(ROUTE_NAMES),
			"gate_N": int(gate.sum()),
			"gate_known_supported_N": int(known.sum()),
			"gate_cases": int(y[known].sum()),
			"calibration_N": int(cal.sum()),
			"target": "error residualized within frozen predicted-risk strata",
			"query_outcomes_used": False,
			"test_outcomes_used": False,
			"risk_only_is_not_a_predictability_route": True,
			"existing_core_selectors_retained": True,
		},
		consensus_models=tuple(consensus_names),
	)
	model.final_predictability = ext
	model.training_audit["final_predictability"] = ext.audit
	for key in ("primary_predictability_selector", "predictability_target_model", "predictability_target_definition"):
		model.training_audit[key] = ext.audit[key]
	print(f"[LE8 C1 LOCK] primary_predictability_selector={primary}; "
	      f"predictability_target_model={predictability_target_model}", flush=True)
	return model


def apply_predictability_extension(model, prediction, X, ids, groups):
	"""Add all deployable route masks and replace the primary selected mask.

    The original core selected mask is preserved for sensitivity analysis.
    """
	ext = getattr(model, "final_predictability", None)
	if ext is None:
		return prediction
	if (
		not np.array_equal(np.asarray(ids, str), np.asarray(prediction["eid"], str))
		or not np.array_equal(np.asarray(groups, str), np.asarray(prediction["group"], str))
		or len(X) != len(ids)
	):
		raise ValueError("Predictability extension must use the original prediction roster")
	p = dict(prediction)
	probs, moe_cluster = _augment_alt_probability(model, X, p["probability"])
	p["probability"] = probs
	diag = p["diagnostics"]
	if moe_cluster is not None:
		diag = dict(diag)
		diag["mixture_expert_cluster"] = np.asarray(moe_cluster, int)
		p["diagnostics"] = diag
	scores = ext.scores(model, probs, diag)
	masks = dict(p.get("masks", {}))
	core_scores = dict(p.get("selector_scores", {}))
	p["selected_core"] = np.asarray(p["selected"], bool).copy()
	p["core_selector"] = model.config.selector

	eligible = np.asarray(diag["technical_qc"], bool) & np.asarray(diag["supported"], bool)
	for name, values in scores.items():
		core_scores[name] = values
		for q in model.config.coverages:
			masks[f"{name}:{q:g}"] = ext.rules[(name, float(q))].apply(
				values, eligible, np.asarray(ids).astype(str)
			)

	selected = masks[f"{ext.primary}:{model.config.primary_coverage:g}"]
	p["selected"] = selected
	p["masks"] = masks
	p["selector_scores"] = core_scores
	p["primary_predictability_selector"] = ext.primary
	p["predictability_target_model"] = ext.audit.get("predictability_target_model", "unspecified")
	p["predictability_routes"] = list(ROUTE_NAMES)
	p["selection_semantics"] = (
		"Primary selector predicts development error relative to participants at similar predicted risk; "
		"all thresholds frozen on calibration_fit; query/test outcomes are unavailable to selection."
	)
	probs = dict(probs)
	probs["fallback_policy"] = np.where(
		np.asarray(diag["technical_qc"], bool),
		np.where(selected, probs[model.config.primary], probs["elasticnet_full"]),
		np.nan,
	)
	p["probability"] = probs
	p["reason"] = np.where(
		~np.asarray(diag["technical_qc"], bool),
		prediction["reason"],
		np.where(
			~np.asarray(diag["supported"], bool),
			"no_adequate_twins",
			np.where(selected, "predictability_selected", "predictability_deferred"),
		),
	)
	p["estimated_error"] = dict(p.get("estimated_error", {}))
	p["estimated_error"]["risk_residual_squared_error"] = -scores["risk_resid_error"]
	p["estimated_error"]["risk_residual_logloss"] = -scores["logloss_resid"]
	return p


def uno_c(time, event, probability, km):
	"""Horizon-truncated Uno C with strict earlier event times and half risk ties.

    Marginal censoring assumption matches the fixed-horizon IPCW implementation.
    Equal event times are not ordered. Returns comparable-pair diagnostics.
    """
	time, event, prob = map(np.asarray, (time, event, probability))
	levels = np.unique(prob)
	rank = np.searchsorted(levels, prob) + 1
	tree = np.zeros(len(levels) + 1)

	def add(k):
		while k < len(tree):
			tree[k] += 1
			k += k & -k

	def prefix(k):
		total = 0.0
		while k > 0:
			total += tree[k]
			k -= k & -k
		return total

	order = np.argsort(-time, kind="stable")
	num = den = 0.0
	later = pairs = start = 0
	while start < len(order):
		stop = start + 1
		while stop < len(order) and time[order[stop]] == time[order[start]]:
			stop += 1
		for i in order[start:stop]:
			if event[i] == 1 and time[i] <= km.horizon and later:
				g = float(km.at(np.array([time[i]]), left=True)[0])
				if g < km.minimum:
					return np.nan, pairs
				weight = 1 / g**2
				low = prefix(int(rank[i] - 1))
				equal = prefix(int(rank[i])) - low
				num += weight * (low + 0.5 * equal)
				den += weight * later
				pairs += later
		for i in order[start:stop]:
			add(int(rank[i]))
			later += 1
		start = stop
	return (num / den if den else np.nan), pairs


def harrell_c(time, event, score, sample_weight=None):
	"""Strict-time Harrell C in O(N log N); weights are family multiplicities.

	Higher scores indicate risk. Equal times are never compared. Pair weights
	are the product of subject multiplicities, never censoring IPCW.
	"""
	time, event, score = map(np.asarray, (time, event, score))
	weight = np.ones(len(time)) if sample_weight is None else np.asarray(sample_weight, float)
	if not (time.shape == event.shape == score.shape == weight.shape) or time.ndim != 1:
		raise ValueError("Harrell C requires equal one-dimensional arrays")
	if not all(np.isfinite(x).all() for x in (time, event, score, weight)) or np.any(weight < 0):
		raise ValueError("Harrell C requires finite inputs and nonnegative multiplicities")
	rank = np.searchsorted(np.unique(score), score) + 1
	tree = np.zeros(len(np.unique(score)) + 1)
	def prefix(k):
		total = 0.
		while k:
			total += tree[k]
			k -= k & -k
		return total
	order = np.argsort(-time, kind="stable")
	num = den = later = 0.
	start = 0
	while start < len(order):
		stop = start + 1
		while stop < len(order) and time[order[stop]] == time[order[start]]:
			stop += 1
		for i in order[start:stop]:
			if event[i] == 1:
				low = prefix(int(rank[i] - 1))
				equal = prefix(int(rank[i])) - low
				num += weight[i] * (low + .5 * equal)
				den += weight[i] * later
		for i in order[start:stop]:
			k = int(rank[i])
			while k < len(tree):
				tree[k] += weight[i]
				k += k & -k
			later += weight[i]
		start = stop
	return (num / den if den else np.nan), den


def _harrell_bootstrap_plan(time, event, score):
	"""Cache rank-tree sweeps for repeated family-weighted Harrell evaluation.

	Each distinct-score pair belongs to exactly one binary rank split. Within
	each split, cumulative lower-score weights count strictly later subjects.
	Risk ties are counted separately at half credit. No N-by-N pair matrix is
	built; storage is O(N log N), evaluation O(N log N), with NumPy sweeps.
	"""
	time, event, score = map(np.asarray, (time, event, score))
	order = np.argsort(-time, kind="stable")
	t, e = time[order], event[order]
	n = len(t)
	if not n:
		return lambda weight: (np.nan, 0.)
	rank = np.unique(score[order], return_inverse=True)[1]
	idx = np.arange(n)
	plans = []
	# Denominator, score ties, then each distinct-score binary split.
	groups = [(np.zeros(n, int), np.ones(n, bool), e == 1),
	          (rank, np.ones(n, bool), e == 1)]
	for bit in range(int(rank.max()).bit_length()):
		high = (rank & (1 << bit)) != 0
		groups.append((rank >> (bit + 1), ~high, (e == 1) & high))
	for group, supply, query in groups:
		ix = np.argsort(group, kind="stable")
		g, tt = group[ix], t[ix]
		new_group = np.r_[True, g[1:] != g[:-1]]
		new_time = new_group | np.r_[True, tt[1:] != tt[:-1]]
		group_start = np.maximum.accumulate(np.where(new_group, idx, 0))
		time_start = np.maximum.accumulate(np.where(new_time, idx, 0))
		plans.append((order[ix], supply[ix], query[ix], time_start, group_start))

	def evaluate(sample_weight):
		weight = np.asarray(sample_weight, float)
		if weight.shape != (n,) or not np.isfinite(weight).all() or np.any(weight < 0):
			raise ValueError("Invalid family multiplicities")
		counts = []
		for ix, supply, query, ts, gs in plans:
			w = weight[ix]
			prefix = np.r_[0., np.cumsum(w * supply)]
			counts.append(float(np.sum(w * query * (prefix[ts] - prefix[gs]))))
		den = counts[0]
		num = .5 * counts[1] + sum(counts[2:])
		return (num / den if den else np.nan), den
	return evaluate


def ipcw_information(y, w):
	"""Weight information, not a replacement for actual cases/families/pairs."""
	y, w = np.asarray(y), np.asarray(w, float)
	if y.shape != w.shape:
		raise ValueError("IPCW labels and weights must align")
	known = np.isfinite(w) & (w > 0)
	def ess(mask):
		v = w[mask & known]
		return float(v.sum() ** 2 / np.dot(v, v)) if len(v) else 0.0
	return dict(IPCW_weight_sum=float(w[known].sum()), IPCW_effective_N=ess(known),
		IPCW_case_effective_N=ess(y == 1), IPCW_control_effective_N=ess(y == 0))


def metrics(time, event, probability, km, *, include_harrell=True, include_uno=True):
	n = len(time)
	if not n:
		return dict(
			N=0,
			known_N=0,
			Harrell_C=np.nan,
			Harrell_comparable_pairs=0,
			Harrell_time_scope="full_follow_up_strict_time",
			AUC_IPCW=np.nan,
			Uno_C_horizon=np.nan,
			comparable_pairs=0,
			Brier_IPCW=np.nan,
			LogLoss_IPCW=np.nan,
			Brier_normalization="Hajek: sum(IPCW*loss)/sum(IPCW)",
			cases=0,
			controls=0,
			status="empty_subset",
			**ipcw_information(np.array([]), np.array([])),
		)
	time, event, probability = map(np.asarray, (time, event, probability))
	if not np.isfinite(probability).all():
		raise ValueError("Nonfinite probabilities inside evaluation mask")
	y, w = km.labels_weights(time, event)
	ok = w > 0
	row = dict(
		N=n,
		known_N=int(ok.sum()),
		cases=int(y.sum()),
		controls=int(np.sum(ok & (y == 0))),
		status="completed",
		Brier_normalization="Hajek: sum(IPCW*loss)/sum(IPCW)",
		Harrell_time_scope="full_follow_up_strict_time",
		**ipcw_information(y, w),
	)
	row["AUC_IPCW"] = (
		float(roc_auc_score(y[ok], probability[ok], sample_weight=w[ok])) if len(set(y[ok])) == 2 else np.nan
	)
	# Bootstrap callers can reuse the cached rank plan instead of sorting again.
	row["Harrell_C"], row["Harrell_comparable_pairs"] = (
		harrell_c(time, event, probability) if include_harrell else (np.nan, 0)
	)
	row["Uno_C_horizon"], row["comparable_pairs"] = (
		uno_c(time, event, probability, km) if include_uno else (np.nan, 0)
	)
	row["Brier_IPCW"] = (
		float(np.average((y[ok] - probability[ok]) ** 2, weights=w[ok])) if ok.any() else np.nan
	)
	row["LogLoss_IPCW"] = (
		float(np.average(logloss_rows(y[ok], probability[ok]), weights=w[ok])) if ok.any() else np.nan
	)
	if not np.isfinite(row["AUC_IPCW"]):
		row["status"] = "insufficient_case_control_information"
	return row


def risk_matched_mask(selected, eligible, risk, ids, risk_edges, seed):
	"""X/prediction-only negative control: no outcomes enter its selection.

    This is a descriptive, coverage-composition-matched control on current query X,
    not the primary deployable policy. Its membership can depend on batch makeup.
    """
	bins = np.searchsorted(risk_edges, risk, side="right")
	mask = np.zeros(len(ids), bool)
	for label in np.unique(bins):
		pool = np.flatnonzero(eligible & (bins == label))
		order = sorted(pool, key=lambda j: tie_hash(str(ids[j]), seed))
		count = int(np.sum(selected & (bins == label)))
		mask[np.asarray(order[:count], int)] = True
	return mask


def evaluate(model, prediction, time, event, truth_ids, bootstrap=None, full_curves=True, fit_calibration=True):
	"""AUC(selected) vs AUC(all) uses joint FAMILY resampling of the full cohort.

    Evaluates selection effects and same-subset training effects separately.
    No claim of individual certainty, external validation or guaranteed coverage.
    """
	a = model.config
	ids, groups = prediction["eid"], prediction["group"]
	if not np.array_equal(ids.astype(str), np.asarray(truth_ids).astype(str)):
		raise ValueError("Truth IDs/order must exactly match frozen prediction roster")
	time, event = np.asarray(time, float), np.asarray(event)
	if (
		len(time) != len(ids)
		or len(event) != len(ids)
		or not np.isfinite(time).all()
		or np.any(time <= 0)
		or not np.isin(event, [0, 1]).all()
	):
		raise ValueError("Invalid evaluation outcomes")
	if prediction["frozen_model_hash"] != model_hash(model):
		raise ValueError("Predictions do not match the supplied frozen model")
	p = {name: np.asarray(prob, float) for name, prob in prediction["probability"].items()}
	if not p or any(prob.shape != (len(ids),) for prob in p.values()):
		raise ValueError("Every model must provide one prediction per frozen participant")
	raw_qc = np.asarray(prediction["diagnostics"]["technical_qc"], bool)
	raw_selected = np.asarray(prediction["selected"], bool)
	if raw_qc.shape != (len(ids),) or raw_selected.shape != (len(ids),):
		raise ValueError("Invalid evaluation selection mask")
	requested_masks = dict(all=raw_qc, selected=raw_selected & raw_qc, rejected=raw_qc & ~raw_selected)
	if full_curves:
		for name, mask in prediction["masks"].items():
			mask = np.asarray(mask, bool)
			if mask.shape != (len(ids),):
				raise ValueError("Invalid frozen coverage mask: " + name)
			requested_masks[name] = mask & raw_qc
	# Freeze a common complete-prediction cohort once. Model-wise deletion would
	# turn a nominal same-subset contrast into a comparison of different people.
	common_finite = np.logical_and.reduce([np.isfinite(prob) for prob in p.values()])
	qc = raw_qc & common_finite
	masks = {name: mask & common_finite for name, mask in requested_masks.items()}
	selected = masks["selected"]
	masks["risk_matched_random"] = risk_matched_mask(
		selected, qc, p["elasticnet_full"], ids, model.risk_edges, a.seed + 412
	)
	pool = np.flatnonzero(qc)
	conf = np.abs(p["elasticnet_full"] - 0.5)
	order = sorted(pool, key=lambda j: (-conf[j], tie_hash(str(ids[j]), a.seed + 413)))
	masks["confidence_count_matched"] = np.zeros(len(ids), bool)
	masks["confidence_count_matched"][np.asarray(order[: int(selected.sum())], int)] = True
	# Count-matched controls are constructed on the common cohort, preserving
	# their intended selected/control counts after a finite-prediction exclusion.
	for name in ("risk_matched_random", "confidence_count_matched"):
		requested_masks[name] = masks[name].copy()
	y, w = model.km.labels_weights(time, event)
	rows = []
	for subset, mask in masks.items():
		for name, prob in p.items():
			valid = mask
			row = metrics(time[valid], event[valid], prob[valid], model.km)
			row.update(
				model=name,
				subset=subset,
				requested_N=int(requested_masks[subset].sum()),
				evaluated_N=int(valid.sum()),
				excluded_N=int(requested_masks[subset].sum() - valid.sum()),
				excluded_nonfinite_N=int(np.sum(requested_masks[subset] & ~common_finite)),
				requested_roster_hash=fingerprint(ids[requested_masks[subset]]),
				evaluation_roster_hash=fingerprint(ids[valid]),
				same_people_within_subset=True,
				coverage=float(valid.mean()),
				case_coverage=float(y[valid].sum() / y.sum()) if y.sum() else np.nan,
				control_coverage=(
					float(np.sum((w[valid] > 0) & (y[valid] == 0)) / np.sum((w > 0) & (y == 0)))
					if np.sum((w > 0) & (y == 0))
					else np.nan
				),
				IPCW_case_coverage=(
					float(np.sum(w[valid] * y[valid]) / np.sum(w * y)) if np.sum(w * y) > 0 else np.nan
				),
				family_N=len(set(groups[valid])),
				horizon=a.horizon,
			)
			rows.append(row)
	lookup = {(r["model"], r["subset"]): r for r in rows}
	primary, baseline = a.primary, "elasticnet_full"
	pairs = [
		("selection_effect", baseline, "selected", baseline, "all"),
		("training_effect_same_subset", primary, "selected", baseline, "selected"),
		("combined_vs_bulk", primary, "selected", baseline, "all"),
		("all_cohort_training_effect", primary, "all", baseline, "all"),
		("fallback_policy", "fallback_policy", "all", baseline, "all"),
		("twin_library_cleaning", "twin_clean", "selected", "twin_full", "selected"),
		("twin_integrated_library_cleaning", "twin_integrated_clean", "selected", "twin_integrated_full", "selected"),
		("molecular_vs_integrated_twin", "twin_clean", "selected", "twin_integrated_clean", "selected"),
		("beyond_risk_composition", baseline, "selected", baseline, "risk_matched_random"),
		("beyond_confidence_selection", baseline, "selected", baseline, "confidence_count_matched"),
	]
	if "elasticnet_random" in p:
		pairs.append(("pruning_specificity", "elasticnet_clean", "selected", "elasticnet_random", "selected"))
	for clean, random, label in (
		("twin_clean", "twin_random", "twin_pruning_specificity"),
		("twin_integrated_clean", "twin_integrated_random", "twin_integrated_pruning_specificity"),
	):
		if clean in p and random in p:
			pairs.append((label, clean, "selected", random, "selected"))
			pairs.append((label + "_all", clean, "all", random, "all"))
	if "lightgbm_full" in p:
		pairs.extend(
			[
				("selected_vs_full_lightgbm", primary, "selected", "lightgbm_full", "all"),
				("same_subset_vs_lightgbm", primary, "selected", "lightgbm_full", "selected"),
			]
		)
	metric_names = ["Harrell_C", "AUC_IPCW", "Uno_C_horizon", "Brier_IPCW", "LogLoss_IPCW"]
	B = a.bootstrap if bootstrap is None else int(bootstrap)
	if B < 1:
		raise ValueError("bootstrap must be positive")
	u, inverse = np.unique(groups, return_inverse=True)
	order = np.argsort(inverse, kind="stable")
	positions = np.split(order, np.cumsum(np.bincount(inverse))[:-1])
	rng = np.random.default_rng(a.seed + 301)
	draws = {(name, m): [] for name, *_ in pairs for m in metric_names}
	needed = set((mn, sn) for _, m1, s1, m2, s2 in pairs for mn, sn in [(m1, s1), (m2, s2)])
	harrell_plans = {
		(mn, sn): _harrell_bootstrap_plan(time[masks[sn]], event[masks[sn]], p[mn][masks[sn]])
		for mn, sn in needed
	}
	for b in range(B):
		family_draw = rng.integers(0, len(u), len(u))
		ix = np.concatenate([positions[j] for j in family_draw])
		multiplicity = np.bincount(family_draw, minlength=len(u))[inverse]
		batch = {}
		for mn, sn in needed:
			jj = ix[masks[sn][ix]]
			batch[(mn, sn)] = metrics(time[jj], event[jj], p[mn][jj], model.km, include_harrell=False)
			batch[(mn, sn)]["Harrell_C"], _ = harrell_plans[(mn, sn)](multiplicity[masks[sn]])
		for label, m1, s1, m2, s2 in pairs:
			for metric in metric_names:
				draws[(label, metric)].append(batch[(m1, s1)][metric] - batch[(m2, s2)][metric])
	contrasts = []
	for label, m1, s1, m2, s2 in pairs:
		for metric in metric_names:
			values = np.asarray(draws[(label, metric)])
			values = values[np.isfinite(values)]
			minimum_valid = max(20, math.ceil(0.8 * B))
			enough = len(values) >= minimum_valid
			bounds = np.quantile(values, [0.025, 0.975]) if enough else [np.nan, np.nan]
			delta = lookup[(m1, s1)][metric] - lookup[(m2, s2)][metric]
			inference_status = ("insufficient_information" if not np.isfinite(delta) else
				"evaluated" if enough else "insufficient_bootstrap")
			contrasts.append(
				dict(
					contrast=label,
					metric=metric,
					model=m1,
					subset=s1,
					reference=m2,
					reference_subset=s2,
					delta=delta,
					same_people=s1 == s2,
					comparison_design="same_subset_model_comparison" if s1 == s2 else "selection_composition_comparison",
					model_evaluated_N=lookup[(m1, s1)]["N"],
					reference_evaluated_N=lookup[(m2, s2)]["N"],
					model_roster_hash=lookup[(m1, s1)]["evaluation_roster_hash"],
					reference_roster_hash=lookup[(m2, s2)]["evaluation_roster_hash"],
					lower=float(bounds[0]),
					upper=float(bounds[1]),
					valid_boot=len(values),
					requested_boot=B,
					minimum_valid_boot=minimum_valid,
					status=inference_status,
					inference_status=inference_status,
					point_status="evaluated" if np.isfinite(delta) else "insufficient_information",
					bootstrap_unit=a.grouping_unit,
					uncertainty="conditional on frozen models, selection rules and build-only censoring model",
					multiplicity="predeclared primary contrast only; other intervals exploratory",
				)
			)
	c = {(r["contrast"], r["metric"]): r for r in contrasts}
	sel = lookup[(primary, "selected")]
	info = sel["cases"] >= a.audit_min_cases and sel["controls"] >= a.audit_min_controls

	def supported(label, metric="AUC_IPCW", lower=True):
		rr = c[(label, metric)]
		return bool(
			info
			and rr["valid_boot"] >= max(200, math.ceil(0.8 * B))
			and (rr["lower"] > 0 if lower else rr["upper"] < 0)
		)

	status = dict(
		information_adequate=bool(info),
		selected_cases=sel["cases"],
		selected_controls=sel["controls"],
		original_N=len(ids),
		technical_QC_excluded=int((~raw_qc).sum()),
		technical_QC_N=int(raw_qc.sum()),
		evaluation_N=int(qc.sum()),
		excluded_nonfinite_N=int(np.sum(raw_qc & ~common_finite)),
		evaluation_roster_hash=fingerprint(ids[qc]),
		requested_selected_N=int(requested_masks["selected"].sum()),
		selected_N=int(selected.sum()),
		actual_coverage=float(selected.mean()),
		no_adequate_twins_N=int((qc & ~prediction["diagnostics"]["supported"]).sum()),
		selective_discrimination_supported=supported("combined_vs_bulk"),
		selection_only_supported=supported("selection_effect"),
		training_gain_supported=supported("training_effect_same_subset"),
		absolute_error_improvement_supported=supported("combined_vs_bulk", "Brier_IPCW", False),
		beyond_risk_composition_supported=supported("beyond_risk_composition"),
		beyond_confidence_selection_supported=supported("beyond_confidence_selection"),
		clinical_release=False,
		external_validation=False,
		interpretation="Independent audit can validate a subset policy even without same-subset algorithm superiority; no individual guarantee",
	)
	# Transparent test calibration diagnostics only; never used to update predictions.
	calibration = []
	for subset in ("all", "selected", "rejected"):
		for name in (primary, baseline):
			ok = masks[subset] & (w > 0)
			yy, pp, ww = y[ok], p[name][ok], w[ok]
			row = dict(model=name, subset=subset, N=int(ok.sum()), evaluation_only=True)
			denominator = float(np.sum(ww * pp))
			row["observed_expected"] = float(np.sum(ww * yy) / denominator) if denominator > 0 else np.nan
			row["diagnostic_regression_fitted"] = False
			if len(set(yy)) == 2 and fit_calibration:
				cc = Calibrator().fit(pp, yy, ww)
				row.update(
					intercept=float(cc.model.intercept_[0]),
					slope=float(cc.model.coef_[0, 0]),
					diagnostic_regression_fitted=True,
					status="evaluation_only_no_prediction_update",
				)
			elif not fit_calibration:
				row["status"] = "diagnostic_regression_not_fitted_report_only"
			else:
				row["status"] = "insufficient_information"
			calibration.append(row)
	strata_rows = []
	strata = {
		"baseline_risk_decile": np.searchsorted(model.risk_edges, p[baseline], side="right"),
		"frozen_quality_group": np.searchsorted(
			model.quality_edges, prediction["selector_scores"][a.selector], side="right"
		),
	}
	for kind, labels in strata.items():
		for label in np.unique(labels):
			for subset, mask in [("all", qc), ("selected", selected)]:
				use = mask & (labels == label)
				if not use.any():
					continue
				pr = metrics(time[use], event[use], p[primary][use], model.km)
				br = metrics(time[use], event[use], p[baseline][use], model.km)
				known = use & (w > 0)
				strata_rows.append(
					dict(
						stratum_type=kind,
						stratum=int(label),
						subset=subset,
						N=int(use.sum()),
						cases=int(y[use].sum()),
						cutpoints_source="calibration_fit",
						candidate_brier=pr["Brier_IPCW"],
						baseline_brier=br["Brier_IPCW"],
						brier_gain=br["Brier_IPCW"] - pr["Brier_IPCW"],
						candidate_AUC=pr["AUC_IPCW"],
						baseline_AUC=br["AUC_IPCW"],
						mean_predicted_squared_error=(
							float(
								np.average(
									prediction["estimated_error"]["absolute_error"][known], weights=w[known]
								)
							)
							if known.any()
							else np.nan
						),
					)
				)
	confidence_rows = []
	confidence = prediction["diagnostics"]["pre_prediction_confidence"]
	# Fixed score intervals: no outcome-driven binning and no score calibration on test.
	for lo in range(0, 100, 10):
		mask = qc & (confidence >= lo) & ((confidence < lo + 10) if lo < 90 else confidence <= 100)
		if mask.sum() < 10:
			continue
		for name in (primary, baseline):
			rr = metrics(time[mask], event[mask], p[name][mask], model.km)
			rr.update(
				model=name,
				confidence_lo=lo,
				confidence_hi=lo + 10,
				mean_confidence=float(np.mean(confidence[mask])),
				supported_fraction=float(np.mean(prediction["diagnostics"]["supported"][mask])),
			)
			confidence_rows.append(rr)
	return dict(
		metrics=rows,
		contrasts=contrasts,
		prediction_finiteness=[dict(
			model=name, requested_N=int(raw_qc.sum()),
			model_nonfinite_N=int(np.sum(raw_qc & ~np.isfinite(prob))),
			evaluated_N=int(qc.sum()), excluded_nonfinite_N=int(np.sum(raw_qc & ~common_finite)),
			evaluation_roster_hash=fingerprint(ids[qc]),
		) for name, prob in p.items()],
		status=status,
		calibration=calibration,
		grouping_unit=a.grouping_unit,
		confidence_validation=confidence_rows,
		strata_diagnostics=strata_rows,
		roster_hash=fingerprint(ids),
		version=prediction["version"],
		layer=a.layer,
		trait=a.trait,
		run_id=prediction["frozen_model_hash"],
		evaluation_scope=prediction.get("evaluation_scope", "unspecified"),
		note="Censoring weights are marginal, and inference is internal/conditional; refit-bootstrap and external validation are additional analyses",
	)


# 🚩 run abm twin
# Separate fit -> project -> evaluate commands. Only fit sees development labels.
import argparse
from threadpoolctl import threadpool_limits

# LE8_FINAL::2026-10-15.final-multiroute-v1


def read_npz(path):
	with np.load(path, allow_pickle=False) as f:
		return {k: f[k] for k in f.files}


def parser():
	p = argparse.ArgumentParser(description=__doc__)
	s = p.add_subparsers(dest="command", required=True)
	a = s.add_parser("fit")
	a.add_argument("--development", required=True)
	a.add_argument("--output", required=True)
	a.add_argument("--config")
	a = s.add_parser("project")
	a.add_argument("--model", required=True)
	a.add_argument("--features", required=True)
	a.add_argument("--output", required=True)
	a = s.add_parser("evaluate")
	a.add_argument("--model", required=True)
	a.add_argument("--predictions", required=True)
	a.add_argument("--truth", required=True)
	a.add_argument("--output", required=True)
	a.add_argument("--bootstrap", type=int)
	a = s.add_parser("demo")
	a.add_argument("--output", required=True)
	a.add_argument("--seed", type=int, default=12345)
	return p


def load_config(path):
	values = json.loads(Path(path).read_text()) if path else {}
	for key in ("coverages", "c_grid", "always_include_indices"):
		if key in values:
			values[key] = tuple(values[key])
	return Config(**values)


def write_result(out, result):
	out.mkdir(parents=True, exist_ok=True)
	for name, value in result.items():
		json_write(out / (name + ".json"), value)
	import pandas as pd
	from le8_common import write_tables

	tables = {
		k: pd.DataFrame(v)
		for k, v in result.items()
		if isinstance(v, list) and (not v or isinstance(v[0], dict))
	}
	tables["status"] = pd.DataFrame([result["status"]])
	for d in tables.values():
		d["evaluation_scope"] = result.get("evaluation_scope", "unspecified")
		d["grouping_unit"] = result.get("grouping_unit", "family")
	write_tables(
		out,
		tables,
		"ABM twin",
		result.get("layer", "prot"),
		result.get("trait", "cvd_cad"),
		run_id=result.get("run_id", ""),
		workbook="ABM_twin.results.xlsx",
	)
	plot_abm(tables, out / "figures")


def fit_command(development, output, config):
	out = Path(output)
	if (out / "model_LOCAL_ONLY.joblib").exists():
		raise FileExistsError("Use a new output directory; old fits are never overwritten")
	data = read_npz(development)
	required = {"X", "time", "event", "eid", "group", "role", "feature_names"}
	if not required <= set(data):
		raise ValueError("Missing development fields: " + str(required - set(data)))
	model = TwinABM(config).fit(
		data["X"],
		data["time"],
		data["event"],
		data["eid"],
		data["group"],
		data["role"],
		data["feature_names"],
	)
	model = fit_predictability_extension(model, data)
	out.mkdir(parents=True, exist_ok=True)
	joblib.dump(model, out / "model_LOCAL_ONLY.joblib", compress=3)
	json_write(
		out / "MODEL_FROZEN.json",
		dict(version=VERSION, hash=model_hash(model), config=asdict(config), test_outcomes_used=False),
	)
	np.savez_compressed(out / "donor_quality_LOCAL_ONLY.npz", **model.donor_table)
	json_write(out / "training_audit.json", model.training_audit)
	json_write(out / "model_tuning.json", model.model_tuning)
	json_write(
		out / "selected_panel.json", dict(features=model.panel_names, input_indices=model.panel_indices)
	)
	import pandas as pd
	from le8_common import write_tables

	audit = model.training_audit
	summary = {k: v for k, v in audit.items() if not isinstance(v, (list, dict))}
	write_tables(
		out / "aggregate",
		{
			"training_summary": pd.DataFrame([summary]),
			"training_folds": pd.DataFrame(audit["oof"]),
			"pruning_audit": pd.DataFrame(audit["pruning"]),
			"selected_panel": pd.DataFrame(
				{
					"feature": model.panel_names,
					"input_index": model.panel_indices,
					"fixed_covariate": [j in config.always_include_indices for j in model.panel_indices],
				}
			),
			"model_tuning": pd.DataFrame(model.model_tuning),
			"role_counts": pd.DataFrame(
				[
					dict(role=r, N=v["N"], cases=v["cases"], technical_QC_excluded=v["technical_QC_excluded"])
					for r, v in audit["roles"].items()
				]
			),
		},
		"ABM twin",
		config.layer,
		config.trait,
		model_hash(model),
	)
	ix = data["role"] == "calibration_audit"
	before = model_hash(model)
	pred = model.predict(
		data["X"][ix], data["eid"][ix], data["group"][ix], data["feature_names"], allow_development_audit=True
	)
	pred = apply_predictability_extension(model, pred, data["X"][ix], data["eid"][ix], data["group"][ix])
	result = evaluate(model, pred, data["time"][ix], data["event"][ix], data["eid"][ix], full_curves=False)
	assert before == model_hash(model), "Audit changed the frozen model"
	write_result(out / "audit", result)
	return model


def demo_command(out, seed):
	"""Synthetic mechanics demonstration only; no clinical performance claim."""
	out = Path(out)
	out.mkdir(parents=True, exist_ok=True)
	rng = np.random.default_rng(seed)
	n = 2400
	p = 12
	X = rng.normal(size=(n, p))
	# Ordinary stochastic survival data; deliberately not engineered to guarantee a pruning gain.
	eta = -2.3 + 1.1 * X[:, 0] - 0.7 * X[:, 1] + 0.4 * X[:, 2]
	rate = np.exp(eta) / 10
	failure = rng.exponential(1 / rate)
	censor = rng.uniform(11, 18, n)
	time = np.minimum(failure, censor)
	event = (failure <= censor).astype(int)
	ids = np.array([f"synthetic_{i}" for i in range(n)])
	groups = ids.copy()
	names = np.array([f"X{j}" for j in range(p)])
	# Split is outcome-blind. This demo is not a UKB result and is not tuned for gains.
	order = np.argsort([tie_hash(v, seed) for v in ids])
	dev = order[: n // 2]
	test = order[n // 2 :]
	counts = np.array([0.6, 0.1, 0.1, 0.1, 0.1])
	cut = np.cumsum(np.r_[0, np.round(counts * len(dev)).astype(int)])
	cut[-1] = len(dev)
	roles = np.empty(len(dev), dtype="U25")
	for j, r in enumerate(ROLES):
		roles[cut[j] : cut[j + 1]] = r
	np.savez_compressed(
		out / "development.npz",
		X=X[dev],
		time=time[dev],
		event=event[dev],
		eid=ids[dev],
		group=groups[dev],
		role=roles,
		feature_names=names,
	)
	np.savez_compressed(out / "test_X.npz", X=X[test], eid=ids[test], group=groups[test], feature_names=names)
	np.savez_compressed(out / "test_truth.npz", time=time[test], event=event[test], eid=ids[test])
	cfg = Config(
		seed=seed,
		folds=3,
		repeats=1,
		c_grid=(0.1,),
		min_events=3,
		k=32,
		min_twins=8,
		min_family_ess=6,
		dimensions=8,
		gate_iterations=20,
		bootstrap=25,
		audit_min_cases=3,
		audit_min_controls=3,
		max_iter=3000,
	)
	model = fit_command(out / "development.npz", out / "fit", cfg)
	pred = model.predict(X[test], ids[test], groups[test], names)
	pred = apply_predictability_extension(model, pred, X[test], ids[test], groups[test])
	joblib.dump(pred, out / "test_predictions_LOCAL_ONLY.joblib")
	result = evaluate(model, pred, time[test], event[test], ids[test], bootstrap=25)
	write_result(out / "test_report", result)
	print(
		json.dumps(
			{
				"synthetic_only": True,
				"selected_n": int(pred["selected"].sum()),
				"test_n": len(test),
				"status": result["status"],
			},
			indent=2,
		)
	)


def run_main():
	a = parser().parse_args()
	if a.command == "fit":
		config = load_config(a.config)
		with threadpool_limits(config.cores):
			fit_command(a.development, a.output, config)
	elif a.command == "demo":
		with threadpool_limits(2):
			demo_command(a.output, a.seed)
	else:
		model = load_model(a.model)  # Only load a trusted model created by this package.
		with threadpool_limits(model.config.cores):
			if a.command == "project":
				data = read_npz(a.features)
				permitted = {"X", "eid", "group", "feature_names"}
				if set(data) != permitted:
					raise ValueError(
						"Projection file must contain exactly X, eid, group, feature_names; no outcomes"
					)
				pred = model.predict(data["X"], data["eid"], data["group"], data["feature_names"])
				pred = apply_predictability_extension(model, pred, data["X"], data["eid"], data["group"])
				output = Path(a.output)
				if output.exists():
					raise FileExistsError(output)
				output.parent.mkdir(parents=True, exist_ok=True)
				joblib.dump(pred, output, compress=3)
				import pandas as pd

				table = pd.DataFrame(
					{
						"eid": pred["eid"],
						"selected": pred["selected"],
						"reason": pred["reason"],
						"candidate_risk": pred["probability"][model.config.primary],
						"fallback_risk": pred["probability"]["elasticnet_full"],
						"policy_risk": pred["probability"]["fallback_policy"],
						**pred["diagnostics"],
					}
				)
				table.to_csv(output.with_name(output.stem + ".confidence.tsv.gz"), sep="\t", index=False)
				json_write(
					output.with_suffix(".manifest.json"),
					dict(
						model_hash=model_hash(model),
						version=VERSION,
						outcomes_used=False,
						selected_N=int(pred["selected"].sum()),
					),
				)
			else:
				pred = joblib.load(a.predictions)
				truth = read_npz(a.truth)
				before = model_hash(model)
				result = evaluate(
					model, pred, truth["time"], truth["event"], truth["eid"], bootstrap=a.bootstrap
				)
				assert before == model_hash(model)
				write_result(Path(a.output), result)


from dataclasses import replace
import pandas as pd
from le8_common import id_table, sha256, write_tables


def prepare_main():
	a = argparse.ArgumentParser(description=__doc__)
	a.add_argument("--cohort-file", required=True)
	a.add_argument("--omics-file", required=True)
	a.add_argument("--config", required=True)
	a.add_argument("--output", required=True)
	a.add_argument("--time-col", default="time")
	a.add_argument("--event-col", default="event")
	a.add_argument("--prevalent-col", default="known_cad_at_index")
	a.add_argument("--clinical", default="")
	a.add_argument(
		"--group-col", help="Optional family/group column; absent: unverified individual ID groups"
	)
	a.add_argument(
		"--assume-unrelated",
		action="store_true",
		help="Optional user declaration, not proof of independent families; unnecessary for running",
	)
	a.add_argument("--role-col")
	x = a.parse_args()
	cfg = load_config(x.config)
	d = id_table(x.cohort_file)
	m = id_table(x.omics_file)
	out = Path(x.output)
	if out.exists() and any(out.iterdir()):
		raise FileExistsError("Prepared output already exists")
	if not {x.time_col, x.event_col, x.prevalent_col} <= set(d):
		raise ValueError(
			"Provide already-derived survival time/event/prevalence from the host outcome code; absent fields are not inferred"
		)
	ids = d.index[d.index.isin(m.index)]
	d = d.loc[ids]
	m = m.loc[ids].drop(columns=["eid"])
	time = pd.to_numeric(d[x.time_col], errors="coerce").to_numpy()
	event = pd.to_numeric(d[x.event_col], errors="coerce").to_numpy()
	prev = pd.to_numeric(d[x.prevalent_col], errors="coerce").to_numpy()
	eligible = (prev == 0) & np.isfinite(time) & (time > 0) & np.isin(event, [0, 1])
	cls = [v.strip() for v in x.clinical.split(",") if v.strip()]
	if set(cls) & {x.time_col, x.event_col, x.prevalent_col}:
		raise ValueError("Outcome leakage into clinical inputs")
	C = d[cls].apply(pd.to_numeric, errors="raise").to_numpy(float) if cls else np.empty((len(d), 0))
	X = np.c_[m.apply(pd.to_numeric, errors="raise").to_numpy(float), C]
	names = list(m.columns) + ["clinical::" + v for v in cls]
	if len(set(names)) != len(names):
		raise ValueError("Duplicate predictor names")
	if x.group_col:
		if x.group_col not in d or d[x.group_col].isna().any():
			raise ValueError("Missing family-group IDs")
		groups = d[x.group_col].astype(str).to_numpy()
	else:
		groups = np.asarray(ids, str)
	support = dict(
		family_status="verified_supplied_groups" if x.group_col else "unavailable_individual_ID_groups_only",
		independent_family_support_verified=bool(x.group_col),
		grouping_unit="family" if x.group_col else "individual_unverified",
		assume_unrelated=bool(x.assume_unrelated),
		validation_status="exploratory",
	)
	if x.role_col:
		role = d[x.role_col].astype(str).to_numpy()
	else:
		v = np.array([tie_hash(g, cfg.seed) / 2**64 for g in groups])
		role = np.full(len(v), "test", dtype="U25")
		edges = [0, 0.50, 0.60, 0.68, 0.76, 0.84]
		for j, r in enumerate(ROLES):
			role[(v >= edges[j]) & (v < edges[j + 1])] = r
	if not set(role) <= set(ROLES) | {"test"}:
		raise ValueError(
			"Invalid roles; use exact explicit build/tune_model/tune_gate/calibration_fit/calibration_audit/test"
		)
	rg = pd.DataFrame({"g": groups, "r": role})
	if rg.groupby("g").r.nunique().max() > 1:
		raise ValueError("Family crosses roles")
	groups = np.asarray(groups, str)
	role = np.asarray(role, str)
	X, time, event, groups, role, ids = [
		v[eligible] for v in (X, time, event, groups, role, np.asarray(ids, str))
	]
	event = event.astype(int)
	dev = role != "test"
	if set(role) != set(ROLES) | {"test"}:
		raise ValueError("Each development role and test must have eligible people")
	cfg = replace(
		cfg,
		always_include_indices=tuple(range(len(m.columns), X.shape[1])),
		grouping_unit=support["grouping_unit"],
	)
	out.mkdir(parents=True)
	local = out / "_LOCAL_ONLY"
	local.mkdir(mode=0o700)
	np.savez_compressed(
		local / "development.npz",
		X=X[dev],
		time=time[dev],
		event=event[dev],
		eid=ids[dev],
		group=groups[dev],
		role=role[dev],
		feature_names=np.asarray(names),
	)
	np.savez_compressed(
		local / "test_X.npz", X=X[~dev], eid=ids[~dev], group=groups[~dev], feature_names=np.asarray(names)
	)
	np.savez_compressed(local / "test_truth.npz", time=time[~dev], event=event[~dev], eid=ids[~dev])
	json_write(out / "config.json", asdict(cfg))
	json_write(
		out / "input_provenance.json",
		dict(
			cohort_sha256=sha256(x.cohort_file),
			omics_sha256=sha256(x.omics_file),
			time_column=x.time_col,
			event_column=x.event_col,
			prevalent_column=x.prevalent_col,
			group_assumption=x.group_col
			or ("declared_unrelated" if x.assume_unrelated else "individual_ID_only_unverified"),
			**support,
			outcome_definition="unchanged supplied host survival fields; death handled according to supplied event/censoring encoding",
		),
	)
	counts = pd.DataFrame({"role": role}).groupby("role").size().reset_index(name="N")
	write_tables(
		out / "aggregate",
		{"cohort_roles": counts, "support_status": pd.DataFrame([support])},
		"ABM twin",
		cfg.layer,
		cfg.trait,
		"prepared",
	)
	print(
		json.dumps(
			dict(
				output=str(out),
				eligible_N=len(X),
				development_N=int(dev.sum()),
				test_N=int((~dev).sum()),
				input_features=X.shape[1],
				molecular_panel_budget=cfg.panel_size,
			)
		)
	)


def run_prepared_numeric_cohort(
	*,
	X,
	time,
	event,
	eid,
	family_group,
	role,
	feature_names,
	clinical_feature_names,
	output,
	config: Config = Config(),
	confirmatory_test=False,
):
	"""Keep the existing roster; accept exactly five development roles plus test.

    clinical_feature_names must be explicit so a protein-only experiment cannot
    silently be described as a clinical-plus-protein baseline comparison.
    confirmatory_test is metadata ONLY; this function cannot verify that the
    user has never inspected that test's labels/results in earlier experiments.
    """
	out = Path(output)
	if out.exists() and any(out.iterdir()):
		raise FileExistsError("Use a fresh versioned directory; never overwrite previous ABM results")
	X = np.asarray(X, float)
	ids = np.asarray(eid).astype(str)
	groups = np.asarray(family_group).astype(str)
	roles = np.asarray(role).astype(str)
	names = list(map(str, feature_names))
	if set(roles) != set(ROLES) | {"test"}:
		raise ValueError("Missing/unknown roles; the adapter never creates a new test split")
	if not clinical_feature_names or not set(clinical_feature_names) <= set(names):
		raise ValueError("Declare included numeric baseline clinical features; no silently weaker comparator")
	if len(set(names)) != len(names):
		raise ValueError("Duplicate feature names")
	forbidden = {
		"time",
		"event",
		"outcome",
		"target",
		"prevalent_status",
		"role",
		"split",
		"diagnosis_date",
		"t2e",
		"Yt2e",
	}
	if forbidden & set(names) or any(name.lower().endswith((".t2e", ".yt2e", ".bi2e")) for name in names):
		raise ValueError("Outcome or split field in X")
	config = replace(
		config, always_include_indices=tuple(names.index(str(n)) for n in clinical_feature_names)
	)
	dev = roles != "test"
	test = ~dev
	if set(groups[dev]) & set(groups[test]):
		raise ValueError("Family overlap with test")
	out.mkdir(parents=True, exist_ok=True)
	np.savez_compressed(
		out / "development_LOCAL_ONLY.npz",
		X=X[dev],
		time=np.asarray(time)[dev],
		event=np.asarray(event)[dev],
		eid=ids[dev],
		group=groups[dev],
		role=roles[dev],
		feature_names=np.asarray(names),
	)
	# Test outcomes are in a separate file and never passed to fit/project.
	np.savez_compressed(
		out / "test_truth_LOCAL_ONLY.npz",
		time=np.asarray(time)[test],
		event=np.asarray(event)[test],
		eid=ids[test],
	)
	np.savez_compressed(
		out / "test_X_LOCAL_ONLY.npz",
		X=X[test],
		eid=ids[test],
		group=groups[test],
		feature_names=np.asarray(names),
	)
	model = fit_command(out / "development_LOCAL_ONLY.npz", out / "fit", config)
	pred = model.predict(X[test], ids[test], groups[test], names)
	pred = apply_predictability_extension(model, pred, X[test], ids[test], groups[test])
	import joblib

	joblib.dump(pred, out / "test_predictions_LOCAL_ONLY.joblib", compress=3)
	import pandas as pd

	pd.DataFrame(
		{
			"eid": pred["eid"],
			"selected": pred["selected"],
			"reason": pred["reason"],
			"candidate_risk": pred["probability"][model.config.primary],
			"fallback_risk": pred["probability"]["elasticnet_full"],
			"policy_risk": pred["probability"]["fallback_policy"],
			**pred["diagnostics"],
		}
	).to_csv(out / "test_predictions_LOCAL_ONLY.confidence.tsv.gz", sep="\t", index=False)
	result = evaluate(model, pred, np.asarray(time)[test], np.asarray(event)[test], ids[test])
	result["status"]["test_provenance"] = (
		"user_declared_untouched" if confirmatory_test else "previously_inspected_or_unverified_exploratory"
	)
	result["status"]["clinical_release"] = False
	write_result(out / "test_report", result)
	json_write(
		out / "host_contract.json",
		dict(
			original_roster_preserved=True,
			clinical_feature_names=list(clinical_feature_names),
			numeric_schema=names,
			all_ids_hash=fingerprint(ids),
			test_ids_hash=fingerprint(ids[test]),
			UKB_phenotype_definition="supplied by host; not reconstructed by this adapter",
			preprocessing="numeric reference core; production categorical/technical/tail adapter requires validation",
			config=asdict(config),
			confirmed_independent_external_validation=False,
		),
	)
	return result


def plot_abm(tables, out):
	import matplotlib

	matplotlib.use("Agg")
	import matplotlib.pyplot as plt

	out = Path(out)
	out.mkdir(parents=True, exist_ok=True)
	t = tables.get("metrics", pd.DataFrame())
	if len(t):
		fig, ax = plt.subplots(figsize=(8, 5))
		d = t[t.subset.str.startswith("support_quality:", na=False)]
		for name, s in d.groupby("model", sort=False):
			if name in ("elasticnet_full", "elasticnet_clean", "twin_clean"):
				s = s.sort_values("coverage")
				ax.plot(s.coverage, s.AUC_IPCW, "o-", label=name)
		ax.set(
			xlabel="Actual held-out population coverage",
			ylabel="IPCW AUC",
			title="Frozen twin-support thresholds: same-mask comparisons",
		)
		if len(d):
			ax.legend(fontsize=8)
		save_figure(fig, out / "Fig1_twin_coverage_AUC.png")
	t = tables.get("confidence_validation", pd.DataFrame())
	if len(t):
		fig, ax = plt.subplots(figsize=(8, 5))
		for name, d in t.groupby("model", sort=False):
			ax.plot(d.mean_confidence, d.Brier_IPCW, "o-", label=name)
		ax.set(
			xlabel="Pre-outcome confidence score (0–100)",
			ylabel="Observed IPCW Brier error",
			title="Empirical validation of twin confidence, not individual certainty",
		)
		ax.legend()
		save_figure(fig, out / "Fig2_confidence_validation.png")
	t = tables.get("contrasts", pd.DataFrame())
	if len(t):
		d = t[
			(t.metric == "AUC_IPCW")
			& t.contrast.isin(
				[
					"selection_effect",
					"training_effect_same_subset",
					"combined_vs_bulk",
					"beyond_confidence_selection",
				]
			)
		]
		fig, ax = plt.subplots(figsize=(10, 5))
		y = np.arange(len(d))
		# Percentile intervals need not contain the point in finite samples.
		ax.hlines(y, d.lower, d.upper)
		ax.plot(d.delta, y, "o")
		ax.set_yticks(y, d.contrast)
		ax.axvline(0, ls="--")
		ax.set(xlabel="Paired ΔAUC", title="Selection versus training-filter effect")
		save_figure(fig, out / "Fig3_effect_decomposition.png")


def main():
	if sys.argv[1:2] == ["prepare"]:
		sys.argv.pop(1)
		return prepare_main()
	return run_main()


if __name__ == "__main__":
	# Models written from the CLI have the same importable types as library use.
	import abm

	abm.main()
