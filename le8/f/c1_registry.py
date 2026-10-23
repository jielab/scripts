"""Scientific C1 methods and an independent, explicit execution-device policy.

This module has no model-fitting imports at import time. Plans and --help therefore
work on machines without R, LightGBM, PyTorch or CUDA. A method is never renamed
because it ran on another device. Parameters, input contracts and code versions
belong in scientific cache keys; the resolved device is separate provenance.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict, fields, replace
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any, Mapping, Sequence
import warnings

VERSION = "LE8-C1-methods-2026-10-22-v1"


@dataclass(frozen=True)
class MethodSpec:
    name: str
    engine: str
    description: str
    abm_overrides: Mapping[str, Any] = field(default_factory=dict)
    experimental: bool = False
    uses_prevalent: bool = False
    cuda_supported: bool = False


_PRESERVE = dict(pruning="case_preserving", keep_fraction=1.0,
                 reference_radius_policy="full_fixed", min_case_retention=1.0,
                 panel_selection="pilot", geometry="coefficient_panel",
                 technical_qc_scope="assay", lightgbm=False, primary="twin_clean",
                 final_multiroute_release="", reference_experiment="unified_case_preserving_20261022",
                 radius_k_ref=64)
METHODS = {
    "glmnet": MethodSpec("glmnet", "r_glmnet", "Penalized logistic incidence model on the shared horizon/roster"),
    "lightgbm": MethodSpec("lightgbm", "lightgbm", "Tree incidence model on the shared horizon/roster"),
    "abm_legacy": MethodSpec("abm_legacy", "twin", "Historical loss-based donor cleaning; explicit sensitivity only",
        dict(pruning="profile_concordance", keep_fraction=.95,
             reference_radius_policy="arm_calibrated", min_case_retention=0.0,
             min_retained_cases=0, geometry="coefficient_panel", panel_selection="pilot",
             technical_qc_scope="legacy_panel", lightgbm=False, primary="elasticnet_clean",
             radius_k_ref=0), True),
    "abm_case_preserving": MethodSpec("abm_case_preserving", "twin",
        "Retain every technically eligible incident case and every control", _PRESERVE),
    "abm_lenient": MethodSpec("abm_lenient", "twin",
        "Retain all incident cases; remove at most 2% of eligible controls by development loss",
        dict(_PRESERVE, keep_fraction=.98), True),
    "abm_profile": MethodSpec("abm_profile", "twin",
        "All molecular features in a build-fitted PCA matching space, retaining cases",
        dict(_PRESERVE, panel_selection="all", geometry="all_molecular_pca")),
    "abm_prevalent": MethodSpec("abm_prevalent", "prevalent_transfer",
        "Prevalent-case similarity learned in build families; incidence fit/calibration remain incident-only",
        dict(_PRESERVE, panel_selection="all", geometry="all_molecular_pca"), True, True),
    "abm_deep": MethodSpec("abm_deep", "torch_transformer",
        "Feature-token Transformer and incident-reference support; same model on CPU or CUDA",
        dict(_PRESERVE, panel_selection="all", geometry="all_molecular_pca"), True, False, True),
    "molecular_deviation": MethodSpec("molecular_deviation", "deviation",
        "Existing molecular deviation analysis; optional within C1", experimental=True),
}
DEFAULT_METHODS = ("glmnet", "lightgbm", "abm_case_preserving", "abm_profile", "abm_prevalent",
                   "abm_legacy", "molecular_deviation", "abm_deep")
METHOD_ALIASES = {"pradeep": "glmnet", "lgbm": "lightgbm", "abm_multiview": "abm_profile"}
LEGACY_MODULE_METHODS = {"c1_abm": ("abm_legacy", "molecular_deviation"), "c1_abm_GPU": ("abm_deep",)}


def get_method(name: str) -> MethodSpec:
    canonical = METHOD_ALIASES.get(str(name).strip(), str(name).strip())
    if canonical not in METHODS:
        raise ValueError(f"Unknown C1 method {name!r}; choose from {', '.join(METHODS)}")
    return METHODS[canonical]


def parse_methods(value: str | Sequence[str] | None = None, *, skip_abm: bool = False) -> tuple[str, ...]:
    raw = DEFAULT_METHODS if value is None else value.split(",") if isinstance(value, str) else value
    if not raw or any(not str(x).strip() for x in raw):
        raise ValueError("A nonempty C1 method list is required")
    result = tuple(dict.fromkeys(get_method(x).name for x in raw))
    if skip_abm:
        result = tuple(x for x in result if not x.startswith("abm_") and x != "molecular_deviation")
    if not result:
        raise ValueError("No C1 methods remain after --skip-abm")
    return result


def normalize_modules(modules: Sequence[str], methods=None, *, skip_abm=False, run_abm=False):
    """Collapse old module aliases once, without duplicating C1 execution."""
    requested = list(dict.fromkeys(modules))
    old = [m for m in requested if m in LEGACY_MODULE_METHODS]
    has_c1 = "c1_correlate" in requested or bool(old) or run_abm
    if methods is not None and not has_c1:
        raise ValueError("--methods applies to c1_correlate; add that module")
    if old:
        warnings.warn("Deprecated C1 module alias: " + ", ".join(old) +
                      "; use c1_correlate --methods METHOD[,METHOD]", FutureWarning, stacklevel=2)
    if methods is None and "c1_correlate" not in requested and old:
        methods = [method for module in old for method in LEGACY_MODULE_METHODS[module]]
    elif methods is None and "c1_correlate" in requested and old:
        methods = list(DEFAULT_METHODS) + [method for module in old for method in LEGACY_MODULE_METHODS[module]]
    canonical = [m for m in requested if m not in LEGACY_MODULE_METHODS and m != "c1_correlate"]
    if has_c1:
        canonical.insert(0, "c1_correlate")
    return canonical, parse_methods(methods, skip_abm=skip_abm) if has_c1 else ()


def _torch_cuda_available() -> bool:
    if importlib.util.find_spec("torch") is None:
        return False
    import torch
    return bool(torch.cuda.is_available())


def resolve_device(method: str | MethodSpec, requested="auto", *, check_available=True):
    spec = get_method(method) if isinstance(method, str) else method
    if requested not in {"auto", "cpu", "cuda"}:
        raise ValueError("--device requires auto, cpu or cuda")
    if not spec.cuda_supported:
        return dict(requested=requested, actual="cpu", backend=spec.engine,
                    reason="This engine executes on CPU; --device controls acceleration-capable methods")
    if requested == "cpu":
        return dict(requested=requested, actual="cpu", backend=spec.engine, reason="Explicit CPU execution")
    if not check_available:
        return dict(requested=requested, actual="cuda" if requested == "cuda" else "auto",
                    backend=spec.engine, reason="Device capability will be checked before fitting")
    available = _torch_cuda_available()
    if requested == "cuda" and not available:
        raise RuntimeError(f"{spec.name}: --device cuda requested, but PyTorch CUDA is unavailable")
    actual = "cuda" if available else "cpu"
    return dict(requested=requested, actual=actual, backend=spec.engine,
                reason="CUDA available" if available else "CUDA unavailable; auto uses the same model on CPU")


def build_abm_config(method, base_config=None, *, overrides=None, **identity):
    """Apply a scientific preset, explicit per-method overrides, then cohort identity.

    Named case-preserving methods cannot silently become case-deleting methods.
    Legacy raw-loss selection must remain explicitly named abm_legacy.
    """
    import abm
    spec = get_method(method)
    if base_config is None:
        cfg = abm.Config()
    elif isinstance(base_config, (str, Path)):
        cfg = abm.load_config(base_config)
    elif isinstance(base_config, Mapping):
        names = {f.name for f in fields(abm.Config)}
        unknown = set(base_config) - names
        if unknown:
            raise ValueError("Unknown ABM parameter(s): " + ", ".join(sorted(unknown)))
        cfg = abm.Config(**base_config)
    else:
        cfg = base_config
    params = dict(spec.abm_overrides)
    params.update(overrides or {})
    params.update(identity)
    if spec.name.startswith("abm_") and spec.name != "abm_legacy":
        if params.get("pruning") not in {"case_preserving", "none"}:
            raise ValueError("This method must preserve cases; select abm_legacy for loss-based deletion")
        if float(params.get("min_case_retention", 1)) != 1:
            raise ValueError("Case-preserving methods require min_case_retention=1")
    cfg = replace(cfg, **params)
    cfg.validate()
    return cfg


def load_method_plan(path=None, *, methods=None):
    """Read strict JSON {methods: [...], overrides: {method: {parameter: value}}}."""
    if path is None:
        return dict(version=VERSION, methods=list(parse_methods(methods)), overrides={})
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict) or set(data) - {"version", "methods", "overrides", "description"}:
        raise ValueError("Method plan accepts only version, methods, overrides and description")
    chosen = parse_methods(methods if methods is not None else data.get("methods"))
    overrides = data.get("overrides", {})
    if not isinstance(overrides, dict):
        raise ValueError("overrides must map method IDs to parameter dictionaries")
    normalized = {}
    for name, params in overrides.items():
        name = get_method(name).name
        if name in normalized or not isinstance(params, dict):
            raise ValueError("Duplicate method alias or invalid override dictionary")
        if {"device", "cuda", "gpu"} & set(params):
            raise ValueError("Execution device belongs in --device, not scientific parameters")
        normalized[name] = params
    return dict(version=VERSION, methods=list(chosen), overrides=normalized)


def scientific_fingerprint(method, parameters, input_contract=None):
    payload = dict(registry_version=VERSION, method=get_method(method).name,
                   parameters=dict(parameters), input_contract=input_contract or {})
    for key in ("device", "requested_device", "actual_device", "cores", "parallel_workers"):
        payload["parameters"].pop(key, None)
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False,
                                     default=str).encode()).hexdigest()


def method_specs():
    return [asdict(spec) for spec in METHODS.values()]
