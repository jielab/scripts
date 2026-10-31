#!/usr/bin/env python3
"""C1 scientific methods, execution policy, family rosters and cohort preparation."""
from __future__ import annotations


# Scientific method definitions and execution policy
from dataclasses import dataclass, field, asdict, fields, replace
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any, Mapping, Sequence
import warnings
from analysis_layout import cohort_directory, cohort_lock_path, resolve_cohort_path

VERSION = "LE8-C1-methods-2026-10-22-v1"
INPUT_CONTRACT_VERSION = "LE8-C1-input-2026-10-22-v1"
MODEL_CONTRACT_VERSION = "LE8-C1-risk-models-2026-10-22-v1"
RELIABILITY_CONTRACT_VERSION = "LE8-C1-reliability-2026-10-26-v2"
LEGACY_MODEL_CONTRACT_VERSION = "LE8-C1-risk-models-2026-10-22-v1"


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
DEFAULT_METHODS = ("glmnet", "lightgbm", "abm_case_preserving", "abm_profile", "abm_prevalent")
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
    payload = dict(model_contract_version=MODEL_CONTRACT_VERSION, method=get_method(method).name,
                   parameters=dict(parameters), input_contract=input_contract or {})
    for key in ("device", "requested_device", "actual_device", "cores", "parallel_workers"):
        payload["parameters"].pop(key, None)
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False,
                                     default=str).encode()).hexdigest()


def method_specs():
    return [asdict(spec) for spec in METHODS.values()]


# Frozen family rosters and cohort preparation
import argparse, fcntl, hashlib, json, os, re, shutil, subprocess, tempfile
from decimal import Decimal
from pathlib import Path
import numpy as np
import pandas as pd

RELEASE="2026-10-17.c1-lockpaper-v1"
PRIMARY_SELECTOR = "risk_resid_error"
COVERAGES = (.2, .4, .5, .6, .8, 1.0)
SEED = 2026
PREDICTABILITY_TARGET_MODEL = "elasticnet_full"
DEFAULT_UKB_REL_FILE = os.getenv("LE8_UKB_REL_FILE", "/mnt/f/gen/ukb/37/typ/ukb_rel.dat")
KINSHIP_PROTOCOL = "2026-10-17.ukb-king-components-v1"
GROUPING_SOURCE = "UKB KING kinship connected components"
COHORT_CACHE_ROOT = Path("/tmp") / f"le8-cohort-{os.getuid()}"
KINSHIP_ALGORITHM = ("union-find connected components on full released relatedness graph; "
                     "then subset to all.rds; nonrelated participants assigned singleton groups")
KINSHIP_SEMANTICS = "UKB KING close-relatedness file; all released rows used as graph edges"
GROUPING_AUDIT_FIELDS = (
    "relatedness_source_file", "relatedness_source_sha256", "group_map_sha256",
    "pair_rows", "related_participant_N", "singleton_participant_N", "component_N",
    "largest_component_size", "grouping_unit", "grouping_source",
    "family_disjoint_outer_split", "family_disjoint_ABM_roles",
)


def grouping_audit(record):
    """Allowlist for public tables: never include IDs or component membership."""
    return {k: record[k] for k in GROUPING_AUDIT_FIELDS if k in record}


def _protocol_path(path):
    return Path(str(resolve_cohort_path(path)) + ".protocol.json")


def _validate_ids(values, context, unique=False):
    s = pd.Series(values, dtype="string")
    if s.empty or s.isna().any() or s.str.contains(r"\s|^$", regex=True).any():
        raise ValueError(f"Missing/invalid IDs in {context}")
    if unique and s.duplicated().any():
        raise ValueError(f"Duplicate eid in {context}")
    return s.astype(str)


def validate_related_pairs(pairs):
    required = {"ID1", "ID2", "HetHet", "IBS0", "Kinship"}
    if not required <= set(pairs):
        raise ValueError("UKB KING file requires ID1 ID2 HetHet IBS0 Kinship")
    d = pairs.copy()
    for c in ("ID1", "ID2"):
        d[c] = _validate_ids(d[c], "UKB KING source").to_numpy()
    if d.ID1.eq(d.ID2).any():
        raise ValueError("UKB KING source contains self pairs (ID1 == ID2)")
    canonical = pd.Series([tuple(sorted(p)) for p in zip(d.ID1, d.ID2)])
    if canonical.duplicated().any():
        raise ValueError("UKB KING source contains duplicated canonical pairs")
    for c in ("HetHet", "IBS0", "Kinship"):
        d[c] = pd.to_numeric(d[c], errors="raise")
        if not np.isfinite(d[c]).all():
            raise ValueError(f"UKB KING {c} must be finite")
    stats = dict(pair_rows=len(d), unique_related_ids_in_source=len(set(d.ID1) | set(d.ID2)),
                 kinship_min=float(d.Kinship.min()), kinship_median=float(d.Kinship.median()),
                 kinship_max=float(d.Kinship.max()),
                 kinship_below_0_0442_N=int((d.Kinship < .0442).sum()))
    # Conservative source-semantics guard, NEVER an edge-selection threshold.
    # Even a small unexpected tail requires source review; no rows are discarded.
    if stats["kinship_below_0_0442_N"]:
        raise ValueError("Unexpected UKB KING distribution below 0.0442; review source, "
                         "do not filter edges: " + json.dumps(stats))
    return d, stats


def build_kinship_groups(pairs, analysis_ids):
    """Build the FULL source graph before projecting onto the analytical cohort."""
    d, stats = validate_related_pairs(pairs)
    ids = _validate_ids(analysis_ids, "all.rds", unique=True)
    parent, size = {}, {}

    def find(x):
        if x not in parent:
            parent[x], size[x] = x, 1
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(d.ID1, d.ID2):
        a, b = find(a), find(b)
        if a != b:
            if (size[a], a) < (size[b], b):
                a, b = b, a
            parent[b] = a
            size[a] += size[b]
    # Numeric EIDs sort numerically; lexical fallback supports synthetic IDs.
    def key(x):
        return (0, int(x), x) if x.isdecimal() else (1, x, x)

    minima = {}
    for x in parent:
        root = find(x)
        if root not in minima or key(x) < key(minima[root]):
            minima[root] = x
    membership = {x: "KIN_" + minima[find(x)] for x in parent}
    result = pd.DataFrame(dict(eid=ids, group=[membership.get(x, "IND_" + x) for x in ids]))
    result = result.sort_values("eid").reset_index(drop=True)
    lookup = result.set_index("eid").group
    both = d.ID1.isin(lookup.index) & d.ID2.isin(lookup.index)
    if not d.loc[both, "ID1"].map(lookup).eq(d.loc[both, "ID2"].map(lookup)).all():
        raise AssertionError("Related pair assigned to different groups")
    sizes = result.groupby("group").size()
    related = result.group.str.startswith("KIN_")
    stats.update(all_rds_N=len(result), all_rds_eid_hash=eid_hash(ids),
                 all_rds_related_N=int(related.sum()), all_rds_singleton_N=int((~related).sum()),
                 component_N=len(sizes), largest_analysis_component_size=int(sizes.max()),
                 connected_components_with_at_least_2_analysis_participants=int((sizes >= 2).sum()),
                 component_size_1_N=int((sizes == 1).sum()), component_size_2_N=int((sizes == 2).sum()),
                 component_size_3_N=int((sizes == 3).sum()), component_size_4plus_N=int((sizes >= 4).sum()))
    return result, stats


def _assert_frozen_cohort(current_ids, frozen_ids):
    """Preserve every real participant and frozen role; allow only retired duplicate aliases."""
    current = set(_validate_ids(current_ids, "current all.rds", unique=True))
    frozen = set(_validate_ids(frozen_ids, "frozen protocol", unique=True))
    added, removed = current - frozen, frozen - current
    aliases = set()
    for value in removed:
        # Exact decimal conversion: float rounding must never merge participant IDs.
        if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?[eE][+-]?[0-9]{1,3}", value):
            continue
        number = Decimal(value)
        if number > 0 and number == number.to_integral_value():
            canonical = str(int(number))
            if canonical in current and canonical in frozen:
                aliases.add(value)
    if added or removed - aliases:
        raise RuntimeError("Frozen protocol cohort changed: "
                           f"added={len(added)}, removed_nonalias={len(removed - aliases)}; "
                           "review the input cohort explicitly; the frozen split was not modified")
    return dict(current_N=len(current), frozen_N=len(frozen), retired_duplicate_alias_N=len(aliases))


def _cohort_projection(allrds, expected_sha, r_bin, group_col=""):
    """Cache only IDs/group columns, bound to exact file bytes; avoid repeated large RDS loads."""
    key = hashlib.sha256(json.dumps([1, str(allrds.resolve()), expected_sha, group_col]).encode()).hexdigest()
    root = COHORT_CACHE_ROOT
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    target, meta = root / f"{key}.tsv", root / f"{key}.json"
    with open(root / f"{key}.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            record = json.loads(meta.read_text())
            valid = record == dict(source_sha256=expected_sha, group_col=group_col, projection_sha256=sha(target))
        except (OSError, ValueError):
            valid = False
        if not valid:
            with tempfile.TemporaryDirectory(dir=root) as td:
                staged = Path(td) / "cohort.tsv"
                _rds_export(allrds, staged, r_bin, group_col)
                if sha(allrds) != expected_sha:
                    raise RuntimeError("Phenotype source changed during cohort validation")
                staged.chmod(0o600)
                record = dict(source_sha256=expected_sha, group_col=group_col, projection_sha256=sha(staged))
                staged_meta = Path(td) / "cohort.json"
                staged_meta.write_text(json.dumps(record) + "\n")
                staged_meta.chmod(0o600)
                os.replace(staged, target)
                os.replace(staged_meta, meta)
        data = pd.read_csv(target, sep="\t", dtype=str)
        expected = ["eid", "group"] if group_col else ["eid"]
        if list(data) != expected:
            raise RuntimeError("Invalid cohort projection columns")
        _validate_ids(data.eid, "current all.rds", unique=True)
        return data


def _report_cohort_reuse(kind, audit):
    print(f"[LE8] SKIP {kind}: updated all.rds verified; current_N={audit['current_N']}; "
          f"frozen_N={audit['frozen_N']}; retired_duplicate_alias_N={audit['retired_duplicate_alias_N']}; "
          "original family groups and outer roles retained", flush=True)


def ensure_kinship_groups(ukb_phe, rel_file, output, r_bin="Rscript", force=False):
    """Validate graph/cohort identity; phenotype-only updates retain the frozen protocol."""
    output, rel_file = resolve_cohort_path(output), Path(rel_file)
    allrds = Path(ukb_phe) / "Rdata/all.rds"
    for source in (rel_file, allrds):
        if not source.is_file() or not source.stat().st_size:
            raise ValueError(f"Required nonempty kinship protocol source: {source}")
    output.parent.mkdir(parents=True, exist_ok=True)
    lock_path = cohort_lock_path(output)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        meta = _protocol_path(output)
        requested = dict(release=RELEASE, kinship_protocol=KINSHIP_PROTOCOL,
                         source_file=str(rel_file.resolve()), source_sha256=sha(rel_file),
                         all_rds_file=str(allrds.resolve()), all_rds_sha256=sha(allrds),
                         algorithm=KINSHIP_ALGORITHM, source_semantics=KINSHIP_SEMANTICS)
        reuse = output.is_file() and meta.is_file() and not force
        if output.exists() != meta.exists() and not force:
            raise RuntimeError("Kinship group map or protocol sidecar is missing")
        if reuse:
            old = json.loads(meta.read_text())
            if not re.fullmatch(r"[0-9a-f]{64}", str(old.get("all_rds_sha256", ""))):
                raise RuntimeError("Kinship protocol is missing its original all.rds SHA256")
            for k, v in requested.items():
                if k != "all_rds_sha256" and old.get(k) != v:
                    raise RuntimeError(f"Kinship group map protocol mismatch for {k}; archive and regenerate explicitly")
            if old.get("group_map_sha256") != sha(output):
                raise RuntimeError("Kinship group map SHA256 mismatch")
            saved = pd.read_csv(output, sep="\t", dtype=str)
            if list(saved) != ["eid", "group"]:
                raise ValueError("Kinship group map must contain only eid,group")
            if old.get("all_rds_sha256") != requested["all_rds_sha256"]:
                current = _cohort_projection(allrds, requested["all_rds_sha256"], r_bin)
                cohort_audit = _assert_frozen_cohort(current.eid, saved.eid)
            # Reconstruct the original graph projection, including provenance for
            # any retired alias. Never rewrite the frozen group map or its hash.
            ids = saved.eid
        else:
            with tempfile.TemporaryDirectory() as td:
                ph = Path(td) / "phenotype.tsv"
                _rds_export(allrds, ph, r_bin)
                ids = pd.read_csv(ph, sep="\t", dtype=str).eid
        pairs = pd.read_csv(rel_file, sep=r"\s+", dtype={"ID1": str, "ID2": str})
        groups, stats = build_kinship_groups(pairs, ids)
        record = {**requested, **stats, "grouping_unit": "family", "grouping_source": GROUPING_SOURCE}
        if reuse:
            if not saved.equals(groups):
                raise RuntimeError("Kinship group map disagrees with full source graph")
            for k, v in record.items():
                if k != "all_rds_sha256" and old.get(k) != v:
                    raise RuntimeError(f"Kinship group map audit mismatch for {k}")
        # Refuse a source that changed while the map was being built/validated.
        if sha(rel_file) != requested["source_sha256"] or sha(allrds) != requested["all_rds_sha256"]:
            raise RuntimeError("Kinship protocol source changed during validation")
        if reuse and old.get("all_rds_sha256") != requested["all_rds_sha256"]:
            _report_cohort_reuse("kinship groups", cohort_audit)
        if not reuse:
            with tempfile.TemporaryDirectory(dir=output.parent) as td:
                staged = Path(td) / output.name
                groups.to_csv(staged, sep="\t", index=False)
                staged.chmod(0o600)
                record["group_map_sha256"] = sha(staged)
                sidecar = _protocol_path(staged)
                sidecar.write_text(json.dumps(record, indent=2) + "\n")
                os.replace(staged, output)
                os.replace(sidecar, meta)
        return output


def kinship_record(group_file, *, ukb_phe=None, r_bin="Rscript"):
    """Recognize generated maps even when supplied explicitly with --group-file."""
    if not group_file:
        return {}
    meta = _protocol_path(group_file)
    if not meta.is_file():
        if Path(group_file).name == "ukb_kinship_groups.tsv":
            raise RuntimeError("Kinship group map protocol sidecar is missing")
        return {}
    record = json.loads(meta.read_text())
    if "kinship_protocol" not in record:
        if Path(group_file).name == "ukb_kinship_groups.tsv":
            raise RuntimeError("Kinship group map protocol marker is missing")
        return {}
    if record.get("kinship_protocol") != KINSHIP_PROTOCOL or record.get("release") != RELEASE:
        raise RuntimeError("Kinship group map release/protocol mismatch")
    if ukb_phe is not None:
        ensure_kinship_groups(ukb_phe, record["source_file"], group_file, r_bin=r_bin)
    if record.get("group_map_sha256") != sha(group_file) or record.get("source_sha256") != sha(record["source_file"]):
        raise RuntimeError("Kinship group map/source SHA256 mismatch")
    return record


def _kinship_audit(record):
    if not record:
        return {}
    return dict(relatedness_source_file=record["source_file"], relatedness_source_sha256=record["source_sha256"],
                group_map_sha256=record["group_map_sha256"], pair_rows=record["pair_rows"],
                related_participant_N=record["all_rds_related_N"], singleton_participant_N=record["all_rds_singleton_N"],
                component_N=record["component_N"], largest_component_size=record["largest_analysis_component_size"],
                grouping_unit="family", grouping_source=GROUPING_SOURCE)


def validate_seed(seed):
    if int(seed) != SEED or float(seed) != SEED:
        raise ValueError("C1 LOCKPAPER fixes the analysis seed at 2026; "
                         "do not change the outer split after inspecting prior results.")


def eid_hash(ids):
    """Language-independent SHA256 of sorted unique IDs, one per line."""
    return hashlib.sha256("".join(x + "\n" for x in sorted(set(map(str, ids)))).encode()).hexdigest()


def roster_record(path):
    path = resolve_cohort_path(path)
    meta = path.with_suffix(path.suffix + ".protocol.json")
    record = json.loads(meta.read_text())
    validate_seed(record.get("seed", -1))
    if record.get("release") != RELEASE or record.get("roster_sha256") != sha(path):
        raise RuntimeError("Frozen outer roster integrity/protocol mismatch: " + str(path))
    if record.get("group_file") and record.get("group_file_sha256") != sha(record["group_file"]):
        raise RuntimeError("Frozen outer roster group-file SHA256 mismatch")
    if record.get("relatedness_source_file"):
        current = kinship_record(record["group_file"])
        if not current or any(record.get(k) != v for k, v in _kinship_audit(current).items()):
            raise RuntimeError("Frozen outer roster kinship protocol mismatch")
    return record


def configure_roster(args, *, settings=None, dry_run=False):
    """Resolve shared CLI/environment settings before inspecting omics or Y."""
    validate_seed(args.seed)
    if settings is None:
        from le8_common import shared_analysis_settings
        settings, _, _ = shared_analysis_settings(args)
    for key, value in settings.items():
        setattr(args, key, value)
    validate_seed(args.seed)
    # Native result transactions use a scratch analysis root. Protocol inputs
    # always live in the published root and are never packaged as result tables.
    root = Path(os.getenv("LE8_PUBLISHED_ROOT") or args.analysis_root)
    roster = resolve_cohort_path(args.outer_roster) if args.outer_roster else (
        cohort_directory(root) / f"c1_outer_roster_seed{args.seed}.tsv"
    )
    if args.group_file:
        args.group_file = str(resolve_cohort_path(args.group_file))
    phe = args.ukb_phe or os.getenv("UKB_PHE", "/mnt/d/data/ukb")
    if not args.group_file and not args.group_col:
        source = Path(getattr(args, "ukb_rel_file", None) or os.getenv("LE8_UKB_REL_FILE", DEFAULT_UKB_REL_FILE))
        # stat() propagates permission/I/O errors; only genuine absence permits fallback.
        try:
            source.stat()
        except FileNotFoundError:
            if source.is_symlink():
                raise RuntimeError(f"Broken UKB KING source symlink: {source}")
            print(f"[LE8 C1 LOCK] WARNING: UKB KING source does not exist: {source}; "
                  "grouping_unit=individual_unverified; close-relative leakage is unverified", flush=True)
        else:
            group_map = cohort_directory(root) / "ukb_kinship_groups.tsv"
            if not dry_run:
                ensure_kinship_groups(phe, source, group_map, r_bin=args.r_bin)
            else:
                print(f"[LE8 C1 LOCK] planned UKB KING groups: {source} -> {group_map}", flush=True)
            args.group_file, args.group_col = str(group_map), ""
    # Native C1 and ABM children must inherit exactly this mapping.
    os.environ["LE8_GROUP_FILE"] = str(args.group_file or "")
    os.environ["LE8_GROUP_COLUMN"] = str(args.group_col or "")
    if not dry_run:
        ensure_roster(ukb_phe=phe,
                      output=roster, seed=args.seed, group_file=args.group_file,
                      group_col=args.group_col, r_bin=args.r_bin)
        print(f"[LE8 C1 LOCK] shared outer roster: {roster}; SHA256={sha(roster)}", flush=True)
        record = roster_record(roster)
        unit = "family" if record["group_file"] or record["group_col"] else "individual_unverified"
        print(f"[LE8 C1 LOCK] release={RELEASE}; seed={record['seed']}; "
              f"outer_test_N={record['test_N']}; test_eid_hash={record['test_eid_hash']}; "
              f"grouping_unit={unit}; independent_family_support_verified={unit == 'family'}", flush=True)
        if record.get("relatedness_source_file"):
            print("C1 LOCKPAPER relatedness:\n"
                  f"  source: {record['relatedness_source_file']}\n"
                  f"  source pairs: {record['pair_rows']}\n"
                  f"  source SHA256: {record['relatedness_source_sha256']}\n"
                  f"  grouping: {record['grouping_source']}\n"
                  f"  analysis N: {record['N']}\n"
                  f"  related analysis participants: {record['related_participant_N']}\n"
                  f"  singleton analysis participants: {record['singleton_participant_N']}\n"
                  f"  connected components (including singletons): {record['component_N']}\n"
                  f"  largest component: {record['largest_component_size']}\n"
                  f"  group-map SHA256: {record['group_map_sha256']}\n"
                  f"  seed: {record['seed']}\n"
                  f"  outer test fraction: {record['test_fraction']:.2f}\n"
                  "  family-disjoint split: TRUE\n  bootstrap unit: family", flush=True)
    else:
        print(f"[LE8 C1 LOCK] planned shared outer roster: {roster}", flush=True)
    args.outer_roster = str(roster)
    os.environ["LE8_OUTER_ROSTER"] = str(roster)
    return roster

def sha(path, block=1<<20):
    path = resolve_cohort_path(path)
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(block),b""): h.update(b)
    return h.hexdigest()

def stable_hash(x,seed):
    return int.from_bytes(hashlib.sha256(f"{seed}:{x}".encode()).digest()[:8],"big")

def _rds_export(rds,out,r_bin,group_col=""):
    code = r'''
args<-commandArgs(trailingOnly=TRUE)
x<-readRDS(args[1])
if(!is.data.frame(x)||!"eid"%in%names(x))stop("all.rds must be a data.frame with eid")
g<-args[3]
z<-data.frame(eid=as.character(x$eid),stringsAsFactors=FALSE)
if(nzchar(g)){
 if(!g%in%names(x))stop("Group column missing: ",g)
 z$group<-as.character(x[[g]])
}
data.table::fwrite(z,args[2],sep="\t",quote=FALSE,na="NA")
'''
    subprocess.run([r_bin,"--vanilla","-e",code,str(rds),str(out),group_col],check=True)

def _read_group_file(path,r_bin):
    path=resolve_cohort_path(path)
    if path.suffix.lower()==".rds":
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/"group.tsv"
            code = r'''
args<-commandArgs(trailingOnly=TRUE)
x<-readRDS(args[1]);if(!is.data.frame(x))stop("group RDS must be data.frame")
if(!all(c("eid","group")%in%names(x)))stop("group file requires eid,group")
data.table::fwrite(data.frame(eid=as.character(x$eid),group=as.character(x$group)),
                   args[2],sep="\t",quote=FALSE,na="NA")
'''
            subprocess.run([r_bin,"--vanilla","-e",code,str(path),str(out)],check=True)
            return pd.read_csv(out,sep="\t",dtype=str)
    sep="\t" if path.suffix.lower() in {".tsv",".txt"} else ","
    return pd.read_csv(path,sep=sep,dtype=str)

def ensure_roster(*,ukb_phe,output,seed=2026,test_fraction=.20,group_file="",group_col="",r_bin="Rscript",force=False):
    validate_seed(seed)
    output = resolve_cohort_path(output)
    group_file = str(resolve_cohort_path(group_file)) if group_file else ""
    output.parent.mkdir(parents=True, exist_ok=True)
    lock_path = cohort_lock_path(output)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _ensure_roster(ukb_phe=ukb_phe, output=output, seed=seed,
                              test_fraction=test_fraction, group_file=group_file,
                              group_col=group_col, r_bin=r_bin, force=force)


def _ensure_roster(*,ukb_phe,output,seed,test_fraction,group_file,group_col,r_bin,force):
    if test_fraction != .20:
        raise ValueError("C1 LOCKPAPER fixes the outer test fraction at 0.20")
    ukb_phe=Path(ukb_phe); output=Path(output)
    allrds=ukb_phe/"Rdata"/"all.rds"
    if not allrds.is_file(): raise FileNotFoundError(allrds)
    kin = kinship_record(group_file, ukb_phe=ukb_phe, r_bin=r_bin)
    meta=output.with_suffix(output.suffix+".protocol.json")
    requested=dict(release=RELEASE,seed=int(seed),test_fraction=float(test_fraction),
                   source=str(allrds.resolve()),source_sha256=sha(allrds),
                   group_file=str(Path(group_file).resolve()) if group_file else "",
                   group_file_sha256=sha(group_file) if group_file else "",
                   group_col=str(group_col or ""))
    requested.update(_kinship_audit(kin))
    # --replace (model fitting) and --force must never rerandomize a frozen roster.
    if force and (output.exists() or meta.exists()):
        raise RuntimeError("Frozen outer roster cannot be overwritten; archive the legacy individual protocol explicitly")
    if output.exists() != meta.exists() and not force:
        raise RuntimeError("Frozen outer roster or its protocol sidecar is missing")
    if output.is_file() and meta.is_file() and not force:
        old=roster_record(output)
        if not re.fullmatch(r"[0-9a-f]{64}", str(old.get("source_sha256", ""))):
            raise RuntimeError("Frozen outer roster is missing its original source SHA256")
        for k,v in requested.items():
            if k == "group_file" and old.get(k) and v and resolve_cohort_path(old[k]).resolve() == Path(v).resolve():
                continue
            if k != "source_sha256" and old.get(k)!=v:
                raise RuntimeError(f"Frozen outer roster protocol mismatch for {k}; review protocol inputs explicitly")
        if old.get("source_sha256") != requested["source_sha256"]:
            # A phenotype-defined group needs an independently frozen map to
            # prove that its memberships have not changed with all.rds.
            if group_col and not group_file:
                raise RuntimeError("Changed all.rds with a phenotype group column requires explicit grouping review")
            current = _cohort_projection(allrds, requested["source_sha256"], r_bin, group_col)
            frozen = pd.read_csv(output, sep="\t", dtype=str)
            audit = _assert_frozen_cohort(current.eid, frozen.eid)
            if group_col:
                groups = _read_group_file(group_file, r_bin)
                supplied = current.eid.map(groups.set_index("eid").group)
                if current.group.isna().any() or not current.group.eq(supplied).all():
                    raise RuntimeError("Phenotype family definitions changed from the frozen group map")
        if group_file:
            _assert_outer_groups(output, _read_group_file(group_file, r_bin))
        if sha(allrds) != requested["source_sha256"] or (group_file and sha(group_file) != requested["group_file_sha256"]):
            raise RuntimeError("Frozen roster source changed during validation")
        if old.get("source_sha256") != requested["source_sha256"]:
            _report_cohort_reuse("outer roster", audit)
        return output
    output.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as td:
        td=Path(td); ph=td/"phenotype.tsv"
        _rds_export(allrds,ph,r_bin,group_col)
        d=pd.read_csv(ph,sep="\t",dtype=str)
        if d.eid.isna().any() or d.eid.duplicated().any() or d.eid.str.strip().eq("").any(): raise ValueError("Invalid/duplicate eid in all.rds")
        if group_file:
            g=_read_group_file(group_file,r_bin)
            if not {"eid","group"}<=set(g): raise ValueError("Group file requires eid,group")
            if g.eid.isna().any() or g.eid.duplicated().any() or g.group.isna().any(): raise ValueError("Invalid group file")
            if "group" in d:
                supplied=d.eid.map(g.set_index("eid").group)
                if not d.group.eq(supplied).all(): raise ValueError("Conflicting family definitions in file and phenotype column")
                d=d.drop(columns="group")
            d=d.merge(g[["eid","group"]],on="eid",how="left",validate="one_to_one")
            if d.group.isna().any(): raise ValueError("Group file does not cover all phenotype IDs")
        elif "group" not in d:
            d["group"]=d.eid
        if d.group.isna().any() or d.group.astype(str).str.strip().eq("").any(): raise ValueError("Incomplete group mapping")
        groups=sorted(d.group.astype(str).unique(),key=lambda x:(stable_hash(x,seed),x))
        if len(groups)<2: raise ValueError("At least two independent groups are required")
        ntest=max(1,min(len(groups)-1,round(len(groups)*test_fraction)))
        test=set(groups[:ntest])
        d["role"]=d.group.astype(str).map(lambda g:"test" if g in test else "training")
        if d.groupby("group").role.nunique().max()!=1: raise AssertionError("Family crosses outer roles")
        outtmp=td/"outer.tsv"; d.sort_values("eid")[["eid","role"]].to_csv(outtmp,sep="\t",index=False); os.replace(outtmp,output)
    record={**requested,"N":int(len(d)),"groups":int(d.group.nunique()),
            "test_N":int((d.role=="test").sum()),"training_N":int((d.role=="training").sum()),
            "test_groups":int(len(test)),"roster_sha256":sha(output),
            "group_definition":group_col or ("external eid,group mapping" if group_file else "individual ID; family mapping unavailable"),
            "grouping_unit":"family" if group_file or group_col else "individual_unverified",
            "grouping_source":GROUPING_SOURCE if kin else (group_col or ("supplied eid,group mapping" if group_file else "individual IDs; kinship unverified")),
            "family_disjoint_outer_split":bool(group_file or group_col),
            "independent_family_support_verified":bool(group_file or group_col),
            "test_eid_hash":eid_hash(d.loc[d.role=="test","eid"]),
            "rule":"sort groups by SHA256(seed:group); exact prespecified group count assigned to test; outcome/omics blind"}
    tmp=Path(str(meta)+".tmp");tmp.write_text(json.dumps(record,indent=2),encoding="utf-8");os.replace(tmp,meta)
    return output


def _assert_outer_groups(roster, groups):
    d = pd.read_csv(roster, sep="\t", dtype=str)
    _validate_ids(d.eid, "outer roster", unique=True)
    _validate_ids(groups.eid, "group map", unique=True)
    if groups.group.isna().any() or groups.group.astype(str).str.strip().eq("").any():
        raise ValueError("Incomplete group mapping")
    if set(d.eid) - set(groups.eid) or set(d.role) != {"training", "test"}:
        raise ValueError("Invalid/incomplete outer roster")
    d = d.merge(groups[["eid", "group"]], on="eid", validate="one_to_one")
    if not d.groupby("group").role.nunique().eq(1).all():
        raise ValueError("Kinship component crosses outer roles")


def archive_pre_kinship_roster(output):
    """Explicit one-time migration only; never archive an already grouped roster."""
    output = resolve_cohort_path(output)
    # Native transactions hold this shared analysis lock while fitting. Refuse
    # migration under a live run, whose already-loaded code/settings cannot upgrade.
    root = output.parents[2] if output.parent.name == "freeze" and output.parent.parent.name == "logs" else output.parent.parent
    root_key = hashlib.sha256(str(root.resolve()).encode()).hexdigest()[:20]
    lock_path = cohort_lock_path(output)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(f"/tmp/le8-{root_key}.lock", "a") as analysis_lock, open(lock_path, "a") as lock:
        try:
            fcntl.flock(analysis_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("C1 analysis is running; stop the old run before archiving its outer roster") from exc
        fcntl.flock(lock, fcntl.LOCK_EX)
        record = roster_record(output)
        # The original LOCKPAPER sidecar predates grouping_unit. Missing that
        # field is only legacy evidence when every grouping field agrees.
        legacy_individual = (
            record.get("grouping_unit") in (None, "", "individual_unverified")
            and not record.get("group_file")
            and not record.get("group_col")
            and not record.get("relatedness_source_file")
            and (record.get("group_definition") or "") in (
                "", "individual ID; family mapping unavailable",
            )
            and (record.get("grouping_source") or "") in (
                "", "individual IDs; kinship unverified",
            )
            and not any(record.get(key) for key in (
                "kinship_protocol", "group_map_sha256",
                "family_disjoint_outer_split", "family_disjoint_ABM_roles",
                "independent_family_support_verified",
            ))
        )
        if not legacy_individual:
            raise RuntimeError("Only a verified legacy individual/unverified roster can be archived by this migration")
        history = output.parent / "_history/pre_kinship_20261017"
        history.mkdir(parents=True, exist_ok=True)
        paths = (output, _protocol_path(output))
        if any((history / p.name).exists() for p in paths):
            raise FileExistsError("Pre-kinship archive already exists; do not repeat the protocol upgrade")
        # Copy and verify both before removing either original, retaining recovery
        # data even if an interrupted migration needs manual repair.
        for p in paths:
            shutil.copy2(p, history / p.name)
            if sha(p) != sha(history / p.name):
                raise RuntimeError("Pre-kinship archive verification failed")
        for p in paths:
            p.unlink()
        return history

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument("--ukb-phe",required=True);p.add_argument("--output",required=True)
    p.add_argument("--seed",type=int,default=2026);p.add_argument("--test-fraction",type=float,default=.20)
    p.add_argument("--group-file",default="");p.add_argument("--group-col",default="");p.add_argument("--r-bin",default="Rscript")
    p.add_argument("--ukb-rel-file", default=os.getenv("LE8_UKB_REL_FILE", DEFAULT_UKB_REL_FILE))
    p.add_argument("--archive-pre-kinship", action="store_true",
                   help="One-time upgrade: archive only an existing individual_unverified roster")
    a=p.parse_args(argv)
    validate_seed(a.seed)
    if a.test_fraction != .20:
        p.error("C1 LOCKPAPER fixes the outer test fraction at 0.20")
    output = resolve_cohort_path(a.output)
    a.analysis_root = output.parents[2] if output.parent.name == "freeze" and output.parent.parent.name == "logs" else output.parent.parent
    a.outer_roster = a.output
    if a.archive_pre_kinship:
        if a.group_file or a.group_col:
            p.error("Pre-kinship migration uses the UKB KING source")
        # Validate/build the new source before touching the previous frozen split.
        ensure_kinship_groups(a.ukb_phe, a.ukb_rel_file,
                              cohort_directory(a.analysis_root) / "ukb_kinship_groups.tsv", r_bin=a.r_bin)
        print("Archived legacy protocol:", archive_pre_kinship_roster(a.output), flush=True)
    configure_roster(a, settings={})


if __name__ == "__main__":
    import sys
    sys.exit(main())
