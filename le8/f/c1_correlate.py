#!/usr/bin/env python3
"""Unified C1: association plus named prediction methods on one frozen cohort.

All selected prediction methods finish fit and outcome-blind projection before
test truth is opened by the common evaluator. Current models and evaluation
live under the canonical c1_correlate scope.
"""
from __future__ import annotations

from file_lock import flock as shared_flock

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
# Keep project modules (notably abm.py) ahead of shared helpers.
_paths_sys.path.append(str(_paths_root / "0f"))
from workflow_paths import configured_path
from workflow_jobs import atomic_json
import sys
sys.dont_write_bytecode = True

import argparse
import copy
from contextlib import nullcontext
from dataclasses import asdict, replace
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import fcntl
import hashlib
import importlib.util
import json
import os
import re
import shlex
import shutil
import subprocess
import tempfile
import traceback

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import abm
from c1_protocol import (VERSION as REGISTRY_VERSION, get_method, load_method_plan,
                        parse_methods, build_abm_config, resolve_device, INPUT_CONTRACT_VERSION,
                        MODEL_CONTRACT_VERSION)
from c1_models import (InsufficientDevelopmentSupport, check_development_support, context, fit_glmnet,
                               fit_lightgbm, fit_twin, fit_prevalent, fit_deep)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
VERSION = "LE8-C1-FINAL-2026-11-02-count"
BASELINE_COMMIT = "40a884324f30219f865e8a7a49ac9708a913e742"
INSUFFICIENT_SUPPORT = "skipped_insufficient_development_support"
INSUFFICIENT_REFERENCE = "skipped_insufficient_reference_support"
PARTIAL_COMPLETION = "completed_with_unavailable_methods"
PARTIAL_RELIABILITY = "completed_with_unavailable_reliability"
INSUFFICIENT_SCOPE_SUPPORT = "skipped_insufficient_scope_support"
PREPARED_FILES = ("development_LOCAL_ONLY.npz", "test_X_LOCAL_ONLY.npz",
                  "test_truth_LOCAL_ONLY.npz", "host_contract.json")
OPTIONAL_PREPARED_FILES = ("prevalent_source_LOCAL_ONLY.joblib", "deviation_LOCAL_ONLY.joblib")


def runtime_source_commit():
    if os.getenv("LE8_SOURCE_COMMIT"):
        return os.environ["LE8_SOURCE_COMMIT"]
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def replace_stage(args, stage):
    value = getattr(args, "replace_" + stage, None)
    return (args.replace if value is None else value) == "TRUE"


def load_module(name, path):
    if name in sys.modules:
        return sys.modules[name]
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec); sys.modules[name] = module
    try: spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(name, None); raise
    return module


def file_size(path):
    """A lightweight output inventory; scientific validity is checked on load."""
    return Path(path).stat().st_size


def stage_directory():
    """A unique plan/scope handoff; never infer the current run from old results."""
    value = os.getenv("LE8_C1_BATCH_DIR")
    return Path(value) if value else None


def prepared_input_contract(source, trait, layer):
    """Check the small completion receipt before accepting a prepared bundle.

    No test outcomes or participant arrays are read by this integrity check.
    Scientific IDs, roles, shapes and values are validated by the model loader.
    """
    source = Path(source)
    marker = json.loads((source / "INPUTS_COMPLETE.json").read_text())
    if not isinstance(marker, dict) or marker.get("status") != "complete":
        raise ValueError("Prepared inputs lack a completed receipt: " + str(source))
    files = marker.get("files")
    if (not isinstance(files, dict) or not set(PREPARED_FILES) <= files.keys()
            or files.keys() - set(PREPARED_FILES + OPTIONAL_PREPARED_FILES)):
        raise ValueError("Invalid prepared input inventory: " + str(source))
    for name, size in files.items():
        path = source / name
        if (not isinstance(size, int) or isinstance(size, bool) or size < 1
                or not path.is_file() or path.stat().st_size != size):
            raise ValueError("Missing or incomplete prepared input: " + str(path))
    contract = json.loads((source / "host_contract.json").read_text())
    if not isinstance(contract, dict):
        raise ValueError("Prepared host contract must be a JSON object")
    identity = contract.get("config", contract.get("input_config", {}))
    if not isinstance(identity, dict):
        raise ValueError("Prepared host input configuration must be a JSON object")
    recorded_trait = contract.get("trait", identity.get("trait"))
    recorded_layer = contract.get("layer", identity.get("layer", identity.get("biom")))
    if recorded_trait != trait or recorded_layer != layer:
        raise ValueError("Prepared inputs belong to a different or unspecified trait/layer")
    return contract, files


def freeze_prepared_inputs(args, source):
    directory = stage_directory()
    if directory is None:
        return
    _, files = prepared_input_contract(source, args.Y, args.biom)
    atomic_json(directory / "prepared.json", dict(trait=args.Y, layer=args.biom,
        source=str(source.resolve()), files=files))


def bind_prepared_inputs(args):
    directory = stage_directory()
    if directory is None:
        return
    saved = json.loads((directory / "prepared.json").read_text())
    if saved.get("trait") != args.Y or saved.get("layer") != args.biom:
        raise ValueError("Stage input receipt belongs to a different trait/layer")
    source = Path(saved["source"])
    if args.prepared_source and args.prepared_source.resolve() != source.resolve():
        raise ValueError("Prepared source differs from this plan's frozen input receipt")
    _, files = prepared_input_contract(source, args.Y, args.biom)
    if saved.get("files") != files:
        raise RuntimeError("Prepared inputs changed after this plan's preparation stage")
    args.prepared_source = source


def frozen_stage_models(args, methods, source):
    """Read exact primary/deep receipts before the CPU evaluator opens test truth."""
    directory = stage_directory()
    if directory is None or args.stage != "evaluation":
        return None
    parts = os.environ.get("LE8_C1_MODEL_PARTS", "").split(",")
    if not parts or any(p not in {"primary", "deep"} for p in parts) or len(set(parts)) != len(parts):
        raise ValueError("Invalid model-stage receipt list")
    result = dict(methods={}, unavailable_methods={})
    seen = set()
    for part in parts:
        saved = json.loads((directory / f"models.{part}.json").read_text())
        if saved.get("trait") != args.Y or saved.get("layer") != args.biom or saved.get("source") != source:
            raise ValueError("Model-stage receipt differs from the frozen input contract")
        if saved.get("failed_methods"):
            raise RuntimeError("A model stage failed; evaluation cannot use older results")
        for category in result:
            entries = saved.get(category, {})
            if seen & entries.keys():
                raise ValueError("Duplicate method in primary/deep stage receipts")
            seen.update(entries)
            result[category].update(entries)
    allowed = set(methods) | ({"abm_prevalent_control"} if "abm_prevalent" in methods else set())
    if not set(methods) <= seen or seen - allowed:
        raise ValueError("primary/deep stage receipts do not cover the requested methods")
    return result


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument("--Y", default="cvd_cad,ra")
    p.add_argument("--stage", choices=["all", "prepare", "association", "models", "evaluation"], default="all",
                   help="Restartable internal stage; evaluation refuses to fit missing models")
    p.add_argument("--biom", default="prot,met")
    p.add_argument("--analysis-root", type=Path,
        default=Path(os.getenv("LE8_ANALYSIS_ROOT", configured_path('ANALYSIS_ROOT', 'le8'))))
    p.add_argument("--methods", help="Comma-separated method IDs; defaults to the methods in --config")
    p.add_argument("--config", type=Path, default=ROOT/"config/c1_methods.json",
                   help="Scientific method-plan JSON")
    p.add_argument("--abm-config", type=Path, default=ROOT/"config/abm.json")
    p.add_argument("--deep-config", type=Path, default=ROOT/"config/abm_deep.json")
    p.set_defaults(history_config=None)
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="cpu")
    p.add_argument("--association", choices=["auto", "run", "skip"], default="run")
    p.add_argument("--pgs", choices=["TRUE", "FALSE"],
                   type=lambda value: {"RUN": "TRUE", "SKIP": "FALSE"}.get(value.upper(), value.upper()),
                   default=os.getenv("C1_PGS_MODE", "FALSE"),
                   help="Enable molecular PGS association scan; default FALSE (not required by primary ABM)")
    p.add_argument("--prepared-source", type=Path, help="Reuse a frozen prepared input directory (one trait/layer only)")
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--min-cases-with-biom", type=int,
                   default=int(os.getenv("C1_MIN_CASES_WITH_BIOM", "50")),
                   help="Minimum known development cases for prediction only; association remains independent [50]")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--cores", type=int, default=8)
    p.add_argument("--replace", choices=["TRUE", "FALSE"], default="TRUE")
    for stage in ("association", "models", "evaluation", "reliability"):
        p.add_argument("--replace-" + stage, choices=["TRUE", "FALSE"], default=None,
                       help="Override --replace for this C1 stage only")
    p.add_argument("--reliability-protocol", type=Path, default=ROOT/"config/c1_reliability_protocol.json")
    p.add_argument("--r-bin", default="Rscript")
    p.add_argument("--host-args", default="")
    p.add_argument("--preflight", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    for flag in ("ukb-phe", "ukb-rel-file", "group-file", "group-col", "outer-roster",
                 "diagnosis-col", "shared-covariates", "deviation-config", "objective-file", "pgs-mapping"):
        p.add_argument("--"+flag)
    return p


def plan_and_base(args):
    base = json.loads(args.abm_config.read_text())
    config = args.config
    plan = load_method_plan(config, methods=args.methods)
    base["seed"] = args.seed
    if args.bootstrap < 1 or args.cores < 1 or args.min_cases_with_biom < 1:
        raise ValueError("Positive bootstrap, core and minimum-case counts are required")
    return plan, base


def public_plan(args, plan, base):
    methods = []
    for name in plan["methods"]:
        spec = get_method(name)
        cfg = build_abm_config(name if spec.engine in {"twin", "twin_count", "prevalent_transfer", "torch_transformer"}
            else "abm_case_preserving", base, overrides=plan["overrides"].get(name, {}),
            seed=args.seed, cores=args.cores, parallel_workers=min(args.cores, 4),
            bootstrap=args.bootstrap, final_multiroute_release="")
        template = asdict(cfg)
        for field in ("trait", "layer", "always_include_indices"):
            template.pop(field, None)
        methods.append(dict(method=name, description=spec.description, experimental=spec.experimental,
            configuration=template, cohort_identity="trait/layer from CLI; clinical schema from prepared inputs",
            device=resolve_device(name, args.device, check_available=False)))
    return dict(version=VERSION, baseline_commit=BASELINE_COMMIT, source_commit=runtime_source_commit(), traits=args.Y.split(","),
        layers=args.biom.split(","), methods=methods, association=args.association, pgs=args.pgs,
        min_cases_with_biom=args.min_cases_with_biom,
        replace_stages={stage:replace_stage(args,stage) for stage in ("association","models","evaluation","reliability")},
        execution_order="trait, then layer: association -> per-method support -> methods -> evaluation; finish each scope before the next",
        common_evaluation="same test people, outcome horizon, build censoring model and paired family draws",
        test_provenance="previously_inspected_or_unverified_exploratory",
        proposed_outputs="<analysis-root>/<Y>/<biom>/c1_correlate/{*.png,*.xlsx,qc,tmp/models}")


def dependencies(args, plan):
    required = set(plan["methods"])
    from c1_evaluation import load_reliability_protocol,DEFAULT_COVERAGES,PRIMARY_COVERAGE,RISK_BINS
    protocol,_=load_reliability_protocol(args.reliability_protocol)
    if (protocol["coverage_grid"]!=list(DEFAULT_COVERAGES) or protocol["primary_coverage"]!=PRIMARY_COVERAGE or protocol["risk_bins"]!=RISK_BINS):
        raise ValueError("Integrated C1 requires the declared FINAL reliability design")
    if "glmnet" in required:
        subprocess.run([args.r_bin, "-e", "p<-c('glmnet','jsonlite','brglm2');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste(m,collapse=', '))"], check=True)
    if "lightgbm" in required:
        try:
            import lightgbm  # Import dependencies too, rather than only finding its top-level spec.
        except ImportError as error:
            raise RuntimeError("LightGBM or one of its dependencies is unavailable; install the LE8 environment") from error
    if "abm_deep" in required:
        if importlib.util.find_spec("torch") is None:
            raise RuntimeError("Missing PyTorch for abm_deep (CPU and CUDA use the same model)")
        resolve_device("abm_deep", args.device)
    if args.history_config:
        from c1_models import PrevalentConfig
        PrevalentConfig(**json.loads(args.history_config.read_text())).validate()


def native_association(args):
    demo = "--demo" in shlex.split(args.host_args)
    if args.association == "skip" or (args.association == "auto" and (args.prepared_source or demo)):
        return {"status": "skipped_by_scope", "existing_outputs_preserved": True}
    cmd = [sys.executable, str(HERE/"0.common.py"), "dispatch", "c1_correlate", "--Y", args.Y, "--biom", args.biom,
           "--analysis-root", str(args.analysis_root), "--r-bin", args.r_bin,
           "--cores", str(args.cores), "--seed", str(args.seed), "--replace", "TRUE" if replace_stage(args,"association") else "FALSE"]
    if args.ukb_phe: cmd += ["--ukb-phe", args.ukb_phe]
    cmd += ["--memory-limit-gb", str(int(float(os.getenv("LE8_MEMORY_LIMIT_GB", "32"))))]
    env = os.environ.copy()
    env.pop("C1_MIN_CASES_WITH_BIOM", None)  # Prediction admission must not invalidate R association caches.
    env["LE8_C1_STAGE"] = "c1_correlate" # only association; prediction handled below
    env["C1_PGS_MODE"] = {"TRUE": "run", "FALSE": "skip"}[args.pgs]
    for key, var in (("diagnosis_col", "LE8_Y_DATE"),
                     ("shared_covariates", "LE8_VARS_ADJ"), ("outer_roster", "LE8_OUTER_ROSTER"),
                     ("group_file", "LE8_GROUP_FILE"), ("group_col", "LE8_GROUP_COLUMN")):
        value = getattr(args, key, None)
        if value: env[var] = str(value)
    print("[LE8 C1] Association: " + shlex.join(cmd), flush=True)
    subprocess.run(cmd, env=env, check=True)
    return {"status": "completed_or_reused", "engine": "existing native R association", "prediction_refit": False}


def save_prepared(incident, deviation, audit, target):
    """Write outcomes apart from features, with private permissions inherited."""
    target.mkdir(parents=True, exist_ok=True)
    roles = np.asarray(incident["role"]).astype(str)
    dev, test = roles != "test", roles == "test"
    ids = np.asarray(incident["eid"]).astype(str)
    groups = np.asarray(incident["family_group"]).astype(str)
    np.savez_compressed(target/"development_LOCAL_ONLY.npz", X=incident["X"][dev],
        time=incident["time"][dev], event=incident["event"][dev], eid=ids[dev], group=groups[dev],
        role=roles[dev], feature_names=np.asarray(incident["feature_names"]))
    np.savez_compressed(target/"test_X_LOCAL_ONLY.npz", X=incident["X"][test], eid=ids[test],
        group=groups[test], feature_names=np.asarray(incident["feature_names"]))
    np.savez_compressed(target/"test_truth_LOCAL_ONLY.npz", time=incident["time"][test],
        event=incident["event"][test], eid=ids[test])
    history = deviation.get("prevalent_aux")
    if history is not None:
        order = np.r_[np.flatnonzero(dev), np.flatnonzero(test)]
        history = dict(history)
        history["incident_metadata"] = {k:np.asarray(v)[order] for k,v in history["incident_metadata"].items()}
        joblib.dump(history, target/"prevalent_source_LOCAL_ONLY.joblib", compress=3)
    abm.json_write(target/"host_contract.json", {**audit,
        "trait": audit["input_config"]["trait"], "layer": audit["input_config"]["biom"],
        "clinical_feature_names": incident["clinical_feature_names"], "numeric_schema": incident["feature_names"],
        "incident_test_ids_hash": abm.fingerprint(ids[test]), "prevalent_in_incident_risk_set": False})


from outputs import named_directory
from outputs import staging_directory


def prepare(args, plan, scope):
    # An input bundle may be requested by two independently submitted plans.
    # Lock only this scope's preparation; all other scopes remain independent.
    directory = scope.parent / "logs/locks"
    directory.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(str(scope.resolve()).encode()).hexdigest()[:20]
    with (directory / f"prepare-{key}.lock").open("a") as lock:
        shared_flock(lock, fcntl.LOCK_EX)
        return _prepare(args, plan, scope)


def _prepare(args, plan, scope):
    host = load_module("le8_unified_host", HERE/"0.host.py")
    source = args.prepared_source.resolve() if args.prepared_source else scope/"tmp/inputs/current"
    if args.prepared_source or args.replace == "FALSE":
        try:
            contract, _ = prepared_input_contract(source, args.Y, args.biom)
        except (OSError, ValueError, TypeError) as error:
            if args.prepared_source:
                raise ValueError("Incomplete explicit prepared input directory: " + str(source)) from error
            print("[LE8 C1] Rebuilding missing/incomplete prepared inputs: " + str(source), flush=True)
        else:
            print("[LE8 C1] Reusing complete prepared inputs; IDs, schema and roles are checked on load", flush=True)
            return source, contract, None
    args.include_prevalent = False
    source.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".prepare-", dir=source.parent) as directory:
        work = Path(directory)
        incident, deviation, audit = host.prepare_host(args, work)
        staged = work/"input_bundle"
        save_prepared(incident, deviation, audit, staged)
        names = PREPARED_FILES + tuple(n for n in OPTIONAL_PREPARED_FILES if (staged / n).is_file())
        abm.json_write(staged/"INPUTS_COMPLETE.json", dict(status="complete",
            files={name:file_size(staged / name) for name in names}))
        prepared_input_contract(staged, args.Y, args.biom)
        host.publish(staged, source)
    return source, json.loads((source/"host_contract.json").read_text()), None


def read_inputs(source, contract):
    d, x = abm.read_npz(source/"development_LOCAL_ONLY.npz"), abm.read_npz(source/"test_X_LOCAL_ONLY.npz")
    clinical = contract.get("clinical_feature_names")
    if not clinical or not set(clinical) <= set(d["feature_names"]):
        raise ValueError("Explicit matching clinical feature schema is required")
    names = list(map(str,d["feature_names"]))
    indices = tuple(names.index(str(v)) for v in clinical)
    history_path = source/"prevalent_source_LOCAL_ONLY.joblib"
    history = abm.load_model(history_path) if history_path.is_file() else None
    grouping = contract.get("independent_family_support_verified", contract.get("config", {}).get("grouping_unit")=="family")
    return d, x, indices, history, "family" if grouping else "individual_unverified"


def completion(path, signature=None):
    """Reuse complete current outputs; code/data checksums are not a cache gate."""
    try:
        value = json.loads((path/"COMPLETE.json").read_text())
        files = value.get("files", {})
        return (isinstance(files, dict) and bool(files) and all(
            not Path(n).is_absolute() and ".." not in Path(n).parts
            and isinstance(size, int) and size > 0
            and (path/n).is_file() and (path/n).stat().st_size == size
            for n, size in files.items()))
    except (OSError, ValueError, KeyError, TypeError):
        return False


def reusable_method(scope, name, signature, default, *, method_cfg=None, context=None, extra=None, source=None):
    return default, signature


def model_matches_cohort(path, c, method):
    """A cached fit must still address these people and this feature schema."""
    try:
        manifest = json.loads((path/"MODEL_FROZEN.json").read_text())
        if manifest.get("method") != method:
            return False
        fitted = abm.load_model(path/"model_LOCAL_ONLY.joblib")
        if list(map(str, fitted.get("feature_names", []))) != list(map(str, c.names)):
            return False
        if fitted.get("config", {}).get("horizon") != c.cfg.horizon:
            return False
        for role in ("calibration_audit", "test"):
            saved = abm.load_model(path/(role+"_predictions_LOCAL_ONLY.joblib"))
            ix = c.parts[role]
            if not np.array_equal(np.asarray(saved["eid"]).astype(str), np.asarray(c.ids[ix]).astype(str)):
                return False
            if not np.array_equal(np.asarray(saved["group"]).astype(str), np.asarray(c.groups[ix]).astype(str)):
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError, AttributeError, EOFError, ImportError):
        return False


def mark_complete(path, signature, metadata):
    # Actual neural/R weights and support donors are part of cache integrity,
    # not just a root-level pointer to them. Optional deviation has nested files.
    anchors = [p for p in path.rglob("*") if p.is_file() and
        (p.name.endswith(("LOCAL_ONLY.joblib", "LOCAL_ONLY.pt", "LOCAL_ONLY.rds")) or
         p.name in {"MODEL_FROZEN.json", "model.joblib", "table_manifest.json", "host_contract.json"})]
    if not anchors:
        raise ValueError("Cannot mark a method complete without a frozen model or verified result anchor")
    abm.json_write(path/"COMPLETE.json", dict(signature=signature, metadata=metadata,
        files={str(p.relative_to(path)):file_size(p) for p in anchors}))


def unavailable_result(path, signature):
    """Unavailable is a verified status, never a completed fitted model."""
    try:
        result = json.loads((path/"UNAVAILABLE.json").read_text())
        if (result.get("status") in {INSUFFICIENT_REFERENCE, INSUFFICIENT_SUPPORT}
                and result.get("files") and not (path/"COMPLETE.json").exists()
                and all((path/n).is_file() and (path/n).stat().st_size == size
                        for n, size in result["files"].items())):
            return result
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return None


def mark_unavailable(path, signature, method, error):
    if isinstance(error, InsufficientDevelopmentSupport):
        pd.DataFrame(error.audit).to_csv(path/"development_support.tsv", sep="\t", index=False)
        result = dict(signature=signature, method=method, status=INSUFFICIENT_SUPPORT,
            reason=str(error), development_support=error.audit, prediction_fitted=False,
            test_truth_read=False, calibration_audit_outcomes_used=False,
            files={"development_support.tsv":file_size(path/"development_support.tsv")})
        abm.json_write(path/"UNAVAILABLE.json", result)
        return result
    pd.DataFrame(error.audit).to_csv(path/"reference_support.tsv", sep="\t", index=False)
    abm.json_write(path/"reference_pruning.json", error.composition)
    result = dict(signature=signature, method=method, status=INSUFFICIENT_REFERENCE,
        reason=str(error), reference_support=error.audit, pruning_composition=error.composition,
        prediction_fitted=False, test_truth_read=False,
        files={n:file_size(path/n) for n in ("reference_support.tsv", "reference_pruning.json")})
    abm.json_write(path/"UNAVAILABLE.json", result)
    return result


def publish(stage, target):
    host = load_module("le8_unified_host", HERE/"0.host.py")
    target.parent.mkdir(parents=True, exist_ok=True)
    host.publish(stage, target)


def fit_deviation(args, source, stage, grouping, contract):
    """Keep the older CAD-specific diagnostic endpoint separate from incident risk."""
    host = load_module("le8_unified_host", HERE/"0.host.py")
    if args.Y != "cvd_cad":
        stage.mkdir()
        host.write_tables(stage/"aggregate",
            {"status":pd.DataFrame([dict(status="not_applicable_non_CAD_trait", trait=args.Y)])},
            "deviation", args.biom, args.Y, "not_applicable")
        host.finish(stage, contract)
        return
    path = source/"deviation_LOCAL_ONLY.joblib"
    if not path.is_file():
        raise ValueError("Molecular deviation needs the all-person prepared source; rerun preparation with this method selected")
    config_path = Path(args.deviation_config) if args.deviation_config else ROOT/f"config/deviation_{args.biom}.json"
    dc = host.deviation.read_config(config_path)
    dc = replace(dc, trait=args.Y, layer=args.biom, seed=args.seed, grouping_unit=grouping)
    dc.validate()
    # The historical runner owns creation of its staging directory.
    host.run_deviation(args, abm.load_model(path), stage, dc, contract)


def final_index(scope, destination, run_receipt):
    unavailable = run_receipt.get("unavailable_methods", {})
    abm.json_write(scope/"qc/C1_FINAL_INDEX.json", dict(version=VERSION,
        evaluation=str(destination) if destination is not None else None,
        methods=run_receipt["methods"], unavailable_methods=unavailable,
        failed_methods=run_receipt.get("failed_methods", {}), baseline_commit=BASELINE_COMMIT, source_commit=runtime_source_commit(),
        reliability=run_receipt.get("reliability"),
        status="failed_methods" if run_receipt.get("failed_methods") else PARTIAL_COMPLETION if unavailable else
            "completed_with_unavailable_reliability" if run_receipt.get("reliability",{}).get("status","completed_primary_P1_vs_P0") != "completed_primary_P1_vs_P0" else "completed",
        reliability_status=run_receipt.get("reliability",{}).get("status"),
        completion_policy="current outputs; --replace TRUE overwrites, FALSE reuses complete results"))
    from outputs import publish_flat_results
    publish_flat_results(scope)
    # Short public aggregate names; signed model/input directories remain intact.
    # Every call with a destination follows verified completion or atomic publish.
    from workflow_layout import publish_result
    publish_result(scope/"qc/C1_FINAL_INDEX.json", scope.parent/"qc/c1.json")
    if destination is not None:
        selected = {"C1_unified_evaluation.xlsx": "c1.xlsx", "performance.tsv": "c1.tsv",
                    "contrasts.tsv": "c1.contrasts.tsv", "coverage.tsv": "c1.coverage.tsv"}
        for source_name, public_name in selected.items():
            source_path = Path(destination)/source_name
            if source_path.is_file():
                publish_result(source_path, scope.parent/public_name)


def record_insufficient_support(args, plan, cfg, source, scope, error):
    """Publish an explicit unavailable status, without fitting or opening held-out truth."""
    receipt = dict(version=VERSION, status=INSUFFICIENT_SUPPORT, trait=args.Y, layer=args.biom,
        reason=str(error), requested_methods=plan["methods"], methods={}, evaluation=None,
        horizon=cfg.horizon, minimum_per_class=cfg.min_events, parameters=asdict(cfg),
        method_plan=plan, source={name:file_size(source/name) for name in
            ("development_LOCAL_ONLY.npz", "test_X_LOCAL_ONLY.npz", "host_contract.json")},
        development_support=error.audit,
        test_truth_read=False, calibration_audit_outcomes_used=False,
        prediction_fitted=False, baseline_commit=BASELINE_COMMIT, source_commit=runtime_source_commit())
    signature = hashlib.sha256(json.dumps(receipt, sort_keys=True, default=str).encode()).hexdigest()
    destination = scope/"qc/prediction_status/current"
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".status-", dir=staging_directory(scope)) as temp:
        stage = Path(temp)/"result"; stage.mkdir()
        pd.DataFrame(error.audit).to_csv(stage/"development_support.tsv", sep="\t", index=False)
        abm.json_write(stage/"status.json", receipt)
        publish(stage, destination)
    abm.json_write(scope/"qc/C1_FINAL_INDEX.json", dict(receipt, prediction_status=str(destination)))
    print(f"[LE8] SKIP C1 prediction Y={args.Y} biom={args.biom}: {error}", flush=True)
    print(f"[LE8 C1] Development support audit: {destination/'development_support.tsv'}", flush=True)
    return {"scope_status": INSUFFICIENT_SUPPORT}


def read_public_tables(destination):
    tables = {}
    for path in sorted(destination.glob("*.tsv")):
        try:
            tables[path.stem] = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            # Support-only fallbacks legitimately have no fitted error-calibration
            # table. The manifest/status documents this; resume must still work.
            tables[path.stem] = pd.DataFrame()
    return tables


def evaluate_all(args, c, predictions, source, scope, run_receipt):
    from c1_evaluation import evaluate_methods
    from le8_common import write_tables
    # Every model is frozen before this boundary. Truth bytes identify evaluation
    # reuse only; they are deliberately absent from model-training signatures.
    run_receipt = dict(run_receipt, test_truth_file_size=file_size(source/"test_truth_LOCAL_ONLY.npz"))
    sig = "current"
    destination = scope/"qc/evaluation"
    destination.parent.mkdir(parents=True, exist_ok=True)
    model_times = [(Path(v["path"])/"COMPLETE.json").stat().st_mtime_ns for v in run_receipt.get("methods", {}).values()
                   if (Path(v["path"])/"COMPLETE.json").is_file()]
    evaluation_current = ((destination/"COMPLETE.json").is_file() and
                          (destination/"COMPLETE.json").stat().st_mtime_ns >= max(model_times, default=0))
    if not replace_stage(args,"evaluation") and evaluation_current and completion(destination, sig):
        print("[LE8] DONE C1 REUSE EVALUATION: verified truth, models, settings and aggregates", flush=True)
        final_index(scope, destination, run_receipt)
        return read_public_tables(destination)
    print(f"[LE8] START C1 RECOMPUTE EVALUATION {args.Y}/{args.biom}", flush=True)
    # First semantic test-truth read, after all selected methods were frozen.
    truth = abm.read_npz(source/"test_truth_LOCAL_ONLY.npz")
    tables = evaluate_methods(predictions, truth, c.km, bootstrap=args.bootstrap, seed=args.seed,
        grouping_unit=c.cfg.grouping_unit, primary_coverage=c.cfg.primary_coverage)
    unavailable = {**run_receipt.get("unavailable_methods", {}), **run_receipt.get("failed_methods", {})}
    if unavailable:
        rows = [dict(model=name, component="model_fit", status=item["status"], reason=item["reason"])
                for name, item in unavailable.items()]
        tables["evaluation_status"] = pd.concat([tables["evaluation_status"], pd.DataFrame(rows)], ignore_index=True)
    with tempfile.TemporaryDirectory(prefix=".evaluation-", dir=staging_directory(scope)) as temp:
        stage = Path(temp)/"result"; stage.mkdir()
        # Aggregates only; person-level predictions remain private in methods/.
        for name, frame in tables.items():
            if any(str(k).lower() in {"eid","person_id","participant_id"} for k in frame.columns):
                raise AssertionError("Identifier leaked into public evaluation")
            frame.to_csv(stage/(name+".tsv"), sep="\t", index=False)
        write_tables(stage, tables, "c1_correlate", args.biom, args.Y, sig[:20], workbook="C1_unified_evaluation.xlsx")
        abm.json_write(stage/"evaluation_contract.json", dict(run_receipt, version=VERSION,
            bootstrap=args.bootstrap, evaluation_km="shared build-only reverse Kaplan-Meier",
            primary_coverage=c.cfg.primary_coverage, all_selected_models_frozen_before_test_truth=not bool(unavailable),
            all_available_models_frozen_before_test_truth=True,
            test_provenance="previously_inspected_or_unverified_exploratory",
            conditional_inference="bootstrap conditional on fitted models; no full retraining uncertainty",
            baseline_commit=BASELINE_COMMIT, source_commit=runtime_source_commit()))
        abm.json_write(stage/"COMPLETE.json", dict(signature=sig, files={p.name:file_size(p) for p in stage.iterdir() if p.is_file()}))
        publish(stage, destination)
    final_index(scope, destination, run_receipt)
    print("[LE8 C1] Complete: "+str(destination/"C1_unified_evaluation.xlsx"), flush=True)
    return tables


def run_scope(args, plan, base):
    scope = args.analysis_root/args.Y/args.biom/"c1_correlate"
    scope.mkdir(parents=True, exist_ok=True)
    (scope/"qc").mkdir(exist_ok=True)
    key = hashlib.sha256(str(scope.parent.resolve()).encode()).hexdigest()[:20]
    lock_dir = args.analysis_root / "logs" / "locks"
    lock_dir.mkdir(parents=True, exist_ok=True)
    with open(lock_dir/("le8-scope-"+key+".lock"), "a") as lock:
        shared_flock(lock, fcntl.LOCK_EX)
        args._c1_stage = "prepare"
        source, contract, _ = getattr(args, "_prepared_inputs", None) or prepare(args, plan, scope)
        args._c1_stage = "development_support"
        d, x, clinical, history, grouping = read_inputs(source, contract)
        cfg = build_abm_config("abm_case_preserving", base, trait=args.Y, layer=args.biom,
            seed=args.seed, cores=args.cores, parallel_workers=min(4,args.cores), bootstrap=args.bootstrap,
            grouping_unit=grouping, always_include_indices=clinical, final_multiroute_release="")
        # Shared preprocessing is outcome-blind apart from build-only censoring.
        # Each method checks its OWN required roles below; a missing calibration
        # class must not exclude a build-only count model or an association scan.
        c = context(d, x, cfg, required_roles=set())
        source_fingerprint = {"directory": str(source.resolve()), "trait": args.Y, "layer": args.biom}
        if history is not None: source_fingerprint["history"] = file_size(source/"prevalent_source_LOCAL_ONLY.joblib")
        frozen = frozen_stage_models(args, plan["methods"], source_fingerprint)
        predictions, receipt = {}, dict(methods={}, source=source_fingerprint,
            trait=args.Y, layer=args.biom, horizon=cfg.horizon, grouping_unit=grouping)
        history_config = json.loads(args.history_config.read_text()) if args.history_config else None
        for name in plan["methods"]:
            args._c1_stage = name
            try:
                spec = get_method(name)
                method_cfg = build_abm_config(name if spec.engine in {"twin","twin_count","prevalent_transfer","torch_transformer"}
                    else "abm_case_preserving", base, overrides=plan["overrides"].get(name, {}),
                    trait=args.Y, layer=args.biom, seed=args.seed, cores=args.cores,
                    parallel_workers=min(4,args.cores), bootstrap=args.bootstrap, grouping_unit=grouping,
                    always_include_indices=clinical, final_multiroute_release="")
                # Clinical roles, censoring horizon and assay QC are a study contract,
                # not parameters that may quietly differ between comparison arms.
                for field in ("horizon","feature_missing","sample_missing"):
                    if getattr(method_cfg,field) != getattr(cfg,field):
                        raise ValueError(f"{name}: {field} differs from common evaluation/input contract; use a separate study run")
                device = resolve_device(name, args.device, check_available=args.stage != "evaluation")
                method_context = copy.copy(c)
                method_context.cfg = method_cfg
                method_context.common_support = None
                method_context.common_geometry = None
                try:
                    check_method_support(method_context, spec.engine, args.min_cases_with_biom)
                except InsufficientDevelopmentSupport as error:
                    # Keep support-only decisions apart from frozen fitted models.
                    # Recompute the counts on every run, including older exclusion logs.
                    support_key = hashlib.sha256(json.dumps(dict(source=source_fingerprint,
                        parameters=asdict(method_cfg), minimum=args.min_cases_with_biom),
                        sort_keys=True, default=str).encode()).hexdigest()
                    target = scope/"qc"/"prediction_status"/name/"current"
                    target.mkdir(parents=True, exist_ok=True)
                    unavailable = mark_unavailable(target, support_key, name, error)
                    receipt.setdefault("unavailable_methods", {})[name] = dict(path=str(target),
                        signature=support_key, status=unavailable["status"], reason=str(error))
                    print(f"[LE8] SKIP C1 method={name} Y={args.Y} biom={args.biom}: {error}", flush=True)
                    continue
                extra = {}
                if spec.uses_prevalent: extra["history_config"] = history_config
                if name == "abm_deep": extra["deep_config"] = json.loads(args.deep_config.read_text())
                if name == "molecular_deviation":
                    host = load_module("le8_unified_host", HERE/"0.host.py")
                    config_path = Path(args.deviation_config) if args.deviation_config else ROOT/f"config/deviation_{args.biom}.json"
                    extra["deviation_config"] = json.loads(config_path.read_text())
                    for label, path in {"pgs_source":host.pgs_path(args), "pgs_mapping":args.pgs_mapping,
                                        "deviation_input":source/"deviation_LOCAL_ONLY.joblib"}.items():
                        if path:
                            path = Path(path)
                            extra[label] = {"path":str(path), "size_bytes":file_size(path) if path.is_file() else "missing"}
                signature = name
                target = scope/"tmp/models"/name/"current"
                target.parent.mkdir(parents=True, exist_ok=True)
                stored_signature = signature
                if frozen is not None:
                    entry = frozen["methods"].get(name, frozen["unavailable_methods"].get(name))
                    if entry is None:
                        raise ValueError(f"{name}: missing model-stage result")
                    target = Path(entry["path"]).resolve()
                    allowed = {(scope/"tmp/models"/name).resolve()}
                    if name == "abm_deep":
                        primary_root = Path(os.getenv("LE8_PRIMARY_ANALYSIS_ROOT", str(args.analysis_root))).resolve()
                        allowed.add((primary_root/"exploratory/abm_deep"/args.Y/args.biom/"c1_correlate/tmp/models"/name).resolve())
                    if target.parent not in allowed:
                        raise ValueError(f"{name}: model-stage path is outside its current method directory")
                    stored_signature = entry.get("stored_signature", signature)
                    valid = (completion(target, stored_signature) if name in frozen["methods"] else
                             unavailable_result(target, stored_signature))
                    if not valid or (name in frozen["methods"] and not model_matches_cohort(target, method_context, name)):
                        raise RuntimeError(f"{name}: model-stage output is incomplete or belongs to another cohort")
                    device = entry.get("device", device)
                elif not replace_stage(args,"models"):
                    method_context.history_config_bytes = args.history_config.read_bytes() if args.history_config else None
                    target, stored_signature = reusable_method(scope, name, signature, target,
                        method_cfg=method_cfg, context=method_context, extra=extra, source=source_fingerprint)
                unavailable = unavailable_result(target, stored_signature) if not replace_stage(args,"models") else None
                if unavailable and unavailable["status"] == INSUFFICIENT_SUPPORT:
                    # Admission may have changed; counts above are authoritative.
                    unavailable = None
                if unavailable:
                    print(f"[LE8 C1] Reuse unavailable status for {name}", flush=True)
                elif not replace_stage(args,"models") and completion(target,stored_signature) and model_matches_cohort(target,method_context,name):
                    device = json.loads((target/"COMPLETE.json").read_text()).get("metadata", {}).get("device", device)
                    print(f"[LE8] DONE C1 REUSE MODEL {name}: verified frozen model and projections", flush=True)
                else:
                    if getattr(args, "stage", "all") == "evaluation":
                        raise RuntimeError(f"Evaluation requires a verified frozen {name} model; run models first")
                    print(f"[LE8] START C1 REFIT MODEL method={name} Y={args.Y} biom={args.biom}; device={device['actual']}; all shared test outcomes withheld",flush=True)
                    with tempfile.TemporaryDirectory(prefix=".fit-",dir=target.parent) as temp:
                        stage = Path(temp)/"result"
                        if spec.engine != "deviation": stage.mkdir()
                        try:
                            if name == "abm":
                                from c1_count import fit_count
                                fit_count(method_context, stage)
                            elif name == "glmnet": fit_glmnet(method_context,stage,args.r_bin)
                            elif name == "lightgbm": fit_lightgbm(method_context,stage)
                            elif spec.engine == "twin":
                                fit_twin(method_context,method_cfg,stage,name)
                            elif name == "abm_prevalent": fit_prevalent(method_context,history,stage,history_config)
                            elif name == "abm_deep":
                                resources = load_module("le8_c1_resources",HERE/"0.resources.py")
                                with resources.gpu_lease() if device["actual"]=="cuda" else nullcontext():
                                    fit_deep(method_context,stage,extra["deep_config"],device["actual"])
                            elif name == "molecular_deviation":
                                fit_deviation(args, source, stage, grouping, contract)
                            else: raise AssertionError(name)
                        except (InsufficientDevelopmentSupport, abm.InsufficientReferenceSupport) as error:
                            stage.mkdir(exist_ok=True)
                            unavailable = mark_unavailable(stage, signature, name, error)
                        if unavailable is None:
                            mark_complete(stage,signature,{"method":name,"device":device,"parameters":asdict(method_cfg),
                                "extra":extra,"source":source_fingerprint,"model_contract_version":MODEL_CONTRACT_VERSION})
                        publish(stage,target)
                if unavailable:
                    receipt.setdefault("unavailable_methods", {})[name] = dict(path=str(target),
                        signature=signature, stored_signature=stored_signature,
                        status=unavailable["status"], reason=unavailable["reason"])
                    print(f"[LE8] SKIP C1 method={name} Y={args.Y} biom={args.biom}: {unavailable['reason']}", flush=True)
                    continue
                receipt["methods"][name] = dict(path=str(target), signature=signature,
                                              stored_signature=stored_signature, device=device)
                if spec.engine != "deviation":
                    predictions[name] = abm.load_model(target/"test_predictions_LOCAL_ONLY.joblib")
                if name == "abm_prevalent":
                    predictions["abm_prevalent_control"] = abm.load_model(target/"nohistory_control/test_predictions_LOCAL_ONLY.joblib")
                    receipt["methods"]["abm_prevalent_control"] = dict(path=str(target/"nohistory_control"),
                        signature=signature,device=device,automatic_matched_ablation=True)
            except Exception as error:
                receipt.setdefault("failed_methods", {})[name] = dict(status="failed", reason=str(error))
                print(f"[LE8] ERROR C1 method={name} Y={args.Y} biom={args.biom}: {error}; continue other methods", flush=True)
                traceback.print_exc()
        if frozen is not None and receipt.get("failed_methods"):
            raise RuntimeError("Frozen primary/deep outputs failed verification; unified evaluation was not run")
        if getattr(args, "stage", "all") == "models":
            directory = stage_directory()
            marker = dict(version=VERSION, **receipt, test_outcomes_used=False, status="models_frozen")
            if directory is not None:
                part = os.environ.get("LE8_C1_MODEL_PART")
                if part not in {"primary", "deep"}:
                    raise ValueError("A split model stage needs a primary/deep receipt label")
                atomic_json(directory / f"models.{part}.json", marker)
                return {"failed_methods": receipt.get("failed_methods", {})}
            abm.json_write(scope/"qc/C1_MODELS_FROZEN.json", marker)
            # Standalone reliability must use this run's available models only.
            final_index(scope, None, receipt)
            return {"failed_methods": receipt.get("failed_methods", {})}
        if frozen is not None:
            atomic_json(scope/"qc/C1_MODELS_FROZEN.json", dict(version=VERSION, **receipt,
                test_outcomes_used=False, status="models_frozen"))
        if predictions:
            # Freeze the new reliability heads before either evaluator opens test truth.
            if "glmnet" in predictions and getattr(args, "stage", "all") == "all":
                args._c1_stage = "reliability"
                from c1_evaluation import evaluate_scope, DEFAULT_COVERAGES, RISK_BINS, PRIMARY_COVERAGE
                receipt["reliability"] = evaluate_scope(scope, args.bootstrap, args.seed,
                    DEFAULT_COVERAGES, PRIMARY_COVERAGE, RISK_BINS, replace_stage(args,"reliability"),
                    method_paths={k:Path(v["path"]) for k,v in receipt["methods"].items()},
                    prepared_source=source, protocol_path=args.reliability_protocol)
                print(f"[LE8 C1] FINAL reliability: {receipt['reliability']['status']}", flush=True)
            args._c1_stage = "evaluation"
            result = evaluate_all(args,c,predictions,source,scope,receipt)
        elif receipt.get("unavailable_methods") or receipt.get("failed_methods"):
            final_index(scope, None, receipt)
            result = {}
        else:
            abm.json_write(scope/"qc/C1_FINAL_INDEX.json", dict(version=VERSION,status="completed",**receipt))
            result = {}
        if receipt.get("unavailable_methods"):
            result = dict(result, scope_status=PARTIAL_COMPLETION,
                unavailable_methods=list(receipt["unavailable_methods"]), evaluation_performed=bool(predictions))
        if receipt.get("reliability",{}).get("status","completed_primary_P1_vs_P0") != "completed_primary_P1_vs_P0":
            result = dict(result,scope_status=PARTIAL_RELIABILITY,
                reliability_status=receipt["reliability"]["status"])
        if receipt.get("failed_methods"):
            result = dict(result, failed_methods=receipt["failed_methods"])
        return result


def check_method_support(c, engine, min_cases_with_biom):
    """Apply prediction admission without imposing it on any other analysis."""
    if engine == "deviation":
        return
    required = {"build"} if engine == "twin_count" else {"build", "tune_model", "calibration_fit"}
    if engine in {"twin", "torch_transformer"}:
        required.add("tune_gate")
    audit = check_development_support(c.parts, c.y, c.w, c.qc, c.cfg, required_roles=required)
    # Include all development roles, but never held-out test outcomes.
    d = c.development
    qc = c.qc[:len(d["eid"])]
    cases = int(np.sum(qc & (np.asarray(d["event"]) == 1) & (np.asarray(d["time"]) <= c.cfg.horizon)))
    audit.append(dict(role="all_development_with_biom", known_cases=cases,
        known_controls=int(np.sum(qc & (np.asarray(d["time"]) > c.cfg.horizon))),
        minimum_per_class=min_cases_with_biom, required=True, sufficient=cases >= min_cases_with_biom))
    if cases < min_cases_with_biom:
        raise InsufficientDevelopmentSupport(audit, c.cfg)


def scope_support_preflight(args, plan, base):
    """Prepare inputs and record counts; never exclude a trait/layer or association."""
    scope = args.analysis_root / args.Y / args.biom / "c1_correlate"
    scope.mkdir(parents=True, exist_ok=True)
    (scope/"qc").mkdir(exist_ok=True)
    source, contract, extra = prepare(args, plan, scope)
    args._prepared_inputs = (source, contract, extra)
    cfg = build_abm_config("abm_case_preserving", base, trait=args.Y, layer=args.biom,
                          seed=args.seed, cores=args.cores, final_multiroute_release="")
    with np.load(source / "development_LOCAL_ONLY.npz", allow_pickle=False) as data:
        roles = data["role"].astype(str)
        time, event = np.asarray(data["time"], float), data["event"]
        if not set(roles) <= set(abm.ROLES) or time.shape != roles.shape or event.shape != roles.shape:
            raise ValueError("Invalid development roles or outcome lengths")
        if not np.isfinite(time).all() or (time <= 0).any() or not np.isin(event, [0, 1]).all():
            raise ValueError("Invalid development survival outcomes")
        rows = [dict(role=role, N=int(np.sum(roles == role)),
                     known_cases=int(np.sum((roles == role) & (event == 1) & (time <= cfg.horizon))),
                     known_controls=int(np.sum((roles == role) & (time > cfg.horizon))),
                     phase="before_technical_qc_counts_only") for role in abm.ROLES]
    receipt = dict(status="prepared", trait=args.Y, layer=args.biom,
        development_support=rows, prepared_source=str(source), test_truth_read=False,
        admission="per_method_after_technical_qc", association="independent", pgs="independent")
    pd.DataFrame(rows).to_csv(scope / "qc/support.tsv", sep="\t", index=False)
    abm.json_write(scope / "qc/C1_SCOPE_PREFLIGHT.json", receipt)
    # Supersede historical whole-scope exclusions; they are no longer admission rules.
    log = args.analysis_root / args.Y / args.biom / "qc/exclusions.json"
    if log.is_file():
        abm.json_write(log, dict(receipt, status="superseded_by_per_analysis_support"))
    return receipt


def main(argv=None):
    os.umask(0o077)
    args = parser().parse_args(argv)
    if not args.Y or any(not re.fullmatch(r"[A-Za-z0-9_]+",s) for s in args.Y.split(",")):
        raise ValueError("Invalid trait names")
    if not args.biom or not set(args.biom.split(",")) <= {"prot","met"}:
        raise ValueError("Omic layers must be prot and/or met")
    if args.prepared_source and ("," in args.Y or "," in args.biom):
        raise ValueError("A prepared source belongs to exactly one trait/layer")
    if args.stage == "evaluation":
        args.replace_models = "FALSE"
    plan, base = plan_and_base(args)
    if args.dry_run:
        print(json.dumps(public_plan(args,plan,base),indent=2)); return
    if args.stage in {"all", "models"} or args.preflight:
        dependencies(args,plan)
    if args.preflight:
        if args.prepared_source:
            for name in ("development_LOCAL_ONLY.npz","test_X_LOCAL_ONLY.npz","test_truth_LOCAL_ONLY.npz","host_contract.json"):
                if not (args.prepared_source/name).is_file(): raise FileNotFoundError(args.prepared_source/name)
        print(json.dumps(public_plan(args,plan,base),indent=2))
        print("[LE8 C1] Dependencies and method plan valid; host inputs validated at cohort preparation")
        return
    skipped, failed, partial = [], [], []
    for trait in args.Y.split(","):
        for layer in args.biom.split(","):
            scoped = argparse.Namespace(**vars(args)); scoped.Y,scoped.biom = trait,layer
            # Complete and publish one trait/layer before starting the next.
            # The native engine still validates/reuses its existing results;
            # run_scope independently reuses verified inputs, models and evaluation.
            print(f"[LE8] START C1 Y={trait} biom={layer}: association -> per-method support -> methods -> evaluation", flush=True)
            stage = "scope_support_preflight"
            try:
                # Association does not depend on prediction preparation or role counts.
                if args.stage in {"all", "association"}:
                    try:
                        native_association(scoped)
                    except Exception as error:
                        failed.append(f"{trait}/{layer}/association")
                        print(f"[LE8] ERROR C1 association {trait}/{layer}: {error}; continue prediction", flush=True)
                        traceback.print_exc()
                if args.stage == "association":
                    continue
                if args.stage in {"models", "evaluation"}:
                    bind_prepared_inputs(scoped)
                with threadpool_limits(args.cores):
                    scope_support_preflight(scoped, plan, base)
                if args.stage == "prepare":
                    freeze_prepared_inputs(scoped, scoped._prepared_inputs[0])
                    print(f"[LE8] DONE C1 prepare Y={trait} biom={layer}; no test truth read", flush=True)
                    continue
                stage = "methods/evaluation"
                with threadpool_limits(args.cores): result = run_scope(scoped,plan,base)
            except Exception as error:
                if trait == "depression" and getattr(error, "missing_outcome_column", None) == "fod_icd10_depression":
                    reason = str(error)
                    record = dict(status="SKIPPED", reason=reason, scope=f"{trait}.{layer}",
                                  missing_column="fod_icd10_depression", exit_code=0)
                    atomic_json(Path(args.analysis_root)/trait/layer/"qc/C1_SCOPE_SKIP.json", record)
                    marker = os.getenv("LE8_SCOPE_SKIP_FILE")
                    if marker:
                        atomic_json(Path(marker), record)
                    print(f"[LE8] SKIP C1 Y={trait} biom={layer}: {reason}; exit=0", flush=True)
                    continue
                # A trait failure must not discard earlier results or prevent
                # the remaining traits from running. Ctrl+C still propagates.
                failed.append(f"{trait}/{layer}")
                stage = getattr(scoped, "_c1_stage", stage)
                # SKIP is forwarded immediately by the public console filter.
                print(f"[LE8] SKIP C1 Y={trait} biom={layer} ERROR stage={stage}: "
                      f"{type(error).__name__}: {error}; continue to next trait/layer; saved results retained", flush=True)
                traceback.print_exc()
                continue
            if result and result.get("failed_methods"):
                failed.extend(f"{trait}/{layer}/{name}" for name in result["failed_methods"])
                print(f"[LE8] PARTIAL C1 Y={trait} biom={layer}: failed methods recorded; available results retained", flush=True)
                continue
            if result and result.get("scope_status") == INSUFFICIENT_SUPPORT:
                skipped.append(f"{trait}/{layer}")
                continue
            if result and result.get("scope_status") == PARTIAL_RELIABILITY:
                reason=result["reliability_status"]
                partial.append(f"{trait}/{layer}: {reason}")
                print(f"[LE8] PARTIAL C1 Y={trait} biom={layer}: frozen risk methods retained; primary reliability unavailable ({reason})",flush=True)
                continue
            if result and result.get("scope_status") == PARTIAL_COMPLETION:
                methods = ", ".join(result["unavailable_methods"])
                partial.append(f"{trait}/{layer}: {methods}")
                detail = "available methods evaluated" if result.get("evaluation_performed") else "no common evaluation"
                print(f"[LE8] PARTIAL C1 Y={trait} biom={layer}: {detail}; unavailable: {methods}", flush=True)
                continue
            print(f"[LE8] DONE C1 Y={trait} biom={layer}: unified C1 completed", flush=True)
    if skipped:
        print(f"[LE8] PARTIAL C1 SUMMARY: insufficient sample support for {len(skipped)} scope(s): "
              f"{', '.join(skipped)}. Association/PGS use independent support checks; "
              "see qc/C1_SCOPE_PREFLIGHT.json and qc/C1_FINAL_INDEX.json. Existing results retained.", flush=True)
    if partial:
        print("[LE8] PARTIAL C1 SUMMARY: unavailable methods or primary reliability: " + "; ".join(partial), flush=True)
    if failed:
        print("[LE8] PARTIAL C1 SUMMARY: failed scopes: " + ", ".join(failed) +
              ". All remaining trait/layer scopes were attempted; resume with --replace FALSE.", flush=True)
    if os.getenv("WORKFLOW_WORKER") != "1":
        from c1_report import write_reliability_cross_trait_summary, write_pgs_cross_trait_summary
        for layer in args.biom.split(","):
            write_reliability_cross_trait_summary(args.analysis_root, layer)
            write_pgs_cross_trait_summary(args.analysis_root, layer)
    # Signal failures only after the entire batch, so callers can detect them.
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
