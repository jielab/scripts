#!/usr/bin/env python3
"""Unified C1: association plus named prediction methods on one frozen cohort.

All selected prediction methods finish fit and outcome-blind projection before
test truth is opened by the common evaluator. Old ABM directories are read-only
historical results; new results live under c1_correlate/methods and evaluation.
"""
from __future__ import annotations
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
from analysis_layout import normalize_cohort_request
from c1_protocol import (VERSION as REGISTRY_VERSION, get_method, load_method_plan,
                        parse_methods, build_abm_config, resolve_device, scientific_fingerprint, INPUT_CONTRACT_VERSION,
                        MODEL_CONTRACT_VERSION, LEGACY_MODEL_CONTRACT_VERSION)
from c1_models import (InsufficientDevelopmentSupport, context, fit_glmnet,
                               fit_lightgbm, fit_twin, fit_prevalent, fit_deep)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
VERSION = "LE8-C1-FINAL-2026-10-26"
BASELINE_COMMIT = "40a884324f30219f865e8a7a49ac9708a913e742"
INSUFFICIENT_SUPPORT = "skipped_insufficient_development_support"
INSUFFICIENT_REFERENCE = "skipped_insufficient_reference_support"
PARTIAL_COMPLETION = "completed_with_unavailable_methods"
PARTIAL_RELIABILITY = "completed_with_unavailable_reliability"
INSUFFICIENT_SCOPE_SUPPORT = "skipped_insufficient_scope_support"


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


def sha_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""): h.update(block)
    return h.hexdigest()


def parser():
    p = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    p.add_argument("--Y", default="cvd_cad,ra")
    p.add_argument("--biom", default="prot,met")
    p.add_argument("--analysis-root", "--dir-analysis-root", type=Path,
        default=Path(os.getenv("LE8_ANALYSIS_ROOT", "/mnt/d/analysis/le8")))
    p.add_argument("--methods", help="Comma-separated method IDs; defaults to the methods in --config")
    p.add_argument("--config", type=Path, default=ROOT/"config/c1_methods.json",
                   help="Scientific method-plan JSON; a legacy flat ABM JSON is also accepted")
    p.add_argument("--abm-config", type=Path, default=ROOT/"config/abm.json")
    p.add_argument("--deep-config", type=Path, default=ROOT/"config/abm_GPU.json")
    p.add_argument("--history-config", type=Path, default=ROOT/"config/c1_history.json",
                   help="PrevalentConfig JSON; weights are explicit sensitivity assumptions")
    p.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    p.add_argument("--association", choices=["auto", "run", "skip"], default="run")
    p.add_argument("--pgs", choices=["TRUE", "FALSE"],
                   type=lambda value: {"RUN": "TRUE", "SKIP": "FALSE"}.get(value.upper(), value.upper()),
                   default=os.getenv("C1_PGS_MODE", "FALSE"),
                   help="Enable molecular PGS association scan; default FALSE (not required by primary ABM)")
    p.add_argument("--prepared-source", type=Path, help="Reuse a frozen prepared input directory (one trait/layer only)")
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--min-cases-with-biom", "--min_cases_with_biom", type=int,
                   default=int(os.getenv("C1_MIN_CASES_WITH_BIOM", "50")),
                   help="Minimum known incident cases with omics in the development cohort [50]")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--cores", type=int, default=8)
    p.add_argument("--replace", choices=["TRUE", "FALSE"], default="FALSE")
    for stage in ("association", "models", "evaluation", "reliability"):
        p.add_argument("--replace-" + stage, choices=["TRUE", "FALSE"], default=None,
                       help="Override --replace for this C1 stage only")
    p.add_argument("--reliability-protocol", type=Path, default=ROOT/"config/c1_reliability_protocol.json")
    p.add_argument("--r-bin", default="Rscript")
    p.add_argument("--host-args", default="")
    p.add_argument("--preflight", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    for flag in ("ukb-phe", "ukb-rel-file", "group-file", "group-col", "outer-roster", "end-date",
                 "diagnosis-col", "shared-covariates", "deviation-config", "objective-file", "pgs-mapping"):
        p.add_argument("--"+flag)
    return p


def plan_and_base(args):
    base = json.loads(args.abm_config.read_text())
    config = args.config
    if config:
        obj = json.loads(config.read_text())
        if "methods" not in obj and "overrides" not in obj:
            base, config = obj, None
    plan = load_method_plan(config, methods=args.methods)
    # The old base config contains a frozen historic release. New scientific
    # presets replace it; a selected abm_legacy keeps legacy scientific choices.
    base["seed"] = args.seed
    if args.bootstrap < 1 or args.cores < 1 or args.min_cases_with_biom < 1:
        raise ValueError("Positive bootstrap, core and minimum-case counts are required")
    return plan, base


def public_plan(args, plan, base):
    methods = []
    for name in plan["methods"]:
        spec = get_method(name)
        cfg = build_abm_config(name if spec.engine in {"twin", "prevalent_transfer", "torch_transformer"}
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
        execution_order="trait, then layer: sample support -> association -> methods -> evaluation; finish each scope before the next",
        common_evaluation="same test people, outcome horizon, build censoring model and paired family draws",
        test_provenance="previously_inspected_or_unverified_exploratory",
        proposed_outputs="<analysis-root>/<Y>/<biom>/c1_correlate/{methods,evaluation/FINAL,reliability/FINAL}")


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
    cmd = ["bash", str(HERE/"0.engine.sh"), "c1_correlate", "--Y", args.Y, "--biom", args.biom,
           "--analysis-root", str(args.analysis_root), "--r-bin", args.r_bin,
           "--cores", str(args.cores), "--seed", str(args.seed), "--replace", "TRUE" if replace_stage(args,"association") else "FALSE"]
    if args.ukb_phe: cmd += ["--ukb-phe", args.ukb_phe]
    env = os.environ.copy()
    env.pop("C1_MIN_CASES_WITH_BIOM", None)  # Scope admission must not invalidate R model caches.
    env["LE8_C1_STAGE"] = "c1_correlate" # only association; prediction handled below
    env["C1_PGS_MODE"] = {"TRUE": "run", "FALSE": "skip"}[args.pgs]
    for key, var in (("end_date", "DATE_FOLLOW_END"), ("diagnosis_col", "LE8_Y_DATE"),
                     ("shared_covariates", "LE8_VARS_ADJ"), ("outer_roster", "LE8_OUTER_ROSTER"),
                     ("group_file", "LE8_GROUP_FILE"), ("group_col", "LE8_GROUP_COLUMN")):
        value = getattr(args, key, None)
        if value: env[var] = str(value)
    print("[LE8 C1] Association: " + shlex.join(cmd), flush=True)
    subprocess.run(cmd, env=env, check=True)
    return {"status": "completed_or_reused", "engine": "existing native R association", "prediction_refit": False}


def input_source_signature(host_args, extra_files=None):
    """Bind input reuse to explicit source bytes and the frozen split contract."""
    fields = vars(host_args).copy()
    for key in ("cores", "shared_config_sources"):
        fields.pop(key, None)
    files = {}
    if not host_args.demo:
        for key in ("phe_file", "omics_file", "met_map", "group_file", "outer_roster", "ukb_rel_file"):
            value = getattr(host_args, key, None)
            if value and Path(value).is_file():
                files[key] = {"path": str(Path(value).resolve()), "sha256": sha_file(value)}
    for key, value in (extra_files or {}).items():
        if value:
            path = Path(value).resolve()
            files[key] = {"path":str(path), "sha256":sha_file(path) if path.is_file() else "missing"}
    payload = dict(input_contract_version=INPUT_CONTRACT_VERSION, fields=fields, files=files)
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()
    return digest, payload


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


def prepare(args, plan, scope):
    host = load_module("le8_unified_host", HERE/"0.host.py")
    if args.prepared_source:
        source = args.prepared_source.resolve()
        for name in ("development_LOCAL_ONLY.npz", "test_X_LOCAL_ONLY.npz", "test_truth_LOCAL_ONLY.npz", "host_contract.json"):
            if not (source/name).is_file(): raise FileNotFoundError(source/name)
        contract = json.loads((source/"host_contract.json").read_text())
        identity = contract.get("config", contract.get("input_config", {}))
        trait, layer = contract.get("trait", identity.get("trait")), contract.get("layer", identity.get("layer", identity.get("biom")))
        if trait and trait != args.Y or layer and layer != args.biom:
            raise ValueError("Prepared source trait/layer differs from the requested analysis")
        if "abm_prevalent" in plan["methods"] and not (source/"prevalent_source_LOCAL_ONLY.joblib").is_file():
            raise ValueError("This legacy incident-only cache has no historical source. Run C1 from host inputs to prepare history.")
        return source, contract, None
    args.include_prevalent = "abm_prevalent" in plan["methods"]
    h = host.host_config(args)
    extra_files = {"objective_file":args.objective_file} if "molecular_deviation" in plan["methods"] else {}
    signature, request = input_source_signature(h, extra_files)
    if "molecular_deviation" in plan["methods"]:
        signature = hashlib.sha256((signature+"|all_person_deviation_bundle").encode()).hexdigest()
    source = scope/"_LOCAL_ONLY"/"inputs"/signature[:20]
    # Prepared cohorts depend on data and split settings, never Python source bytes.
    # Read older input requests without their obsolete code/version fields.
    candidates = [source, *sorted(source.parent.glob("*"), reverse=True)]
    required = {"development_LOCAL_ONLY.npz", "test_X_LOCAL_ONLY.npz",
                "test_truth_LOCAL_ONLY.npz", "host_contract.json"}
    if "abm_prevalent" in plan["methods"]:
        required.add("prevalent_source_LOCAL_ONLY.joblib")
    if "molecular_deviation" in plan["methods"]:
        required.add("deviation_LOCAL_ONLY.joblib")
    for candidate in dict.fromkeys(candidates):
        marker = candidate/"INPUTS_COMPLETE.json"
        if not marker.is_file():
            continue
        saved = json.loads(marker.read_text())
        old_request = {k:v for k,v in saved.get("request", {}).items()
                       if k not in {"version", "preparation_code"}}
        old_request.setdefault("input_contract_version", "LE8-C1-input-2026-10-22-v1")
        if normalize_cohort_request(old_request) != normalize_cohort_request(json.loads(json.dumps(request, default=str))):
            continue
        files = saved.get("files", {})
        if required <= files.keys() and all((candidate/n).is_file() and sha_file(candidate/n)==v
                                           for n,v in files.items()):
            print("[LE8 C1] Reusing verified prepared inputs", flush=True)
            return candidate, json.loads((candidate/"host_contract.json").read_text()), None
        raise RuntimeError("Prepared input cache integrity failed; keep evidence and repair the input bundle")
    source.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".prepare-", dir=source.parent) as temp:
        temp = Path(temp)
        incident, deviation, audit = host.prepare_host(args, temp)
        staged = temp/"input_bundle"
        save_prepared(incident, deviation, audit, staged)
        if "molecular_deviation" in plan["methods"]:
            joblib.dump(deviation, staged/"deviation_LOCAL_ONLY.joblib", compress=3)
        after, _ = input_source_signature(h, extra_files)
        if "molecular_deviation" in plan["methods"]:
            after = hashlib.sha256((after+"|all_person_deviation_bundle").encode()).hexdigest()
        if after != signature: raise RuntimeError("Inputs changed during preparation")
        abm.json_write(staged/"INPUTS_COMPLETE.json", dict(source_signature=signature, request=request,
            files={p.name:sha_file(p) for p in staged.iterdir() if p.is_file()}))
        staged.rename(source)
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


def completion(path, signature):
    try:
        value = json.loads((path/"COMPLETE.json").read_text())
        return value.get("signature") == signature and bool(value.get("files")) and all(
            (path/n).is_file() and sha_file(path/n)==s for n,s in value["files"].items())
    except (OSError, ValueError, KeyError): return False


NONSCIENTIFIC_CONFIG_KEYS = {"cores", "parallel_workers", "device", "requested_device", "actual_device"}


def _canonical_model_config(value):
    return json.loads(json.dumps({k:v for k,v in dict(value or {}).items()
        if k not in NONSCIENTIFIC_CONFIG_KEYS}, sort_keys=True, default=str))


def _current_role_contract(c):
    return {role:dict(N=int(len(ix)), id_hash=abm.fingerprint(np.asarray(c.ids[ix]).astype(str)))
            for role,ix in c.parts.items()}


def _verified_stored_complete(path):
    try:
        obj = json.loads((path/"COMPLETE.json").read_text())
        return obj if obj.get("signature") and completion(path,obj["signature"]) else None
    except (OSError, ValueError, TypeError):
        return None


def _legacy_method_semantically_matches(path, name, cfg, c, extra=None, source=None):
    """Read-only semantic check of a frozen model; never rewrite its receipt."""
    marker = _verified_stored_complete(path)
    if marker is None:
        return None
    meta = marker.get("metadata", {})
    if meta.get("model_contract_version", LEGACY_MODEL_CONTRACT_VERSION) != MODEL_CONTRACT_VERSION:
        return None
    wanted = _canonical_model_config(asdict(cfg))
    if meta.get("method") != name or _canonical_model_config(meta.get("parameters")) != wanted:
        return None
    extra = extra or {}
    if meta.get("source") is not None and source != meta["source"]:
        return None
    try:
        manifest = json.loads((path/"MODEL_FROZEN.json").read_text())
        if manifest.get("method") != name or _canonical_model_config(manifest.get("config")) != wanted:
            return None
        roles = manifest.get("roles", {})
        if roles != _current_role_contract(c):
            return None
        # Both development audit and test rosters must match, including family IDs.
        for role in ("calibration_audit", "test"):
            pred = abm.load_model(path/(role+"_predictions_LOCAL_ONLY.joblib"))
            ix = c.parts[role]
            if not np.array_equal(np.asarray(pred["eid"]).astype(str), np.asarray(c.ids[ix]).astype(str)):
                return None
            if not np.array_equal(np.asarray(pred["group"]).astype(str), np.asarray(c.groups[ix]).astype(str)):
                return None
        fitted = abm.load_model(path/"model_LOCAL_ONLY.joblib")
        if not isinstance(fitted, dict) or list(map(str,fitted.get("feature_names",[]))) != list(map(str,c.names)):
            return None
        if fitted.get("config") is not None and _canonical_model_config(fitted["config"]) != wanted:
            return None
        evidence = dict(checks=["complete_checksums", "parameters", "all_role_rosters", "test_and_audit_families", "feature_schema"])
        if meta.get("extra") is not None:
            if meta["extra"] != extra:
                return None
        elif name == "abm_prevalent":
            from c1_models import PrevalentConfig
            wanted_history = asdict(PrevalentConfig(**{"seed":cfg.seed, **(extra.get("history_config") or {})}))
            old_history = manifest.get("source_transform",{}).get("config")
            if old_history is not None:
                if _canonical_model_config(old_history) != _canonical_model_config(wanted_history):
                    return None
                evidence["history_verification"] = "recorded_source_transform_config"
            else:
                # Exact published bytes, not a checksum of executable code.
                known = b'{\n  "n_components": 12,\n  "n_prototypes": 6,\n  "min_prevalent": 20,\n  "age_band_years": 10.0,\n  "match_sex": true,\n  "weighting_policy": "duration_treatment",\n  "duration_half_weight_years": 5.0,\n  "unknown_duration_weight": 0.5,\n  "medication_weight": 0.5,\n  "unknown_medication_weight": 0.75,\n  "max_neighbours": 64\n}\n'
                current = getattr(c,"history_config_bytes",None)
                if current != known or extra.get("history_config") != json.loads(known):
                    return None
                evidence["history_verification"] = "exception_exact_published_2026_10_22_history_bytes"
        elif name == "abm_deep" and set(extra)=={"deep_config"}:
            old_deep=json.loads((path/"transformer/MODEL_FROZEN.json").read_text()).get("config")
            if old_deep != extra["deep_config"]:return None
            evidence["deep_verification"]="recorded_transformer_config"
        elif extra:
            # Additional deep/deviation settings must be recorded, not inferred from a name.
            return None
        return marker["signature"], evidence
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return None


def reusable_method(scope, name, signature, default, *, method_cfg=None, context=None, extra=None, source=None):
    if completion(default, signature):
        return default, signature
    ledger_path = scope/"_LOCAL_ONLY"/"model_reuse.json"
    try:
        ledger = json.loads(ledger_path.read_text())
        if not isinstance(ledger,dict): ledger = {}
    except (OSError,ValueError):
        ledger = {}
    root = (scope/"methods"/name).resolve()
    accepted_signatures={signature}
    if method_cfg is not None and source is not None and MODEL_CONTRACT_VERSION==LEGACY_MODEL_CONTRACT_VERSION:
        # The immediately preceding semantic key had no explicit contract version.
        # It already bound exact data and settings, with executable code excluded.
        parameters={**asdict(method_cfg),**(extra or {})}
        for key in NONSCIENTIFIC_CONFIG_KEYS:parameters.pop(key,None)
        previous=dict(registry_version=REGISTRY_VERSION,method=name,parameters=parameters,input_contract=source)
        accepted_signatures.add(hashlib.sha256(json.dumps(previous,sort_keys=True,separators=(",",":"),allow_nan=False,default=str).encode()).hexdigest())
    for relative, entry in ledger.items():
        path = (scope/relative).resolve()
        if not isinstance(entry,dict) or path.parent != root or entry.get("signature") not in accepted_signatures:
            continue
        stored = entry.get("stored_signature")
        if stored and (completion(path,stored) or unavailable_result(path,stored)):
            return path,stored
    if method_cfg is not None and context is not None and root.is_dir():
        matches=[]
        for path in root.iterdir():
            if not path.is_dir() or path.is_symlink() or path.name.startswith("."): continue
            result = _legacy_method_semantically_matches(path,name,method_cfg,context,extra,source)
            if result: matches.append((path,*result))
        if matches:
            path,stored,evidence=max(matches,key=lambda row:row[0].stat().st_mtime_ns)
            ledger[str(path.relative_to(scope.resolve()))] = dict(signature=signature,stored_signature=stored,
                migration="semantic_contract_match",model_contract_version=MODEL_CONTRACT_VERSION,**evidence)
            ledger_path.parent.mkdir(parents=True,exist_ok=True)
            with tempfile.NamedTemporaryFile(mode="w",dir=ledger_path.parent,prefix=".reuse-",delete=False) as handle:
                json.dump(ledger,handle,indent=2); handle.write("\n"); temp=Path(handle.name)
            temp.replace(ledger_path)
            return path,stored
    return default,signature


def mark_complete(path, signature, metadata):
    # Actual neural/R weights and support donors are part of cache integrity,
    # not just a root-level pointer to them. Optional deviation has nested files.
    anchors = [p for p in path.rglob("*") if p.is_file() and
        (p.name.endswith(("LOCAL_ONLY.joblib", "LOCAL_ONLY.pt", "LOCAL_ONLY.rds")) or
         p.name in {"MODEL_FROZEN.json", "model.joblib", "table_manifest.json", "host_contract.json"})]
    if not anchors:
        raise ValueError("Cannot mark a method complete without a frozen model or verified result anchor")
    abm.json_write(path/"COMPLETE.json", dict(signature=signature, metadata=metadata,
        files={str(p.relative_to(path)):sha_file(p) for p in anchors}))


def unavailable_result(path, signature):
    """Unavailable is a verified status, never a completed fitted model."""
    try:
        result = json.loads((path/"UNAVAILABLE.json").read_text())
        if (result.get("signature") == signature and result.get("status") in {INSUFFICIENT_REFERENCE, INSUFFICIENT_SUPPORT}
                and result.get("files") and not (path/"COMPLETE.json").exists()
                and all((path/n).is_file() and sha_file(path/n) == digest
                        for n, digest in result["files"].items())):
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
            files={"development_support.tsv":sha_file(path/"development_support.tsv")})
        abm.json_write(path/"UNAVAILABLE.json", result)
        return result
    pd.DataFrame(error.audit).to_csv(path/"reference_support.tsv", sep="\t", index=False)
    abm.json_write(path/"reference_pruning.json", error.composition)
    result = dict(signature=signature, method=method, status=INSUFFICIENT_REFERENCE,
        reason=str(error), reference_support=error.audit, pruning_composition=error.composition,
        prediction_fitted=False, test_truth_read=False,
        files={n:sha_file(path/n) for n in ("reference_support.tsv", "reference_pruning.json")})
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
    abm.json_write(scope/"C1_FINAL_INDEX.json", dict(version=VERSION,
        evaluation=str(destination) if destination is not None else None,
        methods=run_receipt["methods"], unavailable_methods=unavailable, baseline_commit=BASELINE_COMMIT, source_commit=runtime_source_commit(),
        reliability=run_receipt.get("reliability"),
        status=PARTIAL_COMPLETION if unavailable else
            "completed_with_unavailable_reliability" if run_receipt.get("reliability",{}).get("status","completed_primary_P1_vs_P0") != "completed_primary_P1_vs_P0" else "completed",
        reliability_status=run_receipt.get("reliability",{}).get("status"),
        supplementary_historical_results="existing c1_abm and c1_abm_GPU directories retained read-only"))


def record_insufficient_support(args, plan, cfg, source, scope, error):
    """Publish an explicit unavailable status, without fitting or opening held-out truth."""
    receipt = dict(version=VERSION, status=INSUFFICIENT_SUPPORT, trait=args.Y, layer=args.biom,
        reason=str(error), requested_methods=plan["methods"], methods={}, evaluation=None,
        horizon=cfg.horizon, minimum_per_class=cfg.min_events, parameters=asdict(cfg),
        method_plan=plan, source={name:sha_file(source/name) for name in
            ("development_LOCAL_ONLY.npz", "test_X_LOCAL_ONLY.npz", "host_contract.json")},
        development_support=error.audit,
        test_truth_read=False, calibration_audit_outcomes_used=False,
        prediction_fitted=False, baseline_commit=BASELINE_COMMIT, source_commit=runtime_source_commit())
    signature = hashlib.sha256(json.dumps(receipt, sort_keys=True, default=str).encode()).hexdigest()
    destination = scope/"prediction_status"/signature[:20]
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".status-", dir=destination.parent) as temp:
        stage = Path(temp)/"result"; stage.mkdir()
        pd.DataFrame(error.audit).to_csv(stage/"development_support.tsv", sep="\t", index=False)
        abm.json_write(stage/"status.json", receipt)
        publish(stage, destination)
    abm.json_write(scope/"C1_FINAL_INDEX.json", dict(receipt, prediction_status=str(destination)))
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
    run_receipt = dict(run_receipt, test_truth_sha256=sha_file(source/"test_truth_LOCAL_ONLY.npz"))
    sig = hashlib.sha256(json.dumps(dict({k:v for k,v in run_receipt.items() if k != "reliability"}, bootstrap=args.bootstrap,
        version=VERSION), sort_keys=True, default=str).encode()).hexdigest()
    destination = scope/"evaluation"/"FINAL"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not replace_stage(args,"evaluation") and completion(destination, sig):
        print("[LE8] DONE C1 REUSE EVALUATION: verified truth, models, settings and aggregates", flush=True)
        final_index(scope, destination, run_receipt)
        return read_public_tables(destination)
    print(f"[LE8] START C1 RECOMPUTE EVALUATION {args.Y}/{args.biom}", flush=True)
    # First semantic test-truth read, after all selected methods were frozen.
    truth = abm.read_npz(source/"test_truth_LOCAL_ONLY.npz")
    tables = evaluate_methods(predictions, truth, c.km, bootstrap=args.bootstrap, seed=args.seed,
        grouping_unit=c.cfg.grouping_unit, primary_coverage=c.cfg.primary_coverage)
    unavailable = run_receipt.get("unavailable_methods", {})
    if unavailable:
        rows = [dict(model=name, component="model_fit", status=item["status"], reason=item["reason"])
                for name, item in unavailable.items()]
        tables["evaluation_status"] = pd.concat([tables["evaluation_status"], pd.DataFrame(rows)], ignore_index=True)
    with tempfile.TemporaryDirectory(prefix=".evaluation-", dir=destination.parent) as temp:
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
        abm.json_write(stage/"COMPLETE.json", dict(signature=sig, files={p.name:sha_file(p) for p in stage.iterdir() if p.is_file()}))
        publish(stage, destination)
    final_index(scope, destination, run_receipt)
    print("[LE8 C1] Complete: "+str(destination/"C1_unified_evaluation.xlsx"), flush=True)
    return tables


def run_scope(args, plan, base):
    scope = args.analysis_root/args.Y/args.biom/"c1_correlate"
    scope.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(str(scope.resolve()).encode()).hexdigest()[:20]
    with open(Path("/tmp")/("le8-unified-c1-"+key+".lock"), "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        args._c1_stage = "prepare"
        source, contract, _ = getattr(args, "_prepared_inputs", None) or prepare(args, plan, scope)
        args._c1_stage = "development_support"
        d, x, clinical, history, grouping = read_inputs(source, contract)
        cfg = build_abm_config("abm_case_preserving", base, trait=args.Y, layer=args.biom,
            seed=args.seed, cores=args.cores, parallel_workers=min(4,args.cores), bootstrap=args.bootstrap,
            grouping_unit=grouping, always_include_indices=clinical, final_multiroute_release="")
        try:
            c = context(d, x, cfg)
        except InsufficientDevelopmentSupport as error:
            return record_insufficient_support(args, plan, cfg, source, scope, error)
        source_fingerprint = {name:sha_file(source/name) for name in
            ("development_LOCAL_ONLY.npz", "test_X_LOCAL_ONLY.npz", "host_contract.json")}
        if history is not None: source_fingerprint["history"] = sha_file(source/"prevalent_source_LOCAL_ONLY.joblib")
        predictions, receipt = {}, dict(methods={}, source=source_fingerprint,
            trait=args.Y, layer=args.biom, horizon=cfg.horizon, grouping_unit=grouping)
        history_config = json.loads(args.history_config.read_text()) if args.history_config else None
        for name in plan["methods"]:
            args._c1_stage = name
            spec = get_method(name)
            method_cfg = build_abm_config(name if spec.engine in {"twin","prevalent_transfer","torch_transformer"}
                else "abm_case_preserving", base, overrides=plan["overrides"].get(name, {}),
                trait=args.Y, layer=args.biom, seed=args.seed, cores=args.cores,
                parallel_workers=min(4,args.cores), bootstrap=args.bootstrap, grouping_unit=grouping,
                always_include_indices=clinical, final_multiroute_release="")
            # Clinical roles, censoring horizon and assay QC are a study contract,
            # not parameters that may quietly differ between comparison arms.
            for field in ("horizon","feature_missing","sample_missing"):
                if getattr(method_cfg,field) != getattr(cfg,field):
                    raise ValueError(f"{name}: {field} differs from common evaluation/input contract; use a separate study run")
            device = resolve_device(name,args.device)
            method_context = copy.copy(c)
            method_context.cfg = method_cfg
            method_context.common_support = None
            method_context.common_geometry = None
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
                        extra[label] = {"path":str(path), "sha256":sha_file(path) if path.is_file() else "missing"}
            signature = scientific_fingerprint(name, {**asdict(method_cfg), **extra}, source_fingerprint)
            target = scope/"methods"/name/(signature[:16]+"_"+device["actual"])
            target.parent.mkdir(parents=True, exist_ok=True)
            stored_signature = signature
            if not replace_stage(args,"models"):
                method_context.history_config_bytes = args.history_config.read_bytes() if args.history_config else None
                target, stored_signature = reusable_method(scope, name, signature, target,
                    method_cfg=method_cfg, context=method_context, extra=extra, source=source_fingerprint)
            unavailable = unavailable_result(target, stored_signature) if not replace_stage(args,"models") else None
            if unavailable:
                print(f"[LE8 C1] Reuse unavailable status for {name}", flush=True)
            elif not replace_stage(args,"models") and completion(target,stored_signature):
                print(f"[LE8] DONE C1 REUSE MODEL {name}: verified frozen model and projections", flush=True)
            else:
                print(f"[LE8] START C1 REFIT MODEL method={name} Y={args.Y} biom={args.biom}; device={device['actual']}; all shared test outcomes withheld",flush=True)
                with tempfile.TemporaryDirectory(prefix=".fit-",dir=target.parent) as temp:
                    stage = Path(temp)/"result"
                    if spec.engine != "deviation": stage.mkdir()
                    if name == "glmnet": fit_glmnet(method_context,stage,args.r_bin)
                    elif name == "lightgbm": fit_lightgbm(method_context,stage)
                    elif spec.engine == "twin":
                        try:
                            fit_twin(method_context,method_cfg,stage,name)
                        except InsufficientDevelopmentSupport as error:
                            unavailable = mark_unavailable(stage, signature, name, error)
                        except abm.InsufficientReferenceSupport as error:
                            if name != "abm_legacy":
                                raise
                            unavailable = mark_unavailable(stage, signature, name, error)
                    elif name == "abm_prevalent": fit_prevalent(method_context,history,stage,history_config)
                    elif name == "abm_deep":
                        resources = load_module("le8_c1_resources",HERE/"0.resources.py")
                        with resources.gpu_lease() if device["actual"]=="cuda" else nullcontext():
                            fit_deep(method_context,stage,extra["deep_config"],device["actual"])
                    elif name == "molecular_deviation":
                        fit_deviation(args, source, stage, grouping, contract)
                    else: raise AssertionError(name)
                    if unavailable is None:
                        mark_complete(stage,signature,{"method":name,"device":device,"parameters":asdict(method_cfg),
                            "extra":extra,"source":source_fingerprint,"model_contract_version":MODEL_CONTRACT_VERSION})
                    publish(stage,target)
            if unavailable:
                receipt.setdefault("unavailable_methods", {})[name] = dict(path=str(target),
                    signature=signature, status=unavailable["status"], reason=unavailable["reason"])
                print(f"[LE8] SKIP C1 method={name} Y={args.Y} biom={args.biom}: {unavailable['reason']}", flush=True)
                continue
            receipt["methods"][name] = dict(path=str(target),signature=signature,device=device)
            if spec.engine != "deviation":
                predictions[name] = abm.load_model(target/"test_predictions_LOCAL_ONLY.joblib")
            if name == "abm_prevalent":
                predictions["abm_prevalent_control"] = abm.load_model(target/"nohistory_control/test_predictions_LOCAL_ONLY.joblib")
                receipt["methods"]["abm_prevalent_control"] = dict(path=str(target/"nohistory_control"),
                    signature=signature,device=device,automatic_matched_ablation=True)
        if predictions:
            # Freeze the new reliability heads before either evaluator opens test truth.
            if "glmnet" in predictions:
                args._c1_stage = "reliability"
                from c1_evaluation import evaluate_scope, DEFAULT_COVERAGES, RISK_BINS, PRIMARY_COVERAGE
                receipt["reliability"] = evaluate_scope(scope, args.bootstrap, args.seed,
                    DEFAULT_COVERAGES, PRIMARY_COVERAGE, RISK_BINS, replace_stage(args,"reliability"),
                    method_paths={k:Path(v["path"]) for k,v in receipt["methods"].items()},
                    prepared_source=source, protocol_path=args.reliability_protocol)
                print(f"[LE8 C1] FINAL reliability: {receipt['reliability']['status']}", flush=True)
            args._c1_stage = "evaluation"
            result = evaluate_all(args,c,predictions,source,scope,receipt)
        elif receipt.get("unavailable_methods"):
            final_index(scope, None, receipt)
            result = {}
        else:
            abm.json_write(scope/"C1_FINAL_INDEX.json", dict(version=VERSION,status="completed",**receipt))
            result = {}
        if receipt.get("unavailable_methods"):
            result = dict(result, scope_status=PARTIAL_COMPLETION,
                unavailable_methods=list(receipt["unavailable_methods"]), evaluation_performed=bool(predictions))
        if receipt.get("reliability",{}).get("status","completed_primary_P1_vs_P0") != "completed_primary_P1_vs_P0":
            result = dict(result,scope_status=PARTIAL_RELIABILITY,
                reliability_status=receipt["reliability"]["status"])
        return result


def scope_support_preflight(args, plan, base):
    """Reject unsupported whole scopes before associations/PGS; never open test truth."""
    scope = args.analysis_root / args.Y / args.biom / "c1_correlate"
    args._c1_stage = "scope_support_preflight"
    configs = []
    for name in plan["methods"]:
        spec = get_method(name)
        configs.append(build_abm_config(name if spec.engine in {"twin", "prevalent_transfer", "torch_transformer"}
            else "abm_case_preserving", base, overrides=plan["overrides"].get(name, {}),
            trait=args.Y, layer=args.biom, seed=args.seed, cores=args.cores, final_multiroute_release=""))
    cfg = configs[0]
    if any((c.horizon, c.feature_missing, c.sample_missing) !=
           (cfg.horizon, cfg.feature_missing, cfg.sample_missing) for c in configs):
        raise ValueError("Methods must share the horizon and input QC contract")
    minimum = max(c.min_events for c in configs)
    required = {"build", "tune_model", "calibration_fit"}
    if any(get_method(name).engine == "twin" for name in plan["methods"]):
        required.add("tune_gate")
    if "glmnet" in plan["methods"]:
        required.add("calibration_audit")
    criteria = dict(version=1, horizon=cfg.horizon, feature_missing=cfg.feature_missing,
        sample_missing=cfg.sample_missing, technical_qc_scope=cfg.technical_qc_scope,
        minimum_per_class=minimum, required_roles=sorted(required), min_cases_with_biom=args.min_cases_with_biom)
    log = args.analysis_root / "logs/exclusions" / f"{args.Y}.{args.biom}.json"
    if log.is_file() and not args.prepared_source:
        saved = json.loads(log.read_text())
        if saved.get("criteria") == criteria and saved.get("status") == INSUFFICIENT_SCOPE_SUPPORT:
            args.include_prevalent = "abm_prevalent" in plan["methods"]
            host = load_module("le8_unified_host", HERE / "0.host.py")
            h = host.host_config(args)
            extra_files = {"objective_file": args.objective_file} if "molecular_deviation" in plan["methods"] else {}
            _, request = input_source_signature(h, extra_files)
            if normalize_cohort_request(saved.get("input_request", {})) == normalize_cohort_request(json.loads(json.dumps(request, default=str))):
                print(f"[LE8] SKIP C1 ENTIRE SCOPE Y={args.Y} biom={args.biom}: verified exclusion log; "
                      f"{saved['reason']}; no association/PGS/model fitting; continue", flush=True)
                return saved
    scope.mkdir(parents=True, exist_ok=True)
    source, contract, extra = prepare(args, plan, scope)
    args._prepared_inputs = (source, contract, extra)
    with np.load(source / "development_LOCAL_ONLY.npz", allow_pickle=False) as data:
        roles = data["role"].astype(str)
        time, event = np.asarray(data["time"], float), data["event"]
        if set(roles) != set(abm.ROLES) or time.shape != roles.shape or event.shape != roles.shape:
            raise ValueError("Invalid development roles or outcome lengths")
        if not np.isfinite(time).all() or (time <= 0).any() or not np.isin(event, [0, 1]).all():
            raise ValueError("Invalid development survival outcomes")

        def audit(qc, phase):
            rows = []
            for role in abm.ROLES:
                ix = (roles == role) & qc
                cases = int(np.sum(ix & (event == 1) & (time <= cfg.horizon)))
                controls = int(np.sum(ix & (time > cfg.horizon)))
                min_class = max(minimum, 20) if role == "calibration_audit" else minimum
                min_total = 100 if role == "calibration_audit" else 2 * min_class
                rows.append(dict(role=role, N=int(np.sum(roles == role)), technical_qc_pass=int(ix.sum()),
                    known_cases=cases, known_controls=controls, unknown_after_qc=int(ix.sum()) - cases - controls,
                    minimum_per_class=min_class, minimum_known_total=min_total, required=role in required,
                    sufficient=min(cases, controls) >= min_class and cases + controls >= min_total,
                    phase=phase))
            total_cases = int(np.sum(qc & (event == 1) & (time <= cfg.horizon)))
            rows.append(dict(role="all_development_with_biom", N=len(roles), technical_qc_pass=int(qc.sum()),
                known_cases=total_cases, known_controls=int(np.sum(qc & (time > cfg.horizon))),
                minimum_per_class=args.min_cases_with_biom, minimum_known_total=0,
                required=True, sufficient=total_cases >= args.min_cases_with_biom, phase=phase))
            return rows

        # Counts before QC are an upper bound: a failing scope needs no wide X load.
        rows = audit(np.ones(len(roles), bool), "before_technical_qc_upper_bound")
        if all(row["sufficient"] for row in rows if row["required"]):
            names = list(map(str, data["feature_names"]))
            clinical = contract.get("clinical_feature_names")
            if not clinical or not set(clinical) <= set(names):
                raise ValueError("Explicit matching clinical feature schema is required")
            cfg = replace(cfg, always_include_indices=tuple(names.index(v) for v in clinical))
            X = data["X"]
            prep = abm.NumericPrep().fit(X[roles == "build"], cfg)
            rows = audit(prep.technical_qc(X), "after_technical_qc")
    failures = [row for row in rows if row["required"] and not row["sufficient"]]
    receipt = dict(status=INSUFFICIENT_SCOPE_SUPPORT if failures else "supported", trait=args.Y, layer=args.biom,
        horizon=cfg.horizon, requested_methods=plan["methods"], development_support=rows,
        checked_before_association=True, test_truth_read=False,
        calibration_audit_used_for="support_counts_only", existing_outputs_preserved=True,
        prepared_source=str(source), criteria=criteria)
    pd.DataFrame(rows).to_csv(scope / "scope_support.tsv", sep="\t", index=False)
    if failures:
        reason = "; ".join(
            f"{row['role']}: cases={row['known_cases']}, minimum cases={row['minimum_per_class']}"
            if row["role"] == "all_development_with_biom" else
            f"{row['role']}: cases={row['known_cases']}, controls={row['known_controls']}, "
            f"minimum={row['minimum_per_class']}/class, known N>={row['minimum_known_total']}" for row in failures)
        receipt.update(reason=reason, association="skipped", pgs="skipped", prediction="skipped", evaluation="skipped")
        marker = source / "INPUTS_COMPLETE.json"
        if marker.is_file():
            request = json.loads(marker.read_text()).get("request", {})
            request = {k: v for k, v in request.items() if k not in {"version", "preparation_code"}}
            request.setdefault("input_contract_version", INPUT_CONTRACT_VERSION)
            receipt["input_request"] = request
        log.parent.mkdir(parents=True, exist_ok=True)
        abm.json_write(log, receipt)
        index = scope / "C1_FINAL_INDEX.json"
        previous = json.loads(index.read_text()) if index.is_file() else {}
        previous.update(status=INSUFFICIENT_SCOPE_SUPPORT, scope_preflight=str(scope / "C1_SCOPE_PREFLIGHT.json"),
                        reason=reason, existing_outputs_preserved=True)
        abm.json_write(index, previous)
        print(f"[LE8] SKIP C1 ENTIRE SCOPE Y={args.Y} biom={args.biom}: {reason}; "
              "association/PGS/models/evaluation skipped; continue to next trait/layer", flush=True)
    abm.json_write(scope / "C1_SCOPE_PREFLIGHT.json", receipt)
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
    plan, base = plan_and_base(args)
    if args.dry_run:
        print(json.dumps(public_plan(args,plan,base),indent=2)); return
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
            print(f"[LE8] START C1 Y={trait} biom={layer}: sample support -> association -> methods -> evaluation", flush=True)
            stage = "scope_support_preflight"
            try:
                with threadpool_limits(args.cores): support = scope_support_preflight(scoped, plan, base)
                if support["status"] == INSUFFICIENT_SCOPE_SUPPORT:
                    skipped.append(f"{trait}/{layer}")
                    continue
                stage = "association"
                scoped._c1_stage = stage
                native_association(scoped)
                stage = "methods/evaluation"
                with threadpool_limits(args.cores): result = run_scope(scoped,plan,base)
            except Exception as error:
                # A trait failure must not discard earlier results or prevent
                # the remaining traits from running. Ctrl+C still propagates.
                failed.append(f"{trait}/{layer}")
                stage = getattr(scoped, "_c1_stage", stage)
                # SKIP is forwarded immediately by the public console filter.
                print(f"[LE8] SKIP C1 Y={trait} biom={layer} ERROR stage={stage}: "
                      f"{type(error).__name__}: {error}; continue to next trait/layer; saved results retained", flush=True)
                traceback.print_exc()
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
              f"{', '.join(skipped)}. Early scope skips include association/PGS; "
              "see C1_SCOPE_PREFLIGHT.json and C1_FINAL_INDEX.json. Existing results retained.", flush=True)
    if partial:
        print("[LE8] PARTIAL C1 SUMMARY: unavailable methods or primary reliability: " + "; ".join(partial), flush=True)
    if failed:
        print("[LE8] PARTIAL C1 SUMMARY: failed scopes: " + ", ".join(failed) +
              ". All remaining trait/layer scopes were attempted; resume with --replace FALSE.", flush=True)
    from c1_report import write_reliability_cross_trait_summary, write_pgs_cross_trait_summary
    for layer in args.biom.split(","):
        write_reliability_cross_trait_summary(args.analysis_root, layer)
        write_pgs_cross_trait_summary(args.analysis_root, layer)
    # Signal failures only after the entire batch, so callers can detect them.
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
