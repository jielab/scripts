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

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import abm
from c1_registry import (VERSION as REGISTRY_VERSION, get_method, load_method_plan,
                        parse_methods, build_abm_config, resolve_device, scientific_fingerprint)
from c1_unified_models import context, fit_glmnet, fit_lightgbm, fit_twin, fit_prevalent, fit_deep

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
VERSION = "LE8-C1-FINAL-2026-10-22-v1"
UPSTREAM_COMMIT = "40a884324f30219f865e8a7a49ac9708a913e742"


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


def code_fingerprint():
    paths = [HERE/n for n in ("abm.py", "0.host.py", "c1_registry.py", "c1_correlate.py",
        "c1_unified_models.py", "c1_prevalent.py", "c1_unified_glmnet.R", "c1_abm_gpu.py",
        "0.common.py", "c1_lock_roster.py")]
    return hashlib.sha256(json.dumps({p.name: sha_file(p) for p in paths}, sort_keys=True).encode()).hexdigest()


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
    p.add_argument("--prepared-source", type=Path, help="Reuse a frozen prepared input directory (one trait/layer only)")
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--cores", type=int, default=8)
    p.add_argument("--replace", choices=["TRUE", "FALSE"], default="FALSE")
    p.add_argument("--r-bin", default="Rscript")
    p.add_argument("--host-args", default="")
    p.add_argument("--preflight", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    for flag in ("ukb-phe", "ukb-rel-file", "group-file", "group-col", "outer-roster", "end-date",
                 "diagnosis-col", "shared-covariates", "deviation-config", "objective-file", "pgs-root", "pgs-mapping"):
        p.add_argument("--"+flag)
    p.add_argument("--remerge-pgs", action="store_true")
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
    if args.bootstrap < 1 or args.cores < 1:
        raise ValueError("Positive bootstrap count and core count are required")
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
    return dict(version=VERSION, upstream_commit=UPSTREAM_COMMIT, traits=args.Y.split(","),
        layers=args.biom.split(","), methods=methods, association=args.association,
        common_evaluation="same test people, outcome horizon, build censoring model and paired family draws",
        test_provenance="previously_inspected_or_unverified_exploratory",
        proposed_outputs="<analysis-root>/<Y>/<biom>/c1_correlate/{methods,evaluation}")


def dependencies(args, plan):
    required = set(plan["methods"])
    if "glmnet" in required:
        subprocess.run([args.r_bin, "-e", "p<-c('glmnet','jsonlite');m<-p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)];if(length(m))stop(paste(m,collapse=', '))"], check=True)
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
        from c1_prevalent import PrevalentConfig
        PrevalentConfig(**json.loads(args.history_config.read_text())).validate()


def native_association(args):
    demo = "--demo" in shlex.split(args.host_args)
    if args.association == "skip" or (args.association == "auto" and (args.prepared_source or demo)):
        return {"status": "skipped_by_scope", "existing_outputs_preserved": True}
    cmd = ["bash", str(HERE/"0.engine.sh"), "c1_correlate", "--Y", args.Y, "--biom", args.biom,
           "--analysis-root", str(args.analysis_root), "--r-bin", args.r_bin,
           "--cores", str(args.cores), "--seed", str(args.seed), "--replace", args.replace]
    if args.ukb_phe: cmd += ["--ukb-phe", args.ukb_phe]
    env = os.environ.copy()
    env["LE8_C1_STAGE"] = "c1_correlate" # only association; prediction handled below
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
    payload = dict(version=VERSION, fields=fields, files=files,
        preparation_code={name:sha_file(HERE/name) for name in ("0.host.py", "0.common.py", "c1_lock_roster.py")})
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
    marker = source/"INPUTS_COMPLETE.json"
    if marker.is_file():
        saved = json.loads(marker.read_text())
        if saved.get("source_signature") == signature and all((source/n).is_file() and sha_file(source/n)==s for n,s in saved.get("files",{}).items()):
            print("[LE8 C1] Reusing verified prepared inputs", flush=True)
            return source, json.loads((source/"host_contract.json").read_text()), None
        raise RuntimeError("Prepared input cache integrity failed; keep evidence and rebuild with a new source contract")
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
    history = joblib.load(history_path) if history_path.is_file() else None
    grouping = contract.get("independent_family_support_verified", contract.get("config", {}).get("grouping_unit")=="family")
    return d, x, indices, history, "family" if grouping else "individual_unverified"


def completion(path, signature):
    try:
        value = json.loads((path/"COMPLETE.json").read_text())
        return value.get("signature") == signature and bool(value.get("files")) and all(
            (path/n).is_file() and sha_file(path/n)==s for n,s in value["files"].items())
    except (OSError, ValueError, KeyError): return False


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
    host.run_deviation(args, joblib.load(path), stage, dc, contract)


def final_index(scope, destination, run_receipt):
    abm.json_write(scope/"C1_FINAL_INDEX.json", dict(version=VERSION, evaluation=str(destination),
        methods=run_receipt["methods"], source_commit=UPSTREAM_COMMIT, status="completed",
        supplementary_historical_results="existing c1_abm and c1_abm_GPU directories retained read-only"))


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
    from c1_unified_evaluation import evaluate_methods
    from le8_common import write_tables
    # Every model is frozen before this boundary. Truth bytes identify evaluation
    # reuse only; they are deliberately absent from model-training signatures.
    run_receipt = dict(run_receipt, test_truth_sha256=sha_file(source/"test_truth_LOCAL_ONLY.npz"))
    sig = hashlib.sha256(json.dumps(dict(run_receipt, bootstrap=args.bootstrap,
        evaluator=sha_file(HERE/"c1_unified_evaluation.py")), sort_keys=True, default=str).encode()).hexdigest()
    destination = scope/"evaluation"/sig[:20]
    destination.parent.mkdir(parents=True, exist_ok=True)
    if args.replace != "TRUE" and completion(destination, sig):
        print("[LE8 C1] Reuse common evaluation: verified truth, models, evaluator and aggregates", flush=True)
        final_index(scope, destination, run_receipt)
        return read_public_tables(destination)
    # First semantic test-truth read, after all selected methods were frozen.
    truth = abm.read_npz(source/"test_truth_LOCAL_ONLY.npz")
    tables = evaluate_methods(predictions, truth, c.km, bootstrap=args.bootstrap, seed=args.seed,
        grouping_unit=c.cfg.grouping_unit, primary_coverage=c.cfg.primary_coverage)
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
            primary_coverage=c.cfg.primary_coverage, all_selected_models_frozen_before_test_truth=True,
            test_provenance="previously_inspected_or_unverified_exploratory",
            conditional_inference="bootstrap conditional on fitted models; no full retraining uncertainty",
            source_commit=UPSTREAM_COMMIT))
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
        source, contract, _ = prepare(args, plan, scope)
        d, x, clinical, history, grouping = read_inputs(source, contract)
        cfg = build_abm_config("abm_case_preserving", base, trait=args.Y, layer=args.biom,
            seed=args.seed, cores=args.cores, parallel_workers=min(4,args.cores), bootstrap=args.bootstrap,
            grouping_unit=grouping, always_include_indices=clinical, final_multiroute_release="")
        c = context(d, x, cfg)
        source_fingerprint = {name:sha_file(source/name) for name in
            ("development_LOCAL_ONLY.npz", "test_X_LOCAL_ONLY.npz", "host_contract.json")}
        if history is not None: source_fingerprint["history"] = sha_file(source/"prevalent_source_LOCAL_ONLY.joblib")
        code = code_fingerprint()
        predictions, receipt = {}, dict(methods={}, source=source_fingerprint, code=code,
            trait=args.Y, layer=args.biom, horizon=cfg.horizon, grouping_unit=grouping)
        history_config = json.loads(args.history_config.read_text()) if args.history_config else None
        for name in plan["methods"]:
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
                extra["deviation_code"] = sha_file(HERE/"deviation.py")
                for label, path in {"pgs_source":host.pgs_path(args), "pgs_mapping":args.pgs_mapping,
                                    "deviation_input":source/"deviation_LOCAL_ONLY.joblib"}.items():
                    if path:
                        path = Path(path)
                        extra[label] = {"path":str(path), "sha256":sha_file(path) if path.is_file() else "missing"}
                extra["remerge_pgs"] = bool(args.remerge_pgs)
            signature = scientific_fingerprint(name, {**asdict(method_cfg), **extra}, {**source_fingerprint, "code":code})
            target = scope/"methods"/name/(signature[:16]+"_"+device["actual"])
            target.parent.mkdir(parents=True, exist_ok=True)
            if args.replace != "TRUE" and completion(target,signature):
                print(f"[LE8 C1] Reuse {name}: verified frozen model and projections",flush=True)
            else:
                print(f"[LE8 C1] Fit {name}; device={device['actual']}; all shared test outcomes withheld",flush=True)
                with tempfile.TemporaryDirectory(prefix=".fit-",dir=target.parent) as temp:
                    stage = Path(temp)/"result"
                    if spec.engine != "deviation": stage.mkdir()
                    if name == "glmnet": fit_glmnet(method_context,stage,args.r_bin)
                    elif name == "lightgbm": fit_lightgbm(method_context,stage)
                    elif spec.engine == "twin": fit_twin(method_context,method_cfg,stage,name)
                    elif name == "abm_prevalent": fit_prevalent(method_context,history,stage,history_config)
                    elif name == "abm_deep":
                        resources = load_module("le8_c1_resources",HERE/"0.resources.py")
                        with resources.gpu_lease() if device["actual"]=="cuda" else nullcontext():
                            fit_deep(method_context,stage,extra["deep_config"],device["actual"])
                    elif name == "molecular_deviation":
                        fit_deviation(args, source, stage, grouping, contract)
                    else: raise AssertionError(name)
                    if code_fingerprint() != code:
                        raise RuntimeError("Source code changed during fitting; mixed-code output will not be published")
                    mark_complete(stage,signature,{"method":name,"device":device,"parameters":asdict(method_cfg)})
                    publish(stage,target)
            receipt["methods"][name] = dict(path=str(target),signature=signature,device=device)
            if spec.engine != "deviation":
                predictions[name] = joblib.load(target/"test_predictions_LOCAL_ONLY.joblib")
            if name == "abm_prevalent":
                predictions["abm_prevalent_control"] = joblib.load(target/"nohistory_control/test_predictions_LOCAL_ONLY.joblib")
                receipt["methods"]["abm_prevalent_control"] = dict(path=str(target/"nohistory_control"),
                    signature=signature,device=device,automatic_matched_ablation=True)
        if code_fingerprint() != code:
            raise RuntimeError("Source code changed during the run; evaluation requires frozen code")
        if predictions: return evaluate_all(args,c,predictions,source,scope,receipt)
        abm.json_write(scope/"C1_FINAL_INDEX.json",dict(version=VERSION,status="completed",**receipt))
        return {}


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
    native_association(args)
    for trait in args.Y.split(","):
        for layer in args.biom.split(","):
            scoped = argparse.Namespace(**vars(args)); scoped.Y,scoped.biom = trait,layer
            with threadpool_limits(args.cores): run_scope(scoped,plan,base)


if __name__ == "__main__":
    main()
