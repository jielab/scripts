#!/usr/bin/env python3
"""Optional ABM-internal learning curve with nested, group-preserving build budgets.

Default: inspect development roles and print a plan; no fitting or test-truth read.
  python f/c1_training_budget.py --source /path/to/abm_twin --fractions .5,.75,1
Explicit execution, always in a fresh output directory:
  python f/c1_training_budget.py --source /path/to/abm_twin --fractions .5,.75,1 \
      --out /path/to/new_budget_experiment --run

Prepared files are development_LOCAL_ONLY.npz, test_X_LOCAL_ONLY.npz and
test_truth_LOCAL_ONLY.npz. The source's fit/MODEL_FROZEN.json supplies the resolved
configuration. With --source, --config changes algorithm settings while preserving
the source trait, layer, horizon, grouping and mandatory clinical columns. Empty
template clinical indices inherit the source schema and are recorded in the plan;
conflicting nonempty indices or explicit non-template contracts are rejected.
Individual paths can replace --source when their resolved config is provided.

Each budget refits the existing pipeline, including build-only preprocessing and
feature selection. Non-build roles and configuration stay fixed. All predictions
are saved before any test outcome is read. Evaluation is all-only on one common
QC/finite roster, with one predeclared largest-budget build censoring model.
Existing development-only gates are fitted by the reused pipeline; this script
does not tune a selector or use selected/fallback policies in its learning curve.

This controls the offered build sample for internal elastic-net/LightGBM/twin
arms. Effective eligible/pruned training sizes are reported separately. It does
not establish equal clinical features or training budgets with external C1 Cox
models. Metrics are descriptive; no C-index division by sample size is performed.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
import gc
import hashlib
import json
import math
import os
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

import abm
from c1_lock_roster import eid_hash

VERSION = "2026-10-21.abm-training-budget-v1"
DEVELOPMENT_FIELDS = {"X", "time", "event", "eid", "group", "role", "feature_names"}
FEATURE_FIELDS = {"X", "eid", "group", "feature_names"}
INPUT_CONTRACT_FIELDS = ("trait", "layer", "horizon", "grouping_unit", "always_include_indices")


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_fractions(values):
    fractions = [float(v) for v in values.split(",")] if isinstance(values, str) else list(map(float, values))
    if (not fractions or len(set(fractions)) != len(fractions) or
            any(not np.isfinite(v) or not 0 < v <= 1 for v in fractions)):
        raise ValueError("Fractions must be distinct finite numbers in (0, 1]")
    return tuple(sorted(fractions))


def _config_values(path):
    document = json.loads(Path(path).read_text())
    if not isinstance(document, dict) or not isinstance(document.get("config", document), dict):
        raise ValueError("Configuration must be a JSON object")
    return dict(document.get("config", document)), "config" in document


def resolve_config_contract(path, source=None):
    """A prepared source owns input identity; an algorithm template cannot reset it."""
    path = Path(path)
    override, frozen_override = _config_values(path)
    values = dict(override)
    audit = dict(override_config_path=str(path), override_config_sha256=sha256_file(path),
                 override_kind="standalone_resolved_config", inherited_fields=[], template_placeholders_replaced={})
    if source is not None:
        source_path = Path(source) / "fit/MODEL_FROZEN.json"
        original, _ = _config_values(source_path)
        missing = set(INPUT_CONTRACT_FIELDS) - set(original)
        if missing:
            raise ValueError("Source frozen configuration lacks input-contract fields: " + ", ".join(sorted(missing)))
        packaged_templates = {Path(__file__).resolve().parents[1] / "config" / name
                              for name in ("abm.json", "abm_balanced_reference.json")}
        template = (not frozen_override and
                    (path.resolve() in packaged_templates or
                     override.get("reference_experiment") == "balanced_reference_20261021"))
        same_source = path.resolve() == source_path.resolve()
        audit.update(source_config_path=str(source_path), source_config_sha256=sha256_file(source_path),
                     override_kind="source_frozen_configuration" if same_source else
                                   "algorithm_template" if template else "explicit_configuration")
        # These are placeholders in the distributed templates, never identities
        # of a prepared source. Non-placeholder conflicting values still fail.
        placeholders = dict(trait="cvd_cad", layer="prot", horizon=10., grouping_unit="family")
        values = {**original, **override}
        for field in INPUT_CONTRACT_FIELDS:
            inherited = field not in override
            if field == "always_include_indices":
                source_value = tuple(original[field])
                override_value = tuple(override.get(field, ()))
                if override_value and override_value != source_value:
                    raise ValueError("always_include_indices conflicts with the prepared source clinical schema")
                inherited = inherited or not override_value
                if inherited and source_value:
                    audit["clinical_schema_resolution"] = "empty_or_missing_indices_inherited_from_source"
            else:
                source_value = original[field]
                override_value = override.get(field, source_value)
                if override_value != source_value:
                    if not template or override_value != placeholders[field]:
                        raise ValueError(field + " conflicts with the prepared source input contract; --config only overrides algorithm settings")
                    inherited = True
                    audit["template_placeholders_replaced"][field] = dict(template=override_value, source=source_value)
            values[field] = source_value
            if inherited:
                audit["inherited_fields"].append(field)
        audit["algorithm_overrides"] = sorted(key for key, value in override.items()
                                               if key not in INPUT_CONTRACT_FIELDS and value != original.get(key))
    for key in ("coverages", "c_grid", "always_include_indices"):
        if key in values:
            values[key] = tuple(values[key])
    config = abm.Config(**values)
    config.validate()
    audit["resolved"] = {field: list(getattr(config, field)) if field == "always_include_indices" else getattr(config, field)
                         for field in INPUT_CONTRACT_FIELDS}
    return config, audit


def resolved_config(path, source=None):
    return resolve_config_contract(path, source)[0]


def load_development(path, *, include_features=False):
    # NpzFile is lazy: a dry run does not decompress the potentially large X.
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != DEVELOPMENT_FIELDS:
            raise ValueError("Development NPZ must contain exactly " + ", ".join(sorted(DEVELOPMENT_FIELDS)))
        out = {key: data[key] for key in data.files if include_features or key != "X"}
    ids, groups, roles = (np.asarray(out[key], str) for key in ("eid", "group", "role"))
    n = len(ids)
    invalid = {"", "nan", "none", "na", "<na>"}
    if (not n or len(set(ids)) != n or any(len(out[key]) != n for key in ("time", "event", "group", "role")) or
            any(v.strip().lower() in invalid for v in np.r_[ids, groups])):
        raise ValueError("Invalid development identifiers or lengths")
    if (set(roles) != set(abm.ROLES) or not np.isfinite(out["time"]).all() or
            np.any(out["time"] <= 0) or not np.isin(out["event"], [0, 1]).all()):
        raise ValueError("Development must contain valid survival outcomes and exactly the five original roles")
    role_by_group = {}
    for group, role in zip(groups, roles):
        if group in role_by_group and role_by_group[group] != role:
            raise ValueError("Family/group crosses development roles")
        role_by_group[group] = role
    names = np.asarray(out["feature_names"], str)
    if not len(names) or len(set(names)) != len(names):
        raise ValueError("Invalid feature schema")
    if include_features and (out["X"].shape != (n, len(names)) or np.isinf(out["X"]).any()):
        raise ValueError("Invalid development feature matrix")
    out.update(eid=ids, group=groups, role=roles, feature_names=names)
    return out


def budget_plan(development, fractions, seed, horizon):
    """Select outcome-blind hash prefixes of whole build groups, preserving rows."""
    fractions = parse_fractions(fractions)
    roles, groups, ids = (development[key] for key in ("role", "group", "eid"))
    build = roles == "build"
    families, counts = np.unique(groups[build], return_counts=True)
    order = sorted(range(len(families)), key=lambda j: (
        hashlib.sha256(f"{VERSION}:{seed}:{families[j]}".encode()).digest(), families[j]))
    ordered_groups = families[order]
    cumulative = np.cumsum(counts[order])
    role_rows = []
    for role in abm.ROLES:
        mask = roles == role
        role_rows.append(dict(role=role, N=int(mask.sum()), group_N=len(set(groups[mask])), eid_hash=eid_hash(ids[mask])))
    budgets, masks = [], {}
    for i, fraction in enumerate(fractions):
        target = int(math.ceil(fraction * int(build.sum())))
        group_count = int(np.searchsorted(cumulative, target, side="left") + 1)
        chosen = np.isin(groups, ordered_groups[:group_count]) & build
        keep = ~build | chosen
        key = f"b{i + 1:02d}_f{fraction:.12g}".replace(".", "p")
        known = (development["time"] > horizon) | ((development["event"] == 1) & (development["time"] <= horizon))
        cases = (development["event"] == 1) & (development["time"] <= horizon)
        budgets.append(dict(budget=key, target_fraction=fraction, target_build_N=target,
            offered_build_N=int(chosen.sum()), offered_build_fraction=float(chosen.sum() / build.sum()),
            offered_build_group_N=group_count, offered_build_known_N=int(np.sum(chosen & known)),
            offered_build_cases=int(np.sum(chosen & cases)), offered_build_eid_hash=eid_hash(ids[chosen]),
            extra_development_N=int((~build).sum()), total_development_N=int(keep.sum()),
            extra_fitting_role_N=int(np.sum((roles != "build") & (roles != "calibration_audit"))),
            extra_development_eid_hash=eid_hash(ids[~build]),
            **{role + "_N": int(np.sum(roles == role)) for role in abm.ROLES if role != "build"}))
        masks[key] = keep
    return dict(version=VERSION, budget_seed=int(seed), horizon=float(horizon),
        roles=role_rows, budgets=budgets,
        sampling="fixed hash order of whole build groups; shortest prefix reaching ceil(fraction * original build N)",
        group_size_rounding="actual N may exceed target N; groups are never split",
        evaluation_censoring_budget=budgets[-1]["budget"],
        unchanged_roles=list(abm.ROLES[1:]),
        scope="ABM-internal equal offered build budget and input schema; not an external C1 clinical-feature match",
        test_truth_read=False, models_fitted=False), masks


def subset_development(development, keep):
    return {key: value.copy() if key == "feature_names" else value[keep].copy()
            for key, value in development.items()}


def load_features(path, development):
    with np.load(path, allow_pickle=False) as source:
        if set(source.files) != FEATURE_FIELDS:
            raise ValueError("Projection NPZ must contain only X, eid, group, feature_names; no outcomes")
        features = {key: source[key] for key in source.files}
    X, ids, groups = abm.TwinABM.validate_features(features["X"], features["eid"], features["group"])
    if not np.array_equal(features["feature_names"].astype(str), development["feature_names"]):
        raise ValueError("Development/test feature schemas differ")
    if X.shape[1] != len(development["feature_names"]):
        raise ValueError("Test feature column count differs")
    if set(ids) & set(development["eid"]) or set(groups) & set(development["group"]):
        raise ValueError("Development/test identifiers or families overlap")
    features.update(X=X, eid=ids, group=groups)
    return features


def read_truth(path, ids):
    with np.load(path, allow_pickle=False) as source:
        if set(source.files) != {"eid", "time", "event"}:
            raise ValueError("Test truth NPZ must contain exactly eid, time, event")
        truth = {key: source[key] for key in source.files}
    if (not np.array_equal(np.asarray(truth["eid"], str), ids) or
            np.asarray(truth["time"]).shape != ids.shape or np.asarray(truth["event"]).shape != ids.shape or
            not np.isfinite(truth["time"]).all() or np.any(truth["time"] <= 0) or
            not np.isin(truth["event"], [0, 1]).all()):
        raise ValueError("Invalid test-truth roster or outcomes")
    return truth


def risk_arms(prediction, config):
    arms = tuple(sorted(name for name in prediction["probability"]
                        if name != "fallback_policy" and not name.startswith("twin_molecular_")))
    required = {"elasticnet_full", "twin_full", "twin_clean"}
    if config.lightgbm:
        required.add("lightgbm_full")
    if not required <= set(arms):
        raise ValueError("A required internal risk arm is missing")
    return arms


def training_arm_counts(model, development, arms):
    """Audit effective fitting counts; do not equate pruned and offered sizes."""
    build = development["role"] == "build"
    y, w = model.km.labels_weights(development["time"][build], development["event"][build])
    base = (w > 0) & model.prep.technical_qc(development["X"][build])
    donor = model.donor_table
    if not np.array_equal(np.asarray(donor["eid"], str), development["eid"][build]):
        raise ValueError("Donor/build roster mismatch")
    eligible = np.asarray(donor["donor_eligible"], bool)
    masks = dict(full=base, clean=base & donor["retained"], soft=base & (donor["soft_weight"] > 0),
                 random=base & donor["random_retained"])
    rows = {}
    for arm in arms:
        if arm.startswith(("elasticnet_", "lightgbm_")):
            mask = masks[arm.rsplit("_", 1)[1]]
            scope = "global model rows with positive training weight"
        elif arm.startswith("twin_"):
            suffix = arm.rsplit("_", 1)[1]
            mask = eligible if suffix == "full" else eligible & donor["retained" if suffix == "clean" else "random_retained"]
            scope = "reference-library donors"
        elif arm == "knn_unweighted_clean":
            mask, scope = eligible & donor["retained"], "reference-library donors"
        else:
            mask, scope = eligible, "eligible build rows for alternative risk route"
        rows[arm] = dict(effective_training_N=int(mask.sum()), effective_training_cases=int(y[mask].sum()),
            effective_training_controls=int(np.sum(mask & (y == 0))), training_count_scope=scope)
    return rows


def evaluate_frozen(entries, truth, features, km, out, plan):
    """All budgets and arms share one outcome-blind QC/finite evaluation mask."""
    ids = features["eid"]
    contract = plan.get("input_contract", {}).get("resolved", {})
    common_qc, common_finite = np.ones(len(ids), bool), np.ones(len(ids), bool)
    predictions = {}
    expected_arms = tuple(entries[0]["arms"])
    for entry in entries:
        path = out / entry["prediction_file"]
        if sha256_file(path) != entry["prediction_sha256"]:
            raise ValueError("Frozen prediction file changed before evaluation")
        pred = joblib.load(path)
        if (not np.array_equal(np.asarray(pred["eid"], str), ids) or
                not np.array_equal(np.asarray(pred["group"], str), features["group"]) or
                pred["frozen_model_hash"] != entry["model_hash"] or tuple(entry["arms"]) != expected_arms):
            raise ValueError("Budget prediction roster, model hash, or risk arms differ")
        qc = np.asarray(pred["diagnostics"]["technical_qc"], bool)
        if qc.shape != ids.shape:
            raise ValueError("Invalid saved QC mask")
        common_qc &= qc
        for arm in expected_arms:
            p = np.asarray(pred["probability"][arm], float)
            if p.shape != ids.shape or np.any((p[np.isfinite(p)] < 0) | (p[np.isfinite(p)] > 1)):
                raise ValueError("Invalid saved probability")
            common_finite &= np.isfinite(p)
        predictions[entry["budget"]] = pred
    common = common_qc & common_finite
    rows = []
    for entry in entries:
        budget = next(row for row in plan["budgets"] if row["budget"] == entry["budget"])
        for arm in expected_arms:
            result = abm.metrics(truth["time"][common], truth["event"][common],
                                 predictions[entry["budget"]]["probability"][arm][common], km)
            rows.append(dict(**budget, **entry["training_counts"][arm], **result, model=arm, subset="all_common",
                trait=contract.get("trait", ""), layer=contract.get("layer", ""),
                grouping_unit=contract.get("grouping_unit", "individual_unverified"),
                requested_test_N=len(ids), evaluated_N=int(common.sum()), excluded_N=int((~common).sum()),
                excluded_QC_N=int((~common_qc).sum()), excluded_nonfinite_after_QC_N=int(np.sum(common_qc & ~common_finite)),
                requested_test_eid_hash=eid_hash(ids), comparison_eid_hash=eid_hash(ids[common]),
                same_people=True, selection_policy_used=False,
                evaluation_scope="previously_inspected_or_unverified_exploratory",
                evaluation_censoring_budget=plan["evaluation_censoring_budget"],
                evaluation_weights="same frozen largest-budget build-derived censoring KM for every budget and arm",
                IPCW_horizon=km.horizon, inference_status="descriptive_no_bootstrap",
                feature_selection="refit within each budget build; input schema and mandatory clinical columns fixed",
                selected_panel_N=entry["selected_panel_N"], selected_panel_sha256=entry["selected_panel_sha256"]))
    frame = pd.DataFrame(rows)
    frame.to_csv(out / "training_budget_performance.tsv", sep="\t", index=False)
    return frame


def run_experiment(paths, output, config, plan, masks):
    output = Path(output)
    if output.exists():
        raise FileExistsError("Use a fresh output directory; budget experiments are never overwritten")
    development = load_development(paths["development"], include_features=True)
    current_plan, current_masks = budget_plan(development, [row["target_fraction"] for row in plan["budgets"]],
                                              plan["budget_seed"], config.horizon)
    if any(plan.get(key) != value for key, value in current_plan.items()) or any(not np.array_equal(masks[key], current_masks[key]) for key in current_masks):
        raise ValueError("Development inputs changed after budget planning")
    features = load_features(paths["test_features"], development)
    # This hashes inputs needed for fitting/projection only, never test truth.
    source_keys = ("development", "test_features", "config") + (("source_config",) if "source_config" in paths else ())
    sources = {key: dict(path=str(paths[key]), sha256=sha256_file(paths[key])) for key in source_keys}
    if "input_contract" in plan:
        audit = plan["input_contract"]
        actual = {field: list(getattr(config, field)) if field == "always_include_indices" else getattr(config, field)
                  for field in INPUT_CONTRACT_FIELDS}
        if actual != audit["resolved"]:
            raise ValueError("Resolved input contract changed after planning")
        if sources["config"]["sha256"] != audit["override_config_sha256"]:
            raise ValueError("Algorithm configuration changed after planning")
        if "source_config" in sources and sources["source_config"]["sha256"] != audit["source_config_sha256"]:
            raise ValueError("Source input contract changed after planning")
    output.mkdir(parents=True, exist_ok=False)
    frozen_plan = dict(plan, config=asdict(config), sources=sources,
        test_feature_eid_hash=eid_hash(features["eid"]),
        mandatory_clinical_features=[str(development["feature_names"][i]) for i in config.always_include_indices],
        code_sha256={name: sha256_file(Path(__file__).with_name(name)) for name in ("abm.py", Path(__file__).name)},
        environment={key: os.environ.get(key) for key in ("ABM_FINAL_SELECTOR", "ABM_MOE_CLUSTERS")})
    abm.json_write(output / "PLAN_FROZEN.json", frozen_plan)
    pd.DataFrame(plan["budgets"]).to_csv(output / "training_budget_plan.tsv", sep="\t", index=False)
    entries, evaluation_km = [], None
    for budget in plan["budgets"]:
        directory = output / budget["budget"]
        directory.mkdir()
        subset = subset_development(development, masks[budget["budget"]])
        data_path = directory / "development_LOCAL_ONLY.npz"
        np.savez_compressed(data_path, **subset)
        with threadpool_limits(config.cores):
            # A fresh fit sees only this budget's build rows. No full-build
            # normalizer, selected panel, coefficients or fitted model is reused.
            model = abm.fit_command(data_path, directory / "fit", config)
            model_digest = abm.model_hash(model)
            prediction = model.predict(features["X"], features["eid"], features["group"], features["feature_names"])
            prediction = abm.apply_predictability_extension(model, prediction, features["X"], features["eid"], features["group"])
        if abm.model_hash(model) != model_digest:
            raise ValueError("Projection changed the frozen fitted model")
        arms = risk_arms(prediction, config)
        if entries and arms != tuple(entries[0]["arms"]):
            raise ValueError("Risk-arm availability differs across budgets; evaluation aborted")
        prediction_path = directory / "test_predictions_LOCAL_ONLY.joblib"
        joblib.dump(prediction, prediction_path, compress=3)
        panel = list(map(str, model.panel_names))
        entry = dict(budget=budget["budget"], model_hash=model_digest, arms=list(arms),
            prediction_file=str(prediction_path.relative_to(output)), prediction_sha256=sha256_file(prediction_path),
            training_counts=training_arm_counts(model, subset, arms),
            selected_panel_N=len(panel), selected_panel_sha256=hashlib.sha256("\n".join(panel).encode()).hexdigest())
        entries.append(entry)
        abm.json_write(directory / "PREDICTIONS_FROZEN.json", entry)
        if budget["budget"] == plan["evaluation_censoring_budget"]:
            evaluation_km = model.km
        del model, prediction, subset
        gc.collect()
    # This marker precedes the first open of test truth. Failed/incomplete budget
    # runs never evaluate a subset of the requested experiments.
    abm.json_write(output / "ALL_PREDICTIONS_FROZEN.json", dict(version=VERSION, entries=entries, test_truth_read=False))
    truth = read_truth(paths["test_truth"], features["eid"])
    truth_hash = sha256_file(paths["test_truth"])
    frame = evaluate_frozen(entries, truth, features, evaluation_km, output, plan)
    abm.json_write(output / "COMPLETE.json", dict(version=VERSION, status="completed", budget_N=len(entries),
        all_predictions_frozen_before_test_truth=True, test_truth_sha256=truth_hash,
        performance_sha256=sha256_file(output / "training_budget_performance.tsv"),
        scope=plan["scope"], evaluation="all-only common QC/finite roster; fixed shared censoring model",
        test_provenance="previously_inspected_or_unverified_exploratory",
        external_validation=False, clinical_release=False, confidence_intervals="not computed; descriptive learning curve"))
    return frame


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--source", type=Path, help="Existing abm_twin directory containing prepared NPZ files")
    p.add_argument("--development", type=Path)
    p.add_argument("--test-features", type=Path)
    p.add_argument("--test-truth", type=Path)
    p.add_argument("--config", type=Path, help="Algorithm template with --source; otherwise an explicit resolved input config")
    p.add_argument("--fractions", default=".5,.75,1")
    p.add_argument("--budget-seed", type=int, help="Outcome-blind family sampling seed; defaults to fitting config seed")
    p.add_argument("--out", type=Path, help="Fresh experiment directory, required only with --run")
    action = p.add_mutually_exclusive_group()
    action.add_argument("--dry-run", action="store_true", help="Default: print plan without fitting or reading test truth")
    action.add_argument("--run", action="store_true", help="Explicitly perform all requested fresh fits and final evaluation")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    names = dict(development="development_LOCAL_ONLY.npz", test_features="test_X_LOCAL_ONLY.npz",
                 test_truth="test_truth_LOCAL_ONLY.npz", config="fit/MODEL_FROZEN.json")
    paths = {key: getattr(args, key) or (args.source / name if args.source else None) for key, name in names.items()}
    if paths["development"] is None or paths["config"] is None:
        p.error("Provide --source or both --development and --config")
    if args.run and (args.out is None or any(path is None for path in paths.values())):
        p.error("--run requires --out, test features and test truth paths")
    config, contract = resolve_config_contract(paths["config"], args.source)
    development = load_development(paths["development"])
    if args.source is not None:
        paths["source_config"] = args.source / "fit/MODEL_FROZEN.json"
        with np.load(args.source / "development_LOCAL_ONLY.npz", allow_pickle=False) as original:
            if not np.array_equal(np.asarray(original["feature_names"], str), development["feature_names"]):
                raise ValueError("Prepared feature schema differs from the source clinical input contract")
    if any(i < 0 or i >= len(development["feature_names"]) for i in config.always_include_indices):
        raise ValueError("Mandatory clinical column indices do not match the prepared schema")
    plan, masks = budget_plan(development, args.fractions, config.seed if args.budget_seed is None else args.budget_seed, config.horizon)
    contract["mandatory_clinical_features"] = [str(development["feature_names"][i]) for i in config.always_include_indices]
    plan["input_contract"] = contract
    if not args.run:
        print(json.dumps(plan, ensure_ascii=False, indent=2, allow_nan=False))
        return plan
    return run_experiment(paths, args.out, config, plan, masks)


if __name__ == "__main__":
    main()
