#!/usr/bin/env python3
"""Re-evaluate trusted frozen LE8 predictions in a fresh, aggregate-only directory.

No fitting, threshold optimization, phenotype rebuilding, or replacement of old
results occurs. An optional projection of the predeclared calibration_audit role
provides development-only bin boundaries for the support audit. Test predictions
are always reused verbatim. Missing private inputs produce an availability
report, never a substitute training-screening AUC.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import gzip
import json
from pathlib import Path
import re

import joblib
import numpy as np
import pandas as pd

import abm
from c1_support_audit import build_support_audit
from le8_common import sha256, write_tables

REPORT_VERSION = "2026-10-21.fair-evaluation-v1"
REQUIRED = (
    "fit/model_LOCAL_ONLY.joblib", "fit/MODEL_FROZEN.json",
    "test_predictions_LOCAL_ONLY.joblib", "test_truth_LOCAL_ONLY.npz",
    "host_contract.json",
)
INTERNAL_AGGREGATES = ("metrics", "contrasts", "prediction_finiteness", "calibration",
                       "confidence_validation", "strata_diagnostics")


def _read_json(path):
    return json.loads(Path(path).read_text())


def _aggregate_tables(tables):
    """Reject participant-level columns/cells before creating any report files."""
    checked = {}
    for name, values in tables.items():
        frame = pd.DataFrame(values)
        for column in frame.columns:
            label = str(column).lower()
            aggregate_hash = label.endswith(("_hash", "_sha256"))
            if not aggregate_hash and (label in {"id", "ids", "group", "groups", "family"} or re.search(
                    r"(?:^|_)(?:eids?|participant_ids?|sample_ids?|subject_ids?|family_ids?|group_ids?)(?:$|_)", label)):
                raise ValueError(f"Participant identifiers cannot be published: {name}/{column}")
            if frame[column].map(lambda x: isinstance(x, (dict, list, tuple, set, np.ndarray))).any():
                raise ValueError(f"Only scalar aggregate cells can be published: {name}/{column}")
        checked[name] = frame
    return checked


def _record_sources(sources, additions):
    """Never replace an earlier checksum when a component reloads a source."""
    for filename, digest in additions.items():
        key = str(Path(filename).resolve())
        if key in sources and sources[key] != digest:
            raise ValueError(f"Source changed between report components: {filename}")
        sources[key] = digest


def _publish(out, tables, metadata):
    """Only publish to a fresh directory, with a final completion receipt."""
    out = Path(out)
    if out.exists():
        raise FileExistsError(f"Choose a new report directory; will not replace: {out}")
    tables = _aggregate_tables(tables)
    for filename, digest in metadata.get("sources", {}).items():
        if sha256(Path(filename)) != digest:
            raise ValueError(f"Source changed during reporting: {filename}")
    out.mkdir(parents=True)
    write_tables(out, tables, "Frozen fair evaluation", metadata["layer"],
                 metadata["trait"], workbook="fair_evaluation.xlsx")
    abm.json_write(out / "report.json", metadata)
    (out / "README.txt").write_text(
        "This report reuses frozen predictions; no model was fitted.\n"
        "Harrell C uses strict ordered follow-up times; IPCW AUC uses the frozen horizon and KM.\n"
        "Compare models within a table/subset with identical N and roster hash.\n"
        "All versus selected changes the evaluated population, and is not an algorithm comparison.\n"
        "Risk-matched controls have identical counts within risk bins, not exact continuous-risk balance.\n"
        "Support count is not disease risk or a probability of correctness.\n"
        "Different reference sizes, geometries and radii must be reported when comparing counts.\n"
        "Confidence intervals condition on the frozen fit, selectors and censoring model.\n"
        "Unavailable components are not validated findings. No participant IDs are exported.\n",
        encoding="utf-8",
    )
    for filename, digest in metadata.get("sources", {}).items():
        if sha256(Path(filename)) != digest:
            raise ValueError(f"Source changed during reporting: {filename}")
    files = {p.name: sha256(p) for p in out.iterdir() if p.is_file()}
    abm.json_write(out / "COMPLETE.json", dict(status="completed", report_status=metadata["status"], files=files))
    return out


def _reference_prediction(model, development):
    """Project only the allowed audit role; do not relax model overlap checks."""
    if not development.is_file():
        return None, "unavailable_development_file"
    with np.load(development, allow_pickle=False) as data:
        required = {"role", "eid", "group", "X", "feature_names"}
        if not required <= set(data.files):
            raise ValueError("Development file lacks audit projection fields")
        ix = np.asarray(data["role"], str) == "calibration_audit"
        if not ix.any():
            return None, "unavailable_calibration_audit_role"
        # time/event are deliberately never accessed by this function.
        X, ids, groups, names = data["X"][ix], data["eid"][ix], data["group"][ix], data["feature_names"]
    before = abm.model_hash(model)
    pred = model.predict(X, ids, groups, names, allow_development_audit=True)
    pred = abm.apply_predictability_extension(model, pred, X, ids, groups)
    if (pred.get("evaluation_scope") != "calibration_audit" or abm.model_hash(model) != before
            or pred.get("frozen_model_hash") != before):
        raise ValueError("Audit projection changed the frozen model or used an invalid role")
    if (not np.array_equal(np.asarray(pred["eid"], str), np.asarray(ids, str)) or
            not np.array_equal(np.asarray(pred["group"], str), np.asarray(groups, str))):
        raise ValueError("Audit projection changed the requested participant/family roster")
    return pred, "projected_frozen_calibration_audit_without_outcomes"


def _training_tables(model, contract):
    audit = getattr(model, "training_audit", {})
    roles = [dict(role=k, **v) for k, v in audit.get("roles", {}).items()]
    composition = audit.get("pruning_composition")
    if not composition:
        known_cases, retained_cases = audit.get("known_cases"), audit.get("retained_cases")
        composition = dict(status="legacy_known_case_summary_eligibility_denominator_not_reconstructed",
                           known_build_cases=known_cases, retained_cases=retained_cases,
                           eligible_reference_N=audit.get("known_donors"), retained_reference_N=audit.get("retained_donors"),
                           retention_relative_to_known_build_cases=(retained_cases / known_cases
                               if known_cases and retained_cases is not None else np.nan),
                           interpretation="Known build cases may include QC-ineligible cases; this ratio is not attributed exclusively to pruning")
    return {
        "training_roles": pd.DataFrame(roles),
        "training_pruning": pd.DataFrame(audit.get("pruning", [])),
        "training_case_retention": pd.DataFrame([composition]),
        "training_comparability": pd.DataFrame([dict(
            status="information_and_training_budget_not_proven_matched_to_C1",
            same_heldout_people_does_not_imply_same_training_budget=True,
            ABM_molecular_budget=model.config.panel_size,
            ABM_selected_panel_names=json.dumps(getattr(model, "panel_names", [])),
            ABM_selected_molecular_names=json.dumps(audit.get("molecular_panel_features", [])),
            clinical_feature_names=json.dumps(contract.get("clinical_feature_names", [])),
            reference_experiment=getattr(model.config, "reference_experiment", "legacy"),
            radius_policy=getattr(model.config, "reference_radius_policy", "arm_calibrated"),
            twin_random_control=getattr(model.config, "twin_random_control", False),
            action="Use c1_training_budget.py for ABM internal common-input build-budget curves; broad C1 comparators remain separately labelled",
        )]),
    }


def run_report(root, trait, layer, output, *, source=None, bootstrap=1000,
               count_bootstrap=0, reference=True, include_gpu=False):
    if (isinstance(bootstrap, bool) or isinstance(count_bootstrap, bool) or
            not np.isfinite(bootstrap) or not np.isfinite(count_bootstrap) or
            int(bootstrap) != bootstrap or int(count_bootstrap) != count_bootstrap or
            bootstrap < 1 or count_bootstrap < 0):
        raise ValueError("bootstrap must be a positive integer; count_bootstrap must be a nonnegative integer")
    bootstrap, count_bootstrap = int(bootstrap), int(count_bootstrap)
    root, out = Path(root), Path(output)
    scope = root / trait / layer
    twin = Path(source) if source is not None else scope / "c1_abm" / "abm_twin"
    if out.exists():
        raise FileExistsError(f"Choose a new report directory; will not replace: {out}")
    # Keep new reports out of the frozen fit tree and its completion receipts.
    frozen_roots = (twin, scope / "c1_correlate", scope / "c1_abm_GPU")
    if any(out.resolve().is_relative_to(path.resolve()) for path in frozen_roots):
        raise ValueError("Report output must be outside frozen ABM/C1/GPU source directories")
    metadata = dict(version=REPORT_VERSION, trait=trait, layer=layer,
                    bootstrap=bootstrap, count_bootstrap=count_bootstrap,
                    training_performed=False, test_predictions_recomputed=False,
                    exploratory=True, sources={},
                    co_primary_metrics="Harrell_C full follow-up strict-time; AUC_IPCW frozen horizon",
                    secondary_metrics="Brier_IPCW; LogLoss_IPCW; Uno_C_horizon",
                    same_training_information_with_C1=False,
                    censoring="reuse frozen marginal training KM, including its original tie convention")
    missing = [name for name in REQUIRED if not (twin / name).is_file()]
    if missing:
        metadata.update(status="not_available", missing=missing)
        tables = {"availability": pd.DataFrame([dict(
            component="frozen_ABM", status="not_available", missing="; ".join(missing),
            explanation="Public aggregates or a preparation workbook cannot replace private frozen predictions")])}
        return _publish(out, tables, metadata)

    sources = {str((twin / name).resolve()): sha256(twin / name) for name in REQUIRED}
    marker, contract = _read_json(twin / REQUIRED[1]), _read_json(twin / "host_contract.json")
    model = abm.load_model(twin / REQUIRED[0])
    pred = joblib.load(twin / REQUIRED[2])  # trusted local LE8 artifact only
    truth = abm.read_npz(twin / REQUIRED[3])
    frozen_hash = abm.model_hash(model)
    if marker.get("hash") != frozen_hash or pred.get("frozen_model_hash") != frozen_hash:
        raise ValueError("Frozen model, marker and prediction hashes disagree")
    if model.config.trait != trait or model.config.layer != layer:
        raise ValueError("Requested trait/layer differs from frozen model")
    if not np.array_equal(np.asarray(pred["eid"], str), np.asarray(truth["eid"], str)):
        raise ValueError("Frozen prediction/truth IDs or their order differ")
    if len(set(map(str, pred["eid"]))) != len(pred["eid"]):
        raise ValueError("Duplicate frozen participant IDs")
    grouping = getattr(model.config, "grouping_unit", "individual_unverified")
    if contract.get("grouping_unit", grouping) != grouping:
        raise ValueError("Host and frozen model grouping provenance differ")
    print("[LE8 fair evaluation] Frozen fit/roster verified; evaluating saved predictions", flush=True)
    result = abm.evaluate(model, pred, truth["time"], truth["event"], truth["eid"],
                          bootstrap=bootstrap, full_curves=True, fit_calibration=False)
    tables = {"internal_" + name: pd.DataFrame(result[name]) for name in INTERNAL_AGGREGATES if name in result}
    tables["internal_status"] = pd.DataFrame([result["status"]])
    tables.update(_training_tables(model, contract))

    ref, ref_status = None, "disabled_by_user"
    if reference:
        development = twin / "development_LOCAL_ONLY.npz"
        if development.is_file():
            _record_sources(sources, {str(development): sha256(development)})
        print("[LE8 fair evaluation] Building development-only support reference; no refit", flush=True)
        ref, ref_status = _reference_prediction(model, development)
    tables.update(build_support_audit(model, pred, truth["time"], truth["event"], truth["eid"],
                                     reference_prediction=ref, bootstrap=count_bootstrap))
    availability = [dict(component="frozen_ABM", status="evaluated"),
                    dict(component="support_development_reference", status=ref_status)]
    benchmark_path = scope / "c1_correlate/prediction/_LOCAL_ONLY/benchmark.json.gz"
    if benchmark_path.is_file():
        from c1_abm_compare import compare
        _record_sources(sources, {str(benchmark_path): sha256(benchmark_path)})
        with gzip.open(benchmark_path, "rt") as stream:
            benchmark = json.load(stream)
        outer_hash = contract.get("outer_roster_sha256")
        if not outer_hash or benchmark.get("lockpaper_protocol", {}).get("outer_roster_sha256") != outer_hash:
            raise ValueError("C1/ABM outer-roster provenance differs or is unavailable")
        comparison = compare(benchmark, pred, truth, model.km, replace(model.config, bootstrap=bootstrap))
        tables.update({"C1_" + key: value for key, value in comparison.items()})
        availability.append(dict(component="C1_shared_test_comparison", status=str(comparison["status"].iloc[0]["status"])))
    else:
        availability.append(dict(component="C1_shared_test_comparison", status="not_available"))
    if include_gpu:
        from c1_abm_gpu_compare import load_inputs, paired_cohort, performance, primary_tables
        if source is not None and twin.resolve() != (scope / "c1_abm/abm_twin").resolve():
            raise ValueError("GPU comparison requires its original CPU source; do not mix a new experimental fit")
        if not (scope / "c1_abm_GPU/COMPLETE.json").is_file():
            availability.append(dict(component="GPU_shared_test_comparison", status="not_available"))
        else:
            gpu_root = scope / "c1_abm_GPU"
            transformer = gpu_root / "transformer" if (gpu_root / "transformer").is_dir() else gpu_root
            required_gpu = (twin / "COMPLETE.json", gpu_root / "host_contract.json",
                            transformer / "test_predictions_LOCAL_ONLY.tsv.gz", transformer / "MODEL_FROZEN.json")
            missing_gpu = [str(path.relative_to(scope)) for path in required_gpu if not path.is_file()]
            loaded = None
            if not missing_gpu:
                try:
                    loaded = load_inputs(scope)
                except FileNotFoundError as exc:
                    missing_gpu = [str(exc)]
            if missing_gpu:
                availability.append(dict(component="GPU_shared_test_comparison", status="not_available",
                                         detail="; ".join(missing_gpu)))
            else:
                gm, gp, gt, gpu, gc, _, gf, conventional = loaded
                if abm.model_hash(gm) != frozen_hash:
                    raise ValueError("GPU comparison loaded a different frozen CPU model")
                _record_sources(sources, gf)
                data = paired_cohort(gp, gt, gpu, gm.km, conventional)
                data.update(trait=trait, layer=layer, grouping_unit=grouping,
                            outer_roster_sha256=gc["outer_roster_sha256"])
                performance_rows, differences = performance(data, gm.km, bootstrap, model.config.seed)
                tables.update({"GPU_" + key: value for key, value in primary_tables(performance_rows, differences).items()})
                tables["GPU_cohort"] = pd.DataFrame([data["audit"]])
                availability.append(dict(component="GPU_shared_test_comparison", status="evaluated"))
    if abm.model_hash(model) != frozen_hash:
        raise ValueError("Reporting changed the frozen model")
    tables["availability"] = pd.DataFrame(availability)
    code_files = ("abm.py", "c1_fair_evaluation.py", "c1_support_audit.py", "c1_abm_compare.py",
                  "c1_abm_gpu_compare.py", "c1_selector_inference.py")
    metadata.update(status="evaluated", frozen_model_hash=frozen_hash, grouping_unit=grouping,
                    source_scope=str(twin), sources=sources, support_reference=ref_status,
                    components=availability, code_sha256={name: sha256(Path(__file__).with_name(name)) for name in code_files})
    return _publish(out, tables, metadata)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--trait", "--Y", dest="trait", required=True)
    parser.add_argument("--layer", "--biom", dest="layer", choices=("prot", "met"), required=True)
    parser.add_argument("--source", type=Path, help="Alternative trusted abm_twin source directory")
    parser.add_argument("--output", type=Path, required=True, help="Fresh directory outside frozen source")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--count-bootstrap", type=int, default=0,
                        help="Optional count-only coverage inference; 0 leaves descriptive count tables")
    parser.add_argument("--skip-reference", action="store_true",
                        help="Skip frozen development audit projection; non-fixed ordering bins then unavailable")
    parser.add_argument("--include-gpu", action="store_true", help="Also verify/re-evaluate saved GPU predictions")
    args = parser.parse_args()
    out = run_report(args.analysis_root, args.trait, args.layer, args.output, source=args.source,
                     bootstrap=args.bootstrap, count_bootstrap=args.count_bootstrap,
                     reference=not args.skip_reference, include_gpu=args.include_gpu)
    print(f"[LE8 fair evaluation] DONE: {out}", flush=True)


if __name__ == "__main__":
    main()
