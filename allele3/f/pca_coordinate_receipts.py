"""Publish and verify common-panel PCA coordinate receipts.

This module reads metadata and existing scores only; it never fits PCA or
changes scores. Source-basis identity and scored-coordinate-frame identity
remain separate, and reference metadata retains its existing provenance rule.
"""
from __future__ import annotations
import csv
import hashlib
import io
import json
from pathlib import Path
import re
from types import SimpleNamespace

import pca_project as pp
from analysis_integrity import identity

RECEIPT_SCHEMA = "allele3-pca-coordinate-receipt-1"
FRAME_SCHEMA = "allele3-pca-coordinate-frame-1"


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def artifact(path):
    value = identity(path, full=True)
    if value.get("status") == "missing":
        raise ValueError(f"Missing coordinate evidence: {path}")
    return {key: value[key] for key in ("path", "size", "sha256")}


def verify_artifact(saved):
    current = artifact(saved["path"])
    if current != saved:
        raise ValueError("Coordinate evidence changed: " + saved["path"])
    return Path(saved["path"])


def metadata_request(args, basis_plan):
    """Bind inputs without promoting an unverified historical reference table.

    A target metadata file is a new explicit input snapshot.  A reference table
    remains subject to the pre-existing frozen-basis/provenance contract.
    Include this return value in the projection plan signature.
    """
    from pca_contract import historical_identity
    original = basis_plan.get("inputs", {}).get("metadata")
    explicit_reference = getattr(args, "reference_metadata", None)
    reference = explicit_reference or (original or {}).get("path")
    provenance = getattr(args, "provenance_manifest", None)
    if reference and Path(reference).is_file():
        history = historical_identity(reference, provenance, original)
        verified = history["historical_metadata_verification"] == "verified"
        reference_binding = {
            "status": "verified_frozen_reference_metadata" if verified else "historically_unverified",
            "file": artifact(reference),
            # Only a prior verified digest can certify historical identity.
            "recorded": {"path": str(Path(reference).resolve()), "sha256": history["expected_sha256"]}
                        if verified else original,
            "historical_evidence": history,
        }
    elif explicit_reference:
        raise ValueError("Explicit --reference-metadata does not exist")
    else:
        reference_binding = {"status": "reference_metadata_not_available_here",
                             "file": None, "recorded": original}
    target = getattr(args, "target_metadata", None)
    target_snapshot = artifact(target) if target else None
    target_binding = ({"status": "explicit_new_target_metadata_snapshot", "file": target_snapshot,
                       "recorded": target_snapshot} if target else
                      {"status": "not_provided", "file": None, "recorded": None})
    return {"reference": reference_binding, "target": target_binding}


def _by_id(path, required):
    result = {}
    for row in pp.rows(path, required):
        if row["ID"] in result:
            raise ValueError("Duplicate prepared variant: " + row["ID"])
        result[row["ID"]] = row
    return result


def _weights(path, pcs):
    result = {}
    for row in pp.rows(path, {"ID", "A1", *pcs}):
        key = row["ID"], row["A1"]
        if key in result:
            raise ValueError("Duplicate scored allele")
        result[key] = [pp.numeric(row[pc], pc) for pc in pcs]
    return result


def _number_string(value):
    return format(float(value), ".17g")


def frame_spec(plan):
    """Reconstruct the actual common transform from existing prepared files.

    Target native IDs/orientation do not enter the frame identity.  Their
    equivalence to the reference transform is checked explicitly instead.
    """
    root, source_plan, _ = pp._completed_outputs(plan["basis_dir"])
    if source_plan["signature"] != plan["basis_signature"]:
        raise ValueError("Projection source basis signature changed")
    if plan["build"] != "GRCh37" or plan["pcs"] != [f"PC{i+1}" for i in range(len(plan["pcs"]))]:
        raise ValueError("Unexpected coordinate build/PC order")
    for path, expected in plan["prepared_sha256"].items():
        if pp.sha(path) != expected:
            raise ValueError("Prepared projection input changed: " + path)
    selected = [row for row in pp.rows(Path(plan["out"]) / "marker_audit.tsv.gz")
                if row["status"].startswith("matched_")]
    selected.sort(key=lambda r: (int(r["chrom"]), int(r["pos"]), r["ref"], r["alt"], r["id"]))
    marker_ids = {r["id"] for r in selected}
    coordinates = {(r["chrom"], r["pos"]) for r in selected}
    if (not selected or len(marker_ids) != len(selected) or len(coordinates) != len(selected)
            or len(selected) != plan["n_shared_markers"]):
        raise ValueError("Common marker identities/count disagree")
    reference_stages = [s for s in plan["stages"] if not s["target"]]
    if len(reference_stages) != 1:
        raise ValueError("Expected one common-panel reference score stage")
    reference = reference_stages[0]
    ref_weights = _weights(reference["prepared"]["weights"], plan["pcs"])
    ref_counts = _by_id(reference["prepared"]["frequencies"], {"ID", "CHROM", "REF", "ALT", "ALT_CTS", "OBS_CT"})
    source_weights = _weights(root / "ref.eigenvec.allele", plan["pcs"])
    source_counts = _by_id(root / "ref.acount", {"ID", "CHROM", "REF", "ALT", "ALT_CTS", "OBS_CT"})
    if ref_weights != {key: value for key, value in source_weights.items() if key[0] in marker_ids}:
        raise ValueError("Prepared loadings differ from the frozen source axes")
    weights_by_id = {}
    for (vid, allele), value in ref_weights.items():
        weights_by_id.setdefault(vid, {})[allele] = value
    score_args = SimpleNamespace(**plan["settings"], components=len(plan["pcs"]))
    for stage in plan["stages"]:
        if stage["command"] != pp.score_command(score_args, stage["files"], stage["prepared"], stage["output"]):
            raise ValueError("Executed score command does not implement the declared transform")
    ids = Path(reference["prepared"]["extract"]).read_text().splitlines()
    if set(ids) != marker_ids or len(ids) != len(selected) or set(ref_counts) != marker_ids:
        raise ValueError("Prepared reference panel differs from selected coordinates")
    weights, frequencies = [], []
    for marker in selected:
        vid, ref, alt = marker["id"], marker["ref"], marker["alt"]
        row = ref_counts[vid]
        if (row["CHROM"], row["REF"], row["ALT"]) != (marker["chrom"], ref, alt):
            raise ValueError("Prepared reference frequency identity differs")
        original_count = source_counts.get(vid)
        if original_count is None or any(row[key] != original_count[key] for key in ("CHROM", "REF", "ALT")):
            raise ValueError("Prepared reference contrast differs from frozen source")
        ac, an = pp.numeric(row["ALT_CTS"], "AC"), pp.numeric(row["OBS_CT"], "AN")
        if ac != pp.numeric(original_count["ALT_CTS"], "AC") or an != pp.numeric(original_count["OBS_CT"], "AN"):
            raise ValueError("Prepared reference counts differ from frozen source")
        if not 0 < ac < an:
            raise ValueError("Reference variance must be positive")
        alleles = sorted(weights_by_id.get(vid, {}))
        if not alleles or not set(alleles) <= {ref, alt}:
            raise ValueError("Missing/invalid common-panel allele weights")
        weights.extend([[vid, allele, *[_number_string(v) for v in ref_weights[vid, allele]]]
                        for allele in alleles])
        frequencies.append([vid, ref, alt, _number_string(ac), _number_string(an)])
    if {v for v, _ in ref_weights} != marker_ids:
        raise ValueError("Unexpected reference weight variant")
    visited = set()
    for stage in [s for s in plan["stages"] if s["target"]]:
        match = re.fullmatch(r"target_(\d+)", stage["label"])
        if not match:
            raise ValueError("Unrecognized current projection target stage")
        index = int(match[1]) - 1
        subset = [r for r in selected if int(r["target_index"]) == index]
        native_ids = [r["target_id"] for r in subset]
        actual_ids = Path(stage["prepared"]["extract"]).read_text().splitlines()
        if (set(actual_ids) != set(native_ids) or len(actual_ids) != len(native_ids)
                or stage["n_markers"] != len(subset)):
            raise ValueError("Target prepared marker subset differs")
        tw = _weights(stage["prepared"]["weights"], plan["pcs"])
        tc = _by_id(stage["prepared"]["frequencies"], {"ID", "CHROM", "REF", "ALT", "ALT_CTS", "OBS_CT"})
        if set(tc) != set(native_ids):
            raise ValueError("Unexpected target frequency variants")
        expected_tw = {}
        for marker in subset:
            vid, native = marker["id"], marker["target_id"]
            if vid in visited:
                raise ValueError("Common marker scored in multiple target stages")
            visited.add(vid)
            if {marker["target_ref"], marker["target_alt"]} != {marker["ref"], marker["alt"]}:
                raise ValueError("Target allele contrast differs")
            reference_count = ref_counts[vid]
            ac, an = float(reference_count["ALT_CTS"]), float(reference_count["OBS_CT"])
            expected_ac = ac if marker["target_alt"] == marker["alt"] else an - ac
            row = tc[native]
            if ((row["CHROM"], row["REF"], row["ALT"]) !=
                    (marker["chrom"], marker["target_ref"], marker["target_alt"])
                    or float(row["ALT_CTS"]) != expected_ac or float(row["OBS_CT"]) != an):
                raise ValueError("Target does not use the same reference frequencies")
            expected_tw.update({(native, allele): value for allele, value in weights_by_id[vid].items()})
        if tw != expected_tw:
            raise ValueError("Target does not use the same allele loadings")
    if visited != marker_ids:
        raise ValueError("Target stages do not partition the common panel")
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(["chrom", "pos", "id", "ref", "alt"])
    writer.writerows([[r[k] for k in ("chrom", "pos", "id", "ref", "alt")] for r in selected])
    panel_text = buffer.getvalue()
    source_artifacts = {name: artifact(root / name) for name in
                        ("pca_plan.json", "PCA_COMPLETE", "completed.sha256",
                         "ref.eigenvec.allele", "ref.acount", "ref.eigenval")}
    transform = {
        "schema": FRAME_SCHEMA, "build": plan["build"], "source_basis_id": plan["basis_signature"],
        "source_axes_sha256": source_artifacts["ref.eigenvec.allele"]["sha256"],
        "source_frequencies_sha256": source_artifacts["ref.acount"]["sha256"],
        "source_eigenvalues_sha256": source_artifacts["ref.eigenval"]["sha256"],
        "pcs": plan["pcs"], "n_markers": len(selected),
        "marker_sha256": hashlib.sha256(panel_text.encode()).hexdigest(),
        "used_loadings_sha256": digest_json(weights), "used_reference_frequencies_sha256": digest_json(frequencies),
        "scoring": {"centering": "reference_frequency", "variance_standardize": True,
                    "missing": "reference_mean_imputation", "coordinate_kind": "projected_average",
                    "denominator": 2 * len(selected), "chromosomes": "autosomes_1_22",
                    "aggregation": "sum_scores_and_denominators_before_division"},
    }
    return {"schema": FRAME_SCHEMA, "coordinate_frame_id": "sha256:" + digest_json(transform),
            "transform": transform, "source_basis_dir": str(root), "source_artifacts": source_artifacts}, panel_text


def _score_identity(path, pcs, expected_samples, n_markers):
    required = {"FID", "IID", "DENOM", *[pc + "_AVG" for pc in pcs]}
    samples = []
    for row in pp.rows(path, required):
        samples.append((row["FID"], row["IID"]))
        if pp.numeric(row["DENOM"], "DENOM") != 2 * n_markers:
            raise ValueError("Coordinate denominator differs from common panel")
        for pc in pcs:
            pp.numeric(row[pc + "_AVG"], pc)
    if samples != expected_samples or len({sid for _, sid in samples}) != len(samples):
        raise ValueError("Coordinate samples differ from original PSAM")
    return {"file": artifact(path), "sample_ids_sha256": digest_json(samples), "n_samples": len(samples)}


def _check_metadata(binding, sample_ids):
    if not binding.get("file"):
        return
    path = verify_artifact(binding["file"])
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "sample_id" not in (reader.fieldnames or []):
            raise ValueError("Coordinate metadata requires sample_id")
        ids = [r["sample_id"] for r in reader]
    if len(ids) != len(set(ids)) or not set(sample_ids) <= set(ids):
        raise ValueError("Coordinate metadata has duplicates or lacks scored samples")


def publish_projection_receipts(plan, stage_receipts):
    """Call after score aggregation, before publishing PROJECTION_COMPLETE.

    Return every newly generated path; add them to completion.outputs.  A
    consumer rejects these receipts until that final completion exists.
    """
    out = Path(plan["out"])
    frame, panel_text = frame_spec(plan)
    for stage in plan["stages"]:
        saved = stage_receipts[stage["label"]]
        if saved["status"] != "complete" or saved["identity"]["signature"] != plan["signature"]:
            raise ValueError("Scoring stage not verified for this projection")
        pp.check_output_state(saved["outputs"])
        expected = Path(stage["prepared"]["extract"]).read_text().splitlines()
        used = Path(stage["output"] + ".sscore.vars").read_text().splitlines()
        if set(used) != set(expected) or len(used) != len(expected):
            raise ValueError("Actually scored stage markers differ")
    panel_path, frame_path = out / "common_panel.tsv", out / "coordinate_frame.json"
    pp.atomic_text(panel_path, panel_text)
    frame["common_panel"] = artifact(panel_path)
    frame["projection_plan"] = artifact(out / "projection_plan.json")
    pp.save_json(frame_path, frame)
    scores, samples = {}, {}
    for role in ("reference", "target"):
        stage = next(s for s in plan["stages"] if s["target"] == (role == "target"))
        samples[role] = pp.sample_ids(stage["files"]["psam"])
        scores[role] = _score_identity(out / f"{role}.sscore", plan["pcs"], samples[role], plan["n_shared_markers"])
    result = [panel_path, frame_path]
    bindings = plan["coordinate_metadata"]
    for role in ("reference", "target"):
        _check_metadata(bindings[role], [sid for _, sid in samples[role]])
        proof = {
            "schema": RECEIPT_SCHEMA, "status": "verified", "coordinate_kind": "projected_average",
            "role": role, "coordinates": scores[role]["file"], "coordinates_sha256": scores[role]["file"]["sha256"],
            "paired_scores": scores, "source_basis_id": plan["basis_signature"],
            "coordinate_frame_id": frame["coordinate_frame_id"], "basis_id": frame["coordinate_frame_id"],
            "marker_sha256": frame["transform"]["marker_sha256"], "frame": artifact(frame_path),
            "metadata": bindings[role], "projection_signature": plan["signature"],
            "completion_path": str(out / "PROJECTION_COMPLETE.json"),
        }
        path = out / f"{role}.coordinate_receipt.json"
        pp.save_json(path, proof)
        result.append(path)
    return result


def verify_projection_receipt(coordinates, receipt_path, pcs, basis_dir=None):
    """Validate a producer receipt, its complete run, both scores and transform."""
    proof = json.loads(Path(receipt_path).read_text())
    if (proof.get("schema") != RECEIPT_SCHEMA or proof.get("status") != "verified"
            or proof.get("coordinate_kind") != "projected_average" or proof.get("role") not in {"reference", "target"}):
        raise ValueError("Invalid projection coordinate receipt")
    if set(proof.get("paired_scores", {})) != {"reference", "target"}:
        raise ValueError("Coordinate receipt must bind both cohort score files")
    if Path(coordinates).resolve() != Path(proof["coordinates"]["path"]).resolve():
        raise ValueError("Coordinate receipt belongs to a different score file/role")
    completion = json.loads(Path(proof["completion_path"]).read_text())
    if completion["signature"] != proof["projection_signature"]:
        raise ValueError("Coordinate completion signature mismatch")
    pp.check_output_state(completion["outputs"])
    if str(Path(receipt_path).resolve()) not in completion["outputs"]:
        raise ValueError("Coordinate receipt is not bound by the completed projection")
    frame_path = verify_artifact(proof["frame"])
    frame = json.loads(frame_path.read_text())
    plan_path = verify_artifact(frame["projection_plan"])
    plan = json.loads(plan_path.read_text())
    if plan["signature"] != proof["projection_signature"]:
        raise ValueError("Coordinate projection plan mismatch")
    rebuilt, panel_text = frame_spec(plan)
    if (rebuilt["transform"] != frame["transform"] or rebuilt["coordinate_frame_id"] != frame["coordinate_frame_id"]
            or rebuilt["source_artifacts"] != frame["source_artifacts"]
            or rebuilt["source_basis_dir"] != frame["source_basis_dir"]):
        raise ValueError("Recorded coordinate frame differs from actual shared transform")
    for saved in frame["source_artifacts"].values():
        verify_artifact(saved)
    panel_path = verify_artifact(frame["common_panel"])
    if panel_path.read_text() != panel_text:
        raise ValueError("Recorded common panel differs from actual scored markers")
    if (proof["coordinate_frame_id"] != frame["coordinate_frame_id"]
            or proof["basis_id"] != frame["coordinate_frame_id"]
            or proof["source_basis_id"] != frame["transform"]["source_basis_id"]
            or proof["marker_sha256"] != frame["transform"]["marker_sha256"]
            or frame["transform"]["pcs"][:len(pcs)] != [p.removesuffix("_AVG") for p in pcs]):
        raise ValueError("Coordinate frame aliases or PC order disagree")
    if basis_dir and Path(basis_dir).resolve() != Path(frame["source_basis_dir"]).resolve():
        raise ValueError("Explicit basis differs from coordinate receipt source")
    for role, score in proof["paired_scores"].items():
        path = verify_artifact(score["file"])
        stage = next(s for s in plan["stages"] if s["target"] == (role == "target"))
        samples = pp.sample_ids(stage["files"]["psam"])
        if _score_identity(path, plan["pcs"], samples, plan["n_shared_markers"]) != score:
            raise ValueError("Coordinate sample identity differs")
    if proof["coordinates"] != proof["paired_scores"][proof["role"]]["file"]:
        raise ValueError("Role-specific score identity differs")
    if proof["coordinates_sha256"] != proof["coordinates"]["sha256"]:
        raise ValueError("Role-specific score SHA alias differs")
    if proof["metadata"] != plan["coordinate_metadata"][proof["role"]]:
        raise ValueError("Coordinate metadata differs from planned binding")
    score_ids = [row["IID"] for row in pp.rows(coordinates, {"IID"})]
    _check_metadata(proof["metadata"], score_ids)
    audit = {"status": "verified", "basis_id": proof["coordinate_frame_id"],
             "coordinate_frame_id": proof["coordinate_frame_id"], "source_basis_id": proof["source_basis_id"],
             "source_basis_dir": frame["source_basis_dir"], "coordinates": identity(coordinates, full=True),
             "supplied_receipt": identity(receipt_path, full=True), "marker_sha256": proof["marker_sha256"],
             "metadata_recorded": proof["metadata"].get("recorded"),
             "metadata_binding_status": proof["metadata"]["status"], "projection_role": proof["role"]}
    if proof["role"] == "reference":
        # Original reference sample training membership remains meaningful.
        audit["basis_dir"] = frame["source_basis_dir"]
    return audit


def require_same_coordinate_frame(*contracts):
    """Use before comparing/concatenating independently verified coordinates."""
    if not contracts or any(c.get("status") != "verified" or not c.get("coordinate_frame_id") for c in contracts):
        raise ValueError("Comparison requires verified coordinate frame identities")
    if len({c["coordinate_frame_id"] for c in contracts}) != 1:
        raise ValueError("Different marker/weight/frequency/scale frames; re-score on one common panel")
    return contracts[0]["coordinate_frame_id"]
