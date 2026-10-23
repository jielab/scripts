"""Allele-specific, descriptive trajectories from targeted ancient PGEN rows.

This module does not estimate selection coefficients, reconstruct unrepresented
alleles, infer REF ancestry, or claim to reproduce the published selection GLMM.
The aggregation layer needs only the standard library.  PGEN access lazily
imports numpy and the official pgenlib Python API.
"""
from __future__ import annotations

from bisect import bisect_right
from collections import Counter
from contextlib import ExitStack
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import tempfile


SCHEMA_VERSION = "allele3-trajectory-1"
DEFAULT_BINS = (0, 2000, 4000, 6000, 8000, 10000)
STATUSES = {"assayed_pass", "not_in_high_quality", "not_represented_in_panel",
            "represented_multiallelic_unsupported"}
KINDS = {"hardcall", "dosage"}
META_REQUIRED = {"sample_id", "cohort", "age_bp", "include", "stratum"}
SELECTED_REQUIRED = {"variant_id", "pvar_index", "chrom", "pos", "ref", "alt", "status"}
BASE_FIELDS = ["variant_id", "build", "chrom", "pos", "ref", "alt", "pvar_index"]
TIME_FIELDS = ["cohort", "time_bin", "start_bp", "end_bp", "right_inclusive"]
BIN_FIELDS = BASE_FIELDS + TIME_FIELDS + [
    "group_type", "stratum", "count_kind", "dosage_evidence", "n_eligible",
    "n_nonmissing", "call_rate", "ac", "expected_ac", "an", "af", "bin_status",
]
STD_FIELDS = BASE_FIELDS + TIME_FIELDS + [
    "count_kind", "dosage_evidence", "n_reference_strata", "n_supported_strata",
    "supported_weight", "missing_strata", "af_standardized", "standardization_status",
]
AVAIL_FIELDS = BASE_FIELDS + [
    "status", "count_kind", "dosage_evidence", "n_eligible", "n_nonmissing",
    "trajectory_available", "availability_note",
]
QPADM_COMPONENTS = {"3way": ("whg", "eef", "steppe"),
                    "4way": ("anf", "ehg", "whg", "icr")}
OUTPUT_NAMES = ("trajectory_bins.tsv", "trajectory_standardized.tsv",
                "sample_bin_composition.tsv", "sample_inclusion.tsv",
                "allele_availability.tsv", "trajectory_manifest.json")


def _open_text(path):
    path = Path(path)
    opener = gzip.open if path.suffix.lower() in {".gz", ".bgz", ".bgzip"} else open
    return opener(path, "rt", encoding="utf-8-sig", newline="")


def _rows(path, required):
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)) or not required.issubset(fields):
            raise ValueError(f"{path}: requires unique TSV columns {sorted(required)}")
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{line}: wrong number of TSV fields")
            yield {key: value.strip() for key, value in row.items()}


def _missing(value):
    return value is None or (isinstance(value, str) and value.strip().lower() in
                             {"", ".", "na", "nan", "none"})


def _number(value, label):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label}: expected a finite number, got {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label}: expected a finite number, got {value!r}")
    return result


def _integer(value, label, minimum=0):
    result = _number(value, label)
    if result != int(result) or result < minimum:
        raise ValueError(f"{label}: expected an integer >= {minimum}, got {value!r}")
    return int(result)


def _sample_ids(values):
    ids = list(values)
    if not ids or any(not isinstance(x, str) or _missing(x) or x != x.strip() for x in ids):
        raise ValueError("sample_ids must be a nonempty sequence of exact, nonmissing strings")
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate sample_id; IID alone must uniquely identify every PSAM row")
    return ids


def read_psam_ids(psam_path):
    """Read unique IID values in original PSAM order; never guess IID from FID."""
    ids, fields = [], None
    with _open_text(psam_path) as handle:
        for line_no, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("##"):
                continue
            parts = line.split()
            if fields is None:
                fields = [x.lstrip("#") for x in parts]
                if "IID" not in fields or len(fields) != len(set(fields)):
                    raise ValueError(f"{psam_path}: PSAM must have a unique #IID or #FID ... IID header")
                iid_col = fields.index("IID")
                continue
            if line.startswith("#") or len(parts) != len(fields):
                raise ValueError(f"{psam_path}:{line_no}: malformed PSAM row")
            ids.append(parts[iid_col])
    return _sample_ids(ids)


def _variant(row):
    result = dict(row)
    missing = SELECTED_REQUIRED - result.keys()
    if missing:
        raise ValueError(f"Variant is missing fields: {sorted(missing)}")
    vid = result["variant_id"]
    if not isinstance(vid, str) or _missing(vid):
        raise ValueError("Each target needs a nonmissing variant_id")
    result["build"] = result.get("build", "GRCh37")
    if result["build"] != "GRCh37":
        raise ValueError(f"{vid}: this PGEN trajectory adapter requires explicitly harmonized GRCh37 targets")
    chrom = str(result["chrom"]).removeprefix("chr")
    if chrom not in {str(x) for x in range(1, 23)}:
        raise ValueError(f"{vid}: only autosomes 1-22 are supported, got {chrom!r}")
    result["chrom"] = chrom
    result["pos"] = _integer(result["pos"], f"{vid} pos", 1)
    for key in ("ref", "alt"):
        result[key] = str(result[key]).upper()
        if result[key] not in {"A", "C", "G", "T"}:
            raise ValueError(f"{vid}: only single A/C/G/T alleles are supported")
    if result["ref"] == result["alt"]:
        raise ValueError(f"{vid}: REF and ALT must differ")
    status = result["status"]
    if status not in STATUSES:
        raise ValueError(f"{vid}: unsupported availability status {status!r}")
    index = result["pvar_index"]
    result["pvar_index"] = None if _missing(index) else _integer(index, f"{vid} pvar_index")
    if status in {"assayed_pass", "not_in_high_quality"} and result["pvar_index"] is None:
        raise ValueError(f"{vid}: a represented panel row needs its original pvar_index")
    if status in {"not_represented_in_panel", "represented_multiallelic_unsupported"} and result["pvar_index"] is not None:
        raise ValueError(f"{vid}: an unrepresented/unsupported allele must not be assigned a readable PVAR row")
    return result


def _load_pgen_backend():
    try:
        import numpy as np
        import pgenlib
    except ImportError as exc:
        raise RuntimeError("Targeted PGEN reading needs numpy and pgenlib. Install into the active "
                           "environment with: python -m pip install numpy pgenlib") from exc
    return np, pgenlib


def _file_info(path, digest=False):
    path = Path(path).expanduser().resolve()
    info = {"path": str(path), "size": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns}
    if digest:
        hasher = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1048576), b""):
                hasher.update(chunk)
        info["sha256"] = hasher.hexdigest()
    return info


def iter_pgen_variants(pgen_path, psam_path, selected_variants_path, count_kind="hardcall"):
    """Return (sample_ids, iterator) for selected original zero-based PVAR rows.

    Selected TSV: variant_id,pvar_index,chrom,pos,ref,alt,status.  The caller
    must have verified those REF/ALT labels against the original, biallelic
    PVAR; PGEN does not contain the nucleotide labels.  No filtered/reindexed
    PVAR may be substituted.  Read allele_idx=1 explicitly.  Every yielded
    dosages list is independently owned and is in the complete PSAM order.

    ``dosage`` explicitly selects read_dosages.  That API can use hardcalls
    when explicit dosages are absent, so dosage_evidence is always recorded;
    integer-only outputs do not establish that posterior dosages were stored.
    """
    if count_kind not in KINDS:
        raise ValueError("count_kind must be hardcall or dosage")
    sample_ids = read_psam_ids(psam_path)
    selected, seen_ids, seen_indexes = [], set(), set()
    for raw in _rows(selected_variants_path, SELECTED_REQUIRED):
        row = _variant(raw)
        if row["variant_id"] in seen_ids:
            raise ValueError(f"Duplicate selected variant_id: {row['variant_id']}")
        seen_ids.add(row["variant_id"])
        index = row["pvar_index"]
        if index is not None:
            if index in seen_indexes:
                raise ValueError(f"Original PVAR row {index} was assigned to more than one target; "
                                 "this adapter accepts one biallelic ALT per original row")
            seen_indexes.add(index)
        selected.append(row)

    def generate():
        reader = None
        source = {"reader": "pgenlib", "psam": _file_info(psam_path, True),
                  "selected_variants": _file_info(selected_variants_path, True),
                  "requested_count_kind": count_kind, "allele_idx": 1,
                  "index_basis": "zero_based_original_unfiltered_PVAR_data_rows"}
        try:
            for selected_row in selected:
                row = dict(selected_row, count_kind=count_kind, dosages=None,
                           dosage_evidence="not_read", source=source)
                if row["status"] != "assayed_pass":
                    yield row
                    continue
                if reader is None:
                    np, pgenlib = _load_pgen_backend()
                    source["pgen"] = _file_info(pgen_path)
                    try:
                        source["pgenlib_version"] = importlib.metadata.version("pgenlib")
                    except importlib.metadata.PackageNotFoundError:
                        source["pgenlib_version"] = getattr(pgenlib, "__version__", "unknown")
                    reader = pgenlib.PgenReader(os.fsencode(Path(pgen_path).expanduser()),
                                               raw_sample_ct=len(sample_ids))
                    if reader.get_raw_sample_ct() != len(sample_ids):
                        raise ValueError("PGEN sample count disagrees with the full PSAM")
                    if count_kind == "dosage" and not callable(getattr(reader, "read_dosages", None)):
                        raise RuntimeError("This pgenlib has no read_dosages API; upgrade pgenlib or "
                                           "explicitly choose hardcall. No automatic fallback was used.")
                    buffer = np.empty(len(sample_ids), dtype=np.int8 if count_kind == "hardcall" else np.float64)
                    variant_ct = reader.get_variant_ct()
                index = row["pvar_index"]
                if index >= variant_ct:
                    raise ValueError(f"{row['variant_id']}: pvar_index {index} is outside PGEN's "
                                     f"{variant_ct} original rows")
                try:
                    if count_kind == "hardcall":
                        reader.read(index, buffer, allele_idx=1)
                    else:
                        reader.read_dosages(index, buffer, allele_idx=1)
                except Exception as exc:
                    raise RuntimeError(f"PGEN read failed for {row['variant_id']} at original row {index}: {exc}") from exc
                row["dosages"] = [None if x == -9 else (int(x) if count_kind == "hardcall" else float(x))
                                   for x in buffer]
                if count_kind == "dosage":
                    present = [x for x in row["dosages"] if x is not None]
                    row["dosage_evidence"] = ("fractional_values_present" if any(x != int(x) for x in present)
                                              else "integer_only_storage_unverified" if present else "all_missing")
                else:
                    row["dosage_evidence"] = "hardcalls_only"
                yield row
        finally:
            if reader is not None:
                reader.close()

    return sample_ids, generate()


def _metadata(path, sample_ids, bins):
    by_id = {}
    for row in _rows(path, META_REQUIRED):
        sid = row["sample_id"]
        if _missing(sid) or sid in by_id:
            raise ValueError(f"Duplicate or missing metadata sample_id: {sid!r}")
        if row["include"] not in {"0", "1"}:
            raise ValueError(f"{sid}: include must be exactly 0 or 1")
        if row["cohort"] not in {"ancient", "modern"}:
            raise ValueError(f"{sid}: cohort must be ancient or modern")
        row["included"] = row["include"] == "1"
        row["age"] = None
        if not _missing(row["age_bp"]):
            row["age"] = _number(row["age_bp"], f"{sid} age_bp")
        if row["included"] and _missing(row["stratum"]):
            raise ValueError(f"{sid}: included samples need a predeclared stratum")
        if row["included"] and row["cohort"] == "ancient" and (row["age"] is None or row["age"] < 0):
            raise ValueError(f"{sid}: included ancient sample needs a finite nonnegative age_bp")
        row["bin_index"] = None
        if not row["included"]:
            row["inclusion_status"] = "metadata_excluded"
        elif row["cohort"] == "modern":
            row["bin_index"] = -1
            row["inclusion_status"] = "included"
        elif not bins[0] <= row["age"] <= bins[-1]:
            row["inclusion_status"] = "outside_time_range"
        else:
            row["bin_index"] = min(bisect_right(bins, row["age"]) - 1, len(bins) - 2)
            row["inclusion_status"] = "included"
        by_id[sid] = row
    expected, found = set(sample_ids), set(by_id)
    if expected != found:
        raise ValueError("Metadata IDs must exactly match PSAM/sample_ids: "
                         f"missing={sorted(expected-found)[:5]}, extra={sorted(found-expected)[:5]}")
    return [by_id[sid] for sid in sample_ids]


def _weights(path):
    if path is None:
        return None
    weights = {}
    for row in _rows(path, {"stratum", "weight"}):
        key = row["stratum"]
        if _missing(key) or key in weights:
            raise ValueError(f"Duplicate or missing weight stratum: {key!r}")
        weight = _number(row["weight"], f"{key} weight")
        if weight < 0:
            raise ValueError("Reference weights must be nonnegative")
        weights[key] = weight
    total = math.fsum(weights.values())
    if not math.isclose(total, 1.0, rel_tol=0, abs_tol=1e-8):
        raise ValueError(f"Reference weights must sum to 1, got {total}")
    return {key: weight/total for key, weight in weights.items() if weight > 0}


def _time_fields(cohort, index, bins):
    if cohort == "modern":
        return dict(cohort=cohort, time_bin="modern", start_bp=".", end_bp=".", right_inclusive=".")
    start, end = bins[index:index+2]
    last = index == len(bins) - 2
    return dict(cohort=cohort, time_bin=f"[{start:g},{end:g}{']' if last else ')'}",
                start_bp=start, end_bp=end, right_inclusive=int(last))


def _writer(stack, stage, name, fields):
    handle = stack.enter_context((stage / name).open("w", newline="", encoding="utf-8"))
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    return writer


def _write(writer, row):
    writer.writerow({key: "." if value is None else value for key, value in row.items()})


def _composition(metadata, members):
    rows = [metadata[i] for i in members]
    ages = [r["age"] for r in rows if r["age"] is not None and r["cohort"] == "ancient"]
    result = {"n_eligible": len(rows), "mean_age_bp": math.fsum(ages)/len(ages) if ages else None}
    for source in ("region", "batch", "data_source", "data_types"):
        result[source + "_counts_json"] = json.dumps(dict(sorted(Counter(r.get(source, ".") or "."
                                                               for r in rows).items())), ensure_ascii=False)
    for model, components in QPADM_COMPONENTS.items():
        valid = [r for r in rows if r.get("qpadm_valid") == "1" and r.get("qpadm_model") == model]
        values = []
        for row in valid:
            vector = [_number(row.get("qpadm_" + c), f"{row['sample_id']} qpadm_{c}") for c in components]
            if any(not 0 <= x <= 1 for x in vector) or not math.isclose(sum(vector), 1, abs_tol=0.02):
                raise ValueError(f"{row['sample_id']}: qpadm_valid=1 contradicts invalid {model} proportions")
            values.append(vector)
        result[f"n_qpadm_{model}_valid"] = len(values)
        for j, component in enumerate(components):
            result[f"mean_qpadm_{model}_{component}"] = math.fsum(v[j] for v in values)/len(values) if values else None
    return result


def _call(value, kind, label):
    if value is None:
        return None
    # The generic interface uses None for missingness; pgenlib's -9 is
    # translated by the reader.  Reject NaN and all other sentinels here.
    value = _number(value, label)
    if not 0 <= value <= 2 or (kind == "hardcall" and value != int(value)):
        raise ValueError(f"{label}: {kind} needs {'0, 1, or 2' if kind == 'hardcall' else 'a value in [0,2]'} or None")
    return value


def aggregate_trajectories(variants, sample_ids, metadata_path, out_dir,
                           bins=DEFAULT_BINS, min_samples=10, weights_path=None):
    """Write streaming, descriptive ALT trajectories and return a manifest.

    ``variants`` yields dictionaries from iter_pgen_variants (or the same
    protocol).  Low-sample raw bins retain counts and descriptive AF, with
    bin_status=below_min_samples.  Standardized AF requires every positively
    weighted stratum to have >=min_samples nonmissing calls.  Missing or
    unrepresented alleles never become a synthetic all-zero trajectory.
    All output files are committed only after successful validation.
    """
    sample_ids = _sample_ids(sample_ids)
    edges = tuple(_number(x, "bin edge") for x in bins)
    if len(edges) < 2 or edges[0] < 0 or any(a >= b for a, b in zip(edges, edges[1:])):
        raise ValueError("bins must contain at least two strictly increasing nonnegative boundaries")
    min_samples = _integer(min_samples, "min_samples", 1)
    metadata = _metadata(metadata_path, sample_ids, edges)
    weights = _weights(weights_path)
    time_keys = [("ancient", index) for index in range(len(edges)-1)] + [("modern", -1)]
    strata = {cohort: sorted({r["stratum"] for r in metadata
                             if r["cohort"] == cohort and r["inclusion_status"] == "included"})
              for cohort in ("ancient", "modern")}
    groups = [(cohort, index, None) for cohort, index in time_keys]
    groups += [(cohort, index, stratum) for cohort, index in time_keys for stratum in strata[cohort]]
    members = {key: [] for key in groups}
    sample_groups = []
    for index, row in enumerate(metadata):
        keys = []
        if row["inclusion_status"] == "included":
            keys = [(row["cohort"], row["bin_index"], None),
                    (row["cohort"], row["bin_index"], row["stratum"])]
            for key in keys:
                members[key].append(index)
        sample_groups.append(keys)
    active = sum(bool(keys) for keys in sample_groups)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    existing = [name for name in OUTPUT_NAMES if (out / name).exists()]
    if existing:
        raise ValueError(f"Trajectory output already exists ({existing[0]}); use a fresh output directory")
    stage = Path(tempfile.mkdtemp(prefix=".trajectory-", dir=out))
    counts = Counter({key: 0 for key in ("n_targets", "n_bin_rows", "n_standardized_rows", *sorted(STATUSES))})
    kinds, sources = Counter(), {}
    seen_ids, seen_keys = set(), set()
    iterator = iter(variants)
    try:
        with ExitStack() as stack:
            bins_writer = _writer(stack, stage, "trajectory_bins.tsv", BIN_FIELDS)
            std_writer = _writer(stack, stage, "trajectory_standardized.tsv", STD_FIELDS)
            avail_writer = _writer(stack, stage, "allele_availability.tsv", AVAIL_FIELDS)
            composition_fields = TIME_FIELDS + ["group_type", "stratum", "n_eligible", "mean_age_bp"]
            composition_fields += [x + "_counts_json" for x in ("region", "batch", "data_source", "data_types")]
            for model, components in QPADM_COMPONENTS.items():
                composition_fields += [f"n_qpadm_{model}_valid"] + [f"mean_qpadm_{model}_{c}" for c in components]
            comp_writer = _writer(stack, stage, "sample_bin_composition.tsv", composition_fields)
            for cohort, index, stratum in groups:
                _write(comp_writer, dict(_time_fields(cohort, index, edges),
                                         group_type="raw_all" if stratum is None else "by_stratum",
                                         stratum=stratum, **_composition(metadata, members[(cohort, index, stratum)])))
            inclusion = _writer(stack, stage, "sample_inclusion.tsv",
                                ["sample_id", "cohort", "age_bp", "include", "stratum", "inclusion_status"])
            for row in metadata:
                _write(inclusion, {key: row[key] for key in inclusion.fieldnames})
            for raw in iterator:
                row = _variant(raw)
                vid = row["variant_id"]
                key = tuple(row[x] for x in ("build", "chrom", "pos", "ref", "alt"))
                if vid in seen_ids or key in seen_keys:
                    raise ValueError(f"Duplicate trajectory target ID or exact allele: {vid}")
                seen_ids.add(vid)
                seen_keys.add(key)
                kind = row.get("count_kind")
                if kind not in KINDS:
                    raise ValueError(f"{vid}: count_kind must explicitly be hardcall or dosage")
                evidence = row.get("dosage_evidence", "hardcalls_only" if kind == "hardcall" else "external_dosage_provenance_unverified")
                counts["n_targets"] += 1
                counts[row["status"]] += 1
                kinds[kind] += 1
                if row.get("source"):
                    encoded = json.dumps(row["source"], sort_keys=True)
                    sources[encoded] = dict(row["source"])
                base = {x: row[x] for x in BASE_FIELDS}
                availability = dict(base, status=row["status"], count_kind=kind,
                                    dosage_evidence=evidence, n_eligible=active)
                if row["status"] != "assayed_pass":
                    if row.get("dosages") is not None:
                        raise ValueError(f"{vid}: unavailable/failed-QC allele must have dosages=None")
                    note = {"not_represented_in_panel": "unrepresented_allele_not_measured_not_zero",
                            "not_in_high_quality": "excluded_by_variant_quality_control",
                            "represented_multiallelic_unsupported": "represented_but_multiallelic_reader_not_supported"}[row["status"]]
                    _write(avail_writer, dict(availability, n_nonmissing=None, trajectory_available=0,
                                              availability_note=note))
                    continue
                calls = row.get("dosages")
                if calls is None or len(calls) != len(sample_ids):
                    raise ValueError(f"{vid}: dosages length must equal the complete sample_ids length")
                totals = {group: [0, 0.0] for group in groups}
                nonmissing = 0
                for sample_index, value in enumerate(calls):
                    value = _call(value, kind, f"{vid}/{sample_ids[sample_index]}")
                    if value is None or not sample_groups[sample_index]:
                        continue
                    nonmissing += 1
                    for group in sample_groups[sample_index]:
                        totals[group][0] += 1
                        totals[group][1] += value
                _write(avail_writer, dict(availability, n_nonmissing=nonmissing,
                                          trajectory_available=int(nonmissing > 0),
                                          availability_note="descriptive_biallelic_panel_estimate" if nonmissing else "no_nonmissing_included_calls"))
                for cohort, index, stratum in groups:
                    group = (cohort, index, stratum)
                    n, allele_sum = totals[group]
                    eligible = len(members[group])
                    bin_status = ("no_eligible_samples" if not eligible else "all_calls_missing" if not n
                                  else "below_min_samples" if n < min_samples else "ok")
                    record = dict(base, **_time_fields(cohort, index, edges),
                                  group_type="raw_all" if stratum is None else "by_stratum", stratum=stratum,
                                  count_kind=kind, dosage_evidence=evidence, n_eligible=eligible, n_nonmissing=n,
                                  call_rate=n/eligible if eligible else None,
                                  ac=int(allele_sum) if kind == "hardcall" else None,
                                  expected_ac=allele_sum if kind == "dosage" else None,
                                  an=2*n, af=allele_sum/(2*n) if n else None, bin_status=bin_status)
                    _write(bins_writer, record)
                    counts["n_bin_rows"] += 1
                if weights is not None:
                    for cohort, index in time_keys:
                        supported, missing, contributions = [], [], []
                        for stratum, weight in weights.items():
                            n, allele_sum = totals.get((cohort, index, stratum), (0, 0))
                            if n < min_samples:
                                missing.append(stratum)
                            else:
                                supported.append(weight)
                                contributions.append(weight*allele_sum/(2*n))
                        _write(std_writer, dict(base, **_time_fields(cohort, index, edges), count_kind=kind,
                                                dosage_evidence=evidence, n_reference_strata=len(weights),
                                                n_supported_strata=len(supported), supported_weight=math.fsum(supported),
                                                missing_strata=json.dumps(sorted(missing), ensure_ascii=False),
                                                af_standardized=math.fsum(contributions) if not missing else None,
                                                standardization_status="complete_support" if not missing else "insufficient_stratum_support"))
                        counts["n_standardized_rows"] += 1
        manifest = {
            "schema_version": SCHEMA_VERSION, "status": "complete",
            "created_utc": datetime.now(timezone.utc).isoformat(), "build": "GRCh37",
            "scope": "autosomal_diploid_descriptive_allele_specific_trajectories",
            "metadata": _file_info(metadata_path, True),
            "sample_order_sha256": hashlib.sha256("\n".join(sample_ids).encode()).hexdigest(),
            "n_samples": len(sample_ids), "n_included_in_time_scope": active,
            "sample_inclusion_counts": dict(Counter(r["inclusion_status"] for r in metadata)),
            "bins_bp": list(edges), "bin_convention": "left_closed_right_open_except_final_right_closed",
            "modern_endpoint": "separate_cohort_never_pooled_into_youngest_ancient_bin",
            "min_samples": min_samples, "raw_low_sample_af": "retained_with_below_min_samples_status",
            "weights_input": _file_info(weights_path, True) if weights_path else None,
            "reference_weights": weights,
            "standardization": "fixed_weights_complete_support_only" if weights else "not_requested",
            "counts": dict(counts), "count_kind_counts": dict(kinds), "sources": list(sources.values()),
            "outputs": {name: str((out / name).resolve()) for name in OUTPUT_NAMES},
            "limitations": [
                "No selection coefficients or inferential P values are estimated.",
                "Unrepresented alleles are unmeasured, never zero-frequency trajectories.",
                "AN is twice nonmissing diploid calls within the panel's allele model, not independent verification of all bases.",
                "Dosage reading can include hardcall fallback; integer values alone do not verify stored posterior dosages.",
                "PGEN does not preserve genotype probability triplets or per-call read/quality information.",
                "Fixed stratum weights control only the specified observed composition; they do not reproduce the published GLMM.",
                "Point sample dates are used; date and imputation uncertainty are not propagated.",
                "A published biallelic selection coefficient must not be transferred to another ALT at the same coordinate.",
                "Directional selection does not establish balancing selection or pleiotropic maintenance.",
            ],
        }
        (stage / "trajectory_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
        for name in OUTPUT_NAMES:
            os.replace(stage / name, out / name)
        return manifest
    finally:
        close = getattr(iterator, "close", None)
        if callable(close):
            close()
        shutil.rmtree(stage, ignore_errors=True)
