#!/usr/bin/env python3
"""Build an offline UKB disease-date audit; preserve all patient/event rows.

Inputs are emitted by date_check.R. This renderer performs no disease recoding,
date repair, or participant filtering. It uses only the Python standard library.
"""
from __future__ import annotations

import argparse
import base64
from collections import Counter, defaultdict
import csv
from datetime import date, datetime, timezone
from functools import lru_cache
import gzip
import hashlib
import html
import json
from pathlib import Path
import statistics


SOURCES = ("icd10", "icd9", "opcs4", "srd", "gp2", "gp3", "gp", "death", "me", "local_min",
           "reconstructed_ref", "ref", "saved_ref", "official_icd10", "official_icd9", "official_srd",
           "official_gp2", "official_gp3", "official_death")
SPECIAL_DATES = {"1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"}
MISSING = {"", "NA", "NaN", "NULL", "."}
DEFAULT_LIMITATIONS = [
    "This is an audit of the supplied files and disease definitions, not an exact reproduction of every UKB release.",
    "Hospital summary fields hold the earliest date per code; their aod is all retained code-dates, not all hospital episodes.",
    "Repeated GP record dates may be available only when record-level GP inputs were supplied. No missing history is invented.",
    "A blank date means no usable date in the supplied evidence. It does not establish absence of disease or GP coverage.",
    "Coding 819 sentinel dates are kept as raw evidence and excluded from numerical date comparisons.",
    "The audit includes every supplied participant and event row; no patient or event sampling is used.",
]


def read_table(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        try:
            fields = next(reader)
        except StopIteration:
            raise ValueError(f"Empty table: {path}") from None
        if len(set(fields)) != len(fields):
            raise ValueError(f"Duplicate column names in {path}")
        rows = []
        for number, row in enumerate(reader, 2):
            if len(row) != len(fields):
                raise ValueError(f"{path}:{number}: {len(row)} cells, expected {len(fields)}")
            rows.append(["" if value in MISSING else value for value in row])
    return fields, rows


def read_metadata(path: Path):
    if not path.exists():
        return {}
    if path.suffix == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError("Metadata must be a JSON object")
        return data
    fields, rows = read_table(path)
    if not {"key", "value"}.issubset(fields):
        raise ValueError("metadata.tsv needs key and value columns")
    key, value = fields.index("key"), fields.index("value")
    data = {}
    for row in rows:
        if row[key] in data:
            raise ValueError(f"Duplicate metadata key: {row[key]}")
        try:
            data[row[key]] = json.loads(row[value])
        except (ValueError, TypeError):
            data[row[key]] = row[value]
    return data


@lru_cache(maxsize=65536)
def ordinal(value):
    """Only genuine ISO calendar dates; never interpret UKB special values."""
    if not value or value in SPECIAL_DATES:
        return None
    try:
        parsed = date.fromisoformat(value)
    except (TypeError, ValueError):
        return None
    return parsed.toordinal() if parsed.isoformat() == value else None


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"file": path.name, "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def pair_metrics(a, b):
    diffs = [x - y for x, y in zip(a, b) if x is not None and y is not None]
    only_a = sum(x is not None and y is None for x, y in zip(a, b))
    only_b = sum(x is None and y is not None for x, y in zip(a, b))
    exact = sum(x == 0 for x in diffs)
    return {
        "both": len(diffs), "same_day": exact,
        "same_day_pct": round(100 * exact / len(diffs), 2) if diffs else None,
        "a_earlier": sum(x < 0 for x in diffs), "b_earlier": sum(x > 0 for x in diffs),
        "a_only": only_a, "b_only": only_b,
        "neither_date": len(a) - len(diffs) - only_a - only_b,
        "median_days": statistics.median(diffs) if diffs else None,
        "mean_days": round(statistics.mean(diffs), 2) if diffs else None,
        "min_days": min(diffs) if diffs else None,
        "max_days": max(diffs) if diffs else None,
        "within_30d": sum(abs(x) <= 30 for x in diffs),
        "within_365d": sum(abs(x) <= 365 for x in diffs),
    }


def make_summary(patient_fields, patients, event_fields, events, metadata):
    pi = {name: i for i, name in enumerate(patient_fields)}
    ei = {name: i for i, name in enumerate(event_fields)}
    if "eid" not in pi or not {"eid", "source", "date"}.issubset(ei):
        raise ValueError("patients needs eid; events needs eid, source and date")
    ids = [r[pi["eid"]] for r in patients]
    if any(not x for x in ids) or len(set(ids)) != len(ids):
        raise ValueError("patients contains missing or duplicate eid values")
    id_set = set(ids)
    unknown_ids = {r[ei["eid"]] for r in events} - id_set
    if unknown_ids:
        raise ValueError(f"Event eids absent from patients: {len(unknown_ids)}")
    sources = [name for name in SOURCES if name in pi]
    sources += [name for name in patient_fields if name.startswith("computed_")]
    dates = {name: [ordinal(row[pi[name]]) for row in patients] for name in sources}
    event_counts, undated, event_people, granularities = Counter(), Counter(), defaultdict(set), defaultdict(Counter)
    statuses = Counter()
    for row in events:
        source = row[ei["source"]]
        event_counts[source] += 1
        event_people[source].add(row[ei["eid"]])
        if ordinal(row[ei["date"]]) is None:
            undated[source] += 1
        if "granularity" in ei:
            granularities[source][row[ei["granularity"]]] += 1
        if "date_status" in ei:
            statuses[row[ei["date_status"]]] += 1
    source_stats = []
    for source in list(dict.fromkeys([*sources, *event_counts])):
        valid = [x for x in dates.get(source, []) if x is not None]
        source_stats.append({"source": source, "dated_patients": len(valid),
                             "patients_with_event_rows": len(event_people[source]),
                             "event_rows": event_counts[source], "undated_event_rows": undated[source],
                             "date_min": date.fromordinal(min(valid)).isoformat() if valid else "",
                             "date_max": date.fromordinal(max(valid)).isoformat() if valid else "",
                             "granularity": dict(granularities[source])})
    pairs = {a: {} for a in sources}
    for i, a in enumerate(sources):
        for b in sources[i:]:
            metric = pair_metrics(dates[a], dates[b])
            pairs[a][b] = metric
            if a != b:
                reverse = dict(metric)
                for x, y in (("a_earlier", "b_earlier"), ("a_only", "b_only")):
                    reverse[x], reverse[y] = metric[y], metric[x]
                for key in ("median_days", "mean_days"):
                    reverse[key] = -metric[key] if metric[key] is not None else None
                reverse["min_days"] = -metric["max_days"] if metric["max_days"] is not None else None
                reverse["max_days"] = -metric["min_days"] if metric["min_days"] is not None else None
                pairs[b][a] = reverse
    ref_dates = dates.get("ref", [None] * len(patients))
    candidates = []
    # Transparent hypotheses only; the official source mapping/release may differ.
    accumulated = [None] * len(patients)
    gp_sources = ["gp"] if "gp" in dates else ["gp2", "gp3"]
    stages = [("ICD10", ["icd10"]), ("+ ICD9", ["icd9"]), ("+ SRD", ["srd"]),
              ("+ GP", gp_sources), ("+ death", ["death"])]
    used = []
    for label, additions in stages:
        available = [x for x in additions if x in dates]
        used.extend(available)
        for name in available:
            accumulated = [y if x is None else x if y is None else min(x, y)
                           for x, y in zip(accumulated, dates[name])]
        candidates.append({"label": label, "sources": list(used),
                           "dated_patients": sum(x is not None for x in accumulated),
                           **pair_metrics(accumulated, ref_dates)})
    if "reconstructed_ref" in dates:
        candidates.append({"label": "Published UKB maps (candidate)",
                           "sources": [x for x in sources if x.startswith("official_")],
                           "dated_patients": sum(x is not None for x in dates["reconstructed_ref"]),
                           **pair_metrics(dates["reconstructed_ref"], ref_dates)})
    explanations = Counter(row[pi["ref_explanation"]] or "unclassified" for row in patients) if "ref_explanation" in pi else Counter()
    ref_statuses = Counter(row[pi["ref_status"]] or "missing" for row in patients) if "ref_status" in pi else Counter()
    ref_sources = Counter(row[pi["ref_source"]] or "unavailable" for row in patients) if "ref_source" in pi else Counter()
    limitations = list(DEFAULT_LIMITATIONS)
    for key in ("limitations", "audit_limitations"):
        supplied = metadata.get(key, [])
        if isinstance(supplied, str):
            supplied = [supplied]
        if isinstance(supplied, list):
            limitations.extend(str(x) for x in supplied if str(x) not in limitations)
    return {"Y": metadata.get("Y", metadata.get("trait", "disease")),
            "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "n_patients": len(patients), "n_events": len(events),
            "sources": sources, "source_stats": source_stats, "pairs": pairs,
            "reconstruction_candidates": candidates, "ref_explanations": dict(explanations),
            "ref_statuses": dict(ref_statuses), "ref_sources": dict(ref_sources),
            "event_date_statuses": dict(statuses), "metadata": metadata,
            "limitations": limitations,
            "difference_definition": "date(A) - date(B), in calendar days; negative means A earlier",
            "agreement_denominator": "Participants with usable dates in BOTH sources",
            "population_definition": "Exactly the supplied patients table; consult metadata for cohort and case-union counts"}


def build(input_dir, output=None, trait=None, metadata_path=None):
    input_dir = Path(input_dir)
    patient_path, event_path = input_dir / "patients.tsv.gz", input_dir / "events.tsv.gz"
    pf, patients = read_table(patient_path)
    ef, events = read_table(event_path)
    metadata_path = Path(metadata_path) if metadata_path else input_dir / "metadata.tsv"
    metadata = read_metadata(metadata_path)
    if trait:
        metadata["Y"] = trait
    summary = make_summary(pf, patients, ef, events, metadata)
    summary["input_files"] = [fingerprint(p) for p in (patient_path, event_path, metadata_path) if p.exists()]
    payload = {"summary": summary, "patient_fields": pf, "patients": patients,
               "event_fields": ef, "events": events}
    raw = json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")
    packed = base64.b64encode(gzip.compress(raw, compresslevel=6, mtime=0)).decode("ascii")
    template = Path(__file__).with_name("date_check.html").read_text(encoding="utf-8")
    document = template.replace("__TITLE__", html.escape(str(summary["Y"]))).replace("__PAYLOAD_GZIP__", packed)
    output = Path(output) if output else input_dir / "date-check.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")
    summary_path = output.with_name("summary.json")
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Report: {output}\nPatients: {len(patients):,}; event rows: {len(events):,}; no sampling")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--Y")
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()
    build(args.input_dir, args.output, args.Y, args.metadata)


if __name__ == "__main__":
    main()
