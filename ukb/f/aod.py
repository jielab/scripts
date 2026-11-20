#!/usr/bin/env python3
"""AOD analysis: dates, recorded activity and unsupervised record patterns.

Only NumPy is required beyond Python's standard library. Run via ../aod.sh.
Input is the read-only export from aod.R. No UKB data leave this process.
Calendar dates are never winsorized; clipping is confined to model features.
The patient table contains identifiers; the report payload is aggregate only.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import date, datetime, timezone
from functools import lru_cache
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("AOD requires numpy in UKB_PYTHON. Install numpy in that environment.") from exc

VERSION = "2.0.0"
DAYS_YEAR = 365.25
SPECIAL = {"1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"}
MISSING = {"", "NA", "NaN", "nan", ".", "NULL", "null", "N/A"}
SOURCES = {"icd10": "fod.icd10", "ref": "fod.ref", "icd9": "fod.icd9", "gp": "fod.gp",
           "srd": "fod.srd", "opcs4": "fod.opcs4", "min1": "min1 · \u672c\u5730\u6700\u65e9", "min2": "min2 · \u6620\u5c04\u5019\u9009"}
GROUP = {"icd10": "hospital", "icd9": "hospital", "gp": "gp", "gp2": "gp", "gp3": "gp"}
LEVEL_LABEL = {0: "first_date_only", 1: "first_per_code", 2: "cached_record_days", 3: "linked_records"}
LAG_EDGES = [-3653, -366, -31, -1, 0, 30, 365, 3652]
LAG_LABELS = ["≤−10\u5e74", "−10\u81f3−1\u5e74", "−1\u5e74\u81f3−31\u5929", "−30\u81f3−1\u5929", "\u540c\u65e5", "1\u81f330\u5929", "31\u5929\u81f31\u5e74", "1\u81f310\u5e74", ">10\u5e74"]
csv.field_size_limit(1 << 28)


def rows(path):
    path = Path(path)
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt", encoding="utf-8-sig", newline="") as inp:
        reader = csv.DictReader(inp, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"Missing/duplicate header in {path}")
        for n, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Malformed row {n} in {path}")
            yield row


def write_tsv(path, fields, data):
    path = Path(path)
    with (gzip.open if path.suffix == ".gz" else open)(path, "wt", encoding="utf-8", newline="") as out:
        wr = csv.DictWriter(out, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        wr.writeheader()
        for row in data:
            wr.writerow({k: "" if row.get(k) is None else row.get(k, "") for k in fields})


def read_meta(path):
    out = {}
    for row in rows(path):
        key = row["key"]
        if key in out:
            raise ValueError(f"Duplicate metadata key: {key}")
        try:
            out[key] = json.loads(row["value"])
        except (ValueError, TypeError):
            out[key] = row["value"]
    return out


def truth(x):
    return x is True or str(x).lower() in {"1", "true", "t", "yes"}


@lru_cache(maxsize=100000)
def ordinal(x):
    if not x or x in MISSING or x in SPECIAL or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", str(x)):
        return None
    try:
        return date.fromisoformat(x).toordinal()
    except ValueError:
        return None


def iso(x):
    return date.fromordinal(int(x)).isoformat() if x is not None else None


def num(x):
    try:
        f = float(x)
        return f if math.isfinite(f) else None
    except (ValueError, TypeError):
        return None


def clean_json(x):
    if isinstance(x, dict):
        return {str(k): clean_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [clean_json(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return round(float(x), 8) if math.isfinite(x) else None
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def save_json(path, value):
    Path(path).write_text(json.dumps(clean_json(value), ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n", encoding="utf-8")


def percent(a, b):
    return 100.0 * a / b if b else None


def med(vals):
    vals = [v for v in vals if v is not None and math.isfinite(v)]
    return float(np.median(vals)) if vals else None


def hist(vals, edges, labels):
    counts = [0] * (len(edges) + 1)
    for val in vals:
        if val is not None:
            counts[int(np.searchsorted(edges, val, side="left"))] += 1
    return [{"label": label, "n": n} for label, n in zip(labels, counts)]


def pair_metrics(a, b, baseline):
    diffs = [x - y for x, y in zip(a, b) if x is not None and y is not None]
    exact = sum(v == 0 for v in diffs)
    ap = bp = nb = 0
    for x, y, z in zip(a, b, baseline):
        if x is not None and y is not None and z is not None:
            nb += 1
            ap += x <= z < y
            bp += y <= z < x
    return {"both": len(diffs), "same_day": exact, "same_day_pct": percent(exact, len(diffs)),
            "within_30d_pct": percent(sum(abs(d) <= 30 for d in diffs), len(diffs)),
            "within_365d_pct": percent(sum(abs(d) <= 365 for d in diffs), len(diffs)),
            "a_only": sum(x is not None and y is None for x, y in zip(a, b)),
            "b_only": sum(y is not None and x is None for x, y in zip(a, b)),
            "neither": sum(x is None and y is None for x, y in zip(a, b)),
            "median_days": med(diffs), "a_earlier": sum(d < 0 for d in diffs), "b_earlier": sum(d > 0 for d in diffs),
            "histogram": hist(diffs, LAG_EDGES, LAG_LABELS),
            "baseline": {"both_baseline": nb, "n_discordant": ap + bp, "a_pre_b_post": ap, "b_pre_a_post": bp}}


def common_comparison(min1, min2, ref):
    """Compare both candidates against ref in exactly the same people."""
    triples = [(a, b, r) for a, b, r in zip(min1, min2, ref)
               if a is not None and b is not None and r is not None]
    n = len(triples)
    old = [abs(a - r) for a, _, r in triples]
    new = [abs(b - r) for _, b, r in triples]
    result = {"n": n,
              "gained": sum(a != 0 and b == 0 for a, b in zip(old, new)),
              "lost": sum(a == 0 and b != 0 for a, b in zip(old, new)),
              "both_same": sum(a == 0 and b == 0 for a, b in zip(old, new)),
              "neither_same": sum(a != 0 and b != 0 for a, b in zip(old, new))}
    for source, diffs in (("min1", old), ("min2", new)):
        result[source + "_same_day"] = sum(d == 0 for d in diffs)
        for label, limit in (("same_day", 0), ("within_30d", 30), ("within_365d", 365)):
            result[source + "_" + label + "_pct"] = percent(sum(d <= limit for d in diffs), n)
    result["delta_pct_points"] = percent(result["gained"] - result["lost"], n)
    return result


def reference_audit(meta):
    """Keep extraction completeness visible without inferring clinical validity."""
    def fields(key):
        value = meta.get(key, [])
        if isinstance(value, str):
            value = [value] if value else []
        return [str(v) for v in value if v is not None and str(v)] if isinstance(value, list) else []

    expected = fields("reference_fields_expected")
    missing = fields("reference_fields_missing")
    mode = meta.get("reference_input_mode", "unknown")
    if mode not in {"raw", "saved_only", "unavailable"}:
        mode = "unknown"
    complete = "unavailable" if mode == "unavailable" else ("partial" if missing else
                ("complete" if mode == "raw" and expected else "unverified"))
    notes = []
    if missing:
        notes.append(f"ref \u8f93\u5165\u7f3a\u5c11 {len(missing)} \u4e2a\u9884\u671f\u5b57\u6bb5；\u5f53\u524d\u4e00\u81f4\u7387\u4ec5\u57fa\u4e8e\u5df2\u63d0\u4f9b\u5b57\u6bb5，\u4e0d\u80fd\u4ee3\u8868\u5b8c\u6574\u5b9a\u4e49。")
    if mode == "saved_only":
        notes.append("ref \u4ec5\u4f7f\u7528\u65e2\u5b58\u7f13\u5b58；\u66f4\u65b0\u5b9a\u4e49\u540e\u9700\u91cd\u65b0\u63d0\u53d6\u539f\u59cb\u5b57\u6bb5\u5e76\u751f\u6210 fod_ref，\u624d\u80fd\u6838\u5bf9\u65b0\u8303\u56f4。")
    elif mode == "unavailable":
        notes.append("\u5f53\u524d\u8f93\u5165\u6ca1\u6709\u53ef\u7528 ref \u5b57\u6bb5\u6216\u7f13\u5b58，\u4e0d\u80fd\u8ba1\u7b97\u53c2\u8003\u65e5\u671f\u4e00\u81f4\u7387。")
    elif mode == "unknown":
        notes.append("\u5bfc\u51fa\u672a\u8bb0\u5f55 ref \u8f93\u5165\u6a21\u5f0f；\u65e0\u6cd5\u6838\u5b9e\u53c2\u8003\u65e5\u671f\u6765\u81ea\u539f\u59cb\u5b57\u6bb5\u8fd8\u662f\u7f13\u5b58。")
    elif not expected:
        notes.append("ref \u4f7f\u7528\u539f\u59cb\u5b57\u6bb5，\u4f46\u5bfc\u51fa\u672a\u8bb0\u5f55\u9884\u671f\u5b57\u6bb5\u5217\u8868，\u5b8c\u6574\u6027\u5c1a\u672a\u6838\u5b9e。")
    return {"input_mode": mode, "expected_fields": expected, "missing_fields": missing,
            "completeness": complete, "notes": notes}


def resolve_asof(meta, requested=None):
    """An already truncated export cannot be extended by changing a CLI date."""
    exported = meta.get("as_of")
    if not isinstance(exported, str) or ordinal(exported) is None:
        raise ValueError("Extracted metadata needs a valid as_of cutoff; re-extract with aod.sh from the original RDS/records")
    chosen = requested or exported
    if not isinstance(chosen, str) or ordinal(chosen) is None:
        raise ValueError("Invalid as-of date")
    if ordinal(chosen) > ordinal(exported):
        raise ValueError(f"Analysis --as-of {chosen} is later than export cutoff {exported}; re-extract the original RDS/records for a later cutoff")
    return chosen


def read_patients(path, asof, mode, followup_end=None):
    patients = {}
    for p in rows(path):
        eid = p.get("eid", "")
        if eid in MISSING or eid in patients:
            raise ValueError("patients.tsv.gz requires unique, nonmissing eid")
        p["birth"] = ordinal(p.get("date_birth", ""))
        p["baseline"] = ordinal(p.get("date_assess", ""))
        p["death"] = ordinal(p.get("date_death", ""))
        p["lost"] = ordinal(p.get("date_lost", ""))
        # FOD agreement is an audit of all valid dates. The baseline restriction
        # applies only to event-derived features, never to the dates being audited.
        p["audit_cutoff"] = min([asof] + [p[k] for k in ("death", "lost") if p[k] is not None])
        if followup_end is None:
            raise ValueError("Missing shared date_follow_end for clinical features")
        p["clinical_end"] = min(p["audit_cutoff"], followup_end)
        p["cutoff"] = p["clinical_end"]
        if mode == "baseline":
            p["cutoff"] = min(p["cutoff"], p["baseline"]) if p["baseline"] is not None else None
        for key in SOURCES:
            val = ordinal(p.get(key, ""))
            if val is not None and (val > p["audit_cutoff"] or (p["birth"] is not None and val < p["birth"])):
                val = None
            p[key] = val
        p["analysis_case"] = any(p[k] is not None and p[k] <= p["clinical_end"] for k in SOURCES)
        patients[eid] = p
    if not patients:
        raise ValueError("No patient evidence in input")
    return dict(sorted(patients.items()))


def read_coverage(path):
    """Merge only overlapping/adjacent explicit observation intervals."""
    coverage = defaultdict(list)
    if not path:
        return coverage
    path = Path(path)
    with (gzip.open if path.suffix == ".gz" else open)(path, "rt", encoding="utf-8-sig", newline="") as inp:
        fields = next(csv.reader(inp, delimiter="\t"), [])
    if not {"eid", "source", "start", "end"} <= set(fields):
        raise ValueError("Coverage TSV needs eid,source,start,end, including when it has no data rows")
    for row in rows(path):
        if not {"eid", "source", "start", "end"} <= set(row):
            raise ValueError("Coverage TSV needs eid,source,start,end; multiple disjoint intervals are allowed")
        eid = row["eid"]
        src = GROUP.get(row["source"], row["source"])
        a, b = ordinal(row["start"]), ordinal(row["end"])
        if eid in MISSING or src not in {"hospital", "gp"} or a is None or b is None or b < a:
            raise ValueError("Invalid coverage row; use exact participant IDs, hospital/gp, and ISO start <= end")
        coverage[eid, src].append((a, b))
    for key, intervals in coverage.items():
        merged = []
        for a, b in sorted(intervals):
            if merged and a <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else:
                merged.append((a, b))
        coverage[key] = merged
    return coverage


def coverage_days(intervals, start, end):
    return sum(max(0, min(b, end) - max(a, start) + 1) for a, b in intervals) if end >= start else 0


def event_level(source, granularity):
    if granularity in {"record", "episode"}:
        return 3
    if source.startswith("gp") and granularity in {"record_day_cached", "cached_aod"}:
        return 2
    if granularity in {"cached_aod", "code_first_aod", "summary_first_per_code"}:
        return 1
    if granularity == "fod_only":
        return 0
    return -1


def read_events(path, patients, meta):
    banks = {s: {i: defaultdict(set) for i in range(4)} for s in ("hospital", "gp")}
    present = {s: set() for s in banks}
    qc = Counter()
    procedure = defaultdict(set)
    all_local = defaultdict(set)
    for row in rows(path):
        eid = row.get("eid", "")
        if eid not in patients:
            raise ValueError(f"Event participant absent from patients table: {eid}")
        p = patients[eid]
        src = row.get("source", "")
        level = event_level(src, row.get("granularity", ""))
        local = truth(row.get("match_local", ""))
        group = GROUP.get(src)
        raw = row.get("date_raw", row.get("date", ""))
        d = ordinal(row.get("date", ""))
        status = row.get("date_status", "valid")
        if status != "valid" or d is None:
            qc[status if status != "valid" else "invalid_date"] += 1
            continue
        if raw in SPECIAL:
            qc["placeholder"] += 1
            continue
        if p["birth"] is not None and d < p["birth"]:
            qc["before_birth"] += 1
            continue
        if p["death"] is not None and d > p["death"]:
            qc["after_death"] += 1
            continue
        if d > p["clinical_end"]:
            qc["after_analysis_cutoff"] += 1
            continue
        if not local:
            continue
        if src != "opcs4_context":
            p["analysis_case"] = True
        if group and level >= 0:
            present[group].add(level)
        if "date_field=epiend" in row.get("origin", "") or "date_field=disdate" in row.get("origin", ""):
            qc["hospital_end_date_proxy"] += 1
        if src in {"opcs4", "opcs4_context"}:
            procedure[eid].add(d)
        if src in GROUP or src in {"srd", "opcs4", "death"}:
            all_local[eid].add(d)
        if group and level >= 0:
            banks[group][level][eid].add(d)
    # Global source precedence avoids mixing first-per-code caches into people
    # without a record in an otherwise available linked-record extract.
    definitions = meta.get("definitions", {}).get("sources", {})
    hospital_defined = not definitions or any(truth(definitions.get(s, {}).get("local_defined")) for s in ("icd10", "icd9"))
    gp_defined = not definitions or any(truth(definitions.get(s, {}).get("local_defined")) for s in ("gp2", "gp3"))
    if hospital_defined and truth(meta.get("raw_hospital_episode_dates_available")) and any(truth(meta.get("raw_" + s + "_available")) for s in ("icd10", "icd9")):
        present["hospital"].add(3)
    if gp_defined and any(truth(meta.get("raw_" + s + "_available")) for s in ("gp2", "gp3")):
        present["gp"].add(3)
    levels = {s: max(v) if v else -1 for s, v in present.items()}
    selected = {s: banks[s][level] if level >= 0 else {} for s, level in levels.items()}
    return selected, levels, qc, procedure, all_local


def bouts(dates, gap):
    out = []
    for d in dates:
        if not out or d - out[-1][-1] > gap:
            out.append([d])
        else:
            out[-1].append(d)
    return out


PATTERN_LABELS = {
    "quiet_window": "\u7a97\u5185\u672a\u89c1\u540e\u7eed\u8bb0\u5f55",
    "single_bout_quiet": "\u5355\u6bb5\u8bb0\u5f55\u540e\u957f\u671f\u672a\u518d\u8bb0\u5f55",
    "early_limited": "\u65e9\u671f\u96c6\u4e2d\u540e\u51cf\u5c11",
    "increasing": "\u540e\u671f\u8bb0\u5f55\u589e\u52a0",
    "persistent": "\u8de8\u671f\u53cd\u590d\u8bb0\u5f55",
    "intermittent": "\u95f4\u6b47\u518d\u73b0",
    "recent_limited": "\u8fd1\u671f\u96c6\u4e2d\u8bb0\u5f55",
}


def split_name(eid, seed):
    # Stable across input order, case additions and outcome updates. Identifier
    # hashes are used only internally and never written to aggregate reports.
    h = hashlib.blake2b(f"{seed}:{eid}".encode(), digest_size=8).digest()
    return "test" if int.from_bytes(h, "big") % 4 == 0 else "train"


def pattern_for(f, gap=None):
    """An explicit record-course rule, not a disease-remission classifier."""
    if f["cluster_reason"] != "eligible" and f.get("pattern_reason") != "classified":
        return None
    n = f["n_days"]
    if n == 0:
        return "quiet_window"
    nb = len(bouts(f["_days"], gap)) if gap is not None else f["n_bouts"]
    quiet = f["time_since_last_days"]
    # Require an appreciable quiet tail; one late observation is not recovered.
    quiet_limit = min(365, max(90, int((f["_end"] - f["_start"] + 1) * .5)))
    if nb == 1 and quiet is not None and quiet >= quiet_limit:
        return "single_bout_quiet"
    delta = f["late_minus_early"]
    if delta is not None and delta >= .25 and f["late_active_fraction"] >= .5:
        return "increasing"
    if (f["post_index_active_fraction"] is not None and f["post_index_active_fraction"] >= .5
            and f["late_active_fraction"] >= .5 and n >= 2):
        return "persistent"
    if quiet is not None and quiet >= quiet_limit:
        return "early_limited"
    if nb >= 2:
        return "intermittent"
    return "recent_limited"


def feature_rows(patients, selected, levels, coverage, args):
    window = int(round((args.lookback_years if args.mode == "baseline" else args.landmark_years) * DAYS_YEAR))
    nbins = max(4, min(40, int(round(window / DAYS_YEAR * 4))))
    requested = ["hospital", "gp"] if args.trajectory_source == "combined" else [args.trajectory_source]
    used_sources = [s for s in requested if levels.get(s, -1) >= 1]
    partial_sources = [s for s in getattr(args, "partial_sources", []) if s in requested]
    can_measure_encounters = (set(used_sources) == set(requested)
                              and all(levels[s] >= 2 for s in used_sources) and not partial_sources)
    feats = []
    for eid, p in patients.items():
        full = {s: sorted(d for d in selected[s].get(eid, ()) if p["cutoff"] is not None and d <= p["cutoff"])
                if levels[s] >= 1 else [] for s in selected}
        observed = sorted(set().union(*(set(full[s]) for s in used_sources))) if used_sources else []
        anchor = p["baseline"] if args.mode == "baseline" else (observed[0] if observed else None)
        start = anchor - window + 1 if args.mode == "baseline" and anchor is not None else anchor
        end = anchor if args.mode == "baseline" else (start + window - 1 if start is not None else None)
        fod_candidates = [p[k] for k in ("min1", "ref", "icd10", "icd9", "gp", "srd")
                          if p.get(k) is not None and end is not None and p[k] <= min(end, p["clinical_end"])]
        first_fod = min(fod_candidates) if fod_candidates else (observed[0] if observed else None)
        f = {"eid": eid, "window_start": iso(start), "window_end": iso(end), "index_date": iso(observed[0]) if observed else None,
             "n_days_total": len(observed), "record_scope": "+".join(used_sources), "hospital_mode": LEVEL_LABEL.get(levels["hospital"], "unavailable"),
             "gp_mode": LEVEL_LABEL.get(levels["gp"], "unavailable"), "cluster": None, "cluster_reason": "no_clinical_aod",
             "aod": None, "aod_reason": "no_clinical_aod", "adu": num(p.get("adu")), "_start": start, "_end": end,
             "_trajectory": np.full(nbins, np.nan), "_observed_bins": np.zeros(nbins, dtype=bool),
             "_hospital_trajectory": None, "_gp_trajectory": None, "_days": [],
             "pattern": None, "pattern_evidence": "unavailable", "pattern_reason": "no_clinical_aod",
             "assignment_margin": None, "split": "apply" if args.apply_model else split_name(eid, args.seed),
             "clinical_end": iso(p["clinical_end"]), "analysis_case": p["analysis_case"],
             "age_at_index": (observed[0] - p["birth"]) / DAYS_YEAR if observed and p["birth"] is not None else None,
             "sex": p.get("sex") if p.get("sex") not in MISSING else None, "landmark_date": iso(end),
             "age_at_landmark": (end - p["birth"]) / DAYS_YEAR if end is not None and p["birth"] is not None else None,
             "source_duration_years": (end - observed[0]) / DAYS_YEAR if end is not None and observed else None,
             "fod_for_prediction": iso(first_fod),
             "disease_duration_years": (end - first_fod) / DAYS_YEAR if end is not None and first_fod is not None else None,
             "horizon_complete": False, "coverage_complete": False, "recurrence_capable": can_measure_encounters}
        for key in ("n_days", "n_bouts", "span_years", "median_gap_days", "max_gap_days", "active_fraction", "hospital_fraction", "gp_fraction",
                    "event_rate", "episode_rate", "hospital_rate", "gp_rate", "recent_rate", "time_since_last_days", "late_minus_early",
                    "early_active_fraction", "late_active_fraction", "post_index_active_fraction", "early_isolate", "n_days_without_isolate",
                    "n_hospital_days", "n_gp_days", "hospital_observed_years", "gp_observed_years"):
            f[key] = None
        if not used_sources or not observed or start is None or end is None or p["cutoff"] is None:
            reason = ("missing_baseline" if args.mode == "baseline" and p["baseline"] is None else
                      "post_cutoff_or_undated_evidence_only" if not p["analysis_case"] else
                      "no_selected_source_aod_before_anchor")
            f["cluster_reason"] = f["aod_reason"] = f["pattern_reason"] = reason
            feats.append(f)
            continue
        f["horizon_complete"] = p["cutoff"] >= end
        win = {s: [d for d in ds if start <= d <= min(end, p["cutoff"])] for s, ds in full.items()}
        days = sorted(set().union(*(set(win[s]) for s in used_sources)))
        f["_days"] = days
        f["n_days"] = len(days)
        f["n_bouts"] = len(bouts(days, args.episode_gap))
        f["span_years"] = (days[-1] - days[0]) / DAYS_YEAR if days else 0.0
        gaps = np.diff(days).tolist()
        f["median_gap_days"] = med(gaps)
        f["max_gap_days"] = max(gaps) if gaps else None
        # Last observed date can precede the baseline lookback window.
        previous = [d for d in observed if d <= min(end, p["cutoff"])]
        f["time_since_last_days"] = min(end, p["cutoff"]) - previous[-1] if previous else None
        f["n_hospital_days"] = len(win["hospital"]) if levels["hospital"] >= 1 else None
        f["n_gp_days"] = len(win["gp"]) if levels["gp"] >= 1 else None
        source_days = sum(len(win[s]) for s in used_sources)
        both = set(used_sources) == {"hospital", "gp"}
        f["hospital_fraction"] = len(win["hospital"]) / source_days if source_days and both else None
        f["gp_fraction"] = len(win["gp"]) / source_days if source_days and both else None
        # Plot only complete temporal bins; explicit coverage gaps stay missing.
        bounds = [start + int(math.ceil(j * window / nbins)) for j in range(nbins + 1)]
        obs_mask = np.array([bounds[j + 1] - 1 <= p["cutoff"] for j in range(nbins)])
        if args.coverage:
            obs_mask &= np.array([all(coverage_days(coverage.get((eid, src), []), bounds[j], bounds[j+1]-1)
                                     == bounds[j+1]-bounds[j] for src in requested) for j in range(nbins)])
        raw_traj = np.zeros(nbins)
        for d in days:
            raw_traj[min(nbins - 1, int((d - start) / window * nbins))] = 1
        f["_trajectory"] = np.where(obs_mask, raw_traj, np.nan)
        f["_observed_bins"] = obs_mask
        for src in ("hospital", "gp"):
            if levels[src] < 1:
                continue
            raw = np.zeros(nbins)
            for d in win[src]:
                raw[min(nbins - 1, int((d - start) / window * nbins))] = 1
            smask = np.array([bounds[j+1]-1 <= p["cutoff"] for j in range(nbins)])
            if args.coverage:
                smask &= np.array([coverage_days(coverage.get((eid, src), []), bounds[j], bounds[j+1]-1)
                                   == bounds[j+1]-bounds[j] for j in range(nbins)])
            f["_" + src + "_trajectory"] = np.where(smask, raw, np.nan)
        confirmed = set(used_sources) == set(requested)
        for src in used_sources:
            intervals = coverage.get((eid, src), [])
            obsdays = coverage_days(intervals, start, min(end, p["cutoff"]))
            f[src + "_observed_years"] = obsdays / DAYS_YEAR if intervals else None
            if intervals and obsdays and levels[src] >= 2 and src not in partial_sources:
                n_observed = sum(any(a <= d <= b for a, b in intervals) for d in win[src])
                f[src + "_rate"] = n_observed / (obsdays / DAYS_YEAR)
            confirmed = confirmed and bool(intervals) and obsdays == window
        f["coverage_complete"] = confirmed and f["horizon_complete"]
        complete_for_pattern = f["horizon_complete"] and (not args.coverage or f["coverage_complete"])
        if complete_for_pattern:
            half = nbins // 2
            skip = 1 if args.mode == "history" else 0
            f["active_fraction"] = float(raw_traj.mean())
            f["post_index_active_fraction"] = float(raw_traj[skip:].mean())
            f["early_active_fraction"] = float(raw_traj[skip:half].mean())
            f["late_active_fraction"] = float(raw_traj[half:].mean())
            f["late_minus_early"] = f["late_active_fraction"] - f["early_active_fraction"]
        else:
            f["time_since_last_days"] = None
        first = observed[0]
        corroborated = all(any(abs(d - first) <= 365 for d in full[src]) for src in ("hospital", "gp"))
        suspect = len(observed) > 1 and observed[1] - first > 5 * DAYS_YEAR and not corroborated
        f["early_isolate"] = int(suspect)
        altered = observed[1:] if suspect else observed
        # Sensitivity never borrows dates after the clinical cutoff.
        if args.mode == "history" and altered:
            f["n_days_without_isolate"] = sum(altered[0] <= d < altered[0] + window for d in altered)
        else:
            f["n_days_without_isolate"] = sum(start <= d <= end for d in altered)
        if f["coverage_complete"] and can_measure_encounters:
            f["event_rate"] = len(days) / (window / DAYS_YEAR)
            f["episode_rate"] = f["n_bouts"] / (window / DAYS_YEAR)
            last_year = min(window, int(round(DAYS_YEAR)))
            f["recent_rate"] = sum(d > end - last_year for d in days) / (last_year / DAYS_YEAR)
            f["aod_reason"] = "eligible"
        else:
            f["aod_reason"] = ("partial_code_source_extract" if partial_sources else
                               "first_per_code_not_complete_records" if not can_measure_encounters else
                               "unverified_or_incomplete_coverage")
        if not f["horizon_complete"]:
            f["cluster_reason"] = ("death_before_window_end" if p["death"] is not None and p["death"] == p["clinical_end"] else
                                   "loss_before_window_end" if p["lost"] is not None and p["lost"] == p["clinical_end"] else
                                   "administrative_cutoff_before_window_end")
        elif args.coverage and not f["coverage_complete"]:
            f["cluster_reason"] = "incomplete_source_coverage"
        elif not days and not f["coverage_complete"]:
            f["cluster_reason"] = "no_dates_in_window_without_confirmed_coverage"
        else:
            f["cluster_reason"] = "eligible"
        f["pattern"] = pattern_for(f)
        f["pattern_reason"] = "classified" if f["pattern"] else f["cluster_reason"]
        f["pattern_evidence"] = ("repeated_records_confirmed_coverage" if can_measure_encounters and f["coverage_complete"] else
                                 "repeated_records_unknown_coverage" if can_measure_encounters else
                                 "first_per_code_history" if any(levels[src] == 1 for src in used_sources) else "partial_source_history")
        feats.append(f)
    return feats, used_sources, nbins, window


def percentile_values(values):
    """Average-rank percentiles, ties equal. Not a clinical reference scale."""
    a = np.asarray(values, dtype=float)
    order = np.argsort(a, kind="stable")
    ranked = np.empty(len(a))
    i = 0
    while i < len(a):
        j = i + 1
        while j < len(a) and a[order[j]] == a[order[i]]:
            j += 1
        ranked[order[i:j]] = (i + j) / (2 * len(a)) * 100
        i = j
    return ranked


def score_burden(feats, window, saved=None):
    eligible = [f for f in feats if f["aod_reason"] == "eligible"]
    training = [f for f in eligible if f["split"] == "train"]
    model = {"status": "unavailable", "n": 0, "components": ["log1p_event_rate", "post_index_active_fraction", "span_fraction"],
             "weights": [1 / 3] * 3, "method": "mean of three percentiles against frozen training-case distributions; exploratory recorded activity, not severity"}
    make = lambda fs: np.array([[math.log1p(f["event_rate"]), f["post_index_active_fraction"],
                                f["span_years"] * DAYS_YEAR / window] for f in fs], dtype=float).reshape((-1, 3)).round(8)
    if saved:
        refs = saved.get("reference_values", [])
        if not refs:
            for f in eligible:
                f["aod_reason"] = "training_score_unavailable"
            return model
    elif len(training) >= 20:
        xfit = make(training)
        refs = [np.sort(xfit[:, j]).tolist() for j in range(3)]
    else:
        for f in eligible:
            f["aod_reason"] = "fewer_than_20_training_cases"
        return model
    if len(refs) != 3 or any(not r or not np.all(np.isfinite(r)) for r in refs):
        raise ValueError("Invalid frozen burden reference distributions")
    x = make(eligible)
    if len(x):
        ranked = np.column_stack([(np.searchsorted(r, x[:, j], "left") + np.searchsorted(r, x[:, j], "right"))
                                  / (2 * len(r)) * 100 for j, r in enumerate(refs)])
        for f, value in zip(eligible, ranked.mean(axis=1)):
            f["aod"] = float(value)
            f["aod_reason"] = "scored"
    model.update(status="exploratory", reference_values=refs, reference_n=len(refs[0]),
                 n=len(eligible), scope=eligible[0]["record_scope"] if eligible else saved.get("scope", "unavailable"),
                 fit_partition="train", applied_frozen=bool(saved))
    return model


def assign(x, centers):
    labels = np.empty(len(x), dtype=int)
    dist = np.empty(len(x))
    c2 = np.sum(centers * centers, axis=1)
    for start in range(0, len(x), 4096):
        block = x[start:start + 4096]
        ds = np.maximum(0, np.sum(block * block, axis=1)[:, None] + c2[None, :] - 2 * block @ centers.T)
        labels[start:start + len(block)] = ds.argmin(axis=1)
        dist[start:start + len(block)] = ds.min(axis=1)
    return labels, dist


def kmeans(x, k, rng, starts=4, maxiter=80):
    best = None
    for _ in range(starts):
        centers = [x[rng.integers(len(x))].copy()]
        nearest = np.sum((x - centers[0]) ** 2, axis=1)
        for _ in range(1, k):
            if nearest.sum() < 1e-12:
                break
            point = x[rng.choice(len(x), p=nearest / nearest.sum())].copy()
            centers.append(point)
            nearest = np.minimum(nearest, np.sum((x - point) ** 2, axis=1))
        if len(centers) != k:
            continue
        centers = np.asarray(centers)
        for _ in range(maxiter):
            label, ds = assign(x, centers)
            new = np.array([x[label == j].mean(axis=0) if np.any(label == j) else x[int(ds.argmax())] for j in range(k)])
            if np.max(np.abs(new - centers)) < 1e-6:
                centers = new
                break
            centers = new
        label, ds = assign(x, centers)
        score = float(ds.sum())
        if len(set(label)) == k and (best is None or score < best[0]):
            best = (score, centers)
    return None if best is None else best[1]


def silhouette(x, labels):
    groups = sorted(set(labels.tolist()))
    if len(groups) < 2 or len(x) < 3:
        return None
    total = 0.0
    for start in range(0, len(x), 128):
        b = x[start:start + 128]
        ds = np.sqrt(np.maximum(0, np.sum(b * b, axis=1)[:, None] + np.sum(x * x, axis=1)[None, :] - 2 * b @ x.T))
        for j, current in enumerate(labels[start:start + len(b)]):
            same = labels == current
            count = int(same.sum())
            if count <= 1:
                continue
            ds[j, start + j] = 0
            a = ds[j, same].sum() / (count - 1)
            bmean = min(ds[j, labels == other].mean() for other in groups if other != current)
            total += (bmean - a) / max(a, bmean) if max(a, bmean) else 0
    return float(total / len(x))


def adjusted_rand(a, b):
    a = np.asarray(a, dtype=int)
    b = np.asarray(b, dtype=int)
    if len(a) < 2:
        return 1.0
    counts = Counter(zip(a.tolist(), b.tolist()))
    ca, cb = Counter(a.tolist()), Counter(b.tolist())
    choose2 = lambda n: n * (n - 1) / 2
    index = sum(choose2(n) for n in counts.values())
    aa, bb = sum(choose2(n) for n in ca.values()), sum(choose2(n) for n in cb.values())
    expected = aa * bb / choose2(len(a))
    den = (aa + bb) / 2 - expected
    return float((index - expected) / den) if abs(den) > 1e-12 else 1.0


def cluster_features(feats, window, include_source=False):
    names = ["log1p_days", "log1p_bouts", "post_index_active_fraction", "span_fraction", "gap_fraction",
             "late_minus_early", "recency_fraction", "max_gap_fraction"]
    labels = ["\u8bb0\u5f55\u65e5", "\u8bb0\u5f55\u6bb5", "\u540e\u7eed\u6d3b\u8dc3\u5b63\u5ea6", "\u9996\u672b\u8de8\u5ea6", "\u76f8\u90bb\u95f4\u9694", "\u540e\u671f−\u65e9\u671f", "\u8ddd\u672b\u6b21\u8bb0\u5f55", "\u6700\u957f\u95f4\u9694"]
    if include_source:
        names.append("hospital_share")
        labels.append("\u4f4f\u9662\u8bb0\u5f55\u5360\u6bd4")
    x = []
    for f in feats:
        v = [math.log1p(f["n_days"]), math.log1p(f["n_bouts"]), f["post_index_active_fraction"],
             f["span_years"] * DAYS_YEAR / window,
             (f["median_gap_days"] if f["median_gap_days"] is not None else window) / window,
             f["late_minus_early"], min(f["time_since_last_days"], window) / window,
             (f["max_gap_days"] if f["max_gap_days"] is not None else window) / window]
        if include_source:
            v.append(f["hospital_fraction"] if f["hospital_fraction"] is not None else .5)
        x.append(v)
    return np.asarray(x, dtype=float).reshape((-1, len(names))), names, labels


def trajectory_mean(feats, field):
    values = [f[field] for f in feats if f.get(field) is not None]
    if not values:
        return None
    x = np.asarray(values, dtype=float)
    n = np.isfinite(x).sum(axis=0)
    return np.divide(np.nansum(x, axis=0) * 100, n, out=np.full(x.shape[1], np.nan), where=n > 0)


def trajectory_denominators(feats, field="_trajectory"):
    values = [f[field] for f in feats if f.get(field) is not None]
    return np.isfinite(np.asarray(values, dtype=float)).sum(axis=0).tolist() if values else []


def descriptive_label(fs):
    dn, bn = med([f["n_days"] for f in fs]), med([f["n_bouts"] for f in fs])
    delta = med([f["late_minus_early"] for f in fs])
    active = med([f["post_index_active_fraction"] for f in fs])
    if dn == 0:
        return "\u7a97\u5185\u672a\u89c1\u8bb0\u5f55"
    if dn <= 1:
        return "\u5355\u4e00\u65e5\u671f\u4e3a\u4e3b"
    if delta is not None and delta >= .2 and med([f["late_active_fraction"] for f in fs]) >= .5:
        return "\u53cd\u590d\u8bb0\u5f55 · \u540e\u671f\u589e\u52a0"
    if delta is not None and delta <= -.2 and med([f["early_active_fraction"] for f in fs]) >= .5:
        return "\u53cd\u590d\u8bb0\u5f55 · \u540e\u671f\u51cf\u5c11"
    if active is not None and active >= .5:
        return "\u8de8\u671f\u6301\u7eed\u8bb0\u5f55"
    if bn <= 1:
        return "\u5355\u6bb5\u96c6\u4e2d\u8bb0\u5f55"
    return "\u95f4\u6b47\u518d\u73b0 · \u4f4e\u5bc6\u5ea6" if dn <= 4 else "\u95f4\u6b47\u518d\u73b0 · \u591a\u65e5\u671f"


def cluster_profiles(eligible, lab, model, x, names, labels, window):
    profiles = []
    k = int(model["k"])
    label_order = model["label_order"]
    train = np.array([f["split"] == "train" for f in eligible])
    # Labels are descriptions learned in training; apply mode reuses them.
    learned_labels = model.get("cluster_labels", [])
    if not learned_labels:
        learned_labels = [descriptive_label([eligible[i] for i in np.flatnonzero((lab == j) & train)])
                          if np.any((lab == j) & train) else "\u8bb0\u5f55\u6a21\u5f0f" for j in range(k)]
        if k == 1:
            learned_labels = ["\u672a\u53d1\u73b0\u7a33\u5b9a\u5206\u578b"]
        model["cluster_labels"] = learned_labels
    for j in range(k):
        idx = np.flatnonzero(lab == j)
        fs = [eligible[i] for i in idx]
        count_patterns = Counter(f["pattern"] for f in fs if f["pattern"])
        early = med([f["early_active_fraction"] for f in fs])
        late = med([f["late_active_fraction"] for f in fs])
        profile = {"id": label_order[j], "label": learned_labels[j], "n": len(fs), "pct": percent(len(fs), len(eligible)),
                   "train_n": sum(f["split"] == "train" for f in fs), "test_n": sum(f["split"] == "test" for f in fs),
                   "days_median": med([f["n_days"] for f in fs]), "bouts_median": med([f["n_bouts"] for f in fs]),
                   "span_years_median": med([f["span_years"] for f in fs]),
                   "active_fraction_median": med([f["active_fraction"] for f in fs]),
                   "post_index_active_fraction_median": med([f["post_index_active_fraction"] for f in fs]),
                   "early_active_fraction_median": early, "late_active_fraction_median": late,
                   "late_minus_early_median": med([f["late_minus_early"] for f in fs]),
                   "recency_days_median": med([f["time_since_last_days"] for f in fs]),
                   "max_gap_days_median": med([f["max_gap_days"] for f in fs]),
                   "hospital_fraction_median": med([f["hospital_fraction"] for f in fs]),
                   "margin_median": med([f["assignment_margin"] for f in fs]),
                   "trajectory": trajectory_mean(fs, "_trajectory"), "at_risk": trajectory_denominators(fs),
                   "hospital_trajectory": trajectory_mean(fs, "_hospital_trajectory"),
                   "gp_trajectory": trajectory_mean(fs, "_gp_trajectory"),
                   "pattern_mix": [{"id": p, "label": PATTERN_LABELS[p], "n": n} for p, n in count_patterns.most_common()],
                   "features": []}
        # Scale-free display against frozen train percentiles, not test ranks.
        refs = model.get("feature_reference", [])
        if len(idx) and refs:
            vals = []
            for q, feature in enumerate(model["feature_names"]):
                values = np.round(x[idx, names.index(feature)], 8)
                r = np.round(refs[q], 8)
                z = (np.searchsorted(r, values, "left") + np.searchsorted(r, values, "right")) / (2 * len(r))
                vals.append(float(np.median(z)))
            profile["features"] = vals
        profiles.append(profile)
    return profiles


def cluster(feats, used_sources, nbins, window, args):
    eligible = [f for f in feats if f["cluster_reason"] == "eligible"]
    excluded = Counter(f["cluster_reason"] for f in feats if f["cluster_reason"] != "eligible")
    include_source = bool(args.coverage) and set(used_sources) == {"hospital", "gp"} and not getattr(args, "partial_sources", [])
    scope = "\u6765\u6e90\u8986\u76d6\u5df2\u786e\u8ba4\u7684\u56fa\u5b9a\u7a97\u53e3" if args.coverage else "\u56fa\u5b9a\u7a97\u53e3\u53ef\u89c1\u65e5\u671f\u6a21\u5f0f；\u6765\u6e90\u8986\u76d6\u672a\u77e5"
    result = {"status": "insufficient_data", "scope": scope, "k": 0, "n": len(eligible), "silhouette": None, "stability_ari": None,
              "candidates": [], "profiles": [], "excluded": [{"reason": k, "n": n} for k, n in sorted(excluded.items())],
              "feature_labels": [], "bin_labels": [], "projection": {"points": [], "explained": []}, "notes": [],
              "source_availability": {s: s in used_sources for s in ("hospital", "gp")},
              "train_n": sum(f["split"] == "train" for f in eligible), "test_n": sum(f["split"] == "test" for f in eligible),
              "test_silhouette": None, "fit_partition": "train", "applied_frozen": bool(args.apply_model)}
    model = {"version": VERSION, "feature_schema": "aod_trajectory_v2", "fit_scope": scope, "seed": args.seed,
             "source_mix_feature_used": include_source, "status": "insufficient_data", "k": 0,
             "fit_partition": "train", "split_method": "blake2b(seed:eid) modulo4: 0=test, otherwise train",
             "settings": {"winsor_quantiles": [.01, .99], "min_silhouette": .15, "min_stability_ari": .60,
                          "fit_sample_limit": args.cluster_sample, "minimum_cluster_n": args.min_cluster_n,
                          "stability_repeats": args.stability_repeats}}
    axis_offset = -window / DAYS_YEAR if args.mode == "baseline" else 0.0
    axis_edges = [axis_offset + i * window / DAYS_YEAR / nbins for i in range(nbins + 1)]
    if args.mode == "baseline":
        axis_edges[-1] = 0.0
    result["bin_labels"] = [f"{axis_edges[i]:.2f}–{axis_edges[i+1]:.2f}\u5e74" for i in range(nbins)]
    result["notes"] = ["C1…CK\u662f\u540d\u4e49\u7c7b\u522b，\u4e0d\u662f\u4e25\u91cd\u7a0b\u5ea6\u5206\u671f\u6216\u8f6c\u79fb\u987a\u5e8f；K\u4ec5\u7531\u8bad\u7ec3\u96c6\u786e\u5b9a。",
                       "\u805a\u7c7b\u6807\u7b7e\u4e0d\u4f7f\u7528\u672a\u6765\u6b7b\u4ea1\u6216\u540e\u7eed\u8bb0\u5f55。\u6d4b\u8bd5\u96c6\u4ec5\u7528\u4e8e\u5206\u79bb\u5ea6\u4e0e\u4e0b\u6e38\u9884\u540e\u9a8c\u8bc1。",
                       ("\u9996\u4e2a\u5b63\u5ea6\u56e0\u9996\u6b21\u65e5\u671f\u5bf9\u9f50\u5fc5\u7136\u6709\u8bb0\u5f55；\u6301\u7eed\u6027\u548c\u65e9\u665a\u8d8b\u52bf\u7279\u5f81\u5254\u9664\u4e86\u8fd9\u4e2a\u673a\u68b0\u6027\u5cf0\u503c。"
                        if args.mode == "history" else
                        "\u65f6\u95f4\u8f74\u4ee5\u57fa\u7ebf\u8bc4\u4f30\u4e3a0；\u6301\u7eed\u6027\u548c\u65e9\u665a\u8d8b\u52bf\u4f7f\u7528\u57fa\u7ebf\u524d\u7a97\u53e3\u7684\u5168\u90e8\u5b63\u5ea6。") ]
    x, names, labels = cluster_features(eligible, window, include_source)
    if not len(eligible):
        return result, model
    if not np.all(np.isfinite(x)):
        raise ValueError("Eligible clustering features must be finite")
    rng = np.random.default_rng(args.seed)
    train_idx = np.array([i for i, f in enumerate(eligible) if f["split"] == "train"], dtype=int)
    test_idx = np.array([i for i, f in enumerate(eligible) if f["split"] == "test"], dtype=int)
    saved = getattr(args, "apply_model_payload", None)
    if saved:
        expected = {"trait": args.Y, "mode": args.mode, "trajectory_source": args.trajectory_source,
                    "window_days": window, "episode_gap": args.episode_gap, "feature_schema": "aod_trajectory_v2"}
        for key, value in expected.items():
            if saved.get(key) != value:
                raise ValueError(f"Frozen model incompatible for {key}: expected {value}, found {saved.get(key)}")
        if saved.get("source_mix_feature_used") != include_source or saved.get("coverage_supplied") != bool(args.coverage):
            raise ValueError("Frozen model requires the same coverage/source feature definition")
        current_levels = {k: args.source_levels.get(k, -1) for k in used_sources}
        if any(saved.get("source_levels", {}).get(k) != v for k, v in current_levels.items()):
            raise ValueError("Frozen model requires the same source record granularity")
        if int(saved.get("k", 0)) < 1 or not saved.get("centers"):
            raise ValueError("Frozen model has no fitted partition")
        model = dict(saved)
        model.update(applied_frozen=True, application_n=len(eligible))
        keep_names = model["feature_names"]
        if any(v not in names for v in keep_names):
            raise ValueError("Frozen model contains unknown features")
        keep = np.array([names.index(n) for n in keep_names], dtype=int)
        lo, hi, center, scale = [np.asarray(model[n], dtype=float) for n in ("winsor_low", "winsor_high", "center", "scale")]
        z = (np.clip(x[:, keep], lo, hi) - center) / scale if len(keep) else np.zeros((len(x), 0))
        c = np.asarray(model["centers"], dtype=float)
        if c.shape != (model["k"], len(keep)) or (len(scale) and np.any(scale <= 0)):
            raise ValueError("Invalid frozen cluster model dimensions/scales")
        lab, _ = assign(z, c)
        result.update(status="applied_frozen", k=model["k"], silhouette=model.get("silhouette"), stability_ari=model.get("stability_ari"))
    else:
        minimum = max(2 * args.min_cluster_n, 40)
        if len(train_idx) < minimum:
            result["notes"].append(f"\u8bad\u7ec3\u96c6\u4e0d\u8db3{minimum}\u4eba；\u4fdd\u7559\u8bb0\u5f55\u6a21\u5f0f，\u4e0d\u5f3a\u5236\u5206\u7fa4。")
            for f in eligible:
                f["cluster_reason"] = "too_few_training_cases"
            return result, model
        fit_idx = np.sort(rng.choice(train_idx, min(len(train_idx), args.cluster_sample), replace=False))
        rawfit = x[fit_idx]
        lo0, hi0 = np.quantile(rawfit, [.01, .99], axis=0)
        clipped = np.clip(rawfit, lo0, hi0)
        center0 = np.median(clipped, axis=0)
        scale0 = np.quantile(clipped, .75, axis=0) - np.quantile(clipped, .25, axis=0)
        scale0 = np.where(scale0 > 1e-8, scale0, np.std(clipped, axis=0))
        keep = np.flatnonzero(scale0 > 1e-8)
        lo, hi, center, scale = lo0[keep], hi0[keep], center0[keep], scale0[keep]
        z = (np.clip(x[:, keep], lo, hi) - center) / scale if len(keep) else np.zeros((len(x), 0))
        fit = z[fit_idx]
        model.update(feature_names=[names[j] for j in keep], winsor_low=lo, winsor_high=hi, center=center, scale=scale,
                     fit_n=len(fit), train_n=len(train_idx), test_n=len(test_idx), n=len(x),
                     feature_reference=[np.sort(rawfit[:, j]).tolist() for j in keep],
                     clipped_values=int(((x[:, keep] < lo) | (x[:, keep] > hi)).sum()))
        accepted = []
        ev_idx = np.sort(rng.choice(train_idx, min(800, len(train_idx)), replace=False))
        ev = z[ev_idx]
        maxk = min(8, len(np.unique(fit, axis=0)), len(train_idx) // args.min_cluster_n)
        candidates = list(range(2, min(6, maxk) + 1)) if args.k == "auto" else ([int(args.k)] if int(args.k) <= maxk else [])
        for k in candidates if len(keep) else []:
            c = kmeans(fit, k, rng)
            if c is None:
                continue
            train_lab, _ = assign(z[train_idx], c)
            ev_lab, _ = assign(ev, c)
            cnt = np.bincount(train_lab, minlength=k)
            sil = silhouette(ev, ev_lab)
            aris = []
            for _ in range(args.stability_repeats):
                ii = rng.choice(len(fit), max(k, int(.8 * len(fit))), replace=False)
                cc = kmeans(fit[ii], k, rng, starts=3)
                aris.append(adjusted_rand(ev_lab, assign(ev, cc)[0]) if cc is not None else 0.0)
            ari = float(np.mean(aris))
            ok = int(cnt.min()) >= max(args.min_cluster_n, int(math.ceil(.01 * len(train_idx)))) and sil is not None and sil >= .15 and ari >= .60
            result["candidates"].append({"k": k, "silhouette": sil, "stability_ari": ari, "min_n": int(cnt.min()), "accepted": ok})
            if ok:
                accepted.append((float(sil), k, ari, c))
        if accepted:
            top = max(v[0] for v in accepted)
            sil, k, ari, c = min((v for v in accepted if v[0] >= top-.02), key=lambda q: q[1])
            tlab = assign(z[train_idx], c)[0]
            order = sorted(range(k), key=lambda j: (float(np.median([eligible[train_idx[i]]["n_days"]
                                                                      for i in np.flatnonzero(tlab == j)])), j))
            c = c[order]
            result.update(status="ok", k=k, silhouette=sil, stability_ari=ari)
            model.update(status="ok", k=k, silhouette=sil, stability_ari=ari)
        else:
            k = 1
            c = np.mean(fit, axis=0)[None, :]
            status = "no_separation" if not len(keep) else "no_stable_partition"
            result.update(status=status, k=1)
            model.update(status=status, k=1)
            result["notes"].append("\u8bad\u7ec3\u96c6\u672a\u652f\u6301\u591a\u4e2a\u7a33\u5b9a\u5206\u578b；\u4fdd\u7559\u4e00\u4e2a\u672a\u5206\u578b\u7ec4。")
        model.update(centers=c, label_order=[f"C{j+1}" for j in range(k)], candidates=result["candidates"])
        lab, _ = assign(z, c)
        if len(keep):
            pmean = fit.mean(axis=0)
            _, singular, vectors = np.linalg.svd(fit-pmean, full_matrices=False)
            v = vectors[:min(2, len(keep))].copy()
            for axis in range(len(v)):
                if v[axis, np.argmax(np.abs(v[axis]))] < 0:
                    v[axis] *= -1
            exp = singular[:len(v)]**2 / max(float((singular**2).sum()), 1e-12) * 100
            model.update(pca_mean=pmean, pca_loadings=v, pca_explained=exp)
        else:
            model.update(pca_mean=[], pca_loadings=[], pca_explained=[])
    k = int(model["k"])
    if len(test_idx) >= 3 and k > 1:
        ii = np.sort(rng.choice(test_idx, min(800, len(test_idx)), replace=False))
        result["test_silhouette"] = silhouette(z[ii], lab[ii])
    if saved and len(z) >= 3 and k > 1:
        ii = np.sort(rng.choice(len(z), min(800, len(z)), replace=False))
        result["application_silhouette"] = silhouette(z[ii], lab[ii])
    for begin in range(0, len(z), 4096):
        zz = z[begin:begin+4096]
        ds = np.maximum(0, np.sum(zz*zz, axis=1)[:, None] + np.sum(c*c, axis=1)[None, :] - 2*zz@c.T)
        ordered = np.sort(ds, axis=1)
        for j, f in enumerate(eligible[begin:begin+4096]):
            f["cluster"] = model["label_order"][int(lab[begin+j])]
            f["cluster_reason"] = "assigned" if k > 1 else "no_stable_partition"
            f["assignment_margin"] = float((ordered[j, 1]-ordered[j, 0]) / max(ordered[j, 1], 1e-12)) if k > 1 else None
    result["feature_labels"] = [labels[names.index(n)] for n in model["feature_names"]]
    if args.mode == "baseline":
        result["feature_labels"] = ["\u6d3b\u8dc3\u5b63\u5ea6" if label == "\u540e\u7eed\u6d3b\u8dc3\u5b63\u5ea6" else label for label in result["feature_labels"]]
    result["profiles"] = cluster_profiles(eligible, lab, model, x, names, labels, window)
    result["low_margin_n"] = sum(f["assignment_margin"] is not None and f["assignment_margin"] < .1 for f in eligible)
    if len(model.get("pca_loadings", [])):
        xy = (z-np.asarray(model["pca_mean"])) @ np.asarray(model["pca_loadings"]).T
        if xy.shape[1] == 1:
            xy = np.column_stack([xy[:, 0], np.zeros(len(xy))])
        limits = np.quantile(xy, [.005, .995], axis=0)
        width = np.maximum(limits[1]-limits[0], 1e-9)
        points = Counter()
        for row, label in zip(xy, lab):
            bins = np.clip(((row-limits[0])/width*30).astype(int), 0, 29)
            points[int(bins[0]), int(bins[1]), int(label)] += 1
        result["projection"] = {"points": [{"x": float(limits[0, 0]+(bx+.5)/30*width[0]),
                                             "y": float(limits[0, 1]+(by+.5)/30*width[1]), "n": n,
                                             "cluster": model["label_order"][cl]} for (bx, by, cl), n in sorted(points.items())],
                                "explained": clean_json(model["pca_explained"])}
    return result, model


def aggregate(patients, feats, levels, selected, qc, clusters, score, model, meta, args, nbins, window):
    pp = list(patients.values())
    baseline = [p["baseline"] for p in pp]
    reference = reference_audit(meta)
    common = common_comparison([p["min1"] for p in pp], [p["min2"] for p in pp], [p["ref"] for p in pp])
    sources = [{"id": key, "label": label, "n": sum(p[key] is not None for p in pp)} for key, label in SOURCES.items()]
    pairs = []
    for i, a in enumerate(SOURCES):
        for b in list(SOURCES)[i + 1:]:
            pairs.append(dict(a=a, b=b, **pair_metrics([p[a] for p in pp], [p[b] for p in pp], baseline)))
    default = next(p for p in pairs if p["a"] == "icd10" and p["b"] == "ref")
    if default["both"]:
        headline = f"\u5171\u540c\u6709\u65e5\u671f\u8005\u4e2d，ICD10 \u4e0e ref \u540c\u65e5\u7387 {default['same_day_pct']:.1f}%；\u5dee\u503c\u53cd\u6620\u6765\u6e90\u4e0e\u5b9a\u4e49\u5dee\u5f02。"
    else:
        headline = "\u5f53\u524d\u8f93\u5165\u6ca1\u6709\u53ef\u6bd4\u8f83\u7684 ICD10–ref \u914d\u5bf9；\u53ef\u5207\u6362\u6765\u6e90，\u7f3a\u5931\u4e0d\u8ba1\u4f5c\u4e00\u81f4。"
    withdays = [f for f in feats if f["n_days"] is not None and (f["n_days"] > 0 or f["coverage_complete"])]
    reasons = {"quiet": "\u7a97\u5185\u672a\u89c1\u8bb0\u5f55", "single": "\u5355\u4e00\u65e5\u671f", "compact": "\u77ed\u671f\u96c6\u4e2d", "intermittent": "\u8de8\u671f\u95f4\u6b47", "persistent": "\u6301\u7eed\u88ab\u8bb0\u5f55"}
    timeline_groups = defaultdict(list)
    timeline_excluded = Counter()
    for f in withdays:
        if not f["horizon_complete"]:
            timeline_excluded["death_loss_or_cutoff_before_window_end"] += 1
            continue
        if args.coverage and not f["coverage_complete"]:
            timeline_excluded["incomplete_source_coverage"] += 1
            continue
        group = "quiet" if f["n_days"] == 0 else ("single" if f["n_days"] == 1 else ("persistent" if f["active_fraction"] >= .5 else ("intermittent" if f["span_years"] >= 1 else "compact")))
        timeline_groups[group].append(f)
    axis_offset = -window / DAYS_YEAR if args.mode == "baseline" else 0.0
    labels = [f"{axis_offset + i * window / DAYS_YEAR / nbins:.1f}" for i in range(nbins)]
    tl = [{"id": key, "label": reasons[key], "n": len(fs), "values": trajectory_mean(fs, "_trajectory")} for key, fs in timeline_groups.items()]
    timeline_n = sum(len(fs) for fs in timeline_groups.values())
    notes = ["\u4e0b\u65b9\u65f6\u95f4\u5f62\u6001\u662f\u89e3\u91ca\u7528\u89c4\u5219\u5206\u7ec4；\u7b2c\u4e09\u90e8\u5206\u624d\u662f\u6570\u636e\u5b66\u4e60\u5f97\u5230\u7684\u65e0\u76d1\u7763 cluster。",
             "\u540c\u4e00\u4eba\u540c\u4e00\u5929\u7684\u591a\u4ee3\u7801\u4e0eGP\u53cc\u7f16\u7801\u53ea\u8ba1\u4e00\u4e2a\u8bb0\u5f55\u65e5；\u76f8\u90bb\u95f4\u9694\u4e0d\u8d85\u8fc7\u9608\u503c\u7684\u8bb0\u5f55\u65e5\u5408\u5e76\u4e3a\u8bb0\u5f55\u7c07，\u5e76\u975e\u4e34\u5e8a\u4f4f\u9662\u6b21\u6570\u6216\u65b0\u53d1\u4e8b\u4ef6。",
             "\u5408\u6cd5\u7684\u6781\u65e9\u65e5\u671f\u53ea\u6807\u8bb0\u5e76\u505a\u654f\u611f\u6027\u6bd4\u8f83；\u6a21\u578b\u7279\u5f81\u6309\u8bad\u7ec3\u6837\u672c1%/99%\u7f29\u5c3e，\u539f\u59cb\u65e5\u671f、\u4eba\u6570\u548c\u7edf\u8ba1\u503c\u4e0d\u88ab\u5220\u9664。",
             "Y.adu\u8868\u793a\u5e74\u9f84\u76f8\u5173liability；Y.aod\u662f\u75c5\u4f8b\u5185\u90e8\u63a2\u7d22\u6027\u7684\u53ef\u89c2\u5bdf\u8bb0\u5f55\u8d1f\u62c5，\u4e0d\u66ff\u6362ADuLT，\u4e5f\u4e0d\u80fd\u8bc1\u660e\u6cbb\u6108\u6216\u4e25\u91cd\u7a0b\u5ea6。"]
    if not args.coverage:
        notes.insert(0, "\u672a\u63d0\u4f9b\u6765\u6e90\u89c2\u5bdf\u533a\u95f4：\u5c55\u793a\u53ef\u89c1\u65e5\u671f\u6a21\u5f0f，\u5e74\u5ea6\u4e8b\u4ef6\u7387\u4e0e Y.aod \u7559\u7a7a；\u4e0d\u80fd\u628a\u672a\u8bb0\u5f55\u65f6\u95f4\u89e3\u91ca\u4e3a\u5065\u5eb7\u65f6\u95f4。")
    else:
        notes.insert(0, "\u56fa\u5b9a\u7a97\u53e3\u8f68\u8ff9\u4ec5\u7eb3\u5165\u6240\u6709\u5df2\u7528\u6765\u6e90\u8986\u76d6\u5b8c\u6574\u8005；\u8986\u76d6\u4e0d\u5b8c\u6574\u8005\u4fdd\u7559\u53ef\u89c1\u8bb0\u5f55\u7edf\u8ba1，\u4f46\u4e0d\u5c06\u672a\u89c2\u5bdf\u65f6\u95f4\u7ed8\u6210\u96f6\u6d3b\u52a8。")
    if any(level == 1 for level in levels.values()):
        notes.insert(0, "\u533b\u9662\u6c47\u603bAOD\u5305\u542b\u5404\u8bca\u65ad\u4ee3\u7801\u7684\u9996\u6b21\u65e5\u671f，\u9057\u6f0f\u540c\u4e00\u4ee3\u7801\u540e\u7eed\u4f4f\u9662；\u5bf9\u5e94\u5206\u578b\u662f\u4fdd\u7559\u65e5\u671f\u6a21\u5f0f。")
    if args.mode == "baseline":
        notes.append("baseline\u7279\u5f81\u548c\u805a\u7c7b\u4ec5\u4f7f\u7528\u4e8b\u4ef6\u65e5\u671f\u4e0d\u665a\u4e8e\u9996\u6b21\u8bc4\u4f30\u7684\u8bb0\u5f55；\u7b2c\u4e00\u90e8\u5206\u65e5\u671f\u5ba1\u8ba1\u4fdd\u7559\u5b8c\u6574\u6709\u6548FOD。\u7f3a\u5c11\u8bb0\u5f55\u5f55\u5165\u65f6\u95f4\u65f6，\u4ecd\u4e0d\u80fd\u4fdd\u8bc1\u56de\u6eaf\u8bb0\u5f55\u5728\u57fa\u7ebf\u5f53\u65f6\u5df2\u7ecf\u53ef\u83b7\u5f97。")
    warnings = list(reference["notes"])
    if meta.get("reference_scope_warning"):
        warnings.append(str(meta["reference_scope_warning"]))
    if meta.get("official_maps_available") is False:
        warnings.append("\u672a\u63d0\u4f9b\u516c\u5f00\u6620\u5c04\u8868：\u4fdd\u5b58\u7684FOD\u4ecd\u53ef\u6bd4\u8f83，min2\u4ec5\u662f\u5f53\u524d\u53ef\u91cd\u5efa\u90e8\u5206\u7684\u5019\u9009\u65e5\u671f。")
    if not args.coverage:
        warnings.append("\u6765\u6e90\u89c2\u5bdf\u671f\u672a\u77e5：\u805a\u7c7b\u4e3a\u63a2\u7d22\u6027\u65e5\u671f\u6a21\u5f0f；Y.aod\u4e0e\u8054\u5408\u4e8b\u4ef6\u7387\u672a\u8ba1\u7b97。")
    elif timeline_excluded["incomplete_source_coverage"]:
        warnings.append(f"\u5df2\u63d0\u4f9b\u89c2\u5bdf\u671f，\u4f46\u6709 {timeline_excluded['incomplete_source_coverage']} \u4eba\u7684\u56fa\u5b9a\u7a97\u53e3\u8986\u76d6\u4e0d\u5b8c\u6574；\u672a\u7eb3\u5165\u5b8c\u6574\u7a97\u53e3\u8f68\u8ff9，\u672a\u89c2\u5bdf\u65f6\u6bb5\u4e0d\u8865\u96f6。")
    if any(level == 1 for level in levels.values()):
        warnings.append("\u81f3\u5c11\u4e00\u4e2a\u6765\u6e90\u53ea\u6709\u6bcf\u4ee3\u7801\u9996\u6b21\u65e5\u671f，\u4e0d\u80fd\u89e3\u91ca\u6210\u5b8c\u6574\u91cd\u590d\u5c31\u8bca\u5386\u53f2。")
    if getattr(args, "partial_sources", []):
        warnings.append("\u539f\u59cb\u8868\u672a\u5bfc\u51fa\u75be\u75c5\u5b9a\u4e49\u9700\u8981\u7684\u5168\u90e8\u7ec6\u6765\u6e90\u4ee3\u7801\u5217：" + ", ".join(args.partial_sources) + "。\u4fdd\u7559\u8be5\u63d0\u53d6\u8303\u56f4\u7684\u65e5\u671f\u6a21\u5f0f，\u7efc\u5408\u8d1f\u62c5\u5206\u6570\u4e0d\u8ba1\u7b97。")
    eligible_score = [f["aod"] for f in feats if f["aod"] is not None]
    suspect = [f for f in feats if f["early_isolate"]]
    qc_labels = {"after_death": "\u6b7b\u4ea1\u4e4b\u540e", "before_birth": "\u51fa\u751f\u4e4b\u524d", "after_as_of": "\u665a\u4e8e\u622a\u70b9", "after_analysis_cutoff": "\u5206\u6790\u7a97\u5916", "unparseable": "\u65e0\u6cd5\u89e3\u6790", "unknown": "\u65e5\u671f\u672a\u77e5", "hospital_end_date_proxy": "\u4f4f\u9662\u7ed3\u675f\u65e5\u4ee3\u7406"}
    return {"schema_version": 1, "trait": args.Y, "generated_at": datetime.now(timezone.utc).isoformat(), "demo": truth(meta.get("demo")), "mode": args.mode,
            "cohort": {"cases": len(feats), "with_clinical_dates": sum(f["n_days"] > 0 for f in withdays), "clustered": sum(f["cluster"] is not None for f in feats), "score_eligible": len(eligible_score)},
            "availability": [{"source": s, "label": "\u4f4f\u9662" if s == "hospital" else "GP", "mode": LEVEL_LABEL.get(level, "unavailable"),
                              "n_people": sum(bool(v) for v in selected[s].values()), "n_days": sum(len(v) for v in selected[s].values())} for s, level in levels.items()],
            "input_date_ranges": meta.get("input_date_ranges", []),
            "warnings": warnings, "consistency": {"sources": sources, "default_pair": ["icd10", "ref"], "pairs": pairs, "headline": headline,
                       "reference": reference, "common_sample": common,
                       "audit_scope": "\u65e5\u671f\u6bd4\u8f83\u4f7f\u7528\u51fa\u751f、\u6b7b\u4ea1、\u5931\u8bbf\u53ca\u5206\u6790\u622a\u6b62\u65e5\u671f\u7ea6\u675f\u5185\u7684\u5168\u90e8\u6709\u6548FOD；baseline\u6a21\u5f0f\u53ea\u622a\u65ad\u8d1f\u62c5\u4e0e\u805a\u7c7b\u7279\u5f81，\u4e0d\u622a\u65ad\u65e5\u671f\u5ba1\u8ba1。"},
            "burden": {"window_label": (f"\u8bc4\u4f30\u524d {args.lookback_years:g} \u5e74" if args.mode == "baseline" else f"\u9996\u6b21\u53ef\u89c1\u4e34\u5e8a\u65e5\u671f\u8d77 {args.landmark_years:g} \u5e74"),
                       "scope": "+".join(s for s, l in levels.items() if l >= 1), "headline": "\u628a\u8d77\u75c5\u65e9\u665a\u4e0e\u540e\u7eed\u8bb0\u5f55\u5bc6\u5ea6、\u95f4\u9694\u548c\u6301\u7eed\u6027\u5206\u5f00\u770b。",
                       "feature_definitions": [{"key": k, "label": l, "description": d} for k, l, d in [
                         ("n_days", "\u8bb0\u5f55\u65e5", "\u4e0d\u540c\u6709\u6548\u65e5\u671f\u6570；\u540c\u65e5\u4e0d\u540c\u4ee3\u7801\u4e0d\u91cd\u590d\u8ba1\u6570"), ("n_bouts", "\u8bb0\u5f55\u7c07", f"\u76f8\u90bb\u8bb0\u5f55\u95f4\u9694≤{args.episode_gap}\u5929\u4e32\u8054；\u4e0d\u7b49\u540c\u4f4f\u9662\u6b21\u6570"),
                         ("active_fraction", "\u6d3b\u8dc3\u65f6\u95f4\u6bb5", "\u56fa\u5b9a\u7a97\u53e3\u5185\u6709\u8bb0\u5f55\u7684\u5b63\u5ea6\u6bd4\u4f8b；\u65e0\u8986\u76d6\u4fe1\u606f\u65f6\u4ec5\u4f5c\u8bb0\u5f55\u5bc6\u5ea6"), ("median_gap_days", "\u76f8\u90bb\u95f4\u9694", "\u4e2a\u4eba\u53bb\u91cd\u65e5\u671f\u95f4\u9694\u7684\u4e2d\u4f4d\u6570；\u5355\u6b21\u8005\u7f3a\u5931"),
                         ("event_rate", "\u5e74\u5ea6\u8bb0\u5f55\u7387", "\u4ec5\u5728\u5b8c\u6574\u6765\u6e90\u89c2\u5bdf\u671f\u548c\u9010\u6b21\u8bb0\u5f55\u6761\u4ef6\u6210\u7acb\u65f6\u8ba1\u7b97"), ("aod", "Y.aod", "\u4e09\u4e2a\u75c5\u4f8b\u5185\u767e\u5206\u4f4d\u7684\u7b49\u6743\u5747\u503c：\u8bb0\u5f55\u7387、\u6d3b\u8dc3\u5b63\u5ea6、\u9996\u672b\u8de8\u5ea6；\u63a2\u7d22\u6027")]],
                       "distributions": {"days": hist([f["n_days"] for f in withdays], [0, 1, 2, 5, 10, 20, 50], ["0（\u5df2\u786e\u8ba4\u89c2\u5bdf）", "1", "2", "3–5", "6–10", "11–20", "21–50", ">50"]),
                           "bouts": hist([f["n_bouts"] for f in withdays], [0, 1, 2, 5, 10, 20], ["0", "1", "2", "3–5", "6–10", "11–20", ">20"]),
                           "gap": hist([f["median_gap_days"] for f in withdays], [30, 90, 365, 730, 1826], ["≤30\u5929", "31–90\u5929", "91\u5929–1\u5e74", "1–2\u5e74", "2–5\u5e74", ">5\u5e74"]),
                           "score": hist(eligible_score, [20, 40, 60, 80], ["0–20", "20–40", "40–60", "60–80", "80–100"])},
                       "timeline": {"labels": labels, "rows": tl, "eligible_n": timeline_n, "candidate_n": len(withdays),
                                    "coverage_confirmed": bool(args.coverage),
                                    "scope": "\u5b8c\u6574\u6765\u6e90\u8986\u76d6\u7684\u56fa\u5b9a\u7a97\u53e3" if args.coverage else "\u53ef\u89c1\u65e5\u671f\u6a21\u5f0f；\u6765\u6e90\u8986\u76d6\u671f\u672a\u786e\u8ba4",
                                    "excluded": [{"reason": key, "n": value} for key, value in sorted(timeline_excluded.items()) if value]},
                       "score": {"n": len(eligible_score), "method": score["method"], "scope": score.get("scope", "unavailable"), "status": score["status"]},
                       "qc": {"counts": [{"label": qc_labels.get(k, k), "n": n} for k, n in sorted(qc.items())], "sensitivity": [
                           {"label": "\u5b64\u7acb\u6781\u65e9\u65e5\u671f\u8005（\u4e3b\u5206\u6790\u4fdd\u7559）", "value": len(suspect)},
                           {"label": "\u6682\u4e0d\u91c7\u7528\u5b64\u7acb\u65e5\u671f\u540e，\u7a97\u53e3\u8bb0\u5f55\u65e5\u6570\u4e2d\u4f4d\u53d8\u5316", "value": med([f["n_days_without_isolate"] - f["n_days"] for f in suspect])},
                           {"label": "\u4ec5\u5728\u5efa\u6a21\u7279\u5f81\u4e2d\u7f29\u5c3e\u7684\u503c", "value": model.get("clipped_values", 0)}]}, "notes": notes},
            "clusters": clusters, "methods": [{"title": "\u89c2\u5bdf\u5206\u6bcd", "text": "\u75c5\u4f8b\u8bc1\u636e\u5e76\u96c6\u4e0d\u7b49\u4e8e\u5171\u540c\u6709\u65e5\u671f\u4eba\u6570。\u6b7b\u4ea1、\u5931\u8bbf\u6216\u622a\u70b9\u524d\u672a\u5b8c\u6210\u7a97\u53e3\u8005\u4e0d\u8865\u96f6，\u5206\u7fa4\u6392\u9664\u539f\u56e0\u53ef\u67e5。"},
                         {"title": "\u53ef\u91cd\u590d\u805a\u7c7b", "text": "NumPy\u591a\u8d77\u70b9k-means；\u8bad\u7ec3\u6837\u672c1%/99%\u7f29\u5c3e\u4e0e\u4e2d\u4f4d\u6570/IQR\u6807\u51c6\u5316；K=2–6\u6216\u6307\u5b9a2–8；\u62bd\u6837silhouette、80%\u91cd\u62bd\u6837ARI\u53ca\u6700\u5c0f\u7c07\u4eba\u6570\u5171\u540c\u7b5b\u9009，\u5141\u8bb8\u4e0d\u5206\u7fa4。"},
                         {"title": "\u5206\u6790\u7528\u9014", "text": "\u6a21\u578b\u5728\u672c\u6b21\u75c5\u4f8b\u96c6\u62df\u5408，\u4ec5\u4f5c\u63cf\u8ff0\u6027\u8868\u578b\u63a2\u7d22。\u5916\u90e8\u9a8c\u8bc1、GWAS\u548c\u9884\u6d4b\u9700\u8981\u72ec\u7acb\u8bad\u7ec3/\u9a8c\u8bc1\u8bbe\u8ba1；\u4e0d\u80fd\u7528\u672a\u6765\u8f68\u8ff9\u4f5c\u4e3abaseline\u9884\u6d4b\u53d8\u91cf。"}]}



def add_prognosis(feats, patients, selected, coverage, used_sources, levels, args):
    """Endpoints occur strictly after the feature window and never enter clustering.

    Death is all-cause follow-up to the common administrative end. The disease
    endpoint is another *record episode*, not clinically adjudicated recurrence.
    Death ends its cause-specific follow-up; no cumulative-incidence claim is made.
    """
    for f in feats:
        p = patients[f["eid"]]
        end = f["_end"]
        f["feature_mode"] = args.mode
        f.update(death_followup_years=None, death_event=None, future_record_time_years=None,
                 future_record_event=None, future_record_competing_death=None,
                 later_record_days=None, later_record_audit_followup_years=None,
                 prognosis_reason="incomplete_feature_window", record_prognosis_reason="unavailable")
        if end is None or not f["horizon_complete"]:
            continue
        if p["clinical_end"] <= end:
            f["prognosis_reason"] = "no_positive_followup"
            continue
        # This audit may include pre-recruitment windows and is purely retrospective.
        ds = sorted(set().union(*(set(selected[src].get(f["eid"], ())) for src in used_sources)))
        future = [d for d in ds if end < d <= p["clinical_end"]]
        f["later_record_days"] = len(future)
        f["later_record_audit_followup_years"] = (p["clinical_end"]-end) / DAYS_YEAR
        # UKB participants must survive to recruitment. Starting a mortality
        # clock earlier would introduce guaranteed survival/left truncation.
        if p["baseline"] is None:
            f["prognosis_reason"] = f["record_prognosis_reason"] = "missing_recruitment_date"
            continue
        if end < p["baseline"]:
            f["prognosis_reason"] = f["record_prognosis_reason"] = "landmark_before_recruitment"
            continue
        # A supplied coverage gap invalidates the historical phenotype, whereas
        # absent coverage is explicitly labelled exploratory in the report.
        if args.coverage and not f["coverage_complete"]:
            f["prognosis_reason"] = "incomplete_feature_coverage"
            continue
        f["prognosis_reason"] = "eligible"
        f["death_followup_years"] = (p["clinical_end"]-end) / DAYS_YEAR
        f["death_event"] = int(p["death"] is not None and end < p["death"] <= p["clinical_end"])
        if not f["recurrence_capable"]:
            f["record_prognosis_reason"] = "first_per_code_or_incomplete_record_source"
            continue
        if not f["coverage_complete"]:
            f["record_prognosis_reason"] = "unverified_feature_coverage"
            continue
        # Stop at the first coverage break, never treat an unobserved gap as zero.
        coverage_ends = []
        for src in used_sources:
            intervals = coverage.get((f["eid"], src), [])
            containing = [b for a, b in intervals if a <= end+1 <= b]
            if not containing:
                coverage_ends = []
                break
            coverage_ends.append(max(containing))
        if len(coverage_ends) != len(used_sources):
            f["record_prognosis_reason"] = "no_continuous_post_landmark_coverage"
            continue
        observed_end = min([p["clinical_end"]] + coverage_ends)
        if observed_end <= end:
            f["record_prognosis_reason"] = "no_positive_record_followup"
            continue
        episode_starts = [b[0] for b in bouts([d for d in ds if d <= observed_end], args.episode_gap)]
        next_dates = [d for d in episode_starts if end < d <= observed_end]
        nxt = next_dates[0] if next_dates else None
        endpoint = nxt if nxt is not None else observed_end
        f["future_record_time_years"] = (endpoint-end) / DAYS_YEAR
        f["future_record_event"] = int(nxt is not None)
        f["future_record_competing_death"] = int(nxt is None and p["death"] is not None and p["death"] == observed_end)
        f["record_prognosis_reason"] = "eligible_cause_specific_record_episode"


def distribution(values):
    vals = [v for v in values if v is not None and math.isfinite(v)]
    return {"n": len(vals), "median": med(vals), "q25": float(np.quantile(vals, .25)) if vals else None,
            "q75": float(np.quantile(vals, .75)) if vals else None}


def rank_correlation(feats, key):
    pairs = [(f["adu"], f.get(key)) for f in feats if f["adu"] is not None and f.get(key) is not None]
    if len(pairs) < 3:
        return {"feature": key, "n": len(pairs), "rho": None}
    x, y = np.asarray(pairs, dtype=float).T
    rx, ry = percentile_values(x), percentile_values(y)
    rho = float(np.corrcoef(rx, ry)[0, 1]) if np.std(rx) > 0 and np.std(ry) > 0 else None
    return {"feature": key, "n": len(pairs), "rho": rho}


def augment_summary(summary, feats, patients, levels, used_sources, model, meta, args, window):
    eligible = [f for f in feats if f["pattern"] is not None]
    pattern_rows = []
    for key, label in PATTERN_LABELS.items():
        fs = [f for f in eligible if f["pattern"] == key]
        if not fs:
            continue
        later = [f for f in fs if f.get("later_record_days") is not None]
        pattern_rows.append({"id": key, "label": label, "n": len(fs), "pct": percent(len(fs), len(eligible)),
                             "days_median": med([f["n_days"] for f in fs]), "bouts_median": med([f["n_bouts"] for f in fs]),
                             "recency_days_median": med([f["time_since_last_days"] for f in fs]),
                             "trajectory": trajectory_mean(fs, "_trajectory"), "at_risk": trajectory_denominators(fs),
                             "later_observable_n": len(later), "later_with_records_n": sum(f["later_record_days"] > 0 for f in later),
                             "later_audit_followup_years_median": med([f["later_record_audit_followup_years"] for f in later])})
    definitions = [
        {"id": "quiet_window", "text": "\u65e2\u5f80\u6709\u8be5\u6765\u6e90\u8bb0\u5f55；\u5f53\u524d\u5b8c\u6574\u89c2\u5bdf\u7a97\u53e3\u5185\u65e0\u8be5\u75c5\u8bb0\u5f55。"},
        {"id": "single_bout_quiet", "text": "\u7a97\u53e3\u5185\u4ec5\u4e00\u6bb5\u8bb0\u5f55，\u672b\u6b21\u8bb0\u5f55\u8ddd\u7a97\u672b\u81f3\u5c11min(1\u5e74,\u7a97\u53e3\u4e00\u534a)，\u6700\u4f4e90\u5929；\u662f\u56fa\u5b9a\u7a97\u5185\u672a\u518d\u8bb0\u5f55，\u4e0d\u4ee3\u8868\u6c38\u4e45\u6cbb\u6108。"},
        {"id": "increasing", "text": ("\u5254\u9664\u9996\u4e2a\u5bf9\u9f50\u5b63\u5ea6\u540e，" if args.mode == "history" else "\u4f7f\u7528\u57fa\u7ebf\u524d\u5168\u90e8\u5b63\u5ea6，") +
                                      "\u540e\u534a\u7a97\u6d3b\u8dc3\u5b63\u5ea6\u6bd4\u4f8b\u6bd4\u524d\u534a\u7a97\u9ad8\u81f3\u5c110.25，\u4e14\u540e\u534a\u7a97\u81f3\u5c110.5。"},
        {"id": "persistent", "text": "\u81f3\u5c11\u4e24\u4e2a\u8bb0\u5f55\u65e5；\u540e\u7eed\u5b63\u5ea6\u81f3\u5c11\u534a\u6570\u53ca\u540e\u534a\u7a97\u81f3\u5c11\u534a\u6570\u6709\u8bb0\u5f55，\u5141\u8bb8\u5c5e\u4e8e\u540c\u4e00\u6301\u7eed\u8bb0\u5f55\u6bb5。"},
        {"id": "early_limited", "text": "\u4e0d\u6ee1\u8db3\u524d\u8ff0\u89c4\u5219，\u4e14\u672b\u6b21\u8bb0\u5f55\u540e\u5b58\u5728\u81f3\u5c11\u89c4\u5b9a\u957f\u5ea6\u7684\u65e0\u8bb0\u5f55\u5c3e\u90e8。"},
        {"id": "intermittent", "text": "\u5269\u4f59\u5177\u6709\u81f3\u5c11\u4e24\u6bb5\u8bb0\u5f55\u8005；\u8bb0\u5f55\u540e\u518d\u6b21\u51fa\u73b0。"},
        {"id": "recent_limited", "text": "\u5176\u4f59\u4ec5\u4e00\u6bb5、\u9760\u8fd1\u7a97\u53e3\u672b\u7aef\u7684\u8bb0\u5f55。"}]
    sensitivity = []
    for gap in (7, 30, 90):
        counts = Counter(pattern_for(f, gap) for f in eligible)
        sensitivity.append({"episode_gap_days": gap, "n": len(eligible),
                            "median_bouts": med([len(bouts(f["_days"], gap)) for f in eligible]),
                            "changed_from_primary_n": sum(pattern_for(f, gap) != f["pattern"] for f in eligible),
                            "patterns": [{"id": k, "label": PATTERN_LABELS[k], "n": n} for k, n in sorted(counts.items()) if k]})
    summary["settings"] = {"date_follow_end": args.followup_end, "audit_as_of": args.resolved_asof,
                           "feature_cutoff": iso(min(ordinal(args.followup_end), ordinal(args.resolved_asof))),
                           "trajectory_source": args.trajectory_source, "window_days": window,
                           "feature_mode": args.mode, "episode_gap_days": args.episode_gap}
    summary["cohort"].update(analytic_cases=sum(f["analysis_case"] for f in feats), pattern_eligible=len(eligible),
                             feature_window_incomplete=sum(f["cluster_reason"] in {"death_before_window_end", "loss_before_window_end", "administrative_cutoff_before_window_end"} for f in feats),
                             train_n=sum(f["split"] == "train" for f in eligible), test_n=sum(f["split"] == "test" for f in eligible))
    summary["patterns"] = {"n": len(eligible), "rows": pattern_rows, "definitions": definitions, "gap_sensitivity": sensitivity,
                           "evidence": [{"label": k, "n": n} for k, n in Counter(f["pattern_evidence"] for f in eligible).items()],
                           "window_note": "\u6240\u6709\u6a21\u5f0f\u53ea\u63cf\u8ff0\u6307\u5b9a\u7a97\u53e3；\u540e\u7eed\u8bb0\u5f55\u7528\u4e8e\u72ec\u7acb\u5ba1\u8ba1，\u4e0d\u56de\u5199\u5206\u7c7b。",
                           "later_record_note": "later\u5217\u4f7f\u7528\u5206\u7c7b\u7a97\u53e3\u4e4b\u540e、\u884c\u653f\u622a\u6b62\u4e4b\u524d\u7684\u8bb0\u5f55，\u4ec5\u4f5c\u4e8b\u540e\u5ba1\u8ba1；\u968f\u8bbf\u957f\u77ed\u4e0d\u540c，\u4e0d\u80fd\u5f53\u4f5c\u5df2\u9a8c\u8bc1\u7684\u590d\u53d1\u7387。",
                           "index_bin_note": ("history\u7b2c\u4e00\u5b63\u5ea6\u542b\u5bf9\u9f50\u9996\u65e5；\u6a21\u5f0f\u4e0e\u805a\u7c7b\u7684\u6301\u7eed\u6027/\u65e9\u665a\u8d8b\u52bf\u7279\u5f81\u6392\u9664\u6b64\u5b63\u5ea6。"
                                              if args.mode == "history" else
                                              "baseline\u65f6\u95f4\u8f74\u4ee5\u8bc4\u4f30\u65e5\u4e3a0；\u6a21\u5f0f\u4e0e\u805a\u7c7b\u4f7f\u7528\u57fa\u7ebf\u524d\u7a97\u53e3\u7684\u5168\u90e8\u5b63\u5ea6，\u4e0d\u6392\u9664\u9996\u5b63\u5ea6。")}
    # Replace the old generic timeline with the explicit nominal patterns.
    summary["burden"]["timeline"].update(rows=[{"id": x["id"], "label": x["label"], "n": x["n"],
                                                "values": x["trajectory"], "at_risk": x["at_risk"]} for x in pattern_rows],
                                            eligible_n=len(eligible), index_bin_excluded_from_features=args.mode == "history")
    summary["burden"]["scope"] = "+".join(used_sources)
    summary["burden"]["observed_timeline"] = {"values": trajectory_mean(feats, "_trajectory"), "at_risk": trajectory_denominators(feats),
                                                "denominator_type": "complete_source_covered_bins" if args.coverage else "complete_calendar_bins_coverage_unknown",
                                                "notes": "\u672a\u5230\u8fbe\u8be5\u65f6\u95f4\u6bb5\u6216\u5df2\u77e5\u89c2\u5bdf\u4e2d\u65ad\u8005，\u4e0d\u8fdb\u5165\u8be5\u6bb5\u5206\u6bcd；\u672a\u77e5\u6765\u6e90\u8986\u76d6\u4ecd\u6807\u8bb0\u4e3a\u63a2\u7d22\u6027。"}
    summary["burden"]["feature_definitions"] += [
        {"key": "post_index_active_fraction", "label": "\u540e\u7eed\u6d3b\u8dc3\u5b63\u5ea6" if args.mode == "history" else "\u6d3b\u8dc3\u5b63\u5ea6",
         "description": "history\u5254\u9664\u9996\u4e2a\u5bf9\u9f50\u5b63\u5ea6；baseline\u4f7f\u7528\u5168\u90e8\u5b63\u5ea6"},
        {"key": "episode_rate", "label": "\u5e74\u5ea6\u8bb0\u5f55\u6bb5\u7387", "description": "\u53ea\u6709\u9010\u6b21\u8bb0\u5f55\u53ca\u5b8c\u6574\u6765\u6e90\u8986\u76d6\u65f6\u8ba1\u7b97；\u8bb0\u5f55\u6bb5\u4e0d\u7b49\u4e8e\u4e34\u5e8a\u590d\u53d1"}]
    for definition in summary["burden"]["feature_definitions"]:
        if definition["key"] == "aod":
            definition["description"] = "\u8bb0\u5f55\u7387、\u540e\u7eed\u6d3b\u8dc3\u5b63\u5ea6、\u9996\u672b\u8de8\u5ea6\u76f8\u5bf9\u51bb\u7ed3\u8bad\u7ec3\u75c5\u4f8b\u5206\u5e03\u7684\u4e09\u4e2a\u767e\u5206\u4f4d\u7b49\u6743\u5747\u503c"
    summary["consistency"]["audit_scope"] = "FOD\u5ba1\u8ba1\u4fdd\u7559\u622a\u81f3audit_as_of\u7684\u6709\u6548\u65e5\u671f；\u4e34\u5e8a\u7279\u5f81\u53e6\u6309\u516c\u5171date_follow_end\u53ca\u6b7b\u4ea1/\u5931\u8bbf\u622a\u65ad。ref\u662f\u7efc\u5408\u6765\u6e90\u53c2\u8003，\u4e0d\u662f\u4eba\u5de5\u88c1\u5b9a\u91d1\u6807\u51c6。"
    adu_values = [f["adu"] for f in feats]
    summary["adu"] = {"status": "available" if any(v is not None for v in adu_values) else "not_available",
                      **distribution(adu_values), "cohort_n": len(feats),
                      "source": meta.get("adu_source", meta.get("existing_adu_source", "existing phenotype export")),
                      "by_cluster": [{"id": row["id"], **distribution([f["adu"] for f in feats if f["cluster"] == row["id"]])}
                                     for row in summary["clusters"]["profiles"]],
                      "by_pattern": [{"id": row["id"], "label": row["label"], **distribution([f["adu"] for f in eligible if f["pattern"] == row["id"]])}
                                     for row in pattern_rows],
                      "correlations": [rank_correlation(feats, k) for k in ("age_at_index", "n_days", "n_bouts", "post_index_active_fraction", "aod")],
                      "interpretation": "ADuLT\u662f\u4f9d\u8d56\u60a3\u75c5\u5e74\u9f84\u53ca\u4eba\u7fa4\u53d1\u75c5\u98ce\u9669\u7684liability\u8868\u578b；\u6b64\u5904\u5bfc\u5165\u65e2\u6709Y.adu，\u4e0d\u91cd\u7b97，\u4e0d\u89e3\u91ca\u4e3a\u4e25\u91cd\u7a0b\u5ea6\u6216\u590d\u53d1。\u76f8\u5173\u6027\u662f\u63cf\u8ff0\u6027，\u4e14\u53ef\u80fd\u90e8\u5206\u6765\u81ea\u5e74\u9f84\u6784\u9020。"}
    valid = [f for f in feats if f["prognosis_reason"] == "eligible" and f["cluster"] is not None]
    summary["prognosis"] = {"status": "ready_for_validation" if valid else "no_eligible_landmark_cases",
                            "landmark_definition": "\u56fa\u5b9a\u7279\u5f81\u7a97\u6700\u540e\u4e00\u65e5；\u7ed3\u5c40\u4e25\u683c\u53d1\u751f\u4e8e\u968f\u540e\u65e5\u671f",
                            "n": len(valid), "deaths": sum(f["death_event"] for f in valid),
                            "train_n": sum(f["split"] == "train" for f in valid), "test_n": sum(f["split"] == "test" for f in valid),
                            "record_endpoint_n": sum(f["future_record_event"] is not None for f in valid),
                            "by_cluster": [{"id": row["id"], "n": sum(f["cluster"] == row["id"] for f in valid),
                                            "deaths": sum(f["death_event"] for f in valid if f["cluster"] == row["id"]),
                                            "person_years": sum(f["death_followup_years"] for f in valid if f["cluster"] == row["id"])}
                                           for row in summary["clusters"]["profiles"]],
                            "exclusions": [{"reason": key, "n": count} for key, count in
                                           sorted(Counter(f["prognosis_reason"] for f in feats if f["prognosis_reason"] != "eligible").items())],
                            "notes": ["\u4ec5landmark\u4e0d\u65e9\u4e8eUKB\u5165\u7ec4\u7684\u5b58\u6d3b\u8005\u6709\u9884\u540e\u7aef\u70b9；\u907f\u514d\u5165\u7ec4\u524d\u5fc5\u987b\u5b58\u6d3b\u9020\u6210\u7684\u504f\u501a。",
                                      "\u8fd9\u662flandmark\u5b58\u6d3b\u8005\u4e2d\u7684\u540e\u7eed\u9884\u540e；\u4e0d\u80fd\u5916\u63a8\u5230\u7a97\u53e3\u5185\u65e9\u671f\u6b7b\u4ea1\u8005。",
                                      "\u6b7b\u4ea1\u4eba\u6570\u53ca\u4eba\u5e74\u662f\u63cf\u8ff0\u7edf\u8ba1；\u6a21\u578b\u589e\u76ca\u987b\u5728\u76f8\u540c\u6d4b\u8bd5\u4eba\u7fa4\u4e0e\u76f8\u540c\u7ed3\u5c40/\u65f6\u70b9\u4e2d\u6bd4\u8f83。",
                                      "\u540e\u7eed\u75be\u75c5\u8bb0\u5f55\u7aef\u70b9\u9700\u8981\u5b8c\u6574\u6765\u6e90\u8986\u76d6；\u6b7b\u4ea1\u662f\u7ade\u4e89\u4e8b\u4ef6，\u5f53\u524d\u4ec5\u5bfc\u51facause-specific\u8bb0\u5f55\u4e8b\u4ef6。",
                                      "\u65e0\u57fa\u7ebfY\u7684\u53d1\u75c5\u98ce\u9669\u9884\u6d4b\u4e0d\u80fd\u4f7f\u7528\u672a\u6765Y.cluster；\u53ef\u628a\u672a\u6765\u6a21\u5f0f\u4f5c\u4e3a\u7ed3\u5c40，\u6216\u4f7f\u7528\u5176\u4ed6\u57fa\u7ebf\u75c5\u53f2。"]}
    summary["patterns"]["clinical_note"] = "\u540c\u4e00\u75be\u75c5\u53cd\u590d\u7f16\u7801\u53ef\u80fd\u6765\u81ea\u6162\u75c5\u6301\u7eed\u5b58\u5728、\u968f\u8bca\u6216\u5176\u4ed6\u4f4f\u9662\u65f6\u7684\u4f34\u968f\u8bca\u65ad；\u5e76\u4e0d\u81ea\u52a8\u7b49\u4e8e\u591a\u6b21\u6025\u6027\u53d1\u4f5c。"
    summary["methods"] = [
        {"title": "\u4e24\u6761\u65f6\u95f4\u8fb9\u754c", "text": "raw/date\u5ba1\u8ba1\u4f7f\u7528as_of；\u4e34\u5e8a\u7279\u5f81、landmark\u53ca\u7ed3\u5c40\u4f7f\u75280f/phenotype.R::date_follow_end\u5e76\u6309\u6b7b\u4ea1/\u5931\u8bbf\u622a\u65ad。"},
        {"title": "\u8bb0\u5f55\u5c42\u6b21", "text": "\u540c\u65e5\u53bb\u91cd；7/30/90\u5929\u8bb0\u5f55\u6bb5\u654f\u611f\u6027；first-per-code\u5386\u53f2\u53ea\u4f5c\u63a2\u7d22\u6027\u6a21\u5f0f，\u4e0d\u4f30\u8ba1\u771f\u5b9e\u590d\u53d1\u6b21\u6570。"},
        {"title": "\u8bad\u7ec3\u4e0e\u6d4b\u8bd5", "text": "\u6309\u56fa\u5b9aseed\u548cID\u54c8\u5e0c\u5206\u8bad\u7ec375%/\u6d4b\u8bd525%；\u7f29\u5c3e、\u6807\u51c6\u5316、K、\u4e2d\u5fc3\u53ca\u767e\u5206\u4f4d\u53c2\u8003\u4ec5\u8bad\u7ec3，\u6d4b\u8bd5\u51bb\u7ed3\u5e94\u7528。"},
        {"title": "\u672a\u6765\u5e94\u7528", "text": "--apply-model\u590d\u7528\u6a21\u578b；\u540c\u75be\u75c5/\u7a97\u53e3/\u6765\u6e90/\u7c92\u5ea6/\u7279\u5f81\u5b9a\u4e49\u624d\u53ef\u8d4b\u503c。\u7c07ID\u662f\u540d\u4e49\u7c7b\u522b，\u4e0d\u662f\u4e34\u5e8astage。"}]
    if any(p[k] is not None and p[k] > p["clinical_end"] for p in patients.values() for k in SOURCES):
        summary["warnings"].insert(0, "\u90e8\u5206FOD\u5ba1\u8ba1\u65e5\u671f\u665a\u4e8e\u4e34\u5e8a\u968f\u8bbf\u622a\u6b62；\u5b83\u4eec\u4fdd\u7559\u7528\u4e8e\u65e5\u671f\u6765\u6e90\u6838\u67e5，\u4e0d\u8fdb\u5165\u4e34\u5e8a\u7279\u5f81\u6216\u9884\u540e\u7aef\u70b9。")
    return summary


def make_demo(work, trait, seed):
    """Synthetic fixture, never a substitute for real CAD results."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    ps, es, cov = [], [], []
    for i in range(480):
        eid = "DEMO_" + str(i + 1).zfill(4)
        group = i // 120
        first = date(2009 + i % 3, 3, 1).toordinal() + int(rng.integers(0, 90))
        if group == 0:
            ds = [first]
        elif group == 1:
            ds = [first + int(v) for v in sorted(set([0] + rng.integers(1, 100, 6).tolist()))]
        elif group == 2:
            ds = [first + int(v) for v in [0, 365, 730, 1095, 1460]]
        else:
            ds = [first + int(v) for v in range(0, 1800, 40 + i % 20)]
        gp, hospital = [], []
        for j, d in enumerate(ds):
            src = "icd10"
            (gp if src == "gp2" else hospital).append(d)
            if j % 2 == 0:
                gp.append(d)
                es.append(dict(eid=eid, source="gp2", code="G3", date=iso(d), date_raw=iso(d), date_status="valid", origin="SYNTHETIC", granularity="record", match_local="TRUE", match_official="TRUE"))
            es.append(dict(eid=eid, source=src, code="G3" if src == "gp2" else "I259", date=iso(d), date_raw=iso(d), date_status="valid", origin="SYNTHETIC", granularity="record" if src == "gp2" else "episode", match_local="TRUE", match_official="TRUE"))
        hs = min(hospital) if hospital else None
        gs = min(gp) if gp else None
        ref = first + (0 if i % 5 else 365)
        ps.append(dict(eid=eid, date_birth=f"{1940+i%30}-06-15", date_assess="2010-07-01", sex=i%2,
                       date_death=iso(first+int((4+i%7)*DAYS_YEAR)) if i%4 == 0 else "", date_lost="", icd10=iso(hs), gp=iso(gs), ref=iso(ref), min1=iso(first), min2=iso(ref), adu=round(float(rng.normal()), 4)))
        for src in ("hospital", "gp"):
            cov.append(dict(eid=eid, source=src, start="1995-01-01", end="2025-12-31"))
    write_tsv(work / "patients.tsv.gz", ["eid", "date_birth", "date_assess", "date_death", "date_lost", "sex", "icd10", "gp", "ref", "min1", "min2", "adu"], ps)
    write_tsv(work / "events.tsv.gz", ["eid", "source", "code", "date", "date_raw", "date_status", "origin", "granularity", "match_local", "match_official"], es)
    write_tsv(work / "coverage.tsv", ["eid", "source", "start", "end"], cov)
    md = {"Y": trait, "demo": True, "as_of": "2025-12-31", "date_follow_end": "2023-04-01", "raw_hospital_episode_dates_available": True, "raw_icd10_available": True, "raw_gp2_available": True}
    write_tsv(work / "metadata.tsv", ["key", "value"], [{"key": k, "value": json.dumps(v)} for k, v in md.items()])
    return str(work / "coverage.tsv")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--Y", required=True)
    p.add_argument("--coverage")
    p.add_argument("--as-of")
    p.add_argument("--followup-end", help="Internal resolved 0f/phenotype.R::date_follow_end")
    p.add_argument("--trajectory-source", choices=["hospital", "gp", "combined"], default="hospital")
    p.add_argument("--apply-model", help="Assign with frozen model.json; never refit")
    p.add_argument("--mode", choices=["history", "baseline"], default="history")
    p.add_argument("--lookback-years", type=float, default=5)
    p.add_argument("--landmark-years", type=float, default=2)
    p.add_argument("--episode-gap", type=int, default=30)
    p.add_argument("--k", default="auto", choices=["auto"] + list("2345678"))
    p.add_argument("--seed", type=int, default=20261111)
    p.add_argument("--min-cluster-n", type=int, default=20)
    p.add_argument("--cluster-sample", type=int, default=5000)
    p.add_argument("--stability-repeats", type=int, default=8)
    p.add_argument("--demo", action="store_true")
    a = p.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9_]+", a.Y):
        p.error("Y must contain letters, numbers and underscores only")
    if not a.demo and not a.input_dir:
        p.error("--input-dir is required unless --demo")
    if not all(math.isfinite(x) and 1 / 365 <= x <= 100 for x in (a.lookback_years, a.landmark_years)):
        p.error("Year windows must be finite, positive, and at most 100")
    if a.episode_gap < 0 or a.seed < 0 or a.min_cluster_n < 2 or a.cluster_sample < 40 or a.stability_repeats < 2:
        p.error("Require gap>=0, seed>=0, min-cluster-n>=2, cluster-sample>=40, stability-repeats>=2")
    if a.as_of and ordinal(a.as_of) is None:
        p.error("--as-of requires a valid ISO date, not a UKB placeholder")
    if a.followup_end and ordinal(a.followup_end) is None:
        p.error("--followup-end requires the valid resolved shared date")
    return a


def main(argv=None):
    args = parse_args(argv)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if args.demo:
        args.input_dir = str(out)
        args.coverage = make_demo(out, args.Y, args.seed)
    root = Path(args.input_dir)
    meta = read_meta(root / "metadata.tsv")
    if meta.get("Y") and meta["Y"] != args.Y:
        raise ValueError("Requested Y does not match extracted metadata; refusing mislabeled analysis")
    asof_text = resolve_asof(meta, args.as_of)
    asof = ordinal(asof_text)
    args.resolved_asof = asof_text
    args.followup_end = args.followup_end or meta.get("date_follow_end")
    if ordinal(args.followup_end) is None:
        raise ValueError("Missing date_follow_end: use aod.sh to resolve 0f/phenotype.R or re-extract the bundle")
    args.apply_model_payload = json.loads(Path(args.apply_model).read_text(encoding="utf-8")) if args.apply_model else None
    if args.apply_model and (not isinstance(args.apply_model_payload, dict)
                            or args.apply_model_payload.get("feature_schema") != "aod_trajectory_v2"
                            or not isinstance(args.apply_model_payload.get("k"), int)
                            or args.apply_model_payload["k"] < 1):
        raise ValueError("--apply-model requires a fitted AOD v2 model.json")
    patients = read_patients(root / "patients.tsv.gz", asof, args.mode, ordinal(args.followup_end))
    selected, levels, qc, procedure, local = read_events(root / "events.tsv.gz", patients, meta)
    definitions = meta.get("definitions", {}).get("sources", {})
    args.partial_sources = [group for group, members in (("hospital", ("icd10", "icd9")), ("gp", ("gp2", "gp3")))
                            if levels[group] == 3 and any(truth(definitions.get(s, {}).get("local_defined")) and
                                                         not truth(meta.get("raw_" + s + "_available")) for s in members)]
    args.source_levels = levels
    coverage = read_coverage(args.coverage)
    feats, used_sources, nbins, window = feature_rows(patients, selected, levels, coverage, args)
    score = score_burden(feats, window, args.apply_model_payload.get("score_model") if args.apply_model_payload else None)
    cl, model = cluster(feats, used_sources, nbins, window, args)
    add_prognosis(feats, patients, selected, coverage, used_sources, levels, args)
    model.update(trait=args.Y, mode=args.mode, as_of=asof_text, landmark_years=args.landmark_years, lookback_years=args.lookback_years,
                 date_follow_end=args.followup_end, trajectory_source=args.trajectory_source, window_days=window,
                 episode_gap=args.episode_gap, source_levels=levels, score_model=score, coverage_supplied=bool(args.coverage),
                 source_base_commit="fe651ba4de3b5586d442c67b92f541fde5d14aac", numpy_version=np.__version__,
                 analysis_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), export_as_of=meta["as_of"],
                 raw_code_sources={s: truth(meta.get("raw_" + s + "_available")) for s in ("icd10", "icd9", "gp2", "gp3")},
                 partial_sources=args.partial_sources, input_provenance=meta.get("inputs", []) + meta.get("r_inputs", []),
                 limitations=["Descriptive model fitted within this disease-evidence cohort; not a validated severity or causal phenotype.",
                              "Future data must not be used to predict a prior baseline outcome. Baseline event-date truncation cannot remove retrospectively backdated ascertainment.",
                              "All scores are case-conditional; excluded/undated/no-coverage cases and noncases are not zero-scored controls."],
                 input_files=[{"role": name, "bytes": (root / name).stat().st_size} for name in ("patients.tsv.gz", "events.tsv.gz", "metadata.tsv")])
    if args.coverage:
        h = hashlib.sha256()
        with open(args.coverage, "rb") as fp:
            for block in iter(lambda: fp.read(1024 * 1024), b""):
                h.update(block)
        model["coverage_sha256"] = h.hexdigest()
    summary = aggregate(patients, feats, levels, selected, qc, cl, score, model, meta, args, nbins, window)
    summary = augment_summary(summary, feats, patients, levels, used_sources, model, meta, args, window)
    public_fields = ["eid", "window_start", "window_end", "index_date", "n_days_total", "n_days", "n_bouts", "n_hospital_days", "n_gp_days", "span_years",
                     "median_gap_days", "max_gap_days", "active_fraction", "late_minus_early", "time_since_last_days", "hospital_fraction", "gp_fraction",
                     "hospital_observed_years", "gp_observed_years", "event_rate", "hospital_rate", "gp_rate", "recent_rate", "early_isolate", "n_days_without_isolate",
                     "horizon_complete", "coverage_complete", "record_scope", "hospital_mode", "gp_mode", "aod_reason", "cluster_reason",
                     "feature_mode", "clinical_end", "analysis_case", "pattern_reason", "pattern_evidence", "recurrence_capable",
                     "post_index_active_fraction", "early_active_fraction", "late_active_fraction", "episode_rate",
                     "split", "assignment_margin", "landmark_date", "age_at_index", "age_at_landmark", "sex", "disease_duration_years", "source_duration_years", "fod_for_prediction",
                     "prognosis_reason", "death_followup_years", "death_event", "record_prognosis_reason",
                     "future_record_time_years", "future_record_event", "future_record_competing_death",
                     "later_record_days", "later_record_audit_followup_years"]
    output_rows = []
    for f in feats:
        p = patients[f["eid"]]
        row = {k: clean_json(f.get(k)) for k in public_fields}
        row.update({"fod_" + k: iso(p[k]) for k in SOURCES})
        row.update({args.Y + ".aod": f["aod"], args.Y + ".cluster": f["cluster"], args.Y + ".adu": f["adu"], args.Y + ".pattern": f["pattern"]})
        output_rows.append(row)
    public_fields += ["fod_" + k for k in SOURCES] + [args.Y + ".aod", args.Y + ".cluster", args.Y + ".pattern", args.Y + ".adu"]
    write_tsv(out / "aod.tsv.gz", public_fields, output_rows)
    profile_fields = ["id", "label", "n", "pct", "train_n", "test_n", "days_median", "bouts_median", "span_years_median",
                      "active_fraction_median", "post_index_active_fraction_median", "early_active_fraction_median",
                      "late_active_fraction_median", "late_minus_early_median", "recency_days_median", "max_gap_days_median",
                      "hospital_fraction_median", "margin_median"]
    write_tsv(out / "clusters.tsv", profile_fields, cl["profiles"])
    save_json(out / "model.json", model)
    save_json(out / "summary.json", summary)
    print(f"AOD: {len(patients):,} evidence cases; {sum(f['n_days'] is not None and f['n_days'] > 0 for f in feats):,} with dates in window; "
          f"cluster status={cl['status']}, K={cl['k']}; scored={score['n']}. Report data contain aggregates only.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError) as exc:
        raise SystemExit("AOD ERROR: " + str(exc)) from exc

