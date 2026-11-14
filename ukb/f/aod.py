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

VERSION = "1.1.0"
DAYS_YEAR = 365.25
SPECIAL = {"1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"}
MISSING = {"", "NA", "NaN", "nan", ".", "NULL", "null", "N/A"}
SOURCES = {"icd10": "fod.icd10", "ref": "fod.ref", "icd9": "fod.icd9", "gp": "fod.gp",
           "srd": "fod.srd", "opcs4": "fod.opcs4", "min1": "min1 · 本地最早", "min2": "min2 · 映射候选"}
GROUP = {"icd10": "hospital", "icd9": "hospital", "gp": "gp", "gp2": "gp", "gp3": "gp"}
LEVEL_LABEL = {0: "first_date_only", 1: "first_per_code", 2: "cached_record_days", 3: "linked_records"}
LAG_EDGES = [-3653, -366, -31, -1, 0, 30, 365, 3652]
LAG_LABELS = ["≤−10年", "−10至−1年", "−1年至−31天", "−30至−1天", "同日", "1至30天", "31天至1年", "1至10年", ">10年"]
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
        notes.append(f"ref 输入缺少 {len(missing)} 个预期字段；当前一致率仅基于已提供字段，不能代表完整定义。")
    if mode == "saved_only":
        notes.append("ref 仅使用既存缓存；更新定义后需重新提取原始字段并生成 fod_ref，才能核对新范围。")
    elif mode == "unavailable":
        notes.append("当前输入没有可用 ref 字段或缓存，不能计算参考日期一致率。")
    elif mode == "unknown":
        notes.append("导出未记录 ref 输入模式；无法核实参考日期来自原始字段还是缓存。")
    elif not expected:
        notes.append("ref 使用原始字段，但导出未记录预期字段列表，完整性尚未核实。")
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


def read_patients(path, asof, mode):
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
        p["cutoff"] = p["audit_cutoff"]
        if mode == "baseline":
            p["cutoff"] = min(p["cutoff"], p["baseline"]) if p["baseline"] is not None else None
        for key in SOURCES:
            val = ordinal(p.get(key, ""))
            if val is not None and (val > p["audit_cutoff"] or (p["birth"] is not None and val < p["birth"])):
                val = None
            p[key] = val
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
        if p["cutoff"] is None or d > p["cutoff"]:
            qc["after_analysis_cutoff"] += 1
            continue
        if not local:
            continue
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


def feature_rows(patients, selected, levels, coverage, args):
    window = int(round((args.lookback_years if args.mode == "baseline" else args.landmark_years) * DAYS_YEAR))
    nbins = max(4, min(40, int(round(window / DAYS_YEAR * 4))))
    used_sources = [s for s in selected if levels[s] >= 1]
    partial_sources = getattr(args, "partial_sources", [])
    can_measure_encounters = bool(used_sources) and all(levels[s] >= 2 for s in used_sources) and not partial_sources
    feats = []
    for eid, p in patients.items():
        full = {s: sorted(selected[s].get(eid, ())) if levels[s] >= 1 else [] for s in selected}
        observed = sorted(set(full["hospital"]) | set(full["gp"]))
        anchor = p["baseline"] if args.mode == "baseline" else (observed[0] if observed else None)
        start = anchor - window + 1 if args.mode == "baseline" and anchor else anchor
        end = anchor if args.mode == "baseline" else (start + window - 1 if start else None)
        f = {"eid": eid, "window_start": iso(start), "window_end": iso(end), "index_date": iso(observed[0]) if observed else None,
             "n_days_total": len(observed), "record_scope": "+".join(used_sources), "hospital_mode": LEVEL_LABEL.get(levels["hospital"], "unavailable"),
             "gp_mode": LEVEL_LABEL.get(levels["gp"], "unavailable"), "cluster": None, "cluster_reason": "no_clinical_aod",
             "aod": None, "aod_reason": "no_clinical_aod", "adu": num(p.get("adu")), "_start": start, "_end": end,
             "_trajectory": np.zeros(nbins),
             "_hospital_trajectory": np.zeros(nbins) if "hospital" in used_sources else None,
             "_gp_trajectory": np.zeros(nbins) if "gp" in used_sources else None}
        for key in ("n_days", "n_bouts", "span_years", "median_gap_days", "max_gap_days", "active_fraction", "hospital_fraction", "gp_fraction",
                    "event_rate", "hospital_rate", "gp_rate", "recent_rate", "time_since_last_days", "late_minus_early", "early_isolate", "n_days_without_isolate"):
            f[key] = None
        f["coverage_complete"] = False
        f["horizon_complete"] = False
        if not observed or start is None or end is None or p["cutoff"] is None:
            if args.mode == "baseline" and p["baseline"] is None:
                f["cluster_reason"] = f["aod_reason"] = "missing_baseline"
            feats.append(f)
            continue
        f["horizon_complete"] = p["cutoff"] >= end
        win = {s: [d for d in ds if start <= d <= min(end, p["cutoff"])] for s, ds in full.items()}
        days = sorted(set(win["hospital"]) | set(win["gp"]))
        f["n_days"] = len(days)
        f["n_bouts"] = len(bouts(days, args.episode_gap))
        f["span_years"] = (days[-1] - days[0]) / DAYS_YEAR if days else 0.0
        gaps = np.diff(days).tolist()
        f["median_gap_days"] = med(gaps)
        f["max_gap_days"] = max(gaps) if gaps else None
        f["time_since_last_days"] = min(end, p["cutoff"]) - days[-1] if days else None
        f["n_hospital_days"] = len(win["hospital"]) if "hospital" in used_sources else None
        f["n_gp_days"] = len(win["gp"]) if "gp" in used_sources else None
        # Source shares use source-days, so a same-day hospital + GP contact
        # contributes one to each source, but only one to the combined count.
        source_days = len(win["hospital"]) + len(win["gp"])
        both_sources = set(used_sources) == {"hospital", "gp"}
        f["hospital_fraction"] = len(win["hospital"]) / source_days if source_days and both_sources else None
        f["gp_fraction"] = len(win["gp"]) / source_days if source_days and both_sources else None
        for name, ds in (("_trajectory", days), ("_hospital_trajectory", win["hospital"]), ("_gp_trajectory", win["gp"])):
            for d in ds:
                f[name][min(nbins - 1, int((d - start) / window * nbins))] = 1
        f["active_fraction"] = float(f["_trajectory"].mean())
        half = nbins // 2
        f["late_minus_early"] = float(f["_trajectory"][half:].mean() - f["_trajectory"][:half].mean())
        confirmed = True
        for s in used_sources:
            intervals = coverage.get((eid, s), [])
            obsdays = coverage_days(intervals, start, min(end, p["cutoff"]))
            f[s + "_observed_years"] = obsdays / DAYS_YEAR if intervals else None
            if intervals and obsdays and levels[s] >= 2 and s not in partial_sources:
                n_observed = sum(any(a <= d <= b for a, b in intervals) for d in win[s])
                f[s + "_rate"] = n_observed / (obsdays / DAYS_YEAR)
            confirmed = confirmed and bool(intervals) and obsdays == window
        f["coverage_complete"] = confirmed and f["horizon_complete"]
        # A known gap is not an inactive quarter. Visible record counts and
        # intervals remain descriptive, but fixed-window activity/recency need
        # a complete horizon and, when supplied, complete source coverage.
        if not f["horizon_complete"] or (args.coverage and not f["coverage_complete"]):
            for key in ("active_fraction", "late_minus_early", "time_since_last_days"):
                f[key] = None
        # Long-gap singleton is flagged, never deleted from the main result.
        first = observed[0]
        corroborated = all(any(abs(d - first) <= 365 for d in full[s]) for s in ("hospital", "gp"))
        suspect = len(observed) > 1 and observed[1] - first > 5 * DAYS_YEAR and not corroborated
        f["early_isolate"] = int(suspect)
        altered = observed[1:] if suspect else observed
        if args.mode == "history" and altered:
            f["n_days_without_isolate"] = sum(altered[0] <= d < altered[0] + window for d in altered)
        else:
            f["n_days_without_isolate"] = sum(start <= d <= end for d in altered)
        if f["coverage_complete"] and can_measure_encounters:
            f["event_rate"] = len(days) / (window / DAYS_YEAR)
            last_year = min(window, int(round(DAYS_YEAR)))
            f["recent_rate"] = sum(d > end - last_year for d in days) / (last_year / DAYS_YEAR)
            f["aod_reason"] = "eligible"
        else:
            f["aod_reason"] = "partial_code_source_extract" if partial_sources else ("first_per_code_not_complete_records" if not can_measure_encounters else "unverified_or_incomplete_coverage")
        if not days and not (args.mode == "baseline" and f["coverage_complete"] and can_measure_encounters):
            f["cluster_reason"] = "no_dates_in_window"
        elif not f["horizon_complete"]:
            f["cluster_reason"] = "death_loss_or_cutoff_before_window_end"
        elif args.coverage and not f["coverage_complete"]:
            f["cluster_reason"] = "incomplete_source_coverage"
        else:
            f["cluster_reason"] = "eligible"
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


def score_burden(feats, window):
    eligible = [f for f in feats if f["aod_reason"] == "eligible"]
    model = {"status": "unavailable", "n": len(eligible), "components": ["log1p_event_rate", "active_fraction", "span_fraction"],
             "weights": [1 / 3] * 3, "method": "mean of within-analysis-cohort average-rank percentiles; cases only; exploratory recorded activity, not severity"}
    if len(eligible) < 20:
        for f in eligible:
            f["aod_reason"] = "fewer_than_20_comparable_cases"
        model["n"] = 0
        return model
    x = np.array([[math.log1p(f["event_rate"]), f["active_fraction"], f["span_years"] * DAYS_YEAR / window] for f in eligible])
    scores = np.mean(np.column_stack([percentile_values(x[:, j]) for j in range(3)]), axis=1)
    for f, s in zip(eligible, scores):
        f["aod"] = float(s)
        f["aod_reason"] = "scored"
    model.update(status="exploratory", reference_values=[np.sort(x[:, j]).tolist() for j in range(3)],
                 reference_n=len(eligible), n=len(eligible), scope=eligible[0]["record_scope"])
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
    names = ["log1p_days", "active_fraction", "span_fraction", "gap_fraction", "late_minus_early"]
    labels = ["日期数", "活跃时间段", "首末跨度", "相邻间隔", "后期−早期活动"]
    if include_source:
        names.append("hospital_share")
        labels.append("住院记录占比")
    x = []
    for f in feats:
        v = [math.log1p(f["n_days"]), f["active_fraction"], f["span_years"] * DAYS_YEAR / window,
             (f["median_gap_days"] if f["median_gap_days"] is not None else window) / window, f["late_minus_early"]]
        if include_source:
            v.append(f["hospital_fraction"] if f["hospital_fraction"] is not None else .5)
        x.append(v)
    return np.asarray(x, dtype=float), names, labels


def trajectory_mean(feats, field):
    values = [f[field] for f in feats if f.get(field) is not None]
    return np.mean(values, axis=0) * 100 if values else None


def cluster(feats, used_sources, nbins, window, args):
    eligible = [f for f in feats if f["cluster_reason"] == "eligible"]
    excluded = Counter(f["cluster_reason"] for f in feats if f["cluster_reason"] != "eligible")
    include_source = bool(args.coverage) and set(used_sources) == {"hospital", "gp"} and not getattr(args, "partial_sources", [])
    scope = "覆盖期已确认的固定窗口" if args.coverage else "可见日期模式 · 来源覆盖期未确认"
    result = {"status": "insufficient_data", "scope": scope, "k": 0, "n": len(eligible), "silhouette": None, "stability_ari": None,
              "candidates": [], "profiles": [], "excluded": [{"reason": k, "n": n} for k, n in sorted(excluded.items())],
              "feature_labels": [], "bin_labels": [], "projection": {"points": [], "explained": []}, "notes": [],
              "source_availability": {s: s in used_sources for s in ("hospital", "gp")}}
    model = {"version": VERSION, "fit_scope": scope, "seed": args.seed, "source_mix_feature_used": include_source,
             "source_share_zero_total_rule": "neutral 0.5 only for the model when total contacts are zero; raw source shares remain missing",
             "status": "insufficient_data", "k": 0, "settings": {"winsor_quantiles": [.01, .99], "min_silhouette": .15, "min_stability_ari": .60,
             "fit_sample_limit": args.cluster_sample, "minimum_cluster_n": args.min_cluster_n, "stability_repeats": args.stability_repeats}}
    years = window / DAYS_YEAR
    result["bin_labels"] = [f"{i * years / nbins:.1f}–{(i + 1) * years / nbins:.1f}年" for i in range(nbins)]
    if not include_source:
        result["notes"].append("未确认两类来源覆盖时，来源比例只作分群后描述，不作为聚类特征；GP无记录不被当作不看GP。")
    minimum = max(2 * args.min_cluster_n, 40)
    if len(eligible) < minimum:
        result["notes"].append(f"可用人数不足 {minimum}，保留描述统计，不强制分群。")
        for f in eligible:
            f["cluster_reason"] = "too_few_comparable_cases"
        return result, model
    rng = np.random.default_rng(args.seed)
    x, names, labels = cluster_features(eligible, window, include_source)
    fit_idx = np.sort(rng.choice(len(x), min(len(x), args.cluster_sample), replace=False))
    rawfit = x[fit_idx]
    lo, hi = np.quantile(rawfit, [.01, .99], axis=0)
    clipped = np.clip(rawfit, lo, hi)
    center = np.median(clipped, axis=0)
    scale = np.quantile(clipped, .75, axis=0) - np.quantile(clipped, .25, axis=0)
    sd = np.std(clipped, axis=0)
    scale = np.where(scale > 1e-8, scale, sd)
    keep = scale > 1e-8
    if not np.any(keep):
        result.update(status="no_separation", k=1)
        result["notes"].append("特征没有足够变异；没有证据支持多个稳定亚群。")
        for f in eligible:
            f["cluster"] = "C1"
            f["cluster_reason"] = "no_stable_partition"
        result["profiles"] = [{"id": "C1", "label": "未发现稳定分型", "n": len(eligible), "pct": 100.0,
            "days_median": med([f["n_days"] for f in eligible]), "bouts_median": med([f["n_bouts"] for f in eligible]),
            "span_years_median": med([f["span_years"] for f in eligible]), "active_fraction_median": med([f["active_fraction"] for f in eligible]),
            "hospital_fraction_median": med([f["hospital_fraction"] for f in eligible]),
            "trajectory": trajectory_mean(eligible, "_trajectory"),
            "hospital_trajectory": trajectory_mean(eligible, "_hospital_trajectory"),
            "gp_trajectory": trajectory_mean(eligible, "_gp_trajectory"),
            "features": []}]
        model.update(status="no_separation", k=1, constant_model=True, all_feature_names=names,
                     constant_values=center, feature_names=[], center=[], scale=[], winsor_low=[], winsor_high=[],
                     centers=[[]], label_order=["C1"], n=len(eligible), fit_n=len(fit_idx))
        return result, model
    z = (np.clip(x[:, keep], lo[keep], hi[keep]) - center[keep]) / scale[keep]
    fit = z[fit_idx]
    model.update(feature_names=[n for n, k in zip(names, keep) if k], winsor_low=lo[keep], winsor_high=hi[keep],
                 center=center[keep], scale=scale[keep], fit_n=len(fit), n=len(x),
                 clipped_values=int(((x[:, keep] < lo[keep]) | (x[:, keep] > hi[keep])).sum()))
    result["feature_labels"] = [n for n, k in zip(labels, keep) if k]
    validation_idx = np.sort(rng.choice(len(z), min(800, len(z)), replace=False))
    ev = z[validation_idx]
    maxk = min(8, len(np.unique(fit, axis=0)), len(eligible) // args.min_cluster_n)
    candidates = list(range(2, min(6, maxk) + 1)) if args.k == "auto" else ([int(args.k)] if int(args.k) <= maxk else [])
    accepted = []
    for k in candidates:
        c = kmeans(fit, k, rng)
        if c is None:
            continue
        lab, _ = assign(z, c)
        cnt = np.bincount(lab, minlength=k)
        s = silhouette(ev, lab[validation_idx])
        aris = []
        for _ in range(args.stability_repeats):
            ii = rng.choice(len(fit), max(k, int(.8 * len(fit))), replace=False)
            cc = kmeans(fit[ii], k, rng, starts=3)
            aris.append(adjusted_rand(lab[validation_idx], assign(ev, cc)[0]) if cc is not None else 0.0)
        ari = float(np.mean(aris))
        ok = int(cnt.min()) >= max(args.min_cluster_n, int(math.ceil(.01 * len(z)))) and s is not None and s >= .15 and ari >= .60
        info = {"k": k, "silhouette": s, "stability_ari": ari, "min_n": int(cnt.min()), "accepted": ok}
        result["candidates"].append(info)
        if ok:
            accepted.append((float(s), k, ari, c, lab))
    if accepted:
        top = max(v[0] for v in accepted)
        # Prefer a simpler partition when silhouette differs by at most .02.
        chosen = min((v for v in accepted if v[0] >= top - .02), key=lambda a: a[1])
        sil, k, ari, c, lab = chosen
        order = sorted(range(k), key=lambda j: (float(np.median([eligible[i]["n_days"] for i in np.flatnonzero(lab == j)])), j))
        inverse = {old: new for new, old in enumerate(order)}
        lab = np.array([inverse[int(v)] for v in lab])
        c = c[order]
        result.update(status="ok", k=k, silhouette=sil, stability_ari=ari)
        model.update(status="ok", k=k, centers=c, label_order=["C" + str(i + 1) for i in range(k)], silhouette=sil, stability_ari=ari)
    else:
        k = 1
        lab = np.zeros(len(z), dtype=int)
        c = np.mean(fit, axis=0)[None, :]
        result.update(status="no_stable_partition", k=1)
        model.update(status="no_stable_partition", k=1, centers=c, label_order=["C1"])
        result["notes"].append("候选分群未同时通过分离度、稳定性和最小组人数要求，因此保留一个未分型组。")
    for f, l in zip(eligible, lab):
        f["cluster"] = "C" + str(l + 1)
        f["cluster_reason"] = "assigned" if k > 1 else "no_stable_partition"
    feature_percentiles = np.column_stack([percentile_values(z[:, j]) / 100 for j in range(z.shape[1])])
    for j in range(k):
        idx = np.flatnonzero(lab == j)
        fs = [eligible[i] for i in idx]
        dn = med([f["n_days"] for f in fs])
        span = med([f["span_years"] for f in fs])
        act = med([f["active_fraction"] for f in fs])
        gp = med([f["gp_fraction"] for f in fs])
        name = "窗内未见记录" if dn == 0 else ("单一日期为主" if dn <= 1 else ("跨年持续记录" if act >= .45 else ("跨年间歇记录" if span >= 1 else "记录集中")))
        if include_source and gp is not None and gp >= .8:
            name += " · 可见GP为主"
        if k == 1:
            name = "未发现稳定分型"
        result["profiles"].append({"id": "C" + str(j + 1), "label": name, "n": len(fs), "pct": percent(len(fs), len(eligible)),
             "days_median": dn, "bouts_median": med([f["n_bouts"] for f in fs]), "span_years_median": span,
             "active_fraction_median": act, "hospital_fraction_median": med([f["hospital_fraction"] for f in fs]),
             "trajectory": trajectory_mean(fs, "_trajectory"),
             "hospital_trajectory": trajectory_mean(fs, "_hospital_trajectory"),
             "gp_trajectory": trajectory_mean(fs, "_gp_trajectory"),
             "features": np.median(feature_percentiles[idx], axis=0)})
    pca_mean = fit.mean(axis=0)
    _, singular, vectors = np.linalg.svd(fit - pca_mean, full_matrices=False)
    dims = min(2, vectors.shape[0])
    v = vectors[:dims].copy()
    for axis in range(dims):
        if v[axis, np.argmax(np.abs(v[axis]))] < 0:
            v[axis] *= -1
    xy = (z - pca_mean) @ v.T
    if dims == 1:
        xy = np.column_stack([xy[:, 0], np.zeros(len(xy))])
    explained = singular[:dims] ** 2 / max(float((singular ** 2).sum()), 1e-12) * 100
    limits = np.quantile(xy, [.005, .995], axis=0)
    width = np.maximum(limits[1] - limits[0], 1e-9)
    points = Counter()
    for row, label in zip(xy, lab):
        bins = np.minimum(29, np.maximum(0, ((row - limits[0]) / width * 30).astype(int)))
        points[int(bins[0]), int(bins[1]), int(label)] += 1
    # Only occupied aggregate cells, never per-participant points or identifiers.
    result["projection"] = {"points": [{"x": float(limits[0, 0] + (bx + .5) / 30 * width[0]),
          "y": float(limits[0, 1] + (by + .5) / 30 * width[1]), "n": n, "cluster": "C" + str(cl + 1)}
          for (bx, by, cl), n in sorted(points.items())], "explained": explained.tolist()}
    model.update(pca_mean=pca_mean, pca_loadings=v, pca_explained=explained, candidates=result["candidates"])
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
        headline = f"共同有日期者中，ICD10 与 ref 同日率 {default['same_day_pct']:.1f}%；差值反映来源与定义差异。"
    else:
        headline = "当前输入没有可比较的 ICD10–ref 配对；可切换来源，缺失不计作一致。"
    withdays = [f for f in feats if f["n_days"] is not None and (f["n_days"] > 0 or f["coverage_complete"])]
    reasons = {"quiet": "窗内未见记录", "single": "单一日期", "compact": "短期集中", "intermittent": "跨期间歇", "persistent": "持续被记录"}
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
    labels = [f"{i * window / DAYS_YEAR / nbins:.1f}" for i in range(nbins)]
    tl = [{"id": key, "label": reasons[key], "n": len(fs), "values": trajectory_mean(fs, "_trajectory")} for key, fs in timeline_groups.items()]
    timeline_n = sum(len(fs) for fs in timeline_groups.values())
    notes = ["下方时间形态是解释用规则分组；第三部分才是数据学习得到的无监督 cluster。",
             "同一人同一天的多代码与GP双编码只计一个记录日；相邻间隔不超过阈值的记录日合并为记录簇，并非临床住院次数或新发事件。",
             "合法的极早日期只标记并做敏感性比较；模型特征按训练样本1%/99%缩尾，原始日期、人数和统计值不被删除。",
             "Y.adu表示年龄相关liability；Y.aod是病例内部探索性的可观察记录负担，不替换ADuLT，也不能证明治愈或严重程度。"]
    if not args.coverage:
        notes.insert(0, "未提供来源观察区间：展示可见日期模式，年度事件率与 Y.aod 留空；不能把未记录时间解释为健康时间。")
    else:
        notes.insert(0, "固定窗口轨迹仅纳入所有已用来源覆盖完整者；覆盖不完整者保留可见记录统计，但不将未观察时间绘成零活动。")
    if any(level == 1 for level in levels.values()):
        notes.insert(0, "医院汇总AOD包含各诊断代码的首次日期，遗漏同一代码后续住院；对应分型是保留日期模式。")
    if args.mode == "baseline":
        notes.append("baseline特征和聚类仅使用事件日期不晚于首次评估的记录；第一部分日期审计保留完整有效FOD。缺少记录录入时间时，仍不能保证回溯记录在基线当时已经可获得。")
    warnings = list(reference["notes"])
    if meta.get("reference_scope_warning"):
        warnings.append(str(meta["reference_scope_warning"]))
    if meta.get("official_maps_available") is False:
        warnings.append("未提供公开映射表：保存的FOD仍可比较，min2仅是当前可重建部分的候选日期。")
    if not args.coverage:
        warnings.append("来源观察期未知：聚类为探索性日期模式；Y.aod与联合事件率未计算。")
    elif timeline_excluded["incomplete_source_coverage"]:
        warnings.append(f"已提供观察期，但有 {timeline_excluded['incomplete_source_coverage']} 人的固定窗口覆盖不完整；未纳入完整窗口轨迹，未观察时段不补零。")
    if any(level == 1 for level in levels.values()):
        warnings.append("至少一个来源只有每代码首次日期，不能解释成完整重复就诊历史。")
    if getattr(args, "partial_sources", []):
        warnings.append("原始表未导出疾病定义需要的全部细来源代码列：" + ", ".join(args.partial_sources) + "。保留该提取范围的日期模式，综合负担分数不计算。")
    eligible_score = [f["aod"] for f in feats if f["aod"] is not None]
    suspect = [f for f in feats if f["early_isolate"]]
    qc_labels = {"after_death": "死亡之后", "before_birth": "出生之前", "after_as_of": "晚于截点", "after_analysis_cutoff": "分析窗外", "unparseable": "无法解析", "unknown": "日期未知", "hospital_end_date_proxy": "住院结束日代理"}
    return {"schema_version": 1, "trait": args.Y, "generated_at": datetime.now(timezone.utc).isoformat(), "demo": truth(meta.get("demo")), "mode": args.mode,
            "cohort": {"cases": len(feats), "with_clinical_dates": sum(f["n_days"] > 0 for f in withdays), "clustered": sum(f["cluster"] is not None for f in feats), "score_eligible": len(eligible_score)},
            "availability": [{"source": s, "label": "住院" if s == "hospital" else "GP", "mode": LEVEL_LABEL.get(level, "unavailable"),
                              "n_people": sum(bool(v) for v in selected[s].values()), "n_days": sum(len(v) for v in selected[s].values())} for s, level in levels.items()],
            "input_date_ranges": meta.get("input_date_ranges", []),
            "warnings": warnings, "consistency": {"sources": sources, "default_pair": ["icd10", "ref"], "pairs": pairs, "headline": headline,
                       "reference": reference, "common_sample": common,
                       "audit_scope": "日期比较使用出生、死亡、失访及分析截止日期约束内的全部有效FOD；baseline模式只截断负担与聚类特征，不截断日期审计。"},
            "burden": {"window_label": (f"评估前 {args.lookback_years:g} 年" if args.mode == "baseline" else f"首次可见临床日期起 {args.landmark_years:g} 年"),
                       "scope": "+".join(s for s, l in levels.items() if l >= 1), "headline": "把起病早晚与后续记录密度、间隔和持续性分开看。",
                       "feature_definitions": [{"key": k, "label": l, "description": d} for k, l, d in [
                         ("n_days", "记录日", "不同有效日期数；同日不同代码不重复计数"), ("n_bouts", "记录簇", f"相邻记录间隔≤{args.episode_gap}天串联；不等同住院次数"),
                         ("active_fraction", "活跃时间段", "固定窗口内有记录的季度比例；无覆盖信息时仅作记录密度"), ("median_gap_days", "相邻间隔", "个人去重日期间隔的中位数；单次者缺失"),
                         ("event_rate", "年度记录率", "仅在完整来源观察期和逐次记录条件成立时计算"), ("aod", "Y.aod", "三个病例内百分位的等权均值：记录率、活跃季度、首末跨度；探索性")]],
                       "distributions": {"days": hist([f["n_days"] for f in withdays], [0, 1, 2, 5, 10, 20, 50], ["0（已确认观察）", "1", "2", "3–5", "6–10", "11–20", "21–50", ">50"]),
                           "bouts": hist([f["n_bouts"] for f in withdays], [0, 1, 2, 5, 10, 20], ["0", "1", "2", "3–5", "6–10", "11–20", ">20"]),
                           "gap": hist([f["median_gap_days"] for f in withdays], [30, 90, 365, 730, 1826], ["≤30天", "31–90天", "91天–1年", "1–2年", "2–5年", ">5年"]),
                           "score": hist(eligible_score, [20, 40, 60, 80], ["0–20", "20–40", "40–60", "60–80", "80–100"])},
                       "timeline": {"labels": labels, "rows": tl, "eligible_n": timeline_n, "candidate_n": len(withdays),
                                    "coverage_confirmed": bool(args.coverage),
                                    "scope": "完整来源覆盖的固定窗口" if args.coverage else "可见日期模式；来源覆盖期未确认",
                                    "excluded": [{"reason": key, "n": value} for key, value in sorted(timeline_excluded.items()) if value]},
                       "score": {"n": len(eligible_score), "method": score["method"], "scope": score.get("scope", "unavailable"), "status": score["status"]},
                       "qc": {"counts": [{"label": qc_labels.get(k, k), "n": n} for k, n in sorted(qc.items())], "sensitivity": [
                           {"label": "孤立极早日期者（主分析保留）", "value": len(suspect)},
                           {"label": "暂不采用孤立日期后，窗口记录日数中位变化", "value": med([f["n_days_without_isolate"] - f["n_days"] for f in suspect])},
                           {"label": "仅在建模特征中缩尾的值", "value": model.get("clipped_values", 0)}]}, "notes": notes},
            "clusters": clusters, "methods": [{"title": "观察分母", "text": "病例证据并集不等于共同有日期人数。死亡、失访或截点前未完成窗口者不补零，分群排除原因可查。"},
                         {"title": "可重复聚类", "text": "NumPy多起点k-means；训练样本1%/99%缩尾与中位数/IQR标准化；K=2–6或指定2–8；抽样silhouette、80%重抽样ARI及最小簇人数共同筛选，允许不分群。"},
                         {"title": "分析用途", "text": "模型在本次病例集拟合，仅作描述性表型探索。外部验证、GWAS和预测需要独立训练/验证设计；不能用未来轨迹作为baseline预测变量。"}]}


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
            src = "gp2" if (group == 3 or j % 3 == 1) else "icd10"
            (gp if src == "gp2" else hospital).append(d)
            es.append(dict(eid=eid, source=src, code="G3" if src == "gp2" else "I259", date=iso(d), date_raw=iso(d), date_status="valid", origin="SYNTHETIC", granularity="record" if src == "gp2" else "episode", match_local="TRUE", match_official="TRUE"))
        hs = min(hospital) if hospital else None
        gs = min(gp) if gp else None
        ref = first + (0 if i % 5 else 365)
        ps.append(dict(eid=eid, date_birth="1955-06-15", date_assess="2010-07-01", date_death="", date_lost="", icd10=iso(hs), gp=iso(gs), ref=iso(ref), min1=iso(first), min2=iso(ref), adu=round(float(rng.normal()), 4)))
        for src in ("hospital", "gp"):
            cov.append(dict(eid=eid, source=src, start="1995-01-01", end="2025-12-31"))
    write_tsv(work / "patients.tsv.gz", ["eid", "date_birth", "date_assess", "date_death", "date_lost", "icd10", "gp", "ref", "min1", "min2", "adu"], ps)
    write_tsv(work / "events.tsv.gz", ["eid", "source", "code", "date", "date_raw", "date_status", "origin", "granularity", "match_local", "match_official"], es)
    write_tsv(work / "coverage.tsv", ["eid", "source", "start", "end"], cov)
    md = {"Y": trait, "demo": True, "as_of": "2025-12-31", "raw_hospital_episode_dates_available": True, "raw_icd10_available": True, "raw_gp2_available": True}
    write_tsv(work / "metadata.tsv", ["key", "value"], [{"key": k, "value": json.dumps(v)} for k, v in md.items()])
    return str(work / "coverage.tsv")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--Y", required=True)
    p.add_argument("--coverage")
    p.add_argument("--as-of")
    p.add_argument("--mode", choices=["history", "baseline"], default="history")
    p.add_argument("--lookback-years", type=float, default=5)
    p.add_argument("--landmark-years", type=float, default=5)
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
    patients = read_patients(root / "patients.tsv.gz", asof, args.mode)
    selected, levels, qc, procedure, local = read_events(root / "events.tsv.gz", patients, meta)
    definitions = meta.get("definitions", {}).get("sources", {})
    args.partial_sources = [group for group, members in (("hospital", ("icd10", "icd9")), ("gp", ("gp2", "gp3")))
                            if levels[group] == 3 and any(truth(definitions.get(s, {}).get("local_defined")) and
                                                         not truth(meta.get("raw_" + s + "_available")) for s in members)]
    coverage = read_coverage(args.coverage)
    feats, used_sources, nbins, window = feature_rows(patients, selected, levels, coverage, args)
    score = score_burden(feats, window)
    cl, model = cluster(feats, used_sources, nbins, window, args)
    model.update(trait=args.Y, mode=args.mode, as_of=asof_text, landmark_years=args.landmark_years, lookback_years=args.lookback_years,
                 episode_gap=args.episode_gap, source_levels=levels, score_model=score, coverage_supplied=bool(args.coverage),
                 source_base_commit="51409ca789ab0d01011e1982816ebe8f11fc5ea2", numpy_version=np.__version__,
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
    public_fields = ["eid", "window_start", "window_end", "index_date", "n_days_total", "n_days", "n_bouts", "n_hospital_days", "n_gp_days", "span_years",
                     "median_gap_days", "max_gap_days", "active_fraction", "late_minus_early", "time_since_last_days", "hospital_fraction", "gp_fraction",
                     "hospital_observed_years", "gp_observed_years", "event_rate", "hospital_rate", "gp_rate", "recent_rate", "early_isolate", "n_days_without_isolate",
                     "horizon_complete", "coverage_complete", "record_scope", "hospital_mode", "gp_mode", "aod_reason", "cluster_reason"]
    output_rows = []
    for f in feats:
        p = patients[f["eid"]]
        row = {k: clean_json(f.get(k)) for k in public_fields}
        row.update({"fod_" + k: iso(p[k]) for k in SOURCES})
        row.update({args.Y + ".aod": f["aod"], args.Y + ".cluster": f["cluster"], args.Y + ".adu": f["adu"]})
        output_rows.append(row)
    public_fields += ["fod_" + k for k in SOURCES] + [args.Y + ".aod", args.Y + ".cluster"]
    if any(f["adu"] is not None for f in feats):
        public_fields.append(args.Y + ".adu")
    write_tsv(out / "aod.tsv.gz", public_fields, output_rows)
    profile_fields = ["id", "label", "n", "pct", "days_median", "bouts_median", "span_years_median", "active_fraction_median", "hospital_fraction_median"]
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
