#!/usr/bin/env python3
"""Build an aggregate-only UKB date report using the Python standard library.

Identifiers are used for joins in memory only, never embedded in the report.
Recorded dates describe available evidence, not incidence, independent episodes,
disease severity, or a causal change at recruitment.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
from datetime import date, datetime, timezone
from functools import lru_cache
import gzip
import html
import json
import math
from pathlib import Path
import re
import statistics

csv.field_size_limit(1 << 28)
SPECIAL_DATES = {"1900-01-01", "1901-01-01", "1902-02-02",
                 "1903-03-03", "1909-09-09", "2037-07-07"}
MISSING = {"", "NA", "N/A", "NaN", "nan", "NULL", "null", "."}
SOURCE_LABELS = {"min1": "min1", "min2": "min2", "ref": "ref",
                 "icd10": "住院 ICD10", "icd9": "住院 ICD9", "gp": "GP",
                 "srd": "自报", "opcs4": "表型 OPCS4", "death": "疾病死因"}
SOURCE_ORDER = tuple(SOURCE_LABELS)
LOCAL_SOURCES = {"icd10", "icd9", "gp", "gp2", "gp3", "srd", "opcs4"}
OFFICIAL_SOURCES = {"icd10", "icd9", "gp2", "gp3", "srd", "death"}
GROUP_FOR_SOURCE = {"icd10": "hospital", "icd9": "hospital",
                    "gp": "gp", "gp2": "gp", "gp3": "gp", "srd": "srd",
                    "opcs4": "opcs4", "opcs4_context": "opcs4", "death": "death"}
GROUP_LABELS = {"clinical": "住院 + GP", "hospital": "住院诊断", "gp": "GP",
                "srd": "自报", "opcs4": "OPCS4 操作", "death": "疾病死因"}
BINS = [{"lo": -10 + i / 2, "hi": -10 + (i + 1) / 2,
         "mid": -10 + (i + .5) / 2} for i in range(40)]
LAG_LIMITS = (-3653, -1827, -366, -31, -1, 0, 30, 365, 1826, 3652)
LAG_LABELS = ("≤−3653 天", "−3652 至 −1827 天", "−1826 至 −366 天",
              "−365 至 −31 天", "−30 至 −1 天", "同日", "+1 至 +30 天",
              "+31 至 +365 天", "+366 至 +1826 天", "+1827 至 +3652 天", "≥+3653 天")


@contextmanager
def table_rows(path):
    """Validate and stream rows, without materialising the event file."""
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        fields = next(reader, [])
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"Missing or duplicate table header: {path.name}")

        def rows():
            for number, row in enumerate(reader, 2):
                if len(row) != len(fields):
                    raise ValueError(f"{path.name}:{number}: {len(row)} cells, expected {len(fields)}")
                yield ["" if x in MISSING else x for x in row]
        yield fields, rows()


def read_table(path):
    with table_rows(path) as (fields, rows):
        return fields, list(rows)


def read_metadata(path):
    path = Path(path)
    if not path.exists():
        return {}
    if path.suffix == ".json":
        result = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(result, dict):
            raise ValueError("Metadata must be a JSON object")
        return result
    fields, rows = read_table(path)
    if not {"key", "value"} <= set(fields):
        raise ValueError("metadata.tsv needs key and value columns")
    ki, vi, result = fields.index("key"), fields.index("value"), {}
    for row in rows:
        key = row[ki]
        if key in result:
            raise ValueError(f"Duplicate metadata key: {key}")
        try:
            result[key] = json.loads(row[vi])
        except (ValueError, TypeError):
            result[key] = row[vi]
    return result


@lru_cache(maxsize=65536)
def ordinal(value):
    if not value or value in SPECIAL_DATES:
        return None
    try:
        parsed = date.fromisoformat(value)
    except (TypeError, ValueError):
        return None
    return parsed.toordinal() if parsed.isoformat() == value else None


def truth(value):
    return value is True or str(value).strip().lower() in {"true", "t", "1", "yes"}


def positive(value):
    try:
        return float(value) > 0
    except (TypeError, ValueError):
        return False


def pct(n, denominator):
    return round(100 * n / denominator, 4) if denominator else None


def median(values):
    return statistics.median(values) if values else None


def histogram(values, limits, labels):
    counts = [0] * len(labels)
    for value in values:
        index = next((i for i, upper in enumerate(limits) if value <= upper), len(limits))
        counts[index] += 1
    return [{"label": label, "n": n} for label, n in zip(labels, counts)]


def pair_metrics(a, b):
    if len(a) != len(b):
        raise ValueError("Paired date vectors have unequal lengths")
    diffs = [x - y for x, y in zip(a, b) if x is not None and y is not None]
    only_a = sum(x is not None and y is None for x, y in zip(a, b))
    only_b = sum(x is None and y is not None for x, y in zip(a, b))
    exact = sum(d == 0 for d in diffs)
    within30 = sum(abs(d) <= 30 for d in diffs)
    within365 = sum(abs(d) <= 365 for d in diffs)
    return {"both": len(diffs), "same_day": exact, "same_day_pct": pct(exact, len(diffs)),
            "a_earlier": sum(d < 0 for d in diffs), "b_earlier": sum(d > 0 for d in diffs),
            "a_only": only_a, "b_only": only_b,
            "neither_date": len(a) - len(diffs) - only_a - only_b,
            "median_days": median(diffs),
            "mean_days": round(statistics.mean(diffs), 2) if diffs else None,
            "min_days": min(diffs) if diffs else None, "max_days": max(diffs) if diffs else None,
            "within_30d": within30, "within_30d_pct": pct(within30, len(diffs)),
            "within_365d": within365, "within_365d_pct": pct(within365, len(diffs)),
            "histogram": histogram(diffs, LAG_LIMITS, LAG_LABELS)}


def reverse_pair(metric):
    reverse = dict(metric)
    for left, right in (("a_earlier", "b_earlier"), ("a_only", "b_only")):
        reverse[left], reverse[right] = metric[right], metric[left]
    for key in ("median_days", "mean_days"):
        reverse[key] = -metric[key] if metric[key] is not None else None
    reverse["min_days"] = -metric["max_days"] if metric["max_days"] is not None else None
    reverse["max_days"] = -metric["min_days"] if metric["min_days"] is not None else None
    reverse["histogram"] = [{"label": label, "n": item["n"]}
                            for label, item in zip(LAG_LABELS, reversed(metric["histogram"]))]
    return reverse


def coverage_summary(a_name, a, ref, a_evidence, ref_evidence):
    metric, total = pair_metrics(a, ref), len(a)
    labels = (("both", "两边都有有效日期"), ("a_only", f"仅 {a_name} 有有效日期"),
              ("b_only", "仅 ref 有有效日期"), ("neither_date", "两边均无有效日期"))
    reasons = Counter()
    reason_labels = {
        "a_only_ref_undated": f"仅 {a_name} 有日期 · ref 有证据但日期不可用",
        "a_only_no_ref": f"仅 {a_name} 有日期 · 当前输入未见 ref 证据",
        "ref_only_a_undated": f"仅 ref 有日期 · 有 {a_name} 匹配证据但日期不可用",
        "ref_only_no_a": f"仅 ref 有日期 · 当前输入未见 {a_name} 匹配证据",
        "neither_both_evidence": "两边均无日期 · 两边均有匹配证据",
        "neither_a_evidence": f"两边均无日期 · 仅有 {a_name} 匹配证据",
        "neither_ref_evidence": "两边均无日期 · 仅有 ref 证据",
        "neither_other_evidence": "两边均无日期 · 仅见其他定义/来源的证据",
    }
    for x, y, ae, re_ in zip(a, ref, a_evidence, ref_evidence):
        if x is not None and y is not None:
            continue
        if x is not None:
            key = "a_only_ref_undated" if re_ else "a_only_no_ref"
        elif y is not None:
            key = "ref_only_a_undated" if ae else "ref_only_no_a"
        elif ae and re_:
            key = "neither_both_evidence"
        elif ae:
            key = "neither_a_evidence"
        elif re_:
            key = "neither_ref_evidence"
        else:
            key = "neither_other_evidence"
        reasons[key] += 1
    result = {"a": a_name, "b": "ref", "total": total,
              "valid_a": metric["both"] + metric["a_only"],
              "valid_b": metric["both"] + metric["b_only"],
              "excluded": total - metric["both"], **metric,
              "four_way": [{"key": key, "label": label, "n": metric[key], "pct": pct(metric[key], total)}
                           for key, label in labels],
              "excluded_reasons": [{"key": key, "label": reason_labels[key], "n": reasons[key],
                                   "pct": pct(reasons[key], total)}
                                  for key in reason_labels if reasons[key]]}
    if sum(x["n"] for x in result["four_way"]) != total or sum(reasons.values()) != result["excluded"]:
        raise AssertionError("Coverage denominators do not balance")
    return result


def common_comparison(min1, min2, ref):
    counts = Counter()
    for a, b, r in zip(min1, min2, ref):
        if a is None or b is None or r is None:
            continue
        counts["n"] += 1
        old, new = a == r, b == r
        counts["min1_same_day"] += old
        counts["min2_same_day"] += new
        counts["both_same" if old and new else "lost" if old else "gained" if new else "neither_same"] += 1
    result = {key: counts[key] for key in ("n", "min1_same_day", "min2_same_day",
                                          "gained", "lost", "both_same", "neither_same")}
    result["min1_same_day_pct"] = pct(counts["min1_same_day"], counts["n"])
    result["min2_same_day_pct"] = pct(counts["min2_same_day"], counts["n"])
    result["delta_pct_points"] = pct(counts["gained"] - counts["lost"], counts["n"])
    return result


def merge_dates(*mappings):
    result = defaultdict(set)
    for mapping in mappings:
        for person, days in mapping.items():
            result[person].update(days)
    return dict(result)


def time_bin(delta_days):
    return math.floor((delta_days / 365.25 + 10) * 2)


def binned_series(key, retained, records, baseline, retained_available, record_available):
    def collect(mapping):
        counts, people, weights = [0] * 40, [0] * 40, [0.0] * 40
        n_people = n_days = below = above = baseline_day_people = 0
        below_weight = above_weight = 0.0
        for person, days in mapping.items():
            base = baseline[person]
            if base is None or not days:
                continue
            n_people += 1
            n_days += len(days)
            baseline_day_people += base in days
            seen = set()
            for day in days:
                i = time_bin(day - base)
                if i < 0:
                    below += 1
                    below_weight += 1 / len(days)
                elif i >= 40:
                    above += 1
                    above_weight += 1 / len(days)
                else:
                    counts[i] += 1
                    weights[i] += 1 / len(days)
                    seen.add(i)
            for i in seen:
                people[i] += 1
        factor = 100 / n_people if n_people else None
        return {"n_people": n_people, "n_days": n_days, "days": counts, "people": people,
                "weights": [round(w * factor, 6) for w in weights] if factor is not None else None,
                "below": below, "above": above,
                "below_weight": round(below_weight * factor, 6) if factor is not None else None,
                "above_weight": round(above_weight * factor, 6) if factor is not None else None,
                "baseline_day_people": baseline_day_people}
    all_, raw = collect(retained), collect(records)
    return {"key": key, "label": GROUP_LABELS[key],
            "retained_input_available": bool(retained_available), "record_input_available": bool(record_available),
            "n_people": len(retained), "n_record_people": len(records),
            "n_baseline_people": all_["n_people"], "n_baseline_days": all_["n_days"],
            "n_baseline_raw_people": raw["n_people"], "n_baseline_raw_days": raw["n_days"],
            "baseline_day_people": all_["baseline_day_people"], "baseline_day_raw_people": raw["baseline_day_people"],
            "retained_days": all_["days"], "raw_record_days": raw["days"],
            "retained_people": all_["people"], "raw_record_people": raw["people"],
            "patient_weighted_pct": all_["weights"], "raw_record_patient_weighted_pct": raw["weights"],
            "underflow_days": all_["below"], "overflow_days": all_["above"],
            "underflow_raw_days": raw["below"], "overflow_raw_days": raw["above"],
            "underflow_weight_pct": all_["below_weight"], "overflow_weight_pct": all_["above_weight"],
            "underflow_raw_weight_pct": raw["below_weight"], "overflow_raw_weight_pct": raw["above_weight"]}


def paired_windows(clinical, records, baseline, ends, record_available):
    # Equal windows in the same participants; as-of/no known death is NOT
    # evidence of complete source-specific medical-record coverage.
    eligible = [i for i, b in enumerate(baseline) if b is not None and ends[i] >= b + 365]
    output = []
    for key, label, mapping, omit, available in (
        ("all_dates", "全部保留日期 · 基线前后 365 天", clinical, 0, bool(clinical)),
        ("away30", "去掉基线前后 30 天", clinical, 30, bool(clinical)),
        ("records_only", "逐条记录 · 基线前后 365 天", records, 0, record_available),
        ("records_away30", "逐条记录 · 去掉基线前后 30 天", records, 30, record_available),
    ):
        before, after, on_baseline = [], [], 0
        for i in eligible:
            b, days = baseline[i], mapping.get(i, ())
            before.append(sum(-365 <= d - b < -omit for d in days))
            after.append(sum(omit < d - b <= 365 for d in days))
            on_baseline += b in days
        def stats(values):
            return {"people": sum(v > 0 for v in values), "days": sum(values), "median": median(values),
                    "mean": round(statistics.mean(values), 4) if values else None}
        output.append({"key": key, "label": label, "available": bool(available),
                       "eligible_n": len(eligible), "window_days": 365, "excluded_near_baseline_days": omit,
                       "eligible_with_history": sum(i in mapping for i in eligible),
                       "before": stats(before), "after": stats(after), "baseline_day_people": on_baseline,
                       "increased": sum(y > x for x, y in zip(before, after)),
                       "equal": sum(y == x for x, y in zip(before, after)),
                       "decreased": sum(y < x for x, y in zip(before, after))})
    return output


def trajectory_groups(clinical, baseline, n):
    labels = {"no_dates": "无住院 / GP 日期", "missing_baseline": "有日期 · 基线不可用",
              "pre_only": "仅基线前（可含当日）", "baseline_only": "仅基线当日",
              "post_only": "仅基线后（可含当日）", "both_sides": "基线前后均有日期"}
    totals, cells = Counter(), {key: [0] * 40 for key in labels}
    for person in range(n):
        days, b = clinical.get(person, ()), baseline[person]
        if not days:
            key = "no_dates"
        elif b is None:
            key = "missing_baseline"
        else:
            pre, post = min(days) < b, max(days) > b
            key = "both_sides" if pre and post else "pre_only" if pre else "post_only" if post else "baseline_only"
        totals[key] += 1
        if b is not None:
            for j in {time_bin(d - b) for d in days}:
                if 0 <= j < 40:
                    cells[key][j] += 1
    return [{"key": key, "label": label, "n": totals[key], "pct": pct(totals[key], n),
             "heatmap": [pct(c, totals[key]) or 0 for c in cells[key]]} for key, label in labels.items()]


def persistence_summary(clinical):
    counts, spans, gaps, repeat_spans = [], [], [], []
    grid = [[0] * 5 for _ in range(4)]
    for days in clinical.values():
        ordered = sorted(days)
        count, span = len(ordered), (ordered[-1] - ordered[0]) / 365.25
        counts.append(count)
        spans.append(span)
        if count > 1:
            repeat_spans.append(span)
            gaps.append(statistics.median(b - a for a, b in zip(ordered, ordered[1:])))
        x = 0 if count == 1 else 1 if count <= 3 else 2 if count <= 9 else 3 if count <= 29 else 4
        y = 0 if span < 1 else 1 if span < 5 else 2 if span < 10 else 3
        grid[y][x] += 1
    return {"n_people": len(counts), "n_repeated": len(repeat_spans), "median_n_dates": median(counts),
            "median_span_years": round(median(repeat_spans), 3) if repeat_spans else None,
            "median_span_all_years": round(median(spans), 3) if spans else None,
            "median_person_gap_days": median(gaps),
            "counts": histogram(counts, (1, 3, 9, 29), ("1", "2–3", "4–9", "10–29", "≥30")),
            "spans": [{"label": label, "n": sum(row)} for label, row in
                      zip(("<1 年", "1–<5 年", "5–<10 年", "≥10 年"), grid)],
            "gaps": [{"label": label, "n": sum(lo <= g < hi for g in gaps)}
                     for label, lo, hi in (("<30 天", 0, 30), ("30–<365 天", 30, 365),
                                          ("365–<1826 天", 365, 1826), ("≥1826 天", 1826, math.inf))],
            "gap_definition": "每位有至少两个日期者的相邻日期中位间隔；每人只贡献一次",
            "count_span": {"x_labels": ["1", "2–3", "4–9", "10–29", "≥30"],
                           "y_labels": ["<1 年", "1–<5 年", "5–<10 年", "≥10 年"], "counts": grid}}


def procedure_summary(procedures, clinical, records, baseline, metadata, trait, source_available):
    relative_baseline, relative_diagnosis = [], []
    for person, days in procedures.items():
        first = min(days)
        if baseline[person] is not None:
            relative_baseline.append(first - baseline[person])
        if clinical.get(person):
            relative_diagnosis.append(first - min(clinical[person]))
    context = truth(metadata.get("procedure_context_defined"))
    defined = context or truth(metadata.get("procedure_definition_available")) or bool(procedures)
    label = "冠脉搭桥或治疗性冠脉介入" if context and trait == "cvd_cad" else "本表型的 OPCS4 操作"
    return {"available": bool(defined and (source_available or procedures)), "label": label,
            "n_people": len(procedures), "n_record_people": len(records),
            "n_baseline_people": len(relative_baseline), "n_with_diagnosis": len(relative_diagnosis),
            "n_retained_days": sum(map(len, procedures.values())), "n_record_days": sum(map(len, records.values())),
            "timing_baseline": histogram(relative_baseline, (-366, -31, 30, 365),
                                         ("基线前 >1 年", "基线前 31–365 天", "基线 ±30 天",
                                          "基线后 31–365 天", "基线后 >1 年")),
            "timing_diagnosis": histogram(relative_diagnosis, (-366, -1, 0, 365),
                                          ("诊断前 >1 年", "诊断前 1–365 天", "同日",
                                           "诊断后 1–365 天", "诊断后 >1 年")),
            "definition_note": (
                "每人使用最早可见操作日期；诊断锚点为最早可见住院 / GP 日期。"
                "逐条手术记录使用 opdate；OPCS 汇总日期可能是相关住院 episode 起始 / 入院日期。"
                "操作记录不能直接解释为疾病严重度或 ICU 经历。"
                + ("冠脉操作仅作病例内背景，不改变病例并集或 min1。" if context else ""))}


def make_summary(patient_fields, patients, event_fields, events, metadata):
    pi, ei = {name: i for i, name in enumerate(patient_fields)}, {name: i for i, name in enumerate(event_fields)}
    if len(pi) != len(patient_fields) or len(ei) != len(event_fields):
        raise ValueError("Duplicate input columns")
    if "eid" not in pi or not {"eid", "source", "date"} <= set(ei):
        raise ValueError("patients needs eid; events needs eid, source and date")
    if any(len(row) != len(patient_fields) for row in patients):
        raise ValueError("Malformed patient row")
    identifiers = [str(row[pi["eid"]]) for row in patients]
    if any(not x or x in MISSING for x in identifiers) or len(set(identifiers)) != len(identifiers):
        raise ValueError("Missing or duplicate patient identifiers")
    join, n = {value: i for i, value in enumerate(identifiers)}, len(patients)
    trait = str(metadata.get("Y", metadata.get("trait", "disease")))
    if not re.fullmatch(r"[A-Za-z0-9_]+", trait):
        raise ValueError("Trait name must contain only letters, digits, or underscores")
    asof = ordinal(metadata.get("as_of") or date.today().isoformat())
    if asof is None:
        raise ValueError("as_of must be an ordinary YYYY-MM-DD date")
    def value(row, key, fallback=""):
        return row[pi[key]] if key in pi else fallback
    def context_date(row, names):
        key = next((key for key in names if key in pi), None)
        d = ordinal(value(row, key)) if key else None
        return d if d is not None and d <= asof else None
    birth = [context_date(row, ("date_birth", "birth_date")) for row in patients]
    deaths = [context_date(row, ("date_death",)) for row in patients]
    losses = [context_date(row, ("date_lost",)) for row in patients]
    ends = [min(d for d in (asof, deaths[i], losses[i]) if d is not None) for i in range(n)]
    baseline, invalid_baseline = [], 0
    for i, row in enumerate(patients):
        raw = value(row, "date_assess", value(row, "date_attend"))
        d = ordinal(raw)
        if d is not None and (d > ends[i] or (birth[i] is not None and d < birth[i])):
            d = None
        if raw and raw not in MISSING and d is None:
            invalid_baseline += 1
        baseline.append(d)
    def dates_for(*keys):
        key = next((key for key in keys if key in pi), None)
        if key is None:
            return [None] * n
        dates = [ordinal(row[pi[key]]) for row in patients]
        return [d if d is not None and d <= asof and (birth[i] is None or d >= birth[i])
                else None for i, d in enumerate(dates)]
    dates = {name: dates_for(name) for name in (*SOURCE_ORDER, "gp2", "gp3")}
    dates["min1"] = dates_for("min1", "me")
    dates["min2"] = dates_for("min2", "reconstructed_ref")
    if "gp" not in pi:
        dates["gp"] = [min(d for d in (a, b) if d is not None) if a is not None or b is not None else None
                       for a, b in zip(dates["gp2"], dates["gp3"])]
    evidence = {key: [d is not None for d in dates[key]] for key in ("min1", "min2", "ref")}
    for key, column in {"min1": "has_local_evidence", "min2": "has_official_map_evidence",
                        "ref": "has_ref_evidence"}.items():
        for i, row in enumerate(patients):
            evidence[key][i] |= truth(value(row, column))
    for i, row in enumerate(patients):
        evidence["min1"][i] |= any(positive(value(row, "cnt_" + source)) for source in LOCAL_SOURCES)
        evidence["ref"][i] |= bool(value(row, "ref_raw")) or positive(value(row, "cnt_ref"))

    # Retain original-source sets until the FOD fallback decision. Otherwise
    # ICD10 AOD could incorrectly suppress a person's only available ICD9 date.
    retained_source = {key: defaultdict(set) for key in GROUP_FOR_SOURCE}
    record_source = {key: defaultdict(set) for key in GROUP_FOR_SOURCE}
    fallback_source = {key: defaultdict(set) for key in GROUP_FOR_SOURCE}
    summary_people = {key: set() for key in GROUP_LABELS if key != "clinical"}
    event_sources, source_evidence = Counter(), Counter()
    n_events = invalid_rows = after_end = excluded_nonlocal = 0
    for row in events:
        if len(row) != len(event_fields):
            raise ValueError("Malformed event row")
        if str(row[ei["eid"]]) not in join:
            raise ValueError("An event refers to a participant absent from the supplied case union")
        person, source = join[str(row[ei["eid"]])], row[ei["source"]]
        n_events += 1
        event_sources[source] += 1
        local = truth(row[ei["match_local"]]) if "match_local" in ei else False
        official = truth(row[ei["match_official"]]) if "match_official" in ei else False
        if source == "ref":
            evidence["ref"][person] = True
        if local and source in LOCAL_SOURCES:
            evidence["min1"][person] = True
            source_evidence["gp" if source in ("gp", "gp2", "gp3") else source] += 1
        if official and source in OFFICIAL_SOURCES:
            evidence["min2"][person] = True
        group = GROUP_FOR_SOURCE.get(source)
        if group is None:
            continue
        if not local:
            excluded_nonlocal += 1
            continue
        day = ordinal(row[ei["date"]])
        status = row[ei["date_status"]] if "date_status" in ei else ""
        if (day is None or status not in ("", "valid") or day > asof
                or (birth[person] is not None and day < birth[person])):
            invalid_rows += 1
            continue
        if day > ends[person]:
            after_end += 1
            continue
        granularity = row[ei["granularity"]] if "granularity" in ei else ""
        if granularity in ("official_ref", "official_ref_cache"):
            continue
        if granularity == "fod_only":
            fallback_source[source][person].add(day)
        else:
            retained_source[source][person].add(day)
            if granularity in ("record", "episode") or (source == "death" and granularity == "registration"):
                record_source[source][person].add(day)
            else:
                summary_people[group].add(person)
    for source in GROUP_FOR_SOURCE:
        if source in dates:
            for i, day in enumerate(dates[source]):
                if day is not None and day <= ends[i]:
                    fallback_source[source][i].add(day)
    nonfod_people = {key: set() for key in summary_people}
    supplemental_people = {key: set() for key in summary_people}
    for source, group in GROUP_FOR_SOURCE.items():
        nonfod_people[group].update(retained_source[source])
        for person, days in fallback_source[source].items():
            if person not in retained_source[source]:
                retained_source[source][person].update(days)
                supplemental_people[group].add(person)
    retained, records = {}, {}
    for key in summary_people:
        keys = [source for source, group in GROUP_FOR_SOURCE.items() if group == key]
        retained[key] = merge_dates(*(retained_source[source] for source in keys))
        records[key] = merge_dates(*(record_source[source] for source in keys))
    del retained_source, record_source, fallback_source
    fallback_people = {key: supplemental_people[key] - nonfod_people[key] for key in retained}
    clinical = merge_dates(retained["hospital"], retained["gp"])
    clinical_records = merge_dates(records["hospital"], records["gp"])

    sources = [key for key in SOURCE_ORDER
               if key in ("min1", "min2", "ref") or any(d is not None for d in dates[key]) or source_evidence[key]]
    pairs = {key: {} for key in sources}
    for j, a in enumerate(sources):
        for b in sources[j:]:
            metric = pair_metrics(dates[a], dates[b])
            pairs[a][b] = metric
            if a != b:
                pairs[b][a] = reverse_pair(metric)
    coverage = {key: coverage_summary(key, dates[key], dates["ref"], evidence[key], evidence["ref"])
                for key in ("min1", "min2")}
    def meta_true(*names):
        return any(truth(metadata.get(name)) for name in names)
    procedure_defined = meta_true("procedure_context_defined", "procedure_definition_available") or bool(retained["opcs4"])
    raw_available = {"hospital": meta_true("raw_hospital_episode_dates_available"),
                     "gp": meta_true("raw_gp2_available", "raw_gp3_available"),
                     "opcs4": procedure_defined and meta_true("raw_opcs4_available", "raw_hospital_operation_dates_available"),
                     "srd": False, "death": bool(records["death"])}
    retained_available = {key: bool(retained[key] or event_sources[key]) for key in retained}
    saved = metadata.get("saved_sources", [])
    if isinstance(saved, list):
        for item in saved:
            if isinstance(item, dict) and item.get("status") == "saved_fod":
                group = GROUP_FOR_SOURCE.get(item.get("source"))
                if group:
                    retained_available[group] = True
    for key in raw_available:
        raw_available[key] |= bool(records[key])
        retained_available[key] |= raw_available[key]
    retained_available["hospital"] |= bool(event_sources["icd10"] or event_sources["icd9"])
    retained_available["gp"] |= bool(event_sources["gp2"] or event_sources["gp3"])
    retained_available["opcs4"] |= bool(event_sources["opcs4_context"])
    retained_available["opcs4"] |= procedure_defined and meta_true("procedure_summary_available")
    all_retained, all_records = {"clinical": clinical, **retained}, {"clinical": clinical_records, **records}
    retained_available["clinical"] = retained_available["hospital"] or retained_available["gp"]
    raw_available["clinical"] = raw_available["hospital"] or raw_available["gp"]
    series = [binned_series(key, all_retained[key], all_records[key], baseline,
                            retained_available[key], raw_available[key]) for key in GROUP_LABELS]
    provenance = [{"source": key, "label": GROUP_LABELS[key], "people_with_dates": len(retained[key]),
                   "with_record_dates": len(records[key]), "with_summary_dates": len(summary_people[key]),
                   "fod_fallback_only": len(fallback_people[key]),
                   "with_fod_supplement": len(supplemental_people[key]),
                   "distinct_days": sum(map(len, retained[key].values())),
                   "record_days": sum(map(len, records[key].values()))}
                  for key in retained if retained_available[key]]
    groups = trajectory_groups(clinical, baseline, n)
    if sum(item["n"] for item in groups) != n:
        raise AssertionError("Trajectory groups do not partition the case union")
    alerts = []
    missing_fields = metadata.get("reference_fields_missing", [])
    if isinstance(missing_fields, str):
        missing_fields = [missing_fields] if missing_fields else []
    if isinstance(missing_fields, list) and missing_fields:
        alerts.append(f"ref 输入缺少 {len(missing_fields)} 个预期字段；当前一致率仅基于已提供字段。")
    mode = metadata.get("reference_input_mode")
    if mode == "saved_only":
        alerts.append("ref 仅有既存缓存；更新定义后需重新提取原始字段并生成 fod_ref，才能核对新范围。")
    elif mode == "unavailable":
        alerts.append("当前输入没有可用 ref 字段或缓存，不能计算参考日期一致率。")
    if trait == "cvd_cad" and metadata.get("reference_scope_warning"):
        alerts.append("CAD 的本地定义与 ref 范围不一致；需包含 I20–I25，并重新生成参考日期后核对。")
    if invalid_baseline:
        alerts.append(f"{invalid_baseline:,} 人的基线日期不符合出生、死亡、失访或截至日期约束；保留于总人数，排除时间对齐。")
    result = {
        "schema_version": 2, "Y": trait, "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "as_of": date.fromordinal(asof).isoformat(), "is_demo": meta_true("is_demo", "demo"),
        "n_patients": n, "n_events": n_events,
        "population_definition": "任一本地、ref 或公开映射定义有证据者的并集；不是单一病例定义的确诊人数",
        "sources": sources,
        "source_stats": [{"source": key, "label": SOURCE_LABELS[key],
                          "dated_patients": sum(d is not None for d in dates[key]),
                          "pct_cases": pct(sum(d is not None for d in dates[key]), n)} for key in sources],
        "pairs": pairs, "coverage": coverage,
        "common_comparison": common_comparison(dates["min1"], dates["min2"], dates["ref"]),
        "trajectory": {
            "n_baseline": sum(b is not None for b in baseline), "n_clinical_people": len(clinical),
            "n_baseline_clinical": sum(baseline[i] is not None for i in clinical),
            "n_raw_clinical_people": len(clinical_records), "bins": BINS, "series": series,
            "provenance": provenance,
            "paired_windows": paired_windows(clinical, clinical_records, baseline, ends, raw_available["clinical"]),
            "groups": groups, "persistence": persistence_summary(clinical),
            "unit": "同一人、同一来源组、同一天去重；住院 + GP 再跨来源去重",
            "density_definition": "每人的全部可用日期合计权重为 1，再在有可用基线和日期者中取平均；窗口外权重另列",
            "paired_window_definition": "同一批有可用基线、且截至日/已知死亡/失访允许观察至基线后365天的人；来源覆盖仍未知",
        },
        "procedures": procedure_summary(retained["opcs4"], clinical, records["opcs4"], baseline,
                                        metadata, trait, retained_available["opcs4"]),
        "quality": {
            "invalid_baseline": invalid_baseline, "missing_or_invalid_baseline": sum(b is None for b in baseline),
            "invalid_or_undated_event_rows": invalid_rows, "after_known_end_rows": after_end,
            "nonlocal_event_rows_excluded_from_trajectory": excluded_nonlocal,
            "clinical_fod_fallback_people": len(set(clinical) - set(clinical_records)
                                                - summary_people["hospital"] - summary_people["gp"]),
            "raw_clinical_input_available": raw_available["clinical"],
        },
        "data_alerts": alerts,
        "notes": [
            "min1 为本地 ICD10 / ICD9 / OPCS4 / 自报 / GP 的最小日期；min2 为公开映射下的候选日期，并非独立验证。",
            "同日、30天和365天一致率的分母均为两边有有效日期者；min1 与 min2 的共同样本比较另外列出。",
            "AOD 保留了当前输入可见的日期。住院汇总常为每个代码的首次日期；逐条记录也按人-日去重，均不能直接当作独立发作次数或严重度。",
            "FOD 仅在该人该原始来源无 AOD / 逐条记录日期时作补充；这种证据不能评估复发。",
            "基线前后图描述记录分布；记录密度变化可能来自采集与覆盖差异，不能据此断定入组接受治疗或健康轨迹的因果转折。",
            "没有日期不等于没有疾病。未具备个体、来源特异的观察起止区间，故不计算人年发生率；自报日期与临床记录分别显示。",
        ],
    }
    validate_public_summary(result)
    return result


def validate_public_summary(summary):
    """Forbid raw-table keys or arbitrary metadata in the public schema."""
    forbidden = {"eid", "eids", "identifiers", "patients", "events", "patient_fields",
                 "event_fields", "date_raw", "origin", "metadata", "r_inputs", "ref_raw"}
    def visit(node):
        if isinstance(node, dict):
            if forbidden.intersection(node):
                raise ValueError("Nonaggregate field in public report")
            for value in node.values():
                visit(value)
        elif isinstance(node, list):
            for value in node:
                visit(value)
    visit(summary)
    json.dumps(summary, ensure_ascii=False, allow_nan=False)


def build(input_dir, output=None, trait=None, metadata_path=None):
    input_dir = Path(input_dir)
    fields, patients = read_table(input_dir / "patients.tsv.gz")
    metadata_path = Path(metadata_path) if metadata_path else input_dir / "metadata.tsv"
    metadata = read_metadata(metadata_path)
    if trait:
        metadata["Y"] = trait
    with table_rows(input_dir / "events.tsv.gz") as (event_fields, event_rows):
        summary = make_summary(fields, patients, event_fields, event_rows, metadata)
    packed = json.dumps(summary, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    packed = (packed.replace("&", r"\u0026").replace("<", r"\u003c").replace(">", r"\u003e")
              .replace("\u2028", r"\u2028").replace("\u2029", r"\u2029"))
    template = Path(__file__).with_name("date_check.html").read_text(encoding="utf-8")
    if template.count("__AGGREGATE_JSON__") != 1:
        raise ValueError("date_check.html must contain exactly one aggregate payload placeholder")
    document = template.replace("__TITLE__", html.escape(summary["Y"])).replace("__AGGREGATE_JSON__", packed)
    output = Path(output) if output else input_dir / "date-check.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")
    output.with_name("summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Report: {output}\nEvidence-union participants: {summary['n_patients']:,}; "
          f"input rows: {summary['n_events']:,}; aggregate output only")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--Y")
    parser.add_argument("--metadata", type=Path)
    args = parser.parse_args()
    try:
        build(args.input_dir, args.output, args.Y, args.metadata)
    except (ValueError, OSError) as exc:
        parser.exit(2, f"date_check_report: {exc}\n")


if __name__ == "__main__":
    main()

