#!/usr/bin/env python3
"""Cross-cohort concordance/sensitivity for the union of regional O/E highlights.

The complete region tables supply both cohorts' estimates, including real zero
events and low O/E values. Missing rows or unsupported estimates stay missing.
Matching thresholds do not imply matching called-site sets or cohort independence.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import math
import textwrap
from pathlib import Path

from research_common import atomic_json, chrom_key, chrom_name, html_report, html_table, identity, number, open_text, write_tsv

KEY_FIELDS = ("region_kind", "chrom", "start0", "end0", "window_bp", "annotation", "outcome", "population", "sensitivity")
VALUE_FIELDS = ("observed", "expected", "oe", "ci_low", "ci_high", "p", "fdr")
SUPPORT_FIELDS = ("n_eligible_supported", "n_unsupported", "observed_total", "prior_only_fraction")
NA = {"", ".", "NA", "N/A", "None", "nan", "NaN"}


def _canonical(value):
    """Order-insensitive lists of definitions; dictionary order is immaterial."""
    if isinstance(value, dict): return {str(k): _canonical(v) for k, v in sorted(value.items())}
    if isinstance(value, list): return sorted((_canonical(v) for v in value), key=lambda x: json.dumps(x, sort_keys=True))
    return value


def _thresholds(settings):
    values = settings.get("effective_thresholds")
    if not isinstance(values, list) or not values: return None
    result = []
    for supplied in values:
        if not isinstance(supplied, dict): return None
        kind = str(supplied.get("kind", "")).upper()
        fields = ("min_ac", "min_an") if kind == "AC" else ("min_af",) if kind == "AF" else ()
        if not fields: return None
        normalized = dict(supplied, kind=kind)
        for field in fields:
            value = number(supplied.get(field))
            if value is None or value < 0: return None
            normalized[field] = value
        result.append(normalized)
    return _canonical(result)


def _load_input(directory, side):
    root = Path(directory).expanduser().resolve()
    manifest_path = root / "qc/region_highlights.json"
    regions_path = root / "results/regions.tsv.gz"
    source = dict(side=side, root=str(root), manifest_path=str(manifest_path), regions_path=str(regions_path),
                  status="missing_input", target=root.name or side, threshold=3., manifest={}, reason="")
    missing = [str(p) for p in (manifest_path, regions_path) if not p.is_file()]
    if missing:
        source["reason"] = "Missing " + ", ".join(missing)
        return source
    manifest = json.loads(manifest_path.read_text())
    if not isinstance(manifest, dict): raise ValueError("Expected JSON object: " + str(manifest_path))
    source.update(manifest=manifest, target=str(manifest.get("target", root.name or side)))
    if manifest.get("status") != "complete":
        source.update(status="incomplete_input", reason="Region analysis status is " + str(manifest.get("status", "unknown")))
        return source
    threshold = number(manifest.get("threshold", manifest.get("highlight_fold", 3)))
    if threshold is None or threshold <= 0: raise ValueError("Invalid highlight threshold: " + str(manifest_path))
    source.update(status="available", threshold=threshold,
                  inputs={"manifest": identity(manifest_path), "regions": identity(regions_path)})
    return source


def _comparability(left, right):
    flags = []; differences = []
    settings = [s["manifest"].get("scientific_settings", {}) for s in (left, right)]
    settings = [s if isinstance(s, dict) else {} for s in settings]
    thresholds = [_thresholds(s) for s in settings]
    thresholds_match = all(x is not None for x in thresholds) and thresholds[0] == thresholds[1]
    if not thresholds_match: flags.append("not_comparable_thresholds")
    for field in ("mask_inputs", "context", "ref_definition"):
        if any(field not in s for s in settings) or _canonical(settings[0].get(field)) != _canonical(settings[1].get(field)):
            differences.append(field)
    for field in ("min_background_sites", "observation_rule"):
        if any(field in s for s in settings) and _canonical(settings[0].get(field)) != _canonical(settings[1].get(field)):
            differences.append(field)
    for field in ("ref", "build"):
        values = [s["manifest"].get(field, s["manifest"].get("genome_build") if field == "build" else None) for s in (left, right)]
        if any(v is None for v in values) or values[0] != values[1]: differences.append(field)
    if any(s["status"] != "available" for s in (left, right)): differences.append("input_availability")
    if differences: flags.append("not_comparable_definition")
    if not flags: flags.append("threshold_matched_support_may_differ")
    return dict(status=";".join(flags), flags=flags, thresholds_match=thresholds_match,
                definitions_match=not differences, definition_differences=differences,
                left_effective_thresholds=thresholds[0], right_effective_thresholds=thresholds[1],
                interpretation="Cross-cohort concordance/sensitivity. Called-site sets, allele detection and population coverage may differ; cohort independence is not assumed.")


def _key(row):
    absent = set(KEY_FIELDS) - row.keys()
    if absent: raise ValueError("Region table lacks exact-join columns: " + ",".join(sorted(absent)))
    values = {k: str(row[k]).strip() for k in KEY_FIELDS}
    for field in ("start0", "end0", "window_bp"):
        text = values[field]
        if field == "window_bp" and text in NA:
            values[field] = None
            continue
        try: values[field] = int(text)
        except (TypeError, ValueError) as exc: raise ValueError("Invalid integer region coordinate: " + field) from exc
    if values["start0"] < 0 or values["end0"] <= values["start0"]:
        raise ValueError("Regions require BED0 half-open coordinates with start0 < end0")
    if values["window_bp"] is not None and values["window_bp"] < 0: raise ValueError("Negative window_bp")
    values["chrom"] = chrom_name(values["chrom"])
    if values["annotation"] in NA: values["annotation"] = ""
    return tuple(values[k] for k in KEY_FIELDS)


def _estimate(row):
    values = {k: number(row.get(k, row.get("q") if k == "fdr" else row.get("OE") if k == "oe" else None)) for k in VALUE_FIELDS}
    support = {k: number(row.get(k)) for k in SUPPORT_FIELDS}
    valid = (values["oe"] is not None and values["oe"] >= 0 and values["observed"] is not None and
             values["observed"] >= 0 and values["expected"] is not None and values["expected"] > 0 and
             (support["n_eligible_supported"] is None or support["n_eligible_supported"] > 0))
    # A zero numerator with E>0 is a measured O/E=0, even if an inferential CI
    # was unavailable. Low count/precision flags do not remove a raw highlight.
    if not valid: values = {k: None for k in VALUE_FIELDS}
    token = str(row.get("passes_fold", "")).strip().lower()
    passed = True if token in {"true", "1", "yes"} else False if token in {"false", "0", "no"} else None
    return dict(status="measured" if valid else "unsupported", **values, **support,
                region_id=row.get("region_id", ""), input_passes_fold=passed,
                input_support=row.get("support_status", row.get("inference", row.get("status", ""))))


def _rows(source):
    if source["status"] != "available": return
    with open_text(source["regions_path"]) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not set(KEY_FIELDS) <= set(reader.fieldnames or []):
            raise ValueError("Region table lacks exact-join columns: " + source["regions_path"])
        for row in reader:
            yield _key(row), _estimate(row)


def _is_highlight(estimate, threshold):
    value = estimate["oe"]
    if estimate["status"] != "measured": return False
    # Preserve the producer's decision on unrounded O/E when the TSV has a
    # rounded printed value near the cutoff. This flag uses that input's own
    # manifest threshold; it is not a p/FDR or minimum-count filter.
    if estimate.get("input_passes_fold") is not None: return estimate["input_passes_fold"]
    return value >= threshold or math.isclose(value, threshold, rel_tol=1e-12, abs_tol=0.)


def _sort_key(key):
    kind, chrom, start, end, width, annotation, outcome, population, sensitivity = key
    return kind, chrom_key(chrom), start, end, width if width is not None else -1, annotation, outcome, population, sensitivity


def _missing_estimate(source):
    return dict(status=source["status"] if source["status"] != "available" else "missing_region",
                **{k: None for k in (*VALUE_FIELDS, *SUPPORT_FIELDS)}, region_id="", input_passes_fold=None, input_support="")


def _plot(out, rows, left, right, comparison):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from collections import defaultdict
    pairs = [r for r in rows if r["left_status"] == r["right_status"] == "measured"]
    fig, ax = plt.subplots(figsize=(8, 7), layout="constrained")
    groups = defaultdict(list)
    for r in pairs: groups[(r["outcome"], r["sensitivity"])].append(r)
    colors = plt.get_cmap("tab10")
    for i, (key, rr) in enumerate(sorted(groups.items())):
        ax.scatter([r["left_oe"] for r in rr], [r["right_oe"] for r in rr], s=22, alpha=.6,
                   color=colors(i % 10), label=" · ".join(key), rasterized=True, clip_on=False, zorder=3)
    maximum = max([left["threshold"], right["threshold"], 3., *[max(r["left_oe"], r["right_oe"]) for r in pairs]]) * 1.12
    ax.plot([0, maximum], [0, maximum], color="#aab4bb", lw=1, label="Equal O/E")
    ax.axvline(left["threshold"], color="#805c45", ls="--", lw=1.2)
    ax.axhline(right["threshold"], color="#805c45", ls="--", lw=1.2)
    ax.set_xscale("symlog", linthresh=1); ax.set_yscale("symlog", linthresh=1)
    ax.set(xlim=(0, maximum), ylim=(0, maximum), xlabel=left["target"] + " regional O/E", ylabel=right["target"] + " regional O/E")
    ticks=sorted({0., .5, 1., left["threshold"], right["threshold"],
                  *[10.**i for i in range(1, max(1, int(math.log10(maximum))+1))]})
    ticks=[x for x in ticks if x<=maximum]
    ax.set_xticks(ticks,[f"{x:g}" for x in ticks]);ax.set_yticks(ticks,[f"{x:g}" for x in ticks])
    ax.set_aspect("equal", adjustable="box"); ax.grid(alpha=.18)
    if not pairs: ax.text(.5, .5, "No jointly measured regions in the highlight union\nMissing or unsupported estimates are not zero", ha="center", va="center", transform=ax.transAxes)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    status = "Matched allele thresholds; called-site support may differ" if comparison["status"] == "threshold_matched_support_may_differ" else "Descriptive side-by-side values: " + "; ".join(s.replace("not_comparable_", "different/unknown ").replace("_", " ") for s in comparison["flags"])
    fig.suptitle("Cross-cohort region concordance / sensitivity\n" + textwrap.fill(status, width=83), fontsize=12)
    ax.set_title(f"Either-side highlight union; {len(pairs):,} jointly measured regions\nDashed O/E cutoffs: {left['threshold']:g}× / {right['threshold']:g}×; no significance filter", fontsize=10)
    fig.savefig(out / "figures" / ("concordance.png"), dpi=180)
    plt.close(fig)
    return len(pairs)


def compare_regions(left_dir, right_dir, out_dir):
    """Write an exact-coordinate, evidence-aware comparison; return its manifest.

    Missing inputs produce an explicit partial report. Malformed or conflicting
    duplicate rows raise ValueError; callers may keep valid cohort reports.
    """
    out = Path(out_dir).expanduser().resolve()
    for folder in ("results", "figures", "reports", "qc"): (out / folder).mkdir(parents=True, exist_ok=True)
    left, right = _load_input(left_dir, "left"), _load_input(right_dir, "right")
    comparison = _comparability(left, right)
    union = set(); examined = {}
    # Only the highlight union remains in memory; the full tables are streamed.
    for source in (left, right):
        examined[source["side"]] = 0
        for key, estimate in _rows(source):
            examined[source["side"]] += 1
            if _is_highlight(estimate, source["threshold"]): union.add(key)
    selected = {"left": {}, "right": {}}; duplicate_counts = {"left": 0, "right": 0}
    for source in (left, right):
        destination = selected[source["side"]]
        for key, estimate in _rows(source):
            if key not in union: continue
            if key in destination:
                prior = destination[key]
                if any(prior[k] != estimate[k] for k in ("status", "input_passes_fold", *VALUE_FIELDS, *SUPPORT_FIELDS)):
                    raise ValueError("Ambiguous conflicting region rows for exact key " + repr(key))
                duplicate_counts[source["side"]] += 1
            else: destination[key] = estimate
        if source["status"] == "available":
            for kind, path in (("manifest", source["manifest_path"]), ("regions", source["regions_path"])):
                if identity(path) != source["inputs"][kind]: raise ValueError("Input changed while comparing: " + path)
    results = []
    for key in sorted(union, key=_sort_key):
        row = dict(zip(KEY_FIELDS, key)); row["comparison_status"] = comparison["status"]
        for source in (left, right):
            side = source["side"]; estimate = selected[side].get(key, _missing_estimate(source))
            row.update({side + "_" + k: v for k, v in estimate.items()})
            row[side + "_highlight"] = int(_is_highlight(estimate, source["threshold"])) if estimate["status"] == "measured" else None
            row[side + "_target"] = source["target"]
        results.append(row)
    fields = [*KEY_FIELDS, "left_oe", "right_oe", "left_highlight", "right_highlight", "left_status", "right_status", "comparison_status"]
    fields += [side + "_" + k for k in ("target", "region_id", "observed", "expected", "ci_low", "ci_high", "p", "fdr", *SUPPORT_FIELDS, "input_passes_fold", "input_support") for side in ("left", "right")]
    table = out / "results/concordance.tsv.gz"
    write_tsv(table, fields, results)
    n_pairs = _plot(out, results, left, right, comparison)
    counts = dict(union_regions=len(results), jointly_measured_regions=n_pairs,
        left_highlight_regions=sum(r["left_highlight"] == 1 for r in results),
        right_highlight_regions=sum(r["right_highlight"] == 1 for r in results),
        both_highlight_regions=sum(r["left_highlight"] == r["right_highlight"] == 1 for r in results),
        left_missing_or_unsupported=sum(r["left_status"] != "measured" for r in results),
        right_missing_or_unsupported=sum(r["right_status"] != "measured" for r in results),
        left_rows_read=examined["left"], right_rows_read=examined["right"], collapsed_identical_duplicates=duplicate_counts)
    report = out / "reports/concordance.html"
    inputs = {s["side"]: {k: v for k, v in s.items() if k != "manifest"} for s in (left, right)}
    outputs = dict(table=str(table), png=str(out / "figures/concordance.png"), report=str(report))
    manifest = dict(schema="allele3-cross-cohort-concordance-1", status="complete" if all(s["status"] == "available" for s in (left, right)) else "complete_with_missing_input",
        comparability=comparison["status"], comparability_details=comparison, inputs=inputs, counts=counts,
        join_key=list(KEY_FIELDS), selection="Union at either cohort's own point O/E cutoff; producer passes_fold preserves its unrounded cutoff decision when available; no p/FDR filter",
        scatter="All jointly measured pairs in that union, including true zero events and low-O/E counterparts; not restricted to both highlights or both significant",
        missing="Missing input, missing region and unsupported point estimates remain NA, never zero",
        outputs=outputs, report=str(report), output_dir=str(out))
    atomic_json(out / "qc/concordance.json", manifest)
    labels = left["target"] + " / " + right["target"]
    input_text = "".join("<li>" + html.escape(s["side"] + ": " + s["status"] + (" — " + s["reason"] if s["reason"] else "")) + "</li>" for s in (left, right))
    details = html.escape(json.dumps(comparison, indent=2))
    preview = 500
    html_report(report, "Cross-cohort region concordance / sensitivity — " + labels, [
        ("Definition and support", "<p><strong>" + html.escape(comparison["status"]) + "</strong></p><p>" + html.escape(comparison["interpretation"]) + "</p><ul>" + input_text + "</ul><pre>" + details + "</pre>"),
        ("Measured region pairs", '<img src="../figures/concordance.png" alt="Cross-cohort regional O/E scatter"><p>Every jointly measured pair in the either-side highlight union is plotted, regardless of significance or whether the other cohort reaches its cutoff. Measured means point O/E is available; it does not certify adequate counts or precision. The complete table retains original support labels and intervals. True O/E=0 is retained. Axes use a linear scale near zero and a logarithmic scale above one. No correlation or combined p value is inferred from this selected union.</p>'),
        ("Counts", html_table([counts], [k for k in counts if k != "collapsed_identical_duplicates"])),
        ("Exact region comparison", f"<p>{len(results):,} union regions; the first {min(preview,len(results)):,} are previewed below. The complete table contains every union region. Region coordinates are zero-based half-open; matching region IDs alone never overrides different coordinates, hypotheses or sensitivities.</p>" + html_table(results, [*KEY_FIELDS, "left_oe", "right_oe", "left_status", "right_status", "left_observed", "right_observed", "left_fdr", "right_fdr"], limit=preview)),
        ("Downloads", '<p><a href="../results/concordance.tsv.gz">Complete comparison TSV</a> · <a href="../figures/concordance.png">PNG figure</a> · <a href="../qc/concordance.json">Definitions and support manifest</a></p>')])
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", required=True); parser.add_argument("--right", required=True); parser.add_argument("--out", required=True)
    args = parser.parse_args()
    try:
        result = compare_regions(args.left, args.right, args.out)
        print(json.dumps({k: result[k] for k in ("status", "comparability", "counts", "report")}))
    except (OSError, ValueError) as exc:
        raise SystemExit("highlight-compare: " + str(exc))
