#!/usr/bin/env python3
"""Compute primary HERV enrichment and validate current count-derived results."""
from __future__ import annotations

from collections import defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import tempfile
from analysis_integrity import transactional, InputPath, check_output_paths

MISSING = {"", ".", "NA", "NONE", "NULL"}
ROW_STATES = {
    "complete", "insufficient_loco_training", "insufficient_blocks",
    "zero_expected", "no_eligible_sites", "not_available",
}


def optional_number(value, field):
    if value is None or (isinstance(value, str) and value.strip().upper() in MISSING):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {field}: {value!r}")
    return result


def count_value(value, field):
    result = optional_number(value, field)
    if result is None or result < 0 or result != int(result):
        raise ValueError(f"Expected a nonnegative integer for {field}: {value!r}")
    return int(result)


def validate_primary_row(row):
    """Validate arithmetic, not a desired result or significance level."""
    r = dict(row)
    status = str(r.get("status", ""))
    if status not in ROW_STATES:
        raise ValueError(f"Unknown primary row status: {status!r}")
    if status == "not_available":
        if not r.get("reason"):
            raise ValueError("not_available requires an explicit reason")
        for k in ("observed_multi", "expected_multi", "oe", "ci_low", "ci_high", "n_eligible"):
            if optional_number(r.get(k), k) is not None:
                raise ValueError(f"not_available must not contain an estimate: {k}")
        return r
    n = count_value(r.get("n_eligible"), "n_eligible")
    o = count_value(r.get("observed_multi"), "observed_multi")
    missing = count_value(r.get("n_without_loco", 0), "n_without_loco")
    if o > n or missing > n:
        raise ValueError("Multi-ALT/unsupported site counts exceed eligible sites")
    e = optional_number(r.get("expected_multi"), "expected_multi")
    ratio = optional_number(r.get("oe"), "oe")
    lo = optional_number(r.get("ci_low"), "ci_low")
    hi = optional_number(r.get("ci_high"), "ci_high")
    if (lo is None) != (hi is None):
        raise ValueError("Both CI bounds must be present or both unavailable")
    if lo is not None and (lo < 0 or hi < lo):
        raise ValueError("Invalid O/E interval")
    if missing:
        if status != "insufficient_loco_training" or any(x is not None for x in (e, ratio, lo, hi)):
            raise ValueError("Unsupported full-scope estimate must remain unavailable")
    else:
        if e is None or e < 0 or e > n + max(1e-8, n * 1e-10):
            raise ValueError("Expected count is missing or outside the eligible universe")
        if e == 0:
            if ratio is not None or lo is not None:
                raise ValueError("Zero expected count has no estimable O/E")
        elif ratio is None or ratio < 0 or not math.isclose(ratio, o / e, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("O/E does not equal observed / expected")
    if status == "complete" and (n == 0 or missing or e is None or e <= 0):
        raise ValueError("complete requires an estimable point estimate")
    for k in ("bootstrap_p", "fdr"):
        x = optional_number(r.get(k), k)
        if x is not None and not 0 <= x <= 1:
            raise ValueError(f"{k} must be in [0,1]")
    r.update(observed_multi=o, n_eligible=n, n_without_loco=missing,
             expected_multi=e, oe=ratio, ci_low=lo, ci_high=hi,
             ci_status="available" if lo is not None else "unavailable")
    return r


def select_primary_rows(records, families=("ERV1", "ERVK"), thresholds=(5, 10),
                        sensitivity="exclude_combined"):
    wanted = {(f"family:{f}", int(t), sensitivity) for f in families for t in thresholds}
    selected = {}
    for row in records:
        key = (row["label"], int(row["threshold"]), row["sensitivity"])
        if key not in wanted:
            continue
        if key in selected:
            raise ValueError(f"Duplicate primary row: {key}")
        selected[key] = validate_primary_row(row)
    for label, threshold, sensitivity in wanted - selected.keys():
        selected[(label, threshold, sensitivity)] = {
            "label": label, "threshold": threshold, "sensitivity": sensitivity,
            "observed_multi": None, "expected_multi": None, "oe": None,
            "ci_low": None, "ci_high": None, "n_eligible": None,
            "status": "not_available", "reason": "no_matching_counter_summary_row",
        }
    return [selected[k] for k in sorted(selected)]


def _read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def merge_counter_files(paths, destination):
    """Stream existing counters; do not read or reconstruct raw VCF loci."""
    paths = [Path(p).resolve() for p in paths]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("Counter paths must be nonempty and unique")
    header = None
    with gzip.open(destination, "wt", newline="") as dst:
        writer = None
        for path in paths:
            opener = gzip.open if str(path).endswith(".gz") else open
            with opener(path, "rt", newline="") as src:
                reader = csv.DictReader(src, delimiter="\t")
                if not reader.fieldnames or (header is not None and header != reader.fieldnames):
                    raise ValueError("Counter tables require the same nonempty schema")
                if writer is None:
                    header = reader.fieldnames
                    writer = csv.DictWriter(dst, header, delimiter="\t")
                    writer.writeheader()
                writer.writerows(reader)


@transactional("primary_manifest.json")
def compute_primary_from_counters(modern, out, *, bootstrap=2000, seed=20261028,
                                  families=("ERV1", "ERVK"), minimum=1):
    """Compute the declared primary estimator from complete source counters."""
    import scan
    from herv_modern import architecture_paths, write_rows
    from analysis_integrity import identity
    modern, out = Path(modern).resolve(), Path(out).resolve()
    if bootstrap < 2000 or minimum != 1:
        raise ValueError("Primary estimator requires bootstrap >= 2000 and a LOCO training minimum of 1")
    gp = architecture_paths(modern, "global")
    ep = architecture_paths(modern, "exposure")
    protected = [InputPath(modern, "directory"), *[InputPath(p, "file") for p in (*gp, *ep)]]
    check_output_paths(out, protected)
    input_identities = [identity(p) for p in (*gp, *ep)]
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="counter-inputs-", dir=out) as temp:
        global_file = Path(temp) / "global.tsv.gz"
        exposure_file = Path(temp) / "exposure.tsv.gz"
        merge_counter_files(gp, global_file)
        merge_counter_files(ep, exposure_file)
        # All labels remain in this call so the estimator's FDR family is kept.
        result = scan.summarize_exposure(global_file, exposure_file, out,
                                         reps=bootstrap, seed=seed, min_training_sites=minimum)
    if input_identities != [identity(p) for p in (*gp, *ep)]:
        raise ValueError("Counter inputs changed while computing the primary result")
    selected = select_primary_rows(_read_tsv(out / "enrichment_summary.tsv"), families=families)
    write_rows(out / "primary_enrichment.tsv", selected,
               fields=["label", "threshold", "sensitivity", "observed_multi", "expected_multi", "oe",
                       "ci_low", "ci_high", "ci_status", "bootstrap_p", "fdr", "n_eligible",
                       "n_without_loco", "block_count", "status", "reason"])
    plot_primary(selected, out)
    run_path = modern / "run.json"
    run = scan.load_run_metadata(modern) if run_path.is_file() else {}
    manifest = dict(result, schema="herv-primary-1", status="complete",
                    source_counter_inputs=input_identities,
                    source_run=identity(run_path) if run_path.is_file() else None,
                    source_options=run.get("options", run.get("identity", {}).get("options", {})),
                    minimum_training_sites=minimum,
                    primary_rows=len(selected), estimable_primary_rows=sum(r.get("oe") is not None for r in selected),
                    estimand="multi_ALT_occupancy_among_fixed_eligible_SNV_sites",
                    target_scope="exposure_counter_sites_inside_recorded_EUR_LD_partition",
                    background_scope="all_chromosomes_recorded_in_global_counters",
                    missing_rule="any_unsupported_target_stratum_makes_full_scope_OE_unavailable")
    (out / "primary_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    return manifest


def plot_primary(data, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = list(data)
    fig, ax = plt.subplots(figsize=(9, max(2.5, .65 * len(data) + 1.2)))
    extent = [1.0]
    for i, r in enumerate(data):
        value = optional_number(r.get("oe"), "oe")
        lo, hi = (optional_number(r.get(k), k) for k in ("ci_low", "ci_high"))
        if value is None:
            ax.text(.02, i, str(r.get("status", "not_available")),
                    transform=ax.get_yaxis_transform(), va="center", color=".5")
            continue
        extent.append(value)
        ax.plot(value, i, "o", color="#176B87")
        if lo is not None and hi is not None:
            extent.extend((lo, hi)); ax.plot((lo, hi), (i, i), color="#176B87")
        label = f"{value:.3g}" + (f" ({lo:.3g}–{hi:.3g})" if lo is not None else " (CI unavailable)")
        ax.text(1.02, i, label, transform=ax.get_yaxis_transform(), va="center", fontsize=9)
    low, high = min(extent), max(extent)
    padding = max(.05, .1 * (high - low))
    ax.set_xlim(max(0, low - padding), high + padding)
    ax.axvline(1, color=".5", ls="--")
    ax.set(yticks=list(range(len(data))),
           yticklabels=[f"{r['label'].removeprefix('family:')}  AC≥{r['threshold']}" for r in data],
           xlabel="Observed / expected multi-ALT SNVs", title="ERV family enrichment")
    ax.set_ylim(len(data) - .4, -.6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out) / f"primary_enrichment.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def scoped_counter_totals(exposure_paths, background, families, chromosomes):
    """Return [raw O, partial E, raw N, unsupported O, unsupported N]."""
    selected_chromosomes = {str(c).removeprefix("chr") for c in chromosomes}
    labels = {"family:" + f for f in families}
    totals = defaultdict(lambda: [0, 0., 0, 0, 0])
    for path in exposure_paths:
        for r in _read_tsv(path):
            chrom = str(r["chrom"]).removeprefix("chr")
            if (chrom not in selected_chromosomes or r["population_mode"] != "all"
                or r["label"] not in labels or str(r["threshold"]) not in {"5", "10"}
                or r["sensitivity"] not in {"all", "exclude_combined"}
                or r["block_id"] == "outside_partition"):
                continue
            key = (r["label"][7:], "strict" if r["sensitivity"] == "exclude_combined" else "all", int(r["threshold"]))
            n, o = count_value(r["n_sites"], "n_sites"), count_value(r["n_multi"], "n_multi")
            if o > n: raise ValueError("Counter n_multi exceeds n_sites")
            rate = background.rate(chrom, background.key(r))
            v = totals[key]; v[0] += o; v[2] += n
            if rate is None: v[3] += o; v[4] += n
            else:
                if not 0 <= rate <= 1: raise ValueError("Background rate outside [0,1]")
                v[1] += n * rate
    return dict(totals)


def reconcile_spatial_to_counters(counter_totals, partitions, *, scope="eur_ld_partition"):
    """Reconcile measurements over the same target domain, with no fixed table."""
    targets = {}
    for row in partitions:
        if row.get("partition") != "total" or row.get("scope") != scope: continue
        key = row["family"], row["sensitivity"], int(row["threshold"])
        if key in targets: raise ValueError(f"Duplicate spatial total: {key}")
        targets[key] = row
    audit, failures = [], []
    for key in sorted(set(counter_totals) | set(targets)):
        o, partial_e, n, missing_o, missing_n = counter_totals.get(key, [0, 0., 0, 0, 0])
        target = targets.get(key)
        expected = None if missing_n else partial_e
        ratio = o / expected if expected else None
        expected_values = {"observed_multi": o, "eligible_sites": n, "n_without_loco": missing_n,
                           "multi_without_loco": missing_o, "partial_expected_multi": partial_e,
                           "expected_multi": expected, "oe": ratio}
        mismatched = []
        if target is None:
            if n: mismatched.append("missing_spatial_total")
        else:
            for field, value in expected_values.items():
                actual = optional_number(target.get(field), field)
                matches = actual is value if actual is None or value is None else math.isclose(actual, value, rel_tol=1e-10, abs_tol=1e-8)
                if not matches: mismatched.append(field)
        row = dict(family=key[0], sensitivity=key[1], threshold=key[2], scope=scope,
                   status="mismatch" if mismatched else "no_eligible_sites" if n == 0 else "matched",
                   mismatch_fields=",".join(mismatched), counter_observed=o,
                   counter_expected=expected, counter_eligible=n, counter_without_loco=missing_n)
        audit.append(row)
        if mismatched: failures.append(row)
    return audit, failures


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--modern-dir', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--families', default='ERV1,ERVK')
    p.add_argument('--bootstrap', type=int, default=2000)
    p.add_argument('--seed', type=int, default=20261028)
    a = p.parse_args(argv)
    result = compute_primary_from_counters(a.modern_dir, a.out, bootstrap=a.bootstrap,
        seed=a.seed, families=tuple(a.families.split(',')))
    print(json.dumps(result, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
