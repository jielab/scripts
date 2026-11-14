#!/usr/bin/env python3
"""Consolidated herv workflow; internal stages share one explicit module."""
from __future__ import annotations


# ========================================================================
# herv primary
# ========================================================================
_description_herv_primary = 'Compute primary HERV enrichment and validate current count-derived results.'
#!/usr/bin/env python3
"""Compute primary HERV enrichment and validate current count-derived results."""

from collections import defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import tempfile
from analysis_integrity import transactional, InputPath, check_output_paths
from batch_workflows import boolean_option

PRIMARY_MISSING = {"", ".", "NA", "NONE", "NULL"}
PRIMARY_ROW_STATES = {
    "complete", "insufficient_loco_training", "insufficient_blocks",
    "zero_expected", "no_eligible_sites", "not_available",
}


def primary_optional_number(value, field):
    if value is None or (isinstance(value, str) and value.strip().upper() in PRIMARY_MISSING):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {field}: {value!r}")
    return result


def primary_count_value(value, field):
    result = primary_optional_number(value, field)
    if result is None or result < 0 or result != int(result):
        raise ValueError(f"Expected a nonnegative integer for {field}: {value!r}")
    return int(result)


def primary_validate_primary_row(row):
    """Validate arithmetic, not a desired result or significance level."""
    r = dict(row)
    status = str(r.get("status", ""))
    if status not in PRIMARY_ROW_STATES:
        raise ValueError(f"Unknown primary row status: {status!r}")
    if status == "not_available":
        if not r.get("reason"):
            raise ValueError("not_available requires an explicit reason")
        for k in ("observed_multi", "expected_multi", "oe", "ci_low", "ci_high", "n_eligible"):
            if primary_optional_number(r.get(k), k) is not None:
                raise ValueError(f"not_available must not contain an estimate: {k}")
        return r
    n = primary_count_value(r.get("n_eligible"), "n_eligible")
    o = primary_count_value(r.get("observed_multi"), "observed_multi")
    missing = primary_count_value(r.get("n_without_loco", 0), "n_without_loco")
    if o > n or missing > n:
        raise ValueError("Multi-ALT/unsupported site counts exceed eligible sites")
    e = primary_optional_number(r.get("expected_multi"), "expected_multi")
    ratio = primary_optional_number(r.get("oe"), "oe")
    lo = primary_optional_number(r.get("ci_low"), "ci_low")
    hi = primary_optional_number(r.get("ci_high"), "ci_high")
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
        x = primary_optional_number(r.get(k), k)
        if x is not None and not 0 <= x <= 1:
            raise ValueError(f"{k} must be in [0,1]")
    r.update(observed_multi=o, n_eligible=n, n_without_loco=missing,
             expected_multi=e, oe=ratio, ci_low=lo, ci_high=hi,
             ci_status="available" if lo is not None else "unavailable")
    return r


def primary_select_primary_rows(records, families=("ERV1", "ERVK"), thresholds=(5, 10),
                        sensitivity="exclude_combined"):
    wanted = {(f"family:{f}", int(t), sensitivity) for f in families for t in thresholds}
    selected = {}
    for row in records:
        key = (row["label"], int(row["threshold"]), row["sensitivity"])
        if key not in wanted:
            continue
        if key in selected:
            raise ValueError(f"Duplicate primary row: {key}")
        selected[key] = primary_validate_primary_row(row)
    for label, threshold, sensitivity in wanted - selected.keys():
        selected[(label, threshold, sensitivity)] = {
            "label": label, "threshold": threshold, "sensitivity": sensitivity,
            "observed_multi": None, "expected_multi": None, "oe": None,
            "ci_low": None, "ci_high": None, "n_eligible": None,
            "status": "not_available", "reason": "no_matching_counter_summary_row",
        }
    return [selected[k] for k in sorted(selected)]


def primary_read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def primary_merge_counter_files(paths, destination, *, exclude_ibdmix=False):
    """Stream existing counters; do not read or reconstruct raw VCF loci."""
    paths = [Path(p).resolve() for p in paths]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("Counter paths must be nonempty and unique")
    header = None
    excluded = 0
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
                for row in reader:
                    if exclude_ibdmix and modern_archaic_label_reference(row.get("label", "")) is not None:
                        excluded += 1
                        continue
                    writer.writerow(row)
    return excluded


@transactional("primary_manifest.json")
def primary_compute_primary_from_counters(modern, out, *, bootstrap=2000, seed=20261028,
                                  families=("ERV1", "ERVK"), minimum=1, include_ibdmix=False):
    """Compute the declared primary estimator from complete source counters."""
    import scan
    from analysis_integrity import identity
    modern, out = Path(modern).resolve(), Path(out).resolve()
    if bootstrap < 2000 or minimum != 1:
        raise ValueError("Primary estimator requires bootstrap >= 2000 and a LOCO training minimum of 1")
    gp = modern_architecture_paths(modern, "global")
    ep = modern_architecture_paths(modern, "exposure")
    protected = [InputPath(modern, "directory"), *[InputPath(p, "file") for p in (*gp, *ep)]]
    check_output_paths(out, protected)
    input_identities = [identity(p) for p in (*gp, *ep)]
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="counter-inputs-", dir=out) as temp:
        global_file = Path(temp) / "global.tsv.gz"
        exposure_file = Path(temp) / "exposure.tsv.gz"
        primary_merge_counter_files(gp, global_file)
        excluded = primary_merge_counter_files(ep, exposure_file, exclude_ibdmix=not include_ibdmix)
        # Remove optional evidence before inference, including the BH family.
        result = scan.summarize_exposure(global_file, exposure_file, out,
                                         reps=bootstrap, seed=seed, min_training_sites=minimum)
    if input_identities != [identity(p) for p in (*gp, *ep)]:
        raise ValueError("Counter inputs changed while computing the primary result")
    selected = primary_select_primary_rows(primary_read_tsv(out / "enrichment_summary.tsv"), families=families)
    modern_write_rows(out / "primary_enrichment.tsv", selected,
               fields=["label", "threshold", "sensitivity", "observed_multi", "expected_multi", "oe",
                       "ci_low", "ci_high", "ci_status", "bootstrap_p", "fdr", "n_eligible",
                       "n_without_loco", "block_count", "status", "reason"])
    primary_plot_primary(selected, out)
    run_path = modern / "run.json"
    run = scan.load_run_metadata(modern) if run_path.is_file() else {}
    manifest = dict(result, schema="herv-primary-1", status="complete",
                    ibdmix_enabled=include_ibdmix, excluded_ibdmix_counter_rows=excluded,
                    testing_family="all retained annotation labels after optional IBDmix exclusion",
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


def primary_plot_primary(data, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = list(data)
    fig, ax = plt.subplots(figsize=(9, max(2.5, .65 * len(data) + 1.2)))
    extent = [1.0]
    for i, r in enumerate(data):
        value = primary_optional_number(r.get("oe"), "oe")
        lo, hi = (primary_optional_number(r.get(k), k) for k in ("ci_low", "ci_high"))
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


def primary_scoped_counter_totals(exposure_paths, background, families, chromosomes):
    """Return [raw O, partial E, raw N, unsupported O, unsupported N]."""
    selected_chromosomes = {str(c).removeprefix("chr") for c in chromosomes}
    labels = {"family:" + f for f in families}
    totals = defaultdict(lambda: [0, 0., 0, 0, 0])
    for path in exposure_paths:
        for r in primary_read_tsv(path):
            chrom = str(r["chrom"]).removeprefix("chr")
            if (chrom not in selected_chromosomes or r["population_mode"] != "all"
                or r["label"] not in labels or str(r["threshold"]) not in {"5", "10"}
                or r["sensitivity"] not in {"all", "exclude_combined"}
                or r["block_id"] == "outside_partition"):
                continue
            key = (r["label"][7:], "strict" if r["sensitivity"] == "exclude_combined" else "all", int(r["threshold"]))
            n, o = primary_count_value(r["n_sites"], "n_sites"), primary_count_value(r["n_multi"], "n_multi")
            if o > n: raise ValueError("Counter n_multi exceeds n_sites")
            rate = background.rate(chrom, background.key(r))
            v = totals[key]; v[0] += o; v[2] += n
            if rate is None: v[3] += o; v[4] += n
            else:
                if not 0 <= rate <= 1: raise ValueError("Background rate outside [0,1]")
                v[1] += n * rate
    return dict(totals)


def primary_reconcile_spatial_to_counters(counter_totals, partitions, *, scope="eur_ld_partition"):
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
                actual = primary_optional_number(target.get(field), field)
                matches = actual is value if actual is None or value is None else math.isclose(actual, value, rel_tol=1e-10, abs_tol=1e-8)
                if not matches: mismatched.append(field)
        row = dict(family=key[0], sensitivity=key[1], threshold=key[2], scope=scope,
                   status="mismatch" if mismatched else "no_eligible_sites" if n == 0 else "matched",
                   mismatch_fields=",".join(mismatched), counter_observed=o,
                   counter_expected=expected, counter_eligible=n, counter_without_loco=missing_n)
        audit.append(row)
        if mismatched: failures.append(row)
    return audit, failures


def primary_main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=_description_herv_primary)
    p.add_argument('--modern-dir', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--families', default='ERV1,ERVK')
    p.add_argument('--bootstrap', type=int, default=2000)
    p.add_argument('--seed', type=int, default=20261028)
    p.add_argument('--ibdmix', type=boolean_option, default=False, metavar='TRUE|FALSE')
    a = p.parse_args(argv)
    result = primary_compute_primary_from_counters(a.modern_dir, a.out, bootstrap=a.bootstrap,
        seed=a.seed, families=tuple(a.families.split(',')), include_ibdmix=a.ibdmix)
    print(json.dumps(result, ensure_ascii=False))
    return 0


def _entry_herv_primary():
    raise SystemExit(primary_main())


# ========================================================================
# herv modern
# ========================================================================
_description_herv_modern = 'HERV common-support summaries, body decomposition and held-out edge scans.\n\nMulti-ALT occupancy is conditional on the original eligible SNV universe.\nReference-allele observation and longitudinal persistence are not inferred.'
#!/usr/bin/env python3
"""HERV common-support summaries, body decomposition and held-out edge scans.

Multi-ALT occupancy is conditional on the original eligible SNV universe.
Reference-allele observation and longitudinal persistence are not inferred.
"""
import argparse
from bisect import bisect_left
from collections import Counter
import itertools
import os
import subprocess
import time
import numpy as np
import scan
from allele_common import (FastaContext, IntervalMask, canonical_chrom,
                           dump_json, open_text, _parse_header, _record, af_bin, an_bin)

MODERN_SCHEMA = "herv-focus-modern-1.0"
MODERN_AUTOSOMES = {str(i) for i in range(1, 23)}


def modern_rows(path):
    with open_text(path) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def modern_write_rows(path, data, fields=None):
    data = list(data)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    columns = fields or (list(data[0]) if data else ["status"])
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in data:
            writer.writerow({k: "NA" if v is None else v for k, v in row.items()})


def modern_file_state(path):
    p = Path(path).resolve()
    s = p.stat()
    return {"path": str(p), "bytes": s.st_size, "mtime_ns": s.st_mtime_ns}


def modern_architecture_paths(modern, kind):
    combined = Path(modern) / "summary" / f"architecture.{kind}_strata.tsv.gz"
    if combined.is_file():
        return [combined]
    run = scan.load_run(modern)
    prefix = run["options"].get("architecture_prefix", "architecture")
    files = [Path(s) / f"{prefix}.{kind}_strata.tsv.gz" for s in run["shards"]]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise ValueError(f"Missing existing {kind} counters: {missing[:3]}; original caches are read-only")
    return files


class modern_Background:
    """Fixed base denominator; validation uses only discovery chromosomes."""
    def __init__(self, paths, min_training_sites=100, discovery=None, include_all_chromosomes=False):
        self.counts = defaultdict(lambda: np.zeros(2))
        self.total = defaultdict(lambda: np.zeros(2))
        self.train = defaultdict(lambda: np.zeros(2))
        self.discovery = {canonical_chrom(c) for c in (discovery or [str(i) for i in range(1, 23, 2)])}
        if not self.discovery <= MODERN_AUTOSOMES:
            raise ValueError("Discovery chromosomes must be autosomes")
        self.minimum = min_training_sites
        for path in paths:
            for r in modern_rows(path):
                if r["population_mode"] != "all" or r["threshold"] not in {"5", "10"}:
                    continue
                chrom = canonical_chrom(r["chrom"])
                if not include_all_chromosomes and chrom not in MODERN_AUTOSOMES:
                    continue
                key = self.key(r)
                value = np.array([int(r["n_sites"]), int(r["n_multi"])])
                self.counts[(chrom, key)] += value
                self.total[key] += value
                if chrom in self.discovery:
                    self.train[key] += value

    @staticmethod
    def key(r):
        sensitivity = r["sensitivity"].replace("_no_cpg", "")
        if sensitivity == "strict":
            sensitivity = "exclude_combined"
        return (r["context"], r["af_bin"], r["an_bin"], sensitivity, str(r["threshold"]))

    def rate(self, chrom, key, split="all"):
        if split == "all":
            n, m = self.total[key] - self.counts[(chrom, key)]
        elif split == "discovery":
            n, m = self.train[key] - self.counts[(chrom, key)]
        else:
            n, m = self.train[key]
        return None if n < self.minimum else (m + .5) / (n + 1.)


class modern_PrimaryBackground(modern_Background):
    """Exact point-estimate definition in scan.summarize_exposure.

    Include every recorded global chromosome, train with at least one eligible
    site, and never convert an unsupported full-scope total into common support.
    The exposure audit is collected during the existing legacy-counter read.
    """
    def __init__(self, paths, families):
        super().__init__(paths, min_training_sites=1, include_all_chromosomes=True)
        self.families = set(families)
        self.labels = {"family:"+f for f in self.families}
        self.exposure = defaultdict(lambda: np.zeros(5))
        self.exposure_chromosomes = defaultdict(set)

    def observe_exposure(self, r):
        if (r["population_mode"] != "all" or r["threshold"] not in {"5", "10"}
                or r["sensitivity"] not in {"all", "exclude_combined"}
                or r["label"] not in self.labels
                or r["block_id"] == "outside_partition"):
            return
        c = canonical_chrom(r["chrom"])
        sensitivity = "strict" if r["sensitivity"] == "exclude_combined" else "all"
        key = r["label"][7:], sensitivity, int(r["threshold"])
        n, m = int(r["n_sites"]), int(r["n_multi"])
        value = self.exposure[key]
        # Raw observed, partial expected, raw denominator, unsupported M/N.
        value[0] += m
        value[2] += n
        if n:
            self.exposure_chromosomes[key].add(c)
        rate = self.rate(c, self.key(r))
        if rate is None:
            value[3:] += [m, n]
        else:
            value[1] += n*rate


def modern_archaic_label_reference(label):
    """Read only the exact dataset::archaic:REFERENCE:exposure namespace."""
    marker = "::archaic:"
    if marker not in label:
        return None
    if label.count(marker) != 1:
        raise ValueError(f"Ambiguous archaic exposure label: {label}")
    dataset, encoded = label.split(marker)
    reference, separator, exposure = encoded.partition(":")
    if not dataset or not reference or not separator or not exposure:
        raise ValueError(f"Malformed archaic exposure label: {label}")
    return reference


@transactional("summary_manifest.json")
def modern_summarize_common_support(modern, out, minimum=100, coverage=.99, bootstrap=2000, seed=20261027,
                             archaic_references=None, *, include_ibdmix=False):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    selected_references = scan.track_references(archaic_references) if include_ibdmix else ()
    reference_set = set(selected_references)
    reference_audit = Counter()
    background = modern_Background(modern_architecture_paths(modern, "global"), minimum)
    # Supported observed/expected/denominator, then unsupported obs/denominator.
    blocks = defaultdict(lambda: defaultdict(lambda: np.zeros(5)))
    outside = Counter()
    for path in modern_architecture_paths(modern, "exposure"):
        for r in modern_rows(path):
            if r["population_mode"] != "all" or r["threshold"] not in {"5", "10"}:
                continue
            c = canonical_chrom(r["chrom"])
            if c not in MODERN_AUTOSOMES:
                continue
            reference = modern_archaic_label_reference(r["label"])
            if reference is not None:
                selected = reference in reference_set
                reference_audit[(r["label"], reference, selected)] += 1
                if not selected:
                    continue
            key = r["label"], int(r["threshold"]), r["sensitivity"]
            n, m = int(r["n_sites"]), int(r["n_multi"])
            if r["block_id"] == "outside_partition":
                outside[key] += n
                continue
            v = blocks[key][c + ":" + r["block_id"]]
            rate = background.rate(c, background.key(r))
            if rate is None:
                v[3:] += [m, n]
            else:
                v[:3] += [m, n * rate, n]
    results = []
    for key, vals in sorted(blocks.items()):
        o, e, n, removed_o, removed_n = np.sum(list(vals.values()), axis=0)
        supported = n / (n + removed_n) if n + removed_n else 0.
        usable = {b: v[:2] for b, v in vals.items() if v[2] > 0}
        low, high, p = scan.bootstrap_statistics(usable, bootstrap, seed) if e > 0 else (None, None, None)
        status = "complete" if supported >= coverage and len(usable) >= 20 and e > 0 else "low_support_or_few_blocks"
        results.append(dict(label=key[0], threshold=key[1], sensitivity=key[2],
                            observed_multi=o, expected_multi=e, oe=o/e if e else None,
                            ci_low=low, ci_high=high, bootstrap_p=p, fdr=None,
                            supported_sites=int(n), excluded_sites=int(removed_n), excluded_multi=int(removed_o),
                            supported_fraction=supported, block_count=len(usable),
                            outside_partition_sites=outside[key], status=status))
    tested = [r for r in results if r["status"] == "complete" and r["bootstrap_p"] is not None]
    for r, q in zip(tested, scan.bh_fdr([r["bootstrap_p"] for r in tested])):
        r["fdr"] = q
    modern_write_rows(out / "common_support_enrichment.tsv", results)
    selection_rows = [dict(label=label, reference=reference, decision="retained" if selected else "excluded",
                           reason="ibdmix_disabled" if not include_ibdmix else "selected_reference" if selected else "reference_not_selected", counter_rows=count)
                      for (label, reference, selected), count in sorted(reference_audit.items())]
    modern_write_rows(out / "archaic_reference_selection.tsv", selection_rows,
               fields=["label", "reference", "decision", "reason", "counter_rows"])
    manifest = dict(schema=MODERN_SCHEMA, status="complete", source=str(Path(modern).resolve()),
                    ibdmix_enabled=include_ibdmix,
                    summary_definition_version="common-support-reference-selection-2",
                    selected_archaic_references=list(selected_references),
                    observed_archaic_references=sorted({reference for _, reference, _ in reference_audit}),
                    selected_archaic_references_not_present=sorted(reference_set - {reference for _, reference, _ in reference_audit}),
                    excluded_archaic_labels=sorted({label for label, _, selected in reference_audit if not selected}),
                    excluded_archaic_counter_rows=sum(n for (_, _, selected), n in reference_audit.items() if not selected),
                    reference_selection="Exact dataset::archaic:REFERENCE:exposure matching before aggregation and BH; other annotation labels unchanged",
                    reference_audit_scope="Autosomal all-population AC5/10 counter rows; row counts are not unique genomic bases",
                    minimum_training_sites=minimum, minimum_support=coverage, bootstrap=bootstrap,
                    denominator="original eligible SNVs; unsupported strata removed from BOTH numerator and denominator",
                    adjustment="LOCO canonical trinucleotide, primary ALT AF and AN; base denominator fixed across AC tiers",
                    inference="EUR LD-block bootstrap; background offsets treated as fixed; BH over all complete rows",
                    interpretation="cross-sectional multi-ALT occupancy, not selection or temporal persistence")
    dump_json(out / "summary_manifest.json", manifest)
    modern_plot_families(results, out)
    return manifest


def modern_plot_families(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    families = ["ERV1", "ERVK", "ERVL", "ERVL-MaLR"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for t, color, offset in [(5, "#176B87", -.12), (10, "#D16B38", .12)]:
        labelled = False
        for i, family in enumerate(families):
            rr = [r for r in results if r["label"] == "family:"+family and int(r["threshold"]) == t and r["sensitivity"] == "exclude_combined"]
            if rr and rr[0]["oe"] not in (None, "NA"):
                r = rr[0]
                ax.plot(float(r["oe"]), i+offset, "o", color=color, label=f"AC ≥ {t}" if not labelled else None)
                labelled = True
                if r["ci_low"] not in (None, "NA"):
                    ax.plot([float(r["ci_low"]), float(r["ci_high"])], [i+offset]*2, color=color)
    ax.axvline(1, color=".6", ls="--")
    ax.set(yticks=range(4), yticklabels=families, xlabel="Multi-ALT observed / expected", title="HERV families · common support · strict masks")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"family_enrichment.{ext}", dpi=180)
    plt.close(fig)


class modern_Edges:
    """Nearest edge per family/end; base centres avoid strand off-by-one bias."""
    def __init__(self, elements, families):
        self.points = defaultdict(list)
        self.coords = {}
        self.ambiguous_edges = 0
        for e in elements:
            if e.rep_class == "LTR" and e.family in families:
                for end, point in (("5p", e.five_prime), ("3p", e.three_prime)):
                    self.points[e.family, end].append((point, e.strand, e.name, e.source_id))
        for key, v in self.points.items():
            unique = []
            for point, group in itertools.groupby(sorted(v), key=lambda z:z[0]):
                coincident = list(group)
                strands = {z[1] for z in coincident}
                if len(strands) > 1:
                    # The left/right sides must not silently choose different
                    # strands at the same coordinate and create a false peak.
                    self.ambiguous_edges += 1
                    unique.append((point, "?", "ambiguous_orientation", ""))
                else:
                    unique.append(coincident[0])
            self.points[key] = unique
            self.coords[key] = [z[0] for z in unique]

    def locate(self, pos1, flank):
        p = pos1 - .5
        found = []
        for key, vals in self.points.items():
            i = bisect_left(self.coords[key], p)
            candidates = vals[max(0, i-1):min(len(vals), i+1)]
            if not candidates:
                continue
            point, strand, name, source = min(candidates, key=lambda x:(abs(p-x[0]), x))
            if strand == "?":
                continue
            distance = (p-point) * (1 if strand == "+" else -1)
            if -flank <= distance < flank:
                found.append((*key, distance))
        return found

    def intervals(self, flank, length):
        return scan.merge_intervals([(max(0, p-flank), min(length, p+flank))
                                     for vals in self.points.values() for p, *_ in vals
                                     if max(0, p-flank) < min(length, p+flank)])


def modern_query_regions(vcf, fields, bed, bcftools="bcftools"):
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER" + "".join("\t%INFO/"+f for f in fields) + "\n"
    cmd = [bcftools, "query", "--allow-undef-tags", "--regions-overlap", "0", "-R", str(bed), "-f", fmt, str(vcf)]
    with tempfile.TemporaryFile(mode="w+t") as errors:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errors, text=True)
        try:
            previous = None
            for line in proc.stdout:
                record = _record(line.rstrip("\r\n").split("\t"), fields)
                key = (canonical_chrom(record["chrom"]), record["pos"])
                if previous and key[0] == previous[0] and key[1] < previous[1]:
                    raise ValueError("Indexed query was not coordinate sorted")
                previous = key
                yield record
            code = proc.wait()
            if code:
                errors.seek(0)
                raise RuntimeError(f"bcftools query failed ({code}): {errors.read()[-2000:]}")
        finally:
            proc.stdout.close()
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()


def modern_strict_site_excluded(chrom, pos, flags, repeats, masks):
    """Exactly match ArchitectureAccumulator's BED-first annotation policy."""
    flags = set(flags)
    excluded = {name: masks[name].contains(chrom, pos) if name in masks else name in flags
                for name in ("lcr", "segdup")}
    excluded["low_mappability"] = (masks["low_mappability"].contains(chrom, pos)
                                  if "low_mappability" in masks else False)
    excluded["simple_repeat"] = any(e.rep_class in {"Simple_repeat", "Low_complexity"} for e in repeats)
    return any(excluded.values())


def modern_read_indexed_header(vcf, bcftools):
    result = subprocess.run([bcftools, "view", "--no-version", "-h", str(vcf)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode:
        raise ValueError(f"Cannot read VCF header: {result.stderr[-2000:]}")
    return _parse_header(result.stdout)


def modern_scan_boundary_chrom(vcf, chrom, annotation, out, options, args):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"chr{chrom}.bins.tsv.gz"
    marker = out / f"chr{chrom}.json"
    fasta = args.fasta or options["fasta"]
    mask_paths = {n: options.get(k) for n, k in [("segdup", "architecture_segdup_bed"),
                  ("lcr", "architecture_lcr_bed"), ("low_mappability", "architecture_low_mappability_bed")]}
    for name in mask_paths:
        value = getattr(args, name+"_bed", None)
        if value:
            mask_paths[name] = value
    inputs = [vcf, annotation, fasta, fasta+".fai", *[p for p in mask_paths.values() if p]]
    if options.get("callable_bed"):
        inputs.append(options["callable_bed"])
    for suffix in (".tbi", ".csi"):
        if Path(str(vcf)+suffix).is_file():
            inputs.append(str(vcf)+suffix)
            break
    else:
        raise ValueError(f"Indexed VCF required: {vcf}")
    identity = dict(schema=MODERN_SCHEMA, inputs=[modern_file_state(p) for p in inputs], source_options=options,
                    code=[modern_file_state(p) for p in (Path(__file__), Path(scan.__file__), Path(__file__).with_name("0.common.py"))],
                    families=args.families, flank=args.flank, step=args.step, block_bp=args.block_bp)
    if marker.exists():
        old = json.loads(marker.read_text())
        if old["identity"] != identity:
            raise ValueError(f"Boundary inputs/settings/code changed; use a fresh output: {out}")
        if path.exists() and old.get("output") == modern_file_state(path):
            print(f"[resume] boundary chr{chrom}", flush=True)
            return path
        raise ValueError(f"Completed boundary cache changed: {path}")
    header = modern_read_indexed_header(vcf, args.bcftools)
    scan.validate_count_header(header)
    fields = ["AC", "AN"] + [f for f in ("AF", "AS_FilterStatus", "lcr", "segdup") if f in header["info"]]
    elements = scan.parse_joined_rmsk(annotation, chromosome=chrom)
    edges = modern_Edges(elements, args.families.split(","))
    repeat_sweep = scan.ERVOverlapSweep(elements)
    masks = {n: IntervalMask(p, [chrom]) for n, p in mask_paths.items() if p}
    strict_complete = all(n in masks or n in header["info"] for n in ("lcr", "segdup")) and "low_mappability" in masks
    callable_mask = IntervalMask(options["callable_bed"], [chrom]) if options.get("callable_bed") else None
    qc = Counter()
    partial = Path(str(path)+".partial")
    columns = ["chrom", "block", "family", "boundary", "bin_start", "sensitivity", "threshold",
               "context", "af_bin", "an_bin", "n_sites", "n_multi"]
    bins = defaultdict(lambda: [0, 0])
    last_block = None
    with FastaContext(fasta) as fa, gzip.open(partial, "wt", compresslevel=2, newline="") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(columns)
        def flush():
            writer.writerows([*k, *v] for k, v in sorted(bins.items()))
            bins.clear()
        raw_chroms = [c for c in header["contigs"] if canonical_chrom(c) == chrom]
        if len(raw_chroms) != 1:
            raise ValueError(f"No unique VCF contig for chr{chrom}")
        intervals = edges.intervals(args.flank, fa.lengths[chrom])
        bed = out / f"chr{chrom}.regions.bed"
        with bed.open("w") as h:
            for lo, hi in intervals:
                h.write(f"{raw_chroms[0]}\t{lo}\t{hi}\n")
        records = modern_query_regions(vcf, fields, bed, args.bcftools) if intervals else iter(())
        clock = time.monotonic()
        for (c, pos), group in itertools.groupby(records, key=lambda r:(canonical_chrom(r["chrom"]), r["pos"])):
            if c != chrom:
                raise ValueError("Unexpected chromosome from targeted query")
            if callable_mask and not callable_mask.contains(c, pos):
                continue
            near = edges.locate(pos, args.flank)
            if not near:
                continue
            ref, alts, non_snv, flags = scan.group_alleles(group, header, [], options, qc, None, [])
            observed = sorted((a for a in alts.values() if a["status"] == "observed_pass"),
                              key=lambda a:(-a["af"], -a["ac"], a["alt"]))
            if not observed:
                continue
            ctx = fa.context(c, pos)
            if ctx[1] != ref:
                raise ValueError(f"REF/FASTA mismatch: {c}:{pos}")
            context, _ = scan.canonical_context(ctx, ref)
            if any(b not in "ACGT" for b in context):
                # The original scanner retains these eligible loci. Removing
                # them here would silently change the architecture denominator.
                qc["ambiguous_context_retained"] += 1
            primary = observed[0]
            strata = context, af_bin(primary["af"]), an_bin(primary["an"])
            repeats = repeat_sweep.locate(c, pos)
            excluded = modern_strict_site_excluded(c, pos, flags, repeats, masks)
            sensitivities = ["all"]
            if strict_complete and not excluded:
                sensitivities.append("strict")
            cpg = (ref == "C" and ctx[2] == "G") or (ref == "G" and ctx[0] == "C")
            if not cpg:
                sensitivities += [s+"_no_cpg" for s in list(sensitivities)]
            block = f"{c}:{(pos-1)//args.block_bp}"
            if last_block is not None and last_block != block:
                flush()
            last_block = block
            for family, end, distance in near:
                bin_start = math.floor(distance/args.step)*args.step
                for sensitivity in sensitivities:
                    for threshold in (5, 10):
                        key = (c, block, family, end, bin_start, sensitivity, threshold, *strata)
                        v = bins[key]
                        v[0] += 1
                        v[1] += sum(a["ac"] >= threshold for a in observed) >= 2
            qc["eligible_sites"] += 1
            if time.monotonic()-clock > 60:
                print(f"[boundaries chr{chrom}] {qc['eligible_sites']:,} eligible sites", flush=True)
                clock = time.monotonic()
        flush()
    os.replace(partial, path)
    dump_json(marker, dict(identity=identity, output=modern_file_state(path), qc=dict(qc),
                          strict_masks_complete=strict_complete, ambiguous_oriented_edges=edges.ambiguous_edges,
                          queried_bp=sum(b-a for a,b in intervals)))
    return path


def modern_cluster_estimate(matrix, local=False):
    """Log-ratio cluster sandwich; repeated windows share the same block scores."""
    sums = matrix.sum(axis=0)
    if (sums[:2] <= 0).any() or local and (sums[2:4] <= 0).any():
        return None
    theta = sums[0]/sums[1]
    influence = matrix[:, 0]/sums[0] - matrix[:, 1]/sums[1]
    if local:
        theta /= sums[2]/sums[3]
        influence -= matrix[:, 2]/sums[2] - matrix[:, 3]/sums[3]
    used = (matrix[:, :4 if local else 2].sum(axis=1) > 0).sum()
    correction = math.sqrt(used/max(used-1, 1))
    influence *= correction
    se = float(np.sqrt(np.sum(influence**2)))
    if not se or used < 2:
        return None
    from scipy.stats import norm
    z = math.log(theta)/se
    return dict(estimate=float(theta), ci_low=float(np.exp(math.log(theta)-1.96*se)),
                ci_high=float(np.exp(math.log(theta)+1.96*se)), pointwise_p=float(2*norm.sf(abs(z))),
                z=z, influence=influence/se, n_blocks=int(used))


def modern_profile_counts(paths, background, args):
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: np.zeros(4))))
    for path in paths:
        for r in modern_rows(path):
            c = canonical_chrom(r["chrom"])
            split = "discovery" if c in background.discovery else "validation"
            for scope in ("all", split):
                key = (r["family"], r["boundary"], r["sensitivity"], int(r["threshold"]), scope)
                v = data[key][r["block"]][int(r["bin_start"])]
                n, m = int(r["n_sites"]), int(r["n_multi"])
                rate = background.rate(c, background.key(r), scope)
                v[2] += n
                if rate is None:
                    v[3] += n
                else:
                    v[0] += m
                    v[1] += n*rate
    return data


def modern_infer_profiles(data, args):
    all_blocks = sorted({b for profile in data.values() for b in profile})
    index = {b:i for i,b in enumerate(all_blocks)}
    results, inferential = [], []
    for key, blocks in sorted(data.items()):
        family, end, sensitivity, threshold, split = key
        flank = np.zeros((len(all_blocks), 2))
        flank_n = flank_missing = 0.
        for b, bins in blocks.items():
            for pos, v in bins.items():
                if abs(pos+args.step/2) >= args.flank*.75:
                    flank[index[b]] += v[:2]
                    flank_n += v[2]
                    flank_missing += v[3]
        flank_support = 1-flank_missing/flank_n if flank_n else 0.
        for start in range(-args.flank, args.flank-args.window+1, args.step):
            matrix = np.zeros((len(all_blocks), 4))
            matrix[:, 2:] = flank
            n = missing = 0.
            for b, bins in blocks.items():
                for pos in range(start, start+args.window, args.step):
                    v = bins.get(pos)
                    if v is not None:
                        matrix[index[b], :2] += v[:2]
                        n += v[2]
                        missing += v[3]
            support = 1-missing/n if n else 0.
            # Centre contrasts do not overlap the distal control bins.
            eligible_local = start >= -args.flank*.5 and start+args.window <= args.flank*.5
            global_fit = modern_cluster_estimate(matrix)
            local_fit = modern_cluster_estimate(matrix, True) if eligible_local else None
            enough = support >= args.min_support and matrix[:, 0].sum() >= 20
            local_enough = (flank_support >= args.min_support and matrix[:, 2].sum() >= 20)
            record = dict(family=family, boundary=end, sensitivity=sensitivity, threshold=threshold, split=split,
                          start_bp=start, end_bp=start+args.window, observed_multi=matrix[:, 0].sum(),
                          expected_multi=matrix[:, 1].sum(), eligible_sites=int(n), excluded_sites=int(missing),
                          supported_sites=int(n-missing),
                          supported_fraction=support, flank_supported_fraction=flank_support,
                          oe=None, ci_low=None, ci_high=None, local_ratio=None,
                          local_ci_low=None, local_ci_high=None, local_p=None,
                          discovery_maxT_p=None, n_blocks=0, status="insufficient_support")
            if global_fit:
                record.update(oe=global_fit["estimate"], ci_low=global_fit["ci_low"],
                              ci_high=global_fit["ci_high"], n_blocks=global_fit["n_blocks"])
            if local_fit and local_enough:
                record.update(local_ratio=local_fit["estimate"], local_ci_low=local_fit["ci_low"],
                              local_ci_high=local_fit["ci_high"], local_p=local_fit["pointwise_p"])
            if enough and global_fit and global_fit["n_blocks"] >= args.min_blocks:
                record["status"] = "complete"
            results.append(record)
            if (record["status"] == "complete" and local_fit and local_enough
                    and local_fit["n_blocks"] >= args.min_blocks and sensitivity == "strict"
                    and threshold == 5 and split == "discovery"):
                inferential.append((record, local_fit))
    if inferential:
        influences = np.column_stack([x[1]["influence"] for x in inferential])
        rng = np.random.default_rng(args.seed)
        maximum = []
        for begin in range(0, args.bootstrap, 100):
            multipliers = rng.choice([-1., 1.], size=(min(100, args.bootstrap-begin), len(all_blocks)))
            maximum.extend(np.max(np.abs(multipliers @ influences), axis=1))
        maximum = np.array(maximum)
        for record, fit in inferential:
            record["discovery_maxT_p"] = (1+int(np.sum(maximum >= abs(fit["z"]))))/(args.bootstrap+1)
    peaks = []
    for family in args.families.split(","):
        for end in ("5p", "3p"):
            choices = [r for r,_ in inferential if r["family"] == family and r["boundary"] == end]
            if not choices:
                continue
            # Freeze location on discovery chromosomes before looking at validation.
            best = max(choices, key=lambda r:r["local_ratio"])
            for threshold in (5, 10):
                matches = [r for r in results if r["family"] == family and r["boundary"] == end
                           and r["sensitivity"] == "strict" and r["threshold"] == threshold
                           and r["split"] == "validation" and r["start_bp"] == best["start_bp"]]
                val = matches[0] if matches else None
                peaks.append(dict(family=family, boundary=end, start_bp=best["start_bp"], end_bp=best["end_bp"],
                                  discovery_local_ratio=best["local_ratio"], discovery_maxT_p=best["discovery_maxT_p"],
                                  validation_threshold=threshold, validation_local_ratio=val["local_ratio"] if val else None,
                                  validation_ci_low=val["local_ci_low"] if val else None,
                                  validation_ci_high=val["local_ci_high"] if val else None,
                                  validation_p=val["local_p"] if val and val["status"] == "complete" else None,
                                  validation_fdr=None, status="exploratory_or_unconfirmed"))
    valid = [r for r in peaks if r["validation_threshold"] == 5 and r["validation_p"] is not None]
    for r, q in zip(valid, scan.bh_fdr([r["validation_p"] for r in valid])):
        r["validation_fdr"] = q
        if (r["discovery_maxT_p"] < .05 and q < .05
                and r["discovery_local_ratio"] > 1 and r["validation_local_ratio"] > 1):
            r["status"] = "held_out_chromosome_confirmation"
    return results, peaks


def modern_plot_profiles(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    families = sorted({r["family"] for r in results})
    if not families:
        return
    fig, axes = plt.subplots(len(families), 2, figsize=(11, 3.4*len(families)), squeeze=False)
    for i, family in enumerate(families):
        for j, end in enumerate(("5p", "3p")):
            ax = axes[i, j]
            for split, color in [("discovery", "#176B87"), ("validation", "#D16B38")]:
                rr = sorted([r for r in results if r["family"] == family and r["boundary"] == end
                             and r["sensitivity"] == "strict" and r["threshold"] == 5
                             and r["split"] == split and r["oe"] is not None], key=lambda r:r["start_bp"])
                if rr:
                    x = [(r["start_bp"]+r["end_bp"])/2 for r in rr]
                    ax.plot(x, [r["oe"] for r in rr], color=color, label=split)
                    ax.fill_between(x, [r["ci_low"] for r in rr], [r["ci_high"] for r in rr], color=color, alpha=.13)
            ax.axhline(1, ls="--", color=".6")
            ax.axvline(0, ls=":", color=".5")
            ax.set(title=f"{family} · annotated {end} edge", xlabel="Strand-oriented distance (bp)", ylabel="Multi-ALT O/E")
            if ax.get_legend_handles_labels()[0]:
                ax.legend(frameon=False)
    fig.suptitle("AC ≥ 5 · strict masks · fixed discovery / validation chromosomes", y=1.01)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"boundary_profiles.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def modern_assert_annotation_inputs(options, args):
    """Cached global denominators cannot be paired with a newly filtered scan."""
    mappings = [("fasta", "fasta"), ("erv_joined", "erv_joined"),
                ("segdup_bed", "architecture_segdup_bed"),
                ("lcr_bed", "architecture_lcr_bed"),
                ("low_mappability_bed", "architecture_low_mappability_bed")]
    for argument, option in mappings:
        override = getattr(args, argument, None)
        original = options.get(option)
        if override and (not original or Path(override).expanduser().resolve() != Path(original).expanduser().resolve()):
            raise ValueError(f"--{argument.replace('_', '-')} differs from the original modern scan. "
                             "Use the original input with its cached global denominator, or make a separate modern scan.")


def modern_verify_recorded_file(metadata):
    """Use the old run's metadata receipt; do not introduce full-file hashing."""
    if not isinstance(metadata, dict) or not {"path", "size", "mtime_ns"} <= set(metadata):
        return
    current = modern_file_state(metadata["path"])
    if current["bytes"] != metadata["size"] or current["mtime_ns"] != metadata["mtime_ns"]:
        raise ValueError(f"Original modern input changed since its cached global counts: {metadata['path']}")


def modern_targeted_inputs(args, allowed_chromosomes=None):
    """Resolve an adopted store without changing its scientific identity."""
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    run = scan.load_run(args.modern_dir)
    options = run["options"]
    for argument in ("fasta", "erv_joined", "segdup_bed", "lcr_bed", "low_mappability_bed"):
        if getattr(args, argument, None):
            setattr(args, argument, str(scan.resolve_run_path(run, getattr(args, argument))))
    modern_assert_annotation_inputs(options, args)
    receipt = run.get("identity", {})
    for name in ("fasta", "fasta_index", "fasta_gzi", "callable_bed"):
        if receipt.get(name):
            scan.verify_recorded_input(run, receipt[name])
    annotation_receipts = (options.get("architecture_annotations") or {}).get("files", {})
    for name in ("erv_joined", "architecture_segdup_bed", "architecture_lcr_bed", "architecture_low_mappability_bed", "EUR"):
        if annotation_receipts.get(name):
            scan.verify_recorded_input(run, annotation_receipts[name])
    if str(options.get("build", "38")) not in {"38", "GRCh38", "hg38"}:
        raise ValueError("HERV targeted analysis requires GRCh38 modern input")
    if options.get("region"):
        raise ValueError("Use completed chromosome/genome shards, not a region-restricted modern scan")
    annotation = options.get("architecture_annotation_cache")
    if not annotation:
        raise ValueError("Completed scan lacks its annotation cache; restore the cache rather than rebuild frozen inputs")
    annotation = Path(scan.resolve_run_path(run, annotation))
    scan.validate_existing_joined_cache(args.erv_joined or options["erv_joined"], annotation)
    allowed = set(allowed_chromosomes) if allowed_chromosomes is not None else MODERN_AUTOSOMES
    selected = {canonical_chrom(c.strip()) for c in args.chromosomes.split(",")} if args.chromosomes else allowed
    if not selected <= allowed:
        raise ValueError(f"Unsupported targeted chromosomes: {sorted(selected - allowed)}")
    receipts = {str(Path(scan.resolve_run_path(run, v["path"])).resolve()): v
                for v in receipt.get("inputs", []) if isinstance(v, dict) and "path" in v}
    tasks, used, original = [], set(), set()
    for shard in run["shards"]:
        info = json.loads((Path(shard)/"scan.json").read_text())
        for c in info["scope"]:
            c = canonical_chrom(c)
            if c not in allowed:
                continue
            original.add(c)
            if c not in selected:
                continue
            if c in used:
                raise ValueError("Multiple original shards cover one chromosome; use non-overlapping completed chromosome shards")
            used.add(c)
            source = scan.resolve_run_path(run, info["input"])
            metadata = receipts.get(str(Path(source).resolve()))
            if metadata:
                scan.verify_recorded_input(run, metadata)
            tasks.append((c, source))
    if not tasks:
        raise ValueError("No matching completed chromosome shards")
    if args.chromosomes and selected-used:
        raise ValueError(f"Requested chromosomes absent from completed original shards: {sorted(selected-used)}")
    return run, annotation, tasks, original


@transactional("windows_manifest.json")
def modern_run_windows(args):
    if args.step <= 0 or args.window < args.step or args.window % args.step or args.flank % args.step:
        raise ValueError("Positive window/flank must be integer multiples of step")
    if args.window > args.flank or args.block_bp < 100000 or args.bootstrap < 200:
        raise ValueError("window <= flank, block >=100kb and >=200 multiplier draws required")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    run, annotation, tasks, _ = modern_targeted_inputs(args)
    options = run["options"]
    files, used = [], {c for c, _ in tasks}
    for chrom, source in tasks:
        print(f"[boundaries] indexed chr{chrom} query", flush=True)
        files.append(modern_scan_boundary_chrom(source, chrom, annotation/f"chr{chrom}.tsv.gz", out/"cache", options, args))
    background = modern_Background(modern_architecture_paths(args.modern_dir, "global"), args.min_training_sites,
                            args.discovery_chromosomes.split(","))
    data = modern_profile_counts(files, background, args)
    results, peaks = modern_infer_profiles(data, args)
    modern_write_rows(out/"boundary_profiles.tsv", results)
    modern_write_rows(out/"frozen_peak_validation.tsv", peaks)
    modern_plot_profiles(results, out)
    manifest = dict(schema=MODERN_SCHEMA, status="complete", chromosomes=sorted(used, key=int),
                    families=args.families, window=args.window, step=args.step, flank=args.flank,
                    discovery_chromosomes=sorted(background.discovery, key=int),
                    validation_background="frozen discovery chromosomes only",
                    inference=f"{args.block_bp:g}bp physical-block sandwich CI; Rademacher max-|Z| over strict AC5 discovery family/end/central windows; held-out AC5 BH confirmation",
                    local_control=f"outer 25% of +/-{args.flank} bp; local_ratio = window O/E divided by distal-flank O/E",
                    min_support=args.min_support, minimum_training_sites=args.min_training_sites,
                    caveats=["Observed joined alignment edges are not TSS or intact provirus boundaries",
                             "Nearest edge per family/end; no repeated site in one window",
                             "Coincident same-family/end edges with conflicting strands are excluded from oriented assignment",
                             "Original SNV denominator, not callable-bp mutation rate",
                             "Background offsets treated as fixed; physical blocks approximate dependence",
                             "Local GC, repeat age and recombination not included in this O/E",
                             "Chromosome split is internal validation, not an independent cohort"])
    dump_json(out/"windows_manifest.json", manifest)
    return manifest


MODERN_SPATIAL_SCHEMA = "herv-spatial-1.0"
MODERN_PARTITIONS = ("head", "body", "tail", "short_element", "ambiguous_overlap")


class modern_AlignedBodies:
    """One assignment per family/site; no assumption of an intact provirus.

    Distance follows the strand along the concatenated *aligned* bases of a
    joined RepeatMasker record. Gaps do not acquire body positions. Conflicting
    overlapping annotations retain their counts in an explicit residual.
    """
    def __init__(self, elements, families, end_bp=200, bins=20):
        if end_bp < 1 or bins < 2:
            raise ValueError("Positive end-bp and at least two body bins required")
        self.families = set(families)
        self.end_bp, self.bins = end_bp, bins
        self.elements = [e for e in elements if e.rep_class == "LTR" and e.family in self.families]
        self.geometry = {}
        for e in self.elements:
            segments = tuple(scan.merge_intervals(list(e.segments or ((e.start, e.end),))))
            self.geometry[e] = (segments, sum(hi-lo for lo, hi in segments))

    def intervals(self):
        return scan.merge_intervals([segment for segments, _ in self.geometry.values() for segment in segments])

    def coordinate(self, e, pos1):
        segments, length = self.geometry[e]
        prefix, p = 0, int(pos1)-1
        for lo, hi in segments:
            if lo <= p < hi:
                # Twice the base-centre distance is integer and exactly mirrors
                # on the reverse strand, including at normalized bin boundaries.
                distance2 = 2*prefix + 2*(p-lo) + 1
                if e.strand == "-":
                    distance2 = 2*length-distance2
                if length <= 2*self.end_bp:
                    partition = "short_element"
                elif distance2 < 2*self.end_bp:
                    partition = "head"
                elif distance2 >= 2*(length-self.end_bp):
                    partition = "tail"
                else:
                    partition = "body"
                return partition, distance2*self.bins//(2*length)
            prefix += hi-lo
        raise ValueError("A body coordinate was requested for an unaligned base")

    def assign(self, overlaps, pos1):
        by_family = defaultdict(dict)
        for e in overlaps:
            if e.rep_class != "LTR" or e.family not in self.families:
                continue
            segments, _ = self.geometry[e]
            # Name/source aliases at identical geometry do not create replicates.
            geometry = e.strand, segments
            by_family[e.family][geometry] = self.coordinate(e, pos1)
        result = {}
        for family, geometries in by_family.items():
            strands = {g[0] for g in geometries}
            parts = {v[0] for v in geometries.values()}
            bins = {v[1] for v in geometries.values()}
            part = next(iter(parts)) if len(strands) == len(parts) == 1 else "ambiguous_overlap"
            body_bin = next(iter(bins)) if len(strands) == len(bins) == 1 else -1
            overlap = ("single_geometry" if len(geometries) == 1 else
                       "compatible_overlap" if part != "ambiguous_overlap" and body_bin >= 0 else "ambiguous_overlap")
            result[family] = part, body_bin, overlap
        return result


def modern_spatial_annotation_audit(elements, families):
    """Report exact labels, gaps and overlap, without merging biological aliases."""
    output = []
    groups = defaultdict(list)
    for e in elements:
        if e.rep_class == "LTR":
            groups[e.chrom, e.family].append(e)
    for (chrom, family), ee in sorted(groups.items()):
        geometries = [(e.strand, tuple(scan.merge_intervals(list(e.segments or ((e.start, e.end),))))) for e in ee]
        lengths = [sum(t-s for s, t in segments) for _, segments in geometries]
        aligned = sum(t-s for s, t in scan.merge_intervals([p for _, segments in geometries for p in segments]))
        spans = sum(t-s for s, t in scan.merge_intervals([(e.start, e.end) for e in ee]))
        output.append(dict(chrom=chrom, rep_class="LTR", raw_family=family,
                           counter_label="family:"+(family or "unspecified"), selected_exact=family in families,
                           records=len(ee), identical_geometry_extra_records=len(ee)-len(set(geometries)),
                           aligned_record_sum_bp=sum(lengths), aligned_union_bp=aligned,
                           envelope_union_bp=spans, envelope_only_gap_bp=spans-aligned,
                           overlap_multiplicity_extra_bp=sum(lengths)-aligned,
                           min_aligned_bp=min(lengths), median_aligned_bp=float(np.median(lengths)),
                           max_aligned_bp=max(lengths)))
    return output


def modern_scan_spatial_chrom(vcf, chrom, annotation, out, options, args):
    """One indexed query over the family-aligned union, never a genome rescan."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path, marker = out/f"chr{chrom}.spatial.tsv.gz", out/f"chr{chrom}.spatial.json"
    fasta = str(args.fasta or options["fasta"])
    mask_paths = {n: options.get(k) for n, k in (("segdup", "architecture_segdup_bed"),
                  ("lcr", "architecture_lcr_bed"), ("low_mappability", "architecture_low_mappability_bed"))}
    for name in mask_paths:
        mask_paths[name] = getattr(args, name+"_bed", None) or mask_paths[name]
    ld_path = (Path(options["ld_block_dir"])/f"EUR.{options.get('architecture_build', '38')}.bed"
               if options.get("ld_block_dir") else None)
    inputs = [vcf, annotation, fasta, fasta+".fai", *[p for p in mask_paths.values() if p]]
    if ld_path:
        inputs.append(ld_path)
    if options.get("callable_bed"):
        inputs.append(options["callable_bed"])
    for suffix in (".tbi", ".csi"):
        if Path(str(vcf)+suffix).is_file():
            inputs.append(str(vcf)+suffix)
            break
    else:
        raise ValueError(f"Indexed VCF required: {vcf}")
    identity = dict(schema=MODERN_SPATIAL_SCHEMA, inputs=[modern_file_state(p) for p in inputs], source_options=options,
                    code=[modern_file_state(p) for p in (Path(__file__), Path(scan.__file__), Path(__file__).with_name("0.common.py"))],
                    families=args.families, end_bp=args.end_bp, body_bins=args.body_bins, block_bp=args.block_bp)
    if marker.exists():
        old = json.loads(marker.read_text())
        if old["identity"] != identity:
            raise ValueError(f"Spatial inputs/settings/code changed; use a fresh output: {out}")
        expected = {Path(p).name:state for p, state in old.get("outputs", {}).items()}
        required = [path, out/f"chr{chrom}.annotation_audit.tsv", out/f"chr{chrom}.elements.tsv.gz"]
        if all(p.is_file() and expected.get(p.name) == modern_file_state(p) for p in required):
            print(f"[resume] spatial chr{chrom}", flush=True)
            return path
        raise ValueError(f"Completed spatial cache changed: {path}")
    header = modern_read_indexed_header(vcf, args.bcftools)
    scan.validate_count_header(header)
    fields = ["AC", "AN"]+[f for f in ("AF", "AS_FilterStatus", "lcr", "segdup") if f in header["info"]]
    elements = scan.parse_joined_rmsk(annotation, chromosome=chrom)
    families = args.families.split(",")
    body = modern_AlignedBodies(elements, families, args.end_bp, args.body_bins)
    repeat_sweep = scan.ERVOverlapSweep(elements)
    masks = {n: IntervalMask(p, [chrom]) for n, p in mask_paths.items() if p}
    strict_complete = all(n in masks or n in header["info"] for n in ("lcr", "segdup")) and "low_mappability" in masks
    callable_mask = IntervalMask(options["callable_bed"], [chrom]) if options.get("callable_bed") else None
    eur = scan.BlockSweep(scan.read_blocks(ld_path, "EUR", options.get("architecture_build", "38"))[0]) if ld_path else None
    intervals = body.intervals()
    raw_chroms = [c for c in header["contigs"] if canonical_chrom(c) == chrom]
    if len(raw_chroms) != 1:
        raise ValueError(f"No unique VCF contig for chr{chrom}")
    bed = out/f"chr{chrom}.spatial.regions.bed"
    with bed.open("w") as handle:
        for lo, hi in intervals:
            handle.write(f"{raw_chroms[0]}\t{lo}\t{hi}\n")
    columns = ["chrom", "block", "eur_partition", "family", "partition", "body_bin", "overlap_assignment",
               "sensitivity", "threshold", "context", "af_bin", "an_bin", "n_sites", "n_multi"]
    counts, qc, last_block = defaultdict(lambda: [0, 0]), Counter(), None
    partial = Path(str(path)+".partial")
    with FastaContext(fasta) as fa, gzip.open(partial, "wt", compresslevel=2, newline="") as target:
        writer = csv.writer(target, delimiter="\t")
        writer.writerow(columns)
        def flush():
            writer.writerows([*k, *v] for k, v in sorted(counts.items()))
            counts.clear()
        records = modern_query_regions(vcf, fields, bed, args.bcftools) if intervals else iter(())
        clock = time.monotonic()
        for (c, pos), group in itertools.groupby(records, key=lambda r:(canonical_chrom(r["chrom"]), r["pos"])):
            if c != chrom:
                raise ValueError("Unexpected chromosome from targeted spatial query")
            if callable_mask and not callable_mask.contains(c, pos):
                qc["outside_callable_mask"] += 1
                continue
            overlaps = repeat_sweep.locate(c, pos)
            assignment = body.assign(overlaps, pos)
            if not assignment:
                raise ValueError("Indexed body query returned a position outside the selected aligned family union")
            ref, alts, _, flags = scan.group_alleles(group, header, [], options, qc, None, [])
            observed = sorted((a for a in alts.values() if a["status"] == "observed_pass"),
                              key=lambda a:(-a["af"], -a["ac"], a["alt"]))
            if not observed:
                continue
            ctx = fa.context(c, pos).upper()
            if len(ctx) != 3 or ctx[1] != ref:
                raise ValueError(f"REF/FASTA mismatch: {c}:{pos}")
            context, _ = scan.canonical_context(ctx, ref)
            if any(b not in "ACGT" for b in context):
                qc["ambiguous_context_retained"] += 1
            primary = observed[0]
            strata = context, af_bin(primary["af"]), an_bin(primary["an"])
            sensitivities = ["all"]
            if strict_complete and not modern_strict_site_excluded(c, pos, flags, overlaps, masks):
                sensitivities.append("strict")
            if not ((ref == "C" and ctx[2] == "G") or (ref == "G" and ctx[0] == "C")):
                sensitivities += [s+"_no_cpg" for s in list(sensitivities)]
            covered = "inside" if eur and eur.locate(c, pos) else "outside"
            block = f"{c}:{(pos-1)//args.block_bp}"
            if last_block is not None and last_block != block:
                flush()
            last_block = block
            multi = {t: int(sum(a["ac"] >= t for a in observed) >= 2) for t in (5, 10)}
            for family, (part, body_bin, overlap) in assignment.items():
                qc[f"{family}:{overlap}"] += 1
                for sensitivity in sensitivities:
                    for threshold in (5, 10):
                        key = (c, block, covered, family, part, body_bin, overlap, sensitivity, threshold, *strata)
                        counts[key][0] += 1
                        counts[key][1] += multi[threshold]
            qc["eligible_sites_in_union"] += 1
            if time.monotonic()-clock > 60:
                print(f"[spatial chr{chrom}] {qc['eligible_sites_in_union']:,} eligible sites", flush=True)
                clock = time.monotonic()
        flush()
    os.replace(partial, path)
    audit_path, element_path = out/f"chr{chrom}.annotation_audit.tsv", out/f"chr{chrom}.elements.tsv.gz"
    modern_write_rows(audit_path, modern_spatial_annotation_audit(elements, families))
    modern_write_rows(element_path, [dict(chrom=e.chrom, start=e.start, end=e.end, rep_name=e.name,
                                   rep_class=e.rep_class, rep_family=e.family, strand=e.strand, source_id=e.source_id,
                                   aligned_bp=body.geometry[e][1], envelope_bp=e.end-e.start,
                                   unaligned_gap_bp=e.end-e.start-body.geometry[e][1],
                                   name_ends_int=e.name.endswith("-int"),
                                   short_element=body.geometry[e][1] <= 2*args.end_bp,
                                   aligned_segments=",".join(f"{lo}-{hi}" for lo, hi in body.geometry[e][0])) for e in body.elements],
               fields=["chrom", "start", "end", "rep_name", "rep_class", "rep_family", "strand", "source_id",
                       "aligned_bp", "envelope_bp", "unaligned_gap_bp", "name_ends_int", "short_element", "aligned_segments"])
    dump_json(marker, dict(identity=identity, outputs={str(p):modern_file_state(p) for p in (path, audit_path, element_path)},
                          qc=dict(qc), strict_masks_complete=strict_complete, eur_ld_partition_available=eur is not None,
                          queried_aligned_bp=sum(hi-lo for lo, hi in intervals)))
    return path


def modern_spatial_counts(paths, background, with_eur=True, collect_raw=True):
    """Apply the same LOCO rate to every disjoint part of each stratum.

    Values are supported observed/expected/sites and unsupported observed/sites.
    Unsupported loci are never left in the observed numerator. Totals are
    accumulated independently of either partitioning, for accounting checks.
    """
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: np.zeros(5))))
    raw = defaultdict(lambda: np.zeros(2, dtype=np.int64))
    for path in paths:
        for r in modern_rows(path):
            c, n, m = canonical_chrom(r["chrom"]), int(r["n_sites"]), int(r["n_multi"])
            rate = background.rate(c, background.key(r))
            value = np.array([m, n*rate, n, 0, 0]) if rate is not None else np.array([0, 0, 0, m, n])
            scopes = ["all_family_sites"]
            if with_eur and r["eur_partition"] == "inside":
                scopes.append("eur_ld_partition")
            elif with_eur:
                # Keep an explicitly empty legacy-covered profile rather than
                # silently substituting an all-sites result for an absent map.
                prefix = r["family"], r["sensitivity"], int(r["threshold"]), "eur_ld_partition"
                data[(*prefix, "total")]["total"]
            for scope in scopes:
                prefix = r["family"], r["sensitivity"], int(r["threshold"]), scope
                for view, category in (("partition", r["partition"]), ("normalized", int(r["body_bin"])), ("total", "total")):
                    data[(*prefix, view)][category][r["block"]] += value
                if collect_raw and r["sensitivity"] in {"all", "strict"}:
                    raw[(c, *prefix, r["context"], r["af_bin"], r["an_bin"])] += [n, m]
    return data, raw


def modern_spatial_record(blocks, args):
    vals = np.array(list(blocks.values())) if blocks else np.zeros((0, 5))
    o, e, n, removed_o, removed_n = vals.sum(axis=0)
    eligible = n+removed_n
    fit = modern_cluster_estimate(vals[:, :2]) if e > 0 and o > 0 else None
    count = int(np.sum(vals[:, 2] > 0))
    support = n/eligible if eligible else 0.
    enough = bool(support >= args.min_support and count >= args.min_blocks and o >= 20 and fit)
    return dict(observed_multi=float(o), expected_multi=float(e), oe=float(o/e) if e > 0 else None,
                ci_low=fit["ci_low"] if enough else None, ci_high=fit["ci_high"] if enough else None,
                supported_sites=int(n), excluded_sites=int(removed_n), excluded_multi=int(removed_o),
                eligible_sites=int(eligible), raw_multi=int(o+removed_o), supported_fraction=float(support),
                n_blocks=count, status="complete_descriptive" if enough else "no_sites" if eligible == 0 else "low_support_or_few_events_blocks")


def modern_summarize_spatial(data, args):
    partitions, profiles, identities, contrasts = [], [], [], []
    prefixes = sorted({key[:4] for key in data})
    for prefix in prefixes:
        family, sensitivity, threshold, scope = prefix
        dimensions = dict(family=family, sensitivity=sensitivity, threshold=threshold, scope=scope)
        total_blocks = data.get((*prefix, "total"), {}).get("total", {})
        total = modern_spatial_record(total_blocks, args)
        total_e = total["expected_multi"]
        for view, categories, target in (("partition", MODERN_PARTITIONS, partitions),
                                          ("normalized", tuple(range(args.body_bins))+(-1,), profiles)):
            records = []
            for category in categories:
                blocks = data.get((*prefix, view), {}).get(category, {})
                record = dimensions | modern_spatial_record(blocks, args)
                if view == "partition":
                    record["partition"] = category
                else:
                    record.update(body_bin=category,
                                  start_fraction=category/args.body_bins if category >= 0 else None,
                                  end_fraction=(category+1)/args.body_bins if category >= 0 else None,
                                  position_label="ambiguous_overlap" if category < 0 else "aligned_body")
                record.update(expected_weight=record["expected_multi"]/total_e if total_e else None,
                              oe_contribution=record["observed_multi"]/total_e if total_e else None,
                              excess_contribution=(record["observed_multi"]-record["expected_multi"])/total_e if total_e else None)
                records.append(record)
            target.extend(records)
            # Verify both supported and unsupported accounting, not just ratios.
            fields = ("observed_multi", "expected_multi", "supported_sites", "excluded_sites", "excluded_multi")
            deltas = {field:sum(r[field] for r in records)-total[field] for field in fields}
            if any(not math.isclose(delta, 0., rel_tol=0., abs_tol=max(1e-8, abs(total[field])*1e-10))
                   for field, delta in deltas.items()):
                raise ValueError(f"Spatial accounting identity failed for {prefix}/{view}: {deltas}")
            weighted = sum(r["oe_contribution"] for r in records) if total_e else None
            ratios = [r["oe"] for r in records if r["oe"] is not None]
            identities.append(dimensions | dict(view=view, overall_oe=total["oe"],
                              expected_weighted_oe=weighted,
                              unweighted_mean_of_available_part_oe=float(np.mean(ratios)) if ratios else None,
                              observed_delta=deltas["observed_multi"], expected_delta=deltas["expected_multi"],
                              supported_sites_delta=deltas["supported_sites"], excluded_sites_delta=deltas["excluded_sites"],
                              excluded_multi_delta=deltas["excluded_multi"],
                              summed_excess_contribution=sum(r["excess_contribution"] for r in records) if total_e else None,
                              status="exact_accounting"))
            if view == "partition":
                target.append(dimensions | total | dict(partition="total", expected_weight=1. if total_e else None,
                              oe_contribution=total["oe"], excess_contribution=total["oe"]-1 if total_e else None))
        # Prespecified comparisons among the unambiguous long-record parts.
        part_data = data.get((*prefix, "partition"), {})
        block_names = sorted(total_blocks)
        for part in ("head", "tail"):
            matrix = np.zeros((len(block_names), 4))
            for i, block in enumerate(block_names):
                matrix[i, :2] = part_data.get(part, {}).get(block, np.zeros(5))[:2]
                matrix[i, 2:] = part_data.get("body", {}).get(block, np.zeros(5))[:2]
            fit = modern_cluster_estimate(matrix, local=True)
            a = modern_spatial_record(part_data.get(part, {}), args)
            b = modern_spatial_record(part_data.get("body", {}), args)
            valid = fit and a["status"] == b["status"] == "complete_descriptive"
            contrasts.append(dimensions | dict(contrast=part+"_vs_body", ratio_of_oe=fit["estimate"] if fit else None,
                             ci_low=fit["ci_low"] if valid else None, ci_high=fit["ci_high"] if valid else None,
                             pointwise_p=fit["pointwise_p"] if valid else None, primary_fdr=None,
                             n_blocks=fit["n_blocks"] if fit else 0,
                             status="prespecified_association" if valid else "insufficient_support"))
    primary_scope = "eur_ld_partition" if any(r["scope"] == "eur_ld_partition" for r in contrasts) else "all_family_sites"
    tested = [r for r in contrasts if r["scope"] == primary_scope and r["sensitivity"] == "strict"
              and r["threshold"] == 5 and r["pointwise_p"] is not None]
    for r, q in zip(tested, scan.bh_fdr([r["pointwise_p"] for r in tested])):
        r["primary_fdr"] = q
    return partitions, profiles, identities, contrasts


def modern_primary_spatial_record(blocks, args, total_without_loco):
    record = modern_spatial_record(blocks, args)
    supported_m = record["observed_multi"]
    missing_n, missing_m = record["excluded_sites"], record["excluded_multi"]
    partial_e = record["expected_multi"]
    record.update(model="primary_full_scope", observed_multi=float(record["raw_multi"]),
                  supported_observed_multi=supported_m, partial_expected_multi=partial_e,
                  excluded_sites=0, excluded_multi=0, n_without_loco=missing_n,
                  multi_without_loco=missing_m, total_without_loco=int(total_without_loco))
    if total_without_loco:
        # Even a supported head cannot be presented as a weighted share of an
        # unavailable full-scope total. Raw counts and partial E remain auditable.
        record.update(expected_multi=None, oe=None, ci_low=None, ci_high=None,
                      status="unavailable_primary_loco")
    return record


def modern_summarize_primary_spatial(data, args, scoped_exposure=None):
    """Decompose the primary full-scope estimator without substituting common support."""
    partitions, profiles, identities = [], [], []
    prefixes = sorted({key[:4] for key in data if key[1] in {"all", "strict"}
                       and key[3] == "eur_ld_partition"})
    for prefix in prefixes:
        family, sensitivity, threshold, scope = prefix
        dimensions = dict(family=family, sensitivity=sensitivity, threshold=threshold, scope=scope)
        total_blocks = data.get((*prefix, "total"), {}).get("total", {})
        missing = sum(v[4] for v in total_blocks.values())
        if scoped_exposure is not None:
            # A chromosome subset cannot manufacture a decomposition of an
            # unavailable full-scope total by omitting its unsupported strata.
            missing = max(missing, scoped_exposure.get(prefix[:3], np.zeros(5))[4])
        total = modern_primary_spatial_record(total_blocks, args, missing)
        total_e = total["expected_multi"]
        for view, categories, target in (("partition", MODERN_PARTITIONS, partitions),
                                          ("normalized", tuple(range(args.body_bins))+(-1,), profiles)):
            records = []
            for category in categories:
                blocks = data.get((*prefix, view), {}).get(category, {})
                record = dimensions | modern_primary_spatial_record(blocks, args, missing)
                if view == "partition":
                    record["partition"] = category
                else:
                    record.update(body_bin=category, start_fraction=category/args.body_bins if category >= 0 else None,
                                  end_fraction=(category+1)/args.body_bins if category >= 0 else None,
                                  position_label="ambiguous_overlap" if category < 0 else "aligned_body")
                record.update(expected_weight=record["expected_multi"]/total_e if total_e else None,
                              oe_contribution=record["observed_multi"]/total_e if total_e else None,
                              excess_contribution=(record["observed_multi"]-record["expected_multi"])/total_e if total_e else None)
                records.append(record)
            target.extend(records)
            fields = ("observed_multi", "eligible_sites", "n_without_loco", "multi_without_loco", "partial_expected_multi")
            deltas = {field:sum(r[field] for r in records)-total[field] for field in fields}
            if any(not math.isclose(delta, 0., rel_tol=0., abs_tol=max(1e-8, abs(total[field])*1e-10))
                   for field, delta in deltas.items()):
                raise ValueError(f"Primary spatial accounting failed for {prefix}/{view}: {deltas}")
            identities.append(dimensions | dict(model="primary_full_scope", view=view,
                              overall_oe=total["oe"], expected_weighted_oe=sum(r["oe_contribution"] for r in records) if total_e else None,
                              observed_delta=deltas["observed_multi"], eligible_sites_delta=deltas["eligible_sites"],
                              expected_delta=sum(r["expected_multi"] for r in records)-total_e if total_e is not None else None,
                              partial_expected_delta=deltas["partial_expected_multi"],
                              n_without_loco_delta=deltas["n_without_loco"], total_without_loco=int(missing),
                              summed_excess_contribution=sum(r["excess_contribution"] for r in records) if total_e else None,
                              status="raw_accounting_only_primary_loco_unavailable" if missing else "exact_accounting"))
            if view == "partition":
                target.append(dimensions | total | dict(partition="total", expected_weight=1. if total_e else None,
                              oe_contribution=total["oe"], excess_contribution=total["oe"]-1 if total_e else None))
    return partitions, profiles, identities


def modern_reconcile_spatial_counters(modern, raw, background, chromosomes, families, with_eur=True,
                             comparison_scope="complete_original_autosomes", primary_background=None):
    """Compare raw N/M at every old family stratum before interpreting O/E.

    The old family counters include all aligned body sites, while the old formal
    exposure summary excludes sites outside the EUR partition. Both are checked.
    Expected counts are recomputed under the new common-support settings; this
    does not relabel the primary full-scope O/E as a newly reproduced estimate.
    """
    legacy = defaultdict(lambda: np.zeros(2, dtype=np.int64))
    labels = Counter()
    selected_labels = {"family:"+family for family in families}
    sources = [("erv", "all_family_sites")]
    if with_eur:
        sources.append(("exposure", "eur_ld_partition"))
    for kind, scope in sources:
        for path in modern_architecture_paths(modern, kind):
            for r in modern_rows(path):
                if kind == "exposure" and primary_background is not None:
                    primary_background.observe_exposure(r)
                c = canonical_chrom(r["chrom"])
                if c not in chromosomes or r["population_mode"] != "all":
                    continue
                if kind == "erv" and r["label"].startswith("family:") and r["threshold"] == "5" and r["sensitivity"] == "all":
                    labels[r["label"]] += int(r["n_sites"])
                if r["label"] not in selected_labels or r["threshold"] not in {"5", "10"}:
                    continue
                if r["sensitivity"] not in {"all", "exclude_combined"}:
                    continue
                if kind == "exposure" and r["block_id"] == "outside_partition":
                    continue
                sensitivity = "strict" if r["sensitivity"] == "exclude_combined" else "all"
                key = (c, r["label"][7:], sensitivity, int(r["threshold"]), scope, r["context"], r["af_bin"], r["an_bin"])
                legacy[key] += [int(r["n_sites"]), int(r["n_multi"])]
    sums = defaultdict(lambda: np.zeros((2, 7)))
    mismatched = Counter()
    strata_counts = Counter()
    mismatch_rows = []
    for key in sorted(set(raw) | set(legacy)):
        c, family, sensitivity, threshold, scope, context, af, an = key
        group = family, sensitivity, threshold, scope
        old, new = legacy.get(key, np.zeros(2)), raw.get(key, np.zeros(2))
        strata_counts[group] += 1
        if not np.array_equal(old, new):
            mismatched[group] += 1
            mismatch_rows.append(dict(chrom=c, family=family, sensitivity=sensitivity, threshold=threshold, scope=scope,
                                      context=context, af_bin=af, an_bin=an, counter_sites=int(old[0]), target_sites=int(new[0]),
                                      counter_multi=int(old[1]), target_multi=int(new[1])))
        rate = background.rate(c, background.key(dict(context=context, af_bin=af, an_bin=an,
                                                      sensitivity=sensitivity, threshold=threshold)))
        for i, (n, m) in enumerate((old, new)):
            sums[group][i] += ([n, m, n, m, n*rate, 0, 0] if rate is not None else [n, m, 0, 0, 0, n, m])
    result = []
    for (family, sensitivity, threshold, scope), values in sorted(sums.items()):
        group = family, sensitivity, threshold, scope
        old, new = values
        result.append(dict(family=family, sensitivity=sensitivity, threshold=threshold, scope=scope,
                           comparison_scope=comparison_scope, compared_chromosomes=",".join(sorted(chromosomes, key=int)),
                           counter_raw_sites=int(old[0]), targeted_raw_sites=int(new[0]), sites_delta=int(new[0]-old[0]),
                           counter_raw_multi=int(old[1]), targeted_raw_multi=int(new[1]), multi_delta=int(new[1]-old[1]),
                           counter_supported_sites=int(old[2]), targeted_supported_sites=int(new[2]),
                           counter_supported_multi=int(old[3]), targeted_supported_multi=int(new[3]),
                           counter_recomputed_expected=float(old[4]), targeted_expected=float(new[4]),
                           expected_delta=float(new[4]-old[4]), excluded_sites=int(new[5]), excluded_multi=int(new[6]),
                           compared_strata=strata_counts[group], mismatched_strata=mismatched[group],
                           status="matched" if not mismatched[group] else "mismatch"))
    label_rows = [dict(counter_label=label, exact_family=label[7:], selected_exact=label[7:] in families,
                       eligible_sites_ac5_all=count, alias_policy="exact family strings only; no alias collapsing")
                  for label, count in sorted(labels.items())]
    return result, mismatch_rows, label_rows


def modern_optional_number(value):
    if value is None or str(value).strip().lower() in {"", ".", "na", "none", "nan"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite original enrichment value: {value}")
    return result


def modern_plot_spatial(partitions, profiles, out, end_bp, prefix="", model_label="Common-support model"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scope = "eur_ld_partition" if any(r["scope"] == "eur_ld_partition" for r in partitions) else "all_family_sites"
    primary = [r for r in partitions if r["scope"] == scope and r["sensitivity"] == "strict" and r["threshold"] == 5]
    families = sorted({r["family"] for r in primary})
    if not families:
        return
    fig, axes = plt.subplots(len(families), 2, figsize=(12, 3.8*len(families)), squeeze=False)
    categories = (*MODERN_PARTITIONS, "total")
    names = [f"Head\n{end_bp} bp", "Body", f"Tail\n{end_bp} bp", "Short", "Ambiguous", "Total"]
    for i, family in enumerate(families):
        rr = {r["partition"]:r for r in primary if r["family"] == family}
        for j, part in enumerate(categories):
            r = rr.get(part)
            if not r or r["oe"] is None:
                continue
            color = "#176B87" if part != "total" else ".2"
            axes[i, 0].plot(j, r["oe"], "o", color=color,
                            markerfacecolor=color if r["status"] == "complete_descriptive" else "white")
            if r["ci_low"] is not None:
                axes[i, 0].plot([j, j], [r["ci_low"], r["ci_high"]], color="#176B87")
            axes[i, 1].bar(j, r["excess_contribution"], color="#D16B38" if part != "total" else ".3")
        axes[i, 0].axhline(1, ls="--", color=".6")
        axes[i, 1].axhline(0, color=".6", lw=.8)
        axes[i, 0].set(title=f"{family}: disjoint aligned-body parts", ylabel="Multi-ALT observed / expected")
        axes[i, 1].set(title=f"{family}: weighted contribution to total excess",
                       ylabel="Contribution to overall O/E − 1\n(Observed part − expected part) / expected total")
        if rr.get("total", {}).get("status") == "unavailable_primary_loco":
            for ax in axes[i]:
                ax.text(.5, .5, "Original total unavailable:\nexposed strata lack LOCO training", transform=ax.transAxes, ha="center")
        for ax in axes[i]:
            ax.set(xticks=range(len(categories)), xticklabels=names)
    fig.suptitle(f"{model_label}\nAC ≥ 5 · strict masks · {scope}\nShort ≤ {2*end_bp} aligned bp; physical-block CI; open markers: low support / few events or blocks")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"{prefix}spatial_partitions.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    fig, axes = plt.subplots(len(families), 1, figsize=(9, 3.5*len(families)), squeeze=False)
    for i, family in enumerate(families):
        ax = axes[i, 0]
        for threshold, color in ((5, "#176B87"), (10, "#D16B38")):
            rr = sorted([r for r in profiles if r["family"] == family and r["scope"] == scope
                         and r["sensitivity"] == "strict" and r["threshold"] == threshold and r["body_bin"] >= 0],
                        key=lambda r:r["body_bin"])
            if not rr:
                continue
            x = [(r["start_fraction"]+r["end_fraction"])/2 for r in rr]
            y = [r["oe"] if r["oe"] is not None else np.nan for r in rr]
            low = [r["ci_low"] if r["ci_low"] is not None else np.nan for r in rr]
            high = [r["ci_high"] if r["ci_high"] is not None else np.nan for r in rr]
            ax.plot(x, y, "-", color=color, label=f"AC ≥ {threshold}")
            for xx, yy, r in zip(x, y, rr):
                ax.plot(xx, yy, "o", ms=3, color=color,
                        markerfacecolor=color if r["status"] == "complete_descriptive" else "white")
            ax.fill_between(x, low, high, color=color, alpha=.12)
        ax.axhline(1, ls="--", color=".6")
        ax.set(xlim=(0, 1), title=family, xlabel="Strand-oriented position along concatenated aligned bases (5′ → 3′)",
               ylabel="Multi-ALT observed / expected")
        ax.legend(frameon=False)
        if all(r["oe"] is None for r in profiles if r["family"] == family and r["scope"] == scope and r["sensitivity"] == "strict"):
            ax.text(.5, .5, "Original O/E unavailable or no estimable bins", transform=ax.transAxes, ha="center")
    fig.suptitle(f"{model_label}\nWhole joined-repeat profile · short fragments included · gaps excluded\nAmbiguous assignments retained in table; pointwise physical-block CI\nOpen markers: low support / few events or blocks")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"{prefix}whole_element_profiles.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def modern_merge_spatial_tables(paths, destination):
    """Stream large element tables instead of collecting all records in RAM."""
    writer = None
    with open_text(destination, "wt") as target:
        for path in paths:
            with open_text(path) as source:
                reader = csv.DictReader(source, delimiter="\t")
                if reader.fieldnames == ["status"]:
                    continue
                if writer is None:
                    writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t")
                    writer.writeheader()
                elif reader.fieldnames != writer.fieldnames:
                    raise ValueError(f"Inconsistent spatial table fields: {path}")
                writer.writerows(reader)
        if writer is None:
            target.write("status\n")


@transactional("spatial_manifest.json")
def modern_run_spatial(args):
    if args.end_bp < 1 or args.body_bins < 2 or args.block_bp < 100000:
        raise ValueError("Require end-bp > 0, body-bins >= 2 and physical blocks >= 100 kb")
    families = [family.strip() for family in args.families.split(",")]
    if not all(families) or len(set(families)) != len(families):
        raise ValueError("Families must be distinct nonempty exact RepeatMasker family strings")
    args.families = ",".join(families)
    run, annotation, tasks, source_chromosomes = modern_targeted_inputs(args)
    options, out = run["options"], Path(args.out)
    if not {5, 10} <= {int(t) for t in options.get("thresholds", [])}:
        raise ValueError("The source architecture scan must contain AC5 and AC10 counters for spatial reconciliation")
    if not options.get("ld_block_dir"):
        raise ValueError("Primary spatial decomposition requires the EUR LD partition recorded in the source scan")
    files, used = [], {c for c, _ in tasks}
    for chrom, source in tasks:
        print(f"[spatial] indexed chr{chrom} aligned-body query", flush=True)
        files.append(modern_scan_spatial_chrom(source, chrom, annotation/f"chr{chrom}.tsv.gz", out/"cache", options, args))
    with_eur = bool(options.get("ld_block_dir"))
    global_paths = modern_architecture_paths(args.modern_dir, "global")
    background = modern_Background(global_paths, args.min_training_sites)
    primary_background = modern_PrimaryBackground(global_paths, families)
    data, raw = modern_spatial_counts(files, background, with_eur)
    partitions, profiles, identities, contrasts = modern_summarize_spatial(data, args)
    primary_data, _ = modern_spatial_counts(files, primary_background, with_eur, collect_raw=False)
    comparison_scope = "all_22_autosomes" if used == MODERN_AUTOSOMES else "selected_chromosomes_subset"
    reconciliation, mismatches, labels = modern_reconcile_spatial_counters(args.modern_dir, raw, background, used,
                                                                families, with_eur, comparison_scope, primary_background)
    # All recorded global chromosomes train the primary rates; target totals
    # are restricted to the exact chromosomes actually queried by this stage.
    counter_totals = primary_scoped_counter_totals(modern_architecture_paths(args.modern_dir, "exposure"),
                                          primary_background, families, used)
    primary_parts, primary_profiles, primary_identities = modern_summarize_primary_spatial(primary_data, args, counter_totals)
    primary_reconciliation, primary_failures = primary_reconcile_spatial_to_counters(counter_totals, primary_parts)
    primary_checks = {(r["family"], r["sensitivity"], r["threshold"]): r for r in primary_reconciliation}
    for record in [*primary_parts, *primary_profiles, *primary_identities]:
        check = primary_checks.get((record["family"], record["sensitivity"], record["threshold"]), {})
        record.update(comparison_scope=comparison_scope, counter_reconciliation_status=check.get("status"))
    modern_write_rows(out/"spatial_partitions.tsv", partitions)
    modern_write_rows(out/"whole_element_profiles.tsv", profiles)
    modern_write_rows(out/"spatial_identity.tsv", identities)
    modern_write_rows(out/"spatial_contrasts.tsv", contrasts)
    modern_write_rows(out/"primary_spatial_partitions.tsv", primary_parts)
    modern_write_rows(out/"primary_whole_element_profiles.tsv", primary_profiles)
    modern_write_rows(out/"primary_spatial_identity.tsv", primary_identities)
    modern_write_rows(out/"primary_counter_reconciliation.tsv", primary_reconciliation)
    modern_write_rows(out/"family_counter_reconciliation.tsv", reconciliation)
    modern_write_rows(out/"family_counter_mismatches.tsv.gz", mismatches,
               fields=["chrom", "family", "sensitivity", "threshold", "scope", "context", "af_bin", "an_bin",
                       "counter_sites", "target_sites", "counter_multi", "target_multi"])
    modern_write_rows(out/"family_counter_labels.tsv", labels)
    modern_merge_spatial_tables([out/"cache"/f"chr{c}.annotation_audit.tsv" for c, _ in tasks], out/"spatial_annotation_audit.tsv")
    modern_merge_spatial_tables([out/"cache"/f"chr{c}.elements.tsv.gz" for c, _ in tasks], out/"spatial_elements.tsv.gz")
    if not mismatches and not primary_failures and not getattr(args, "no_plots", False):
        label = "Primary full-scope estimator · LOCO min 1"
        if used != source_chromosomes:
            label += "\nSelected target chromosomes"
        modern_plot_spatial(primary_parts, primary_profiles, out, args.end_bp, prefix="primary_", model_label=label)
        modern_plot_spatial(partitions, profiles, out, args.end_bp,
                     model_label=f"Common-support sensitivity · autosomes · LOCO min {args.min_training_sites}")
    manifest = dict(schema=MODERN_SPATIAL_SCHEMA, status="counter_count_mismatch" if mismatches else "primary_counter_mismatch" if primary_failures else "complete",
                    modern_dir=str(Path(args.modern_dir).resolve()),
                    families=families, exact_family_label_policy="LTR class and exact raw family strings; no alias collapsing",
                    chromosomes=sorted(used, key=int), available_autosomal_chromosomes=sorted(source_chromosomes, key=int),
                    comparison_scope=comparison_scope, full_available_autosomal_scope=used == source_chromosomes,
                    raw_counts_matched=not mismatches, mismatched_strata=len(mismatches),
                    primary_model="primary_full_scope", primary_counter_matches=not primary_failures,
                    primary_counter_mismatches=len(primary_failures), primary_minimum_training_sites=1,
                    primary_background_chromosomes=sorted({c for c, _ in primary_background.counts}, key=scan.chrom_sort_key),
                    primary_missing_rule="Retain raw observed and eligible counts; any unsupported target stratum makes the full-scope total and its weighted decomposition unavailable",
                    models_share_one_count_cache=True, additional_vcf_queries_for_primary_model=0,
                    primary_scope="eur_ld_partition" if with_eur else "all_family_sites",
                    end_bp=args.end_bp, body_bins=args.body_bins, physical_cluster_bp=args.block_bp,
                    denominator="Base-eligible SNV sites, fixed across AC5/AC10",
                    background="LOCO all recorded global chromosomes, minimum 1; canonical trinucleotide / primary ALT AF / AN; (multi+.5)/(sites+1)",
                    common_support_denominator="Unsupported observed and expected counts removed together; autosomal training background",
                    minimum_training_sites=args.min_training_sites, minimum_training_sites_applies_to="common_support_sensitivity",
                    minimum_support=args.min_support,
                    partition_definition={"head":f"first {args.end_bp} aligned bases in strand order",
                                          "tail":f"last {args.end_bp} aligned bases in strand order",
                                          "body":"remaining aligned bases of a long joined record",
                                          "short_element":f"entire record when aligned length <= {2*args.end_bp} bp",
                                          "ambiguous_overlap":"conflicting strand or part assignments; retained in total"},
                    normalized_position="base centres along concatenated merged aligned segments, reversed on minus strand; conflicting bins retained as -1",
                    accounting="overall O/E=sum(E_part/E_total*O_part/E_part); sum excess contributions=overall O/E-1; residual classes retained",
                    counter_reconciliation="N/M checked per chromosome/context/AF/AN/threshold/sensitivity; E recomputed from the same counts and training rule",
                    inference=f"{args.block_bp:g}bp cluster log-ratio sandwich pointwise CI; head/body and tail/body ratios; BH across strict AC5 primary-scope contrasts",
                    caveats=["Joined repeat records are not necessarily complete proviruses, true LTRs, insertion boundaries or TSSs",
                             "The -int suffix is a name audit, not a reconstructed structural assignment",
                             "Sites are counted once within each family; overlapping different families cannot be added",
                             "Normalized bins mix record lengths and are exploratory; short fragments remain separate in the fixed-bp partition",
                             "This is an SNV-denominator statistic, not a mutation rate per callable base",
                             "Background offsets are treated as fixed; local GC/recombination/repeat-age effects are not fitted",
                             "Aligned-family queries exclude external flanks; use the windows stage for flanks",
                             "Spatial physical-block CIs and primary LD-block bootstrap CIs use different uncertainty estimators"],
                    cache_files=[str(p.resolve()) for p in files])
    dump_json(out/"spatial_manifest.json", manifest)
    if mismatches:
        raise ValueError(f"Targeted spatial counts disagree with {len(mismatches)} source counter strata; see {out/'family_counter_mismatches.tsv.gz'}")
    if primary_failures:
        raise ValueError(f"Spatial enrichment disagrees with {len(primary_failures)} same-scope counter totals; see {out/'primary_counter_reconciliation.tsv'}")
    return manifest


def modern_parser():
    p = argparse.ArgumentParser(description=_description_herv_modern)
    p.add_argument("command", choices=["summary", "spatial", "windows"])
    p.add_argument("--modern-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--erv-joined")
    p.add_argument("--fasta")
    p.add_argument("--families", default="ERV1,ERVK")
    p.add_argument("--ibdmix", type=boolean_option, default=False, metavar="TRUE|FALSE")
    p.add_argument("--archaic-references", default=",".join(scan.DEFAULT_ARCHAIC_REFERENCES),
                   help="Common-support summary reference selection; .2013 references require explicit selection; reads existing counters only")
    p.add_argument("--end-bp", type=int, default=200, help="Spatial head/tail length in aligned bases; records <= twice this length form the short-element part")
    p.add_argument("--body-bins", type=int, default=20, help="Exploratory strand-oriented normalized aligned-body bins")
    p.add_argument("--flank", type=int, default=2000)
    p.add_argument("--window", type=int, default=200)
    p.add_argument("--step", type=int, default=50)
    p.add_argument("--block-bp", type=int, default=5000000)
    p.add_argument("--min-blocks", type=int, default=20)
    p.add_argument("--min-training-sites", type=int, default=100)
    p.add_argument("--min-support", type=float, default=.99)
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20261027)
    p.add_argument("--discovery-chromosomes", default="1,3,5,7,9,11,13,15,17,19,21")
    p.add_argument("--chromosomes")
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--no-plots", action="store_true", help="Skip spatial plots")
    for name in ("segdup", "lcr", "low-mappability"):
        p.add_argument("--"+name+"-bed")
    return p


def modern_main(argv=None):
    a = modern_parser().parse_args(argv)
    if not 0 < a.min_support <= 1 or a.min_training_sites < 1:
        raise ValueError("Require 0 < min-support <=1 and positive training count")
    if a.command == "summary":
        return modern_summarize_common_support(a.modern_dir, a.out, a.min_training_sites, a.min_support, a.bootstrap, a.seed,
                                        a.archaic_references, include_ibdmix=a.ibdmix)
    if not a.erv_joined:
        raise ValueError("--erv-joined required for targeted spatial/boundary scan")
    if a.command == "spatial":
        return modern_run_spatial(a)
    return modern_run_windows(a)


def _entry_herv_modern():
    try:
        modern_main()
    except (ValueError, OSError, RuntimeError) as modern_exc:
        raise SystemExit(f"herv-modern: {modern_exc}")



# ========================================================================
# herv candidates
# ========================================================================
_description_herv_candidates = 'HERV allele catalogs and seeded validation panels from completed site caches.\n\nThis is a conditional candidate design: both cases and controls already have\nat least two observed ALT nucleotides with AC>=5. It does not estimate HERV\nenrichment or test differences in the allele frequencies used for matching.'
#!/usr/bin/env python3
"""HERV allele catalogs and seeded validation panels from completed site caches.

This is a conditional candidate design: both cases and controls already have
at least two observed ALT nucleotides with AC>=5. It does not estimate HERV
enrichment or test differences in the allele frequencies used for matching.
"""
import random
import re
from allele_common import chrom_key
from scan_store import MASK_OPTIONS, resolve_mask_options


CANDIDATES_FAMILIES = ("ERV1", "ERVK")
CANDIDATES_AUTOSOMES = {str(i) for i in range(1, 23)}
CANDIDATES_FIELDS = ["chrom", "pos", "ref", "alts", "site_id", "family", "group", "match_id",
          "overlap_family", "case_region", "candidate_families", "context", "canonical_context",
          "partition", "comparison_partition", "body_bin", "spatial_overlap", "family_partitions",
          "primary_alt", "primary_af", "primary_an", "secondary_alt", "secondary_af",
          "primary_af_bin", "an_bin", "secondary_af_bin", "n_alt_ac5", "n_alt_ac10",
          "ac10_supported", "observed_alts_all", "alt_ac_all", "alt_an_all",
          "ref_presence_status", "strict_filter_complete", "missing_strict_filters",
          "match_stratum", "sampling_probability", "sampling_weight"]


def candidates_metadata(path):
    p = Path(path).expanduser().resolve()
    s = p.stat()
    return {"path": str(p), "size": s.st_size, "mtime_ns": s.st_mtime_ns}


def candidates_write_json(path, value):
    p = Path(path); tmp = p.with_name(p.name + ".partial")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(tmp, p)


def candidates_write_tsv(path, fields, rows):
    p = Path(path); tmp = p.with_name(p.name + ".partial")
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, p)


class candidates_Reservoir:
    """Uniform bounded sample; sample size is independent of stream length."""
    def __init__(self, capacity, seed):
        self.capacity = capacity; self.rng = random.Random(seed)
        self.n = 0; self.rows = []

    def add(self, row):
        self.n += 1
        if self.capacity is None or len(self.rows) < self.capacity:
            self.rows.append(row)
        else:
            j = self.rng.randrange(self.n)
            if j < self.capacity:
                self.rows[j] = row


def candidates_stratified_quotas(counts, maximum):
    """Proportional largest-remainder quotas; one per chromosome when feasible."""
    counts = {k: n for k, n in counts.items() if n}
    budget = min(maximum, sum(counts.values()))
    if not budget:
        return {k: 0 for k in counts}
    keys = sorted(counts, key=chrom_key)
    base = int(budget >= len(keys))
    quotas = {k: base for k in keys}
    remaining = budget - sum(quotas.values())
    available = {k: counts[k] - quotas[k] for k in keys}
    total = sum(available.values())
    raw = {k: remaining * available[k] / total if total else 0 for k in keys}
    for k in keys:
        quotas[k] += math.floor(raw[k])
    left = budget - sum(quotas.values())
    order = sorted(keys, key=lambda k: (-(raw[k] - math.floor(raw[k])), chrom_key(k)))
    for k in order[:left]:
        quotas[k] += 1
    return quotas


def candidates_load_sources(modern_dir):
    root = Path(modern_dir).expanduser().resolve()
    run_path = root / "run.json"
    run = scan.load_run_metadata(root)
    if run.get("status") not in {"scanned", "complete"}:
        raise ValueError("Use a completed/scanned modern run.json")
    options = run.get("options", run.get("identity", {}).get("options", {}))
    if options.get("build", "GRCh38") not in {"38", "GRCh38", "hg38"}:
        raise ValueError("Candidates currently require GRCh38 modern shards and annotation")
    shards = []
    for raw in run.get("shards", []):
        d = Path(raw).expanduser()
        if not d.is_absolute(): d = root / d
        d = d.resolve()
        receipt = json.loads((d / "COMPLETE.json").read_text())
        if not run.get("signature") or receipt.get("signature") != run["signature"]:
            raise ValueError(f"Shard completion signature differs from run.json: {d}")
        site = d / "multiallelic_sites.tsv.gz"
        saved = receipt.get("outputs", {}).get(site.name)
        if not site.is_file() or (saved and saved.get("size") != site.stat().st_size):
            raise ValueError(f"Missing or size-changed site cache: {site}")
        info = json.loads((d / "scan.json").read_text())
        if not info.get("complete"):
            raise ValueError(f"Incomplete scan.json: {d}")
        scope = info.get("scope", [])
        if not scope:
            raise ValueError(f"Missing chromosome scope in {d}/scan.json")
        shards.append({"directory": d, "site": site, "flags": set(info.get("fields_read", [])) & {"lcr", "segdup"},
                       "scope": scope, "metadata": candidates_metadata(site), "scan_metadata": candidates_metadata(d / "scan.json")})
    if not shards or len({s["directory"] for s in shards}) != len(shards):
        raise ValueError("No shards or duplicate shard directories in run.json")
    shards.sort(key=lambda s: (min(chrom_key(c) for c in s["scope"]), str(s["directory"])))
    return run, options, shards


def candidates_split_annotation(source, options, cache):
    """Reuse verified-by-metadata chromosome cache read-only, otherwise split once."""
    source_meta = candidates_metadata(source)
    def valid_files(directory, manifest):
        for i in range(1, 23):
            name = f"chr{i}.tsv.gz"; p = directory / name
            saved = manifest.get("files", {}).get(name, {})
            if not p.is_file() or not saved: return False
            actual = candidates_metadata(p)
            if any(saved.get(k) != actual[k] for k in ("size", "mtime_ns")): return False
        return True
    previous = options.get("architecture_annotation_cache")
    if previous:
        marker = Path(previous) / "manifest.json"
        if marker.is_file():
            manifest = json.loads(marker.read_text())
            recorded = manifest.get("source", {})
            same = all(recorded.get(k) == v for k, v in source_meta.items())
            if same and valid_files(Path(previous), manifest):
                return Path(previous), "original_read_only"
    cache.mkdir(parents=True, exist_ok=True)
    marker = cache / "manifest.json"
    if marker.is_file():
        old = json.loads(marker.read_text())
        if old.get("source") == source_meta and valid_files(cache, old):
            return cache, "derived_chromosome_cache"
    handles = {}
    try:
        for i in range(1, 23): handles[str(i)] = gzip.open(cache / f"chr{i}.tsv.gz.partial", "wt", compresslevel=2)
        with open_text(source) as f:
            for ln, line in enumerate(f, 1):
                if not line.strip() or line.startswith("#"): continue
                # The final joined-RMSK field may be empty; preserve trailing tabs.
                a = line.rstrip("\r\n").split("\t") if "\t" in line else line.split()
                if len(a) != 14:
                    raise ValueError(f"{source}:{ln}: expected 14 joined RepeatMasker columns, got {len(a)}")
                chrom = canonical_chrom(a[1])
                if chrom in handles: handles[chrom].write(line if line.endswith("\n") else line + "\n")
    finally:
        for h in handles.values(): h.close()
    for i in range(1, 23): os.replace(cache / f"chr{i}.tsv.gz.partial", cache / f"chr{i}.tsv.gz")
    candidates_write_json(marker, {"source": source_meta, "method": "autosomal split; metadata consistency, no full-file checksums",
                        "files": {f"chr{i}.tsv.gz": candidates_metadata(cache / f"chr{i}.tsv.gz") for i in range(1, 23)}})
    return cache, "derived_chromosome_cache"


def candidates_load_masks(options):
    # Compact merged interval arrays are reused between both passes.
    return {name: IntervalMask(options[key], [str(i) for i in range(1, 23)])
            for name, key in MASK_OPTIONS.items() if options.get(key)}


def candidates_parse_site(row):
    """Preserve observed ALT identities; never infer an observed REF allele."""
    chrom = canonical_chrom(row["chrom"]); pos = int(row["pos"])
    ref = row["ref"].upper()
    alts = row["alts_by_af"].split(",")
    ac = [int(x) for x in row["alt_ac"].split(",")]
    an = [int(x) for x in row["alt_an"].split(",")]
    if len({len(alts), len(ac), len(an)}) != 1 or len(set(alts)) != len(alts):
        raise ValueError(f"Misaligned/duplicated ALT arrays at {chrom}:{pos}")
    if pos < 1 or ref not in "ACGT" or len(ref) != 1 or any(a not in {"A", "C", "G", "T"} or a == ref for a in alts):
        raise ValueError(f"Invalid observed SNV identity at {chrom}:{pos}")
    if any(n <= 0 or c < 1 or c > n for c, n in zip(ac, an)):
        raise ValueError(f"Invalid ALT AC/AN at {chrom}:{pos}")
    n5 = sum(c >= 5 for c in ac); n10 = sum(c >= 10 for c in ac)
    if "n_alt_ac5" in row and int(row["n_alt_ac5"]) != n5:
        raise ValueError(f"Cached AC5 count disagrees with observed ALTs at {chrom}:{pos}")
    if "n_alt_ac10" in row and int(row["n_alt_ac10"]) != n10:
        raise ValueError(f"Cached AC10 count disagrees with observed ALTs at {chrom}:{pos}")
    if n5 < 2: return None
    kept = sorted([(a, c, n, c / n) for a, c, n in zip(alts, ac, an) if c >= 5], key=lambda x: (-x[3], x[0]))
    context = row.get("context", "NA")
    canonical = row.get("canonical_context") or scan.canonical_context(context, ref)[0]
    first, second = kept[:2]
    return {"chrom": chrom, "pos": pos, "ref": ref, "alts": ",".join(x[0] for x in kept),
            "site_id": f"{chrom}:{pos}:{ref}", "context": context, "canonical_context": canonical,
            "primary_alt": first[0], "primary_af": first[3], "primary_an": first[2],
            "secondary_alt": second[0], "secondary_af": second[3],
            "primary_af_bin": af_bin(first[3]), "an_bin": an_bin(first[2]), "secondary_af_bin": af_bin(second[3]),
            "n_alt_ac5": n5, "n_alt_ac10": n10, "ac10_supported": int(n10 >= 2),
            "observed_alts_all": ",".join(alts), "alt_ac_all": ",".join(map(str, ac)),
            "alt_an_all": ",".join(map(str, an)), "ref_presence_status": "not_inferable_from_sites_summary"}


class candidates_ChromAnnotation:
    def __init__(self, cache, chrom, *, end_bp, body_bins):
        self.elements = scan.parse_joined_rmsk(cache / f"chr{chrom}.tsv.gz", chromosome=chrom)
        self.body = scan.ERVOverlapSweep(self.elements)
        self.spatial = modern_AlignedBodies(self.elements, CANDIDATES_FAMILIES, end_bp=end_bp, bins=body_bins)
        self.last_spatial = {}
        self.boundaries = {fam: scan.BoundarySweep([e for e in self.elements if e.rep_class == "LTR" and e.family == fam])
                           for fam in CANDIDATES_FAMILIES}

    def classify(self, chrom, pos):
        overlap = self.body.locate(chrom, pos)
        self.last_spatial = self.spatial.assign(overlap, pos)
        body_families = {e.family for e in overlap if e.rep_class == "LTR" and e.family in CANDIDATES_FAMILIES}
        near200, near2000 = set(), set()
        for fam, sweep in self.boundaries.items():
            for edge in ("5p", "3p"):
                near = sweep.nearest(chrom, pos, edge, max_distance=2000)
                if near:
                    near2000.add(fam)
                    if near["distance"] < 200: near200.add(fam)
        return body_families, near200, bool(near2000), any(e.rep_class == "LTR" for e in overlap), any(
            e.rep_class in {"Simple_repeat", "Low_complexity"} for e in overlap)


def candidates_eligible_rows(shards, annotation_cache, masks, counts, *, end_bp, body_bins):
    last = None; current = None; annotation = None
    for shard in shards:
        print(f"[herv-candidates] read {shard['site']}", flush=True)
        with gzip.open(shard["site"], "rt") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for raw in reader:
                chrom = canonical_chrom(raw["chrom"]); key = (chrom_key(chrom), int(raw["pos"]))
                if last is not None and key <= last:
                    raise ValueError("Site caches must be globally chromosome/position sorted and have no duplicate positions")
                last = key; counts["rows_read"] += 1
                if chrom not in CANDIDATES_AUTOSOMES: continue
                if "n_alt_ac5" in raw and int(raw["n_alt_ac5"]) < 2: continue
                row = candidates_parse_site(raw)
                if row is None: continue
                counts["ac5_eligible"] += 1
                flags = set(re.split("[,;|]", raw.get("region_flags", ".")))
                missing = sorted(set(MASK_OPTIONS) - set(masks) - shard["flags"])
                # Match the original ArchitectureAccumulator: an explicit BED
                # defines that mask; INFO is a fallback only when no BED exists.
                excluded = any(mask.contains(chrom, row["pos"]) for mask in masks.values())
                excluded |= any(name in flags for name in ("lcr", "segdup") if name not in masks)
                if excluded:
                    counts["excluded_difficult_regions"] += 1; continue
                if current != chrom:
                    current = chrom; annotation = candidates_ChromAnnotation(annotation_cache, chrom, end_bp=end_bp, body_bins=body_bins)
                body, near, near_control, any_ltr, simple = annotation.classify(chrom, row["pos"])
                if simple:
                    counts["excluded_simple_low_complexity"] += 1; continue
                row.update(strict_filter_complete=int(not missing), missing_strict_filters=",".join(missing) or ".")
                row["_spatial_parts"] = annotation.last_spatial
                row["family_partitions"] = ";".join(f"{family}:{value[0]}" for family,value in sorted(annotation.last_spatial.items())) or "."
                row["match_stratum"] = "|".join(str(row[k]) for k in ("chrom", "canonical_context", "primary_af_bin", "an_bin", "secondary_af_bin"))
                yield row, body, near, near_control, any_ltr


def candidates_match_key(row):
    # A per-shard missing annotation must not silently weaken just one match arm.
    return row["match_stratum"], row["missing_strict_filters"]


@transactional("candidate_manifest.json")
def candidates_run_candidates(args):
    if args.max_sites < 0 or args.controls_per_case < 0:
        raise ValueError("--max-sites must be nonnegative (0 means all); --controls-per-case must be nonnegative")
    end_bp, body_bins = args.end_bp, args.body_bins
    if end_bp < 1 or body_bins < 2:
        raise ValueError("--end-bp must be positive and --body-bins at least two")
    spatial_definition = {"schema": "aligned-record-partition-1", "end_bp": end_bp, "body_bins": body_bins}
    catalog_all = bool(getattr(args, "all_alleles_catalog", False))
    root = Path(args.modern_dir).expanduser().resolve(); out = Path(args.out).expanduser().resolve()
    run, options, shards = candidates_load_sources(root)
    options = resolve_mask_options(run, {name: getattr(args, name + "_bed", None)
                                         for name in MASK_OPTIONS}, verify=True)
    if out == root or any(out == s["directory"] or out.is_relative_to(s["directory"]) for s in shards):
        raise ValueError("Candidate output must not overwrite the original modern run/shards")
    identity = {"schema": 3, "spatial_definition": spatial_definition,
                "code": {name: candidates_metadata(Path(__file__).with_name(name))
                         for name in ('herv.py', 'herv.py', "scan.py", "scan_store.py", "scan_store.py", "0.common.py")},
                "modern_run": candidates_metadata(root / "run.json"), "source_signature": run.get("signature"),
                "sites": [s["metadata"] for s in shards], "scan_metadata": [s["scan_metadata"] for s in shards],
                "erv_joined": candidates_metadata(args.erv_joined),
                "masks": {k: candidates_metadata(options[v]) for k, v in MASK_OPTIONS.items() if options.get(v)},
                "max_sites_per_family": args.max_sites, "all_alleles_catalog": catalog_all,
                "controls_per_case": args.controls_per_case, "seed": args.seed,
                "method": "body_or_abs_boundary_lt200; controls_noLTR_and_abs_boundary_ge2000; AC5-supported_ALT_AF_matching_v1"}
    receipt = out / "candidate_manifest.json"
    if receipt.is_file():
        old = json.loads(receipt.read_text())
        if old.get("status") == "complete" and old.get("identity") == identity:
            if all((out / p).is_file() and candidates_metadata(out / p) == fp for p, fp in old.get("outputs", {}).items()) and old.get("outputs"):
                print(f"[herv-candidates] reuse {out}", flush=True)
                return old
    out.mkdir(parents=True, exist_ok=True)
    annotation_cache, cache_mode = candidates_split_annotation(args.erv_joined, options, out / ".cache" / "joined")
    masks = candidates_load_masks(options)
    pools = {}; pass1 = Counter()
    for row, body, near, _, _ in candidates_eligible_rows(shards, annotation_cache, masks, pass1, end_bp=end_bp, body_bins=body_bins):
        families = body | near
        if not families: continue
        # A multiply annotated locus belongs to one seeded family for selection;
        # original memberships are retained, preventing duplicated target loci.
        family = random.Random(f"{args.seed}:family:{row['site_id']}").choice(sorted(families))
        part, body_bin, spatial_overlap = row["_spatial_parts"].get(family, ("boundary_only", -1, "outside_aligned_body"))
        row.update(family=family, group="case", overlap_family=",".join(sorted(body)) or "boundary_only",
                   candidate_families=",".join(sorted(families)), case_region=("body_and_boundary" if family in body & near else "body" if family in body else "boundary"),
                   partition=part, comparison_partition=part, body_bin=body_bin, spatial_overlap=spatial_overlap)
        key = family, row["chrom"]
        if key not in pools: pools[key] = candidates_Reservoir(args.max_sites or None, f"{args.seed}:case:{family}:{row['chrom']}")
        pools[key].add(row)
    cases = []; sampling = []
    for family in CANDIDATES_FAMILIES:
        counts = {c: p.n for (f, c), p in pools.items() if f == family}
        quotas = candidates_stratified_quotas(counts, args.max_sites or sum(counts.values()))
        for chrom in sorted(counts, key=chrom_key):
            pool = pools[family, chrom]; quota = quotas[chrom]
            random.Random(f"{args.seed}:quota:{family}:{chrom}").shuffle(pool.rows)
            probability = quota / pool.n
            for row in pool.rows[:quota]:
                row.update(sampling_probability=probability, sampling_weight=1 / probability)
                cases.append(row)
            sampling.append({"family": family, "chrom": chrom, "eligible_cases": pool.n, "selected_cases": quota,
                             "sampling_probability": probability, "sampling_weight": 1 / probability if probability else "NA"})
    cases.sort(key=lambda r: (r["family"], chrom_key(r["chrom"]), r["pos"]))
    demands = defaultdict(list)
    for i, row in enumerate(cases, 1):
        row["match_id"] = f"HERV{i:06d}"
        demands[candidates_match_key(row)].append(row)
    control_pools = {k: candidates_Reservoir(len(v) * args.controls_per_case, f"{args.seed}:control:{k}")
                     for k, v in demands.items() if args.controls_per_case}
    pass2 = Counter()
    if control_pools:
        for row, _, _, near_control, any_ltr in candidates_eligible_rows(shards, annotation_cache, masks, pass2, end_bp=end_bp, body_bins=body_bins):
            if any_ltr or near_control: continue
            key = candidates_match_key(row)
            if key in control_pools: control_pools[key].add(row)
    targets = list(cases); audit = []; used_controls = set()
    for key, group in sorted(demands.items()):
        pool = control_pools.get(key)
        controls = list(pool.rows) if pool else []
        random.Random(f"{args.seed}:allocate:{key}").shuffle(controls)
        matched = Counter()
        for i, control in enumerate(controls):
            case = group[i % len(group)]
            if control["site_id"] in used_controls:
                raise ValueError("Duplicate control locus encountered")
            used_controls.add(control["site_id"]); matched[case["match_id"]] += 1
            control.update(family=case["family"], group="control", match_id=case["match_id"], overlap_family="nonHERV",
                           candidate_families=".", case_region="nonLTR_distant_control", sampling_probability="NA", sampling_weight="NA",
                           partition="nonLTR_control", comparison_partition=case["partition"], body_bin=-1, spatial_overlap="nonLTR_control")
            targets.append(control)
        for case in group:
            n = matched[case["match_id"]]
            audit.append({"match_id": case["match_id"], "family": case["family"], "site_id": case["site_id"],
                          "match_stratum": case["match_stratum"], "missing_strict_filters": case["missing_strict_filters"],
                          "requested_controls": args.controls_per_case, "matched_controls": n,
                          "eligible_controls_in_stratum": pool.n if pool else 0, "cases_in_stratum": len(group),
                          "status": "complete" if n == args.controls_per_case else "insufficient_exact_controls"})
    if (any(candidates_metadata(s["site"]) != s["metadata"] or candidates_metadata(s["directory"] / "scan.json") != s["scan_metadata"] for s in shards)
            or candidates_metadata(root / "run.json") != identity["modern_run"]
            or candidates_metadata(args.erv_joined) != identity["erv_joined"]
            or any(candidates_metadata(options[MASK_OPTIONS[k]]) != saved for k, saved in identity["masks"].items())):
        raise ValueError("Inputs changed during candidate selection; rerun with stable inputs")
    targets.sort(key=lambda r: (r["match_id"], r["group"] != "case", chrom_key(r["chrom"]), r["pos"]))
    candidates_write_tsv(out / "targets.tsv", CANDIDATES_FIELDS, targets)
    candidates_write_tsv(out / "match_audit.tsv", ["match_id", "family", "site_id", "match_stratum", "missing_strict_filters",
              "requested_controls", "matched_controls", "eligible_controls_in_stratum", "cases_in_stratum", "status"], audit)
    candidates_write_tsv(out / "sampling.tsv", ["family", "chrom", "eligible_cases", "selected_cases", "sampling_probability", "sampling_weight"], sampling)
    output_names = ["targets.tsv", "match_audit.tsv", "sampling.tsv"]
    if catalog_all:
        candidates_write_allele_catalogs(out, targets)
        output_names += ["site_catalog.tsv.gz", "allele_catalog.tsv.gz"]
    result = {"status": "complete", "identity": identity, "spatial_definition": spatial_definition,
              "case_boundary_bp": 200, "control_boundary_exclusion_bp": 2000,
              "annotation_cache_mode": cache_mode,
              "cases_by_family": dict(Counter(r["family"] for r in cases)), "n_cases": len(cases), "n_controls": len(used_controls),
              "n_unmatched_cases": sum(r["matched_controls"] == 0 for r in audit), "n_incomplete_matches": sum(r["status"] != "complete" for r in audit),
              "pass1": dict(pass1), "pass2": dict(pass2),
              "all_selected_strict_complete": bool(targets) and all(r["strict_filter_complete"] for r in targets),
              "coordinate_build": "GRCh38", "site_coordinates": "1-based", "boundaries": "observed joined repeat edges, not TSS; abs(distance)<200 bp",
              "case_selection": ("ALL AC5 count-supported multi-ALT HERV cases; no site cap; no ancient statistic or window peak used" if args.max_sites == 0 else
                                 "AC5 count-supported multi-ALT; per-family chromosome-stratified reservoir; no ancient statistic or window peak used"),
              "all_observed_alt_catalog": catalog_all,
              "catalog_scope": "All selected modern loci with all observed ALT bases (including additional ALT AC<5, labelled separately). REF is a label, not an inferred observed count. All cached eligible HERV loci only when max_sites=0; not all single-ALT genomic positions.",
              "overlapping_family_assignment": "one seeded family per locus; memberships retained in candidate_families",
              "controls": "unique loci globally; both arms AC5 multi-ALT; no aligned LTR body and >=2000 bp from ERV1/ERVK observed edges",
              "matching": "chrom, canonical trinucleotide, primary ALT AF bin, primary ALT AN bin, secondary ALT AF bin; same missing-filter signature",
              "af_definition": "highest and second-highest AC/AN among AC>=5 observed ALT alleles; no inferred REF AF",
              "estimand_warning": "Candidate follow-up conditional on contemporary multi-ALT occupancy and AF/AN matching. Not an HERV enrichment estimate, an allele-age/retention test, or an AF difference test.",
              "input_validation": "completion signatures, file sizes and modification metadata; explicit mask relocations checked against original size/edge fingerprints; no added full-file checksums",
              "outputs": {p: candidates_metadata(out / p) for p in output_names}}
    candidates_write_json(receipt, result)
    print(f"[herv-candidates] {len(cases)} cases, {len(used_controls)} unique controls; {out / 'targets.tsv'}", flush=True)
    return result


def candidates_write_allele_catalogs(out, targets):
    """Preserve every observed modern base and count for the complete time panel.

    A locus is selected because at least two ALTs meet AC5. Extra lower-count
    bases remain explicitly lower-count, and reference counts are never inferred
    by subtracting unrelated per-ALT denominators.
    """
    site_fields = [*CANDIDATES_FIELDS, "all_alt_ac", "all_alt_an", "catalog_scope"]
    allele_fields = ["allele_id", "site_id", "chrom", "pos", "ref", "allele", "role",
                     "modern_ac", "modern_an", "modern_af", "ac5_supported", "ac10_supported",
                     "family", "group", "match_id", "partition", "comparison_partition",
                     "candidate_families", "family_partitions", "strict_filter_complete", "missing_strict_filters"]
    site_path, allele_path = Path(out)/"site_catalog.tsv.gz", Path(out)/"allele_catalog.tsv.gz"
    temporary = [p.with_name(p.name+".partial") for p in (site_path, allele_path)]
    with gzip.open(temporary[0], "wt", newline="", compresslevel=3) as sf, gzip.open(temporary[1], "wt", newline="", compresslevel=3) as af:
        sw = csv.DictWriter(sf, site_fields, delimiter="\t", extrasaction="ignore")
        aw = csv.DictWriter(af, allele_fields, delimiter="\t", extrasaction="ignore")
        sw.writeheader(); aw.writeheader()
        for row in targets:
            alts = row["observed_alts_all"].split(",")
            ac, an = row["alt_ac_all"].split(","), row["alt_an_all"].split(",")
            if not len(alts) == len(ac) == len(an):
                raise ValueError(f"Catalog allele/count lengths disagree: {row['site_id']}")
            sw.writerow(dict(row, alts=",".join(alts), all_alt_ac=",".join(ac), all_alt_an=",".join(an),
                             catalog_scope="modern_AC5_multiALT_loci_all_observed_bases"))
            common = {k:row.get(k, ".") for k in allele_fields}
            aw.writerow(common | dict(allele_id=row["site_id"]+":"+row["ref"], allele=row["ref"], role="reference_label",
                        modern_ac="NA", modern_an="NA", modern_af="NA", ac5_supported="NA", ac10_supported="NA"))
            for base, count, denominator in zip(alts, ac, an):
                count, denominator = int(count), int(denominator)
                role = "primary" if base == row["primary_alt"] else "secondary" if base == row["secondary_alt"] else "additional_supported" if count >= 5 else "additional_low_count"
                aw.writerow(common | dict(allele_id=row["site_id"]+":"+base, allele=base, role=role,
                            modern_ac=count, modern_an=denominator, modern_af=count/denominator,
                            ac5_supported=int(count >= 5), ac10_supported=int(count >= 10)))
    for temp, destination in zip(temporary, (site_path, allele_path)):
        os.replace(temp, destination)


def candidates_main(argv=None):
    p = argparse.ArgumentParser(description=_description_herv_candidates)
    p.add_argument("--modern-dir", required=True, help="Completed modern run.json plus original shards")
    p.add_argument("--erv-joined", required=True, help="GRCh38 14-column joined RepeatMasker table")
    p.add_argument("--out", required=True, help="Separate derived candidate output directory")
    for name in MASK_OPTIONS:
        p.add_argument("--" + name.replace("_", "-") + "-bed",
                       help="Relocated copy of the original scan mask; verifies saved size/edge fingerprint")
    p.add_argument("--max-sites", type=int, default=500, help="Maximum CASE sites PER FAMILY; 0 retains ALL qualified HERV cases")
    p.add_argument("--all-alleles-catalog", action="store_true", help="Also write full-base site_catalog.tsv.gz and allele_catalog.tsv.gz for ancient trajectories")
    p.add_argument("--controls-per-case", type=int, default=4)
    p.add_argument("--end-bp", type=int, default=200)
    p.add_argument("--body-bins", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261027)
    args = p.parse_args(argv)
    try:
        candidates_run_candidates(args)
    except (ValueError, OSError, KeyError) as e:
        p.exit(2, f"herv-candidates: {e}\n")
    return 0


def _entry_herv_candidates():
    raise SystemExit(candidates_main())


# ========================================================================
# herv genotypes
# ========================================================================
_description_herv_genotypes = 'Validate predefined HERV targets with actual A/C/G/T genotype counts.\n\nIndexed 1KG concordance in a fixed ancestry-balanced subset, not independent\nreplication unless sample overlap has been ruled out. No genome-wide VCF scan.'
#!/usr/bin/env python3
"""Validate predefined HERV targets with actual A/C/G/T genotype counts.

Indexed 1KG concordance in a fixed ancestry-balanced subset, not independent
replication unless sample overlap has been ruled out. No genome-wide VCF scan.
"""
from allele_measurement import phase_relation


GENOTYPES_GROUPS = ("AFR", "EAS", "EUR", "SAS")
GENOTYPES_BASES = ("A", "C", "G", "T")
GENOTYPES_SITE_FIELDS = ["site_id", "chrom", "pos", "ref", "alts", "family", "group", "match_id", "population",
               "n_selected", "n_called", "AN", "call_rate", "A_AC", "C_AC", "G_AC", "T_AC",
               "encoded_alts", "expected_alts_not_encoded", "represented_status", "record_status",
               "representation", "n_phase_comparable", "n_phase_unknown", "n_phase_conflict", "n_implicit_phase_set", "n_explicit_phase_set", "n_missing_samples", "n_conflicting_samples", "n_bases_ac_ge_min",
               "observed_three_base_in_encoded_calls", "confirmed_multi_allelic", "target_alts_retained",
               "evaluable", "min_ac", "min_call_rate", "interpretation"]
GENOTYPES_SUMMARY_FIELDS = ["family", "group", "population", "n_targets", "n_all_expected_encoded", "n_evaluable",
                  "n_confirmed_multi_allelic", "confirmed_fraction", "confirmed_wilson_low", "confirmed_wilson_high",
                  "n_target_alts_retained", "retained_fraction", "retained_wilson_low", "retained_wilson_high",
                  "min_ac", "min_call_rate", "ci_method", "scope"]
GENOTYPES_BURDEN_FIELDS = ["sample_id", "population", "family", "group", "n_targets_requested", "n_targets_represented",
                 "called_sites", "called_alleles", "nonref_allele_copies", "nonref_per_called_site",
                 "nonref_allele_fraction", "sample_target_call_rate", "measure"]


def genotypes_text_open(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def genotypes_metadata(path):
    p = Path(path).expanduser().resolve(); s = p.stat()
    return {"path": str(p), "size": s.st_size, "mtime_ns": s.st_mtime_ns}


def genotypes_write_json(path, value):
    p = Path(path); tmp = p.with_name(p.name + ".partial")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(tmp, p)


def genotypes_write_tsv(path, fields, rows):
    p = Path(path); tmp = p.with_name(p.name + ".partial")
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, p)


def genotypes_read_panel(path):
    with genotypes_text_open(path) as f:
        rows = [line.split() for line in f if line.strip() and not line.startswith("##")]
    if not rows: raise ValueError("Empty sample panel")
    normalize = lambda s: re.sub("[^a-z0-9]", "", s.lower())
    header = [normalize(x) for x in rows[0]]
    sample_names = {"sample", "sampleid", "iid", "id"}
    pop_names = {"superpop", "superpopulation", "superpopulationcode"}
    if set(header) & sample_names:
        si = next(i for i, v in enumerate(header) if v in sample_names)
        pi = next((i for i, v in enumerate(header) if v in pop_names), None)
        if pi is None:
            pi = next((i for i, v in enumerate(header) if v == "population"), None)
        if pi is None: raise ValueError("Sample panel requires a superpopulation column")
        rows = rows[1:]
    else:
        si, pi = 0, (2 if len(rows[0]) > 2 else 1)
    panel = {}
    for row in rows:
        if len(row) <= max(si, pi): raise ValueError("Malformed sample-panel row")
        if row[si] in panel: raise ValueError(f"Duplicate sample ID in panel: {row[si]}")
        panel[row[si]] = row[pi].upper()
    return panel


def genotypes_read_targets(path):
    targets = []; seen = set()
    with genotypes_text_open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            chrom = r["chrom"].removeprefix("chr"); pos = int(r["pos"]); ref = r["ref"].upper()
            alts = tuple(r["alts"].upper().split(","))
            if chrom not in {str(i) for i in range(1, 23)} or pos < 1:
                raise ValueError("Targets must have autosomal GRCh38 1-based coordinates")
            if ref not in GENOTYPES_BASES or len(alts) < 2 or len(set(alts)) != len(alts) or any(a not in GENOTYPES_BASES or a == ref for a in alts):
                raise ValueError(f"Target must have at least two distinct observed SNV ALTs: {chrom}:{pos}")
            if (chrom, pos) in seen: raise ValueError(f"Duplicate target coordinate: {chrom}:{pos}")
            seen.add((chrom, pos))
            if r["group"] not in {"case", "control"}: raise ValueError("Target group must be case/control")
            targets.append(dict(r, chrom=chrom, pos=pos, ref=ref, alts=",".join(alts),
                                site_id=r.get("site_id") or f"{chrom}:{pos}:{ref}"))
    return sorted(targets, key=lambda r: (int(r["chrom"]), r["pos"]))


def genotypes_balanced_samples(common_samples, panel, seed):
    members = {g: sorted(s for s in common_samples if panel.get(s) == g) for g in GENOTYPES_GROUPS}
    n = min(map(len, members.values()))
    if not n: raise ValueError("At least one of AFR/EAS/EUR/SAS has no panel-labelled sample in the common VCF sample set")
    rng = random.Random(seed)
    selected = {g: sorted(rng.sample(members[g], n)) for g in GENOTYPES_GROUPS}
    return selected, {g: len(v) for g, v in members.items()}


def genotypes_decode_records(records, target_ref, samples):
    """Return complete diploid base-count dictionaries or sample-level missing.

    Split-biallelic REF dosage is 2-sum(ALT dosages) only if every record is
    completely called with the same REF and non-overlapping ALT identities.
    Negative REF dosage is a conflict, never clipped or counted as zero.
    Fully phased GTs without PS share the VCF implicit phase set and therefore
    undergo the same per-haplotype conflict check as equal explicit PS values.
    """
    result = {"record_status": "ok", "representation": "unavailable", "encoded_alts": set(),
              "calls": {}, "missing": set(), "conflicting": set(), "phase_comparable": set(), "phase_unknown": set(), "phase_conflict": set()}
    def fail(status):
        result["record_status"] = status
        return result
    if not records: return fail("record_not_found")
    if any(set(r.filter.keys()) - {"PASS", "."} for r in records): return fail("filtered_record_at_target")
    if any(r.ref != target_ref for r in records): return fail("inconsistent_or_mismatched_ref")
    alleles = [tuple(r.alleles or ()) for r in records]
    if any(len(a) < 2 or len(set(a)) != len(a) or any(x not in GENOTYPES_BASES for x in a) for a in alleles):
        return fail("non_snv_or_invalid_alleles")
    encoded = [a for rec in alleles for a in rec[1:]]
    result["encoded_alts"] = set(encoded)
    if len(set(encoded)) != len(encoded): return fail("duplicated_alt_records")
    if len(records) > 1 and any(len(a) != 2 for a in alleles): return fail("mixed_record_representation")
    result["representation"] = "single_multiallelic" if len(records) == 1 and len(alleles[0]) > 2 else "single_biallelic" if len(records) == 1 else "split_biallelic"
    for sample in samples:
        gts = [r.samples[sample].get("GT") for r in records]
        if any(gt is None or len(gt) != 2 or any(a is None for a in gt) for gt in gts):
            result["missing"].add(sample); continue
        if any(any(not isinstance(a, int) or a < 0 or a >= len(bases) for a in gt) for gt, bases in zip(gts, alleles)):
            result["conflicting"].add(sample); continue
        counts = dict.fromkeys(GENOTYPES_BASES, 0)
        if len(records) == 1:
            for a in gts[0]: counts[alleles[0][a]] += 1
        else:
            comparable, phase_origin = phase_relation(records, sample)
            if comparable:
                result.setdefault(phase_origin, set()).add(sample)
            result['phase_comparable' if comparable else 'phase_unknown'].add(sample)
            if comparable and any(sum(gt[h] == 1 for gt in gts) > 1 for h in (0, 1)):
                result['phase_conflict'].add(sample); result['conflicting'].add(sample); continue
            for gt, bases in zip(gts, alleles): counts[bases[1]] += sum(a == 1 for a in gt)
            total_alt = sum(counts.values())
            if total_alt > 2:
                result["conflicting"].add(sample); continue
            counts[target_ref] = 2 - total_alt
        result["calls"][sample] = counts
    return result


def genotypes_evaluate_population(target, decoded, population, samples, min_ac, min_call_rate):
    row = {k: target.get(k, ".") for k in ("site_id", "chrom", "pos", "ref", "alts", "family", "group", "match_id")}
    expected = set(target["alts"].split(",")); encoded = decoded["encoded_alts"]
    missing_alts = expected - encoded
    usable = decoded["record_status"] == "ok"
    represented = "all_expected_alts_encoded" if usable and not missing_alts else "expected_alt_not_encoded" if usable or decoded["record_status"] == "record_not_found" else "unavailable_record"
    calls = [decoded["calls"][s] for s in samples if s in decoded["calls"]]
    n = len(calls); call_rate = n / len(samples) if samples and usable else "NA"
    ac = {b: sum(c[b] for c in calls) for b in GENOTYPES_BASES}
    n_bases = sum(v >= min_ac for v in ac.values()) if n else "NA"
    evaluable = bool(usable and not missing_alts and n and call_rate >= min_call_rate)
    row.update(population=population, n_selected=len(samples), n_called=n if usable else "NA",
               AN=2*n if usable else "NA", call_rate=call_rate, encoded_alts=",".join(sorted(encoded)) or ".",
               expected_alts_not_encoded=",".join(sorted(missing_alts)) if usable or decoded["record_status"] == "record_not_found" else "unavailable",
               represented_status=represented, record_status=decoded["record_status"], representation=decoded["representation"],
               n_phase_comparable=sum(s in decoded.get('phase_comparable', set()) for s in samples),
               n_phase_unknown=sum(s in decoded.get('phase_unknown', set()) for s in samples),
               n_phase_conflict=sum(s in decoded.get('phase_conflict', set()) for s in samples),
               n_implicit_phase_set=sum(s in decoded.get('implicit_phase_set', set()) for s in samples),
               n_explicit_phase_set=sum(s in decoded.get('explicit_phase_set', set()) for s in samples),
               n_missing_samples=sum(s in decoded["missing"] for s in samples) if usable else "NA",
               n_conflicting_samples=sum(s in decoded["conflicting"] for s in samples) if usable else "NA",
               n_bases_ac_ge_min=n_bases, observed_three_base_in_encoded_calls=int(n_bases >= 3) if n else "NA",
               confirmed_multi_allelic=int(n_bases >= 3) if evaluable else "NA",
               target_alts_retained=int(all(ac[a] >= min_ac for a in expected)) if evaluable else "NA",
               evaluable=int(evaluable), min_ac=min_ac, min_call_rate=min_call_rate,
               interpretation="conditional_target_set_dataset_concordance")
    # Zero is reserved for an encoded nucleotide with genuinely observed calls.
    for b in GENOTYPES_BASES:
        row[b + "_AC"] = ac[b] if n and b in encoded | {target["ref"]} else "NA"
    return row


def genotypes_wilson(k, n):
    if not n: return "NA", "NA"
    z = 1.959963984540054; p = k / n; den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return max(0., center-half), min(1., center+half)


def genotypes_summarize(rows, min_ac, min_call_rate):
    groups = defaultdict(list)
    for row in rows: groups[row["family"], row["group"], row["population"]].append(row)
    result = []
    for (family, group, pop), r in sorted(groups.items()):
        usable = [x for x in r if x["evaluable"]]
        n = len(usable); k = sum(x["confirmed_multi_allelic"] for x in usable)
        kept = sum(x["target_alts_retained"] for x in usable)
        lo, hi = genotypes_wilson(k, n); rlo, rhi = genotypes_wilson(kept, n)
        result.append(dict(family=family, group=group, population=pop, n_targets=len(r),
                           n_all_expected_encoded=sum(x["represented_status"] == "all_expected_alts_encoded" for x in r),
                           n_evaluable=n, n_confirmed_multi_allelic=k, confirmed_fraction=k/n if n else "NA",
                           confirmed_wilson_low=lo, confirmed_wilson_high=hi,
                           n_target_alts_retained=kept, retained_fraction=kept/n if n else "NA",
                           retained_wilson_low=rlo, retained_wilson_high=rhi, min_ac=min_ac, min_call_rate=min_call_rate,
                           ci_method="descriptive_Wilson_not_LD_adjusted",
                           scope="predefined_AF_matched_modern_multi_ALT_targets; not_genome_wide_enrichment"))
    return result


def genotypes_unavailable(status):
    return {"record_status": status, "representation": "unavailable", "encoded_alts": set(),
            "calls": {}, "missing": set(), "conflicting": set(), "phase_comparable": set(), "phase_unknown": set(), "phase_conflict": set()}


def genotypes_plot_concordance(summary, out):
    # A rerun can have no evaluable sites or lose its plotting dependency. Do
    # not leave a figure from a different target set looking like this result.
    figure = out / "genotype_concordance.png"
    figure.unlink(missing_ok=True)
    if not any(r["n_evaluable"] for r in summary): return "not_measured"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return "matplotlib_unavailable"
    families = sorted({r["family"] for r in summary})
    fig, axes = plt.subplots(1, len(families), figsize=(5 * len(families), 4.6), squeeze=False)
    for ax, family in zip(axes[0], families):
        for j, role in enumerate(("case", "control")):
            rr = {(r["population"]): r for r in summary if r["family"] == family and r["group"] == role}
            color = "#c7653d" if role == "case" else "#347b9c"
            for i, pop in enumerate(GENOTYPES_GROUPS):
                r = rr.get(pop); x = i + (j-.5)*.30
                if r and r["n_evaluable"]:
                    y = r["confirmed_fraction"]
                    ax.errorbar(x, y, yerr=[[y-r["confirmed_wilson_low"]], [r["confirmed_wilson_high"]-y]],
                                fmt="o", color=color, capsize=3, label=role if i == 0 else None)
        ax.set_title(family); ax.set_xticks(range(4), GENOTYPES_GROUPS); ax.set_ylim(-.025, 1.025)
        ax.set_ylabel("Fraction with ≥3 observed nucleotide states"); ax.legend(frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Dataset concordance in predefined targets")
    fig.text(.5, .015, "Equal samples per ancestry; descriptive Wilson intervals do not adjust for LD. Missing ALTs are not zeros.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, .94)); fig.savefig(figure, dpi=180); plt.close(fig)
    return "complete"


@transactional("genotypes_manifest.json")
def genotypes_run_genotypes(args):
    if args.min_ac < 1 or not 0 < args.min_call_rate <= 1:
        raise ValueError("--min-ac must be positive and --min-call-rate in (0,1]")
    targets = genotypes_read_targets(args.targets); panel = genotypes_read_panel(args.sample_panel)
    out = Path(args.out).expanduser().resolve(); out.mkdir(parents=True, exist_ok=True)
    sources = {}; missing_sources = {}; index_paths = {}
    for chrom in sorted({r["chrom"] for r in targets}, key=int):
        path = Path(args.vcf_dir).expanduser() / args.vcf_pattern.format(chrom=chrom)
        index = next((Path(str(path)+x) for x in (".tbi", ".csi") if Path(str(path)+x).is_file()), None)
        if not path.is_file(): missing_sources[chrom] = "missing_vcf"
        elif index is None: missing_sources[chrom] = "missing_vcf_index"
        else: sources[chrom] = path; index_paths[chrom] = index
    identity = {"schema": 2, "code": genotypes_metadata(__file__),
                "measurement_code": genotypes_metadata(Path(__file__).with_name("allele_measurement.py")),
                "targets": genotypes_metadata(args.targets), "panel": genotypes_metadata(args.sample_panel),
                "vcfs": {c: genotypes_metadata(p) for c, p in sources.items()}, "indexes": {c: genotypes_metadata(p) for c, p in index_paths.items()},
                "missing_sources": missing_sources, "seed": args.seed, "min_ac": args.min_ac, "min_call_rate": args.min_call_rate,
                "method": "complete_diploid_A_C_G_T_v2; balanced_fixed_samples; split_unique_ALT_residual_REF; implicit_PS_checked"}
    receipt = out / "genotypes_manifest.json"
    if receipt.is_file():
        old = json.loads(receipt.read_text())
        if old.get("status") == "complete" and old.get("identity") == identity and old.get("outputs"):
            if all((out/p).is_file() and genotypes_metadata(out/p) == fp for p, fp in old["outputs"].items()):
                print(f"[herv-genotypes] reuse {out}", flush=True); return old
    # Invalidate the previous success before any refreshed output is written.
    # A failed/interrupted direct invocation must not expose an older completed
    # manifest as evidence for the new target set or input files.
    genotypes_write_json(receipt, {"status": "incomplete", "identity": identity, "outputs": {},
                         "reason": "Refresh started; successful finalization has not occurred"})
    handles = {}; selected = {g: [] for g in GENOTYPES_GROUPS}; eligible_counts = {}; not_run_reason = ""
    sites = []; sample_rows = []; unknown_panel_samples = set()
    try:
        if sources:
            try: import pysam
            except ImportError as e: raise ValueError("pysam is required (included in allele3 env.yml)") from e
            for chrom, path in sources.items():
                handles[chrom] = pysam.VariantFile(str(path))
            common = set.intersection(*(set(v.header.samples) for v in handles.values()))
            unknown_panel_samples = common - set(panel)
            selected, eligible_counts = genotypes_balanced_samples(common, panel, args.seed)
            samples = [s for g in GENOTYPES_GROUPS for s in selected[g]]
            for handle in handles.values(): handle.subset_samples(samples)
        else:
            samples = []; not_run_reason = "no_indexed_target_VCFs" if targets else "empty_targets"
        burden = defaultdict(lambda: [0, 0, 0])  # represented targets, called sites, nonreference copies
        requests = Counter((r["family"], r["group"]) for r in targets)
        for target in targets:
            chrom = target["chrom"]; decoded = None
            if chrom not in handles:
                decoded = genotypes_unavailable(missing_sources.get(chrom, "not_assayed"))
            else:
                vf = handles[chrom]
                contig = chrom if chrom in vf.header.contigs else "chr"+chrom if "chr"+chrom in vf.header.contigs else None
                if contig is None:
                    decoded = genotypes_unavailable("contig_not_in_vcf")
                else:
                    records = [r for r in vf.fetch(contig, target["pos"]-1, target["pos"]) if r.pos == target["pos"]]
                    decoded = genotypes_decode_records(records, target["ref"], samples)
            for g in GENOTYPES_GROUPS:
                sites.append(genotypes_evaluate_population(target, decoded, g, selected[g], args.min_ac, args.min_call_rate))
            represented = decoded["record_status"] == "ok" and set(target["alts"].split(",")) <= decoded["encoded_alts"]
            for g in GENOTYPES_GROUPS:
                for sample in selected[g]:
                    key = sample, g, target["family"], target["group"]
                    z = burden[key]
                    if not represented: continue
                    z[0] += 1
                    if sample in decoded["calls"]:
                        z[1] += 1
                        z[2] += 2 - decoded["calls"][sample][target["ref"]]
        for (sample, pop, family, group), (represented, called, nonref) in sorted(burden.items()):
            sample_rows.append(dict(sample_id=sample, population=pop, family=family, group=group,
                                    n_targets_requested=requests[family, group], n_targets_represented=represented,
                                    called_sites=called, called_alleles=2*called, nonref_allele_copies=nonref if called else "NA",
                                    nonref_per_called_site=nonref/called if called else "NA",
                                    nonref_allele_fraction=nonref/(2*called) if called else "NA",
                                    sample_target_call_rate=called/represented if represented else "NA",
                                    measure="nonreference_diploid_copy_burden_per_called_predefined_target; not_individual_three_alleles"))
    finally:
        for handle in handles.values(): handle.close()
    # A concurrent file update must not leave a falsely reusable result.
    if (genotypes_metadata(args.targets) != identity["targets"]
            or genotypes_metadata(args.sample_panel) != identity["panel"]
            or any(genotypes_metadata(p) != identity["vcfs"][c] for c, p in sources.items())
            or any(genotypes_metadata(p) != identity["indexes"][c] for c, p in index_paths.items())):
        raise ValueError("Genotype inputs changed during analysis")
    summary = genotypes_summarize(sites, args.min_ac, args.min_call_rate)
    balanced = [dict(sample_id=s, population=g, available_common_samples=eligible_counts[g], selected_per_group=len(selected[g]))
                for g in GENOTYPES_GROUPS for s in selected[g]]
    genotypes_write_tsv(out / "genotype_validation.tsv", GENOTYPES_SITE_FIELDS, sites)
    genotypes_write_tsv(out / "cohort_retention.tsv", GENOTYPES_SUMMARY_FIELDS, summary)
    genotypes_write_tsv(out / "sample_burden.tsv", GENOTYPES_BURDEN_FIELDS, sample_rows)
    genotypes_write_tsv(out / "balanced_samples.tsv", ["sample_id", "population", "available_common_samples", "selected_per_group"], balanced)
    plot_status = genotypes_plot_concordance(summary, out)
    names = ["genotype_validation.tsv", "cohort_retention.tsv", "sample_burden.tsv", "balanced_samples.tsv"]
    if plot_status == "complete": names.append("genotype_concordance.png")
    result = {"status": "complete" if sources else "not_run", "not_run_reason": not_run_reason, "identity": identity,
              "n_targets": len(targets), "n_selected_samples": len(balanced), "selected_per_population": {g: len(v) for g, v in selected.items()},
              "available_common_samples_by_population": eligible_counts, "unlabelled_common_samples_excluded": len(unknown_panel_samples),
              "record_status_counts": dict(Counter(r["record_status"] for r in sites)),
              "phase_audit_site_population_counts": {k:sum(r.get(k,0) for r in sites if isinstance(r.get(k),int)) for k in ("n_phase_comparable","n_phase_unknown","n_phase_conflict","n_implicit_phase_set","n_explicit_phase_set")},
              "n_evaluable_site_population_pairs": sum(r["evaluable"] for r in sites),
              "min_ac": args.min_ac, "min_call_rate": args.min_call_rate, "coordinate_build": "GRCh38", "plot_status": plot_status,
              "validation_scope": "dataset_concordance; sample_overlap_with_gnomAD_not_ruled_out; not_independent_replication",
              "threshold_warning": "1KG minAC=2 by default in equal-sized ancestry subsets; not the same threshold/sample size as gnomAD AC5/AC10",
              "confirmed_definition": "all expected ALTs encoded, call rate passes, and at least three ACTUAL nucleotide states each have AC>=min_ac; REF counts must be observed in direct calls or reconstructed from fully called consistent split records",
              "missing_policy": "unencoded expected ALT, absent records and incomplete genotypes are unavailable, never evidence of a zero allele count; nonencoded base AC is NA",
              "split_records": "only unique biallelic SNVs with identical REF; complete diploid GT on every record; REF=2-sumALT only if totalALT<=2; conflicting samples excluded and reported; mixed/duplicate/nonSNV record structures rejected",
              "phase_sets": "All phased calls with no PS share the implicit phase set (VCF specification); equal explicit PS also permits haplotype checks; mixed/mismatched PS or unphased calls remain phase-unknown",
              "uncertainty": "Wilson intervals are descriptive, conditional on tested sites, and do not account for LD; no significance tests",
              "burden": "per-sample nonreference copies normalized by called predefined targets; cannot infer three alleles within a diploid person",
              "outputs": {name: genotypes_metadata(out / name) for name in names}}
    genotypes_write_json(receipt, result)
    print(f"[herv-genotypes] {len(targets)} targets; {len(balanced)} balanced samples; {out / 'genotype_validation.tsv'}", flush=True)
    return result


def genotypes_main(argv=None):
    p = argparse.ArgumentParser(description=_description_herv_genotypes)
    p.add_argument("--targets", required=True)
    p.add_argument("--vcf-dir", required=True)
    p.add_argument("--sample-panel", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--vcf-pattern", default="chr{chrom}.vcf.gz")
    p.add_argument("--seed", type=int, default=20261027)
    p.add_argument("--min-ac", type=int, default=2)
    p.add_argument("--min-call-rate", type=float, default=.9)
    args = p.parse_args(argv)
    try:
        genotypes_run_genotypes(args)
    except (ValueError, OSError, KeyError) as e:
        p.exit(2, f"herv-genotypes: {e}\n")
    return 0


def _entry_herv_genotypes():
    raise SystemExit(genotypes_main())


# ========================================================================
# herv domains
# ========================================================================
_description_herv_domains = "Decompose primary HERV multi-ALT enrichment with external domain annotations.\n\nBED4/6 column 4 must contain reviewed biological labels (e.g. Gag, Pol, Env,\n5LTR, 3LTR). The caller supplies the annotation source and GRCh38 build.\nGeometric repeat ends and '-int' names are never converted to domain labels.\nThe original eligible SNV universe is queried over aligned family intervals,\nincluding single-ALT controls and all unannotated/ambiguous residual sites."
#!/usr/bin/env python3
"""Decompose primary HERV multi-ALT enrichment with external domain annotations.

BED4/6 column 4 must contain reviewed biological labels (e.g. Gag, Pol, Env,
5LTR, 3LTR). The caller supplies the annotation source and GRCh38 build.
Geometric repeat ends and '-int' names are never converted to domain labels.
The original eligible SNV universe is queried over aligned family intervals,
including single-ALT controls and all unannotated/ambiguous residual sites.
"""
from dataclasses import dataclass
import hashlib
import sys
from analysis_integrity import sha256

DOMAINS_SCHEMA = "herv-biological-domains-1.0"
DOMAINS_RESIDUAL = ("unannotated", "ambiguous", "strand_mismatch")
DOMAINS_FIELDS = ["chrom", "block_id", "family", "domain", "sensitivity", "threshold", "context", "af_bin", "an_bin", "n_sites", "n_multi"]


@dataclass(frozen=True)
class domains_Domain:
    chrom: str
    start: int
    end: int
    label: str
    strand: str = "."


def domains_read_domains(path):
    domains, seen, labels = [], set(), set()
    with open_text(path) as f:
        for number, line in enumerate(f, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) not in {4, 6}:
                raise ValueError(f"{path}:{number}: use BED4 or BED6 with a reviewed domain label in column 4")
            chrom, start, end, label = fields[:4]
            strand = fields[5] if len(fields) == 6 else "."
            start, end = int(start), int(end)
            if start < 0 or end <= start or not label or label.lower() in {*DOMAINS_RESIDUAL, "total"} or strand not in {"+", "-", "."}:
                raise ValueError(f"{path}:{number}: invalid domain interval, reserved label, or strand")
            domain = domains_Domain(canonical_chrom(chrom), start, end, label, strand)
            if domain not in seen:
                domains.append(domain)
                seen.add(domain)
                labels.add(label)
    if not domains:
        raise ValueError("No domain annotation intervals")
    return sorted(domains, key=lambda d: (scan.chrom_sort_key(d.chrom), d.start, d.end, d.label, d.strand)), sorted(labels)


class domains_DomainSweep:
    """Return unique domain annotations for sorted 1-based nucleotide positions."""
    def __init__(self, domains):
        self.by_chrom = defaultdict(list)
        for domain in domains:
            self.by_chrom[domain.chrom].append(domain)
        for chrom in self.by_chrom:
            self.by_chrom[chrom].sort(key=lambda d: (d.start, d.end, d.label, d.strand))
        self.chrom, self.index, self.active, self.last = None, 0, [], {}

    def locate(self, chrom, pos1):
        chrom, pos = canonical_chrom(chrom), int(pos1) - 1
        if pos < 0 or pos < self.last.get(chrom, -1):
            raise ValueError("Domain sweep requires sorted positive positions")
        self.last[chrom] = pos
        if self.chrom != chrom:
            self.chrom, self.index, self.active = chrom, 0, []
        arr = self.by_chrom[chrom]
        while self.index < len(arr) and arr[self.index].start <= pos:
            self.active.append(arr[self.index])
            self.index += 1
        self.active = [d for d in self.active if d.end > pos]
        return tuple(self.active)


def domains_assign_domain(hits, repeat_strands):
    """One category per site/family, retaining conflicting evidence explicitly."""
    if not hits:
        return "unannotated"
    usable = [d for d in hits if d.strand == "." or d.strand in repeat_strands]
    if not usable:
        return "strand_mismatch"
    # An overlapping, discordant-strand annotation remains ambiguity evidence.
    if len(usable) != len(hits) or len(repeat_strands) != 1:
        return "ambiguous"
    labels = {d.label for d in usable}
    return next(iter(labels)) if len(labels) == 1 else "ambiguous"


def domains_scan_chrom(vcf, chrom, annotation, options, args, domains):
    """Query all eligible sites in selected aligned-family unions once per chr."""
    out = Path(args.out) / "cache"
    out.mkdir(parents=True, exist_ok=True)
    fasta = str(args.fasta or options["fasta"])
    mask_paths = {name: getattr(args, name + "_bed", None) or options.get(option)
                  for name, option in (("segdup", "architecture_segdup_bed"), ("lcr", "architecture_lcr_bed"),
                                       ("low_mappability", "architecture_low_mappability_bed"))}
    ld_path = Path(options["ld_block_dir"]) / f"EUR.{options.get('architecture_build', '38')}.bed"
    indexes = [str(vcf) + suffix for suffix in (".tbi", ".csi") if Path(str(vcf) + suffix).is_file()]
    if not indexes:
        raise ValueError(f"Indexed VCF required: {vcf}")
    header = modern_read_indexed_header(vcf, args.bcftools)
    scan.validate_count_header(header)
    fields = ["AC", "AN"] + [f for f in ("AF", "AS_FilterStatus", "lcr", "segdup") if f in header["info"]]
    elements = scan.parse_joined_rmsk(annotation, chromosome=chrom)
    families = set(args.families.split(","))
    selected = [e for e in elements if e.rep_class == "LTR" and e.family in families]
    intervals = scan.merge_intervals([s for e in selected for s in (e.segments or ((e.start, e.end),))])
    repeats = scan.ERVOverlapSweep(elements)
    domain_sweep = domains_DomainSweep(domains)
    masks = {name: IntervalMask(path, [chrom]) for name, path in mask_paths.items() if path}
    strict_complete = all(name in masks or name in header["info"] for name in ("lcr", "segdup")) and "low_mappability" in masks
    callable_mask = IntervalMask(options["callable_bed"], [chrom]) if options.get("callable_bed") else None
    eur = scan.BlockSweep(scan.read_blocks(ld_path, "EUR", options.get("architecture_build", "38"))[0])
    contigs = [c for c in header["contigs"] if canonical_chrom(c) == chrom]
    if len(contigs) != 1:
        raise ValueError(f"No unique VCF contig for chromosome {chrom}")
    bed = out / f"chr{chrom}.regions.bed"
    with bed.open("w") as f:
        for start, end in intervals:
            f.write(f"{contigs[0]}\t{start}\t{end}\n")
    counts, qc, previous_block = defaultdict(lambda: [0, 0]), Counter(), None
    path = out / f"chr{chrom}.domain_strata.tsv.gz"
    partial = Path(str(path) + ".partial")
    with FastaContext(fasta) as fa, gzip.open(partial, "wt", newline="", compresslevel=2) as f:
        for d in domains:
            if d.chrom == chrom and (chrom not in fa.lengths or d.end > fa.lengths[chrom]):
                raise ValueError(f"Domain BED interval outside input reference: {d}")
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(DOMAINS_FIELDS)
        def flush():
            writer.writerows([*key, *value] for key, value in sorted(counts.items()))
            counts.clear()
        records = modern_query_regions(vcf, fields, bed, args.bcftools) if intervals else iter(())
        started = time.monotonic()
        for (c, pos), group in itertools.groupby(records, key=lambda r: (canonical_chrom(r["chrom"]), r["pos"])):
            if c != chrom:
                raise ValueError("Indexed domain query returned an unexpected chromosome")
            if callable_mask and not callable_mask.contains(c, pos):
                qc["outside_callable_mask"] += 1
                continue
            overlap = repeats.locate(c, pos)
            by_family = defaultdict(set)
            for e in overlap:
                if e.rep_class == "LTR" and e.family in families:
                    by_family[e.family].add(e.strand)
            if not by_family:
                raise ValueError("Indexed domain query returned a site outside selected aligned family bases")
            ref, alleles, _, flags = scan.group_alleles(group, header, [], options, qc, None, [])
            observed = sorted((a for a in alleles.values() if a["status"] == "observed_pass"), key=lambda a: (-a["af"], -a["ac"], a["alt"]))
            if not observed:
                continue
            ctx = fa.context(c, pos).upper()
            if len(ctx) != 3 or ctx[1] != ref:
                raise ValueError(f"REF/FASTA mismatch: {c}:{pos}")
            context, _ = scan.canonical_context(ctx, ref)
            primary = observed[0]
            strata = context, af_bin(primary["af"]), an_bin(primary["an"])
            sensitivity = ["all"]
            if strict_complete and not modern_strict_site_excluded(c, pos, flags, overlap, masks):
                sensitivity.append("strict")
            block = eur.locate(c, pos)
            block_id = block.block_id if block else "outside_partition"
            if previous_block is not None and block_id != previous_block:
                flush()
            previous_block = block_id
            hits = domain_sweep.locate(c, pos)
            for family, strands in by_family.items():
                label = domains_assign_domain(hits, strands)
                qc[f"{family}:{label}:eligible"] += 1
                for sens in sensitivity:
                    for threshold in (5, 10):
                        key = c, block_id, family, label, sens, threshold, *strata
                        counts[key][0] += 1
                        counts[key][1] += int(sum(a["ac"] >= threshold for a in observed) >= 2)
            qc["eligible_sites_in_family_union"] += 1
            if time.monotonic() - started > 60:
                print(f"[domains chr{chrom}] {qc['eligible_sites_in_family_union']:,} eligible SNVs", flush=True)
                started = time.monotonic()
        flush()
    os.replace(partial, path)
    inputs = [vcf, indexes[0], annotation, fasta, fasta + ".fai", ld_path, *[p for p in mask_paths.values() if p]]
    if options.get("callable_bed"):
        inputs.append(options["callable_bed"])
    dump_json(out / f"chr{chrom}.query.json", dict(status="complete", qc=dict(qc),
        queried_aligned_bp=sum(b-a for a, b in intervals), strict_masks_complete=strict_complete,
        inputs=[modern_file_state(p) for p in inputs], domain_bed_sha256=sha256(args.domain_bed),
        output=modern_file_state(path), query_scope="all eligible SNVs in selected RepeatMasker aligned-family union"))
    return path


def domains_aggregate(paths, background):
    # Vector: observed, partial expected, eligible, unsupported observed, unsupported eligible.
    by_domain = defaultdict(lambda: defaultdict(lambda: np.zeros(5)))
    raw, outside = defaultdict(lambda: np.zeros(2, dtype=np.int64)), Counter()
    for path in paths:
        for r in modern_rows(path):
            chrom, n, m = canonical_chrom(r["chrom"]), int(r["n_sites"]), int(r["n_multi"])
            if n < 0 or m < 0 or m > n:
                raise ValueError("Invalid domain N/M counters")
            prefix = r["family"], r["sensitivity"], int(r["threshold"])
            scopes = ["all_family_sites"]
            if r["block_id"] != "outside_partition":
                scopes.append("eur_ld_partition")
            for scope in scopes:
                raw[chrom, *prefix, scope, r["context"], r["af_bin"], r["an_bin"]] += [n, m]
            if r["block_id"] == "outside_partition":
                outside[prefix] += n
                continue
            rate = background.rate(chrom, background.key(r))
            if rate is not None and not 0 <= rate <= 1:
                raise ValueError("LOCO rate outside [0,1]")
            value = np.array([m, n*rate, n, 0, 0]) if rate is not None else np.array([m, 0, n, m, n])
            by_domain[*prefix, r["domain"]][r["block_id"]] += value
    return by_domain, raw, outside


def domains_summarize_domains(data, outside, labels, families, bootstrap=2000, seed=20261102, min_blocks=20):
    if bootstrap < 2000 or min_blocks < 2:
        raise ValueError("Domain inference requires >=2000 bootstrap draws and >=2 minimum blocks")
    categories = sorted(set(labels)) + list(DOMAINS_RESIDUAL)
    rows, identities, block_rows = [], [], []
    prefixes = sorted({key[:3] for key in data} | set(outside))
    for prefix in prefixes:
        total_blocks = defaultdict(lambda: np.zeros(5))
        for category in categories:
            for block, value in data.get((*prefix, category), {}).items():
                total_blocks[block] += value
        total = sum(total_blocks.values(), np.zeros(5))
        overall_missing = int(total[4])
        total_expected = None if overall_missing else float(total[1])
        group_rows = []
        for category in [*categories, "total"]:
            blocks = total_blocks if category == "total" else data.get((*prefix, category), {})
            observed, partial_expected, eligible, missing_o, missing_n = sum(blocks.values(), np.zeros(5))
            # Suppress every weighted part when full-scope LOCO support is incomplete.
            expected = None if overall_missing else float(partial_expected)
            oe = float(observed / expected) if expected else None
            enough = not overall_missing and len(blocks) >= min_blocks and observed > 0 and expected
            ci = scan.bootstrap_statistics({b: v[:2] for b, v in sorted(blocks.items())}, bootstrap, seed) if enough else (None, None, None)
            status = ("insufficient_loco_training" if overall_missing else "no_eligible_sites" if not eligible else
                      "zero_expected" if not expected else "insufficient_blocks_or_events" if not enough else "complete")
            row = dict(family=prefix[0], sensitivity=prefix[1], threshold=prefix[2], domain=category,
                partition=category, scope="eur_ld_partition", observed_multi=int(observed),
                eligible_sites=int(eligible), n_without_loco=int(missing_n), multi_without_loco=int(missing_o),
                total_without_loco=overall_missing, partial_expected_multi=float(partial_expected),
                expected_multi=expected, oe=oe, ci_low=ci[0], ci_high=ci[1], bootstrap_p=ci[2], fdr=None,
                n_blocks=len(blocks), outside_partition_sites=outside.get(prefix, 0) if category == "total" else None,
                expected_weight=expected / total_expected if expected is not None and total_expected else None,
                oe_contribution=float(observed / total_expected) if total_expected else None,
                excess_contribution=float((observed - expected) / total_expected) if total_expected and expected is not None else None,
                status=status)
            group_rows.append(row)
            for block, value in sorted(blocks.items()):
                block_rows.append(dict(family=prefix[0], sensitivity=prefix[1], threshold=prefix[2], domain=category,
                    block_id=block, observed_multi=int(value[0]), partial_expected_multi=float(value[1]),
                    n_eligible=int(value[2]), multi_without_loco=int(value[3]), n_without_loco=int(value[4])))
        total_row, parts = group_rows[-1], group_rows[:-1]
        fields = ("observed_multi", "eligible_sites", "n_without_loco", "multi_without_loco", "partial_expected_multi")
        deltas = {field: sum(r[field] for r in parts) - total_row[field] for field in fields}
        if any(not math.isclose(delta, 0, abs_tol=max(1e-8, abs(total_row[field])*1e-10)) for field, delta in deltas.items()):
            raise ValueError(f"Domain partition failed exact count accounting: {prefix}")
        identities.append(dict(family=prefix[0], sensitivity=prefix[1], threshold=prefix[2],
            overall_oe=total_row["oe"], expected_weighted_oe=sum(r["oe_contribution"] for r in parts) if total_expected else None,
            expected_weight_sum=sum(r["expected_weight"] for r in parts) if total_expected else None,
            observed_delta=deltas["observed_multi"], eligible_delta=deltas["eligible_sites"], expected_delta=deltas["partial_expected_multi"],
            total_without_loco=overall_missing, status="raw_accounting_only_loco_unavailable" if overall_missing else "exact_accounting"))
        rows.extend(group_rows)
    # Biological labels form a distinct exploratory testing family; total and
    # ambiguity residuals remain reported but do not inflate this label FDR.
    tests = defaultdict(list)
    for row in rows:
        if row["domain"] in labels and row["bootstrap_p"] is not None:
            tests[row["threshold"], row["sensitivity"]].append(row)
    for tested in tests.values():
        for row, q in zip(tested, scan.bh_fdr([r["bootstrap_p"] for r in tested])):
            row["fdr"] = q
    return rows, identities, block_rows


@transactional("domain_manifest.json", input_files=("domain_bed",))
def domains_run(args):
    if args.domain_build not in {"38", "GRCh38", "hg38"} or not args.domain_source.strip():
        raise ValueError("A GRCh38 annotation and explicit domain source are required")
    families = args.families.split(",")
    if not all(families) or len(set(families)) != len(families):
        raise ValueError("Families must be unique nonempty exact RepeatMasker labels")
    domains, labels = domains_read_domains(args.domain_bed)
    original_hash = sha256(args.domain_bed)
    source_run, annotation, tasks, source_chromosomes = modern_targeted_inputs(args)
    options = source_run["options"]
    if not options.get("ld_block_dir") or str(options.get("architecture_build", "38")) != "38":
        raise ValueError("Primary domain decomposition requires the original EUR GRCh38 LD partition")
    if not {5, 10} <= {int(t) for t in options.get("thresholds", [])}:
        raise ValueError("Original counters must contain AC5 and AC10")
    files = [domains_scan_chrom(vcf, chrom, annotation / f"chr{chrom}.tsv.gz", options, args, domains) for chrom, vcf in tasks]
    used = {c for c, _ in tasks}
    background = modern_PrimaryBackground(modern_architecture_paths(args.modern_dir, "global"), families)
    data, raw, outside = domains_aggregate(files, background)
    rows, identities, block_rows = domains_summarize_domains(data, outside, labels, families, args.bootstrap, args.seed, args.min_blocks)
    scope = "all_22_autosomes" if used == MODERN_AUTOSOMES else "selected_chromosomes_subset"
    reconciliation, mismatches, _ = modern_reconcile_spatial_counters(args.modern_dir, raw, background, used, families,
        comparison_scope=scope, primary_background=background)
    counter_totals = primary_scoped_counter_totals(modern_architecture_paths(args.modern_dir, "exposure"), background, families, used)
    total_audit, total_failures = primary_reconcile_spatial_to_counters(counter_totals, rows)
    out = Path(args.out)
    for row in rows:
        row["comparison_scope"] = scope
    modern_write_rows(out / "domain_enrichment.tsv", rows)
    modern_write_rows(out / "domain_accounting.tsv", identities)
    modern_write_rows(out / "domain_block_counts.tsv.gz", block_rows)
    modern_write_rows(out / "domain_counter_reconciliation.tsv", reconciliation)
    modern_write_rows(out / "domain_counter_mismatches.tsv.gz", mismatches)
    modern_write_rows(out / "domain_total_reconciliation.tsv", total_audit)
    modern_write_rows(out / "domain_annotations.tsv.gz", [dict(chrom=d.chrom, start=d.start, end=d.end, domain=d.label, strand=d.strand) for d in domains])
    if sha256(args.domain_bed) != original_hash:
        raise ValueError("Domain annotation changed during analysis")
    if mismatches or total_failures:
        raise ValueError("Domain totals disagree with original same-scope counters; diagnostic staging is retained")
    manifest = dict(status="complete" if rows else "not_estimable", schema=DOMAINS_SCHEMA, domain_source=args.domain_source,
        domain_bed=str(Path(args.domain_bed).resolve()), domain_bed_sha256=original_hash, domain_build="GRCh38", domain_labels=labels,
        source_annotation_intervals=len(domains), queried_chromosomes=sorted(used, key=int),
        annotation_chromosomes_not_queried=sorted({d.chrom for d in domains} - used, key=scan.chrom_sort_key),
        source_chromosomes=sorted(source_chromosomes, key=int), comparison_scope=scope,
        families=families, original_scan=str(Path(args.modern_dir).resolve()),
        source_options=options, bootstrap_replicates=args.bootstrap, seed=args.seed, min_blocks=args.min_blocks,
        estimator="Original count-supported multi-ALT definition; LOCO trinucleotide/primary ALT AF/AN, minimum 1; (M+.5)/(N+1)",
        inference="EUR LD-block percentile bootstrap; background fitted once and treated as fixed; label FDR within AC/sensitivity",
        partition="One domain per eligible site/family; unannotated, conflicting labels, and strand mismatch retained; gaps/flanks excluded",
        accounting="Every domain and residual sums to the same eligible N and multi-ALT O; total O/E is weighted by domain E",
        caveats=["External domain labels are supplied evidence, not validated by their names; no geometric head is called a 5-prime LTR.",
                 "Data describe observed multi-ALT enrichment, not the time or mechanism of HERV insertion or nucleotide mutation.",
                 "Cross-family overlap is retained independently; different families cannot be added.",
                 "X/Y annotations are disclosed as unqueried; existing X/Y allele3 analyses remain separate because the EUR autosomal LD inference does not apply.",
                 "Different HERV age, local GC and recombination distributions remain potential confounders."])
    dump_json(out / "domain_manifest.json", manifest)
    return manifest


def domains_parser():
    p = argparse.ArgumentParser(description=_description_herv_domains)
    p.add_argument("--modern-dir", required=True)
    p.add_argument("--domain-bed", required=True)
    p.add_argument("--domain-source", required=True, help="Publication DOI/version or reviewed annotation provenance")
    p.add_argument("--domain-build", required=True, choices=("38", "GRCh38", "hg38"))
    p.add_argument("--out", required=True)
    p.add_argument("--families", default="ERV1,ERVK")
    p.add_argument("--chromosomes", help="Optional original autosomal chromosome subset; never called a genome-wide result")
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20261102)
    p.add_argument("--min-blocks", type=int, default=20)
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--fasta")
    p.add_argument("--erv-joined")
    for name in ("segdup", "lcr", "low-mappability"):
        p.add_argument("--" + name + "-bed")
    return p


def domains_main(argv=None):
    args = domains_parser().parse_args(argv)
    try:
        result = domains_run(args)
        print(json.dumps({k: result[k] for k in ("status", "queried_chromosomes", "domain_labels")}, indent=2))
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        print(f"herv-domains: {exc}", file=sys.stderr)
        return 2
    return 0


def _entry_herv_domains():
    raise SystemExit(domains_main())


# ========================================================================
# herv archaic
# ========================================================================
_description_herv_archaic = 'New GU/IBDmix cohorts against the completed, read-only gnomAD background.\n\nEach cohort retains its own analyzed-person denominator. Query the union of its\nlifted tracts with the original SNV QC and ArchitectureAccumulator; never infer\nthe eligible denominator from the cached multi-ALT-only table.'
#!/usr/bin/env python3
"""New GU/IBDmix cohorts against the completed, read-only gnomAD background.

Each cohort retains its own analyzed-person denominator. Query the union of its
lifted tracts with the original SNV QC and ArchitectureAccumulator; never infer
the eligible denominator from the cached multi-ALT-only table.
"""

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
archaic_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(archaic_paths_root / "0f"))
from workflow_paths import configured_path
import fcntl
import fnmatch
import shutil
import tarfile
from types import SimpleNamespace

from analysis_integrity import identity, validate_manifest

ARCHAIC_SCHEMA = "allele3-archaic-3"
ARCHAIC_DEFAULT_REFERENCES = ("Altai", "Chagyr", "Denisova", "Denisova25", "Vindija")
ARCHAIC_AUTOSOMES = tuple(str(i) for i in range(1, 23))
ARCHAIC_CHROMOSOMES = (*ARCHAIC_AUTOSOMES, "X", "Y")
ARCHAIC_CODE = ('herv.py', 'herv.py', "scan.py", "scan_store.py", "0.common.py", "analysis_integrity.py")


def archaic_sample_ids(path):
    result = set()
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            value = line.split()[0].lstrip("#")
            if value.lower() in {"sample", "sample_id", "iid", "id", "fid"}:
                continue
            if value in result:
                raise ValueError(f"Duplicate analyzed sample in {path}: {value}")
            result.add(value)
    if not result:
        raise ValueError(f"Empty analyzed sample panel: {path}")
    return result


def archaic_metadata(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        parts = line.replace("\\t", "\t").split("\t", 1)
        if len(parts) == 2:
            result[parts[0]] = parts[1]
    return result


def archaic_archive_complete(path):
    """Older GU publications keep .complete in the leading archive metadata.

    Inspect bounded headers only: never decompress a large tract member during
    discovery. GU's publisher puts completion/provenance before final tables.
    """
    with tarfile.open(path, 'r|gz') as archive:
        consumed = 0
        for member in archive:
            if member.name.removeprefix('./') == '.complete' and member.isfile():
                return True
            consumed += member.size + 512
            if consumed > 1024*1024:
                break
    return False


def archaic_inspect_dataset(path, module="auto", chromosomes=ARCHAIC_AUTOSOMES, references=ARCHAIC_DEFAULT_REFERENCES):
    """Metadata and file inventory only; do not open raw tracts at preflight."""
    from scan import track_references
    selected_references = track_references(references)
    path = Path(path).resolve()
    modules = [p for p in (path/"ibdmix", path/"ibdmix-GPU") if p.is_dir()]
    if module != "auto":
        modules = [path/module] if (path/module).is_dir() else []
    if len(modules) != 1:
        raise ValueError(f"Expected one IBDmix module at {path}/ibdmix or {path}/ibdmix-GPU; found {[str(p) for p in modules]}; select --archaic-module")
    source = modules[0]
    inputs, files, panels = [], [], []
    cohort, contract = None, None
    for chrom in chromosomes:
        directory = source/f"chr{chrom}"
        native = sorted((directory/"final").glob("all_archaic_refs*.segments.tsv.gz"))
        if len(native) > 1:
            raise ValueError(f"Ambiguous finalized tract versions: {directory}")
        tract = native[0] if native else directory/"ibdmix.raw.tar.gz"
        panel, meta = directory/"samples/ALL.txt", directory/"run.meta.tsv"
        required = [tract, panel, meta]
        missing = [str(p) for p in required if not p.is_file()]
        if missing:
            raise ValueError("Missing finalized chromosome inputs: " + ", ".join(missing))
        # Published GU runs retain a publication marker after unpacked .complete
        # moves into their archive. Live workspaces instead retain .complete.
        markers = [p for p in (directory/".published-content", directory/".complete") if p.is_file()]
        if not markers and not (tract.name == 'ibdmix.raw.tar.gz' and archaic_archive_complete(tract)):
            raise ValueError(f"No GU completion/publication marker: {directory}")
        data = archaic_metadata(meta)
        build = data.get("genome_build", "").removeprefix("GRCh").removeprefix("b")
        if build not in {"37", "38"} or data.get("coordinate_system") != "0-based-half-open":
            raise ValueError(f"Unknown genome build or coordinate convention: {meta}")
        refs = tuple(sorted(data.get("refs", "").split()))
        if not refs:
            raise ValueError(f"Missing archaic reference list: {meta}")
        missing_refs = set(selected_references) - set(refs)
        if missing_refs:
            raise ValueError(f"Requested references were not tested in {meta}: {sorted(missing_refs)}")
        lod_cut, len_cut = float(data["lod_cut"]), int(data["len_cut"])
        if not math.isfinite(lod_cut) or lod_cut < 0 or len_cut < 1:
            raise ValueError(f"Invalid native tract thresholds: {meta}")
        if not data.get("pipeline_version"):
            raise ValueError(f"Missing GU caller provenance: {meta}")
        current = (build, refs, lod_cut, len_cut,
                   tuple((key, data.get(key, "not_recorded")) for key in
                         ("pipeline_version", "profile", "background_filter", "denisovan_call_definition",
                          "emit_lod_cut", "minor_allele_count", "archaic_error", "modern_error_max", "modern_error_proportion")))
        if contract is not None and current != contract:
            raise ValueError(f"Build/reference/filter contract differs between chromosomes: {source}")
        contract = current
        people = archaic_sample_ids(panel)
        if cohort is not None and people != cohort:
            raise ValueError(f"Analyzed sample panels differ between chromosomes: {source}")
        cohort = people
        inputs.extend(required + markers)
        if (directory/"cache.meta.tsv").is_file():
            inputs.append(directory/"cache.meta.tsv")
        files.append(str(tract)); panels.append(str(panel))
    return dict(dataset=path.name, source=str(source), status="ready", build=contract[0],
                references=list(contract[1]), selected_references=list(selected_references),
                caller_contract=dict(contract[4]), native_min_lod=contract[2], native_min_length=contract[3],
                n_samples=len(cohort), sample_panel=panels[0], chromosomes=list(chromosomes),
                files=files, inputs=[str(p) for p in inputs])


def archaic_chromosome_groups(selection="auto"):
    """Keep independent sex-chromosome panels and backgrounds out of autosomal LOCO."""
    tokens = [x.strip() for x in selection.split(",")]
    if selection == "auto":
        return [("autosomes", ARCHAIC_AUTOSOMES), ("X", ("X",)), ("Y", ("Y",))]
    if not tokens or any(x not in {"autosomes", "X", "Y"} for x in tokens) or len(set(tokens)) != len(tokens):
        raise ValueError("Chromosomes must be auto or a comma-separated selection of autosomes,X,Y")
    return [(x, ARCHAIC_AUTOSOMES if x == "autosomes" else (x,)) for x in tokens]


def archaic_discover(root, datasets="auto", module="auto", chromosomes="auto", references=ARCHAIC_DEFAULT_REFERENCES,
             source_build="auto"):
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"GU source directory is missing: {root}")
    if source_build not in {"auto", "37", "38"}:
        raise ValueError("source_build must be auto, 37 or 38")
    # QiMing stores IBDmix next to its genotype build. Retain stable cohort
    # labels (1kg, ukb013) for downstream tables and explicit dataset selectors.
    candidates = defaultdict(list)
    builds = ("37", "38") if source_build == "auto" else (source_build,)
    def has_module(path):
        return any((path / m).is_dir() for m in ("ibdmix", "ibdmix-GPU"))
    for path in sorted(root.iterdir()):
        if not path.is_dir() or not (path.name in {"1kg", "ukb"} or re.fullmatch(r"ukb\d+", path.name)):
            continue
        if has_module(path):
            candidates[path.name].append(path)
        for number in builds:
            build = path / number
            if has_module(build):
                candidates[path.name].append(build)
            if path.name == "ukb" and build.is_dir():
                for cohort in sorted(build.glob("cohort*")):
                    if cohort.is_dir() and re.fullmatch(r"cohort\d+", cohort.name) and has_module(cohort):
                        candidates["ukb" + cohort.name[len("cohort"):]].append(cohort)
    names = sorted(candidates)
    auto = datasets == "auto"
    if not auto:
        requested = []
        for pattern in datasets.split(","):
            matches = [n for n in names if fnmatch.fnmatchcase(n, pattern.strip())]
            if not matches:
                raise ValueError(f"No GU dataset matches {pattern!r} under {root}; expected 1kg/37|38/ibdmix or ukb/37|38/cohortNNN/ibdmix; source-build={source_build}")
            requested.extend(matches)
        names = sorted(set(requested))
    records = []
    for name in names:
        paths = candidates[name]
        if len(paths) != 1:
            raise ValueError(f"Ambiguous GU cohort {name}: {[str(p) for p in paths]}; select --source-build 37|38 (pipeline: --archaic-source-build), or one legacy/build layout")
        for scope, selected in archaic_chromosome_groups(chromosomes):
            identifier = name if scope == "autosomes" else f"{name}.chr{scope}"
            try:
                record = archaic_inspect_dataset(paths[0], module, selected, references)
                layout_build = paths[0].name if paths[0].name in {"37", "38"} else paths[0].parent.name
                if layout_build in {"37", "38"} and record["build"] != layout_build:
                    raise ValueError(f"Caller metadata build {record['build']} conflicts with source directory build {layout_build}")
                if source_build != "auto" and record["build"] != source_build:
                    raise ValueError(f"Caller metadata build {record['build']} conflicts with requested source-build {source_build}")
                record.update(dataset=identifier, cohort=name, chromosome_scope=scope)
                records.append(record)
            except (OSError, ValueError, KeyError, tarfile.TarError, EOFError) as exc:
                if not auto and not (chromosomes == "auto" and scope in {"X", "Y"}):
                    raise ValueError(f"Requested dataset {identifier} is not ready: {exc}") from exc
                records.append(dict(dataset=identifier, chromosome_scope=scope, status="not_available", reason=str(exc), inputs=[]))
    if not any(r["status"] == "ready" for r in records):
        raise ValueError(f"No complete selected GU datasets under {root}; " +
                         "; ".join(f"{r['dataset']}: {r['reason']}" for r in records))
    return records


def archaic_modern_inputs(modern):
    return modern_files(modern, "windows")


def archaic_file_identities(paths):
    return [identity(p) for p in dict.fromkeys(str(Path(p).resolve()) for p in paths)]


def archaic_source_identity(dataset, args, common):
    if args.min_lod < dataset["native_min_lod"] or args.min_length < dataset["native_min_length"]:
        raise ValueError(f"{dataset['dataset']}: cannot relax filters below the finalized GU calls")
    extra = [args.chain] if dataset["build"] != "38" else []
    executables = [shutil.which(args.bcftools), shutil.which(args.liftover)]
    return dict(schema=ARCHAIC_SCHEMA, dataset=dataset,
                settings={k: getattr(args, k) for k in ("min_lod", "min_length", "carrier_thresholds", "bootstrap", "seed", "min_training_sites")},
                inputs=archaic_file_identities([*dataset["inputs"], *common, *extra,
                                        *[p for p in executables if p],
                                        *[Path(__file__).with_name(p) for p in ARCHAIC_CODE]]))


def archaic_scan_tract_chrom(vcf, chrom, annotation, track, dataset, out, options, bcftools):
    """Original per-locus QC and exposure hook, restricted to the tract union."""
    import scan
    from allele_common import FastaContext, IntervalMask, canonical_chrom
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    header = modern_read_indexed_header(vcf, bcftools)
    scan.validate_count_header(header)
    fields = ["AC", "AN"] + [f for f in ("AF", "AS_FilterStatus", "lcr", "segdup") if f in header["info"]]
    tracks = scan.read_label_bed(track, chromosome=chrom)
    intervals = scan.merge_intervals([(x.start, x.end) for x in tracks])
    raw = [c for c in header["contigs"] if canonical_chrom(c) == chrom]
    if len(raw) != 1:
        raise ValueError(f"No unique VCF contig for chr{chrom}")
    bed = out/f"chr{chrom}.regions.bed"
    with bed.open("w") as handle:
        for lo, hi in intervals:
            handle.write(f"{raw[0]}\t{lo}\t{hi}\n")
    masks = {name: IntervalMask(options[key], [chrom]) for name, key in
             (("lcr", "architecture_lcr_bed"), ("segdup", "architecture_segdup_bed"),
              ("low_mappability", "architecture_low_mappability_bed")) if options.get(key)}
    callable_mask = IntervalMask(options["callable_bed"], [chrom]) if options.get("callable_bed") else None
    blocks = []
    if chrom in ARCHAIC_AUTOSOMES:
        blocks, _ = scan.read_blocks(Path(options["ld_block_dir"])/f"EUR.{options['architecture_build']}.bed",
                                    "EUR", options["architecture_build"])
    elements = scan.parse_joined_rmsk(annotation, chromosome=chrom)
    prefix = out/f"chr{chrom}"
    counter = scan.ArchitectureAccumulator({"EUR": blocks}, elements,
        thresholds=options["thresholds"], out_prefix=prefix, masks=masks,
        flag_annotations={"lcr", "segdup"} & set(header["info"]),
        tracks={dataset: tracks}, track_only=True)
    qc = Counter()
    try:
        with FastaContext(options["fasta"]) as fa:
            if any(hi > fa.lengths[chrom] for _, hi in intervals):
                raise ValueError(f"Lifted tract outside reference chromosome {chrom}")
            records = modern_query_regions(vcf, fields, bed, bcftools) if intervals else iter(())
            for (c, pos), group in itertools.groupby(records, key=lambda r: (canonical_chrom(r["chrom"]), r["pos"])):
                if c != chrom:
                    raise ValueError("Unexpected chromosome in targeted archaic query")
                ref, alts, _, flags = scan.group_alleles(group, header, [], options, qc, None, [])
                if ref is None or not alts or (callable_mask and not callable_mask.contains(c, pos)):
                    continue
                observed = sorted((a for a in alts.values() if a["status"] == "observed_pass"),
                                  key=lambda a: (-a["af"], -a["ac"], a["alt"]))
                if not observed:
                    continue
                ctx = fa.context(c, pos).upper()
                if len(ctx) != 3 or ctx[1] != ref:
                    raise ValueError(f"REF/FASTA mismatch: {c}:{pos}")
                context, _ = scan.canonical_context(ctx, ref)
                primary = observed[0]
                counter.observe(c, pos, len(observed), context, primary["af"], primary["an"],
                    threshold_counts={t: sum(a["ac"] >= t for a in observed) for t in options["thresholds"]},
                    region_flags=flags)
                qc["eligible_sites"] += 1
            counter.write(prefix)
    finally:
        counter.close()
    audit = dict(chrom=chrom, queried_bp=sum(b-a for a, b in intervals), intervals=len(intervals), qc=dict(qc))
    (out/f"chr{chrom}.audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    return Path(str(prefix)+".exposure_strata.tsv.gz"), audit


@transactional("manifest.json")
def archaic_analyze_dataset(args):
    import scan
    out = Path(args.out)
    dataset = args.dataset
    track_args = SimpleNamespace(dir=dataset["source"], format="ibdmix", out=str(out/"tracks"),
        sample_panel=dataset["sample_panel"], n_samples=None, build=dataset["build"], target_build="38",
        chain=args.chain, liftover=args.liftover, min_lod=args.min_lod, min_length=args.min_length,
        carrier_thresholds=args.carrier_thresholds, coordinates="bed0",
        files=dataset.get("files"), chromosomes=dataset.get("chromosomes"),
        references=dataset["selected_references"])
    track_result = scan.build_tracks(track_args)
    files, audits = [], []
    for chrom, vcf in args.tasks:
        print(f"[archaic {dataset['dataset']}] indexed chr{chrom}", flush=True)
        path, audit = archaic_scan_tract_chrom(vcf, chrom, Path(args.annotation)/f"chr{chrom}.tsv.gz",
            track_result["output"], dataset["dataset"], out/"counts", args.options, args.bcftools)
        files.append(path); audits.append(audit)
    scope = dataset.get("chromosome_scope", "autosomes")
    # Keep autosomal LOCO stable, and use a separate within-chromosome
    # descriptive background for X/Y. Permanent counters remain read-only.
    with tempfile.TemporaryDirectory(prefix="archaic-autosomes-", dir=out) as temporary:
        background = Path(temporary)/"global.tsv.gz"
        with gzip.open(background, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, scan.HEADERS["global"], delimiter="\t")
            writer.writeheader()
            for source in modern_architecture_paths(args.modern_dir, "global"):
                with gzip.open(source, "rt") as rows:
                    for row in csv.DictReader(rows, delimiter="\t"):
                        if scan.canonical_chrom(row["chrom"]) in (ARCHAIC_AUTOSOMES if scope == "autosomes" else (scope,)):
                            writer.writerow(row)
        if scope == "autosomes":
            result = scan.summarize_exposure(background, files, out,
                reps=args.bootstrap, seed=args.seed, min_training_sites=args.min_training_sites)
        else:
            result = archaic_summarize_sex_exposure(background, files, out, scope, args.min_training_sites)
    if archaic_file_identities([r["path"] for r in args.source_identity["inputs"]]) != args.source_identity["inputs"]:
        raise ValueError("Archaic inputs changed during analysis")
    result.update(schema=ARCHAIC_SCHEMA, source_identity=args.source_identity, dataset=dataset["dataset"],
        n_samples=dataset["n_samples"], track=track_result, chromosome_audits=audits, chromosome_scope=scope,
        background_signature=args.background_signature, background_read_only=True,
        denominator="Original eligible SNV positions within lifted tracts; includes eligible single-ALT loci",
        interpretation="Unphased tract carrier prevalence, not archaic allele frequency or individual ancestry")
    return result


def archaic_summarize_sex_exposure(background, files, out, chrom, min_training_sites=1):
    """Descriptive X/Y O/E against matched loci outside each label on that chromosome.

    AC/AN come from the original scanner, never from an assumed diploid 2*N.
    Global counters cannot resolve PAR or sex strata, so no inferential P/CI is
    produced. Include outside_partition: autosomal LD coverage is not required.
    """
    dims = ("context", "af_bin", "an_bin", "sensitivity", "population_mode", "threshold")
    total = defaultdict(lambda: [0, 0])
    exposed = defaultdict(lambda: [0, 0])
    def read(paths):
        for path in paths:
            with gzip.open(path, "rt") as handle:
                yield from csv.DictReader(handle, delimiter="\t")
    def counts(row):
        n, m = int(row["n_sites"]), int(row["n_multi"])
        if not 0 <= m <= n:
            raise ValueError("Invalid eligible/multi-ALT counts")
        return n, m
    for row in read([background]):
        if row["threshold"] not in {"5", "10"} or row["population_mode"] != "all":
            continue
        if row["chrom"].removeprefix("chr") != chrom:
            raise ValueError("Sex-chromosome background contains another chromosome")
        n, m = counts(row)
        key = tuple(row[x] for x in dims)
        total[key][0] += n; total[key][1] += m
    if not total:
        raise ValueError(f"No original chr{chrom} global background counts")
    for row in read(files):
        if row["threshold"] not in {"5", "10"} or row["population_mode"] != "all":
            continue
        if row["chrom"].removeprefix("chr") != chrom:
            raise ValueError("Exposure chromosome differs from selected background")
        n, m = counts(row)
        key = (row["label"], tuple(row[x] for x in dims))
        exposed[key][0] += n; exposed[key][1] += m
    summaries = defaultdict(lambda: [0, 0, 0., 0])
    out = Path(out)
    with gzip.open(out/"sex_chromosome_strata.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chrom", "label", *dims, "exposed_sites", "observed_multi", "control_sites", "control_multi", "expected_multi"])
        for (label, st), (n, m) in sorted(exposed.items()):
            tn, tm = total[st][0] - n, total[st][1] - m
            if not 0 <= tm <= tn:
                raise ValueError(f"chr{chrom} exposure is not a subset of its original background: {label}")
            expected = n * (tm + .5)/(tn + 1) if tn >= min_training_sites else None
            writer.writerow([chrom, label, *st, n, m, tn, tm, expected if expected is not None else "NA"])
            z = summaries[label, st[-1], st[3]]
            z[0] += n; z[1] += m
            if expected is None: z[3] += n
            else: z[2] += expected
    fields = ['label','threshold','sensitivity','observed_multi','expected_multi','oe','ci_low','ci_high','bootstrap_p','fdr','block_count','n_eligible','n_without_loco','outside_partition_sites','status']
    with (out/"enrichment_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t"); writer.writeheader()
        for (label, threshold, sensitivity), (n, m, e, missing) in sorted(summaries.items()):
            writer.writerow(dict(label=label, threshold=threshold, sensitivity=sensitivity,
                observed_multi=m, expected_multi="NA" if missing else e,
                oe=m/e if e and not missing else "NA", ci_low="NA", ci_high="NA",
                bootstrap_p="NA", fdr="NA", block_count=0, n_eligible=n,
                n_without_loco="NA", outside_partition_sites="NA",
                status="insufficient_same_chromosome_control" if missing else "descriptive_same_chromosome"))
    manifest = dict(status="complete", n_tests=0, n_descriptive_rows=len(summaries),
        chromosome_scope=chrom, primary_thresholds=[5, 10],
        inference="Descriptive only; no P, CI or FDR; no autosomal LD-block assumption",
        background="Same chromosome outside each label, matched context/primary-AF/AN/sensitivity; (multi+0.5)/(sites+1)",
        empty_exposure=not bool(summaries),
        limitations=["Original AC/AN retained; no assumption AN=2*N",
            "PAR/non-PAR and sex composition are not resolved by existing counters",
            "Carrier denominator is the chromosome-specific GU analyzed-person panel, including zero-tract people; not allele frequency",
            "X/Y estimates must not be interpreted as autosomal LOCO tests"])
    (out/"enrichment_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def archaic_reusable(directory, expected):
    try:
        result = validate_manifest(directory, require_evidence=True)
        return result if result.get("status") == "complete" and result.get("source_identity") == expected else None
    except (OSError, ValueError, KeyError):
        return None


def archaic_combine_summaries(sources, out):
    """Recompute multiplicity over the selected cohorts, keeping cohort O/E."""
    import scan
    import scan
    data, fields = [], None
    for dataset, path in sources:
        with Path(path).open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if fields is not None and fields != reader.fieldnames:
                raise ValueError("Archaic cohort summary schemas differ")
            fields = reader.fieldnames
            for row in reader:
                row["dataset"] = dataset
                row["fdr"] = "NA"
                data.append(row)
    groups = {}
    for row in data:
        try:
            p = float(row["bootstrap_p"])
        except (ValueError, TypeError):
            continue
        if row["status"] == "complete" and math.isfinite(p):
            groups.setdefault((row["threshold"], row["sensitivity"]), []).append(row)
    for group in groups.values():
        for row, q in zip(group, scan.bh_fdr([float(r["bootstrap_p"]) for r in group])):
            row["fdr"] = q
    with (Path(out)/"archaic_enrichment.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, ["dataset", *(fields or [])], delimiter="\t")
        writer.writeheader(); writer.writerows(data)
    return data


def archaic_plot_results(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    selected = [r for r in rows if r["status"] == "complete" and r["threshold"] == "10"
                and r["sensitivity"] == "exclude_combined" and r["label"].endswith(":any")
                and r["ci_low"] not in {"", "NA", "None"}]
    if not selected:
        return
    fig, ax = plt.subplots(figsize=(10, max(4, .24*len(selected))))
    for i, row in enumerate(selected):
        lo, value, hi = (float(row[k]) for k in ("ci_low", "oe", "ci_high"))
        ax.plot([lo, hi], [i, i], color="#176B87"); ax.plot(value, i, "o", color="#176B87", ms=3)
    ax.set_yticks(range(len(selected)), [r["label"] for r in selected], fontsize=7)
    ax.axvline(1, color="black", lw=1, ls="--")
    ax.set_xlabel("Observed / LOCO expected; AC ≥ 10, strict masks; 95% LD-block bootstrap CI")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(Path(out)/f"archaic_enrichment.{ext}", dpi=160)
    plt.close(fig)


@transactional("manifest.json")
def archaic_run(args):
    selection = getattr(args, "chromosomes", "auto")
    records = archaic_discover(args.gu_root, args.datasets, args.module, selection, getattr(args, "references", ARCHAIC_DEFAULT_REFERENCES),
                       getattr(args, "source_build", "auto"))
    common = archaic_modern_inputs(args.modern_dir)
    # This validates the relocated source and all original QC inputs. New
    # exposure queries use the original options, not fresh scanner defaults.
    query_args = SimpleNamespace(out=args.out, modern_dir=args.modern_dir, fasta=None,
        erv_joined=None, segdup_bed=None, lcr_bed=None, low_mappability_bed=None, chromosomes=None)
    original, annotation, tasks, chromosomes = modern_targeted_inputs(query_args, allowed_chromosomes=ARCHAIC_CHROMOSOMES)
    options = original["options"]
    if not {5, 10} <= set(options.get("thresholds", [])):
        raise ValueError("Archaic analysis requires original AC5/10 counters")
    ready = []
    for record in records:
        if record["status"] != "ready":
            continue
        missing = set(record["chromosomes"]) - chromosomes
        if missing:
            reason = f"Original modern scan lacks chromosomes: {sorted(missing)}"
            if selection != "auto" or record["chromosome_scope"] == "autosomes":
                raise ValueError(reason)
            record.update(status="not_available", reason=reason)
        else:
            ready.append(record)
    if not ready:
        raise ValueError("No selected GU/modern chromosome pairs are ready")
    if any(r["chromosome_scope"] == "autosomes" for r in ready) and not options.get("ld_block_dir"):
        raise ValueError("Autosomal LOCO requires EUR LD blocks")
    summaries, manifests, snapshots = [], {}, []
    for dataset in ready:
        expected = archaic_source_identity(dataset, args, common)
        snapshots.append(expected['inputs'])
        directory = Path(args.cache_dir)/dataset["dataset"]
        result = archaic_reusable(directory, expected)
        if result is None:
            # cache_dir is this child's output container, not a scientific input.
            # Keep actual sources in the transaction's protected input set.
            child = SimpleNamespace(**{k: v for k, v in vars(args).items() if k != 'cache_dir'})
            child.out = str(directory); child.dataset = dataset; child.source_identity = expected
            child.tasks = [(c, path) for c, path in tasks if c in dataset["chromosomes"]]
            child.annotation = str(annotation); child.options = options
            child.background_signature = original["signature"]
            archaic_analyze_dataset(child)
            result = validate_manifest(directory, require_evidence=True)
            dataset["operation"] = "computed"
        else:
            print(f"[reuse] archaic {dataset['dataset']}: inputs and outputs verified", flush=True)
            dataset["operation"] = "reused"
        manifests[dataset["dataset"]] = str(directory/"manifest.json")
        summaries.append((dataset["dataset"], directory/"enrichment_summary.tsv"))
    for saved in snapshots:
        if archaic_file_identities([r['path'] for r in saved]) != saved:
            raise ValueError('Archaic inputs changed while preparing the combined result')
    rows = archaic_combine_summaries(summaries, args.out)
    with (Path(args.out)/"datasets.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, ["dataset", "chromosome_scope", "status", "operation", "source", "n_samples", "build", "reason"],
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(records)
    archaic_plot_results(rows, args.out)
    return dict(status="complete", schema=ARCHAIC_SCHEMA, datasets=records, dataset_manifests=manifests,
        n_tests=sum(r["status"] == "complete" for r in rows), n_summary_rows=len(rows), background_signature=original["signature"], background_read_only=True,
        fdr_scope="All selected archaic datasets within each AC threshold and sensitivity",
        method="Original scanner QC; autosomes: LOCO/EUR LD-block bootstrap; X/Y: separate same-chromosome outside-tract descriptive O/E",
        interpretation="Cohorts remain separate; overlapping people are not independent replication; zero-tract people remain in carrier denominator",
        estimand="gnomAD multi-ALT enrichment inside regions called as tracts in each GU cohort; not an individual-level genotype association",
        limitations=["A tract overlap does not establish that the third/fourth nucleotide was carried by the tract carrier or introduced by archaic gene flow",
            "The gnomAD mutation/count background is shared; it is not the GU cohort's ancestry-specific allele-frequency background",
            "For individual co-occurrence and exact archaic-base support use archaic-alleles with the same persons' genotypes",
            "Context, primary AF and AN matching does not remove all mappability, local ancestry or recombination confounding"])


def archaic_parser():
    p = argparse.ArgumentParser(description=_description_herv_archaic)
    p.add_argument("--gu-root", default=configured_path('ANALYSIS_ROOT', 'gu'))
    p.add_argument("--datasets", default="auto")
    p.add_argument("--source-build", choices=("auto", "37", "38"), default="auto",
                   help="Select the GU source assembly when both build directories exist; caller metadata must agree")
    p.add_argument("--references", default=",".join(ARCHAIC_DEFAULT_REFERENCES), help="Comma-separated tested references; historical .2013 references require explicit selection")
    p.add_argument("--module", choices=("auto", "ibdmix", "ibdmix-GPU"), default="auto")
    p.add_argument("--chromosomes", default="auto", help="auto: autosomes plus available X/Y; or autosomes,X,Y / X,Y / X / Y")
    p.add_argument("--modern-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cache-dir", required=True)
    p.add_argument("--chain", default=os.environ.get('ALLELE_CHAIN37TO38', configured_path('DIR_FILES', 'liftOver/hg19ToHg38.over.chain.gz')))
    p.add_argument("--liftover", default=configured_path('DIR_SOFTWARE', 'bin/liftOver'))
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--min-lod", type=float, default=4.)
    p.add_argument("--min-length", type=int, default=50000)
    p.add_argument("--carrier-thresholds", default="0.01,0.05,0.10,0.20")
    p.add_argument("--bootstrap", type=int, default=2000)
    p.add_argument("--seed", type=int, default=20261024)
    p.add_argument("--min-training-sites", type=int, default=1)
    return p


def archaic_main(argv=None):
    args = archaic_parser().parse_args(argv)
    if args.bootstrap < 2000 or args.min_length < 1 or args.min_lod < 0 or args.min_training_sites < 1:
        raise ValueError("Require ≥2000 bootstrap draws, positive length/training minimum and nonnegative LOD")
    thresholds = [float(x) for x in args.carrier_thresholds.split(",")]
    if not thresholds or any(not 0 < x <= 1 for x in thresholds) or len(set(thresholds)) != len(thresholds):
        raise ValueError("Carrier thresholds must be distinct values in (0,1]")
    protected = [Path(args.gu_root).resolve(), Path(args.modern_dir).resolve()]
    outputs = [Path(args.out).resolve(), Path(args.cache_dir).resolve()]
    for out in outputs:
        for source in protected:
            if out == source or out in source.parents or source in out.parents:
                raise ValueError(f"Generated output overlaps permanent input: {out}")
    if outputs[0] == outputs[1] or outputs[0] in outputs[1].parents or outputs[1] in outputs[0].parents:
        raise ValueError("Summary output and dataset cache must be separate, non-nested directories")
    # Different report directories may share this persistent cohort cache.
    # Keep a stable view while validating, computing and combining its files.
    outputs[1].mkdir(parents=True, exist_ok=True)
    with (outputs[1]/'.archaic.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f'Archaic dataset cache is already in use: {outputs[1]}') from exc
        archaic_run(args)
    return 0


def _entry_herv_archaic():
    try:
        raise SystemExit(archaic_main())
    except (OSError, ValueError, RuntimeError, KeyError) as archaic_exc:
        raise SystemExit(f"archaic: {archaic_exc}")



# ========================================================================
# herv ld
# ========================================================================
_description_herv_ld = 'Compare HERV sequence and LD partitions on the same physical windows.\n\nThis module measures an association between an annotation and existing LD maps.\nIt does not estimate recombination from genotypes or identify a causal driver.\nUse maps made with comparable samples/SNPs/algorithms for biological contrasts.\n\nExample:\n  python f/herv_ld.py --ld-block-dir /mnt/f/LD/block     --erv-joined /mnt/f/annot/allele3/repeats/hg38/rmskJoinedCurrent.txt.gz     --out /path/to/results/ld --chrom-sizes /path/to/hg38.fa.fai\n\nOptional covariates: a tab-delimited table with chrom,start,end (BED0), exactly\nmatching the physical grid, and gc,recombination,mappability,snp_density.\nMissing rows/values remain missing; they are never silently replaced with zero.'
#!/usr/bin/env python3
'Compare HERV sequence and LD partitions on the same physical windows.\n\nThis module measures an association between an annotation and existing LD maps.\nIt does not estimate recombination from genotypes or identify a causal driver.\nUse maps made with comparable samples/SNPs/algorithms for biological contrasts.\n\nExample:\n  python f/herv.py --ld-block-dir /mnt/f/LD/block     --erv-joined /mnt/f/annot/allele3/repeats/hg38/rmskJoinedCurrent.txt.gz     --out /path/to/results/ld --chrom-sizes /path/to/hg38.fa.fai\n\nOptional covariates: a tab-delimited table with chrom,start,end (BED0), exactly\nmatching the physical grid, and gc,recombination,mappability,snp_density.\nMissing rows/values remain missing; they are never silently replaced with zero.\n'
ld_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(ld_paths_root / "0f"))
from bisect import bisect_right

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scan import parse_joined_rmsk, read_blocks  # Existing BED0/joined parser.

LD_VERSION = "herv_ld_1.0"
LD_DEFAULT_ANCESTRIES = ("AFR", "EAS", "EUR", "SAS")
LD_DEFAULT_FAMILIES = ("ERV1", "ERVK", "ERVL", "ERVL-MaLR")
LD_REQUIRED_COVARIATES = ("gc", "recombination", "mappability", "snp_density")
LD_METRICS = ("boundary_density_per_mb", "log2_bp_weighted_median_block_bp")


def ld_merge_intervals(intervals):
    """Union BED0 half-open intervals, including touching intervals."""
    result = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError(f"Invalid interval {start}-{end}")
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def ld_intersect_intervals(left, right):
    """Intersection of two sorted, nonoverlapping interval lists."""
    result, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        if start < end:
            result.append((start, end))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return result


class ld_IntervalIndex:
    def __init__(self, intervals):
        self.intervals = ld_merge_intervals(intervals)
        self.starts = [v[0] for v in self.intervals]
        self.ends = [v[1] for v in self.intervals]

    def clip(self, start, end):
        lo = bisect_right(self.ends, start)
        hi = bisect_left(self.starts, end)
        return [(max(a, start), min(b, end)) for a, b in self.intervals[lo:hi]]

    def bp(self, start, end):
        return sum(b-a for a, b in self.clip(start, end))

    def spans_boundary(self, position):
        # Both bases flanking an LD junction must belong to common coverage.
        i = bisect_right(self.starts, position)-1
        return i >= 0 and self.starts[i] < position < self.ends[i]


def ld_weighted_median(values, weights):
    values, weights = np.asarray(values), np.asarray(weights)
    if not len(values) or not np.any(weights > 0):
        return float("nan")
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[order])
    return float(values[order[np.searchsorted(cumulative, cumulative[-1]/2)]])


def ld_read_sizes(path):
    result = {}
    with open(path) as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"):
                continue
            row = line.split()
            chrom = canonical_chrom(row[0])
            if chrom not in {str(i) for i in range(1, 23)}:
                continue
            length = int(row[1])
            if length <= 0 or chrom in result:
                raise ValueError(f"Invalid/duplicate chromosome size: {line.strip()}")
            result[chrom] = length
    if not result:
        raise ValueError("Chromosome sizes contain no autosomes")
    return result


def ld_build_profiles(block_sets, elements, *, sizes=None, window_bp=1_000_000,
                   min_common_coverage=.95, min_window_fraction=.5,
                   families=LD_DEFAULT_FAMILIES):
    """One row per physical window and ancestry; use identical common bases.

    A boundary is a shared end/start of two adjacent blocks, not an edge next to
    an unobserved gap. Block lengths refer to original blocks; their weights are
    the number of common-covered bases in this window, not clipped lengths.
    """
    if window_bp < 1 or not 0 < min_common_coverage <= 1 or not 0 < min_window_fraction <= 1:
        raise ValueError("Invalid physical-window settings")
    autosomes = {str(i) for i in range(1, 23)}
    ancestry = list(block_sets)
    by_anc, max_ends = {}, defaultdict(int)
    for anc, blocks in block_sets.items():
        by_chrom = defaultdict(list)
        for block in blocks:
            if block.chrom not in autosomes:
                continue
            by_chrom[block.chrom].append(block)
            max_ends[block.chrom] = max(max_ends[block.chrom], block.end)
        by_anc[anc] = by_chrom
    lengths = dict(sizes or max_ends)
    for chrom, end in max_ends.items():
        if chrom not in lengths or end > lengths[chrom]:
            raise ValueError(f"LD intervals exceed chromosome sizes on chr{chrom}")
    if not lengths:
        raise ValueError("No autosomal LD intervals")

    repeats = defaultdict(list)
    for element in elements:
        if element.rep_class == "LTR" and element.family in families and element.chrom in lengths:
            for start, end in element.segments:
                if end > lengths[element.chrom] and sizes is not None:
                    raise ValueError(f"Repeat annotation exceeds chromosome size on chr{element.chrom}")
                repeats[element.chrom, element.family].append((start, end))
    repeat_index = {key: ld_IntervalIndex(value) for key, value in repeats.items()}
    empty = ld_IntervalIndex([])
    rows, windows = [], []
    for chrom in sorted(lengths, key=int):
        coverage, block_index, junctions = {}, {}, {}
        for anc in ancestry:
            blocks = sorted(by_anc[anc].get(chrom, []), key=lambda b: b.start)
            coverage[anc] = ld_IntervalIndex((b.start, b.end) for b in blocks)
            block_index[anc] = (blocks, [b.start for b in blocks], [b.end for b in blocks])
            junctions[anc] = [right.start for left, right in zip(blocks, blocks[1:])
                              if left.end == right.start]
        common = coverage[ancestry[0]].intervals
        for anc in ancestry[1:]:
            common = ld_intersect_intervals(common, coverage[anc].intervals)
        common_index = ld_IntervalIndex(common)
        for start in range(0, lengths[chrom], window_bp):
            end = min(start+window_bp, lengths[chrom])
            shared = common_index.clip(start, end)
            common_bp = sum(b-a for a, b in shared)
            eligible = (common_bp/(end-start) >= min_common_coverage
                        and (end-start)/window_bp >= min_window_fraction)
            base = dict(chrom=chrom, start=start, end=end, window_bp=end-start,
                        common_bp=common_bp, common_coverage_fraction=common_bp/(end-start),
                        eligible=int(eligible))
            for family in families:
                index = repeat_index.get((chrom, family), empty)
                bp = sum(index.bp(a, b) for a, b in shared)
                base[family+"_aligned_bp"] = bp
                base[family+"_fraction"] = bp/common_bp if common_bp else float("nan")
            windows.append(base)
            for anc in ancestry:
                blocks, starts, ends = block_index[anc]
                overlap_weights, original_lengths = [], []
                # A block may cross several common-coverage fragments: count it
                # once and sum its overlap before taking either median.
                i, j = bisect_right(ends, start), bisect_left(starts, end)
                for block in blocks[i:j]:
                    overlap = sum(max(0, min(b, block.end)-max(a, block.start)) for a, b in shared)
                    if overlap:
                        original_lengths.append(block.end-block.start)
                        overlap_weights.append(overlap)
                bounds = junctions[anc]
                lo, hi = bisect_left(bounds, start), bisect_left(bounds, end)
                n_boundary = sum(common_index.spans_boundary(p) for p in bounds[lo:hi])
                median = ld_weighted_median(original_lengths, overlap_weights)
                rows.append(base | dict(
                    ancestry=anc, covered_bp=coverage[anc].bp(start, end),
                    coverage_fraction=coverage[anc].bp(start, end)/(end-start),
                    n_blocks_intersecting_common=len(original_lengths),
                    median_intersecting_block_bp=float(np.median(original_lengths)) if original_lengths else float("nan"),
                    bp_weighted_median_block_bp=median,
                    log2_bp_weighted_median_block_bp=math.log2(median) if median > 0 else float("nan"),
                    internal_boundary_count=n_boundary,
                    boundary_density_per_mb=n_boundary*1_000_000/common_bp if common_bp else float("nan")))
    return rows, windows, lengths


def ld_load_covariates(path, windows, columns=LD_REQUIRED_COVARIATES):
    """Join covariates to exact, canonical BED0 physical-window keys."""
    if not path:
        return {}, [], {"provided": False, "matched_windows": 0}
    opener = gzip.open if str(path).endswith(".gz") else open
    target = {(w["chrom"], w["start"], w["end"]) for w in windows}
    values, ignored = {}, 0
    with opener(path, "rt") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"chrom", "start", "end"}.issubset(reader.fieldnames or []):
            raise ValueError("Covariate TSV requires chrom,start,end (BED0 physical grid)")
        present = [name for name in columns if name in reader.fieldnames]
        if not present:
            raise ValueError("Covariate TSV has none of the requested covariate columns")
        seen = set()
        for row in reader:
            key = (canonical_chrom(row["chrom"]), int(row["start"]), int(row["end"]))
            if key in seen:
                raise ValueError(f"Duplicate covariate window: {key}")
            seen.add(key)
            if key not in target:
                ignored += 1
                continue
            result = {}
            for name in present:
                raw = row[name].strip()
                result[name] = float("nan") if raw in {"", ".", "NA", "NaN", "nan"} else float(raw)
                if math.isinf(result[name]):
                    raise ValueError(f"Infinite covariate {name} at {key}")
            values[key] = result
    if not values:
        raise ValueError("No covariate rows match the physical window grid; check build/window size")
    return values, present, dict(provided=True, matched_windows=len(values),
                                missing_windows=len(target)-len(values), ignored_rows=ignored)


def ld_center_by_chromosome(matrix, weights, chromosomes):
    centered = np.asarray(matrix, float).copy()
    for chrom in np.unique(chromosomes):
        use = chromosomes == chrom
        centered[use] -= np.average(centered[use], weights=weights[use], axis=0)
    return centered


def ld_fit_association(x, y, chromosomes, weights, covariates=None, covariate_names=(), *,
                    bootstrap=1000, seed=20261027, min_chromosomes=8, min_windows=30):
    """Weighted within-chromosome slope, resampling whole chromosomes.

    x is HERV percentage (0--100), so beta is per percentage point. Chromosome
    fixed effects are removed by weighted centering. Paired population-difference
    outcomes use this same routine and the same physical windows.
    """
    x, y, weights = map(lambda v: np.asarray(v, float), (x, y, weights))
    chromosomes = np.asarray(chromosomes, str)
    covariates = np.empty((len(x), 0)) if covariates is None else np.asarray(covariates, float)
    if covariates.ndim != 2 or len(covariates) != len(x):
        raise ValueError("Covariates must have one row per window")
    use = (np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0)
           & np.isfinite(covariates).all(axis=1))
    x, y, chromosomes, weights, covariates = [v[use] for v in (x, y, chromosomes, weights, covariates)]
    groups = np.unique(chromosomes)
    result = dict(n_windows=len(x), n_chromosomes=len(groups), beta_per_percentage_point=None,
                  ci_low=None, ci_high=None, bootstrap_p=None, fdr=None,
                  bootstrap_successful=0, test_status="insufficient_windows_or_chromosomes",
                  covariates_used="", covariates_absorbed_or_redundant="")
    if len(x) < min_windows or len(groups) < min_chromosomes:
        return result
    within = ld_center_by_chromosome(np.column_stack([x, y, covariates]), weights, chromosomes)
    x, yc, cc = within[:, 0], within[:, 1], within[:, 2:]
    if np.sum(weights*yc*yc) <= 1e-20:
        result.update(beta_per_percentage_point=None, test_status="no_within_chromosome_outcome_variation")
        return result
    if np.sum(weights*x*x) <= 1e-15:
        result["test_status"] = "no_within_chromosome_exposure_variation"
        return result
    kept, dropped = [], []
    cdesign = np.empty((len(x), 0))
    # Choose identifiable covariates without considering the exposure. A
    # confounder collinear with HERV must not be dropped to rescue its estimate.
    for i, name in enumerate(covariate_names):
        scale = math.sqrt(np.average(cc[:, i]**2, weights=weights))
        if scale <= 1e-10:
            dropped.append(name)
            continue
        column = cc[:, i]/scale
        trial = np.column_stack([cdesign, column])
        if np.linalg.matrix_rank(trial*np.sqrt(weights)[:, None]) > cdesign.shape[1]:
            kept.append(name)
            cdesign = trial
        else:
            dropped.append(name)
    result["covariates_used"] = ",".join(kept)
    result["covariates_absorbed_or_redundant"] = ",".join(dropped)
    design = np.column_stack([x, cdesign])
    weighted = design*np.sqrt(weights)[:, None]
    if np.linalg.matrix_rank(weighted) != design.shape[1]:
        result["test_status"] = "exposure_not_identifiable_after_adjustment"
        return result
    xx, xy = [], []
    for chrom in groups:
        use = chromosomes == chrom
        xx.append(design[use].T @ (design[use]*weights[use, None]))
        xy.append(design[use].T @ (yc[use]*weights[use]))
    xx, xy = np.asarray(xx), np.asarray(xy)
    beta = float(np.linalg.solve(xx.sum(axis=0), xy.sum(axis=0))[0])
    result.update(beta_per_percentage_point=beta, test_status="estimated_no_bootstrap")
    if bootstrap == 0:
        return result
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(bootstrap):
        counts = rng.multinomial(len(groups), np.full(len(groups), 1/len(groups)))
        gram = np.tensordot(counts, xx, axes=1)
        if np.linalg.matrix_rank(gram) != len(gram):
            continue
        estimate = np.linalg.solve(gram, np.tensordot(counts, xy, axes=1))[0]
        if np.isfinite(estimate):
            estimates.append(estimate)
    result["bootstrap_successful"] = len(estimates)
    if len(estimates) < max(20, .8*bootstrap):
        result["test_status"] = "bootstrap_unstable"
        return result
    estimates = np.asarray(estimates)
    lower, upper = np.quantile(estimates, [.025, .975])
    # Centered two-sided bootstrap test; descriptive association, not a causal
    # or randomized-test P value. BH correction spans all reported slope tests.
    pvalue = (1+np.sum(np.abs(estimates-beta) >= abs(beta)))/(len(estimates)+1)
    result.update(ci_low=float(lower), ci_high=float(upper), bootstrap_p=float(pvalue), test_status="complete")
    return result


def ld_mean_difference(y, chromosomes, weights, bootstrap, seed):
    groups = np.unique(chromosomes)
    sums = np.array([np.sum(y[chromosomes == c]*weights[chromosomes == c]) for c in groups])
    denominators = np.array([np.sum(weights[chromosomes == c]) for c in groups])
    result = dict(mean_difference=float(sums.sum()/denominators.sum()), mean_ci_low=None, mean_ci_high=None)
    if bootstrap and len(groups) >= 3:
        counts = np.random.default_rng(seed).multinomial(len(groups), np.full(len(groups), 1/len(groups)), bootstrap)
        values = counts @ sums/(counts @ denominators)
        result.update(mean_ci_low=float(np.quantile(values, .025)), mean_ci_high=float(np.quantile(values, .975)))
    return result


def ld_analyze(profiles, families, ancestries, covariate_values, covariate_names, *,
            bootstrap=1000, seed=20261027, min_chromosomes=8, min_windows=30):
    index = defaultdict(dict)
    for row in profiles:
        if row["eligible"]:
            index[row["chrom"], row["start"], row["end"]][row["ancestry"]] = row
    keys = sorted(index, key=lambda k: (int(k[0]), k[1]))
    # Every slope and every paired contrast shares a complete-case window set.
    keys = [key for key in keys if set(index[key]) == set(ancestries)
            and all(np.isfinite(index[key][anc][metric]) for anc in ancestries for metric in LD_METRICS)
            and all(np.isfinite(covariate_values.get(key, {}).get(name, np.nan)) for name in covariate_names)]
    if set(LD_REQUIRED_COVARIATES).issubset(covariate_names):
        level = "adjusted_observational_association"
    elif covariate_names:
        level = "descriptive_partial_adjustment"
    else:
        level = "descriptive_without_sequence_recombination_covariates"
    base = [index[key][ancestries[0]] for key in keys]
    chromosomes = np.asarray([key[0] for key in keys])
    weights = np.asarray([row["common_bp"]/1_000_000 for row in base])
    covariates = np.array([[covariate_values[key][name] for name in covariate_names] for key in keys]) if covariate_names else np.empty((len(keys), 0))
    if covariate_names and not keys:
        covariates = np.empty((0, len(covariate_names)))
    association, contrasts = [], []
    for family in families:
        x = np.asarray([row[family+"_fraction"]*100 for row in base])
        for metric in LD_METRICS:
            outcomes = {anc: np.asarray([index[key][anc][metric] for key in keys]) for anc in ancestries}
            kwargs = dict(bootstrap=bootstrap, seed=seed, min_chromosomes=min_chromosomes, min_windows=min_windows)
            for anc in ancestries:
                fit = ld_fit_association(x, outcomes[anc], chromosomes, weights, covariates, covariate_names, **kwargs)
                association.append(dict(family=family, ancestry=anc, metric=metric, evidence_level=level) | fit)
            for first, second in itertools.combinations(ancestries, 2):
                delta = outcomes[first]-outcomes[second]
                fit = ld_fit_association(x, delta, chromosomes, weights, covariates, covariate_names, **kwargs)
                means = ld_mean_difference(delta, chromosomes, weights, bootstrap, seed) if keys else dict(mean_difference=None, mean_ci_low=None, mean_ci_high=None)
                contrasts.append(dict(family=family, ancestry_a=first, ancestry_b=second, metric=metric,
                                      contrast=f"{first}-minus-{second}", evidence_level=level) | means | fit)
    tests = [row for row in association+contrasts if row["bootstrap_p"] is not None]
    ordered = sorted(tests, key=lambda row: row["bootstrap_p"])
    previous = 1.
    for i in range(len(ordered)-1, -1, -1):
        previous = min(previous, ordered[i]["bootstrap_p"]*len(ordered)/(i+1))
        ordered[i]["fdr"] = previous
    return association, contrasts, dict(eligible_physical_windows=len(index), complete_case_windows=len(keys),
                                       dropped_incomplete_windows=len(index)-len(keys), evidence_level=level,
                                       slope_tests=len(association)+len(contrasts), tested_slopes=len(tests),
                                       fdr_scope="all families, both outcomes, all ancestry-specific and paired-difference slope tests")


def ld_write_tsv(path, rows):
    if not rows:
        raise ValueError(f"No output rows for {path}")
    path = Path(path)
    temporary = path.with_name(path.name+".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(temporary, "wt", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if value is None or isinstance(value, (float, np.floating)) and not np.isfinite(value)
                             else value for key, value in row.items()})
    temporary.replace(path)


def ld_plot_results(profiles, association, contrasts, families, ancestries, sizes, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 180, "pdf.fonttype": 42})
    eligible = [row for row in profiles if row["eligible"]]
    if not eligible:
        return []
    offsets, total = {}, 0
    for chrom in sorted(sizes, key=int):
        offsets[chrom] = total
        total += sizes[chrom]
    palette = ("#9C4A75", "#177D8D", "#D18E24", "#547749", "#6966A6", "#626262")
    colors = {anc: palette[i % len(palette)] for i, anc in enumerate(ancestries)}
    files = []
    def save(fig, name):
        for extension in ("png", "pdf"):
            file = out/(name+"."+extension)
            fig.savefig(file, bbox_inches="tight")
            files.append(file.name)
        plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True)
    first = [r for r in eligible if r["ancestry"] == ancestries[0]]
    x = [(offsets[r["chrom"]]+(r["start"]+r["end"])/2)/1e6 for r in first]
    for family in families:
        axes[0].plot(x, [r[family+"_fraction"]*100 for r in first], lw=.7, alpha=.8, label=family)
    axes[0].set_ylabel("Aligned HERV sequence (%)\nof common-covered bases")
    axes[0].legend(ncol=len(families), frameon=False, loc="upper right")
    for anc in ancestries:
        group = [r for r in eligible if r["ancestry"] == anc]
        gx = [(offsets[r["chrom"]]+(r["start"]+r["end"])/2)/1e6 for r in group]
        axes[1].scatter(gx, [r["boundary_density_per_mb"] for r in group], s=3, alpha=.6, label=anc, color=colors[anc])
        axes[2].scatter(gx, [r["bp_weighted_median_block_bp"]/1e3 for r in group], s=3, alpha=.6, label=anc, color=colors[anc])
    axes[1].set_ylabel("Internal LD boundaries\nper common-covered Mb")
    axes[1].legend(ncol=len(ancestries), frameon=False)
    axes[2].set_ylabel("Base-weighted median\nLD block length (kb)")
    axes[2].set_yscale("log")
    ticks = [(offsets[c]+sizes[c]/2)/1e6 for c in sorted(sizes, key=int)]
    axes[2].set_xticks(ticks, ["chr"+c for c in sorted(sizes, key=int)], rotation=45)
    axes[2].set_xlabel("Identical physical windows and common LD-map coverage")
    for ax in axes:
        for offset in offsets.values():
            ax.axvline(offset/1e6, color="#CCCCCC", lw=.4, zorder=0)
    fig.suptitle("HERV sequence and LD-map architecture | descriptive association")
    fig.tight_layout()
    save(fig, "ld_landscape")

    fig, axes = plt.subplots(1, 2, figsize=(12, max(4, len(families)*1.2)), sharey=True)
    for ax, metric in zip(axes, LD_METRICS):
        for ai, anc in enumerate(ancestries):
            group = {r["family"]: r for r in association if r["ancestry"] == anc and r["metric"] == metric}
            for fi, family in enumerate(families):
                row = group[family]
                beta = row["beta_per_percentage_point"]
                if beta is None:
                    continue
                yy = fi+(ai-(len(ancestries)-1)/2)*.13
                ax.plot(beta, yy, "o", ms=4, color=colors[anc], label=anc if fi == 0 else None)
                if row["ci_low"] is not None:
                    ax.plot([row["ci_low"], row["ci_high"]], [yy, yy], color=colors[anc], lw=1.2)
        ax.axvline(0, color="#555555", lw=.8)
        ax.set_xlabel("LD outcome change per +1 percentage point HERV\n(chromosome fixed effects; chromosome bootstrap 95% CI)")
        ax.set_title("Boundary density per Mb" if metric == LD_METRICS[0] else "log2 block length")
        ax.set_yticks(range(len(families)), families)
        ax.grid(axis="x", color="#EEEEEE")
    axes[0].legend(frameon=False, ncol=2)
    fig.suptitle("Within-chromosome association; no causal interpretation")
    fig.tight_layout()
    save(fig, "ld_association")

    pairs = [f"{a}-minus-{b}" for a, b in itertools.combinations(ancestries, 2)]
    fig, axes = plt.subplots(1, 2, figsize=(12, max(4, len(pairs)*.5)))
    for ax, metric in zip(axes, LD_METRICS):
        group = {(r["contrast"], r["family"]): r for r in contrasts if r["metric"] == metric}
        matrix = np.array([[group[pair, family]["beta_per_percentage_point"]
                            if group[pair, family]["beta_per_percentage_point"] is not None else np.nan
                            for family in families] for pair in pairs], dtype=float)
        finite = np.isfinite(matrix)
        maximum = max(float(np.max(np.abs(matrix[finite]))) if finite.any() else 0, 1e-8)
        im = ax.imshow(np.ma.masked_invalid(matrix), vmin=-maximum, vmax=maximum, cmap="RdBu_r", aspect="auto")
        ax.set_xticks(range(len(families)), families, rotation=30)
        ax.set_yticks(range(len(pairs)), [p.replace("-minus-", " − ") for p in pairs])
        ax.set_title("Boundary density contrast" if metric == LD_METRICS[0] else "log2 block-length contrast")
        for i, pair in enumerate(pairs):
            for j, family in enumerate(families):
                row = group[pair, family]
                if row["fdr"] is not None and row["fdr"] < .05:
                    ax.text(j, i, "•", ha="center", va="center", color="black")
        fig.colorbar(im, ax=ax, fraction=.04, label="Paired slope difference per +1 pp HERV")
    fig.suptitle("Population contrasts at the same physical windows | dot: BH FDR < 0.05")
    fig.tight_layout()
    save(fig, "ld_ancestry_contrasts")
    return files


def ld_file_info(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def ld_parser():
    p = argparse.ArgumentParser(description=_description_herv_ld, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ld-block-dir", default=configured_path('DIR_LD', 'block'))
    p.add_argument("--erv-joined", default=configured_path('DIR_ANNOT', 'allele3/repeats/hg38/rmskJoinedCurrent.txt.gz'))
    p.add_argument("--out", required=True)
    p.add_argument("--build", choices=("37", "38"), default="38")
    p.add_argument("--annotation-build", choices=("37", "38"), default="38", help="Declared build of joined annotation; must match --build")
    p.add_argument("--ancestries", default=",".join(LD_DEFAULT_ANCESTRIES))
    p.add_argument("--families", default=",".join(LD_DEFAULT_FAMILIES))
    p.add_argument("--chrom-sizes", help="Two-column chromosome sizes or FASTA .fai; inferred from LD envelope if omitted")
    p.add_argument("--window-bp", type=int, default=1_000_000)
    p.add_argument("--min-common-coverage", type=float, default=.95)
    p.add_argument("--min-window-fraction", type=float, default=.5)
    p.add_argument("--covariates", help="BED0 exact-grid TSV: chrom,start,end plus numerical covariates")
    p.add_argument("--covariate-columns", default=",".join(LD_REQUIRED_COVARIATES))
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261027)
    p.add_argument("--min-chromosomes", type=int, default=8)
    p.add_argument("--min-windows", type=int, default=30)
    p.add_argument("--no-plots", action="store_true")
    return p


@transactional("manifest.json")
def ld_run(args):
    if args.build != args.annotation_build:
        raise ValueError("LD BED and RepeatMasker builds must match; this module does not liftover")
    if args.bootstrap < 0 or args.min_chromosomes < 3 or args.min_windows < 3:
        raise ValueError("Require bootstrap >= 0, min-chromosomes >= 3, min-windows >= 3")
    ancestries = tuple(v.strip() for v in args.ancestries.split(",") if v.strip())
    families = tuple(v.strip() for v in args.families.split(",") if v.strip())
    if len(ancestries) < 2 or len(set(ancestries)) != len(ancestries) or not families or len(set(families)) != len(families):
        raise ValueError("Require >=2 unique ancestries and unique nonempty families")
    paths = {anc: Path(args.ld_block_dir)/f"{anc}.{args.build}.bed" for anc in ancestries}
    sets, block_stats = {}, {}
    for anc, path in paths.items():
        sets[anc], block_stats[anc] = read_blocks(path, anc, args.build)
        if not sets[anc]:
            raise ValueError(f"Empty LD partition: {path}")
    elements = parse_joined_rmsk(args.erv_joined)
    sizes = ld_read_sizes(args.chrom_sizes) if args.chrom_sizes else None
    profiles, windows, lengths = ld_build_profiles(sets, elements, sizes=sizes, window_bp=args.window_bp,
        min_common_coverage=args.min_common_coverage, min_window_fraction=args.min_window_fraction, families=families)
    columns = tuple(v.strip() for v in args.covariate_columns.split(",") if v.strip())
    covariate_values, covariate_names, covariate_qc = ld_load_covariates(args.covariates, windows, columns)
    for row in profiles:
        values = covariate_values.get((row["chrom"], row["start"], row["end"]), {})
        row.update({"cov_"+name: values.get(name, float("nan")) for name in covariate_names})
    association, contrasts, analysis_qc = ld_analyze(profiles, families, ancestries, covariate_values, covariate_names,
        bootstrap=args.bootstrap, seed=args.seed, min_chromosomes=args.min_chromosomes, min_windows=args.min_windows)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ld_write_tsv(out/"window_profiles.tsv.gz", profiles)
    ld_write_tsv(out/"association.tsv", association)
    ld_write_tsv(out/"ancestry_contrasts.tsv", contrasts)
    figures = [] if args.no_plots else ld_plot_results(profiles, association, contrasts, families, ancestries, lengths, out)
    manifest = dict(version=LD_VERSION, status="complete", build=args.build, annotation_build=args.annotation_build,
        inputs={"ld_partitions": {anc: ld_file_info(path) for anc, path in paths.items()},
                "erv_joined": ld_file_info(args.erv_joined),
                "chrom_sizes": ld_file_info(args.chrom_sizes) if args.chrom_sizes else None,
                "covariates": ld_file_info(args.covariates) if args.covariates else None},
        settings=vars(args), block_counts=block_stats,
        coordinate_convention="BED0 half-open; autosomes only",
        genome_extent="chromosome sizes" if sizes else "largest LD interval end per chromosome; not a callable genome mask",
        common_coverage="intersection of interval unions of ALL requested ancestries, identically used for every exposure and outcome",
        boundary_definition="junctions of touching blocks with both flanking bases in common coverage; gap edges are excluded",
        repeat_definition="union of aligned joined RepeatMasker segments per family; consensus tails and nested insertion gaps excluded",
        block_length_definition="length of original block, weighted by common bases overlapped in the physical window",
        model="WLS within chromosome (chromosome fixed effects); per +1 percentage point HERV aligned fraction; complete cases shared across all tests",
        uncertainty="whole-chromosome bootstrap, percentile 95% CI; centered two-sided bootstrap P; global BH for slopes",
        covariates=covariate_qc | {"included": covariate_names, "missing_key_covariates": sorted(set(LD_REQUIRED_COVARIATES)-set(covariate_names))},
        analysis=analysis_qc,
        interpretation="Association between repeat annotation and supplied LD partitions; no genotypic LD recalculation, recombination causality, or viral driving-force inference",
        map_comparability="Sample size, SNP ascertainment, LD method and map construction must be comparable; this module cannot establish that from BED intervals",
        tables=["window_profiles.tsv.gz", "association.tsv", "ancestry_contrasts.tsv"], figures=figures)
    temporary = out/"manifest.json.tmp"
    temporary.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
    temporary.replace(out/"manifest.json")
    return manifest


def ld_main(argv=None):
    args = ld_parser().parse_args(argv)
    try:
        result = ld_run(args)
    except (ValueError, OSError) as error:
        raise SystemExit(f"herv_ld: {error}") from error
    print(json.dumps({"status": result["status"], "out": str(Path(args.out).resolve()), **result["analysis"]}))


def _entry_herv_ld():
    ld_main()


# ========================================================================
# herv cohort
# ========================================================================
_description_herv_cohort = 'Plan a balanced HERV census subset and audit completed individual censuses.\n\nThe existing ``census.py cohort burden`` remains the scalable all-sample\ngenotype analysis. This module never runs RepeatMasker or changes a completed\ncensus. Percentages describe reference-conditioned sequence models, not\nindependent assemblies or the prevalence of polymorphic HERV insertions.'
#!/usr/bin/env python3
"""Plan a balanced HERV census subset and audit completed individual censuses.

The existing ``census.py cohort burden`` remains the scalable all-sample
genotype analysis. This module never runs RepeatMasker or changes a completed
census. Percentages describe reference-conditioned sequence models, not
independent assemblies or the prevalence of polymorphic HERV insertions.
"""
import html
import shlex
import statistics

COHORT_VERSION = "1.1.0"
COHORT_CATEGORIES = ("ERV_without_MaLR", "ERV_with_MaLR", "all_LTR")
COHORT_COUNT_FIELDS = ("union_bp", "genome_bp", "union_non_N_bp", "genome_non_N_bp")
COHORT_PERCENT_FIELDS = ("percent_genome", "percent_non_N")
COHORT_LABELS = ("diploid", "reference_callable_diploid")


def cohort_sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cohort_write_tsv(path, rows, fields):
    with Path(path).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if r.get(k) is None else r.get(k, "") for k in fields} for r in rows)


def cohort_save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def cohort_metadata(path, group_column="auto", sex_column="auto"):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        lines = [line for line in f if line.strip() and not line.startswith("##")]
    if not lines:
        raise ValueError("Empty sample metadata")
    values = list(csv.reader(lines, delimiter="\t")) if "\t" in lines[0] else [line.split() for line in lines]
    header = [h.lower().lstrip("#").strip() for h in values[0]]
    if len(set(header)) != len(header):
        raise ValueError("Duplicate metadata header")
    def column(request, candidates):
        names = candidates if request == "auto" else (request.lower(),)
        return next((name for name in names if name in header), None)
    sample_column = column("auto", ("sample_id", "sample", "iid", "id"))
    group_column = column(group_column, ("super_pop", "superpop", "super_population", "population", "pop"))
    sex_column = column(sex_column, ("sex", "gender"))
    if not sample_column or not group_column or not sex_column:
        raise ValueError("Metadata requires sample/sample_id/IID, a population column, and sex/gender; auto prefers super_pop")
    if any(len(row) != len(header) for row in values[1:]):
        raise ValueError("Metadata row does not match its header")
    rows = [dict(zip(header, row)) for row in values[1:]]
    seen, result = set(), []
    for row in rows:
        sample = row[sample_column].strip()
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            raise ValueError(f"Invalid sample_id: {sample!r}")
        if sample in seen:
            raise ValueError(f"Duplicate sample_id: {sample}")
        seen.add(sample)
        group, sex = row[group_column].strip(), row[sex_column].strip()
        if not group or not sex or group in {"NA", "."} or sex in {"NA", "."}:
            raise ValueError(f"Missing population/sex metadata for {sample}")
        include = row.get("include", "1").strip()
        unrelated = row.get("unrelated", "unknown").strip()
        if include not in {"0", "1"} or unrelated not in {"0", "1", "unknown"}:
            raise ValueError("include must be 0/1; unrelated must be 0/1/unknown")
        result.append(dict(sample_id=sample, population=group, sex=sex,
                           include=include, unrelated=unrelated))
    if not result:
        raise ValueError("Empty sample metadata")
    return result


def cohort_balanced_subset(rows, total=20, per_stratum=2, seed=20261102, always=()):
    if total < 1 or per_stratum < 1:
        raise ValueError("total and per-stratum must be positive")
    eligible = {r["sample_id"]: r for r in rows if r["include"] == "1" and r["unrelated"] != "0"}
    if len(set(always)) != len(always) or set(always) - eligible.keys():
        raise ValueError("Explicit included samples must be unique and eligible in metadata")
    if len(always) > total:
        raise ValueError("Explicit included samples exceed total")
    # Sorting first makes the result independent of input row order.
    grouped = defaultdict(list)
    for sample in sorted(eligible):
        row = eligible[sample]
        grouped[row["population"], row["sex"]].append(row)
    rng = random.Random(seed)
    for group in sorted(grouped):
        rng.shuffle(grouped[group])
    chosen = [eligible[s] for s in always]
    selected = set(always)
    counts = defaultdict(int)
    for row in chosen:
        counts[row["population"], row["sex"]] += 1
    # Round robin gives each population/sex stratum a turn before a second pick.
    for _ in range(per_stratum):
        for key in sorted(grouped):
            if len(chosen) >= total:
                break
            candidates = [r for r in grouped[key] if r["sample_id"] not in selected]
            if candidates and counts[key] < per_stratum:
                row = candidates[0]
                chosen.append(row)
                selected.add(row["sample_id"])
                counts[key] += 1
    if not chosen:
        raise ValueError("No eligible samples")
    support = [dict(population=k[0], sex=k[1], eligible=len(v), selected=counts[k])
               for k, v in sorted(grouped.items())]
    return chosen, support


def cohort_validate_fraction(row):
    parsed = dict(row)
    for field in COHORT_COUNT_FIELDS:
        raw = row[field]
        if isinstance(raw, bool) or isinstance(raw, float) and (not math.isfinite(raw) or not raw.is_integer()):
            raise ValueError(f"Noninteger {field}")
        val = int(raw)
        if val < 0:
            raise ValueError(f"Negative {field}")
        parsed[field] = val
    if (parsed["union_bp"] > parsed["genome_bp"] or
            parsed["union_non_N_bp"] > parsed["genome_non_N_bp"] or
            parsed["genome_non_N_bp"] > parsed["genome_bp"] or
            parsed["union_non_N_bp"] > parsed["union_bp"]):
        raise ValueError("HERV union or non-N counts exceed their denominator")
    for field, numerator, denominator in (
            ("percent_genome", "union_bp", "genome_bp"),
            ("percent_non_N", "union_non_N_bp", "genome_non_N_bp")):
        expected = 100 * parsed[numerator] / parsed[denominator] if parsed[denominator] else None
        raw = row.get(field)
        actual = None if raw in (None, "", "NA") else float(raw)
        if expected is None:
            if actual is not None:
                raise ValueError("A zero denominator requires an unavailable percentage")
        elif actual is None or not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(f"{field} disagrees with its counts")
        parsed[field] = actual
    return parsed


def cohort_inspect_census(root, sample):
    root = Path(root)
    run = json.loads((root / "run.json").read_text())
    done = json.loads((root / "analysis.done.json").read_text())
    result = json.loads((root / "genome_fraction.json").read_text())
    settings = run.get("settings", {})
    if run.get("status") != "complete":
        raise ValueError("Census run is not complete")
    if settings.get("sample") != sample:
        raise ValueError("Census sample differs from metadata")
    if done.get("settings") != settings or result.get("settings") != settings:
        raise ValueError("Summary/completion settings differ from the census run")
    required = {"method", "fasta_sha256", "library_sha256", "repeatmasker_sha256",
                "rmblast_version", "chunk_bp", "halo", "uncertain_flank", "parameters", "code_sha256"}
    if required - settings.keys():
        raise ValueError("Census lacks comparison provenance: " + ", ".join(sorted(required - settings.keys())))
    scope = [str(c).removeprefix("chr") for c in run.get("scope", [])]
    if len(scope) != len(set(scope)) or not scope or any(c not in {str(i) for i in range(1, 23)} for c in scope):
        raise ValueError("This census implementation requires an explicit unique autosomal scope")
    table = root / "genome_fraction.tsv"
    if done.get("results_sha256") != cohort_sha256(table):
        raise ValueError("Census table SHA256 differs from the completion receipt")
    with table.open() as f:
        rows = [cohort_validate_fraction(r) for r in csv.DictReader(f, delimiter="\t")]
    indexed = {(r["label"], r["category"]): r for r in rows}
    expected_rows = [cohort_validate_fraction(r) for r in result.get("categories", [])]
    expected = {(r["label"], r["category"]): r for r in expected_rows}
    if len(indexed) != len(rows) or len(expected) != len(expected_rows) or set(indexed) != set(expected):
        raise ValueError("Duplicate or inconsistent fraction categories")
    for key, row in indexed.items():
        if any(row[f] != expected[key][f] for f in COHORT_COUNT_FIELDS + COHORT_PERCENT_FIELDS):
            raise ValueError("Census JSON and TSV disagree")
    if any((label, cat) not in indexed for label in COHORT_LABELS for cat in COHORT_CATEGORIES):
        raise ValueError("Missing diploid/matched-reference categories")
    for label in COHORT_LABELS:
        categories = [indexed[label, cat] for cat in COHORT_CATEGORIES]
        if len({(r['genome_bp'], r['genome_non_N_bp']) for r in categories}) != 1:
            raise ValueError("HERV categories do not share the same sequence denominator")
        for field in ('union_bp', 'union_non_N_bp'):
            counts = [r[field] for r in categories]
            if counts != sorted(counts):
                raise ValueError("HERV categories must nest: ERV without MaLR <= ERV with MaLR <= all LTR")
    # Code versions remain separate unless proven comparable by a later audit.
    comparison = {k: settings[k] for k in sorted(settings) if k not in {"sample", "library", "repeatmasker"}}
    comparison["chromosomes"] = sorted(scope, key=int)
    comparison_id = hashlib.sha256(json.dumps(comparison, sort_keys=True).encode()).hexdigest()
    output = []
    for cat in COHORT_CATEGORIES:
        model, ref = (indexed[label, cat] for label in COHORT_LABELS)
        row = dict(sample_id=sample, category=cat, comparison_id=comparison_id,
                   scope="autosomes:" + ",".join(comparison["chromosomes"]),
                   model_genome_bp=model["genome_bp"], model_non_N_bp=model["genome_non_N_bp"],
                   callable_model_fraction=model["genome_non_N_bp"] / model["genome_bp"] if model["genome_bp"] else None,
                   model_percent=model["percent_genome"], model_percent_non_N=model["percent_non_N"],
                   matched_reference_percent=ref["percent_genome"], matched_reference_percent_non_N=ref["percent_non_N"],
                   delta_percentage_points=(model["percent_non_N"] - ref["percent_non_N"]
                                            if model["percent_non_N"] is not None and ref["percent_non_N"] is not None else None),
                   source=str(root), results_sha256=done["results_sha256"])
        output.append(row)
    return output, comparison


def cohort_census_inputs(raw):
    """Track the saved summaries actually consulted, without rehashing genomes."""
    if raw.get('metadata') and raw.get('census_root'):
        for row in cohort_metadata(raw['metadata'], raw.get('group_column', 'auto'), raw.get('sex_column', 'auto')):
            if row['include'] != '1':
                continue
            root = Path(raw['census_root']) / row['sample_id']
            for name in ('run.json', 'analysis.done.json', 'genome_fraction.json', 'genome_fraction.tsv'):
                path = root / name
                if path.is_file():
                    yield InputPath(path)


@transactional('plan.json', input_files=('metadata',), input_dirs=('census_root',), resolve_inputs=cohort_census_inputs)
def cohort_plan(args):
    roster = cohort_metadata(args.metadata, args.group_column, args.sex_column)
    chosen, support = cohort_balanced_subset(roster, args.total, args.per_stratum, args.seed,
                                     tuple(s for s in args.include_samples.split(",") if s))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    planned, commands = [], ["#!/usr/bin/env bash", "set -euo pipefail", "# Sequential, resumable individual HERV validation; this file is a plan only."]
    for row in chosen:
        target = Path(args.census_root).resolve() / row["sample_id"]
        action, reason = "run_or_resume", ""
        if (target / "analysis.done.json").exists():
            try:
                cohort_inspect_census(target, row["sample_id"])
                action = "reuse_complete"
            except (ValueError, OSError, KeyError, TypeError) as exc:
                action, reason = "blocked_invalid_completed_result", str(exc)
        command = ""
        if action == "run_or_resume":
            argv = ["bash", str(Path(args.entry).resolve()), "--dat", args.dataset,
                "--build", args.build, "--sample", row["sample_id"], "--out", str(target),
                "--replace", "FALSE", "--no-reuse", "--jobs", str(args.jobs), "--cores", str(args.cores)]
            for name, flag in (("pfile_dir", "--pfile-dir"), ("vcf_dir", "--vcf-dir"),
                               ("fasta", "--fasta"), ("library", "--library"), ("python", "--python")):
                value = getattr(args, name, None)
                if value:
                    resolved = value if name == 'python' and '/' not in value else str(Path(value).expanduser().resolve())
                    argv.extend([flag, resolved])
            command = shlex.join(argv)
            commands.append(command)
        planned.append(dict(row, run_dir=str(target), action=action, reason=reason, command=command))
    cohort_write_tsv(out / "sample_plan.tsv", planned, ["sample_id", "population", "sex", "unrelated", "run_dir", "action", "reason", "command"])
    cohort_write_tsv(out / "selected_metadata.tsv", chosen,
              ["sample_id", "population", "sex", "include", "unrelated"])
    (out / "selected_samples.txt").write_text("\n".join(r["sample_id"] for r in chosen) + "\n")
    cohort_write_tsv(out / "stratum_support.tsv", support, ["population", "sex", "eligible", "selected"])
    (out / "run_selected.sh").write_text("\n".join(commands) + "\n")
    return dict(status="planned", version=COHORT_VERSION, out=str(out), metadata_sha256=cohort_sha256(args.metadata),
        selected=len(planned), seed=args.seed, requested_maximum=args.total, per_stratum=args.per_stratum,
        selected_metadata=str(out / "selected_metadata.tsv"), sample_ids=str(out / "selected_samples.txt"),
        census_root=str(Path(args.census_root).resolve()),
        censuses_completed=sum(r["action"] == "reuse_complete" for r in planned),
        censuses_blocked=sum(r["action"] == "blocked_invalid_completed_result" for r in planned),
        note="No analyses launched. Zero-selected strata are disclosed; this small validation subset is not a population prevalence sample.",
        all_sample_route="census.py cohort burden: shared reference atlas and cohort PGEN; no per-person RepeatMasker required.")


def cohort_census_report(out, rows, audit):
    """Report saved measurements without treating an absent census as a zero."""
    verified = sum(r['status'] == 'verified_summary' for r in audit)
    requested = sum(r['status'] != 'excluded_metadata' for r in audit)
    fields = ('sample_id', 'population', 'sex', 'category', 'model_percent_non_N',
              'matched_reference_percent_non_N', 'delta_percentage_points', 'callable_model_fraction')
    def table(records, columns):
        def value(v):
            return f'{v:.6g}' if isinstance(v, float) else str(v if v is not None else 'NA')
        header = ''.join('<th>' + html.escape(k) + '</th>' for k in columns)
        body = ''.join('<tr>' + ''.join('<td>' + html.escape(value(r.get(k))) + '</td>' for k in columns) + '</tr>' for r in records)
        return '<div class="table"><table><thead><tr>' + header + '</tr></thead><tbody>' + body + '</tbody></table></div>'
    text = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>HERV census: measured sequence models</title>
<style>body{font-family:system-ui,sans-serif;margin:3em;max-width:1400px;color:#18334a}h1,h2{color:#195e72}.table{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border-bottom:1px solid #dbe4ea;padding:8px;text-align:left}th{background:#edf4f7}p{line-height:1.6}.note{background:#edf4f7;padding:1em}</style>
<h1>HERV census: measured sequence models</h1>'''
    text += f'<p><strong>{verified} / {requested} included samples have verified completed summaries.</strong> Missing or invalid censuses remain unavailable; they are not assigned zero HERV content.</p>'
    text += '''<p class="note">These are autosomal VCF/PGEN-conditioned reference sequence models. The percentage is RepeatMasker-aligned bases divided by model sequence bases; the non-N denominator excludes gaps and masked uncertainty. ERV_without_MaLR, ERV_with_MaLR and all_LTR are separate, nested categories. No generic “8%” value is inserted. The matched reference uses the same uncertainty masks; delta is measured in percentage points.</p>
<h2>Individual measurements</h2>'''
    text += table(rows, fields) if rows else '<p>No verified measurements are available.</p>'
    text += '''<h2>Interpretation and next analysis</h2><p>Agreement among a small ancestry/sex-balanced subset validates reproducibility of this reference-conditioned calculation. It does not estimate the prevalence of polymorphic HERV insertions or prove that viral insertion created a third nucleotide allele. Whole-element insertion/deletion presence requires a separate SV or mobile-element genotype assay. For all-sample sequence variation, use the shared HERV atlas plus the existing cohort genotype-burden analysis, retaining missingness and ancestry covariates.</p>
<p>Comparisons across different FASTA, repeat-library, software, scope or calling settings are kept in separate comparison_id groups. A small validation subset is not a representative population sample.</p>
<h2>Availability audit</h2>'''
    text += table(audit, ('sample_id', 'status', 'reason'))
    text += '<p><a href="individual_census.tsv">Individual table</a> · <a href="group_summary.tsv">Stratified descriptive summaries</a> · <a href="sample_audit.tsv">Availability audit</a></p></html>'
    (Path(out) / 'report.html').write_text(text)


@transactional('manifest.json', input_files=('metadata',), input_dirs=('census_root',), resolve_inputs=cohort_census_inputs)
def cohort_summarize(args):
    roster = cohort_metadata(args.metadata, args.group_column, args.sex_column)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows, audit, comparisons = [], [], {}
    for meta in roster:
        root = Path(args.census_root) / meta["sample_id"]
        if meta["include"] != "1":
            audit.append(dict(sample_id=meta["sample_id"], status="excluded_metadata", reason="include=0"))
            continue
        try:
            found, comparison = cohort_inspect_census(root, meta["sample_id"])
            for row in found:
                row.update(population=meta["population"], sex=meta["sex"], unrelated=meta["unrelated"])
                rows.append(row)
                comparisons[row["comparison_id"]] = comparison
            audit.append(dict(sample_id=meta["sample_id"], status="verified_summary", reason=""))
        except (ValueError, OSError, KeyError, TypeError) as exc:
            audit.append(dict(sample_id=meta["sample_id"], status="unavailable", reason=str(exc)))
    fields = ["sample_id", "population", "sex", "unrelated", "category", "comparison_id", "scope",
              "model_genome_bp", "model_non_N_bp", "callable_model_fraction", "model_percent", "model_percent_non_N",
              "matched_reference_percent", "matched_reference_percent_non_N", "delta_percentage_points", "source", "results_sha256"]
    cohort_write_tsv(out / "individual_census.tsv", rows, fields)
    cohort_write_tsv(out / "sample_audit.tsv", audit, ["sample_id", "status", "reason"])
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["comparison_id"], row["population"], row["sex"], row["category"]].append(row)
    summary = []
    for key, rr in sorted(grouped.items()):
        for metric in ("model_percent_non_N", "matched_reference_percent_non_N", "delta_percentage_points", "callable_model_fraction"):
            vals = [r[metric] for r in rr if r[metric] is not None]
            summary.append(dict(comparison_id=key[0], population=key[1], sex=key[2], category=key[3], metric=metric,
                n=len(vals), mean=statistics.mean(vals) if vals else None,
                median=statistics.median(vals) if vals else None,
                sd=statistics.stdev(vals) if len(vals) > 1 else None,
                minimum=min(vals) if vals else None, maximum=max(vals) if vals else None,
                inference="descriptive_validation_subset"))
    cohort_write_tsv(out / "group_summary.tsv", summary, ["comparison_id", "population", "sex", "category", "metric", "n", "mean", "median", "sd", "minimum", "maximum", "inference"])
    unavailable = sum(r["status"] == "unavailable" for r in audit)
    state = "complete_with_exclusions" if rows and unavailable else "complete" if rows else "not_available"
    cohort_census_report(out, rows, audit)
    return dict(status=state, version=COHORT_VERSION, out=str(out), metadata_sha256=cohort_sha256(args.metadata),
        sample_coverage="complete" if not unavailable else "partial" if rows else "none",
        included_samples=sum(r["status"] != "excluded_metadata" for r in audit),
        verified_samples=sum(r["status"] == "verified_summary" for r in audit),
        unavailable_samples=sum(r["status"] == "unavailable" for r in audit), comparisons=comparisons,
        verification="Checks saved summary SHA256, count/percentage consistency, sample, completion and scientific settings; does not rehash raw genome/chunk files.",
        scope_note="The current individual census is autosomal. X/Y are retained in the separate allele3 chromosome analyses; no autosomal estimate is labeled whole-genome.",
        interpretation="Reference-conditioned sequence-model fractions and matched callable-reference contrasts. Not independent assemblies, insertion counts, or proof that HERV insertion created a third nucleotide allele.")


def cohort_main(argv=None):
    p = argparse.ArgumentParser(description=_description_herv_cohort)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("plan", "summarize"):
        q = sub.add_parser(name)
        q.add_argument("--metadata", required=True, help="Existing 1KG panel or TSV with sample/sample_id, population/super_pop, sex/gender")
        q.add_argument("--group-column", default="auto", help="Exact column; auto prefers super_pop then population/pop")
        q.add_argument("--sex-column", default="auto", help="Exact column; auto accepts sex or gender")
        q.add_argument("--census-root", required=True, help="Contains one census directory per sample ID")
        q.add_argument("--out", required=True)
        if name == "plan":
            q.add_argument("--entry", default=str(Path(__file__).resolve().parents[1] / "2.herv_census.sh"))
            q.add_argument("--dataset", default="1kg")
            q.add_argument("--build", choices=("37", "38"), default="38")
            q.add_argument("--total", type=int, default=20)
            q.add_argument("--per-stratum", type=int, default=2)
            q.add_argument("--include-samples", default="", help="Comma-separated eligible IDs counted within --total")
            q.add_argument("--seed", type=int, default=20261102)
            q.add_argument("--jobs", type=int, default=2)
            q.add_argument("--cores", type=int, default=8)
            for name in ("pfile-dir", "vcf-dir", "fasta", "library", "python"):
                q.add_argument("--" + name, help="Optional explicit path forwarded to every selected individual census")
    args = p.parse_args(argv)
    if args.command == "plan" and (not 1 <= args.jobs <= 8 or args.cores < 1):
        p.error("jobs must be 1..8 and cores must be positive")
    try:
        result = cohort_plan(args) if args.command == "plan" else cohort_summarize(args)
        print(json.dumps({key: result[key] for key in ('status', 'out', 'selected', 'verified_samples',
            'unavailable_samples', 'sample_coverage') if key in result}, indent=2))
    except (ValueError, OSError, KeyError, TypeError) as exc:
        print(f"herv-cohort: {exc}", file=sys.stderr)
        return 2
    return 0


def _entry_herv_cohort():
    raise SystemExit(cohort_main())


# ========================================================================
# herv report
# ========================================================================
_description_herv_report = 'Render measured HERV analyses with explicit estimands and availability.'
#!/usr/bin/env python3
"""Render measured HERV analyses with explicit estimands and availability."""


def report_table_html(data, columns, limit=30):
    data = list(data)
    if not data:
        return "<p class='muted'>No estimable rows available.</p>"
    head = "".join(f"<th>{html.escape(c.replace('_',' '))}</th>" for c in columns)
    body = "".join("<tr>"+"".join(f"<td>{html.escape(str(r.get(c,'')))}</td>" for c in columns)+"</tr>" for r in data[:limit])
    return f"<div class='table'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def report_load_json(path):
    return json.loads(Path(path).read_text()) if Path(path).is_file() else {}


def report_recorded_output(base, name, result):
    """Only files from a validated current manifest may be shown."""
    return name in result.get('output_inventory', {}) and (base / name).is_file()


def report_make_report(out, modern, analysis_root, pca_temporal_dir=None, *, include_ibdmix=False):
    module_out, modern, root = Path(out), Path(modern), Path(analysis_root)
    out = module_out / "runtime"
    out.mkdir(parents=True, exist_ok=True)
    sections = ["<h1>ERV1 / ERVK and multi-ALT variation</h1><p class='lead'>Measured results, analysis scope and completion status.</p>"]
    if not include_ibdmix:
        sections.append("<p>IBDmix evidence is disabled for this run. Archaic tract results are excluded from the main analysis.</p>")
    sections.append("<h2>1. Primary family enrichment</h2>")
    primary_dir = out / "primary"
    primary_available = False
    pipeline_state = report_load_json(module_out / "qc/pipeline_status.json")
    primary_state = next((r for r in pipeline_state.get("stages", []) if isinstance(r, dict) and r.get("stage") == "primary"), {})
    try:
        if primary_state.get("status") in {"failed", "not_available", "blocked", "skipped"}:
            raise ValueError("Primary stage did not complete in this invocation")
        primary = validate_manifest(primary_dir, "primary_manifest.json", require_evidence=True)
        if bool(primary.get("ibdmix_enabled", True)) != include_ibdmix:
            raise ValueError("Primary result uses a different IBDmix evidence scope; recompute primary")
        for expected in primary.get("source_counter_inputs", []):
            if identity(expected["path"]) != expected:
                raise ValueError("Primary counter inputs changed; recompute the primary stage")
        selected = primary_select_primary_rows(modern_rows(primary_dir / "enrichment_summary.tsv"))
        primary_available = True
        sections.append(report_table_html(selected, ["label", "threshold", "observed_multi", "expected_multi", "oe", "ci_low", "ci_high", "fdr", "status"]))
        if report_recorded_output(primary_dir, "primary_enrichment.png", primary):
            sections.append("<img src='primary/primary_enrichment.png' alt='Primary ERV family enrichment'>")
        sections.append("<p><a href='primary/enrichment_summary.tsv'>Complete primary result table</a> · <a href='primary/primary_manifest.json'>Analysis definition and input provenance</a></p>")
        details = {k: primary.get(k) for k in ("estimand", "target_scope", "background_scope", "minimum_training_sites", "bootstrap_replicates", "seed", "missing_rule")}
        sections.append("<pre>" + html.escape(json.dumps(details, ensure_ascii=False, indent=2)) + "</pre>")
    except (OSError, ValueError, KeyError) as exc:
        sections.append("<p class='muted'>Primary estimate unavailable: " + html.escape(str(exc)) + ". Other completed analyses remain reportable.</p>")
    sections.append("<p>AC5/AC10 requires two distinct nonreference ALT bases to meet the specified count threshold. Reference-base observation is established separately by joint genotypes. This SNV-denominator statistic is not a mutation rate per callable genomic base.</p>")
    status = report_load_json(module_out/"qc/pipeline_status.json")
    stages = status.get("stages", [])
    if isinstance(stages, dict):
        stages = [dict(stage=k, **v) for k,v in stages.items()]
    stage_status = {r.get("stage"): r for r in stages}
    if stages:
        sections += ["<h2>Execution status</h2>", report_table_html(stages, ["stage","status","result_status","reason"], 30)]
    evidence = [
        ("archaic", "archaic_enrichment.tsv", "Archaic tract enrichment in GU cohorts",
         "Finalized IBDmix cohorts use each chromosome scope's analyzed sample panel including zero-tract people. This estimates gnomAD multi-ALT enrichment inside GU tract regions; it is not a same-person genotype association or evidence that the extra nucleotide arose in an archaic genome. The shared gnomAD background is not the GU cohort's ancestry-specific allele-frequency background. Autosomes retain LOCO and EUR LD-block bootstrap, with BH FDR across selected cohorts. X/Y are separate descriptive O/E comparisons against matched loci outside each label on the same chromosome, including sites without LD blocks. Existing counters do not resolve PAR/non-PAR or sex strata: X/Y have no P, CI or FDR and are not significance tests. Cohorts may overlap; carrier prevalence is not archaic allele frequency. Use archaic-alleles for within-person target-base/tract co-occurrence and exact archaic-base support.", "archaic_enrichment.png"),
        ("spatial", "primary_spatial_partitions.tsv", "2. Spatial distribution within annotated ERVs",
         "The primary estimator is decomposed into strand-aware aligned-base head, body and tail partitions, with short elements and ambiguous overlaps retained. The end length and normalized-bin count are recorded in the spatial manifest. Partitions use expected-count weights and the same target domain as their total. External flanks are separate. Joined repeat records do not necessarily represent complete proviruses or TSS.", "primary_spatial_partitions.png"),
        ("summary", "common_support_enrichment.tsv", "5. Common-support sensitivity analysis",
         "Unsupported background strata are removed from both observed and expected counts. This sensitivity estimand has a different support rule from the primary full-scope estimator.", "family_enrichment.png"),
        ("windows", "frozen_peak_validation.tsv", "6. Family-specific observed-boundary profiles",
         "Strand-aware windows use the width, step and flank recorded in the windows manifest. The central window-to-distal-flank ratio tests positional concentration. Peak location is selected on discovery chromosomes, then evaluated on held-out chromosomes. These alignment edges are not transcription start sites.", "boundary_profiles.png"),
        ("genotypes", "cohort_retention.tsv", "7. Multi-person genotype concordance",
         "Actual A/C/G/T counts establish whether three observed bases coexist in sampled populations. Results are conditional on the modern candidate set and are not a new genome-wide enrichment estimate. Sample-independent replication requires verified non-overlap with gnomAD.", "genotype_concordance.png"),
        ("ancient", "fixed_panel_family_trajectories.tsv", "3. Full ancient cohort: modern-defined alleles over the past 10,000 years",
         "The full gnomAD HERV catalog, with no default 500-site cap, is followed across every original PSAM sample. Primary, secondary and additional ALT bases are audited separately. Primary time models use officially unrelated ancient individuals aged 0–10,000 BP; modern EUR samples remain a separate comparison group. Exact represented nucleotides are counted by their native PGEN allele indices. A reference base cannot be inferred as two minus one ALT when split records have unknown other-ALT encoding. Unmeasurable bases remain unavailable. Shared PCs and available region, data-source and coverage covariates support descriptive two-sided temporal associations, not selection coefficients or proof of historical tri-allelic occupancy.", "family_allele_time.png"),
        ("selection", "published_selection_summary.tsv", "4. Published exact-allele selection evidence",
         "Published Akbari 2026 statistics are attached to their exact original REF/ALT contrasts. The modern catalog is fixed before this join, and no coefficient is transferred to an unrepresented third base. Raw S, SE, calibrated X, P_X, posterior and FDR are preserved. Signed effects are not reoriented under the default raw_only contract. Counts are descriptive base identities, not independent loci or a new enrichment test; the ancient data overlap the source study.", None),
        ("ld", "association.tsv", "8. HERV and population LD structure",
         "Common physical windows and the intersection of the four ancestry maps supply a shared denominator. Chromosome bootstrap preserves paired regional comparisons. Association is not evidence that HERVs drive LD, especially without sequence and recombination covariates.", "ld_association.png"),
    ]
    evidence_order = {name:i for i,name in enumerate(("spatial", "ancient", "selection", "archaic", "summary", "windows", "genotypes", "ld"))}
    evidence.sort(key=lambda item:evidence_order[item[0]])
    for directory, filename, title, description, figure in evidence:
        if directory == "archaic" and not include_ibdmix:
            continue
        stage_name = "time" if directory == "ancient" else directory
        current = stage_status.get(stage_name, {})
        manifest_name = {"summary":"summary_manifest.json", "spatial":"spatial_manifest.json", "windows":"windows_manifest.json",
                         "genotypes":"genotypes_manifest.json", "ancient":"manifest.json", "selection":"manifest.json", "ld":"manifest.json", "archaic":"manifest.json"}[directory]
        # Read the canonical stage output only.
        candidates = [out/"time"] if directory == "ancient" else [out/directory]
        recorded = Path(current["result_manifest"]) if current.get("result_manifest") else None
        base = next((p for p in candidates if recorded is not None and (p/manifest_name).resolve() == recorded.resolve()), None)
        if base is None:
            base = next((p for p in candidates if (p/manifest_name).is_file() or (p/filename).is_file()), candidates[0])
        sections += [f"<h2>{title}</h2><p>{description}</p>"]
        path = base/filename
        if current.get("status") in {"failed", "not_available", "error", "blocked", "skipped"}:
            sections.append("<p class='muted'>This stage did not complete in the current run: "
                            + html.escape(str(current.get("status"))) + ". "
                            + html.escape(str(current.get("reason", "See execution status and stage log."))) + "</p>")
            if path.is_file():
                sections.append("<p class='muted'>Previous output exists but is not displayed as current evidence. "
                                "Its inputs or parameters may differ from this run.</p>")
            continue
        try:
            result = validate_manifest(base, manifest_name, require_evidence=True)
            if directory == "summary" and bool(result.get("ibdmix_enabled", True)) != include_ibdmix:
                raise ValueError("Summary uses a different IBDmix evidence scope; recompute summary")
        except (OSError, ValueError) as exc:
            sections.append("<p class='muted'>Current evidence unavailable: " + html.escape(str(exc)) + "</p>")
            continue
        complete_states = {"complete", "complete_descriptive_imputed_panel", "complete_descriptive_all_allele_panel", "complete_adjustment_unavailable"}
        # A direct module rerun can invalidate its manifest after an earlier
        # successful pipeline run. Either record of incompletion blocks display.
        result_state = next((s for s in (result.get("status"), current.get("result_status")) if s and s not in complete_states), None)
        if recorded is not None and not (base/manifest_name).is_file():
            result_state = "missing_result_manifest"
        if result_state:
            sections.append("<p class='muted'>The stage audit did not produce a usable effect estimate: "
                            + html.escape(str(result_state)) + ". Previous measurement files, if present, are not displayed.</p>")
            continue
        if report_recorded_output(base, filename, result):
            if current.get("status") not in {"complete", "reused"}:
                sections.append("<p class='muted'>Previously saved output. This stage was not validated in this invocation; "
                                "the table below is not a new measurement from the current run.</p>")
            data = list(modern_rows(path))
            if directory == "ancient":
                sections.append("<p>Genetic adjustment: " + html.escape(str(result.get('genetic_adjustment', 'unverified'))) + ". Raw family curves are not PC-adjusted.</p>")
                sections.append("<p>Actual full-catalog execution counts: <code>" + html.escape(json.dumps(
                    {k:result.get(k) for k in ("n_PSAM_samples", "sample_age_scopes", "catalog_sites", "exact_HQ_ALT_contrasts", "n_fixed_contrasts", "n_tested_primary_strict_AC5", "maximum_genotype_chunk_rows") if k in result}, ensure_ascii=False)) + "</code>. See the complete sample and allele-availability audits below.</p>")
            if directory == "summary":
                data = [r for r in data if r.get("label", "").startswith("family:") and r.get("sensitivity") == "exclude_combined"]
            if directory == 'archaic':
                data = [r for r in data if r.get('threshold') == '10' and r.get('sensitivity') == 'exclude_combined' and r.get('label', '').endswith(':any')]
                if report_recorded_output(base, 'datasets.tsv', result):
                    sections.append(f"<p><a href='{str((base/'datasets.tsv').relative_to(out))}'>Cohort inputs, sample counts, and unavailable cohorts</a></p>")
            if directory == "spatial":
                data = [r for r in data if r.get("scope") == "eur_ld_partition" and r.get("sensitivity") in {"exclude_combined", "strict"}]
            if directory == "ld":
                data = data[:20]
            relative = path.relative_to(out)
            sections.append(f"<p><a href='{html.escape(str(relative))}'>{html.escape(filename)}</a></p>")
            if data:
                columns = (["family","threshold","partition","observed_multi","expected_multi","oe","ci_low","ci_high","expected_weight","excess_contribution","status"]
                           if directory == "spatial" else ["dataset", "label", "threshold", "sensitivity", "n_eligible", "observed_multi", "expected_multi", "oe", "status"]
                           if directory == "archaic" else list(data[0])[:10])
                sections.append(report_table_html(data, columns, 30 if directory == "spatial" else 20))
            else:
                sections.append("<p class='muted'>The output contains no estimable rows; this is not a negative biological result.</p>")
            if figure and report_recorded_output(base, figure, result):
                sections.append(f"<img src='{html.escape(str((base/figure).relative_to(out)))}' alt='{html.escape(title)}'>")
            if directory == "selection":
                for name in ("published_selection_alleles.tsv.gz", "availability.tsv", "manifest.json"):
                    if report_recorded_output(base, name, result) or name == manifest_name and (base/name).is_file():
                        sections.append(f"<p><a href='{str((base/name).relative_to(out))}'>{name}</a></p>")
            if directory == "ancient":
                sections.append('<h3>Same-person HERV matched-panel statistics</h3><p>Incomplete personal means are observed-pair conditional. Formal personal models require the complete frozen panel; period comparisons retain fixed pair and composition weights. Person and chromosome uncertainty have separate scopes. These quantities are not ancient gnomAD O/E.</p>')
                for name in ("allele_time_heatmap.png", "allele_panel_availability.png", "F6_joint_measurability.png"):
                    if report_recorded_output(base, name, result):
                        sections.append(f"<img src='{str((base/name).relative_to(out))}' alt='Full ancient allele trajectories and representation'>")
                for name in ("allele_trajectories.tsv.gz", "measurement_plans.json", "adjusted_allele_trends.tsv.gz", "all_allele_audit.tsv.gz", "sample_time_region_coverage.tsv", "catalog_coverage.tsv", "matched_locus_role_trajectories.tsv", "model_covariates.json", "joint_measurability_funnel.tsv", "joint_person_statistics.tsv.gz", "joint_fixed_pair_panel.tsv", "joint_fixed_panel_trajectories.tsv", "joint_standardized_trajectories.tsv", "joint_person_models.tsv", "joint_missingness_models.tsv", "standardized_family_trajectories.tsv", "regional_fixed_panel_family_trajectories.tsv", "joint_manifest.json"):
                    if report_recorded_output(base, name, result):
                        sections.append(f"<p><a href='{str((base/name).relative_to(out))}'>{name}</a></p>")
            if directory == "spatial":
                if report_recorded_output(base, "primary_whole_element_profiles.png", result):
                    sections.append(f"<img src='{str((base/'primary_whole_element_profiles.png').relative_to(out))}' alt='Primary expected-count definition along aligned repeat bases'>")
                for name in ("primary_spatial_identity.tsv", "primary_counter_reconciliation.tsv", "family_counter_reconciliation.tsv", "spatial_annotation_audit.tsv"):
                    if report_recorded_output(base, name, result):
                        sections.append(f"<p><a href='{str((base/name).relative_to(out))}'>{name}</a></p>")
                sections.append("<p>The common-support sensitivity applies the recorded training minimum and removes unsupported observed and expected counts together.</p>")
                for name in ("spatial_partitions.tsv", "whole_element_profiles.tsv", "spatial_identity.tsv", "spatial_contrasts.tsv"):
                    if report_recorded_output(base, name, result):
                        sections.append(f"<p><a href='{str((base/name).relative_to(out))}'>{name}</a></p>")
        else:
            sections.append("<p class='muted'>Not run / inputs unavailable. No effect estimate is inferred.</p>")
    if pca_temporal_dir:
        sections.append('<h2>Genome-wide PCA temporal structure (separate estimand)</h2>')
        try:
            temporal = validate_manifest(pca_temporal_dir, require_evidence=True)
            sections.append('<p>Fixed shared coordinates; raw, fixed-composition and adjusted PC responses. This is not HERV selection evidence.</p>')
            sections.append('<pre>' + html.escape(json.dumps({k: temporal.get(k) for k in ('status', 'sample_counts', 'coordinate_contract', 'limitations')}, ensure_ascii=False, indent=2)) + '</pre>')
            for name in ('temporal_models.tsv', 'standardized_trajectories.tsv', 'report.html'):
                if report_recorded_output(Path(pca_temporal_dir), name, temporal):
                    link = os.path.relpath(Path(pca_temporal_dir) / name, out)
                    sections.append(f'<p><a href="{html.escape(link)}">PCA {name}</a></p>')
        except (OSError, ValueError) as exc:
            sections.append('<p>PCA temporal results unavailable: ' + html.escape(str(exc)) + '</p>')
    cohort = root/"final/herv_cohort/herv_cohort_burden.tsv"
    if not cohort.is_file() and os.environ.get("ALLELE_COHORT_TABLE"):
        cohort = Path(os.environ["ALLELE_COHORT_TABLE"])
    if cohort.is_file():
        n = len({r["sample_id"] for r in modern_rows(cohort)})
        sections.append(f"<h2>Measured cohort context</h2><p>The existing cohort table contains {n:,} people. Its measure is nonreference allele-copy burden within shared reference HERV-related intervals, including indels. It is not a count of HERV insertions or personal multi-allelic loci.</p>")
    sections.append("<h2>Interpretation</h2><p>Regional association and descriptive temporal trends do not establish mutation origin or selection. Independent genotype support, repeat-aware mapping, and an appropriate population design are required for those questions.</p>")
    css = "body{font:16px/1.6 system-ui,sans-serif;color:#20333c;max-width:1160px;margin:40px auto;padding:0 24px;background:#fbfcfc}h1{font-size:34px}h2{margin-top:36px;color:#176b87}.lead{font-size:20px}.muted{color:#687a83}img{max-width:100%;display:block;margin:24px auto;background:white}.table{overflow:auto}table{border-collapse:collapse;font-size:13px;width:100%;background:white}th,td{padding:9px;border-bottom:1px solid #dae4e8;text-align:left;white-space:nowrap}th{background:#eaf2f5}a{color:#146d97}code{font-size:13px}"
    text = "<!doctype html><html lang='en'><meta charset='utf-8'><title>HERV integrated evidence report</title><style>"+css+"</style><body>"+"\n".join(sections)+"</body></html>"
    from analysis_integrity import publish_linked_report
    report = publish_linked_report(module_out, out, text, 'herv_enrichment.html')
    return dict(status="complete", report=str(report), primary_available=primary_available)


def report_main(argv=None):
    p = argparse.ArgumentParser(description=_description_herv_report)
    p.add_argument("command", choices=["report"])
    p.add_argument("--out", required=True)
    p.add_argument("--modern-dir", required=True)
    p.add_argument("--analysis-root", required=True)
    p.add_argument("--pca-temporal-dir")
    p.add_argument("--ibdmix", type=boolean_option, default=False, metavar="TRUE|FALSE")
    a = p.parse_args(argv)
    result = report_make_report(a.out, a.modern_dir, a.analysis_root, a.pca_temporal_dir, include_ibdmix=a.ibdmix)
    print(json.dumps(result, indent=2))


def _entry_herv_report():
    report_main()


# ========================================================================
# herv pipeline
# ========================================================================
_description_herv_pipeline = 'Run current HERV analyses with explicit input and result availability.\n\nSource data are read-only. Missing data are reported without inventing estimates.'
#!/usr/bin/env python3
"""Run current HERV analyses with explicit input and result availability.

Source data are read-only. Missing data are reported without inventing estimates.
"""
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
import datetime as dt
import importlib.util
from analysis_integrity import identity as content_identity
from batch_workflows import STAGES, DEFAULT_STAGES, select_stages

VERSION = "6.0.1"
SOURCE = Path(__file__).resolve().parent
PREFIX = "AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes"
SCAN_STORE = str(Path(os.environ.get('ALLELE_ANALYSIS_ROOT', configured_path('ANALYSIS_ROOT', 'allele3'))) / 'gnomad_scan')


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def state(path):
    return content_identity(path)


def save(path, value):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    temporary = p.with_name(p.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n")
    os.replace(temporary, p)


def first_file(paths):
    return next((Path(p) for p in paths if Path(p).is_file()), None)


def resolve_output_dir(analysis_root, output=None):
    root = Path(analysis_root).expanduser().resolve()
    out = Path(output).expanduser().resolve() if output else root
    if out == root:
        out = root / "herv_enrichment"
    if out in root.parents:
        raise ValueError("Output cannot be a parent of the analysis root")
    return out


def resolve(args):
    root = Path(args.dir_analysis_root).expanduser().resolve()
    out = resolve_output_dir(root, args.dir_out)
    from analysis_integrity import resolve_result_input
    modern = Path(resolve_result_input(args.dir_modern or SCAN_STORE)).expanduser().resolve()
    if out == modern or out in modern.parents or modern in out.parents:
        raise ValueError("Output and source scan must use separate, non-nested directories")
    ancient = Path(args.dir_ancient).expanduser().resolve()
    prefix = Path(args.pgen_prefix).expanduser().resolve() if args.pgen_prefix else ancient / PREFIX
    pvar = first_file([str(prefix) + ".pvar", str(prefix) + ".pvar.zst"]) or Path(str(prefix) + ".pvar")
    old_metadata = first_file([root / "pca3d/tables/sample_metadata.tsv", root / "final/evolution/metadata/sample_metadata.tsv", root / "evolution/metadata/sample_metadata.tsv",
                               root / "pca3d/layered_1000yr/metadata/sample_metadata.tsv"])
    metadata_mode = "explicit" if args.metadata else "existing_metadata" if old_metadata else "prepare_from_original_supplement"
    metadata = Path(args.metadata).expanduser().resolve() if args.metadata else old_metadata or out / "runtime/metadata/sample_metadata.tsv"
    coordinate_candidates = []
    for directory in (root / "pca3d/coordinates", root / "pca3d/layered_1000yr/pca", root / "final/evolution/pca", root / "pca3d/basis", root / "evolution/pca"):
        coordinate_candidates += [directory / "all.sscore", directory / "all.sscore.gz", directory / "shared.eigenvec", directory / "pca.eigenvec"]
    coordinates = Path(args.coordinates).expanduser().resolve() if args.coordinates else first_file(coordinate_candidates)
    return dict(root=root, out=out, modern=modern, ancient=ancient, prefix=prefix, pvar=pvar,
        targets=Path(args.targets).expanduser().resolve() if args.targets else out / "runtime/candidates/targets.tsv",
        catalog=Path(args.allele_catalog).expanduser().resolve() if args.allele_catalog else out / "runtime/alleles/site_catalog.tsv.gz",
        selection=Path(args.selection_stats).expanduser().resolve() if args.selection_stats else ancient / "Selection_Summary_Statistics_01OCT2025.tsv.gz",
        metadata=metadata, metadata_mode=metadata_mode, coordinates=coordinates,
        hq=Path(args.hq_variants).expanduser().resolve() if args.hq_variants else ancient / "HighQualityVariants.tsv",
        supplement=Path(args.supplement_xlsx).expanduser().resolve() if args.supplement_xlsx else ancient / "41586_2026_10358_MOESM3_ESM.xlsx",
        chrom_sizes=Path(args.chrom_sizes).expanduser().resolve() if args.chrom_sizes else Path(args.fasta38 + ".fai"))


def stage_selection(args):
    return select_stages(args.command, args.stages, args.ibdmix)


def modern_files(modern, kind, mask_overrides=None):
    """Audit the files actually read, including merged-counter/shard fallback.

    Keep this metadata-only: preflight must report missing scientific packages
    rather than importing them before the dependency audit can run.
    """
    from scan_store import load_run_metadata, resolve_run_path
    modern = Path(modern)
    sidecars = [modern / name for name in ("relocation.json", "adoption_receipt.json", "input_locations.json") if (modern/name).is_file()]
    run_path = modern / "run.json"
    inputs = [run_path, *sidecars]
    if not run_path.is_file():
        return inputs
    run = load_run_metadata(modern)
    options = run.get("options", {})
    if mask_overrides:
        from scan_store import resolve_mask_options
        options = resolve_mask_options(run, mask_overrides)
    prefix = options.get("architecture_prefix", "architecture")
    targeted = kind in {"windows", "spatial"}
    counter_kinds = {"primary": ("global", "exposure"), "summary": ("global", "exposure"), "windows": ("global",),
                     "spatial": ("global", "erv") + (("exposure",) if options.get("ld_block_dir") else ())}.get(kind, ())
    # Match herv_modern.architecture_paths: each kind independently prefers its
    # merged file. An absent merged file is legitimate when complete shards exist.
    fallback = []
    for counter in counter_kinds:
        merged = modern / "summary" / f"architecture.{counter}_strata.tsv.gz"
        if merged.is_file():
            inputs.append(merged)
        else:
            fallback.append(counter)
    need_shards = targeted or kind == "candidates" or bool(fallback)
    for directory in run.get("shards", []):
        directory = Path(directory)
        if not directory.is_absolute():
            directory = modern / directory
        if need_shards:
            inputs += [directory / "COMPLETE.json", directory / "scan.json"]
        if kind == "candidates":
            inputs.append(directory / "multiallelic_sites.tsv.gz")
        inputs += [directory / f"{prefix}.{counter}_strata.tsv.gz" for counter in fallback]
        if (targeted or bool(fallback)) and (directory / "COMPLETE.json").is_file():
            # Existing scan.load_run verifies all recorded outputs; include these
            # in the cheap missing-input audit instead of reporting a false error.
            receipt = json.loads((directory / "COMPLETE.json").read_text())
            inputs += [directory / x for x in receipt.get("outputs", {})]
        if targeted and (directory / "scan.json").is_file():
            info = json.loads((directory / "scan.json").read_text())
            path = info.get("input")
            if isinstance(path, dict):
                path = path.get("path")
            if path:
                path = resolve_run_path(run, path)
                inputs.append(Path(path))
                index = first_file([str(path) + ".tbi", str(path) + ".csi"])
                inputs.append(index or Path(str(path) + ".tbi"))
    for key in ("architecture_segdup_bed", "architecture_lcr_bed", "architecture_low_mappability_bed"):
        if (kind == "candidates" or targeted) and options.get(key):
            inputs.append(Path(options[key]))
    if targeted:
        # These are verified/read by targeted_inputs and the regional scanners,
        # even though their old fingerprints are also embedded in run.json.
        for key in ("callable_bed", "fasta", "erv_joined"):
            if options.get(key): inputs.append(Path(options[key]))
        if options.get("fasta"):
            inputs.append(Path(options["fasta"] + ".fai"))
            if options["fasta"].endswith((".gz", ".bgz")):
                inputs.append(Path(options["fasta"] + ".gzi"))
        if options.get("ld_block_dir"):
            inputs.append(Path(options["ld_block_dir"]) / f"EUR.{options.get('architecture_build', '38')}.bed")
        identity = run.get("identity", {})
        annotations = (options.get("architecture_annotations") or {}).get("files", {})
        recorded = [identity.get(k) for k in ("fasta", "fasta_index", "fasta_gzi", "callable_bed")]
        recorded += [annotations.get(k) for k in ("erv_joined", "architecture_segdup_bed", "architecture_lcr_bed", "architecture_low_mappability_bed", "EUR")]
        inputs += [resolve_run_path(run, r["path"]) for r in recorded if isinstance(r, dict) and r.get("path")]
    if (targeted or kind == "candidates") and options.get("architecture_annotation_cache"):
        cache = Path(options["architecture_annotation_cache"])
        marker = cache / "manifest.json"
        if targeted and cache.is_dir():
            inputs.append(marker)
        if marker.is_file():
            inputs.append(marker)
            # prepare_joined_cache validates every committed chromosome file.
            # Candidate selection can instead rebuild an invalid optional split
            # cache, so missing optional files are omitted in that branch only.
            cached = json.loads(marker.read_text())
            inputs += [cache / name for name in cached.get("files", {})
                       if targeted or (cache / name).is_file()]
    return list(dict.fromkeys(inputs))


@dataclass
class Stage:
    name: str
    command: list
    inputs: list
    out: Path
    modules: tuple = ()
    packages: tuple = ()
    executables: tuple = ()
    manifest: str | None = None
    details: dict | None = None


def plan(args, paths, names):
    py = sys.executable
    module_out, modern = paths["out"], paths["modern"]
    out = module_out / "runtime"
    report = str(SOURCE / 'herv.py')
    commands = {}
    ibdmix = "TRUE" if args.ibdmix else "FALSE"
    if "archaic" in names and not args.ibdmix:
        raise ValueError("Archaic analysis requires --ibdmix TRUE")
    common_report = ["--analysis-root", str(paths["root"]), "--modern-dir", str(modern)]
    commands["primary"] = Stage("primary", [py, str(SOURCE / 'herv.py'),'primary',"--modern-dir", str(modern), "--out", str(out / "primary"), "--bootstrap", str(args.bootstrap), "--seed", str(args.seed), "--ibdmix", ibdmix],
        modern_files(modern, "primary") if "primary" in names else [], out / "primary",
        ('herv.py', 'herv.py', "scan.py", "scan_store.py", "0.common.py"),
        ("numpy", "matplotlib"), manifest="primary_manifest.json")
    commands["summary"] = Stage("summary", [py, str(SOURCE / 'herv.py'),'modern',"summary", "--modern-dir", str(modern),
        "--out", str(out / "summary"), "--bootstrap", str(args.bootstrap), "--seed", str(args.seed), "--min-training-sites", str(args.min_training_sites), "--min-support", str(args.min_support),
        "--archaic-references", args.archaic_references, "--ibdmix", ibdmix],
        modern_files(modern, "summary") if "summary" in names else [], out / "summary", ('herv.py', "scan.py", "scan_store.py", "0.common.py"), ("numpy", "scipy", "matplotlib"), manifest="summary_manifest.json")
    if "archaic" in names:
        datasets = archaic_discover(args.archaic_root, args.archaic_datasets, args.archaic_module, args.archaic_chromosomes,
                            args.archaic_references, getattr(args, 'archaic_source_build', 'auto'))
        ready = [d for d in datasets if d['status'] == 'ready']
        needs_lift = any(d['build'] != '38' for d in ready)
        cache = args.archaic_cache or str(paths['root'] / 'archaic_overlap')
        archaic_cmd = [py, str(SOURCE / 'herv.py'),'archaic','--gu-root', args.archaic_root,
            '--datasets', args.archaic_datasets, '--module', args.archaic_module,
            '--source-build', getattr(args, 'archaic_source_build', 'auto'),
            '--chromosomes', args.archaic_chromosomes, '--references', args.archaic_references,
            '--modern-dir', str(modern), '--out', str(out/'archaic'), '--cache-dir', cache,
            '--chain', args.chain37to38, '--liftover', args.archaic_liftover, '--bcftools', args.bcftools,
            '--min-lod', str(args.archaic_min_lod), '--min-length', str(args.archaic_min_length),
            '--carrier-thresholds', args.archaic_carrier_thresholds, '--bootstrap', str(args.bootstrap),
            '--seed', str(args.archaic_seed), '--min-training-sites', str(args.archaic_min_training_sites)]
        archaic_inputs = [*archaic_modern_inputs(modern), *(Path(p) for d in ready for p in d['inputs'])]
        if needs_lift: archaic_inputs.append(Path(args.chain37to38))
        commands['archaic'] = Stage('archaic', archaic_cmd, archaic_inputs, out/'archaic',
            ('herv.py', 'herv.py', 'scan.py', 'scan_store.py', '0.common.py'),
            ('numpy', 'matplotlib', 'pysam'), ('sort', args.bcftools) + ((args.archaic_liftover,) if needs_lift else ()),
            manifest='manifest.json', details={'datasets': [{k: v for k, v in d.items() if k not in {'files', 'inputs', 'sample_panel'}} for d in datasets],
                'dataset_cache': cache, 'gnomad_access': 'indexed tract-union queries; completed background remains read-only'})
    mask_overrides = {name: getattr(args, name + "_bed", None)
                      for name in ("segdup", "lcr", "low_mappability")}
    mask_flags = [value for name, path in mask_overrides.items() if path
                  for value in ("--" + name.replace("_", "-") + "-bed", str(path))]
    commands["candidates"] = Stage("candidates", [py, str(SOURCE / 'herv.py'),'candidates',"--modern-dir", str(modern),
        "--erv-joined", args.erv_joined, "--out", str(out / "candidates"), "--max-sites", str(args.max_sites),
        "--controls-per-case", str(args.controls_per_case), "--seed", str(args.seed),
        "--end-bp", str(args.end_bp), "--body-bins", str(args.body_bins), *mask_flags],
        [*(modern_files(modern, "candidates", mask_overrides) if "candidates" in names else []), Path(args.erv_joined)], out / "candidates", ('herv.py', 'herv.py', "scan.py", "scan_store.py", "scan_store.py", "0.common.py"), ("numpy",), manifest="candidate_manifest.json")
    commands["alleles"] = Stage("alleles", [py, str(SOURCE / 'herv.py'),'candidates',"--modern-dir", str(modern),
        "--erv-joined", args.erv_joined, "--out", str(out / "alleles"), "--max-sites", str(args.allele_max_sites),
        "--controls-per-case", str(args.allele_controls_per_case), "--all-alleles-catalog", "--seed", str(args.seed),
        "--end-bp", str(args.end_bp), "--body-bins", str(args.body_bins), *mask_flags],
        [*(modern_files(modern, "candidates", mask_overrides) if "alleles" in names else []), Path(args.erv_joined)], out / "alleles",
        ('herv.py', 'herv.py', "scan.py", "scan_store.py", "scan_store.py", "0.common.py"), ("numpy",), manifest="candidate_manifest.json")
    kg_inputs = [paths["targets"], Path(args.kg_sample_panel)]
    # Only chromosomes represented in frozen targets are required for the targeted fetch.
    chromosomes = []
    if "genotypes" in names and paths["targets"].is_file():
        import csv
        with paths["targets"].open() as handle:
            chromosomes = sorted({r["chrom"].removeprefix("chr") for r in csv.DictReader(handle, delimiter="\t")}, key=lambda x: int(x))
    else:
        chromosomes = [str(x) for x in range(1, 23)]
    for chrom in chromosomes:
        vcf = Path(args.kg_vcf_dir) / args.kg_vcf_pattern.format(chrom=chrom)
        kg_inputs += [vcf, first_file([str(vcf) + ".tbi", str(vcf) + ".csi"]) or Path(str(vcf) + ".tbi")]
    commands["genotypes"] = Stage("genotypes", [py, str(SOURCE / 'herv.py'),'genotypes',"--targets", str(paths["targets"]),
        "--vcf-dir", args.kg_vcf_dir, "--sample-panel", args.kg_sample_panel, "--vcf-pattern", args.kg_vcf_pattern,
        "--out", str(out / "genotypes"), "--seed", str(args.seed), "--min-ac", str(args.kg_min_ac), "--min-call-rate", str(args.min_call_rate)],
        kg_inputs, out / "genotypes", ('herv.py',), ("pysam", "numpy"), manifest="genotypes_manifest.json")
    window_cmd = [py, str(SOURCE / 'herv.py'),'modern',"windows", "--modern-dir", str(modern), "--erv-joined", args.erv_joined,
        "--fasta", args.fasta38, "--out", str(out / "windows"), "--window", str(args.window), "--step", str(args.step), "--flank", str(args.flank),
        "--bootstrap", str(args.bootstrap), "--seed", str(args.seed), "--min-training-sites", str(args.min_training_sites),
        "--min-support", str(args.min_support), "--bcftools", args.bcftools]
    if args.chromosomes:
        window_cmd += ["--chromosomes", args.chromosomes]
    target_inputs = [Path(args.erv_joined), Path(args.fasta38), Path(args.fasta38 + ".fai")]
    if args.fasta38.endswith((".gz", ".bgz")):
        target_inputs.append(Path(args.fasta38 + ".gzi"))
    for key in ("segdup", "lcr", "low_mappability"):
        path = getattr(args, key + "_bed")
        if path:
            window_cmd += ["--" + key.replace("_", "-") + "-bed", path]
            target_inputs.append(Path(path))
    window_inputs = [*(modern_files(modern, "windows") if "windows" in names else []), *target_inputs]
    commands["windows"] = Stage("windows", window_cmd, window_inputs, out / "windows", ('herv.py', "scan.py", "scan_store.py", "0.common.py"),
                                 ("numpy", "scipy", "matplotlib") + (("pysam",) if args.fasta38.endswith((".gz", ".bgz")) else ()), (args.bcftools,), manifest="windows_manifest.json")
    spatial_cmd = [py, str(SOURCE / 'herv.py'),'modern',"spatial", "--modern-dir", str(modern), "--erv-joined", args.erv_joined,
        "--fasta", args.fasta38, "--out", str(out / "spatial"), "--end-bp", str(args.end_bp), "--body-bins", str(args.body_bins),
        "--bootstrap", str(args.bootstrap), "--seed", str(args.seed), "--min-training-sites", str(args.min_training_sites),
        "--min-support", str(args.min_support), "--bcftools", args.bcftools]
    if args.chromosomes:
        spatial_cmd += ["--chromosomes", args.chromosomes]
    for key in ("segdup", "lcr", "low_mappability"):
        if getattr(args, key + "_bed"):
            spatial_cmd += ["--" + key.replace("_", "-") + "-bed", getattr(args, key + "_bed")]
    spatial_inputs = [*modern_files(modern, "spatial"), *target_inputs] if "spatial" in names else []
    # Shared regional inputs plus the ERV/exposure counters used for reconciliation.
    commands["spatial"] = Stage("spatial", spatial_cmd, list(dict.fromkeys(spatial_inputs)), out / "spatial",
        ('herv.py', "scan.py", "scan_store.py", "0.common.py"), ("numpy", "scipy", "matplotlib") + (("pysam",) if args.fasta38.endswith((".gz", ".bgz")) else ()),
        (args.bcftools,), manifest="spatial_manifest.json")
    time_cmd = [py, str(SOURCE / 'chronology.py'),'allele',"run", "--catalog", str(paths["catalog"]), "--pgen-prefix", str(paths["prefix"]),
        "--metadata", str(paths["metadata"]), "--hq-variants", str(paths["hq"]), "--fasta38", args.fasta38, "--fasta37", args.fasta37,
        "--chain38to37", args.chain38to37, "--chain37to38", args.chain37to38, "--out", str(out / "time"),
        "--bin-years", str(args.bin_years), "--min-samples", str(args.min_samples), "--min-call-rate", str(args.min_call_rate),
        "--max-age-bp", str(args.max_age_bp), "--chunk-size", str(args.time_chunk_size), "--pcs", str(args.time_pcs),
        "--bootstrap", str(args.bootstrap), "--seed", str(args.seed)]
    time_inputs = [paths["catalog"], Path(str(paths["prefix"]) + ".pgen"), Path(str(paths["prefix"]) + ".psam"), paths["pvar"],
        paths["hq"], Path(args.fasta38), Path(args.fasta38 + ".fai"), Path(args.fasta37), Path(args.fasta37 + ".fai"), Path(args.chain38to37), Path(args.chain37to38)]
    time_inputs += [Path(path + ".gzi") for path in (args.fasta38, args.fasta37) if path.endswith((".gz", ".bgz"))]
    catalog_manifest = paths["catalog"].parent / "candidate_manifest.json"
    if catalog_manifest.is_file():
        time_inputs.append(catalog_manifest)
    if paths["metadata_mode"] == "prepare_from_original_supplement":
        time_inputs += [paths["supplement"]]
    else:
        time_inputs += [paths["metadata"]]
    if paths["coordinates"]:
        time_cmd += ["--coordinates", str(paths["coordinates"])]
        time_inputs.append(paths["coordinates"])
    for option in ('basis_dir', 'coordinate_receipt', 'provenance_manifest'):
        value = getattr(args, option, None)
        if value:
            time_cmd += ['--' + option.replace('_', '-'), value]
            if option != 'basis_dir': time_inputs.append(Path(value))
    commands["time"] = Stage("time", time_cmd, time_inputs, out / "time", ('chronology.py', 'chronology.py', "evolve.py", "0.common.py"),
        ("numpy", "scipy", "matplotlib", "pgenlib", "pyliftover", "pysam") + (("openpyxl",) if paths["metadata_mode"] == "prepare_from_original_supplement" else ())
        + (("zstandard",) if any(str(p).endswith(".zst") for p in time_inputs) else ()), manifest="manifest.json")
    selection_inputs = [paths["catalog"], paths["selection"], Path(args.fasta38), Path(args.fasta38 + ".fai"),
        Path(args.fasta37), Path(args.fasta37 + ".fai"), Path(args.chain38to37), Path(args.chain37to38)]
    selection_inputs += [Path(p+".gzi") for p in (args.fasta38,args.fasta37) if p.endswith((".gz",".bgz"))]
    commands["selection"] = Stage("selection", [py, str(SOURCE / 'chronology.py'),'selection',"--catalog", str(paths["catalog"]),
        "--selection-stats", str(paths["selection"]), "--fasta38", args.fasta38, "--fasta37", args.fasta37,
        "--chain38to37", args.chain38to37, "--chain37to38", args.chain37to38,
        "--effect-allele", args.selection_effect_allele, "--out", str(out/"selection")],
        selection_inputs, out/"selection", ('chronology.py', 'chronology.py', "0.common.py"), ("pyliftover", "pysam")
        + (("zstandard",) if any(str(p).endswith(".zst") for p in selection_inputs) else ()), manifest="manifest.json")
    ld_cmd = [py, str(SOURCE / 'herv.py'),'ld',"--ld-block-dir", args.ld_block_dir, "--erv-joined", args.erv_joined,
              "--build", "38", "--out", str(out / "ld"), "--chrom-sizes", str(paths["chrom_sizes"]), "--bootstrap", str(args.ld_bootstrap), "--seed", str(args.seed)]
    ld_inputs = [Path(args.ld_block_dir) / f"{ancestry}.38.bed" for ancestry in ("AFR", "EAS", "EUR", "SAS")]
    ld_inputs += [Path(args.erv_joined), paths["chrom_sizes"]]
    if args.ld_covariates:
        ld_cmd += ["--covariates", args.ld_covariates]
        ld_inputs.append(Path(args.ld_covariates))
    commands["ld"] = Stage("ld", ld_cmd, ld_inputs, out / "ld", ('herv.py', "scan.py", "0.common.py"), ("numpy", "scipy", "matplotlib"), manifest="manifest.json")
    # Optional regional association. A request never silently starts VCF counting.
    commands["report"] = Stage("report", [py, report, "report", "report", *common_report, "--out", str(module_out), "--ibdmix", ibdmix], [], module_out, ('herv.py', 'herv.py', "scan.py", "0.common.py"), ("numpy", "matplotlib"))
    if getattr(args, 'pca_temporal_dir', None):
        commands['report'].command += ['--pca-temporal-dir', args.pca_temporal_dir]
    return [commands[name] for name in names]


def availability(stage):
    missing = [str(p) for p in stage.inputs if not Path(p).is_file()]
    missing += [str(SOURCE / m) for m in stage.modules if not (SOURCE / m).is_file()]
    packages = [name for name in stage.packages if importlib.util.find_spec(name) is None]
    programs = [name for name in stage.executables if shutil.which(name) is None]
    return dict(missing_inputs=list(dict.fromkeys(missing)), missing_packages=packages, missing_programs=programs)


def snapshot(stage):
    inputs = list(dict.fromkeys(stage.inputs + [SOURCE / p for p in stage.modules] + [SOURCE / p for p in ("batch_workflows.py", "analysis_integrity.py", 'pca.py', 'pca.py', "allele_measurement.py", 'chronology.py', "allele_measurement.py", 'herv.py')]))
    return dict(schema=VERSION, command=stage.command, inputs=[state(p) for p in inputs])


def outputs(directory):
    return {str(p):state(p) for p in directory.rglob("*") if p.is_file() and not p.name.endswith((".partial", ".lock")) and "logs" not in p.parts}


def completed_result(stage):
    if not stage.manifest:
        if stage.name == "report" and not (stage.out / "reports/herv_enrichment.html").is_file():
            raise RuntimeError("Report command exited without report.html")
        return {}
    path = stage.out / stage.manifest
    if not path.is_file():
        if not path or not path.is_file():
            raise RuntimeError(f"{stage.name} exited successfully without required manifest: {stage.manifest}")
    result = validate_manifest(stage.out, stage.manifest)
    allowed = {
        'primary': {'complete'}, 'summary': {'complete'}, 'spatial': {'complete'}, 'windows': {'complete'},
        'alleles': {'complete'}, 'candidates': {'complete'}, 'selection': {'complete'}, 'ld': {'complete'},
        'genotypes': {'complete', 'not_run'},
        'time': {'complete_descriptive_imputed_panel', 'complete_descriptive_all_allele_panel',
                 'complete_adjustment_unavailable', 'no_usable_exact_high_quality_contrast', 'no_usable_contrast'},
    }
    if result.get('status') not in allowed.get(stage.name, {'complete'}):
        raise RuntimeError(f'{stage.name} has invalid terminal status: {result.get("status")}')
    return {"result_status":result.get("status", "unknown"), "result_manifest":str(path)}


def run_stage(stage, receipts, always=False):
    available = availability(stage)
    if any(available.values()):
        print(f"[not_available] {stage.name}: inputs={len(available['missing_inputs'])}, packages={available['missing_packages']}, programs={available['missing_programs']}", flush=True)
        return dict(stage=stage.name, status="not_available", **available)
    identity = snapshot(stage)
    receipt = receipts / (stage.name + ".json")
    if not always and receipt.is_file():
        old = json.loads(receipt.read_text())
        if old.get("identity") == identity and old.get("outputs") and all(state(path) == saved for path, saved in old["outputs"].items()):
            completed_result(stage)
            print(f"[reuse] {stage.name}: input/parameter/output content identities match", flush=True)
            return dict(stage=stage.name, status="reused", **old.get("result", {}))
    print(f"[stage] {stage.name}: {shlex.join(stage.command)}", flush=True)
    start = now()
    try:
        subprocess.run(stage.command, check=True)
        result = completed_result(stage)
        if snapshot(stage) != identity:
            raise RuntimeError(f"{stage.name}: inputs changed during execution")
        if stage.name != "report":
            produced = outputs(stage.out)
            if not produced:
                raise RuntimeError(f"{stage.name}: no outputs produced")
            save(receipt, dict(identity=identity, outputs=produced, result=result, started=start, completed=now()))
        return dict(stage=stage.name, status="complete", started=start, completed=now(), **result)
    except subprocess.CalledProcessError as exc:
        return dict(stage=stage.name, status="failed", exit_code=exc.returncode, reason="stage command failed; see complete log", started=start, completed=now())
    except (OSError, ValueError, RuntimeError) as exc:
        return dict(stage=stage.name, status="failed", reason=str(exc), started=start, completed=now())


def prepare_metadata_if_needed(paths, receipts):
    if paths["metadata_mode"] != "prepare_from_original_supplement":
        return
    identity = dict(schema=VERSION, input=[state(str(paths["prefix"]) + ".psam"), state(paths["supplement"]), state(SOURCE / "evolve.py")])
    receipt = receipts / "metadata.json"
    if receipt.is_file() and paths["metadata"].is_file():
        old = json.loads(receipt.read_text())
        if old.get("identity") == identity and old.get("output") == state(paths["metadata"]):
            return
    from evolve import prepare_metadata
    print("[metadata] Preparing canonical sample table from the original supplement into the new output", flush=True)
    prepare_metadata(str(paths["prefix"]) + ".psam", paths["supplement"], paths["metadata"].parent)
    save(receipt, dict(identity=identity, output=state(paths["metadata"])))


def overall(records):
    if any(r["status"] == "failed" for r in records):
        return "failed"
    if any(r["status"] == "not_available" for r in records):
        return "partial_missing_inputs"
    return "complete_requested_stages"


def execute(args):
    paths = resolve(args)
    names = stage_selection(args)
    configuration = {k:str(v) if isinstance(v, Path) else v for k,v in paths.items()}
    configuration.update(stages=names, version=VERSION, ibdmix_enabled=args.ibdmix,
                         disabled_stages=[] if args.ibdmix else ["archaic"],
                         reuse="size/mtime and parameters; no added large-file checksum")
    if args.dry_run or args.preflight:
        stages = plan(args, paths, names)
        details = []
        for stage in stages:
            details.append(dict(stage=stage.name, command=shlex.join(stage.command), **availability(stage),
                                **({'details': stage.details} if stage.details is not None else {})))
        print(json.dumps(dict(status="plan_only", configuration=configuration, stages=details,
            notes=["No files written or analyses started", "Run the integrated pipeline after the scan and PCA inputs have completed; primary reads the selected current scan counters", "Full allele catalog and small validation targets generated earlier in this plan may be absent at preflight", "Time uses the full allele catalog and all eligible ancient samples, not the small genotype-validation panel", "Existing shared PCA is reused; absent coordinates do not trigger PCA fitting", "Actual execution through allele3.sh retains the systemd/cgroup memory cap"]), indent=2))
        return 0
    out = paths["out"]
    out.mkdir(parents=True, exist_ok=True)
    receipts = out / "qc/receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (out / "qc/pipeline.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another focused pipeline is using {out}") from exc
        records = []
        status_path = out / "qc/pipeline_status.json"
        for name in names:
            if args.command == "pipeline" and not args.stages and records and records[-1]["status"] not in {"complete", "reused"}:
                records.append(dict(stage=name, status="not_available",
                                    reason="previous module did not complete; all-module execution stopped"))
                save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
                continue
            # Re-resolve generated candidate files and their exact chromosomes.
            try:
                stage = plan(args, paths, [name])[0]
            except (OSError, ValueError, RuntimeError) as exc:
                records.append(dict(stage=name, status="failed", reason="input planning: " + str(exc)))
                save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
                continue
            candidate_status = next((r["status"] for r in records if r["stage"] == "candidates"), None)
            if name == "genotypes" and not args.targets and candidate_status not in {None, "complete", "reused"}:
                records.append(dict(stage=name, status="not_available", reason="candidates stage did not complete; stale candidate output is not reused"))
                save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
                continue
            catalog_status = next((r["status"] for r in records if r["stage"] == "alleles"), None)
            if name in {"time", "selection"} and not args.allele_catalog and catalog_status not in {None, "complete", "reused"}:
                records.append(dict(stage=name, status="not_available", reason="full allele catalog stage did not complete; stale catalog is not reused"))
                save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
                continue
            if name == "time" and not any(availability(stage).values()):
                try:
                    prepare_metadata_if_needed(paths, receipts)
                    if paths["metadata"] not in stage.inputs:
                        stage.inputs.append(paths["metadata"])
                except (OSError, ValueError, RuntimeError) as exc:
                    records.append(dict(stage="time", status="failed", reason="metadata adapter: " + str(exc)))
                    save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
                    continue
            if name == "report":
                save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
            # Archaic checks its persistent per-cohort caches itself, including
            # their inventories, even when the aggregate table is unchanged.
            result = run_stage(stage, receipts, always=name in {"report", "archaic"})
            records.append(result)
            save(status_path, dict(status=overall(records), stages=records, configuration=configuration, updated=now()))
        final = overall(records)
        print(f"[herv-focus] {final}: {status_path}", flush=True)
        return 2 if final == "failed" else 3 if final == "partial_missing_inputs" else 0


def parser():
    p = argparse.ArgumentParser(description=_description_herv_pipeline)
    p.add_argument("command", nargs="?", choices=("pipeline", *STAGES), default="pipeline")
    p.add_argument("--stages", help="Comma-separated stage subset; only requested stages run")
    p.add_argument("--ibdmix", type=boolean_option, default=False, metavar="TRUE|FALSE",
                   help="Enable exploratory IBDmix evidence and archaic analysis [FALSE]")
    p.add_argument("--dir-analysis-root", default=os.environ.get("ALLELE_ANALYSIS_ROOT", configured_path('ANALYSIS_ROOT', 'allele3')))
    p.add_argument("--dir-modern", "--modern-dir", "--scan-dir", default=os.environ.get("ALLELE_SCAN_STORE", SCAN_STORE))
    p.add_argument("--dir-out", "--out", default=os.environ.get("ALLELE_HERV_FOCUS_OUT"),
                   help="Output directory; the analysis root expands to ROOT/herv_enrichment")
    p.add_argument("--dir-ancient", default=os.environ.get("ALLELE_ANCIENT_DIR", configured_path('DIR_GEN', 'ancient')))
    p.add_argument("--erv-joined", default=os.environ.get("ALLELE_ERV_JOINED", configured_path('DIR_ANNOT', 'allele3/repeats/hg38/rmskJoinedCurrent.txt.gz')))
    p.add_argument("--fasta38", default=os.environ.get("ALLELE_FASTA", configured_path('DIR_GEN', '1hgp/38/GRCH38.fasta.gz')))
    p.add_argument("--fasta37", default=os.environ.get("ALLELE_FASTA37", configured_path('DIR_GEN', '1hgp/37/GRCH37.fasta.gz')))
    p.add_argument("--chain38to37", default=os.environ.get("ALLELE_CHAIN38TO37", configured_path('DIR_FILES', 'liftOver/hg38ToHg19.over.chain.gz')))
    p.add_argument("--chain37to38", default=os.environ.get("ALLELE_CHAIN37TO38", configured_path('DIR_FILES', 'liftOver/hg19ToHg38.over.chain.gz')))
    p.add_argument("--ld-block-dir", default=os.environ.get("ALLELE_LD_BLOCK_DIR", configured_path('DIR_LD', 'block')))
    for name in ("pgen-prefix", "hq-variants", "metadata", "coordinates", "supplement-xlsx", "targets", "allele-catalog", "selection-stats", "ld-covariates", "chrom-sizes", "chromosomes"):
        p.add_argument("--" + name)
    for name in ("basis-dir", "coordinate-receipt", "provenance-manifest", "pca-temporal-dir"):
        p.add_argument("--" + name)
    p.add_argument("--selection-effect-allele", choices=("raw_only", "ALT"), default="raw_only", help="Keep published signs raw unless the source's ALT effect contract has been verified")
    for name in ("segdup", "lcr", "low-mappability"):
        p.add_argument("--" + name + "-bed")
    p.add_argument("--kg-vcf-dir", default=os.environ.get("ALLELE_1KG_VCF_DIR", configured_path('DIR_GEN', '1kg/38/vcf')))
    p.add_argument("--kg-sample-panel", default=os.environ.get("ALLELE_1KG_PANEL", configured_path('DIR_GEN', '1kg/38/samples.txt')))
    p.add_argument("--kg-vcf-pattern", default="chr{chrom}.vcf.gz")
    p.add_argument("--kg-min-ac", type=int, default=2)
    p.add_argument('--archaic-root', default=os.environ.get('ALLELE_ARCHAIC_ROOT', str(Path(os.environ.get('DIR_GEN', configured_path('DIR_GEN'))))))
    p.add_argument('--archaic-datasets', default=os.environ.get('ALLELE_ARCHAIC_DATASETS', '1kg'), help='Comma-separated names/globs [default: %(default)s]; auto selects all complete 1kg/ukbNNN datasets')
    p.add_argument('--archaic-references', default='Altai,Chagyr,Denisova,Denisova25,Vindija', help='Explicit tested reference set; .2013 references are opt-in only')
    p.add_argument('--archaic-module', choices=('auto', 'ibdmix', 'ibdmix-GPU'), default='auto')
    p.add_argument('--archaic-source-build', choices=('auto', '37', '38'), default='auto',
                   help='Select GU assembly if both 37/ and 38/ layouts exist; checked against caller metadata')
    p.add_argument('--archaic-chromosomes', default='auto', help='auto includes available X/Y separately; or autosomes,X,Y / X,Y / X / Y')
    p.add_argument('--archaic-cache', default=os.environ.get('ALLELE_ARCHAIC_CACHE'), help='Persistent per-cohort products; default ANALYSIS-ROOT/archaic_overlap')
    p.add_argument('--archaic-liftover', default=os.environ.get('ALLELE_ARCHAIC_LIFTOVER', configured_path('DIR_SOFTWARE', 'bin/liftOver')))
    p.add_argument('--archaic-min-lod', type=float, default=4.)
    p.add_argument('--archaic-min-length', type=int, default=50000)
    p.add_argument('--archaic-carrier-thresholds', default='0.01,0.05,0.10,0.20')
    p.add_argument('--archaic-min-training-sites', type=int, default=1)
    p.add_argument('--archaic-seed', type=int, default=20261024)
    for name, default in (("window",200),("step",50),("flank",2000),("bootstrap",2000),("min-training-sites",100),
                          ("max-sites",500),("controls-per-case",4),("bin-years",1000),("min-samples",20),
                          ("allele-max-sites",0),("allele-controls-per-case",1),("max-age-bp",10000),
                          ("time-chunk-size",64),("time-pcs",10),("ld-bootstrap",1000),("end-bp",200),("body-bins",20),("seed",20261027)):
        p.add_argument("--" + name, type=int, default=default)
    p.add_argument("--min-support", type=float, default=.99)
    p.add_argument("--min-call-rate", type=float, default=.9)
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--preflight", action="store_true")
    p.add_argument("--version", action="version", version=VERSION)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        return execute(args)
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"herv-pipeline: {exc}", file=sys.stderr)
        return 2


def _entry_herv_pipeline():
    raise SystemExit(main())



if __name__ == '__main__':
    routes = {'primary': _entry_herv_primary, 'modern': _entry_herv_modern, 'candidates': _entry_herv_candidates, 'genotypes': _entry_herv_genotypes, 'domains': _entry_herv_domains, 'archaic': _entry_herv_archaic, 'ld': _entry_herv_ld, 'cohort': _entry_herv_cohort, 'report': _entry_herv_report, 'pipeline': _entry_herv_pipeline}
    arguments = list(sys.argv[1:])
    action = arguments.pop(0) if arguments and arguments[0] in routes else 'pipeline'
    if action is None:
        raise SystemExit('Choose a stage: ' + ', '.join(routes))
    sys.argv[1:] = arguments
    routes[action]()
