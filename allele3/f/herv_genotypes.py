#!/usr/bin/env python3
"""Validate predefined HERV targets with actual A/C/G/T genotype counts.

Indexed 1KG concordance in a fixed ancestry-balanced subset, not independent
replication unless sample overlap has been ruled out. No genome-wide VCF scan.
"""
from __future__ import annotations
from analysis_integrity import transactional
from allele_measurement import phase_relation

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
import math
import os
from pathlib import Path
import random
import re


GROUPS = ("AFR", "EAS", "EUR", "SAS")
BASES = ("A", "C", "G", "T")
SITE_FIELDS = ["site_id", "chrom", "pos", "ref", "alts", "family", "group", "match_id", "population",
               "n_selected", "n_called", "AN", "call_rate", "A_AC", "C_AC", "G_AC", "T_AC",
               "encoded_alts", "expected_alts_not_encoded", "represented_status", "record_status",
               "representation", "n_phase_comparable", "n_phase_unknown", "n_phase_conflict", "n_implicit_phase_set", "n_explicit_phase_set", "n_missing_samples", "n_conflicting_samples", "n_bases_ac_ge_min",
               "observed_three_base_in_encoded_calls", "confirmed_multi_allelic", "target_alts_retained",
               "evaluable", "min_ac", "min_call_rate", "interpretation"]
SUMMARY_FIELDS = ["family", "group", "population", "n_targets", "n_all_expected_encoded", "n_evaluable",
                  "n_confirmed_multi_allelic", "confirmed_fraction", "confirmed_wilson_low", "confirmed_wilson_high",
                  "n_target_alts_retained", "retained_fraction", "retained_wilson_low", "retained_wilson_high",
                  "min_ac", "min_call_rate", "ci_method", "scope"]
BURDEN_FIELDS = ["sample_id", "population", "family", "group", "n_targets_requested", "n_targets_represented",
                 "called_sites", "called_alleles", "nonref_allele_copies", "nonref_per_called_site",
                 "nonref_allele_fraction", "sample_target_call_rate", "measure"]


def text_open(path):
    return gzip.open(path, "rt") if str(path).endswith((".gz", ".bgz")) else open(path)


def metadata(path):
    p = Path(path).expanduser().resolve(); s = p.stat()
    return {"path": str(p), "size": s.st_size, "mtime_ns": s.st_mtime_ns}


def write_json(path, value):
    p = Path(path); tmp = p.with_name(p.name + ".partial")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(tmp, p)


def write_tsv(path, fields, rows):
    p = Path(path); tmp = p.with_name(p.name + ".partial")
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, p)


def read_panel(path):
    with text_open(path) as f:
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


def read_targets(path):
    targets = []; seen = set()
    with text_open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            chrom = r["chrom"].removeprefix("chr"); pos = int(r["pos"]); ref = r["ref"].upper()
            alts = tuple(r["alts"].upper().split(","))
            if chrom not in {str(i) for i in range(1, 23)} or pos < 1:
                raise ValueError("Targets must have autosomal GRCh38 1-based coordinates")
            if ref not in BASES or len(alts) < 2 or len(set(alts)) != len(alts) or any(a not in BASES or a == ref for a in alts):
                raise ValueError(f"Target must have at least two distinct observed SNV ALTs: {chrom}:{pos}")
            if (chrom, pos) in seen: raise ValueError(f"Duplicate target coordinate: {chrom}:{pos}")
            seen.add((chrom, pos))
            if r["group"] not in {"case", "control"}: raise ValueError("Target group must be case/control")
            targets.append(dict(r, chrom=chrom, pos=pos, ref=ref, alts=",".join(alts),
                                site_id=r.get("site_id") or f"{chrom}:{pos}:{ref}"))
    return sorted(targets, key=lambda r: (int(r["chrom"]), r["pos"]))


def balanced_samples(common_samples, panel, seed):
    members = {g: sorted(s for s in common_samples if panel.get(s) == g) for g in GROUPS}
    n = min(map(len, members.values()))
    if not n: raise ValueError("At least one of AFR/EAS/EUR/SAS has no panel-labelled sample in the common VCF sample set")
    rng = random.Random(seed)
    selected = {g: sorted(rng.sample(members[g], n)) for g in GROUPS}
    return selected, {g: len(v) for g, v in members.items()}


def decode_records(records, target_ref, samples):
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
    if any(len(a) < 2 or len(set(a)) != len(a) or any(x not in BASES for x in a) for a in alleles):
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
        counts = dict.fromkeys(BASES, 0)
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


def evaluate_population(target, decoded, population, samples, min_ac, min_call_rate):
    row = {k: target.get(k, ".") for k in ("site_id", "chrom", "pos", "ref", "alts", "family", "group", "match_id")}
    expected = set(target["alts"].split(",")); encoded = decoded["encoded_alts"]
    missing_alts = expected - encoded
    usable = decoded["record_status"] == "ok"
    represented = "all_expected_alts_encoded" if usable and not missing_alts else "expected_alt_not_encoded" if usable or decoded["record_status"] == "record_not_found" else "unavailable_record"
    calls = [decoded["calls"][s] for s in samples if s in decoded["calls"]]
    n = len(calls); call_rate = n / len(samples) if samples and usable else "NA"
    ac = {b: sum(c[b] for c in calls) for b in BASES}
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
    for b in BASES:
        row[b + "_AC"] = ac[b] if n and b in encoded | {target["ref"]} else "NA"
    return row


def wilson(k, n):
    if not n: return "NA", "NA"
    z = 1.959963984540054; p = k / n; den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return max(0., center-half), min(1., center+half)


def summarize(rows, min_ac, min_call_rate):
    groups = defaultdict(list)
    for row in rows: groups[row["family"], row["group"], row["population"]].append(row)
    result = []
    for (family, group, pop), r in sorted(groups.items()):
        usable = [x for x in r if x["evaluable"]]
        n = len(usable); k = sum(x["confirmed_multi_allelic"] for x in usable)
        kept = sum(x["target_alts_retained"] for x in usable)
        lo, hi = wilson(k, n); rlo, rhi = wilson(kept, n)
        result.append(dict(family=family, group=group, population=pop, n_targets=len(r),
                           n_all_expected_encoded=sum(x["represented_status"] == "all_expected_alts_encoded" for x in r),
                           n_evaluable=n, n_confirmed_multi_allelic=k, confirmed_fraction=k/n if n else "NA",
                           confirmed_wilson_low=lo, confirmed_wilson_high=hi,
                           n_target_alts_retained=kept, retained_fraction=kept/n if n else "NA",
                           retained_wilson_low=rlo, retained_wilson_high=rhi, min_ac=min_ac, min_call_rate=min_call_rate,
                           ci_method="descriptive_Wilson_not_LD_adjusted",
                           scope="predefined_AF_matched_modern_multi_ALT_targets; not_genome_wide_enrichment"))
    return result


def unavailable(status):
    return {"record_status": status, "representation": "unavailable", "encoded_alts": set(),
            "calls": {}, "missing": set(), "conflicting": set(), "phase_comparable": set(), "phase_unknown": set(), "phase_conflict": set()}


def plot_concordance(summary, out):
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
            for i, pop in enumerate(GROUPS):
                r = rr.get(pop); x = i + (j-.5)*.30
                if r and r["n_evaluable"]:
                    y = r["confirmed_fraction"]
                    ax.errorbar(x, y, yerr=[[y-r["confirmed_wilson_low"]], [r["confirmed_wilson_high"]-y]],
                                fmt="o", color=color, capsize=3, label=role if i == 0 else None)
        ax.set_title(family); ax.set_xticks(range(4), GROUPS); ax.set_ylim(-.025, 1.025)
        ax.set_ylabel("Fraction with ≥3 observed nucleotide states"); ax.legend(frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Dataset concordance in predefined targets")
    fig.text(.5, .015, "Equal samples per ancestry; descriptive Wilson intervals do not adjust for LD. Missing ALTs are not zeros.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, .94)); fig.savefig(figure, dpi=180); plt.close(fig)
    return "complete"


@transactional("genotypes_manifest.json")
def run_genotypes(args):
    if args.min_ac < 1 or not 0 < args.min_call_rate <= 1:
        raise ValueError("--min-ac must be positive and --min-call-rate in (0,1]")
    targets = read_targets(args.targets); panel = read_panel(args.sample_panel)
    out = Path(args.out).expanduser().resolve(); out.mkdir(parents=True, exist_ok=True)
    sources = {}; missing_sources = {}; index_paths = {}
    for chrom in sorted({r["chrom"] for r in targets}, key=int):
        path = Path(args.vcf_dir).expanduser() / args.vcf_pattern.format(chrom=chrom)
        index = next((Path(str(path)+x) for x in (".tbi", ".csi") if Path(str(path)+x).is_file()), None)
        if not path.is_file(): missing_sources[chrom] = "missing_vcf"
        elif index is None: missing_sources[chrom] = "missing_vcf_index"
        else: sources[chrom] = path; index_paths[chrom] = index
    identity = {"schema": 2, "code": metadata(__file__),
                "measurement_code": metadata(Path(__file__).with_name("allele_measurement.py")),
                "targets": metadata(args.targets), "panel": metadata(args.sample_panel),
                "vcfs": {c: metadata(p) for c, p in sources.items()}, "indexes": {c: metadata(p) for c, p in index_paths.items()},
                "missing_sources": missing_sources, "seed": args.seed, "min_ac": args.min_ac, "min_call_rate": args.min_call_rate,
                "method": "complete_diploid_A_C_G_T_v2; balanced_fixed_samples; split_unique_ALT_residual_REF; implicit_PS_checked"}
    receipt = out / "genotypes_manifest.json"
    if receipt.is_file():
        old = json.loads(receipt.read_text())
        if old.get("status") == "complete" and old.get("identity") == identity and old.get("outputs"):
            if all((out/p).is_file() and metadata(out/p) == fp for p, fp in old["outputs"].items()):
                print(f"[herv-genotypes] reuse {out}", flush=True); return old
    # Invalidate the previous success before any refreshed output is written.
    # A failed/interrupted direct invocation must not expose an older completed
    # manifest as evidence for the new target set or input files.
    write_json(receipt, {"status": "incomplete", "identity": identity, "outputs": {},
                         "reason": "Refresh started; successful finalization has not occurred"})
    handles = {}; selected = {g: [] for g in GROUPS}; eligible_counts = {}; not_run_reason = ""
    sites = []; sample_rows = []; unknown_panel_samples = set()
    try:
        if sources:
            try: import pysam
            except ImportError as e: raise ValueError("pysam is required (included in allele3 environment.yml)") from e
            for chrom, path in sources.items():
                handles[chrom] = pysam.VariantFile(str(path))
            common = set.intersection(*(set(v.header.samples) for v in handles.values()))
            unknown_panel_samples = common - set(panel)
            selected, eligible_counts = balanced_samples(common, panel, args.seed)
            samples = [s for g in GROUPS for s in selected[g]]
            for handle in handles.values(): handle.subset_samples(samples)
        else:
            samples = []; not_run_reason = "no_indexed_target_VCFs" if targets else "empty_targets"
        burden = defaultdict(lambda: [0, 0, 0])  # represented targets, called sites, nonreference copies
        requests = Counter((r["family"], r["group"]) for r in targets)
        for target in targets:
            chrom = target["chrom"]; decoded = None
            if chrom not in handles:
                decoded = unavailable(missing_sources.get(chrom, "not_assayed"))
            else:
                vf = handles[chrom]
                contig = chrom if chrom in vf.header.contigs else "chr"+chrom if "chr"+chrom in vf.header.contigs else None
                if contig is None:
                    decoded = unavailable("contig_not_in_vcf")
                else:
                    records = [r for r in vf.fetch(contig, target["pos"]-1, target["pos"]) if r.pos == target["pos"]]
                    decoded = decode_records(records, target["ref"], samples)
            for g in GROUPS:
                sites.append(evaluate_population(target, decoded, g, selected[g], args.min_ac, args.min_call_rate))
            represented = decoded["record_status"] == "ok" and set(target["alts"].split(",")) <= decoded["encoded_alts"]
            for g in GROUPS:
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
    if (metadata(args.targets) != identity["targets"]
            or metadata(args.sample_panel) != identity["panel"]
            or any(metadata(p) != identity["vcfs"][c] for c, p in sources.items())
            or any(metadata(p) != identity["indexes"][c] for c, p in index_paths.items())):
        raise ValueError("Genotype inputs changed during analysis")
    summary = summarize(sites, args.min_ac, args.min_call_rate)
    balanced = [dict(sample_id=s, population=g, available_common_samples=eligible_counts[g], selected_per_group=len(selected[g]))
                for g in GROUPS for s in selected[g]]
    write_tsv(out / "genotype_validation.tsv", SITE_FIELDS, sites)
    write_tsv(out / "cohort_retention.tsv", SUMMARY_FIELDS, summary)
    write_tsv(out / "sample_burden.tsv", BURDEN_FIELDS, sample_rows)
    write_tsv(out / "balanced_samples.tsv", ["sample_id", "population", "available_common_samples", "selected_per_group"], balanced)
    plot_status = plot_concordance(summary, out)
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
              "outputs": {name: metadata(out / name) for name in names}}
    write_json(receipt, result)
    print(f"[herv-genotypes] {len(targets)} targets; {len(balanced)} balanced samples; {out / 'genotype_validation.tsv'}", flush=True)
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
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
        run_genotypes(args)
    except (ValueError, OSError, KeyError) as e:
        p.exit(2, f"herv-genotypes: {e}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
