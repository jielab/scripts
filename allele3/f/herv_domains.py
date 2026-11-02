#!/usr/bin/env python3
"""Decompose primary HERV multi-ALT enrichment with external domain annotations.

BED4/6 column 4 must contain reviewed biological labels (e.g. Gag, Pol, Env,
5LTR, 3LTR). The caller supplies the annotation source and GRCh38 build.
Geometric repeat ends and '-int' names are never converted to domain labels.
The original eligible SNV universe is queried over aligned family intervals,
including single-ALT controls and all unannotated/ambiguous residual sites.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
import gzip
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
import scan
import herv_modern as hm
from analysis_integrity import transactional, sha256
from allele_common import FastaContext, IntervalMask, canonical_chrom, af_bin, an_bin, open_text, dump_json
from herv_primary import scoped_counter_totals, reconcile_spatial_to_counters

SCHEMA = "herv-biological-domains-1.0"
RESIDUAL = ("unannotated", "ambiguous", "strand_mismatch")
FIELDS = ["chrom", "block_id", "family", "domain", "sensitivity", "threshold", "context", "af_bin", "an_bin", "n_sites", "n_multi"]


@dataclass(frozen=True)
class Domain:
    chrom: str
    start: int
    end: int
    label: str
    strand: str = "."


def read_domains(path):
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
            if start < 0 or end <= start or not label or label.lower() in {*RESIDUAL, "total"} or strand not in {"+", "-", "."}:
                raise ValueError(f"{path}:{number}: invalid domain interval, reserved label, or strand")
            domain = Domain(canonical_chrom(chrom), start, end, label, strand)
            if domain not in seen:
                domains.append(domain)
                seen.add(domain)
                labels.add(label)
    if not domains:
        raise ValueError("No domain annotation intervals")
    return sorted(domains, key=lambda d: (scan.chrom_sort_key(d.chrom), d.start, d.end, d.label, d.strand)), sorted(labels)


class DomainSweep:
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


def assign_domain(hits, repeat_strands):
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


def scan_chrom(vcf, chrom, annotation, options, args, domains):
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
    header = hm.read_indexed_header(vcf, args.bcftools)
    scan.validate_count_header(header)
    fields = ["AC", "AN"] + [f for f in ("AF", "AS_FilterStatus", "lcr", "segdup") if f in header["info"]]
    elements = scan.parse_joined_rmsk(annotation, chromosome=chrom)
    families = set(args.families.split(","))
    selected = [e for e in elements if e.rep_class == "LTR" and e.family in families]
    intervals = scan.merge_intervals([s for e in selected for s in (e.segments or ((e.start, e.end),))])
    repeats = scan.ERVOverlapSweep(elements)
    domain_sweep = DomainSweep(domains)
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
        writer.writerow(FIELDS)
        def flush():
            writer.writerows([*key, *value] for key, value in sorted(counts.items()))
            counts.clear()
        records = hm.query_regions(vcf, fields, bed, args.bcftools) if intervals else iter(())
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
            if strict_complete and not hm.strict_site_excluded(c, pos, flags, overlap, masks):
                sensitivity.append("strict")
            block = eur.locate(c, pos)
            block_id = block.block_id if block else "outside_partition"
            if previous_block is not None and block_id != previous_block:
                flush()
            previous_block = block_id
            hits = domain_sweep.locate(c, pos)
            for family, strands in by_family.items():
                label = assign_domain(hits, strands)
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
        inputs=[hm.file_state(p) for p in inputs], domain_bed_sha256=sha256(args.domain_bed),
        output=hm.file_state(path), query_scope="all eligible SNVs in selected RepeatMasker aligned-family union"))
    return path


def aggregate(paths, background):
    # Vector: observed, partial expected, eligible, unsupported observed, unsupported eligible.
    by_domain = defaultdict(lambda: defaultdict(lambda: np.zeros(5)))
    raw, outside = defaultdict(lambda: np.zeros(2, dtype=np.int64)), Counter()
    for path in paths:
        for r in hm.rows(path):
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


def summarize_domains(data, outside, labels, families, bootstrap=2000, seed=20261102, min_blocks=20):
    if bootstrap < 2000 or min_blocks < 2:
        raise ValueError("Domain inference requires >=2000 bootstrap draws and >=2 minimum blocks")
    categories = sorted(set(labels)) + list(RESIDUAL)
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
def run(args):
    if args.domain_build not in {"38", "GRCh38", "hg38"} or not args.domain_source.strip():
        raise ValueError("A GRCh38 annotation and explicit domain source are required")
    families = args.families.split(",")
    if not all(families) or len(set(families)) != len(families):
        raise ValueError("Families must be unique nonempty exact RepeatMasker labels")
    domains, labels = read_domains(args.domain_bed)
    original_hash = sha256(args.domain_bed)
    source_run, annotation, tasks, source_chromosomes = hm.targeted_inputs(args)
    options = source_run["options"]
    if not options.get("ld_block_dir") or str(options.get("architecture_build", "38")) != "38":
        raise ValueError("Primary domain decomposition requires the original EUR GRCh38 LD partition")
    if not {5, 10} <= {int(t) for t in options.get("thresholds", [])}:
        raise ValueError("Original counters must contain AC5 and AC10")
    files = [scan_chrom(vcf, chrom, annotation / f"chr{chrom}.tsv.gz", options, args, domains) for chrom, vcf in tasks]
    used = {c for c, _ in tasks}
    background = hm.PrimaryBackground(hm.architecture_paths(args.modern_dir, "global"), families)
    data, raw, outside = aggregate(files, background)
    rows, identities, block_rows = summarize_domains(data, outside, labels, families, args.bootstrap, args.seed, args.min_blocks)
    scope = "all_22_autosomes" if used == hm.AUTOSOMES else "selected_chromosomes_subset"
    reconciliation, mismatches, _ = hm.reconcile_spatial_counters(args.modern_dir, raw, background, used, families,
        comparison_scope=scope, primary_background=background)
    counter_totals = scoped_counter_totals(hm.architecture_paths(args.modern_dir, "exposure"), background, families, used)
    total_audit, total_failures = reconcile_spatial_to_counters(counter_totals, rows)
    out = Path(args.out)
    for row in rows:
        row["comparison_scope"] = scope
    hm.write_rows(out / "domain_enrichment.tsv", rows)
    hm.write_rows(out / "domain_accounting.tsv", identities)
    hm.write_rows(out / "domain_block_counts.tsv.gz", block_rows)
    hm.write_rows(out / "domain_counter_reconciliation.tsv", reconciliation)
    hm.write_rows(out / "domain_counter_mismatches.tsv.gz", mismatches)
    hm.write_rows(out / "domain_total_reconciliation.tsv", total_audit)
    hm.write_rows(out / "domain_annotations.tsv.gz", [dict(chrom=d.chrom, start=d.start, end=d.end, domain=d.label, strand=d.strand) for d in domains])
    if sha256(args.domain_bed) != original_hash:
        raise ValueError("Domain annotation changed during analysis")
    if mismatches or total_failures:
        raise ValueError("Domain totals disagree with original same-scope counters; diagnostic staging is retained")
    manifest = dict(status="complete" if rows else "not_estimable", schema=SCHEMA, domain_source=args.domain_source,
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


def parser():
    p = argparse.ArgumentParser(description=__doc__)
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


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        result = run(args)
        print(json.dumps({k: result[k] for k in ("status", "queried_chromosomes", "domain_labels")}, indent=2))
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        print(f"herv-domains: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
