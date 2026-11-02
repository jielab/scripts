#!/usr/bin/env python3
"""Summarize actual nucleotide counts from genotype_validation.tsv.

This is conditional on the predefined target and genotype-encoded nucleotide
panel. It neither infers REF counts from split population summaries nor dates
mutations. Rarefaction samples chromosome copies, not unrelated individuals.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path
from analysis_integrity import transactional, identity

BASES = "ACGT"
MISSING = {"", ".", "NA", "nan", "None"}


def integer(value, label):
    x = float(value)
    if not math.isfinite(x) or x < 0 or x != int(x):
        raise ValueError(f"{label} must be a finite nonnegative integer")
    return int(x)


def rarefaction(counts, n):
    """Exact finite-population P(>=3 observed bases), each with >=1 copy.

    Inclusion/exclusion uses integer binomial coefficients, so cancellation
    cannot create a negative tiny probability. No independent-ALT assumption.
    """
    values = [int(c) for c in counts if c > 0]
    N, k = sum(values), len(values)
    if n < 1 or n > N:
        return None, None
    choose = lambda x: math.comb(x, n) if x >= n else 0
    denominator = choose(N)
    expected_states = sum(denominator - choose(N-c) for c in values) / denominator
    if k < 3:
        return 0.0, expected_states
    # Each one-state outcome appears k-1 times in the pair sum.
    pair_sum = sum(choose(values[i]+values[j]) for i in range(k) for j in range(i+1, k))
    numerator = denominator - pair_sum + (k-2)*sum(choose(c) for c in values)
    if not 0 <= numerator <= denominator:
        raise AssertionError("Invalid multivariate rarefaction probability")
    return numerator/denominator, expected_states


def summarize_row(row, minimum, target_an, min_call_rate):
    out = {k: row.get(k, "") for k in ("site_id", "chrom", "pos", "ref", "family", "group", "population")}
    out.update(status="unavailable", reason="", min_ac=minimum, rarefy_an=target_an,
               scope="predefined_targets_genotype_encoded_bases")
    if row.get("evaluable") != "1" or row.get("represented_status") != "all_expected_alts_encoded":
        out["reason"] = "target_not_fully_evaluable"
        return out
    if row.get("record_status") != "ok":
        raise ValueError("evaluable=1 conflicts with record_status")
    ref = row["ref"]
    if ref not in BASES or len(ref) != 1:
        raise ValueError("REF must be one DNA nucleotide")
    encoded = set(row["encoded_alts"].split(",")) - {"", "."}
    if not encoded <= set(BASES)-{ref}:
        raise ValueError("Invalid encoded ALT nucleotides")
    counts = {}
    for b in BASES:
        raw = row.get(b+"_AC", "NA")
        if raw in MISSING:
            if b in encoded | {ref}:
                raise ValueError("Missing count for an encoded nucleotide")
            continue
        counts[b] = integer(raw, b+"_AC")
        if b not in encoded | {ref}:
            raise ValueError("Count supplied for an unencoded nucleotide")
    an = integer(row["AN"], "AN")
    call_rate = float(row["call_rate"])
    if an < 1 or sum(counts.values()) != an or not math.isfinite(call_rate) or not 0 <= call_rate <= 1:
        raise ValueError("Counts, AN or call rate are inconsistent")
    if call_rate < min_call_rate:
        out["reason"] = "low_call_rate"
        return out
    order = sorted(counts, key=lambda b: (-counts[b], b))
    f = [counts[b]/an for b in order]
    n_supported = sum(c >= minimum for c in counts.values())
    p3, expected = rarefaction(counts.values(), target_an)
    out.update(status="complete_descriptive", AN=an, call_rate=call_rate,
               n_observed_bases=sum(c > 0 for c in counts.values()),
               n_bases_ac_ge_min=n_supported, at_least_three_supported=int(n_supported >= 3),
               exactly_three_supported=int(n_supported == 3), four_supported=int(n_supported == 4),
               observed_ref_ac=counts[ref], observed_ref_af=counts[ref]/an,
               rank1_tied=int(len(f)>1 and f[0] == f[1]),
               rank3_tied=int(len(f)>3 and f[2] == f[3]),
               entropy_nats=-sum(x*math.log(x) for x in f if x),
               diversity_one_minus_sum_p2=1-sum(x*x for x in f),
               p_rarefied_ge3_each_ac_ge1=p3, expected_rarefied_states=expected,
               rarefaction_status="complete" if p3 is not None else "AN_below_target",
               rank_meaning="current_count_order_not_mutation_order")
    for i in range(4):
        out[f"rank{i+1}_base"] = order[i] if i<len(order) else None
        out[f"rank{i+1}_ac"] = counts[order[i]] if i<len(order) else None
        out[f"rank{i+1}_af"] = f[i] if i<len(order) else None
    return out


def write_table(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, keys, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows({k: "NA" if r.get(k) is None else r.get(k, "NA") for k in keys} for r in rows)


@transactional("balance_manifest.json")
def run(args):
    if args.min_ac < 1 or not 1 <= args.rarefy_an <= 5000 or not 0 <= args.min_call_rate <= 1:
        raise ValueError("Invalid min-ac, rarefy-an (1..5000) or min-call-rate")
    before = identity(args.genotypes, full=True)
    source = gzip.open(args.genotypes, "rt") if str(args.genotypes).endswith(".gz") else open(args.genotypes)
    found, results = {}, []
    with source as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"site_id", "chrom", "pos", "ref", "family", "group", "population", "AN", "call_rate",
                    "A_AC", "C_AC", "G_AC", "T_AC", "encoded_alts", "evaluable", "represented_status", "record_status"}
        if required - set(reader.fieldnames or []):
            raise ValueError("Input is not genotype_validation.tsv; missing: "+str(sorted(required-set(reader.fieldnames or []))))
        for row in reader:
            result = summarize_row(row, args.min_ac, args.rarefy_an, args.min_call_rate)
            key = tuple(result[k] for k in ("chrom", "pos", "ref", "family", "group", "population"))
            if key in found:
                if result != found[key]:
                    raise ValueError("Conflicting duplicate site/population count")
                continue
            found[key] = result
            results.append(result)
    if not results:
        raise ValueError("Empty genotype table")
    groups = defaultdict(list)
    for row in results:
        groups[row["family"], row["group"], row["population"]].append(row)
    summary=[]
    for (fam, group, pop), rr in sorted(groups.items()):
        ok=[r for r in rr if r["status"]=="complete_descriptive"]
        equal=[r for r in ok if r["p_rarefied_ge3_each_ac_ge1"] is not None]
        summary.append(dict(family=fam, group=group, population=pop, n_targets=len(rr), n_evaluable=len(ok),
                            n_three_supported=sum(r["at_least_three_supported"] for r in ok),
                            n_rarefiable=len(equal), rarefy_an=args.rarefy_an,
                            mean_rarefied_ge3_probability=sum(r["p_rarefied_ge3_each_ac_ge1"] for r in equal)/len(equal) if equal else None,
                            interpretation="conditional_target_summary; no selection or genome-wide enrichment inference"))
    out=Path(args.out); out.mkdir(parents=True, exist_ok=True)
    write_table(out/"allele_balance.tsv", results)
    write_table(out/"balance_summary.tsv", summary)
    if before != identity(args.genotypes, full=True):
        raise ValueError("Genotype input changed during analysis")
    manifest=dict(status="complete", input=before, settings=vars(args),
                  estimand="Nucleotide-count diversity conditional on predefined targets and fully encoded callable genotypes",
                  rarefaction="exact multivariate hypergeometric over observed chromosome copies; >=1 per base, unrelatedness and LD not modeled",
                  caveat="Balanced frequencies or three states alone do not establish balancing selection; ranks do not imply ancestral or temporal order")
    (out/"balance_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genotypes", required=True, help="Existing genotypes/genotype_validation.tsv")
    p.add_argument("--out", required=True)
    p.add_argument("--min-ac", type=int, default=2)
    p.add_argument("--min-call-rate", type=float, default=.9)
    p.add_argument("--rarefy-an", type=int, default=100)
    run(p.parse_args(argv))


if __name__ == "__main__":
    main()
