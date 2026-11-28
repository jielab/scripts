#!/usr/bin/env python3
"""Stream population SNV summaries into reusable GRCh38 chromosome shards.

The default inputs are gnomAD PASS VCF records with AC evidence, or the supplied
All of Us frequency TSV. Available columns, not the cohort's name, determine
whether allele counts and singleton status can be evaluated. AF>0 can establish
multiple observed ALT bases even when counts are not provided.
All alternatives at a coordinate are combined before counting occupied bases.
Split records do not establish a common called sample set merely because their
AN values happen to match.  The detailed historical gnomAD workflow remains
available through ``3.allele3_scan.sh scan --mode detailed --target gnomad``.

Public reader: ``iter_sites(store, chromosomes=None, population='ALL')``.
It streams numerically ordered chromosomes, ascending positions within each
chromosome, and the requested populations within each position.  Numeric NA is
None; ``alts`` and the three ``alt_*s`` columns become aligned Python lists.
``base_af``, ``base_ac``, and ``base_an`` map actual nucleotide identities to
their values.  ``primary_af`` is the largest KNOWN ALT AF, not the major-allele
frequency.  Only a verified joint record (or explicit custom joint declaration)
permits REF reconstruction and strict bi/tri/quad-allelic state assignment.
"""
from __future__ import annotations

from product_manifest import atomic_replace, write_inventory
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import dataclass, field
import gzip
import hashlib
import html
import json
import math
import os
from pathlib import Path
from output_layout import data_path
import re
import sys
import time
from population_patterns import (FIELDS as PATTERN_FIELDS, COUNT_FIELDS as PATTERN_COUNTS,
                                 classify_site, add_counts as add_pattern_counts)

VERSION = "summary-scan-20261115.1"
EVIDENCE_SCHEMA = "allele-evidence-20261115.3"
PATTERN_SCHEMA = "population-patterns-20261115.3"
BASES = "ACGT"
NA = {"", ".", "NA", "N/A", "nan", "NaN", "None", "NULL"}
GNOMAD_POPS = (
    "ALL", "AFR", "AMI", "AMR", "ASJ", "EAS",
    "FIN", "MID", "NFE", "OTH", "REMAINING", "SAS",
)
AOU_POPS = ("ALL", "AFR", "AMR", "EAS", "EUR", "MID", "OTH", "SAS")
BASE_SITE_FIELDS = [
    "chrom", "pos", "ref", "population", "alts", "alt_afs", "alt_acs", "alt_ans",
    "n_alt_observed", "n_alt_qualified", "ref_af", "ref_ac", "AN", "observed_bases",
    "qualified_bases", "observed_state", "qualified_state", "singleton_status",
    "count_mode", "threshold_kind", "threshold", "ref_evidence", "flags", "control",
    "primary_af",
]
EVIDENCE_FIELDS = [
    "observed_alt_bases", "zero_alt_bases", "unknown_alt_bases", "unreported_alt_bases",
    "qualified_alt_bases", "record_bases", "record_allele_count", "record_state",
    "ref_plus_observed_alt_bases", "ref_plus_observed_alt_count", "ref_plus_observed_alt_state",
    "observed_alleles_min", "observed_alleles_max", "possible_observed_states",
]
SITE_FIELDS = BASE_SITE_FIELDS + EVIDENCE_FIELDS
EVIDENCE_COUNTS = [
    "sites", "one_observed_alt_sites", "two_observed_alt_sites", "three_observed_alt_sites",
    "at_least_three_observed_bases_sites", "proven_quadallelic_sites",
    "ref_plus_observed_triallelic_sites", "ref_plus_observed_quadallelic_sites",
    "record_quadallelic_sites", "unresolved_occupancy_sites", "reported_zero_alt_sites",
    "reported_unknown_alt_sites", "af_available_sites", "ac_available_sites",
    "an_available_sites", "singleton_known_sites",
]
COUNT_FIELDS = [
    "sites", "sites_with_observed_alt", "sites_with_qualified_alt", "multi_alt_sites",
    "biallelic_sites", "triallelic_sites", "quadallelic_sites", "monomorphic_sites",
    "ref_unknown_sites", "alt_singleton_sites", "singleton_unknown_sites",
    "alt_observed_total", "alt_qualified_total", "flagged_sites",
]
WINDOW_FIELDS = ["chrom", "start", "end", "population", "threshold_kind", "threshold", *COUNT_FIELDS]


def chromosome(value):
    value = str(value).strip()
    value = re.sub(r"^chr", "", value, flags=re.I).upper()
    if value in {"M", "MT"}:
        return "chrM"
    if value in {"X", "Y"}:
        return "chr" + value
    if value.isdigit() and 1 <= int(value) <= 22:
        return "chr" + str(int(value))
    return None


def chromosome_key(value):
    c = chromosome(value)
    if c is None:
        raise ValueError("Unsupported primary chromosome: " + str(value))
    return int(c[3:]) if c[3:].isdigit() else {"chrX": 23, "chrY": 24, "chrM": 25}[c]


def selected_chromosomes(value):
    if value is None or value == "all":
        return ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    if value == "autosomes":
        return ["chr" + str(i) for i in range(1, 23)]
    out = set()
    parts = value.split(",") if isinstance(value, str) else value
    for part in parts:
        match = re.fullmatch(r"(?:chr)?(\d+)-(?:chr)?(\d+)", str(part), re.I)
        if match:
            a, b = map(int, match.groups())
            if not 1 <= a <= b <= 22:
                raise ValueError("Chromosome ranges must be within 1-22")
            out.update("chr" + str(i) for i in range(a, b + 1))
        else:
            c = chromosome(part)
            if c is None:
                raise ValueError("Unsupported chromosome: " + str(part))
            out.add(c)
    return sorted(out, key=chromosome_key)


def open_text(path, mode="rt"):
    path = Path(path)
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix in {".gz", ".bgz"} else path.open(mode, encoding="utf-8", newline="")


def number(value, integer=False):
    if value is None or str(value).strip() in NA:
        return None
    value = str(value).strip()
    if integer:
        if not re.fullmatch(r"\d+", value):
            raise ValueError("Invalid nonnegative integer: " + value)
        return int(value)
    result = float(value)
    if not math.isfinite(result) or result < 0 or result > 1:
        raise ValueError("AF must be finite and in [0,1]: " + value)
    return result


def fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return format(value, ".12g")
    return str(value)


def vector(value, n, integer=False, broadcast=False):
    if value is None or value in NA:
        return [None] * n
    parts = value.split(",")
    if broadcast and len(parts) == 1:
        parts *= n
    if len(parts) != n:
        raise ValueError(f"Expected {n} ALT-aligned values, got {len(parts)}")
    return [number(x, integer) for x in parts]


def same_number(a, b):
    if a is None or b is None:
        return a is b
    return a == b if isinstance(a, int) and isinstance(b, int) else math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-10)


@dataclass
class Observation:
    ac: int | None = None
    an: int | None = None
    af: float | None = None
    conflict: bool = False

    def values(self):
        return self.ac, self.an, self.af


@dataclass
class Record:
    chrom: str
    pos: int
    ref: str
    # Per-population dictionaries include only valid SNV ALT identities.
    observations: dict
    modes: dict
    passed: bool
    joint_record: bool
    flags: set = field(default_factory=set)


@dataclass
class SiteGroup:
    chrom: str
    pos: int
    ref: str
    observations: dict = field(default_factory=dict)
    modes: dict = field(default_factory=dict)
    flags: set = field(default_factory=set)
    first_signature: tuple | None = None
    distinct_records: bool = False
    joint_record: bool = True
    ref_conflict: bool = False
    raw_records: int = 0
    passing_records: int = 0
    duplicate_values: int = 0
    conflicts: int = 0

    def add(self, record):
        self.raw_records += 1
        self.flags.update(record.flags)
        if record.ref != self.ref:
            self.ref_conflict = True
            self.flags.add("REF_conflict")
        self.joint_record &= record.joint_record
        if not record.passed:
            self.flags.add("filtered_record_at_site")
            self.joint_record = False
            return
        self.passing_records += 1
        signature = tuple((p, tuple((b, o.values()) for b, o in sorted(v.items()))) for p, v in sorted(record.observations.items()))
        if self.first_signature is None:
            self.first_signature = signature
        elif self.first_signature != signature:
            self.distinct_records = True
        for population, alternatives in record.observations.items():
            self.modes[population] = record.modes[population]
            dest = self.observations.setdefault(population, {})
            for base, obs in alternatives.items():
                if base not in dest:
                    dest[base] = Observation(*obs.values(), obs.conflict)
                    continue
                old = dest[base]
                if old.conflict:
                    continue
                values = []
                conflict = obs.conflict
                complementary = False
                for a, b in zip(old.values(), obs.values()):
                    if a is None and b is not None:
                        values.append(b)
                        complementary = True
                    elif b is None:
                        values.append(a)
                        complementary |= a is not None
                    elif not same_number(a, b):
                        conflict = True
                        values.append(None)
                    else:
                        values.append(a)
                if conflict:
                    dest[base] = Observation(conflict=True)
                    self.flags.add("conflicting_ALT_summary")
                    self.conflicts += 1
                else:
                    dest[base] = Observation(*values)
                    self.duplicate_values += 1
                    self.flags.add("duplicate_complemented" if complementary else "duplicate_deduplicated")


def info_values(text, wanted):
    # Do not parse AC_raw/AN_raw as filtered counts. Unrelated annotations are
    # never materialized into a per-variant dictionary.
    found = {}
    for token in text.split(";"):
        key, sep, val = token.partition("=")
        if sep and key in wanted:
            if key in found and found[key] != val:
                raise ValueError("Conflicting duplicate INFO field: " + key)
            found[key] = val
    return found


def population_names(config, available):
    requested = config["populations"]
    if requested == ["all"]:
        return list(available)
    missing = set(requested) - set(available)
    if missing:
        raise ValueError("Requested populations not in input schema/header: " + ",".join(sorted(missing)))
    return requested


def vcf_records(path, config, stats, wanted_chroms):
    declared = set()
    pops = None
    wanted = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##"):
                match = re.match(r"##INFO=<ID=([^,>]+)", line)
                if match:
                    declared.add(match.group(1))
                if line.startswith("##reference=") and re.search(r"GRCh37|hg19", line, re.I):
                    raise ValueError(f"Input declares an incompatible GRCh37 reference: {path}")
                continue
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t", 8)
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF line {path}:{line_number}")
            c = chromosome(fields[0])
            if c is None or c not in wanted_chroms:
                continue
            if pops is None:
                available = ["ALL"]
                available += [p for p in GNOMAD_POPS[1:] if "AC_" + p.lower() in declared]
                # Hand-written/custom VCFs without INFO metadata are legal.
                if not declared:
                    available = list(GNOMAD_POPS)
                pops = population_names(config, available)
                wanted = {prefix + ("" if p == "ALL" else "_" + p.lower()) for p in pops for prefix in ("AC", "AN", "AF")}
            try:
                pos = int(fields[1])
            except ValueError as exc:
                raise ValueError(f"Invalid VCF position {path}:{line_number}") from exc
            if pos < 1:
                raise ValueError(f"VCF positions must be positive: {path}:{line_number}")
            stats["source_records_selected_chromosomes"] += 1
            ref, alts = fields[3].upper(), fields[4].upper().split(",")
            if len(ref) != 1 or ref not in BASES:
                stats["non_snv_records"] += 1
                if ref and ref[0] in BASES:
                    yield Record(c, pos, ref[0], {}, {}, False, False, {"non_SNV_record_at_site"})
                continue
            valid = [i for i, a in enumerate(alts) if a in BASES and len(a) == 1 and a != ref]
            if not valid:
                stats["non_snv_records"] += 1
                yield Record(c, pos, ref, {}, {}, False, False, {"non_SNV_record_at_site"})
                continue
            passed = fields[6] == "PASS"
            if not passed:
                stats["non_PASS_records"] += 1
            flags = set()
            joint = len(valid) == len(alts) and len(set(alts)) == len(alts)
            if not joint:
                flags.add("mixed_or_duplicate_ALT_record")
            try:
                vals = info_values(fields[7], wanted)
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc
            observations = {}
            for p in pops:
                suffix = "" if p == "ALL" else "_" + p.lower()
                try:
                    ac = vector(vals.get("AC" + suffix), len(alts), True)
                    an = vector(vals.get("AN" + suffix), len(alts), True, broadcast=True)
                    af = vector(vals.get("AF" + suffix), len(alts))
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_number}, {p}: {exc}") from exc
                out = {}
                for i in valid:
                    invalid = ac[i] is not None and an[i] is not None and ac[i] > an[i]
                    # Tiny positive AF may round to 0, so do not reject AC>0/AF=0.
                    invalid |= ac[i] == 0 and af[i] is not None and af[i] > 0
                    if invalid:
                        out[alts[i]] = Observation(conflict=True)
                        flags.add("invalid_AC_AN_AF")
                        stats["invalid_count_summaries"] += 1
                    else:
                        out[alts[i]] = Observation(ac[i], an[i], af[i])
                    if ac[i] is None:
                        stats["missing_AC_summaries"] += 1
                observations[p] = out
            yield Record(c, pos, ref, observations, {p: "AC" for p in pops}, passed, joint, flags)


def default_aou_schema():
    return dict(format="tsv", build="GRCh38", coordinate_system="1-based", chrom="#Chr", pos="Start", end="End", ref="Ref", alt="Alt",
                populations={p: {"af": "gvs_" + p.lower() + "_af"} for p in AOU_POPS}, joint_counts=False)


def discover_aou_schema(path):
    """Discover standard gvs columns from one header, without scanning variants."""
    with open_text(path) as stream:
        header = next((line.rstrip("\r\n").split("\t") for line in stream if line.strip() and not line.startswith("##")), None)
    if header is None:
        raise ValueError("Empty All of Us TSV: " + str(path))
    schema = default_aou_schema()
    if {"#Chr","Start","Ref","Alt"} <= set(header):
        if "End" not in header:
            schema.pop("end")
    elif {"contig","position","ref_allele","alt_allele"} <= set(header):
        schema.update(chrom="contig", pos="position", ref="ref_allele", alt="alt_allele")
        schema.pop("end")
    else:
        raise ValueError("Unrecognized All of Us coordinate columns; provide --schema with the actual coordinate system and names")
    schema["populations"] = {}
    for p in AOU_POPS:
        columns = {kind: "gvs_"+p.lower()+"_"+kind for kind in ("af","ac","an")
                   if "gvs_"+p.lower()+"_"+kind in header}
        if "af" in columns or "ac" in columns:
            schema["populations"][p] = columns
    if not schema["populations"]:
        raise ValueError("No standard gvs_{population}_af/ac columns found; provide --schema for this All of Us file")
    return schema


def tsv_records(path, config, stats, wanted_chroms):
    schema = config["schema"]
    header = None
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("##"):
                continue
            fields = line.rstrip("\r\n").split("\t")  # Empty tab fields stay missing.
            if header is None:
                header = {name: i for i, name in enumerate(fields)}
                if len(header) != len(fields):
                    raise ValueError("Duplicate TSV column names: " + str(path))
                required = {schema[x] for x in ("chrom", "pos", "ref", "alt")}
                if schema.get("end"):
                    required.add(schema["end"])
                if schema.get("filter_column"):
                    required.add(schema["filter_column"])
                available = [p for p, columns in schema["populations"].items() if all(v in header for v in columns.values())]
                pops = population_names(config, available)
                for p in pops:
                    required.update(schema["populations"][p].values())
                if required - set(header):
                    raise ValueError("Missing TSV columns: " + ",".join(sorted(required - set(header))))
                modes = {p: "AC" if "ac" in schema["populations"][p] else "AF" for p in pops}
                continue
            if len(fields) != len(header):
                raise ValueError(f"Malformed tab-separated row {path}:{line_number}: expected {len(header)} fields, got {len(fields)}")
            c = chromosome(fields[header[schema["chrom"]]])
            if c is None or c not in wanted_chroms:
                continue
            stats["source_records_selected_chromosomes"] += 1
            pos = int(fields[header[schema["pos"]]])
            zero_based = schema.get("coordinate_system", "1-based") == "0-based-half-open"
            end = int(fields[header[schema["end"]]]) if schema.get("end") else None
            ref = fields[header[schema["ref"]]].upper()
            if end is not None and end != pos + int(zero_based):
                stats["non_snv_records"] += 1
                if ref and ref[0] in BASES:
                    yield Record(c, pos + int(zero_based), ref[0], {}, {}, False, False, {"non_SNV_record_at_site"})
                continue
            pos += int(zero_based)
            if pos < 1:
                raise ValueError(f"Invalid TSV position: {path}:{line_number}")
            alts = fields[header[schema["alt"]]].upper().split(",")
            if len(ref) != 1 or ref not in BASES:
                stats["non_snv_records"] += 1
                if ref and ref[0] in BASES:
                    yield Record(c, pos, ref[0], {}, {}, False, False, {"non_SNV_record_at_site"})
                continue
            valid = [i for i, a in enumerate(alts) if len(a) == 1 and a in BASES and a != ref]
            if not valid:
                stats["non_snv_records"] += 1
                yield Record(c, pos, ref, {}, {}, False, False, {"non_SNV_record_at_site"})
                continue
            passed = not schema.get("filter_column") or fields[header[schema["filter_column"]]] in schema.get("pass_values", ["PASS"])
            if not passed:
                stats["non_PASS_records"] += 1
            observations = {}
            flags = set()
            if len(valid) != len(alts) or len(set(alts)) != len(alts):
                flags.add("mixed_or_duplicate_ALT_record")
            for p in pops:
                columns = schema["populations"][p]
                try:
                    ac = vector(fields[header[columns["ac"]]] if "ac" in columns else None, len(alts), True)
                    an = vector(fields[header[columns["an"]]] if "an" in columns else None, len(alts), True, broadcast=True)
                    af = vector(fields[header[columns["af"]]] if "af" in columns else None, len(alts))
                except ValueError as exc:
                    raise ValueError(f"{path}:{line_number}, {p}: {exc}") from exc
                observations[p] = {}
                for i in valid:
                    invalid = ac[i] is not None and an[i] is not None and ac[i] > an[i]
                    invalid |= ac[i] == 0 and af[i] is not None and af[i] > 0
                    observations[p][alts[i]] = Observation(conflict=True) if invalid else Observation(ac[i], an[i], af[i])
                    if invalid:
                        flags.add("invalid_AC_AN_AF")
                        stats["invalid_count_summaries"] += 1
                    if modes[p] == "AF" and af[i] is None:
                        stats["missing_AF_summaries"] += 1
                    if modes[p] == "AC" and ac[i] is None:
                        stats["missing_AC_summaries"] += 1
            joint = bool(schema.get("joint_counts")) and len(valid) == len(alts) and len(set(alts)) == len(alts)
            yield Record(c, pos, ref, observations, modes, passed, joint, flags)
    if header is None:
        raise ValueError("Empty TSV input: " + str(path))


def site_groups(records):
    current = None
    finished_chromosomes = set()
    for record in records:
        if current is None or (record.chrom, record.pos) != (current.chrom, current.pos):
            if current is not None:
                if record.chrom == current.chrom and record.pos < current.pos:
                    raise ValueError(f"Input is not position-sorted: {record.chrom}:{record.pos} after {current.pos}. Sort once upstream; the scanner does not hold a chromosome in RAM.")
                if record.chrom != current.chrom:
                    finished_chromosomes.add(current.chrom)
                    if record.chrom in finished_chromosomes:
                        raise ValueError("Chromosome appears in multiple noncontiguous blocks: " + record.chrom)
                yield current
            current = SiteGroup(record.chrom, record.pos, record.ref)
        current.add(record)
    if current is not None:
        yield current


def stable_hash_fraction(chrom, pos, ref, seed):
    text = f"{seed}|{chrom}|{pos}|{ref}".encode("ascii")
    return int.from_bytes(hashlib.blake2b(text, digest_size=8).digest(), "big") / 18446744073709551616


def state(bases):
    return {0: "unobserved", 1: "monomorphic", 2: "biallelic", 3: "triallelic", 4: "quadallelic"}[len(bases)]


def qualifies(value, threshold):
    return value is not None and value > 0 and value >= threshold


def allele_evidence(row):
    """Positive-base identities and explicit occupancy bounds from stored evidence.

    Bounds concern reported population evidence, not allele chronology. A zero
    AF is a reported zero (it may reflect rounding); an unreported base remains
    unmeasured unless a coherent joint record already establishes occupancy.
    This O(1) transformation also works on v2 shards without rewriting them.
    """
    alts, ref = row["alts"], row["ref"]
    values = row["alt_afs"] if row["count_mode"] == "AF_only" else row["alt_acs"]
    positive = [b for b, v in zip(alts, values) if v is not None and v > 0]
    zeros = [b for b, v in zip(alts, values) if v == 0]
    unknown = [b for b, v in zip(alts, values) if v is None]
    unreported = [b for b in BASES if b != ref and b not in alts]
    if row["threshold_kind"] == "AC":
        qualification = row["alt_acs"]
    elif row["count_mode"] == "AF_only":
        qualification = row["alt_afs"]
    else:
        qualification = [af if ac is not None and ac > 0 else (0. if ac == 0 else None)
                         for ac, af in zip(row["alt_acs"], row["alt_afs"])]
    qualified = [b for b, v in zip(alts, qualification) if qualifies(v, row["threshold"])]
    record = sorted(set([ref, *alts]))
    ref_plus_positive = sorted(set([ref, *positive]))
    exact = row.get("observed_bases")
    if exact is not None:
        lower = upper = len(exact)
    else:
        lower = len(set(positive + ([ref] if row["ref_evidence"] == "observed" else [])))
        upper = 4 - len(set(zeros + ([ref] if row["ref_evidence"] == "zero" else [])))
    return dict(observed_alt_bases=positive, zero_alt_bases=zeros,
        unknown_alt_bases=unknown, unreported_alt_bases=unreported, qualified_alt_bases=qualified,
        record_bases=record, record_allele_count=len(record), record_state=state(record),
        ref_plus_observed_alt_bases=ref_plus_positive, ref_plus_observed_alt_count=len(ref_plus_positive),
        ref_plus_observed_alt_state=state(ref_plus_positive),
        observed_alleles_min=lower, observed_alleles_max=upper,
        possible_observed_states=[state(BASES[:n]) for n in range(lower, upper+1)])


def add_evidence_counts(counter, row):
    e = row
    counter["sites"] += 1
    for n, name in [(1,"one"),(2,"two"),(3,"three")]:
        counter[name + "_observed_alt_sites"] += int(row["n_alt_observed"] == n)
    counter["at_least_three_observed_bases_sites"] += int(e["observed_alleles_min"] >= 3)
    counter["proven_quadallelic_sites"] += int(e["observed_alleles_min"] == e["observed_alleles_max"] == 4)
    counter["ref_plus_observed_triallelic_sites"] += int(e["ref_plus_observed_alt_count"] == 3)
    counter["ref_plus_observed_quadallelic_sites"] += int(e["ref_plus_observed_alt_count"] == 4)
    counter["record_quadallelic_sites"] += int(e["record_allele_count"] == 4)
    counter["unresolved_occupancy_sites"] += int(e["observed_alleles_min"] != e["observed_alleles_max"])
    counter["reported_zero_alt_sites"] += int(bool(e["zero_alt_bases"]))
    counter["reported_unknown_alt_sites"] += int(bool(e["unknown_alt_bases"]))
    for prefix, field in [("af","alt_afs"),("ac","alt_acs"),("an","alt_ans")]:
        counter[prefix + "_available_sites"] += int(any(x is not None for x in row[field]))
    counter["singleton_known_sites"] += int(row["singleton_status"] != "unknown")


def prepare_population(group, population, config):
    observations = group.observations[population]
    flags = set(group.flags)
    source_mode = group.modes[population]
    mode = config.get("threshold_kind", "auto")
    if mode == "auto":
        mode = source_mode
    if mode == "AC" and source_mode == "AF":
        raise ValueError("Cannot apply AC thresholds to an AF-only population summary")
    alts = sorted(observations)
    ac = [observations[b].ac for b in alts]
    an = [observations[b].an for b in alts]
    af = [observations[b].af for b in alts]
    if group.ref_conflict:
        ac = an = af = [None] * len(alts)
    if any(observations[b].conflict for b in alts):
        flags.add("conflicting_or_invalid_ALT_summary")
    known_an = {v for v in an if v is not None}
    if len(known_an) > 1:
        flags.add("ALT_AN_mismatch")
    complete_ac = bool(ac) and all(x is not None for x in ac)
    coherent_an = len(known_an) == 1 and all(x is not None for x in an) and next(iter(known_an)) > 0
    declared_joint = bool(config.get("schema", {}).get("joint_counts"))
    joint = (group.joint_record and (not group.distinct_records or declared_joint)
             and not group.ref_conflict and not any(o.conflict for o in observations.values())
             and complete_ac and coherent_an)
    ref_ac = ref_af = common_an = None
    count_mode = "AF_only" if source_mode == "AF" else "per_ALT_counts"
    evidence = "unknown"
    if joint:
        common_an = next(iter(known_an))
        if sum(ac) <= common_an:
            ref_ac = common_an - sum(ac)
            ref_af = ref_ac / common_an
            count_mode = "joint_declared_counts" if declared_joint else "joint_VCF_record"
            evidence = "observed" if ref_ac > 0 else "zero"
        else:
            flags.add("ALT_AC_sum_exceeds_AN")
            common_an = None
    elif source_mode == "AC" and group.distinct_records:
        flags.add("split_records_joint_denominator_unverified")
    # Missing reported AF remains NA even when AC/AN exist. Consumers may
    # explicitly derive an AF from a valid count; the source field is not erased.
    observed_values = ac if source_mode == "AC" else af
    observed = [b for b, v in zip(alts, observed_values) if v is not None and v > 0]
    if mode == "AC":
        values = ac
    elif source_mode == "AF":
        values = af
    else:
        # Filtered AC, not raw AC or AF alone, establishes observed gnomAD
        # alleles. A documented AC=0 is known absent; a missing AC is unknown.
        values = [v if a is not None and a > 0 else (0.0 if a == 0 else None)
                  for a, v in zip(ac, af)]
    exact = ref_ac is not None and all(v is not None for v in observed_values)
    observed_bases = sorted(observed + ([group.ref] if ref_ac is not None and ref_ac > 0 else [])) if exact else None
    if source_mode == "AF":
        singleton = "unknown"
    elif any(x == 1 for x in ac):
        singleton = "present"
    elif any(x is None for x in ac):
        singleton = "unknown"
    else:
        singleton = "absent"
    return dict(alts=alts, alt_acs=ac, alt_ans=an, alt_afs=af, values=values,
                ref_ac=ref_ac, ref_af=ref_af, AN=common_an, observed_bases=observed_bases,
                observed_state=state(observed_bases) if observed_bases is not None else "unknown",
                n_alt_observed=len(observed), count_mode=count_mode, threshold_kind=mode,
                ref_evidence=evidence, singleton_status=singleton, flags=";".join(sorted(flags)) or ".",
                primary_af=max((x for x in af if x is not None), default=None))


def population_row(group, population, prepared, threshold, config, hash_value):
    p = prepared
    qualified = [b for b, v in zip(p["alts"], p["values"]) if qualifies(v, threshold)]
    ref_value = p["ref_ac"] if p["threshold_kind"] == "AC" else p["ref_af"]
    exact = p["observed_bases"] is not None and all(v is not None for v in p["values"])
    qualified_bases = sorted(qualified + ([group.ref] if qualifies(ref_value, threshold) else [])) if exact else None
    row = dict(chrom=group.chrom, pos=group.pos, ref=group.ref, population=population,
               **{k: v for k, v in p.items() if k != "values"}, n_alt_qualified=len(qualified),
               qualified_bases=qualified_bases, qualified_state=state(qualified_bases) if exact else "unknown",
               threshold=threshold,
               control=int(len(qualified) == 1 and all(v is not None for v in p["values"]) and hash_value < config["control_rate"] and "conflict" not in p["flags"] and "invalid" not in p["flags"]))
    return row


def serialize_row(row):
    result = {}
    for key in SITE_FIELDS:
        value = row.get(key)
        result[key] = ",".join(fmt(v) for v in value) if isinstance(value, list) else fmt(value)
    return result


def add_counts(counter, row):
    counter["sites"] += 1
    counter["sites_with_observed_alt"] += int(row["n_alt_observed"] > 0)
    counter["sites_with_qualified_alt"] += int(row["n_alt_qualified"] > 0)
    counter["multi_alt_sites"] += int(row["n_alt_qualified"] >= 2)
    key = row["qualified_state"] + "_sites"
    if key in COUNT_FIELDS:
        counter[key] += 1
    counter["ref_unknown_sites"] += int(row["ref_evidence"] == "unknown")
    counter["alt_singleton_sites"] += int(row["singleton_status"] == "present")
    counter["singleton_unknown_sites"] += int(row["singleton_status"] == "unknown")
    counter["alt_observed_total"] += row["n_alt_observed"]
    counter["alt_qualified_total"] += row["n_alt_qualified"]
    counter["flagged_sites"] += int(row["flags"] != ".")


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".partial")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temp, path)


def file_identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def shard_paths(out, chrom):
    out = Path(out)
    return dict(sites=out / "shards" / (chrom + ".sites.tsv.gz"),
                controls=out / "shards" / (chrom + ".controls.tsv.gz"),
                windows=out / "shards" / (chrom + ".windows.tsv.gz"),
                patterns=out / "shards" / (chrom + ".population_patterns.tsv.gz"),
                receipt=data_path(out, 'qc') / (chrom + ".json"))


def usable_receipt(out, chrom, source=None, config=None):
    paths = shard_paths(out, chrom)
    try:
        item = json.loads(paths["receipt"].read_text())
        if item.get("status") != "complete" or item.get("version") != VERSION:
            return None
        if source is not None and item.get("source") != source:
            return None
        if config is not None and scientific_config(item.get("config", {})) != scientific_config(config):
            return None
        for k in ("sites", "controls", "windows", *(["patterns"] if "patterns" in item["outputs"] else [])):
            if not paths[k].is_file() or paths[k].stat().st_size != item["outputs"][k]["size"]:
                return None
        return item
    except (OSError, ValueError, KeyError, TypeError):
        return None


def scientific_config(config):
    # Chromosome assignment and progress cadence may differ across independent
    # jobs that are subsequently merged; biological settings may not.
    return {k: v for k, v in config.items() if k not in {"chromosomes", "progress_every"}}


class ShardWriter:
    def __init__(self, out, chrom, source, config):
        self.out, self.chrom, self.source, self.config = Path(out), chrom, source, config
        self.paths = shard_paths(out, chrom)
        self.temporary = {}
        self.handles = {}
        self.writers = {}
        for key in ("sites", "controls", "windows", "patterns"):
            self.paths[key].parent.mkdir(parents=True, exist_ok=True)
            temp = data_path(self.out, 'tmp') / (self.paths[key].name + ".partial")
            temp.parent.mkdir(parents=True, exist_ok=True)
            self.temporary[key] = temp
            self.handles[key] = gzip.open(temp, "wt", encoding="utf-8", newline="", compresslevel=3)
            fields = WINDOW_FIELDS if key == "windows" else PATTERN_FIELDS if key == "patterns" else SITE_FIELDS
            writer = csv.DictWriter(self.handles[key], fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            self.writers[key] = writer
        self.window = None
        self.window_counts = defaultdict(Counter)
        self.totals = defaultdict(Counter)
        self.evidence_totals = defaultdict(Counter)
        self.pattern_counts = Counter()
        self.pattern_rows = 0
        self.rows = self.controls = self.sites = 0
        self.issues = Counter()

    def flush_window(self):
        if self.window is None:
            return
        for (p, kind, threshold), counts in sorted(self.window_counts.items()):
            self.writers["windows"].writerow(dict(chrom=self.chrom, start=self.window,
                end=self.window + self.config["window_size"], population=p, threshold_kind=kind,
                threshold=fmt(threshold), **{k: counts[k] for k in COUNT_FIELDS}))
        self.window_counts.clear()

    def add(self, group):
        if not group.observations:
            return
        current_window = ((group.pos - 1) // self.config["window_size"]) * self.config["window_size"]
        if self.window != current_window:
            self.flush_window()
            self.window = current_window
        self.sites += 1
        self.issues["duplicate_ALT_summaries"] += group.duplicate_values
        self.issues["conflicting_ALT_summaries"] += group.conflicts
        self.issues["REF_conflict_sites"] += int(group.ref_conflict)
        hv = stable_hash_fraction(group.chrom, group.pos, group.ref, self.config["seed"])
        populations = list(group.observations) if self.config["populations"] == ["all"] else self.config["populations"]
        primary_rows = []
        for p in populations:
            if p not in group.observations:
                continue
            prepared = prepare_population(group, p, self.config)
            kind = prepared["threshold_kind"]
            thresholds = self.config["ac_thresholds"] if kind == "AC" else self.config["af_thresholds"]
            primary = self.config["min_ac"] if kind == "AC" else self.config["min_af"]
            row = population_row(group, p, prepared, primary, self.config, hv)
            row.update(allele_evidence(row))
            primary_rows.append(row)
            add_evidence_counts(self.evidence_totals[p], row)
            serialized = serialize_row(row)
            self.writers["sites"].writerow(serialized)
            self.rows += 1
            if row["control"]:
                self.writers["controls"].writerow(serialized)
                self.controls += 1
            for threshold in thresholds:
                r = row if threshold == primary else population_row(group, p, prepared, threshold, self.config, hv)
                key = p, kind, threshold
                add_counts(self.window_counts[key], r)
                add_counts(self.totals[key], r)
        if primary_rows:
            pattern = classify_site(primary_rows, required_populations=populations)
            add_pattern_counts(self.pattern_counts, pattern)
            if pattern["candidate"]:
                self.writers["patterns"].writerow({k:fmt(pattern.get(k)) for k in PATTERN_FIELDS})
                self.pattern_rows += 1

    def close(self):
        self.flush_window()
        for h in self.handles.values():
            h.close()
        for key in self.temporary:
            atomic_replace(self.temporary[key], self.paths[key])
        totals = [dict(population=p, threshold_kind=k, threshold=t, **{c: v[c] for c in COUNT_FIELDS}) for (p, k, t), v in sorted(self.totals.items())]
        receipt = dict(version=VERSION, status="complete", chrom=self.chrom, source=self.source,
                       config=self.config, sites=self.sites, rows=self.rows, control_rows=self.controls,
                       evidence_schema=EVIDENCE_SCHEMA, site_columns=SITE_FIELDS,
                       population_pattern_schema=PATTERN_SCHEMA, population_pattern_rows=self.pattern_rows,
                       population_pattern_counts={k:self.pattern_counts[k] for k in PATTERN_COUNTS},
                       evidence_totals=[dict(population=p, **{k:n[k] for k in EVIDENCE_COUNTS})
                                        for p,n in sorted(self.evidence_totals.items())],
                       totals=totals, issues=dict(self.issues),
                       outputs={k: dict(path=str(self.paths[k].relative_to(self.out)), size=self.paths[k].stat().st_size) for k in self.temporary})
        write_json(self.paths["receipt"], receipt)
        return receipt

    def abort(self):
        for h in self.handles.values():
            h.close()
        for path in self.temporary.values():
            path.unlink(missing_ok=True)


def filename_chromosome(path):
    match = re.search(r"(?:^|[._-])(chr(?:[1-9]|1\d|2[0-2]|X|Y|M|MT))(?=[._-]|$)", Path(path).name, re.I)
    return chromosome(match.group(1)) if match else None


def scan_source(task):
    path, out, config, resume = task
    source = file_identity(path)
    selected = set(config["chromosomes"])
    skipped = {c for c in selected if resume and usable_receipt(out, c, source, config)}
    hint = filename_chromosome(path)
    if hint in skipped and config["format"] == "vcf":
        return dict(source=source, scanned=[], reused=sorted(skipped, key=chromosome_key), stats={})
    wanted = selected - skipped
    if not wanted:
        return dict(source=source, scanned=[], reused=sorted(skipped, key=chromosome_key), stats={})
    stats = Counter()
    records = vcf_records(path, config, stats, wanted) if config["format"] == "vcf" else tsv_records(path, config, stats, wanted)
    writer = None
    completed = []
    started = time.monotonic()
    try:
        for group in site_groups(records):
            if hint and config["format"] == "vcf" and group.chrom != hint:
                raise ValueError(f"Chromosome-named VCF {path} also contains {group.chrom}; use one chromosome per shard or a neutral filename for a sequential multi-chromosome VCF")
            if writer is None or writer.chrom != group.chrom:
                if writer is not None:
                    completed.append(writer.close()["chrom"])
                prior = usable_receipt(out, group.chrom)
                if prior and prior["source"]["path"] != source["path"]:
                    raise ValueError(f"Two input files claim {group.chrom}; supply one exhaustive file per chromosome")
                writer = ShardWriter(out, group.chrom, source, config)
            writer.add(group)
            if writer.sites and writer.sites % config["progress_every"] == 0:
                print(f"[summary-scan] {Path(path).name} {group.chrom}:{group.pos:,} sites={writer.sites:,} elapsed={time.monotonic()-started:.0f}s", file=sys.stderr, flush=True)
        if writer is not None:
            completed.append(writer.close()["chrom"])
            writer = None
    except BaseException:
        if writer is not None:
            writer.abort()
        raise
    return dict(source=source, scanned=completed, reused=sorted(skipped, key=chromosome_key), stats=dict(stats))


def write_table(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k)) for k in fields})


def svg_bars(rows, key):
    if not rows:
        return '<p>No qualifying chromosome summary is available.</p>'
    width, height, left, bottom = 920, 250, 65, 210
    maximum = max((r.get(key, 0) for r in rows), default=0) or 1
    gap = (width - left - 25) / len(rows)
    parts = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Multi ALT sites by chromosome">']
    for i, r in enumerate(rows):
        h = 175 * r[key] / maximum
        x = left + i * gap
        parts += [f'<rect x="{x:.1f}" y="{bottom-h:.1f}" width="{max(2,gap-5):.1f}" height="{h:.1f}" rx="2" fill="#277da1"><title>{html.escape(r["chrom"])}: {r[key]:,}</title></rect>',
                  f'<text x="{x+gap/2:.1f}" y="230" text-anchor="middle" font-size="11">{html.escape(r["chrom"][3:])}</text>']
    parts += [f'<text x="55" y="38" text-anchor="end" font-size="11">{maximum:,}</text>', '<line x1="65" y1="210" x2="895" y2="210" stroke="#9aa7b2"/>', '</svg>']
    return "".join(parts)


def make_report(out, manifest, chromosome_rows, totals, evidence_rows):
    out = Path(out)
    report = out / "reports" / "index.html"
    report.parent.mkdir(parents=True, exist_ok=True)
    primary = [r for r in chromosome_rows if r["population"] == "ALL" and r["threshold"] == manifest["config"]["min_ac" if r["threshold_kind"] == "AC" else "min_af"]]
    table = []
    for r in totals:
        table.append("<tr>" + "".join(f'<td>{html.escape(fmt(r[k]))}</td>' for k in ["population", "threshold_kind", "threshold", "sites_with_observed_alt", "multi_alt_sites", "triallelic_sites", "quadallelic_sites", "ref_unknown_sites"]) + "</tr>")
    source_stats = Counter()
    for result in manifest.get("source_runs", []):
        source_stats.update(result.get("stats", {}))
    evidence_fields = ["population","two_observed_alt_sites","three_observed_alt_sites",
        "ref_plus_observed_quadallelic_sites","proven_quadallelic_sites","unresolved_occupancy_sites","summary_status"]
    evidence_table = "".join("<tr>"+"".join("<td>"+html.escape(fmt(r.get(k)))+"</td>" for k in evidence_fields)+"</tr>" for r in evidence_rows)
    report.write_text('''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Population allele scan</title>
<style>body{font:16px/1.5 system-ui,sans-serif;max-width:1120px;margin:35px auto;padding:0 22px;color:#183148;background:#f6f8fa}h1,h2{color:#16324f}section{background:white;border-radius:12px;padding:22px;margin:18px 0;box-shadow:0 1px 5px #d9e0e7}table{border-collapse:collapse;font-size:13px;width:100%}td,th{padding:7px;border-bottom:1px solid #dce5eb;text-align:right}td:first-child,th:first-child{text-align:left}a{color:#12678d}code{background:#edf2f5;padding:2px 4px}svg{width:100%}.scroll{overflow-x:auto}.note{color:#526877}</style>
<h1>Population allele scan</h1>'''+f'<p>{html.escape(manifest["target"])} · GRCh38 · {len(manifest["shards"])} chromosomes · {manifest["site_rows"]:,} population rows</p>'+
        '<section><h2>What the scan measures</h2><p>Positive reported AF can identify two or three observed ALT bases and supports multi-ALT analysis. The supplied All of Us gvs_*_af file contains such frequency evidence; the cohort name does not restrict its capabilities. An explicit schema can supply actual AC/AN columns from a count-bearing release.</p><p>Multi ALT means at least two observed non-reference nucleotides above the stated threshold. Site shards additionally distinguish positive, reported-zero, missing and unreported ALT identities, record-defined REF+ALT bases, REF+positive-ALT candidates, and bounds on observed total bases. For T with A/C/G all at positive AF, three ALT bases are observed and REF+positive-ALT gives a four-base candidate; without REF evidence the total observed count is bounded at 3–4.</p><p>Exact singleton means AC=1. AF alone does not identify AC or AN. Strict bi-, tri- and quad-allelic states require joint REF evidence; equal-looking AN in independent split records does not establish one joint called sample set. Missing AF is not zero; reported AF=0 can reflect rounding. REF is a sequence label, not necessarily ancestral or most frequent. Thresholds apply within the stated population.</p><p class="note">Counts are descriptive. They do not by themselves demonstrate a mutation rate, balancing selection, an allele origin date or an enrichment over callable genomic sequence.</p></section>'+
        '<section><h2>Multi ALT sites by chromosome · ALL · primary threshold</h2>'+svg_bars(primary, "multi_alt_sites")+'</section>'+
        '<section><h2>Observed-base evidence · before the primary threshold</h2><div class="scroll"><table><thead><tr>'+''.join('<th>'+html.escape(s)+'</th>' for s in ["Population","Two observed ALTs","Three observed ALTs","REF+positive ALT quad candidate","Proven quad","Unresolved occupancy bound","Summary status"])+"</tr></thead><tbody>"+evidence_table+'</tbody></table></div><p>REF+positive ALT is an explicit reference-based candidate convention, not proof that REF was observed. Record-based counts can also include reported AF=0 alternatives. Old completed shards remain reusable; unavailable new aggregate statistics stay NA, while the streaming reader derives the new site-level evidence from their existing columns.</p></section>'+
        '<section><h2>Within-population coexistence versus pooled diversity</h2><p><a href="../results/population_patterns.tsv">Complete classification counts</a> · <a href="../results/population_candidates.tsv.gz">All primary-threshold candidate sites</a></p><p>ALL is pooled and is not counted as an independent ancestry group. A site can show multiple ALTs within a named ancestry group, or only in the union of groups whose individual summaries each contain at most one positive ALT. Missing observations prevent the latter conclusion. These patterns describe the available summaries and are not evidence that an allele arose recently or is maintained by balancing selection.</p><p>Use <code>--populations all</code> for the available ancestry groups. With only ALL, the structure is <code>pooled_only</code> and population-pattern inference remains unresolved, not negative. Candidate rows also retain pairwise AF-rank reversals across comparable groups. Source estimates may differ in detection power.</p><pre>'+html.escape(json.dumps(manifest["population_pattern_summary"], indent=2))+'</pre></section>'+
        '<section><h2>Threshold sensitivity</h2><div class="scroll"><table><thead><tr>'+''.join('<th>'+s+'</th>' for s in ["Population","Kind","Threshold","Observed sites","Multi ALT","Strict tri","Strict quad","REF unknown"])+"</tr></thead><tbody>"+"".join(table)+"</tbody></table></div></section>"+
        '<section><h2>Reusable results</h2><ul><li><a href="../results/window_summary.tsv.gz">Window summaries</a> (0-based half-open intervals; raw counts, without a callable-base denominator)</li><li><a href="../results/chromosome_summary.tsv">Chromosome counts</a></li><li><a href="../results/population_summary.tsv">Population and threshold counts</a></li><li><a href="../results/allele_evidence_summary.tsv">Observed-base evidence by population</a></li><li><a href="../results/source_capabilities.tsv">Capabilities of the effective input columns</a></li><li><a href="../results/site_catalog.tsv">Shard catalog</a></li></ul><p>Site shards retain exact nucleotide identities and aligned AF/AC/AN arrays. Control shards contain an independent deterministic hash sample of sites with one qualified ALT, with no annotation or functional-score selection.</p></section>'+
        '<section><h2>Input accounting</h2><pre>'+html.escape(json.dumps(dict(source_stats), indent=2))+ '</pre><p class="note">Source accounting is for files read in this invocation; reused chromosome receipts retain their original counters. An ordinary All of Us .idx is not assumed to be Tabix. A single compressed TSV is scanned sequentially once and split into chromosome outputs.</p></section></html>')


def merge_store(out, config=None, target=None, source_runs=None, sources=None):
    out = Path(out)
    receipts = []
    for receipt_path in (data_path(out, 'qc')).glob("chr*.json"):
        c = chromosome(receipt_path.stem)
        if c is None or (config is not None and c not in config["chromosomes"]):
            continue
        receipt = usable_receipt(out, c, config=config)
        if receipt:
            if sources is not None and receipt["source"] not in sources:
                continue
            receipts.append(receipt)
    receipts.sort(key=lambda r: chromosome_key(r["chrom"]))
    if not receipts:
        raise ValueError("No complete chromosome shards were found in " + str(out))
    if config is None:
        config = dict(receipts[0]["config"], chromosomes=[r["chrom"] for r in receipts])
    if any(scientific_config(r["config"]) != scientific_config(config) for r in receipts):
        raise ValueError("Cannot merge chromosome shards produced with different scientific parameters")
    sums = defaultdict(Counter)
    chrom_rows = []
    for r in receipts:
        for total in r["totals"]:
            chrom_rows.append(dict(chrom=r["chrom"], **total))
            key = total["population"], total["threshold_kind"], total["threshold"]
            sums[key].update({c: total[c] for c in COUNT_FIELDS})
    totals = [dict(population=p, threshold_kind=k, threshold=t, **{c: n[c] for c in COUNT_FIELDS}) for (p, k, t), n in sorted(sums.items())]
    results = out / "results"
    fields = ["population", "threshold_kind", "threshold", *COUNT_FIELDS]
    write_table(results / "chromosome_summary.tsv", chrom_rows, ["chrom", *fields])
    write_table(results / "population_summary.tsv", totals, fields)
    write_table(results / "occupancy_summary.tsv", totals, fields)
    write_table(results / "singleton_summary.tsv", totals, ["population", "threshold_kind", "threshold", "sites_with_observed_alt", "alt_singleton_sites", "singleton_unknown_sites"])
    evidence_rows = []
    populations = sorted({r["population"] for r in totals})
    for p in populations:
        summed = Counter(); measured = []; unavailable = []
        for receipt in receipts:
            if not any(x["population"] == p for x in receipt["totals"]):
                continue
            values = next((x for x in receipt.get("evidence_totals", []) if x["population"] == p), None)
            if values is None:
                unavailable.append(receipt["chrom"])
            else:
                summed.update({k: values[k] for k in EVIDENCE_COUNTS}); measured.append(receipt["chrom"])
        evidence_rows.append(dict(population=p, summary_status="not_available_from_prior_receipts" if unavailable else "complete",
            measured_chromosomes=",".join(measured), unavailable_chromosomes=",".join(unavailable),
            **{k: None if unavailable else summed[k] for k in EVIDENCE_COUNTS}))
    write_table(results / "allele_evidence_summary.tsv", evidence_rows,
                ["population","summary_status","measured_chromosomes","unavailable_chromosomes",*EVIDENCE_COUNTS])
    capabilities = []
    for p in populations:
        if config["format"] == "vcf":
            suffix = "" if p == "ALL" else "_"+p.lower()
            columns = {k: k.upper()+suffix for k in ("af","ac","an")}
        else:
            columns = config["schema"]["populations"][p]
        has_counts = "ac" in columns
        capabilities.append(dict(population=p, target=target or config.get("target","custom"),
            format=config["format"], af_column=columns.get("af"), ac_column=columns.get("ac"), an_column=columns.get("an"),
            observation_evidence="positive_filtered_AC" if has_counts else "positive_reported_AF",
            multi_alt_analysis="available", exact_singleton="requires_nonmissing_AC" if has_counts else "unavailable_without_AC",
            joint_ref_counts="requires_coherent_joint_counts" if (config["format"] == "vcf" or config.get("schema",{}).get("joint_counts")) else "not_established_by_split_rows",
            qualification=config["threshold_kind"] if config["threshold_kind"] != "auto" else "AC" if has_counts else "AF"))
    write_table(results / "source_capabilities.tsv", capabilities,
                ["population","target","format","af_column","ac_column","an_column","observation_evidence","multi_alt_analysis","exact_singleton","joint_ref_counts","qualification"])
    pattern_rows = []; pattern_sums = Counter(); unavailable_patterns = []
    for receipt in receipts:
        values = receipt.get("population_pattern_counts")
        available = values is not None and "patterns" in receipt["outputs"]
        if available:
            pattern_sums.update(values)
        else:
            unavailable_patterns.append(receipt["chrom"])
        pattern_rows.append(dict(chrom=receipt["chrom"], summary_status="complete" if available else "not_available_from_prior_receipts",
                                 **{k:values[k] if available else None for k in PATTERN_COUNTS}))
    pattern_status = "not_available_from_prior_receipts" if unavailable_patterns else "complete"
    pattern_rows.append(dict(chrom="all_chromosomes", summary_status=pattern_status,
                             **{k:None if unavailable_patterns else pattern_sums[k] for k in PATTERN_COUNTS}))
    write_table(results / "population_patterns.tsv", pattern_rows, ["chrom","summary_status",*PATTERN_COUNTS])
    def candidates():
        for receipt in receipts:
            if "patterns" not in receipt["outputs"]:
                continue
            with open_text(out / receipt["outputs"]["patterns"]["path"]) as stream:
                yield from csv.DictReader(stream, delimiter="\t")
    write_table(results / "population_candidates.tsv.gz", candidates(), PATTERN_FIELDS)
    def windows():
        for r in receipts:
            with open_text(shard_paths(out, r["chrom"])["windows"]) as f:
                yield from csv.DictReader(f, delimiter="\t")
    write_table(results / "window_summary.tsv.gz", windows(), WINDOW_FIELDS)
    catalog = [dict(chrom=r["chrom"], file=r["outputs"]["sites"]["path"], controls=r["outputs"]["controls"]["path"],
                    population_patterns=r["outputs"].get("patterns",{}).get("path"),
                    sites=r["sites"], rows=r["rows"], control_rows=r["control_rows"]) for r in receipts]
    write_table(results / "site_catalog.tsv", catalog, ["chrom", "file", "controls", "population_patterns", "sites", "rows", "control_rows"])
    old = {}
    try:
        old = json.loads((data_path(out, 'qc') / "manifest.json").read_text())
    except (OSError, ValueError):
        pass
    manifest = dict(version=VERSION, status="complete", build="GRCh38", target=target or config.get("target", old.get("target", "custom")),
                    config=config, shards=catalog, site_rows=sum(r["rows"] for r in receipts),
                    unique_sites=sum(r["sites"] for r in receipts), control_rows=sum(r["control_rows"] for r in receipts),
                    sources=sources if sources is not None else list({r["source"]["path"]: r["source"] for r in receipts}.values()),
                    source_runs=source_runs if source_runs is not None else old.get("source_runs", []),
                    invocation=old.get("invocation"),
                    window_coordinates="0-based-half-open", site_coordinates="1-based", populations=sorted({r["population"] for r in totals}),
                    site_columns=SITE_FIELDS if all(r.get("site_columns") == SITE_FIELDS for r in receipts) else BASE_SITE_FIELDS,
                    derived_site_columns=EVIDENCE_FIELDS, evidence_schema=EVIDENCE_SCHEMA,
                    population_pattern_schema=PATTERN_SCHEMA,
                    population_pattern_summary=dict(status=pattern_status, unavailable_chromosomes=unavailable_patterns,
                        populations="ALL is pooled; ancestry rows are considered separately", counts=pattern_rows[-1]),
                    absent_selected_chromosomes=sorted(set(config["chromosomes"]) - {r["chrom"] for r in receipts}, key=chromosome_key),
                    interpretation={"primary_af": "maximum known ALT AF within the named population; not major AF", "AF_only": "No inferred AC, AN, singleton, or REF count", "gnomad": "FILTER must equal PASS; filtered AC>0 establishes an observed ALT; AC_raw is never substituted", "AF_observation": "AF>0 establishes observation for AF-only summaries; zero is excluded, missing remains unknown", "split_records": "equal AN alone does not prove a common called sample set", "control_sampling": "blake2b on seed/chrom/pos/ref, independent of annotations and functional scores"})
    manifest["interpretation"].update(
        input_capabilities="Effective columns/schema determine AC/AN availability, not the target name; source_capabilities.tsv records the mapping",
        record_state="REF plus every reported SNV identity, including zero/missing AF; a record definition, not population occupancy",
        ref_plus_observed_alt_state="REF plus positive-evidence ALT bases; a conventional candidate definition that does not assert REF observation",
        occupancy_bounds="observed_alleles_min/max retain positive-base evidence while unresolved REF or unreported ALT identities remain uncertain; reported zero AF may reflect rounding",
        old_shards="Existing base-schema shards remain readable; new evidence fields are derived by iter_sites without rewriting or rescanning input")
    make_report(out, manifest, chrom_rows, totals, evidence_rows)
    write_json(data_path(out, 'qc') / "manifest.json", manifest)
    return manifest


def verify_store(out, read_all=False):
    out = Path(out)
    manifest = json.loads((data_path(out, 'qc') / "manifest.json").read_text())
    if manifest.get("status") != "complete" or manifest.get("version") != VERSION:
        raise ValueError("Incomplete or unsupported summary scan: " + str(out))
    required = ['results/'+name for name in (
        'chromosome_summary.tsv','population_summary.tsv','occupancy_summary.tsv',
        'singleton_summary.tsv','allele_evidence_summary.tsv','source_capabilities.tsv',
        'population_patterns.tsv','population_candidates.tsv.gz','window_summary.tsv.gz','site_catalog.tsv')]
    required += ['reports/index.html']
    for name in required:
        if not (out/name).is_file() or (out/name).stat().st_size == 0:
            raise ValueError('Missing merged summary product: '+str(out/name))
    for shard in manifest["shards"]:
        receipt = usable_receipt(out, shard["chrom"], config=manifest["config"])
        if not receipt:
            raise ValueError("Missing or incomplete shard: " + shard["chrom"])
        with open_text(out / shard["file"]) as f:
            if next(f, "").rstrip("\r\n").split("\t") != receipt.get("site_columns", BASE_SITE_FIELDS):
                raise ValueError("Unexpected shard schema: " + shard["file"])
            if read_all:
                count = sum(1 for _ in f)
                if count != shard["rows"]:
                    raise ValueError("Shard row count mismatch: " + shard["file"])
    return dict(status="verified", shards=len(manifest["shards"]), validation="completion receipts, sizes and headers" + ("; full gzip/row-count read" if read_all else ""))


def iter_sites(store, chromosomes=None, population="ALL", controls_only=False):
    """Yield parsed rows without loading a chromosome or summary table in RAM.

    population=None returns all stored populations. ``controls_only`` reads the
    sparse control shards directly. ALTs remain nucleotide-sorted, never sorted
    by frequency. No unknown count or frequency is replaced with zero.
    """
    root = Path(store)
    manifest = json.loads((data_path(root, 'qc') / "manifest.json").read_text())
    if manifest.get("status") != "complete" or manifest.get("version") != VERSION:
        raise ValueError("A completed unified summary scan is required: " + str(root))
    chroms = set(selected_chromosomes(chromosomes)) if chromosomes is not None else None
    for shard in sorted(manifest["shards"], key=lambda s: chromosome_key(s["chrom"])):
        if chroms is not None and shard["chrom"] not in chroms:
            continue
        from study_workflow import shard_path
        with open_text(shard_path(root, shard["controls" if controls_only else "file"])) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if population is not None and row["population"] != population:
                    continue
                row["alts"] = [] if row["alts"] in NA else row["alts"].split(",")
                for key in ("alt_afs", "alt_acs", "alt_ans"):
                    row[key] = [number(v, key != "alt_afs") for v in row[key].split(",")] if row["alts"] else []
                for key in ("pos", "n_alt_observed", "n_alt_qualified", "ref_ac", "AN", "control"):
                    row[key] = number(row[key], True)
                for key in ("ref_af", "primary_af"):
                    row[key] = number(row[key])
                row["threshold"] = float(row["threshold"])
                for key in ("observed_bases", "qualified_bases"):
                    row[key] = [] if row[key] == "" else (None if row[key] in NA else row[key].split(","))
                row["base_af"] = dict(zip(row["alts"], row["alt_afs"]))
                row["base_ac"] = dict(zip(row["alts"], row["alt_acs"]))
                row["base_an"] = dict(zip(row["alts"], row["alt_ans"]))
                row["base_af"][row["ref"]] = row["ref_af"]
                row["base_ac"][row["ref"]] = row["ref_ac"]
                row["base_an"][row["ref"]] = row["AN"]
                row.update(allele_evidence(row))
                yield row


def load_schema(path):
    schema = json.loads(Path(path).read_text())
    schema.setdefault("format", "tsv")
    if schema.get("build", "GRCh38") != "GRCh38":
        raise ValueError("Custom summaries must use GRCh38")
    if schema.get("coordinate_system", "1-based") not in {"1-based", "0-based-half-open"}:
        raise ValueError("Declare a 1-based or 0-based-half-open coordinate system")
    if schema.get("format", "tsv") not in {"tsv", "vcf"}:
        raise ValueError("Custom schema format must be tsv or vcf")
    if schema.get("format", "tsv") == "tsv":
        for key in ("chrom", "pos", "ref", "alt", "populations"):
            if key not in schema:
                raise ValueError("Custom schema lacks " + key)
        schema["populations"] = {p.upper(): spec for p, spec in schema["populations"].items()}
        for p, spec in schema["populations"].items():
            if not isinstance(spec, dict) or not (set(spec) <= {"ac", "an", "af"}) or not ("ac" in spec or "af" in spec):
                raise ValueError("Population schema needs af and/or ac, optionally an: " + p)
        if schema.get("joint_counts") and any("ac" not in s or "an" not in s for s in schema["populations"].values()):
            raise ValueError("joint_counts requires explicit AC and AN columns; AF alone is insufficient")
    return schema


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", nargs="?", default="scan", choices=["scan", "merge", "verify"])
    p.add_argument("--target", choices=["gnomad", "allofus", "custom"], default="gnomad")
    p.add_argument("--input", action="append", default=[], help="Source file; repeat for separate chromosome VCFs")
    p.add_argument("--vcf", action="append", default=[])
    p.add_argument("--vcf-dir", "--dir-vcf", dest="vcf_dir")
    p.add_argument("--schema", help="Explicit JSON column schema; overrides cohort defaults and can provide actual AC/AN columns")
    p.add_argument("--build", choices=["GRCh38", "38", "hg38"], default="GRCh38")
    p.add_argument("--out")
    p.add_argument("--chr", default="autosomes")
    p.add_argument("--populations", default="all", help="all available populations by default; ALL requests only the pooled stratum")
    p.add_argument("--min-ac", type=int, default=None)
    p.add_argument("--min-af", type=float, default=None)
    p.add_argument("--threshold-kind", choices=["auto", "AC", "AF"], default="auto", help="auto uses AC when counts exist; AF for AF-only sources")
    p.add_argument("--ac-thresholds", default=None)
    p.add_argument("--af-thresholds", default=None)
    p.add_argument("--window-size", type=int, default=1000000)
    p.add_argument("--control-rate", type=float, default=0.001)
    p.add_argument("--seed", type=int, default=20261115)
    p.add_argument("--jobs", type=int, default=1, help="Parallel workers for independent chromosome VCFs")
    p.add_argument("--replace", choices=["TRUE", "FALSE"], default="FALSE", help="TRUE recomputes selected shards; FALSE resumes complete matching shards")
    p.add_argument("--progress-every", type=int, default=1000000)
    p.add_argument("--read-all", action="store_true", help="verify only: read complete gzip streams instead of lightweight verification")
    p.add_argument("--no-merge", action="store_true", help="Only write chromosome shards/receipts; use merge after independent chromosome jobs")
    p.add_argument("--plan-only", action="store_true")
    return p


def run(args):
    root = Path(os.environ.get("ALLELE_ANALYSIS_ROOT", str(Path(os.environ.get("ANALYSIS_ROOT", "analysis")) / "allele3" / "38")))
    out = Path(args.out or Path(os.environ.get("ALLELE_SUMMARY_STORE", str(__import__('product_paths').defaults()['summary']))) / args.target).expanduser().resolve()
    if args.command == "verify":
        return verify_store(out, args.read_all)
    if args.command == "merge":
        result = merge_store(out)
        return dict(status=result["status"], out=str(out), shards=len(result["shards"]), site_rows=result["site_rows"])
    if args.window_size < 1 or args.jobs < 1 or args.progress_every < 1 or not 0 <= args.control_rate <= 1:
        raise ValueError("Window size, jobs and progress interval must be positive; control rate must be in [0,1]")
    inputs = args.input + args.vcf
    if not inputs and args.target == "allofus":
        annot = Path(os.environ.get("DIR_ANNOT", str(Path(os.environ.get("DIR_DATA", "data")) / "annot")))
        inputs = [os.environ.get("ALLELE_ALL_OF_US", str(annot / "allofus/38/hg38_allofus.txt.gz"))]
    schema = load_schema(args.schema) if args.schema else (discover_aou_schema(Path(inputs[0]).expanduser()) if args.target == "allofus" else {})
    if args.target == "custom" and not args.schema:
        raise ValueError("A custom summary requires --schema")
    form = schema.get("format", "vcf")
    requested = None if args.populations == "all" else {x.strip().upper() for x in args.populations.split(",")}
    specifications = [spec for p, spec in schema.get("populations", {}).items() if requested is None or p in requested]
    count_source = form == "vcf" or any("ac" in spec for spec in specifications)
    af_only_source = form == "tsv" and specifications and all("ac" not in spec for spec in specifications)
    if af_only_source and (args.min_ac is not None or args.ac_thresholds is not None or args.threshold_kind == "AC"):
        raise ValueError("The selected columns are AF-only: use --min-af/--af-thresholds, or provide actual AC/AN columns with --schema")
    if count_source and args.min_af is not None and args.threshold_kind != "AF":
        raise ValueError("For an AF-threshold analysis with count-bearing columns also select --threshold-kind AF; the default is a separate AC-threshold analysis")
    min_ac, min_af = args.min_ac if args.min_ac is not None else 1, args.min_af if args.min_af is not None else 0.0
    if min_ac < 1 or not 0 <= min_af <= 1:
        raise ValueError("--min-ac must be >=1; --min-af must be in [0,1]")
    ac_thresholds = sorted(set([int(x) for x in (args.ac_thresholds or "1,2,5,10").split(",")] + [min_ac]))
    af_thresholds = sorted(set([float(x) for x in (args.af_thresholds or "0,0.00001,0.0001,0.001").split(",")] + [min_af]))
    if min(ac_thresholds) < 1 or any(not math.isfinite(t) or not 0 <= t <= 1 for t in af_thresholds):
        raise ValueError("Invalid AC or AF sensitivity thresholds")
    populations = ["all"] if args.populations == "all" else list(dict.fromkeys(x.strip().upper() for x in args.populations.split(",")))
    if not populations or any(not x for x in populations):
        raise ValueError("--populations cannot be empty")
    config = dict(target=args.target, format=form, schema=schema, chromosomes=selected_chromosomes(args.chr), populations=populations, threshold_kind=args.threshold_kind,
                  min_ac=min_ac, min_af=min_af, ac_thresholds=ac_thresholds, af_thresholds=af_thresholds,
                  window_size=args.window_size, control_rate=args.control_rate, seed=args.seed, progress_every=args.progress_every)
    if not inputs:
        annot = Path(os.environ.get("DIR_ANNOT", str(Path(os.environ.get("DIR_DATA", "data")) / "annot")))
        if args.target == "gnomad":
            directory = Path(args.vcf_dir or os.environ.get("ALLELE_VCF_DIR", str(annot / "gnomad/38/sites")))
            inputs = [str(p) for pattern in ("chr*.vcf.bgz", "chr*.vcf.gz", "chr*.vcf") for p in directory.glob(pattern)]
            selected = set(config["chromosomes"])
            inputs = [x for x in inputs if filename_chromosome(x) is None or filename_chromosome(x) in selected]
        elif args.target == "allofus":
            inputs = [os.environ.get("ALLELE_ALL_OF_US", str(annot / "allofus/38/hg38_allofus.txt.gz"))]
    if not inputs:
        raise ValueError("No source files found; supply --input or --vcf-dir")
    inputs = sorted({str(Path(p).expanduser().resolve()) for p in inputs}, key=lambda p: (chromosome_key(filename_chromosome(p)) if filename_chromosome(p) else 99, p))
    hints = [filename_chromosome(p) for p in inputs]
    claimed = [c for c in hints if c is not None]
    if len(claimed) != len(set(claimed)):
        raise ValueError("Multiple source files match the same chromosome; select the intended release explicitly")
    for path in inputs:
        if not Path(path).is_file() or not Path(path).stat().st_size:
            raise ValueError("Missing or empty input: " + path)
        if out == Path(path) or out in Path(path).parents:
            raise ValueError("The reusable output store must not contain its source input: " + path)
    if args.plan_only:
        return dict(status="planned", out=str(out), target=args.target, inputs=inputs, config=config,
                    jobs=args.jobs if len(inputs)>1 and all(hints) else 1, ordinary_idx_used=False)
    from output_layout import prepare
    prepare(out)
    out.mkdir(parents=True, exist_ok=True)
    sources = [file_identity(p) for p in inputs]
    previous = {}
    try:
        previous = json.loads((data_path(out, 'qc') / "manifest.json").read_text())
    except (OSError, ValueError):
        pass
    if not args.no_merge and args.replace == "FALSE" and previous.get("config") == config and previous.get("sources") == sources:
        try:
            verified = verify_store(out)
            if previous.get("evidence_schema") != EVIDENCE_SCHEMA or previous.get("population_pattern_schema") != PATTERN_SCHEMA:
                merge_store(out)  # Refresh compact reports without rereading source variants.
            return dict(verified, status="reused", out=str(out), source_reads=0)
        except (OSError, ValueError, KeyError):
            pass
    if not args.no_merge:
        write_json(data_path(out, 'qc') / "manifest.json", dict(version=VERSION, status="running", target=args.target, config=config, sources=sources, invocation=vars(args)))
    tasks = [(path, str(out), config, args.replace == "FALSE") for path in inputs]
    if args.jobs > 1 and len(tasks) > 1 and all(hints) and form == "vcf":
        with ProcessPoolExecutor(max_workers=min(args.jobs, len(tasks))) as pool:
            runs = list(pool.map(scan_source, tasks))
    else:
        runs = [scan_source(task) for task in tasks]
    if args.no_merge:
        return dict(status="shards_complete", out=str(out), source_runs=runs, next_step="merge this store after all chromosome jobs finish")
    manifest = merge_store(out, config, args.target, runs, sources)
    result = dict(status="complete", out=str(out), shards=len(manifest["shards"]), unique_sites=manifest["unique_sites"], site_rows=manifest["site_rows"], control_rows=manifest["control_rows"])
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    result = run(args)
    if result.get('status') in {'complete','reused'} and not args.plan_only:
        from product_paths import defaults
        write_inventory(args.out or defaults()['summary']/args.target)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        raise SystemExit("summary-scan: " + str(exc))
