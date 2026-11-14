#!/usr/bin/env python3
"""Audit and compare an explicitly described All of Us site-summary table.

The default mode previews a bounded prefix only. ``targeted`` makes one full
sequential pass, retaining only requested coordinates. ``normalize`` explicitly
streams all primary-chromosome SNVs and requires chromosome-contiguous,
position-sorted input. An ANNOVAR .idx is never treated as a tabix index.

No sample-size assumption or AF-to-AC conversion is used. All ancestry metrics
are selected by exact column names and remain separate. These are descriptive
summary-file comparisons, not an independent enrichment or selection test.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import math
from pathlib import Path
import re

from analysis_integrity import identity, transactional


BASES = frozenset("ACGT")
PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
MISSING = {"", ".", "na", "nan", "none", "null", "n/a"}
BUILDS = {"37": "GRCh37", "hg19": "GRCh37", "grch37": "GRCh37",
          "38": "GRCh38", "hg38": "GRCh38", "grch38": "GRCh38"}
ANNOVAR_COLUMNS = ["Chr", "Start", "End", "Ref", "Alt", "gvs_all_af", "gvs_max_af",
                   "gvs_afr_af", "gvs_amr_af", "gvs_eas_af", "gvs_eur_af", "gvs_mid_af",
                   "gvs_oth_af", "gvs_sas_af"]
ANNOVAR_POPULATIONS = {name: {"af": "gvs_" + name.lower() + "_af"}
                      for name in ("ALL", "AFR", "AMR", "EAS", "EUR", "MID", "OTH", "SAS")}
ALLELE_FIELDS = ["build", "chrom", "pos", "ref", "alt", "population", "AC", "AN", "AF",
                 "ref_AC", "af_source", "evidence", "count_supported", "min_ac", "reported_max_population_AF"]
SITE_FIELDS = ["build", "chrom", "pos", "ref", "population", "reported_alts", "n_reported_alts",
               "AN", "ref_AC", "n_alts_with_usable_counts", "n_alts_ac_ge_min",
               "count_supported_multi_alt", "n_alts_positive_AF", "frequency_supported_multi_alt",
               "reported_AF_sum", "reported_AF_sum_exceeds_one",
               "n_supported_bases_lower_bound", "at_least_three_bases_demonstrated",
               "evidence", "min_ac", "scope"]
TARGET_ALLELE_FIELDS = ["site_id", "family", "group", "target_alt", "match_status", "source_ref"] + ALLELE_FIELDS
TARGET_SITE_FIELDS = ["site_id", "family", "group", "target_alts", "matched_target_alts",
                      "unreported_target_alts", "additional_reported_alts", "target_all_alts_count_supported",
                      "match_status", "source_ref"] + SITE_FIELDS


def apply_preset(args):
    """Use the user's actual ANNOVAR schema, preserving tab-separated blanks."""
    if getattr(args, "preset", None) is None:
        if not args.build:
            raise ValueError("An explicit --build or the verified ANNOVAR preset is required")
        return args
    if args.build and build_name(args.build) != "GRCh38":
        raise ValueError("annovar-allofus-hg38 requires GRCh38; it cannot relabel GRCh37")
    if args.delimiter != "tab" or args.coordinate_system not in {None, "1-based"}:
        raise ValueError("ANNOVAR preset requires tab fields and 1-based inclusive coordinates")
    if args.population or any((args.af_column, args.ac_column, args.an_column, args.ref_ac_column)):
        raise ValueError("ANNOVAR preset defines eight AF-only populations; do not override metric columns")
    if args.population_map and json.loads(args.population_map) != ANNOVAR_POPULATIONS:
        raise ValueError("ANNOVAR preset population mapping conflicts with the supplied mapping")
    for key, ordinary, exact in (("chrom_column", "chrom", "Chr"), ("pos_column", "pos", "Start"),
                                 ("ref_column", "ref", "Ref"), ("alt_column", "alt", "Alt"),
                                 ("end_column", None, "End")):
        if getattr(args, key) not in {ordinary, exact}:
            raise ValueError("ANNOVAR preset conflicts with --" + key.replace("_", "-"))
        setattr(args, key, exact)
    if args.columns and [x.strip().lstrip("#") for x in args.columns.split(",")] != ANNOVAR_COLUMNS:
        raise ValueError("ANNOVAR preset headerless order must exactly match the supplied 14-column schema")
    args.build, args.coordinate_system = "GRCh38", "1-based"
    args.population_map = json.dumps(ANNOVAR_POPULATIONS)
    args.max_af_column = "gvs_max_af"
    return args


def build_name(value):
    try:
        return BUILDS[str(value).lower()]
    except KeyError as exc:
        raise ValueError("Build must be explicitly 37/GRCh37/hg19 or 38/GRCh38/hg38") from exc


def chrom_name(value):
    c = str(value).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    c = c.upper()
    if c in {"M", "MT"}:
        return "MT"
    if c.isdigit():
        return str(int(c))
    return c


def open_text(path):
    # Inspect only the two-byte magic; gzip and BGZF can both stream here.
    path = Path(path)
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="utf-8-sig") if compressed else path.open(encoding="utf-8-sig")


def integer(raw, label, optional=False):
    raw = str(raw).strip()
    if optional and raw.lower() in MISSING:
        return None
    if raw.isdigit():
        return int(raw)
    try:
        value = Decimal(raw)
        if not value.is_finite() or value < 0 or value != value.to_integral_value():
            raise ValueError()
        return int(value)
    except (ValueError, InvalidOperation) as exc:
        raise ValueError(f"{label} must be a finite nonnegative integer: {raw!r}") from exc


def frequency(raw, label):
    if str(raw).strip().lower() in MISSING:
        return None
    try:
        value = float(raw)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid {label}: {raw!r}") from exc
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{label} must lie between zero and one: {raw!r}")
    return value


def split_line(line, delimiter):
    line = line.rstrip("\r\n")
    return line.split("\t") if delimiter == "tab" else line.split()


def populations(args):
    if args.population_map:
        if args.population or any((args.af_column, args.ac_column, args.an_column, args.ref_ac_column)):
            raise ValueError("Use --population-map OR --population and metric-column options")
        result = json.loads(args.population_map)
        if not isinstance(result, dict) or not result:
            raise ValueError("--population-map must be a nonempty JSON object")
    else:
        if not args.population:
            raise ValueError("Analysis requires --population and exact metric columns, or --population-map")
        result = {args.population: {k: v for k, v in dict(
            af=args.af_column, ac=args.ac_column, an=args.an_column, ref_ac=args.ref_ac_column).items() if v}}
    for name, mapping in result.items():
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or not isinstance(mapping, dict):
            raise ValueError("Invalid population label or mapping")
        if set(mapping) - {"af", "ac", "an", "ref_ac"} or any(not isinstance(v, str) or not v for v in mapping.values()):
            raise ValueError("Population mappings use only af/ac/an/ref_ac with exact column names")
        if bool(mapping.get("ac")) != bool(mapping.get("an")):
            raise ValueError("Count evidence requires both AC and AN columns")
        if not mapping.get("af") and not mapping.get("ac"):
            raise ValueError("Each population requires AF or AC/AN")
        if mapping.get("ref_ac") and not mapping.get("an"):
            raise ValueError("An explicit REF count requires the AC/AN mapping too")
        if len(mapping.values()) != len(set(mapping.values())):
            raise ValueError("Different metrics cannot refer to the same source column")
    return result


class Table:
    """Read a declared header or explicit headerless column sequence."""
    def __init__(self, path, args, required):
        self.path, self.args, self.required = path, args, set(required)
        self.columns = None

    def rows(self):
        explicit = [x.strip().lstrip("#") for x in self.args.columns.split(",")] if self.args.columns else None
        if explicit:
            self.validate_columns(explicit)
            self.columns = explicit
        with open_text(self.path) as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("##"):
                    continue
                comment = line.startswith("#")
                values = split_line(line.lstrip("#") if comment else line, self.args.delimiter)
                if self.columns is None:
                    if comment and not self.required <= set(values):
                        continue
                    if not self.required <= set(values):
                        raise ValueError("No usable source header. Preview first, then supply --columns and exact column mappings")
                    self.validate_columns(values)
                    self.columns = values
                    continue
                if comment:
                    continue
                if values == self.columns:
                    continue
                if len(values) != len(self.columns):
                    raise ValueError(f"Source line {number}: {len(values)} fields, expected {len(self.columns)}")
                yield number, values
        if self.columns is None:
            raise ValueError("No source header or data found")

    def validate_columns(self, names):
        if not names or any(not x for x in names) or len(names) != len(set(names)):
            raise ValueError("Column names must be nonempty and unique")
        missing = self.required - set(names)
        if missing:
            raise ValueError("Missing mapped source columns: " + ",".join(sorted(missing)))


def vector(value, n, label):
    if str(value).strip().lower() in MISSING:
        return [None] * n
    parts = str(value).split(",")
    if len(parts) != n:
        raise ValueError(f"{label}: require one value per ALT, in the supplied ALT order")
    return parts


@dataclass(frozen=True)
class Metric:
    ac: int | None
    an: int | None
    af: float | None
    ref_ac: int | None
    af_source: str = field(compare=False)
    reported_max_af: float | None = None

    @property
    def counted(self):
        return self.ac is not None and self.an is not None and self.an > 0

    @property
    def evidence(self):
        if self.counted:
            return "reported_AC_AN"
        if self.an == 0:
            return "no_callable_alleles"
        return "frequency_only" if self.af is not None else "metrics_missing"


def decode_metrics(row, mappings, n_alts, tolerance, max_af_column=None):
    result = [{} for _ in range(n_alts)]
    max_afs = vector(row[max_af_column], n_alts, "maximum population AF") if max_af_column else [None] * n_alts
    max_afs = [frequency(value, "maximum population AF") if value is not None else None for value in max_afs]
    for population, mapping in mappings.items():
        an = integer(row[mapping["an"]], "AN", True) if mapping.get("an") else None
        ref_ac = integer(row[mapping["ref_ac"]], "REF_AC", True) if mapping.get("ref_ac") else None
        acs = vector(row[mapping["ac"]], n_alts, "AC") if mapping.get("ac") else [None] * n_alts
        afs = vector(row[mapping["af"]], n_alts, "AF") if mapping.get("af") else [None] * n_alts
        if ref_ac is not None and an is not None and ref_ac > an:
            raise ValueError("REF_AC exceeds AN")
        for i in range(n_alts):
            ac = integer(acs[i], "AC", True) if acs[i] is not None else None
            af = frequency(afs[i], "AF") if afs[i] is not None else None
            af_source = "reported" if af is not None else "missing"
            if ac is not None and an is not None:
                if ac > an or (ref_ac is not None and ac + ref_ac > an):
                    raise ValueError("Allele counts exceed AN")
                if an == 0 and af not in {None, 0.0}:
                    raise ValueError("Positive AF with AN=0")
                if an > 0:
                    if af is not None and abs(af - ac / an) > tolerance:
                        raise ValueError("Reported AF disagrees with AC/AN; check mapping or explicitly set rounding tolerance")
                    if af is None:
                        af, af_source = ac / an, "derived_from_reported_AC_AN"
            elif an == 0 and af not in {None, 0.0}:
                raise ValueError("Positive AF with AN=0")
            result[i][population] = Metric(ac, an, af, ref_ac, af_source, max_afs[i])
        # Check the original vector before excluding indels or other ALT
        # identities. Missing AC contributes no known lower bound, and the
        # explicit REF count is added once, not once per ALT.
        known_alt_sum = sum(m[population].ac for m in result if m[population].ac is not None)
        if an is not None and known_alt_sum + (ref_ac if ref_ac is not None else 0) > an:
            raise ValueError("Sum of original ALT counts plus explicit REF count exceeds AN")
    return result


class Site:
    def __init__(self, chrom, pos, ref):
        self.chrom, self.pos, self.ref = chrom, pos, ref
        self.alleles = {}

    def add(self, ref, alt, metrics):
        if ref != self.ref:
            raise ValueError(f"Conflicting REF at {self.chrom}:{self.pos}; no pooling across REF alleles")
        if alt in self.alleles:
            if metrics != self.alleles[alt]:
                raise ValueError(f"Conflicting duplicate ALT metrics at {self.chrom}:{self.pos}:{ref}>{alt}")
            return False
        self.alleles[alt] = metrics
        return True

    def validate(self, population, tolerance):
        mm = [values[population] for values in self.alleles.values()]
        ans = {m.an for m in mm if m.an is not None}
        if len(ans) > 1:
            raise ValueError(f"Inconsistent AN across ALTs at {self.chrom}:{self.pos}, population {population}")
        an = next(iter(ans), None)
        # A missing AN does not license pooling that row's AC/REF_AC with
        # another row's denominator, even when both have the same pop label.
        comparable = [m for m in mm if m.an is not None and m.an == an]
        refs = {m.ref_ac for m in comparable if m.ref_ac is not None}
        if len(refs) > 1:
            raise ValueError(f"Inconsistent explicit REF count at {self.chrom}:{self.pos}, population {population}")
        ref_ac = next(iter(refs), None)
        if an is not None and sum(m.ac or 0 for m in comparable) + (ref_ac or 0) > an:
            raise ValueError(f"Sum of distinct reported allele counts exceeds AN at {self.chrom}:{self.pos}, {population}")
        return mm, an, ref_ac


def allele_row(site, alt, population, build, minimum):
    m = site.alleles[alt][population]
    return dict(build=build, chrom=site.chrom, pos=site.pos, ref=site.ref, alt=alt,
                population=population, AC=m.ac, AN=m.an, AF=m.af, ref_AC=m.ref_ac,
                af_source=m.af_source, evidence=m.evidence,
                count_supported=int(m.ac >= minimum) if m.counted else None, min_ac=minimum,
                reported_max_population_AF=m.reported_max_af)


def site_row(site, population, build, minimum, tolerance):
    mm, an, ref_ac = site.validate(population, tolerance)
    counted = [m for m in mm if m.counted]
    supported = sum(m.ac >= minimum for m in counted)
    af_known = [m for m in mm if m.af is not None and m.an != 0]
    positive = sum(m.af > 0 for m in af_known)
    lower = supported + int(ref_ac is not None and an is not None and an > 0 and ref_ac >= minimum)
    evidence = "reported_AC_AN" if len(counted) == len(mm) else "mixed_counts_and_missing_metrics" if counted else "frequency_only" if af_known else "metrics_unavailable"
    return dict(build=build, chrom=site.chrom, pos=site.pos, ref=site.ref, population=population,
                reported_alts=",".join(sorted(site.alleles)), n_reported_alts=len(mm), AN=an, ref_AC=ref_ac,
                n_alts_with_usable_counts=len(counted), n_alts_ac_ge_min=supported,
                count_supported_multi_alt=1 if supported >= 2 else 0 if len(counted) == len(mm) else None,
                n_alts_positive_AF=positive,
                frequency_supported_multi_alt=1 if positive >= 2 else 0 if len(af_known) == len(mm) else None,
                reported_AF_sum=sum(m.af for m in af_known) if af_known else None,
                reported_AF_sum_exceeds_one=int(sum(m.af for m in af_known) > 1 + len(mm) * tolerance) if af_known else None,
                n_supported_bases_lower_bound=lower if counted else None,
                at_least_three_bases_demonstrated=int(lower >= 3) if counted else None,
                evidence=evidence, min_ac=minimum,
                scope="reported_summary_alleles_only; unreported_not_zero; REF_count_never_inferred")


def read_targets(path, declared_build, source_build):
    if build_name(declared_build) != source_build:
        raise ValueError("Target build differs from source build; no liftover is performed")
    targets = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "pos", "ref", "alts"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Targets require chrom,pos,ref,alts; positions are 1-based")
        for row in reader:
            c, p, ref = chrom_name(row["chrom"]), integer(row["pos"], "target position"), row["ref"].upper()
            alts = tuple(row["alts"].upper().split(","))
            if c not in PRIMARY or p < 1 or ref not in BASES or not alts or len(alts) != len(set(alts)) or not set(alts) <= BASES - {ref}:
                raise ValueError("Targets must specify primary-chromosome 1-based SNVs with distinct DNA ALT bases")
            for name in ("build", "genome_build", "assembly", "target_build"):
                if row.get(name) and build_name(row[name]) != source_build:
                    raise ValueError("Target row build conflicts with the declared build")
            key = c, p
            normalized = dict(row, chrom=c, pos=p, ref=ref, alts=",".join(sorted(alts)),
                              site_id=row.get("site_id") or f"{c}:{p}:{ref}")
            if key in targets and targets[key] != normalized:
                raise ValueError(f"Conflicting/ambiguous duplicate target coordinate: {c}:{p}")
            targets[key] = normalized
    if not targets:
        raise ValueError("Empty targets table")
    return targets


def sort_key(key):
    c, p = key
    return (int(c) if c.isdigit() else {"X": 23, "Y": 24, "MT": 25}.get(c, 26), p)


class Writer:
    def __init__(self, path, fields):
        self.handle = gzip.open(path, "wt", newline="", compresslevel=3) if str(path).endswith(".gz") else open(path, "w", newline="")
        self.fields = fields
        self.writer = csv.DictWriter(self.handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        self.writer.writeheader()

    def write(self, row):
        self.writer.writerow({key: "NA" if row.get(key) is None else row.get(key, "NA") for key in self.fields})

    def close(self):
        self.handle.close()


def preview(args, out):
    with open_text(args.source) as handle:
        prefix = handle.read(args.preview_chars + 1)
    bounded = prefix[:args.preview_chars]
    lines = bounded.splitlines()[:args.preview_lines]
    records = [{"line": i + 1, "field_count": len(split_line(line, args.delimiter)), "text": line}
               for i, line in enumerate(lines)]
    data = dict(status="bounded_preview_only", declared_build=build_name(args.build),
                input=identity(args.source), delimiter=args.delimiter, lines=records,
                max_chars=args.preview_chars, max_lines=args.preview_lines,
                prefix_truncated=len(prefix) > args.preview_chars,
                prefix_sha256=hashlib.sha256(bounded.encode()).hexdigest(),
                schema="Inspect these fields; analysis never guesses AF/AC/AN or population columns",
                index_policy="Sequential reads only; an ANNOVAR .idx is not used as tabix",
                counts_policy="AF alone never supplies AC/AN; missing variants are not zero")
    (out / "preview.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    print("[allofus-preview] bounded prefix read; schema evidence is in preview.json", flush=True)
    return data


def target_outputs(out, targets, sites, mappings, args, build):
    aw, sw = Writer(out / "target_alleles.tsv", TARGET_ALLELE_FIELDS), Writer(out / "target_sites.tsv", TARGET_SITE_FIELDS)
    matches = Counter()
    try:
        for key in sorted(targets, key=sort_key):
            target, site = targets[key], sites.get(key)
            expected = target["alts"].split(",")
            correct_ref = site is not None and site.ref == target["ref"]
            for population in mappings:
                common = {k: target.get(k, ".") for k in ("site_id", "family", "group")}
                common.update(build=build, chrom=key[0], pos=key[1], ref=target["ref"], population=population,
                              min_ac=args.min_ac, source_ref=site.ref if site else None)
                status = "coordinate_matched" if correct_ref else "reference_mismatch" if site else "not_in_summary"
                matched = [a for a in expected if correct_ref and a in site.alleles]
                unmatched = sorted(set(expected) - set(matched))
                summary = site_row(site, population, build, args.min_ac, args.af_tolerance) if correct_ref else {}
                metrics = [site.alleles[a][population] for a in matched] if correct_ref else []
                target_supported = int(all(m.ac >= args.min_ac for m in metrics)) if len(matched) == len(expected) and all(m.counted for m in metrics) else None
                sw.write(dict(common, **{k: v for k, v in summary.items() if k not in common},
                              target_alts=target["alts"], matched_target_alts=",".join(matched) or ".",
                              unreported_target_alts=",".join(unmatched) or ".",
                              additional_reported_alts=",".join(sorted(set(site.alleles) - set(expected))) if correct_ref else None,
                              target_all_alts_count_supported=target_supported, match_status=status))
                for alt in expected:
                    data = allele_row(site, alt, population, build, args.min_ac) if alt in matched else {}
                    match_status = "matched_" + data["evidence"] if data else "reference_mismatch" if site and not correct_ref else "not_in_summary"
                    matches[match_status] += 1
                    aw.write(common | data | dict(alt=alt, target_alt=alt, match_status=match_status))
    finally:
        aw.close(); sw.close()
    return dict(matches)


def analyze(args, out, build):
    mappings = populations(args)
    targets = read_targets(args.targets, args.targets_build, build) if args.mode == "targeted" else None
    coords = {"chrom": args.chrom_column, "pos": args.pos_column, "ref": args.ref_column, "alt": args.alt_column}
    if args.end_column:
        coords["end"] = args.end_column
    if args.build_column:
        coords["build"] = args.build_column
    required = set(coords.values()) | {name for mapping in mappings.values() for name in mapping.values()}
    if getattr(args, "max_af_column", None):
        required.add(args.max_af_column)
    table = Table(args.source, args, required)
    counters, sites, completed_chroms = Counter(), {}, set()
    current = None
    aw = Writer(out / "allofus_alleles.tsv.gz", ALLELE_FIELDS)
    sw = Writer(out / "allofus_sites.tsv.gz", SITE_FIELDS)

    def emit(site):
        if site is None:
            return
        for pop in mappings:
            summary = site_row(site, pop, build, args.min_ac, args.af_tolerance)
            sw.write(summary)
            counters["site_population_pairs"] += 1
            counters["count_supported_multi_alt_pairs"] += summary["count_supported_multi_alt"] == 1
            counters["frequency_supported_multi_alt_pairs"] += summary["frequency_supported_multi_alt"] == 1
            counters["AF_sum_exceeds_one_pairs"] += summary["reported_AF_sum_exceeds_one"] == 1
            for alt in sorted(site.alleles):
                aw.write(allele_row(site, alt, pop, build, args.min_ac))

    try:
        indices = None
        for line_number, values in table.rows():
            counters["source_rows"] += 1
            if args.progress_every and counters["source_rows"] % args.progress_every == 0:
                print(f"[allofus] scanned {counters['source_rows']:,} rows; retained {counters['unique_snv_alleles']:,} SNV alleles", flush=True)
            if indices is None:
                indices = {name: table.columns.index(name) for name in required}
            if args.build_column and build_name(values[indices[args.build_column]]) != build:
                raise ValueError(f"Source line {line_number}: row build conflicts with declared build")
            chrom = chrom_name(values[indices[args.chrom_column]])
            if chrom not in PRIMARY:
                counters["nonprimary_rows"] += 1
                continue
            raw_pos = integer(values[indices[args.pos_column]], f"position on line {line_number}")
            pos = raw_pos + int(args.coordinate_system == "0-based")
            if pos < 1:
                raise ValueError(f"Source line {line_number}: invalid coordinate")
            key = chrom, pos
            if targets is not None and key not in targets:
                continue
            counters["selected_coordinate_rows"] += 1
            ref = values[indices[args.ref_column]].upper()
            alts = values[indices[args.alt_column]].upper().split(",")
            if ref not in BASES or not any(a in BASES for a in alts):
                counters["selected_non_snv_rows"] += 1
                continue
            if len(alts) != len(set(alts)) or ref in alts:
                raise ValueError(f"Source line {line_number}: REF/ALT identity or duplicate ALT")
            if any(a not in BASES for a in alts):
                counters["selected_mixed_alt_rows"] += 1
                counters["excluded_non_snv_alts"] += sum(a not in BASES for a in alts)
            if args.end_column:
                end = integer(values[indices[args.end_column]], "end coordinate")
                if end != pos:
                    raise ValueError(f"Source line {line_number}: SNV end is not the expected one-base interval")
            row = {name: values[index] for name, index in indices.items()}
            decoded = decode_metrics(row, mappings, len(alts), args.af_tolerance, getattr(args, "max_af_column", None))
            if targets is not None:
                site = sites.setdefault(key, Site(chrom, pos, ref))
            else:
                if current is None or (current.chrom, current.pos) != key:
                    if current is not None:
                        if chrom == current.chrom and pos < current.pos:
                            raise ValueError("normalize requires position-sorted input")
                        if chrom != current.chrom:
                            completed_chroms.add(current.chrom)
                            if chrom in completed_chroms:
                                raise ValueError("normalize requires each chromosome to form one contiguous block")
                        emit(current)
                    current = Site(chrom, pos, ref)
                site = current
            for alt, metrics in zip(alts, decoded):
                if alt not in BASES:
                    continue
                added = site.add(ref, alt, metrics)
                counters["unique_snv_alleles" if added else "identical_duplicate_rows"] += 1
        if targets is None:
            emit(current)
        else:
            for key in sorted(sites, key=sort_key):
                emit(sites[key])
    finally:
        aw.close(); sw.close()
    if not counters["source_rows"]:
        raise ValueError("Source contains no data rows")
    match_counts = target_outputs(out, targets, sites, mappings, args, build) if targets is not None else None
    return dict(mode=args.mode, counts=dict(counters), target_match_counts=match_counts,
                n_targets=len(targets) if targets is not None else None,
                raw_output_scope="source SNVs at selected coordinates; target_*.tsv applies exact REF/ALT matches" if targets is not None else "all streamed primary-chromosome SNVs",
                target_reference_checks={
                    "matched_coordinates": sum(key in sites and sites[key].ref == target["ref"] for key, target in targets.items()),
                    "reference_mismatch_coordinates": sum(key in sites and sites[key].ref != target["ref"] for key, target in targets.items()),
                    "coordinates_not_in_summary": sum(key not in sites for key in targets),
                } if targets is not None else None,
                schema=dict(columns=table.columns, coordinates=coords, populations=mappings,
                            delimiter=args.delimiter, coordinate_system=args.coordinate_system))


@transactional("allofus_manifest.json", input_files=("source", "targets"))
def run(args):
    apply_preset(args)
    build = build_name(args.build)
    if args.min_ac < 1 or not math.isfinite(args.af_tolerance) or args.af_tolerance < 0:
        raise ValueError("min-ac must be positive and AF tolerance finite and nonnegative")
    if not 1 <= args.preview_chars <= 1048576 or not 1 <= args.preview_lines <= 1000 or args.progress_every < 0:
        raise ValueError("Invalid bounded preview or progress settings")
    if args.mode != "preview" and args.coordinate_system is None:
        raise ValueError("Analysis requires explicit --coordinate-system 1-based or 0-based")
    if args.mode == "targeted" and (not args.targets or not args.targets_build):
        raise ValueError("targeted requires --targets and --targets-build")
    if args.mode != "targeted" and (args.targets or args.targets_build):
        raise ValueError("Target inputs are only used in targeted mode")
    source_before = identity(args.source)
    target_before = identity(args.targets, full=True) if args.targets else None
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    result = preview(args, out) if args.mode == "preview" else analyze(args, out, build)
    if source_before != identity(args.source) or (args.targets and target_before != identity(args.targets, full=True)):
        raise ValueError("Input changed during analysis")
    manifest = dict(status="complete_descriptive_unverified", source=source_before, target=target_before,
                    build=build, build_verification="explicit declarations and any mapped row build; no liftover or reference FASTA verification",
                    settings=vars(args), analysis=result,
                    source_identity_policy="size and mtime only; the large source is not hashed or scanned a second time",
                    scan_policy="bounded prefix" if args.mode == "preview" else "one sequential full source pass; no idx/tabix random access",
                    inference="Descriptive source-summary evidence; populations remain separate; no gnomAD+AllofUs count pooling",
                    missingness="not_in_summary is unreported, never AC=0 or allele extinction; reporting/callability policy may be unknown",
                    ref_policy="REF count is used only if explicitly mapped; >=2 supported ALTs does not by itself demonstrate >=3 observed bases",
                    count_policy="Only reported, consistent AC/AN can establish AC thresholds; AF alone remains frequency evidence",
                    frequency_policy="Marginal AF values may have unverified callable denominators; sums above one are flagged, never renormalized or interpreted as a joint nucleotide distribution; gvs_max_af is a cross-population indicator, never a ninth population",
                    enrichment="Not estimated. Independent genome-wide enrichment requires an AllofUs eligible-site denominator, matching callability and covariates, and its own background model; gnomAD expected counts are not transferable",
                    independence="Cohort overlap and independent replication are not established")
    (out / "allofus_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    return manifest


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True, help="Local text/gzip summary; default mode reads only a bounded prefix")
    p.add_argument("--build", choices=tuple(BUILDS) + ("GRCh37", "GRCh38"), help="Required without a build-specific preset")
    p.add_argument("--preset", choices=("annovar-allofus-hg38",), help="User-provided 14-column tab-separated AF-only schema; GRCh38 and 1-based inclusive")
    p.add_argument("--out", required=True)
    p.add_argument("--mode", choices=("preview", "targeted", "normalize"), default="preview")
    p.add_argument("--preview-lines", type=int, default=12)
    p.add_argument("--preview-chars", type=int, default=65536)
    p.add_argument("--delimiter", choices=("tab", "whitespace"), default="tab")
    p.add_argument("--columns", help="Exact comma-separated column order for a headerless file; never guessed")
    p.add_argument("--chrom-column", default="chrom")
    p.add_argument("--pos-column", default="pos")
    p.add_argument("--end-column", help="Optional end; one-base SNV interval is checked")
    p.add_argument("--ref-column", default="ref")
    p.add_argument("--alt-column", default="alt")
    p.add_argument("--build-column", help="Optional per-row build; must agree with --build")
    p.add_argument("--coordinate-system", choices=("1-based", "0-based"), help="0-based positions are starts of half-open intervals")
    p.add_argument("--population", help="User-declared label for exactly the selected metric columns")
    p.add_argument("--af-column")
    p.add_argument("--ac-column")
    p.add_argument("--an-column")
    p.add_argument("--ref-ac-column", help="Optional actual REF allele count; never computed as AN minus displayed ALTs")
    p.add_argument("--population-map", help='JSON object, e.g. {"ALL":{"af":"AF"},"AFR":{"af":"AF_afr"}}; no automatic population aliases')
    p.add_argument("--targets", help="Existing TSV with chrom,pos,ref,alts and optional site_id/family/group/build")
    p.add_argument("--targets-build", choices=tuple(BUILDS) + ("GRCh37", "GRCh38"))
    p.add_argument("--min-ac", type=int, default=5)
    p.add_argument("--af-tolerance", type=float, default=1e-6, help="Absolute rounding tolerance for AF vs real AC/AN; never used to infer counts")
    p.add_argument("--progress-every", type=int, default=1000000)
    p.add_argument("--source-version", default="unverified", help="User-declared release; never inferred from filename or sample count")
    p.add_argument("--reporting-policy", default="unknown", help="Document source censoring/QC policy if known")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    run(args)
    print(f"[allofus] {args.mode} complete; output={Path(args.out).resolve()}", flush=True)


if __name__ == "__main__":
    main()
