#!/usr/bin/env python3
"""Consolidated burden workflow; internal stages share one explicit module."""
from __future__ import annotations


# ========================================================================
# assoc common
# ========================================================================
_description_assoc_common = 'Shared, side-effect-free contracts for the optional individual association module.\n\nNo phenotype pipeline is sourced and no scientific inputs are inferred from row\norder.  These helpers deliberately do not import any PCA or window-scan module.'
"""Shared, side-effect-free contracts for the optional individual association module.

No phenotype pipeline is sourced and no scientific inputs are inferred from row
order.  These helpers deliberately do not import any PCA or window-scan module.
"""

import copy
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).resolve().parents[2] / "0f"))
from ukb_followup import resolve_followup, validate_record

import numpy as np
import pandas as pd

VERSION = "allele3-assoc-1.0.0"
DEFAULT_Y = ["cvd_cad", "ra", "depression"]
MISSING = {"", ".", "NA", "NaN", "nan", "None", "null", "NULL", "<NA>"}
DEFAULT_CONFIG = {
    "phenotypes": {
        "id_col": "auto",
        "baseline_col": "date_attend", "birth_col": "birth_date",
        "lost_col": "date_lost", "death_col": "date_death",
        "diagnosis_template": "fod_icd10_{Y}",
        "evidence_template": "icd10Ct_{Y}", "case_template": None,
        "ancestry_column": "ethnic.c", "ancestry_values": ["White", "1"],
        "bbc_columns": [], "outcomes": {},
    },
    "models": {
        "endpoints": ["cox", "logistic"], "logistic_kind": "prevalent",
        "exposures": ["extra_copy_rate", "multi_het_rate"],
        "primary_exposure": "extra_copy_rate",
        "covariates": ["age", "sex", "tdi"] + [f"PC{k}" for k in range(1, 11)],
        "categorical": ["sex"],
        "enhanced_covariates": ["F_ROH", "heterozygosity", "other_rare_burden"],
        "adjustments": ["standard", "enhanced"],
        "technical_covariates": ["call_rate"], "quantiles": 4,
        "min_n": 100, "min_events": 20, "min_controls": 20,
        "min_events_per_parameter": 10, "min_group_n": 20,
        "min_group_events": 5, "min_call_rate": .95,
        "relatedness_policy": "require_unrelated",
        "min_locus_mac": 20, "max_loci": 100,
        "max_dosage_rows": 10000000, "condition_limit": 1e8,
        "max_iter": 100, "ties": "efron",
    },
}


def csv_words(value):
    if value is None:
        return []
    if isinstance(value, str):
        return list(dict.fromkeys(v.strip() for v in value.split(",") if v.strip()))
    return list(dict.fromkeys(str(v) for v in value))


def load_config(path=None):
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    fu = resolve_followup()
    if path:
        supplied = json.loads(Path(path).read_text())
        validate_record(supplied, "phenotypes", "end_date", fu)
        supplied.get("phenotypes", {}).pop("end_date", None)
        for section in ("phenotypes", "models"):
            if section in supplied:
                if not isinstance(supplied[section], dict):
                    raise ValueError(f"{section} config must be an object")
                unknown = set(supplied[section]) - set(cfg[section])
                if unknown:
                    raise ValueError(f"Unknown {section} config keys: {sorted(unknown)}")
                cfg[section].update(supplied[section])
    for section, keys in {
        "phenotypes": ["ancestry_values", "bbc_columns"],
        "models": ["endpoints", "exposures", "covariates", "categorical", "adjustments",
                   "enhanced_covariates", "technical_covariates"],
    }.items():
        for key in keys:
            cfg[section][key] = csv_words(cfg[section][key])
    m = cfg["models"]
    if not m["endpoints"] or set(m["endpoints"]) - {"cox", "logistic"}:
        raise ValueError("models.endpoints must contain cox and/or logistic")
    if m["logistic_kind"] not in {"prevalent", "ever", "explicit"}:
        raise ValueError("logistic_kind must be prevalent, ever or explicit")
    if m["relatedness_policy"] not in {"require_unrelated", "fixed_effect_exploratory"}:
        raise ValueError("Related samples require an unrelated keep list or explicit fixed_effect_exploratory")
    if set(m["adjustments"]) - {"standard", "enhanced"} or not m["adjustments"]:
        raise ValueError("adjustments must be standard and/or enhanced")
    if not m["exposures"] or m["primary_exposure"] not in m["exposures"]:
        raise ValueError("primary_exposure must occur in exposures")
    for key in ("min_n", "min_events", "min_controls", "min_group_n", "min_group_events",
                "min_locus_mac", "max_loci", "max_dosage_rows", "max_iter"):
        if int(m[key]) != m[key] or m[key] < 1:
            raise ValueError(f"{key} must be a positive integer")
    if m["quantiles"] not in (0, 2, 3, 4, 5):
        raise ValueError("quantiles must be 0 (disabled), 2, 3, 4 or 5")
    if not 0 <= m["min_call_rate"] <= 1 or m["min_events_per_parameter"] <= 0:
        raise ValueError("Invalid support thresholds")
    if m["condition_limit"] <= 1 or m["ties"] not in ("efron", "breslow"):
        raise ValueError("Invalid numerical model settings")
    cfg["phenotypes"]["end_date"] = fu["date"]
    cfg["followup_provenance"] = fu
    return cfg


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (pd.Timestamp,)):
        return value.isoformat()
    if value is pd.NA or value is pd.NaT:
        return None
    return value


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(json_safe(value), ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sample_hash(ids):
    return hashlib.sha256("\n".join(sorted(str(v) for v in ids)).encode()).hexdigest()


def file_record(path):
    p = Path(path).resolve()
    return {"path": str(p), "bytes": p.stat().st_size, "sha256": sha256_file(p)}


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header) != len(set(header)):
        raise ValueError(f"Duplicate columns in {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                       na_values=sorted(MISSING), low_memory=False)


def write_tsv(path, data, columns=None):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    frame = data if isinstance(data, pd.DataFrame) else pd.DataFrame(data, columns=columns)
    temporary = p.with_name(p.name + ".tmp")
    frame.to_csv(temporary, sep="\t", index=False, na_rep="NA",
                 compression="gzip" if str(p).endswith(".gz") else None)
    temporary.replace(p)


def validate_ids(values, context, unique=True):
    """Preserve text IDs; normalize genuine safe numeric R IDs, never row indices."""
    ids = []
    for value in values:
        if pd.isna(value):
            raise ValueError(f"Missing sample ID in {context}")
        if isinstance(value, (float, np.floating, int, np.integer)):
            if not np.isfinite(value) or value != int(value) or abs(value) >= 2**53:
                raise ValueError(f"Unsafe numeric sample ID in {context}; export IDs as text")
            value = str(int(value))
        else:
            value = str(value)
        if value in MISSING or value != value.strip() or any(x in value for x in "\t\r\n"):
            raise ValueError(f"Blank/ambiguous sample ID in {context}")
        ids.append(value)
    if unique and len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate sample IDs in {context}")
    return pd.Series(ids, index=getattr(values, "index", None), dtype="string")


def numeric(values, name, allow_missing=True):
    values = pd.Series(values)
    out = pd.to_numeric(values, errors="coerce")
    bad = values.notna() & (out.isna() | ~np.isfinite(out))
    if bad.any() or (not allow_missing and out.isna().any()):
        raise ValueError(f"Invalid numeric values in {name}")
    if ((out.abs() > 0) & (out.abs() < np.finfo(float).tiny)).any():
        raise ValueError(f"Subnormal numeric values in {name}; possible unsupported R integer64 encoding, use the R export adapter")
    return out.astype(float)


def adjust_all_tests(frame):
    """All predeclared attempted rows count; unestimable tests enter as p=1.

    Missing p values retain NA adjusted values in the public table.  Counting
    their hypotheses in the denominator avoids rewarding selective failures.
    """
    frame = frame.copy()
    if frame.empty:
        return frame, pd.DataFrame(columns=["family", "n_planned", "n_estimable", "method"])
    p = pd.to_numeric(frame["p"], errors="coerce")
    if ((p.dropna() < 0) | (p.dropna() > 1)).any():
        raise ValueError("P values outside [0,1]")
    summaries = []

    def correct(indices, name, suffix):
        indices = np.asarray(indices)
        val = p.loc[indices].fillna(1).to_numpy(float)
        order = np.argsort(val, kind="stable")
        n = len(order)
        ranked = val[order] * n / np.arange(1, n + 1)
        bh = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
        by = (bh * np.sum(1 / np.arange(1, n + 1))).clip(0, 1)
        valid = p.loc[indices].notna().to_numpy()[order]
        frame.loc[indices[order], "q_bh" + suffix] = np.where(valid, bh, np.nan)
        frame.loc[indices[order], "q_by" + suffix] = np.where(valid, by, np.nan)
        summaries.append({"family": name, "n_planned": n, "n_estimable": int(valid.sum()),
                          "n_unestimable": int(n - valid.sum()),
                          "method": "BH_and_BY; unestimable planned tests included as p=1"})

    correct(frame.index.to_numpy(), "ALL_ATTEMPTED", "_all")
    for family, group in frame.groupby("family", sort=True, dropna=False):
        correct(group.index.to_numpy(), str(family), "_family")
    return frame, pd.DataFrame(summaries)


# ========================================================================
# assoc catalog
# ========================================================================
_description_assoc_catalog = 'Freeze a phenotype-independent, external three/four-nucleotide catalog.\n\nInput is a common-callset nucleotide-count TSV, including the existing\ngenotype_validation.tsv format. REF is an assembly label, never the major\nallele. Rank2=rank3 ties are excluded. The extra alleles are the externally\ncount-supported rank3 AND rank4 bases; an unmeasured/unsupported fourth base\nis not silently added. These ranks are not a mutation chronology.\n\nBED4 sets use 0-based half-open intervals. ALL means all supplied catalog loci,\nnot every callable genomic base. No disease or case/control label is read.'
#!/usr/bin/env python3
"""Freeze a phenotype-independent, external three/four-nucleotide catalog.

Input is a common-callset nucleotide-count TSV, including the existing
genotype_validation.tsv format. REF is an assembly label, never the major
allele. Rank2=rank3 ties are excluded. The extra alleles are the externally
count-supported rank3 AND rank4 bases; an unmeasured/unsupported fourth base
is not silently added. These ranks are not a mutation chronology.

BED4 sets use 0-based half-open intervals. ALL means all supplied catalog loci,
not every callable genomic base. No disease or case/control label is read.
"""

import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
from contextlib import contextmanager
from decimal import Decimal, InvalidOperation
import sys

from analysis_integrity import InputPath, identity, transactional

CATALOG_VERSION = "assoc-catalog-1.1"
CATALOG_BASES = "ACGT"
CATALOG_MISSING = {"", ".", "NA", "N/A", "nan", "NaN", "None"}
CATALOG_AUTO_SETS = ("ALL", "TRIALLELIC", "QUADALLELIC")
CATALOG_CATALOG_FIELDS = [
    "site_id", "catalog_id", "build", "chrom", "pos", "ref", "source_cohort", "source_scope",
    "population", "AN", "call_rate", "A_AC", "C_AC", "G_AC", "T_AC",
    "n_measured_bases", "n_observed_bases", "n_supported_bases", "supported_bases",
    "top2_bases", "extra_bases", "unmeasured_bases", "unsupported_observed_bases",
    "baseline_base", "baseline_tied", "rank1_base", "rank1_ac", "rank2_base", "rank2_ac",
    "rank3_base", "rank3_ac", "rank4_base", "rank4_ac", "min_ac", "extra_policy",
]
CATALOG_AUDIT_FIELDS = ["site_id", "chrom", "pos", "population", "status", "reason", "n_supported_bases", "AN"]
CATALOG_SET_FIELDS = ["site_id", "catalog_id", "set_id"]


@contextmanager
def catalog_text_open(path, mode="rt"):
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, mode, newline="") as handle:
        yield handle


def catalog_table_rows(path, required=()):
    """Strict TSV reader; the first column may carry a VCF-style # prefix."""
    with catalog_text_open(path) as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("##"))
        reader = csv.DictReader(lines, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        aliases = {"CHROM": "chrom", "POS": "pos", "REF": "ref"}
        fields = [aliases.get(x.lstrip("#"), x.lstrip("#")) for x in reader.fieldnames]
        if len(fields) != len(set(fields)) or not set(required) <= set(fields):
            raise ValueError(f"{path}: need unique columns including {','.join(required)}")
        reader.fieldnames = fields
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Malformed TSV row in {path}")
            yield row


def catalog_build(value):
    aliases = {"37": "GRCh37", "38": "GRCh38", "hg19": "GRCh37", "hg38": "GRCh38"}
    value = aliases.get(str(value), str(value))
    if value not in {"GRCh37", "GRCh38"}:
        raise ValueError("An explicit GRCh37 or GRCh38 assembly is required")
    return value


def catalog_chromosome(value):
    value = str(value).removeprefix("chr")
    if value not in {str(i) for i in range(1, 23)}:
        raise ValueError("Only autosomes 1-22 are supported; X/Y require a separate ploidy contract")
    return value


def catalog_integer(value, name, minimum=0):
    try:
        number = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if not number.is_finite() or number != number.to_integral_value() or number < minimum:
        raise ValueError(f"{name} must be a finite integer >= {minimum}")
    return int(number)


def catalog_token(value, name):
    value = str(value)
    if not value or value in CATALOG_MISSING or any(c in value for c in "\t\r\n"):
        raise ValueError(f"Invalid {name}")
    return value


def catalog_clean(row, fields):
    return {k: "NA" if row.get(k) is None else row.get(k, "NA") for k in fields}


def catalog_reject_preview(row):
    """A preview anywhere in the supplied table invalidates a formal catalog."""
    if "is_preview" in row:
        preview = catalog_integer(row["is_preview"], "is_preview")
        if preview not in (0, 1):
            raise ValueError("is_preview must be 0 or 1 when supplied")
        if preview:
            raise ValueError("Preview count rows (is_preview=1) cannot define a formal association catalog")


@contextmanager
def catalog_writer(path, fields):
    with catalog_text_open(path, "wt") as handle:
        output = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        output.writeheader()
        yield lambda row: output.writerow(catalog_clean(row, fields))


class catalog_BedSets:
    def __init__(self, path=None):
        grouped = defaultdict(list)
        self.names = set()
        if path:
            with catalog_text_open(path) as handle:
                for line_no, line in enumerate(handle, 1):
                    if not line.strip() or line.startswith(("#", "track ", "browser ")):
                        continue
                    fields = line.split()
                    if len(fields) < 4:
                        raise ValueError(f"BED4 required at {path}:{line_no}")
                    c = catalog_chromosome(fields[0]); lo = catalog_integer(fields[1], "BED start"); hi = catalog_integer(fields[2], "BED end", 1)
                    if hi <= lo:
                        raise ValueError("BED interval must be nonempty and half-open")
                    name = catalog_token(fields[3], "BED set_id")
                    if name in CATALOG_AUTO_SETS:
                        raise ValueError(f"BED set name is reserved: {name}")
                    self.names.add(name); grouped[c].append((lo, hi, name))
        self.index = {}
        for c, values in grouped.items():
            values.sort(); maximum = []; current = -1
            for _, hi, _ in values:
                current = max(current, hi); maximum.append(current)
            self.index[c] = ([r[0] for r in values], values, maximum)

    def labels(self, chrom, pos):
        if chrom not in self.index:
            return set()
        starts, intervals, maximum = self.index[chrom]
        point = pos - 1; i = bisect_right(starts, point) - 1; labels = set()
        while i >= 0 and maximum[i] > point:
            lo, hi, name = intervals[i]
            if lo <= point < hi:
                labels.add(name)
            i -= 1
        return labels


def catalog_catalog_row(row, args):
    """Return a frozen catalog row or a scientific exclusion reason."""
    catalog_reject_preview(row)
    for field in ("build", "genome_build", "assembly"):
        if row.get(field) not in CATALOG_MISSING | {None} and catalog_build(row[field]) != catalog_build(args.build):
            raise ValueError("Count table assembly conflicts with --build; no liftover is performed")
    c = catalog_chromosome(row["chrom"]); pos = catalog_integer(row["pos"], "pos", 1); ref = row["ref"].upper()
    if ref not in CATALOG_BASES or len(ref) != 1:
        raise ValueError("REF must be a single A/C/G/T nucleotide")
    sid = catalog_token(row.get("site_id") or f"{c}:{pos}:{ref}", "site_id")
    base = dict(site_id=sid, chrom="chr"+c, pos=pos, population=row["population"], AN=None, n_supported_bases=None)
    for field, expected in (("evaluable", "1"), ("represented_status", "all_expected_alts_encoded"), ("record_status", "ok")):
        if field in row and row[field] != expected:
            return None, dict(base, status="excluded", reason="source_"+field+"_not_eligible")
    an = catalog_integer(row["AN"], "AN")
    base["AN"] = an
    if an == 0:
        return None, dict(base, status="excluded", reason="source_AN_zero")
    if "callable_all_four" in row and row["callable_all_four"] != "1":
        return None, dict(base, status="excluded", reason="source_callable_all_four_not_eligible")
    counts = {b: catalog_integer(row[b+"_AC"], b+"_AC") for b in CATALOG_BASES if row[b+"_AC"] not in CATALOG_MISSING}
    if sum(counts.values()) != an:
        raise ValueError(f"{sid}: measured A/C/G/T counts do not sum to the common AN")
    if "encoded_alts" in row:
        encoded = set(row["encoded_alts"].split(",")) - {"", ".", "NA"}
        if not encoded <= set(CATALOG_BASES)-{ref} or set(counts) != encoded | {ref}:
            raise ValueError(f"{sid}: genotype-validation counts and encoded nucleotide labels disagree")
    if row.get("n_called") not in CATALOG_MISSING | {None} and "encoded_alts" in row:
        if 2 * catalog_integer(row["n_called"], "n_called") != an:
            raise ValueError(f"{sid}: AN differs from the common diploid n_called")
    call_rate = None
    if row.get("call_rate") not in CATALOG_MISSING | {None}:
        call_rate = float(row["call_rate"])
        if not math.isfinite(call_rate) or not 0 <= call_rate <= 1:
            raise ValueError("call_rate must be finite in [0,1]")
        if call_rate < args.min_call_rate:
            return None, dict(base, AN=an, status="excluded", reason="source_call_rate_below_threshold")
    ranked = sorted((b for b in counts if counts[b] >= args.min_ac), key=lambda b: (-counts[b], b))
    base.update(AN=an, n_supported_bases=len(ranked))
    if len(ranked) < 3:
        return None, dict(base, status="excluded", reason="fewer_than_three_count_supported_bases")
    if counts[ranked[1]] == counts[ranked[2]]:
        return None, dict(base, status="excluded", reason="top2_boundary_tie")
    result = dict(base, catalog_id=args.catalog_id, build=catalog_build(args.build), ref=ref,
                  source_cohort=args.source_cohort, source_scope=args.source_scope,
                  call_rate=call_rate, n_measured_bases=len(counts), n_observed_bases=sum(v > 0 for v in counts.values()),
                  supported_bases=",".join(ranked), top2_bases=",".join(ranked[:2]), extra_bases=",".join(ranked[2:]),
                  unmeasured_bases=",".join(b for b in CATALOG_BASES if b not in counts) or ".",
                  unsupported_observed_bases=",".join(b for b in CATALOG_BASES if 0 < counts.get(b, 0) < args.min_ac) or ".",
                  baseline_base=ranked[0], baseline_tied=int(counts[ranked[0]] == counts[ranked[1]]),
                  min_ac=args.min_ac, extra_policy="externally_count_supported_rank3_and_rank4")
    result.update({b+"_AC": counts.get(b) for b in CATALOG_BASES})
    for index in range(4):
        b = ranked[index] if index < len(ranked) else None
        result[f"rank{index+1}_base"] = b
        result[f"rank{index+1}_ac"] = counts[b] if b else None
    return result, dict(base, status="included", reason="three_or_four_count_supported_bases")


def catalog_input_paths(raw):
    return [InputPath(raw[key]) for key in ("counts", "sets_bed") if raw.get(key)]


@transactional("manifest.json", input_prefixes=(), resolve_inputs=catalog_input_paths)
def catalog_run(args):
    if args.min_ac < 1 or not 0 <= args.min_call_rate <= 1:
        raise ValueError("min-ac must be positive; min-call-rate must be in [0,1]")
    for name in ("source_cohort", "population", "catalog_id", "source_scope"):
        catalog_token(getattr(args, name), name)
    declared_build = catalog_build(args.build)
    required = ["chrom", "pos", "ref", "population", "AN", *(b+"_AC" for b in CATALOG_BASES)]
    sets = catalog_BedSets(args.sets_bed); out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    audit = Counter(); chosen = []; seen_sites = set(); seen_ids = set()
    with catalog_writer(out/"catalog_audit.tsv.gz", CATALOG_AUDIT_FIELDS) as emit:
        for row in catalog_table_rows(args.counts, required):
            audit["input_rows"] += 1
            # Check before selecting a population: a preview file must not become
            # eligible simply by selecting another stratum from that file.
            catalog_reject_preview(row)
            if row["population"] != args.population:
                audit["other_population_rows"] += 1; continue
            audit["selected_population_rows"] += 1
            c = catalog_chromosome(row["chrom"]); pos = catalog_integer(row["pos"], "pos", 1)
            if (c, pos) in seen_sites:
                raise ValueError(f"Duplicate selected-population position: chr{c}:{pos}")
            seen_sites.add((c, pos))
            item, status = catalog_catalog_row(row, args); emit(status); audit[status["reason"]] += 1
            if item is not None:
                if item["site_id"] in seen_ids:
                    raise ValueError("Duplicate site_id assigned to different positions")
                seen_ids.add(item["site_id"]); chosen.append(item)
    if not audit["selected_population_rows"]:
        raise ValueError(f"No rows found for explicit population {args.population!r}")
    chosen.sort(key=lambda r: (int(catalog_chromosome(r["chrom"])), r["pos"]))
    set_counts = Counter({s: 0 for s in (*CATALOG_AUTO_SETS, *sorted(sets.names))})
    with catalog_writer(out/"catalog.tsv.gz", CATALOG_CATALOG_FIELDS) as emit, catalog_writer(out/"site_sets.tsv.gz", CATALOG_SET_FIELDS) as membership:
        for row in chosen:
            emit(row)
            names = {"ALL", "TRIALLELIC" if row["n_supported_bases"] == 3 else "QUADALLELIC"}
            names |= sets.labels(catalog_chromosome(row["chrom"]), row["pos"])
            for name in sorted(names):
                membership(dict(site_id=row["site_id"], catalog_id=args.catalog_id, set_id=name)); set_counts[name] += 1
    with catalog_writer(out/"sets.tsv", ["set_id", "n_targets", "scope"]) as emit:
        for name, number in sorted(set_counts.items()):
            emit(dict(set_id=name, n_targets=number, scope="all supplied catalog loci" if name == "ALL" else "subset of supplied catalog loci"))
    return dict(status="complete" if chosen else "no_usable_contrast", schema=CATALOG_VERSION,
                build=declared_build, catalog_id=args.catalog_id, source_cohort=args.source_cohort,
                source_scope=args.source_scope, population=args.population, n_sites=len(chosen),
                chromosomes=sorted({catalog_chromosome(r["chrom"]) for r in chosen}, key=int),
                set_counts=dict(set_counts), audit=dict(audit), min_ac=args.min_ac,
                counts_input=identity(args.counts, full=Path(args.counts).stat().st_size <= 64*1024*1024),
                catalog_sha256=identity(out/"catalog.tsv.gz", full=True)["sha256"],
                site_sets_sha256=identity(out/"site_sets.tsv.gz", full=True)["sha256"],
                ranking="fixed external population count order; rank2=rank3 excluded; top1 ties use alphabetical baseline",
                extra_policy="externally count-supported rank3 AND rank4; unmeasured or unsupported bases remain unlisted",
                count_contract="measured nucleotide counts sum to one common AN; common-callset provenance remains an input-provider assertion",
                optional_input_flags={
                    "is_preview": "when present, must be 0 or 1; any 1 aborts the entire formal catalog",
                    "callable_all_four": "when present, only 1 is eligible; absent fields retain the legacy common-count contract"},
                phenotype_use="none; source case/control columns are not used for ranking, inclusion or BED membership",
                limitations=["ALL is all supplied catalog loci, not genome-wide callable bases.",
                             "Counts do not establish mutation age, selection, or independent source/target cohorts.",
                             "No source sample roster is inferred from site summaries; overlap must be audited separately.",
                             "Unencoded fourth bases are not assigned measured zero; AN coherence does not prove four-base callability.",
                             "A missing source call-rate field is recorded as unknown, not assigned 1.",
                             "AN=0 and explicitly ineligible four-base calls are audited exclusions, never association controls.",
                             "Absent optional eligibility flags retain the legacy input contract and do not prove four-base callability."])


def catalog_parser():
    p = argparse.ArgumentParser(description=_description_assoc_catalog)
    p.add_argument("--counts", required=True)
    p.add_argument("--population", required=True)
    p.add_argument("--build", required=True, choices=("GRCh37", "GRCh38", "37", "38"))
    p.add_argument("--source-cohort", required=True)
    p.add_argument("--source-scope", default="prespecified_catalog")
    p.add_argument("--catalog-id", default="main", help="Independent named catalogs retain distinct identities")
    p.add_argument("--min-ac", type=int, default=2)
    p.add_argument("--min-call-rate", type=float, default=.9)
    p.add_argument("--sets-bed", help="Optional BED4 with nonreserved set_id in column 4")
    p.add_argument("--out", required=True)
    return p


def catalog_main(argv=None):
    args = catalog_parser().parse_args(argv)
    try:
        result = catalog_run(args)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"assoc-catalog: {exc}") from exc
    print(json.dumps({k: result[k] for k in ("status", "schema", "n_sites", "catalog_id")}, ensure_ascii=False))


def _entry_assoc_catalog():
    catalog_main()


# ========================================================================
# assoc exposure
# ========================================================================
_description_assoc_exposure = 'Stream individual exposures to a fixed external three/four-base catalog.\n\nNative multiallelic PGEN or VCF GTs supply actual nucleotide identities. Split,\nduplicate and mixed same-position records are conservatively unavailable.\nNeither a missing record nor an unencoded nucleotide is converted to zero.\nOnly externally count-supported rank3/rank4 bases in extra_bases contribute to\nextra_copy_rate. A genotype containing an unlisted base is jointly missing.\n\nOne genotype locus is held at a time; accumulation is samples x named sets,\nnot samples x loci. Chromosome shards can be resumed separately and merged\nonly with identical catalogs, sample rosters and QC and complete disjoint scope.\nNo phenotype is used in either counting or the saved callability risk set.'
#!/usr/bin/env python3
"""Stream individual exposures to a fixed external three/four-base catalog.

Native multiallelic PGEN or VCF GTs supply actual nucleotide identities. Split,
duplicate and mixed same-position records are conservatively unavailable.
Neither a missing record nor an unencoded nucleotide is converted to zero.
Only externally count-supported rank3/rank4 bases in extra_bases contribute to
extra_copy_rate. A genotype containing an unlisted base is jointly missing.

One genotype locus is held at a time; accumulation is samples x named sets,
not samples x loci. Chromosome shards can be resumed separately and merged
only with identical catalogs, sample rosters and QC and complete disjoint scope.
No phenotype is used in either counting or the saved callability risk set.
"""
from contextlib import ExitStack
import re

from analysis_integrity import validate_manifest
from allele_measurement import open_pgen_with_pvar, panel_pvar

EXPOSURE_VERSION = "assoc-exposure-1.0"
EXPOSURE_EXPOSURE_FIELDS = ["sample_id", "set_id", "catalog_id", "cohort_id", "build", "source_scope",
                   "n_callable", "n_targets", "n_missing", "extra_copies", "extra_carrier_sites",
                   "heterozygous_sites", "extra_copy_rate", "extra_carrier_rate", "multi_het_rate",
                   "call_rate", "missing_fraction", "eligible", "eligibility_reason"]
EXPOSURE_RISK_FIELDS = ["sample_id", "set_id", "n_callable", "n_targets", "n_missing", "call_rate",
               "missing_fraction", "eligible", "eligibility_reason"]
EXPOSURE_SAMPLE_FIELDS = ["sample_id", "cohort_id", "population", "n_callable", "n_targets", "call_rate", "eligible"]
EXPOSURE_DOSAGE_FIELDS = ["sample_id", "site_id", "baseline_base", "allele", "allele_rank", "allele_role",
                 "dosage", "evaluable", "measurement_status", "catalog_id", "cohort_id", "build"]
EXPOSURE_QC_FIELDS = ["site_id", "chrom", "pos", "measurement_status", "representation", "n_records",
             "n_samples", "n_callable", "n_missing_gt", "n_quality_failed", "n_unlisted_base",
             "n_non_diploid", "dp_missing_field", "gq_missing_field", "ad_missing_field", "call_rate"]
EXPOSURE_RAW_FIELDS = ("n_callable", "extra_copies", "extra_carrier_sites", "heterozygous_sites")
EXPOSURE_QC_NAMES = ("min_dp", "min_gq", "min_allele_balance", "min_sample_call_rate", "max_missing_fraction")


def exposure_parse_chromosomes(value, default=()):
    if value is None:
        return sorted(set(default), key=int)
    result = []
    for item in str(value).split(","):
        item = item.strip().removeprefix("chr")
        if re.fullmatch(r"\d+-\d+", item):
            lo, hi = map(int, item.split("-"))
            if lo > hi:
                raise ValueError("Chromosome range must increase")
            result.extend(catalog_chromosome(str(i)) for i in range(lo, hi+1))
        else:
            result.append(catalog_chromosome(item))
    if not result or len(result) != len(set(result)):
        raise ValueError("Specify unique autosomal chromosomes")
    return sorted(result, key=int)


def exposure_load_catalog(path, declared_build):
    rows = []; seen = set(); ids = set(); metadata = set()
    for raw in catalog_table_rows(path, CATALOG_CATALOG_FIELDS):
        if catalog_build(raw["build"]) != catalog_build(declared_build):
            raise ValueError("Catalog and target assembly differ; no implicit liftover")
        c = catalog_chromosome(raw["chrom"]); pos = catalog_integer(raw["pos"], "pos", 1)
        sid = catalog_token(raw["site_id"], "site_id"); ref = raw["ref"]
        if ref not in CATALOG_BASES or len(ref) != 1 or sid in ids or (c, pos) in seen:
            raise ValueError("Catalog needs unique site IDs/positions and canonical REF")
        seen.add((c, pos)); ids.add(sid)
        minimum = catalog_integer(raw["min_ac"], "min_ac", 1)
        counts = {b: catalog_integer(raw[b+"_AC"], b+"_AC") for b in CATALOG_BASES if raw[b+"_AC"] not in CATALOG_MISSING}
        if sum(counts.values()) != catalog_integer(raw["AN"], "AN", 1):
            raise ValueError(f"{sid}: catalog counts do not sum to common AN")
        supported = sorted((b for b in counts if counts[b] >= minimum), key=lambda b: (-counts[b], b))
        if len(supported) < 3 or counts[supported[1]] == counts[supported[2]]:
            raise ValueError(f"{sid}: invalid supported bases or ambiguous top2 boundary")
        for key, values in (("supported_bases", supported), ("top2_bases", supported[:2]), ("extra_bases", supported[2:])):
            if raw[key].split(",") != values:
                raise ValueError(f"{sid}: {key} differs from frozen external count ranking")
        if raw["baseline_base"] != supported[0] or catalog_integer(raw["n_supported_bases"], "n_supported_bases") != len(supported):
            raise ValueError(f"{sid}: baseline/count contract differs")
        for index, b in enumerate(supported, 1):
            if raw[f"rank{index}_base"] != b or catalog_integer(raw[f"rank{index}_ac"], "rank AC") != counts[b]:
                raise ValueError(f"{sid}: rank labels/counts differ")
        metadata.add((raw["catalog_id"], raw["source_cohort"], raw["source_scope"], raw["population"]))
        rows.append(dict(raw, chrom=c, pos=pos, supported=supported, extras=supported[2:]))
    if not rows or len(metadata) != 1:
        raise ValueError("Use one nonempty fixed catalog with one source cohort/population/scope")
    rows.sort(key=lambda r: (int(r["chrom"]), r["pos"]))
    return rows


def exposure_sets_path(raw):
    if raw.get("site_sets"):
        return str(raw["site_sets"])
    path = Path(raw["catalog"]).parent/"site_sets.tsv.gz"
    return str(path) if path.is_file() else None


def exposure_memberships(catalog, path):
    by_id = {r["site_id"]: r for r in catalog}
    result = {r["site_id"]: {"ALL", "TRIALLELIC" if len(r["supported"]) == 3 else "QUADALLELIC"} for r in catalog}
    names = set(CATALOG_AUTO_SETS); seen = set()
    if path:
        for row in catalog_table_rows(path, ("site_id", "set_id")):
            sid = row["site_id"]; name = catalog_token(row["set_id"], "set_id")
            if sid not in by_id or (sid, name) in seen:
                raise ValueError("Unknown/duplicate site-set membership")
            if row.get("catalog_id", by_id[sid]["catalog_id"]) != by_id[sid]["catalog_id"]:
                raise ValueError("Site-set catalog_id differs")
            if name in {"TRIALLELIC", "QUADALLELIC"} and name not in result[sid]:
                raise ValueError("Automatic tri/quad set conflicts with supported nucleotide count")
            seen.add((sid, name)); result[sid].add(name); names.add(name)
    return result, sorted(names)


def exposure_template_path(template, chrom):
    if template.count("{chrom}") != 1 or "{" in template.replace("{chrom}", "") or "}" in template.replace("{chrom}", ""):
        raise ValueError("Source template requires exactly one literal {chrom}")
    return str(Path(template.format(chrom=chrom)).expanduser().resolve())


def exposure_source_paths(raw):
    """Metadata-only path expansion for workflow preflight; raw is a dict.

No genotype file is opened. If the catalog is not yet built, every declared
chromosome is expanded. Existing catalogs limit data reads to targeted scopes;
an empty shard without --samples uses one catalog source as a roster anchor.
"""
    actual = None
    if Path(raw["catalog"]).is_file():
        actual = {catalog_chromosome(r["chrom"]) for r in catalog_table_rows(raw["catalog"], ("chrom",))}
    selected = exposure_parse_chromosomes(raw.get("chromosomes"), actual or [str(i) for i in range(1, 23)])
    needed = sorted(set(selected) & actual, key=int) if actual is not None else selected
    if not needed and not raw.get("samples") and actual:
        needed = [min(actual, key=int)]
    result = []
    if raw.get("pgen_template"):
        for c in needed:
            prefix = exposure_template_path(raw["pgen_template"], c)
            present = [p for p in (prefix+".pvar", prefix+".pvar.zst") if Path(p).is_file()]
            if len(present) > 1:
                raise ValueError("Ambiguous PVAR and PVAR.ZST companions")
            result.append(dict(kind="pgen", chrom=c, pgen=prefix+".pgen", pvar=present[0] if present else prefix+".pvar", psam=prefix+".psam"))
    elif raw.get("vcf_template"):
        for c in needed:
            path = exposure_template_path(raw["vcf_template"], c)
            result.append(dict(kind="vcf", chrom=c, vcf=path,
                               index=next((path+s for s in (".tbi", ".csi") if Path(path+s).is_file()), None)))
    elif raw.get("vcf"):
        path = str(Path(raw["vcf"]).expanduser().resolve())
        if needed or not raw.get("samples"):
            result.append(dict(kind="vcf", chrom="all", vcf=path,
                               index=next((path+s for s in (".tbi", ".csi") if Path(path+s).is_file()), None)))
    else:
        raise ValueError("Supply one PGEN template, VCF template, or VCF")
    return result


def exposure_run_inputs(raw):
    files = [raw["catalog"], exposure_sets_path(raw)]
    files += [raw.get(x) for x in ("samples", "dosage_sites", "upstream_qc")]
    for source in exposure_source_paths(raw):
        files += [source.get(k) for k in ("pgen", "pvar", "psam", "vcf", "index")]
    return [InputPath(x) for x in files if x]


def exposure_sample_selection(path):
    rows = list(catalog_table_rows(path, ("sample_id",))) if path else []
    samples = [catalog_token(r["sample_id"], "sample_id") for r in rows]
    if len(samples) != len(set(samples)) or (path and not samples):
        raise ValueError("Sample roster must be nonempty and have unique exact-string sample_id")
    return samples, {r["sample_id"]: r.get("population", "unknown") for r in rows}


def exposure_read_psam(path):
    ids = []; names = None
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            parts = line.split()
            if names is None:
                names = [x.lstrip("#") for x in parts]
                if len(names) != len(set(names)) or "IID" not in names:
                    raise ValueError("PSAM requires a unique IID column")
                column = names.index("IID"); continue
            if len(parts) != len(names):
                raise ValueError("Malformed PSAM row")
            ids.append(catalog_token(parts[column], "IID"))
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("PSAM IIDs must be nonempty and globally unique; ambiguous FID/IID is not guessed")
    return ids


def exposure_sample_hash(samples):
    return hashlib.sha256(("\n".join(sorted(samples))+"\n").encode()).hexdigest()


def exposure_verify_header_build(line, expected):
    found = {catalog_build(x) for x in re.findall(r"(?:GRCh(?:37|38)|hg(?:19|38))", str(line))}
    if found and found != {catalog_build(expected)}:
        raise ValueError("Source header assembly differs from the declared build")
    return found


def exposure_pvar_targets(path, targets, expected_build):
    """One PVAR pass; retain ALL same-position companions before QC."""
    found = defaultdict(list); names = None; index = 0; headers = set()
    if str(path).endswith(".zst"):
        import io, zstandard
        stack = ExitStack(); binary = stack.enter_context(open(path, "rb"))
        stream = stack.enter_context(zstandard.ZstdDecompressor().stream_reader(binary))
        handle = stack.enter_context(io.TextIOWrapper(stream))
    else:
        stack = ExitStack(); handle = stack.enter_context(open(path))
    with stack:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("##"):
                headers |= exposure_verify_header_build(line, expected_build); continue
            values = line.split()
            if names is None:
                names = [x.lstrip("#") for x in values]
                if len(names) != len(set(names)) or not {"CHROM", "POS", "REF", "ALT"} <= set(names):
                    raise ValueError("PVAR needs unique CHROM/POS/REF/ALT columns")
                continue
            if len(values) != len(names):
                raise ValueError("Malformed PVAR row")
            row = dict(zip(names, values)); c = row["CHROM"].removeprefix("chr"); pos = catalog_integer(row["POS"], "PVAR POS", 1)
            if (c, pos) in targets:
                found[c, pos].append(dict(index=index, ref=row["REF"], alleles=[row["REF"], *row["ALT"].split(",")],
                                          passed=row.get("FILTER", "PASS") in {"PASS", "."}))
            index += 1
    return found, index, "recognized_header" if headers else "declared_only_no_recognized_header"


def exposure_locus_status(records, site):
    if not records:
        return "record_not_found"
    if len(records) != 1:
        return "split_or_duplicate_records_rejected"
    record = records[0]; alleles = record["alleles"]
    if len(set(alleles)) != len(alleles) or not set(alleles) <= set(CATALOG_BASES) or not 2 <= len(alleles) <= 4:
        return "non_snv_or_invalid_record"
    if record["ref"] != site["ref"]:
        return "reference_label_mismatch"
    if len(alleles) < 3:
        return "not_native_multiallelic"
    if not record["passed"]:
        return "record_filter_failed"
    if not set(site["supported"]) <= set(alleles):
        return "catalog_supported_base_not_encoded"
    return "native_multiallelic"


def exposure_missing_call(n, reason, records=0):
    return np.full((n, 2), -9, dtype=np.int8), dict(measurement_status=reason, representation="unavailable",
        n_records=records, n_missing_gt=n, n_quality_failed=0, n_non_diploid=0,
        dp_missing_field=0, gq_missing_field=0, ad_missing_field=0)


def exposure_nucleotide_calls(indices, alleles, site):
    if indices.ndim != 2 or indices.shape[1] != 2:
        raise ValueError("Expected diploid genotype indices")
    if np.any((indices < 0) & (indices != -9)) or np.any(indices >= len(alleles)):
        raise ValueError("Invalid genotype indices; PGEN/PVAR or VCF allele contract differs")
    valid = np.all(indices >= 0, axis=1)
    calls = np.full(indices.shape, -9, dtype=np.int8)
    codes = np.asarray([CATALOG_BASES.index(b) for b in alleles], dtype=np.int8)
    calls[valid] = codes[indices[valid]]
    allowed = np.asarray([CATALOG_BASES.index(b) for b in site["supported"]], dtype=np.int8)
    unlisted = valid & ~np.isin(calls, allowed).all(axis=1)
    calls[unlisted] = -9
    return calls, int(unlisted.sum())


def exposure_vcf_gt(record, site, samples, args):
    """Read one native record; all allele dosages share one callability mask."""
    n = len(samples); pairs = np.full((n, 2), -9, dtype=np.int16); audit = Counter()
    format_fields = set(record.format.keys())
    checks = (("DP", args.min_dp > 0), ("GQ", args.min_gq > 0), ("AD", args.min_allele_balance > 0))
    for field, enabled in checks:
        absent = enabled and field not in format_fields
        audit[field.lower()+"_missing_field"] = int(absent)
        if absent and not args.upstream_qc:
            raise ValueError(f"VCF lacks {field} QC at {site['site_id']}; provide explicit upstream-QC evidence or explicitly disable that criterion")
    for i, sample in enumerate(samples):
        call = record.samples[sample]; gt = call.get("GT")
        if gt is None or any(a is None for a in gt):
            audit["n_missing_gt"] += 1; continue
        if len(gt) != 2:
            audit["n_non_diploid"] += 1; continue
        if any(not isinstance(a, int) or a < 0 or a >= len(record.alleles) for a in gt):
            raise ValueError("Invalid VCF GT allele index")
        failed = False
        # Explicit per-genotype FILTER failures must not be rescued by good
        # DP/GQ. A dot is the VCF missing/unapplied filter convention.
        if "FT" in format_fields:
            ft = call.get("FT")
            values = ft if isinstance(ft, (tuple, list)) else (ft,)
            failed |= any(v not in {None, "", ".", "PASS"} for v in values)
        for field, threshold in (("DP", args.min_dp), ("GQ", args.min_gq)):
            if threshold > 0 and field in format_fields:
                value = call.get(field)
                failed |= value is None or not isinstance(value, (int, float)) or not math.isfinite(value) or value < threshold
        if args.min_allele_balance > 0 and gt[0] != gt[1] and "AD" in format_fields:
            ad = call.get("AD")
            if ad is None or len(ad) != len(record.alleles) or any(v is None or v < 0 for v in ad) or sum(ad) <= 0:
                failed = True
            else:
                fractions = [ad[a]/sum(ad) for a in gt]
                failed |= any(x < args.min_allele_balance or x > 1-args.min_allele_balance for x in fractions)
        if failed:
            audit["n_quality_failed"] += 1; continue
        pairs[i] = gt
    calls, unlisted = exposure_nucleotide_calls(pairs, record.alleles, site)
    audit.update(n_unlisted_base=unlisted)
    return calls, dict(audit, measurement_status="native_multiallelic", representation="native_multiallelic", n_records=1)


def exposure_vcf_records(handle, source, sites):
    """Indexed exact queries, or one bounded-memory pass over a standard VCF."""
    if source.get("index"):
        contigs = list(handle.header.contigs)
        for site in sites:
            choices = [c for c in (site["chrom"], "chr"+site["chrom"]) if c in contigs]
            if len(choices) > 1:
                raise ValueError("VCF contains ambiguous chr/plain chromosome aliases")
            records = [r for r in handle.fetch(choices[0], site["pos"]-1, site["pos"]) if r.pos == site["pos"]] if choices else []
            yield site, records
        return
    wanted = {(s["chrom"], s["pos"]): s for s in sites}; seen = set(); current = None; records = []; previous = None
    def finish(key, group):
        if key in wanted:
            seen.add(key)
            return wanted[key], group
        return None
    for record in handle:
        c = record.contig.removeprefix("chr")
        if c not in {str(i) for i in range(1, 23)}:
            continue
        key = (c, record.pos); order = (int(c), record.pos)
        if previous is not None and order < previous:
            raise ValueError("Unindexed VCF must be autosomally coordinate sorted")
        previous = order
        if current is not None and key != current:
            result = finish(current, records)
            if result: yield result
            records = []
        current = key
        if key in wanted:
            records.append(record)
    if current is not None:
        result = finish(current, records)
        if result: yield result
    for site in sites:
        if (site["chrom"], site["pos"]) not in seen:
            yield site, []


def exposure_allocate(n, names, maximum_mb):
    required = n * len(names) * len(EXPOSURE_RAW_FIELDS) * np.dtype(np.int64).itemsize
    if maximum_mb <= 0 or required > maximum_mb*1024*1024:
        raise ValueError(f"samples x sets aggregation needs {required} bytes; reduce named sets or raise --max-aggregation-mb explicitly")
    return np.zeros((len(names), len(EXPOSURE_RAW_FIELDS), n), dtype=np.int64)


def exposure_qc_parameters(args):
    values = {k: getattr(args, k) for k in EXPOSURE_QC_NAMES}
    if args.min_dp < 0 or args.min_gq < 0 or not 0 <= args.min_allele_balance <= .5:
        raise ValueError("DP/GQ thresholds must be nonnegative and allele balance in [0,.5]")
    if not 0 <= args.min_sample_call_rate <= 1 or not 0 <= args.max_missing_fraction <= 1:
        raise ValueError("Callability thresholds must be in [0,1]")
    for value in values.values():
        if not math.isfinite(value): raise ValueError("QC thresholds must be finite")
    if args.pgen_template and not args.upstream_qc:
        raise ValueError("PGEN has no DP/GQ/read balance: --upstream-qc evidence file is required")
    if args.upstream_qc:
        if not Path(args.upstream_qc).is_file() or not Path(args.upstream_qc).stat().st_size:
            raise ValueError("upstream-qc must be a nonempty QC evidence/declaration file")
        values["upstream_qc_sha256"] = identity(args.upstream_qc, full=True)["sha256"]
    else:
        values["upstream_qc_sha256"] = None
    values["split_policy"] = "reject_all_split_duplicate_or_mixed_records"
    values["unlisted_base_policy"] = "joint_sample_locus_missing"
    return values


def exposure_target_counts(sites, member, names):
    index = {name: i for i, name in enumerate(names)}; result = np.zeros(len(names), dtype=np.int64)
    for site in sites:
        for name in member[site["site_id"]]: result[index[name]] += 1
    return result


def exposure_publish_summaries(out, samples, populations, names, counts, targets, metadata, qc):
    all_index = names.index("ALL"); all_total = int(targets[all_index])
    all_called = counts[all_index, 0]
    global_ok = (all_called >= all_total*qc["min_sample_call_rate"]-1e-12) & ((all_total-all_called) <= all_total*qc["max_missing_fraction"]+1e-12) if all_total else np.zeros(len(samples), dtype=bool)
    with catalog_writer(out/"exposures.tsv.gz", EXPOSURE_EXPOSURE_FIELDS) as emit, catalog_writer(out/"riskset.tsv.gz", EXPOSURE_RISK_FIELDS) as risk, catalog_writer(out/"samples.tsv", EXPOSURE_SAMPLE_FIELDS) as roster:
        for k, sample in enumerate(samples):
            roster(dict(sample_id=sample, cohort_id=metadata["cohort_id"], population=populations.get(sample, "unknown"),
                        n_callable=int(all_called[k]), n_targets=all_total, call_rate=float(all_called[k])/all_total if all_total else None,
                        eligible=int(global_ok[k])))
            for j, name in enumerate(names):
                n, extra, carrier, hetero = map(int, counts[j, :, k]); total = int(targets[j]); missing = total-n
                if not 0 <= n <= total or not 0 <= extra <= 2*n or not 0 <= carrier <= n or not 0 <= hetero <= n or not carrier <= extra <= 2*carrier:
                    raise ValueError("Invalid aggregate numerator/denominator")
                rate = n/total if total else None; fraction = missing/total if total else None
                reasons = []
                if not total: reasons.append("empty_set")
                elif not n: reasons.append("zero_callable_loci")
                if total and (rate+1e-12 < qc["min_sample_call_rate"] or fraction-1e-12 > qc["max_missing_fraction"]):
                    reasons.append("set_callability_below_threshold")
                if not global_ok[k]: reasons.append("ALL_callability_below_threshold")
                row = dict(metadata, sample_id=sample, set_id=name, n_callable=n, n_targets=total, n_missing=missing,
                           extra_copies=extra, extra_carrier_sites=carrier, heterozygous_sites=hetero,
                           extra_copy_rate=extra/n if n else None, extra_carrier_rate=carrier/n if n else None,
                           multi_het_rate=hetero/n if n else None, call_rate=rate, missing_fraction=fraction,
                           eligible=int(not reasons), eligibility_reason=";".join(reasons) or "eligible")
                emit(row); risk(row)


def exposure_shared_metadata(catalog, args):
    first = catalog[0]
    return dict(catalog_id=first["catalog_id"], cohort_id=catalog_token(args.cohort_id, "cohort_id"), build=catalog_build(args.build), source_scope=first["source_scope"])


def exposure_dosage_selection(path, catalog):
    ids = [] if not path else [r["site_id"] for r in catalog_table_rows(path, ("site_id",))]
    available = {r["site_id"] for r in catalog}
    if len(ids) != len(set(ids)) or not set(ids) <= available:
        raise ValueError("Dosage selection has duplicate or unknown site IDs")
    return set(ids)


@transactional("manifest.json", input_prefixes=(), resolve_inputs=exposure_run_inputs)
def exposure_run(args):
    qc = exposure_qc_parameters(args); catalog = exposure_load_catalog(args.catalog, args.build); source = exposure_source_paths(vars(args))
    requested = exposure_parse_chromosomes(args.chromosomes, {s["chrom"] for s in catalog})
    sites = [s for s in catalog if s["chrom"] in requested]
    set_file = exposure_sets_path(vars(args)); member, names = exposure_memberships(catalog, set_file)
    samples, populations = exposure_sample_selection(args.samples); explicit_samples = bool(args.samples)
    metadata = exposure_shared_metadata(catalog, args); out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    dosages = exposure_dosage_selection(args.dosage_sites, catalog)
    handles = {}; source_rosters = {}; source_audits = []; input_evidence = []
    with ExitStack() as stack:
        for i, item in enumerate(source):
            for field in ("pgen", "pvar", "psam", "vcf", "index"):
                if item.get(field):
                    if not Path(item[field]).is_file(): raise ValueError(f"Missing source {field}: {item[field]}")
                    input_evidence.append(identity(item[field]))
            if item["kind"] == "pgen":
                ids = exposure_read_psam(item["psam"])
            else:
                import pysam
                handle = pysam.VariantFile(item["vcf"]); stack.callback(handle.close); handles[i] = handle
                ids = list(handle.header.samples)
                if not ids or len(ids) != len(set(ids)):
                    raise ValueError("VCF needs unique individual sample IDs and GT; a sites-only VCF is insufficient")
                for header in handle.header.records:
                    exposure_verify_header_build(str(header), args.build)
            if not samples: samples = ids.copy()
            if explicit_samples:
                if not set(samples) <= set(ids): raise ValueError("Selected sample IDs are missing from a chromosome source")
            elif set(samples) != set(ids):
                raise ValueError("Chromosome sample rosters differ; supply an explicit common --samples roster")
            source_rosters[i] = ids
        if not samples:
            raise ValueError("No sample roster; an empty chromosome shard needs --samples or a readable catalog source")
        totals = exposure_target_counts(sites, member, names); counts = exposure_allocate(len(samples), names, args.max_aggregation_mb)
        whole_panel_rows = len(samples)*sum(len(s["supported"])-1 for s in catalog if s["site_id"] in dosages)
        if args.max_dosage_rows < 1 or whole_panel_rows > args.max_dosage_rows:
            raise ValueError(f"Requested single-site dosage panel would write {whole_panel_rows} rows, above --max-dosage-rows")
        emit_qc = stack.enter_context(catalog_writer(out/"site_qc.tsv.gz", EXPOSURE_QC_FIELDS))
        emit_dosage = stack.enter_context(catalog_writer(out/"locus_dosages.tsv.gz", EXPOSURE_DOSAGE_FIELDS)) if args.dosage_sites else None
        set_index = {s: i for i, s in enumerate(names)}; audit = Counter(); visited = set()

        def consume(site, calls, detail):
            sid = site["site_id"]
            if sid in visited: raise ValueError("A catalog locus was counted twice")
            visited.add(sid)
            called = np.all(calls >= 0, axis=1)
            extra_codes = [CATALOG_BASES.index(b) for b in site["extras"]]
            copies = np.isin(calls, extra_codes).sum(axis=1); copies[~called] = 0
            carrier = (copies > 0) & called
            hetero = (calls[:, 0] != calls[:, 1]) & called
            for name in member[sid]:
                j = set_index[name]; counts[j, 0] += called; counts[j, 1] += copies
                counts[j, 2] += carrier; counts[j, 3] += hetero
            fields = dict.fromkeys(EXPOSURE_QC_FIELDS, 0)
            fields.update(detail, site_id=sid, chrom="chr"+site["chrom"], pos=site["pos"],
                          n_samples=len(samples), n_callable=int(called.sum()), call_rate=float(called.mean()))
            emit_qc(fields); audit[fields["measurement_status"]] += 1
            if sid in dosages:
                for k, sample in enumerate(samples):
                    for rank, b in enumerate(site["supported"][1:], 2):
                        emit_dosage(dict(metadata, sample_id=sample, site_id=sid, baseline_base=site["baseline_base"], allele=b,
                                         allele_rank=rank, allele_role="extra" if rank >= 3 else "rank2",
                                         dosage=int((calls[k] == CATALOG_BASES.index(b)).sum()) if called[k] else None,
                                         evaluable=int(called[k]), measurement_status="called_native_joint_GT" if called[k] else "joint_sample_locus_missing"))

        for i, item in enumerate(source):
            selected_sites = [s for s in sites if item["chrom"] == "all" or s["chrom"] == item["chrom"]]
            if not selected_sites: continue
            ids = source_rosters[i]; id_index = {sid: k for k, sid in enumerate(ids)}
            indices = np.asarray([id_index[sid] for sid in samples], dtype=np.int64)
            if item["kind"] == "pgen":
                found, n_variants, header_status = exposure_pvar_targets(item["pvar"], {(s["chrom"], s["pos"]) for s in selected_sites}, args.build)
                with open_pgen_with_pvar(item["pgen"], item["pvar"], len(ids)) as reader:
                    if reader.get_raw_sample_ct() != len(ids) or reader.get_variant_ct() != n_variants:
                        raise ValueError("PGEN dimensions differ from actual PSAM/PVAR")
                    buffer = np.empty(2*len(ids), dtype=np.int32)
                    for site in selected_sites:
                        records = found.get((site["chrom"], site["pos"]), []); status = exposure_locus_status(records, site)
                        if status != "native_multiallelic":
                            calls, detail = exposure_missing_call(len(samples), status, len(records))
                        else:
                            record = records[0]; reader.read_alleles(record["index"], buffer)
                            raw = buffer.reshape((-1, 2))[indices]
                            calls, unlisted = exposure_nucleotide_calls(raw, record["alleles"], site)
                            detail = dict(measurement_status=status, representation="native_multiallelic", n_records=1,
                                          n_missing_gt=int((raw < 0).any(axis=1).sum()), n_unlisted_base=unlisted)
                        consume(site, calls, detail)
                source_audits.append(dict(chrom=item["chrom"], kind="pgen", n_source_samples=len(ids), n_selected_samples=len(samples),
                                          pvar_rows=n_variants, build_verification=header_status, genotype_qc="upstream_QC_declaration_only_no_read_fields"))
            else:
                handle = handles[i]; handle.subset_samples(samples)
                for site, records in exposure_vcf_records(handle, item, selected_sites):
                    summaries = [dict(ref=r.ref, alleles=list(r.alleles or ()), passed=not(set(r.filter.keys())-{"PASS", "."})) for r in records]
                    status = exposure_locus_status(summaries, site)
                    if status != "native_multiallelic":
                        calls, detail = exposure_missing_call(len(samples), status, len(records))
                    else:
                        calls, detail = exposure_vcf_gt(records[0], site, samples, args)
                    consume(site, calls, detail)
                source_audits.append(dict(chrom=item["chrom"], kind="vcf", n_source_samples=len(ids), n_selected_samples=len(samples),
                                          mode="indexed_exact_targets" if item.get("index") else "single_coordinate_sorted_stream"))
        if visited != {s["site_id"] for s in sites}:
            raise ValueError("Not every requested catalog locus was explicitly assessed")
    exposure_publish_summaries(out, samples, populations, names, counts, totals, metadata, qc)
    result = dict(status="complete", schema=EXPOSURE_VERSION, mode="chromosome_counts", **metadata,
                  chromosomes=requested, catalog_chromosomes=sorted({s["chrom"] for s in catalog}, key=int),
                  full_catalog_coverage={s["chrom"] for s in catalog} <= set(requested),
                  catalog_sha256=identity(args.catalog, full=True)["sha256"],
                  site_sets_sha256=identity(set_file, full=True)["sha256"] if set_file else None,
                  sample_ids_sha256=exposure_sample_hash(samples), n_samples=len(samples), n_sites=len(sites), sets=names,
                  set_target_counts={name: int(totals[j]) for j, name in enumerate(names)}, qc_parameters=qc,
                  evidence_kind=args.evidence_kind, source_cohort=catalog[0]["source_cohort"], source_population=catalog[0]["population"],
                  source_inputs=input_evidence, source_audits=source_audits, site_qc=dict(audit),
                  dosage_sites_sha256=identity(args.dosage_sites, full=True)["sha256"] if args.dosage_sites else None,
                  planned_full_catalog_dosage_rows=whole_panel_rows, max_dosage_rows=args.max_dosage_rows,
                  genotype_memory="one locus; accumulators samples x sets; no sample x locus matrix",
                  sample_overlap="not_assessed; distinct cohort labels do not establish independence",
                  riskset="saved phenotype-independent per-set AND ALL callability gates; recomputed after shard merge",
                  definitions={"extra_copy_rate": "externally supported rank3+rank4 copies / callable catalog loci; range0..2",
                               "extra_carrier_rate": "loci carrying either supported extra base / callable catalog loci",
                               "multi_het_rate": "heterozygous supported genotypes at catalog loci / callable catalog loci",
                               "ALL": "all supplied catalog loci in this explicit chromosome scope; not genome-wide callable bases"},
                  limitations=["A diploid individual carries at most two copies at a locus; tri/quad describes the population catalog.",
                               "REF is not major/ancestral and rank3/rank4 are current external count ranks, not temporal order.",
                               "Unknown or externally unsupported bases are not silently included; carriers of an unlisted base are missing jointly at that locus.",
                               "Split/duplicate records are rejected rather than reconstructing REF from a single ALT.",
                               "PGEN upstream-QC evidence is supplied provenance, not independently revalidated reads; imputed_hardcall is not direct sequencing evidence.",
                               "Raw counts remain even when callability risk-set gates fail; zero-call rates are NA."])
    return result


def exposure_merge_inputs(raw):
    result = [InputPath(raw["catalog"])]
    if exposure_sets_path(raw): result.append(InputPath(exposure_sets_path(raw)))
    for path in raw["shard_dir"]:
        result.append(InputPath(path, "directory"))
        # Track actual checkpoint files too, so a shard refreshed while the
        # reducer runs cannot produce a falsely complete mixed-version result.
        for name in ("manifest.json", "exposures.tsv.gz", "riskset.tsv.gz", "site_qc.tsv.gz", "samples.tsv", "locus_dosages.tsv.gz"):
            candidate = Path(path)/name
            if candidate.is_file(): result.append(InputPath(candidate))
    return result


@transactional("manifest.json", input_prefixes=(), resolve_inputs=exposure_merge_inputs)
def exposure_merge(args):
    paths = [Path(p).expanduser().resolve() for p in args.shard_dir]
    if not paths or len(paths) != len(set(paths)):
        raise ValueError("Merge requires distinct completed shard directories")
    manifests = [validate_manifest(p, require_evidence=True) for p in paths]
    first = manifests[0]
    if any(m.get("schema") != EXPOSURE_VERSION or m.get("mode") != "chromosome_counts" for m in manifests):
        raise ValueError("Only completed assoc-exposure chromosome count shards can be merged")
    catalog = exposure_load_catalog(args.catalog, first["build"])
    set_file = exposure_sets_path(vars(args)); member, names = exposure_memberships(catalog, set_file)
    expected = exposure_parse_chromosomes(args.chromosomes, {s["chrom"] for s in catalog})
    keys = ("catalog_sha256", "site_sets_sha256", "sample_ids_sha256", "qc_parameters", "catalog_id", "cohort_id", "build",
            "evidence_kind", "source_scope", "sets", "n_samples", "dosage_sites_sha256", "max_dosage_rows")
    if any(any(m.get(k) != first.get(k) for k in keys) for m in manifests[1:]):
        raise ValueError("Shard catalog, sets, sample roster, cohort, dosage panel or QC parameters differ")
    if first["catalog_sha256"] != identity(args.catalog, full=True)["sha256"] or first["site_sets_sha256"] != (identity(set_file, full=True)["sha256"] if set_file else None):
        raise ValueError("Merge catalog/site sets differ from the measured shards")
    if first["sets"] != names: raise ValueError("Shard set names differ from catalog membership")
    present = []
    for manifest in manifests:
        current = exposure_parse_chromosomes(",".join(manifest["chromosomes"]))
        if set(present) & set(current): raise ValueError("Overlapping chromosome shards cannot be summed")
        present += current
    if set(present) != set(expected):
        raise ValueError("Incomplete or unexpected chromosome coverage; a missing shard is never filled with zero")
    roster = list(catalog_table_rows(paths[0]/"samples.tsv", ("sample_id", "population")))
    samples = [r["sample_id"] for r in roster]; populations = {r["sample_id"]: r["population"] for r in roster}
    if len(samples) != len(set(samples)) or len(samples) != first["n_samples"] or exposure_sample_hash(samples) != first["sample_ids_sha256"]:
        raise ValueError("Shard roster identity is invalid")
    counts = exposure_allocate(len(samples), names, args.max_aggregation_mb)
    index = {sample: k for k, sample in enumerate(samples)}; set_index = {name: k for k, name in enumerate(names)}
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    sites = [s for s in catalog if s["chrom"] in expected]; totals = exposure_target_counts(sites, member, names)
    metadata = {k: first[k] for k in ("catalog_id", "cohort_id", "build", "source_scope")}
    with ExitStack() as stack:
        emit_qc = stack.enter_context(catalog_writer(out/"site_qc.tsv.gz", EXPOSURE_QC_FIELDS))
        emit_dose = stack.enter_context(catalog_writer(out/"locus_dosages.tsv.gz", EXPOSURE_DOSAGE_FIELDS)) if first["dosage_sites_sha256"] else None
        for path, manifest in zip(paths, manifests):
            shard_roster = {r["sample_id"]: r["population"] for r in catalog_table_rows(path/"samples.tsv", ("sample_id", "population"))}
            if shard_roster != populations: raise ValueError("Shard participant identities/populations differ")
            selected = [s for s in catalog if s["chrom"] in manifest["chromosomes"]]
            local_targets = exposure_target_counts(selected, member, names)
            seen = np.zeros((len(names), len(samples)), dtype=bool)
            for row in catalog_table_rows(path/"exposures.tsv.gz", EXPOSURE_EXPOSURE_FIELDS):
                if row["sample_id"] not in index or row["set_id"] not in set_index:
                    raise ValueError("Unknown participant/set in shard exposures")
                k = index[row["sample_id"]]; j = set_index[row["set_id"]]
                if seen[j, k]: raise ValueError("Duplicate sample x set in shard exposures")
                seen[j, k] = True
                values = [catalog_integer(row[name], name) for name in EXPOSURE_RAW_FIELDS]
                total = catalog_integer(row["n_targets"], "n_targets"); missing = catalog_integer(row["n_missing"], "n_missing")
                n, extra, carrier, het = values
                if total != local_targets[j] or missing != total-n or not 0 <= n <= total or not 0 <= extra <= 2*n or not 0 <= carrier <= n or not 0 <= het <= n or not carrier <= extra <= 2*carrier:
                    raise ValueError("Invalid shard raw counts or target denominator")
                counts[j, :, k] += values
            if not seen.all(): raise ValueError("Incomplete sample x set shard rows; missing is not zero")
            qc_ids = set()
            for row in catalog_table_rows(path/"site_qc.tsv.gz", EXPOSURE_QC_FIELDS):
                if row["site_id"] in qc_ids: raise ValueError("Duplicate shard site-QC row")
                qc_ids.add(row["site_id"]); emit_qc(row)
            if qc_ids != {s["site_id"] for s in selected}: raise ValueError("Incomplete shard locus audit")
            if emit_dose:
                for row in catalog_table_rows(path/"locus_dosages.tsv.gz", EXPOSURE_DOSAGE_FIELDS): emit_dose(row)
    exposure_publish_summaries(out, samples, populations, names, counts, totals, metadata, first["qc_parameters"])
    return dict(status="complete", schema=EXPOSURE_VERSION, mode="merged_counts", **metadata, chromosomes=expected,
                catalog_chromosomes=sorted({s["chrom"] for s in catalog}, key=int), full_catalog_coverage={s["chrom"] for s in catalog} <= set(expected),
                catalog_sha256=first["catalog_sha256"], site_sets_sha256=first["site_sets_sha256"],
                sample_ids_sha256=first["sample_ids_sha256"], n_samples=len(samples), n_sites=len(sites), sets=names,
                qc_parameters=first["qc_parameters"], evidence_kind=first["evidence_kind"],
                source_cohort=first["source_cohort"], source_population=first["source_population"],
                sample_overlap=first["sample_overlap"],
                set_target_counts={name: int(totals[j]) for j, name in enumerate(names)},
                dosage_sites_sha256=first["dosage_sites_sha256"], definitions=first["definitions"], limitations=first["limitations"],
                source_scope_description="ALL means all supplied catalog loci in the declared complete merge scope",
                merge_rule="sum raw numerators and denominators across disjoint complete chromosome shards; never average rates",
                riskset="ALL and per-set callability gates recomputed after summation; shard eligibility is not treated as final",
                shard_manifests=[identity(p/"manifest.json", full=True) for p in paths])


def exposure_parser():
    p = argparse.ArgumentParser(description=_description_assoc_exposure)
    p.add_argument("--catalog", required=True)
    p.add_argument("--site-sets")
    p.add_argument("--build", required=True, choices=("GRCh37", "GRCh38", "37", "38"))
    p.add_argument("--cohort-id", required=True)
    p.add_argument("--evidence-kind", required=True, choices=("direct_gt", "imputed_hardcall"))
    sources = p.add_mutually_exclusive_group(required=True)
    sources.add_argument("--pgen-template")
    sources.add_argument("--vcf-template")
    sources.add_argument("--vcf")
    p.add_argument("--chromosomes", help="Explicit chromosome shard/scope, e.g. 1 or 1-22")
    p.add_argument("--samples", help="Unique exact-string sample_id TSV; explicit common roster")
    p.add_argument("--upstream-qc", help="Nonempty QC provenance/declaration file; required for PGEN or unavailable VCF read fields")
    p.add_argument("--min-dp", type=float, default=10)
    p.add_argument("--min-gq", type=float, default=20)
    p.add_argument("--min-allele-balance", type=float, default=.2)
    p.add_argument("--min-sample-call-rate", type=float, default=.9)
    p.add_argument("--max-missing-fraction", type=float, default=.1)
    p.add_argument("--dosage-sites", help="Bounded, predefined site_id TSV; all nonbaseline supported bases are emitted jointly")
    p.add_argument("--max-dosage-rows", type=int, default=10000000)
    p.add_argument("--max-aggregation-mb", type=int, default=1024)
    p.add_argument("--out", required=True)
    return p


def exposure_merge_parser():
    p = argparse.ArgumentParser(description="Merge complete disjoint assoc-exposure chromosome checkpoints")
    p.add_argument("--catalog", required=True)
    p.add_argument("--site-sets")
    p.add_argument("--shard-dir", required=True, action="append")
    p.add_argument("--chromosomes", help="Required scope, including deliberately empty chromosome shards if requested")
    p.add_argument("--max-aggregation-mb", type=int, default=1024)
    p.add_argument("--out", required=True)
    return p


def exposure_main(argv=None):
    values = list(sys.argv[1:] if argv is None else argv)
    try:
        if values and values[0] == "merge":
            result = exposure_merge(exposure_merge_parser().parse_args(values[1:]))
        else:
            result = exposure_run(exposure_parser().parse_args(values))
    except (OSError, ValueError, RuntimeError) as exc:
        raise SystemExit(f"assoc-exposure: {exc}") from exc
    print(json.dumps({k: result[k] for k in ("status", "schema", "n_samples", "n_sites", "chromosomes")}, ensure_ascii=False))


def _entry_assoc_exposure():
    exposure_main()


# ========================================================================
# assoc phenotypes
# ========================================================================
_description_assoc_phenotypes = 'Read-only LE8/UKB phenotype bridge and audited incident-outcome construction.\n\nReads DIR/Rdata/all.rds (and explicitly requested BBC columns), or an explicit\nportable TSV.  It does not run LE8, rebuild disease definitions, impute missing\ncovariates, or turn a missing diagnosis column into an all-control outcome.'
#!/usr/bin/env python3
"""Read-only LE8/UKB phenotype bridge and audited incident-outcome construction.

Reads DIR/Rdata/all.rds (and explicitly requested BBC columns), or an explicit
portable TSV.  It does not run LE8, rebuild disease definitions, impute missing
covariates, or turn a missing diagnosis column into an all-control outcome.
"""
import datetime as dt


def phenotype_parse_dates(values, name):
    values = pd.Series(values)
    # Numeric days in an R Date are already converted by pyreadr.  Unlabelled
    # numeric dates could be Unix seconds, days or Excel serials: never guess.
    if pd.api.types.is_numeric_dtype(values) and values.notna().any():
        raise ValueError(f"{name}: numeric dates are ambiguous; export ISO dates using burden_export.R")
    present = values.notna()
    valid_objects = values.map(lambda x: isinstance(x, (dt.date, dt.datetime, pd.Timestamp, np.datetime64)))
    text = values.astype("string")
    iso = text.str.fullmatch(r"\d{4}-\d{2}-\d{2}", na=False)
    if (present & ~valid_objects & ~iso).any():
        raise ValueError(f"{name}: use valid YYYY-MM-DD dates, not partial/locale dates")
    normalized = values.map(lambda x: x.strftime("%Y-%m-%d")
                            if isinstance(x, (dt.date, dt.datetime, pd.Timestamp)) and not pd.isna(x)
                            else x)
    result = pd.to_datetime(normalized, format="%Y-%m-%d", errors="coerce")
    if (present & result.isna()).any():
        raise ValueError(f"Invalid calendar dates in {name}")
    return result


def phenotype_read_source(path):
    path = Path(path)
    if path.suffix.lower() == ".rds":
        try:
            import pyreadr
        except ImportError as exc:
            raise RuntimeError("Reading UKB RDS requires pyreadr; alternatively export a TSV with burden_export.R") from exc
        try:
            objects = pyreadr.read_r(str(path))
        except Exception as exc:
            raise ValueError("RDS cannot be represented safely by pyreadr. Use the read-only burden_export.R exporter; "
                             f"no fallback values were created. Reader error: {exc}") from exc
        if len(objects) != 1 or not isinstance(next(iter(objects.values())), pd.DataFrame):
            raise ValueError("Expected one R data.frame in all.rds")
        data = next(iter(objects.values())).reset_index(drop=True)
        # Normalize R missing numeric payloads while retaining factor labels.
        for col in data:
            if pd.api.types.is_numeric_dtype(data[col]):
                data[col] = data[col].mask(data[col].isna(), np.nan)
        return data
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return read_tsv(path)


def phenotype_id_column(data, setting="auto"):
    if setting != "auto":
        if setting not in data:
            raise ValueError(f"Missing ID column: {setting}")
        return setting
    candidates = [x for x in ("sample_id", "eid") if x in data]
    if not candidates:
        raise ValueError("Need sample_id or eid; participant IDs cannot be inferred from row numbers")
    if len(candidates) > 1:
        a = validate_ids(data[candidates[0]], candidates[0])
        b = validate_ids(data[candidates[1]], candidates[1])
        if not a.equals(b):
            raise ValueError("sample_id and eid disagree")
    return candidates[0]


def phenotype_source_fields(config, outcome):
    p = config["phenotypes"]
    override = p.get("outcomes", {}).get(outcome, {})
    if not isinstance(override, dict):
        raise ValueError(f"outcome mapping must be an object: {outcome}")
    allowed = {"diagnosis_col", "evidence_col", "case_col"}
    if set(override) - allowed:
        raise ValueError(f"Unknown outcome mapping keys for {outcome}")
    fields = {
        "diagnosis_col": p["diagnosis_template"].format(Y=outcome),
        "evidence_col": p["evidence_template"].format(Y=outcome) if p["evidence_template"] else None,
        "case_col": p["case_template"].format(Y=outcome) if p["case_template"] else None,
        "end_date": p["end_date"],
    }
    fields.update(override)
    return fields


def phenotype_normalize_cohort(data, outcomes, config, bbc=None):
    """Select and name existing columns; preserve NA and mark absent outcomes."""
    data = data.copy()
    if data.columns.duplicated().any():
        raise ValueError("Duplicate phenotype columns")
    p, m = config["phenotypes"], config["models"]
    source_id = phenotype_id_column(data, p["id_col"])
    ids = validate_ids(data[source_id], "phenotype source")
    data["sample_id"] = ids
    input_n = len(data)
    if p["ancestry_values"]:
        col = p["ancestry_column"]
        if not col or col not in data:
            raise ValueError("Configured ancestry filter column is absent; set ancestry_values=[] explicitly for another design")
        data = data.loc[data[col].astype("string").isin(p["ancestry_values"])].copy()
    if not len(data):
        raise ValueError("No participants remain after the declared phenotype cohort filter")
    if bbc is not None:
        bid = phenotype_id_column(bbc, p["id_col"])
        bbc = bbc.copy()
        bbc["sample_id"] = validate_ids(bbc[bid], "BBC source")
        requested = [x for x in p["bbc_columns"] if x not in data]
        absent = set(requested) - set(bbc)
        if absent:
            raise ValueError(f"Requested BBC columns absent: {sorted(absent)}")
        data = data.merge(bbc[["sample_id"] + requested], on="sample_id", how="left", validate="one_to_one")
    requested_bbc_missing = set(p["bbc_columns"]) - set(data)
    if requested_bbc_missing:
        raise ValueError(f"Missing requested BBC columns: {sorted(requested_bbc_missing)}")

    mappings = {"date_attend": p["baseline_col"], "birth_date": p["birth_col"],
                "date_lost": p["lost_col"], "date_death": p["death_col"]}
    selected = list(dict.fromkeys(["sample_id", p["ancestry_column"]] + m["covariates"]
                   + m["technical_covariates"] + m["enhanced_covariates"] + p["bbc_columns"]))
    out = data[[v for v in selected if v and v in data]].copy()
    audit = [{"outcome": "ALL", "metric": "input_participants", "value": input_n},
             {"outcome": "ALL", "metric": "after_declared_ancestry_filter", "value": len(data)}]
    for canonical, original in mappings.items():
        if original and original in data:
            out[canonical] = phenotype_parse_dates(data[original], original).dt.strftime("%Y-%m-%d")
        else:
            audit.append({"outcome": "ALL", "metric": "absent_column", "value": canonical})
    for outcome in outcomes:
        f = phenotype_source_fields(config, outcome)
        for key, canonical in (("diagnosis_col", "fod_icd10_" + outcome),
                               ("evidence_col", "icd10Ct_" + outcome),
                               ("case_col", "case_status_" + outcome)):
            source = f[key]
            if source and source in data:
                if key == "diagnosis_col":
                    out[canonical] = phenotype_parse_dates(data[source], source).dt.strftime("%Y-%m-%d")
                else:
                    out[canonical] = numeric(data[source], source)
                    if key == "case_col" and not out[canonical].dropna().isin([0, 1]).all():
                        raise ValueError(f"Explicit case column must be 0/1/NA: {source}")
                audit.append({"outcome": outcome, "metric": key, "value": source})
            else:
                audit.append({"outcome": outcome, "metric": "unavailable_" + key,
                              "value": source or "not_requested"})
        derived, counts = phenotype_derive_outcome(out, outcome, config)
        audit.extend({"outcome": outcome, "metric": k, "value": v} for k, v in counts.items())
    return out.reset_index(drop=True), audit


def phenotype_derive_outcome(data, outcome, config):
    """LE8 first-recorded onset semantics, with conservative undated-case handling.

    Diagnosis on the baseline date is prevalent. An event on the censoring date
    is incident. Death/loss/admin dates censor at their minimum. A diagnosis
    after censoring is not an incident event. Positive undated disease evidence
    is unresolved, not a control. Zero follow-up never enters a Cox risk set.
    """
    p = config["phenotypes"]
    f = phenotype_source_fields(config, outcome)
    end = phenotype_parse_dates(pd.Series([f["end_date"]]), "administrative end").iloc[0]
    if pd.isna(end):
        raise ValueError("Administrative end date must be explicit and valid")
    index = data.index
    out = pd.DataFrame(index=index)
    audit = {"source_n": len(data), "end_date": str(end.date()),
             "followup_provenance": config["followup_provenance"]}
    raw_case = "case_status_" + outcome
    if raw_case in data:
        out["explicit_case"] = numeric(data[raw_case], raw_case)
        if not out["explicit_case"].dropna().isin([0, 1]).all():
            raise ValueError(f"Explicit case status is not binary: {outcome}")
    else:
        out["explicit_case"] = np.nan
    required = ["date_attend", "fod_icd10_" + outcome]
    absent = [c for c in required if c not in data]
    if absent:
        for col in ("event", "time", "prevalent_case", "ever_case"):
            out[col] = np.nan
        out["eligible_cox"] = False
        out["eligible_prevalent"] = False
        out["eligible_ever"] = False
        out["outcome_reason"] = "missing_outcome_columns:" + ",".join(absent)
        audit.update(status="unavailable", reason=out["outcome_reason"].iloc[0] if len(out) else "missing_columns")
        return out, audit
    baseline = phenotype_parse_dates(data["date_attend"], "date_attend")
    onset = phenotype_parse_dates(data["fod_icd10_" + outcome], "fod_icd10_" + outcome)
    censor_inputs = [pd.Series(end, index=index)]
    for col in ("date_lost", "date_death"):
        if col in data:
            censor_inputs.append(phenotype_parse_dates(data[col], col))
        else:
            audit["missing_" + col] = "no_date_record_available; admin_censor_applies"
    censor = pd.concat(censor_inputs, axis=1).min(axis=1)
    before_birth = pd.Series(False, index=index)
    if "birth_date" in data:
        birth = phenotype_parse_dates(data["birth_date"], "birth_date")
        before_birth = (baseline.notna() & birth.notna() & baseline.lt(birth)) | (
            onset.notna() & birth.notna() & onset.lt(birth))
    prevalent = onset.notna() & baseline.notna() & onset.le(baseline)
    evidence_name = "icd10Ct_" + outcome
    positive_evidence = out["explicit_case"].eq(1)
    if evidence_name in data:
        evidence = numeric(data[evidence_name], evidence_name)
        if (evidence.dropna() < 0).any():
            raise ValueError(f"Negative disease evidence count: {evidence_name}")
        positive_evidence = positive_evidence | evidence.gt(0)
    undated = onset.isna() & positive_evidence
    audit["explicit_positive_undated"] = int((onset.isna() & out["explicit_case"].eq(1)).sum())
    audit["explicit_negative_dated"] = int((onset.notna() & out["explicit_case"].eq(0)).sum())
    valid = baseline.notna() & censor.ge(baseline) & ~before_birth
    incident = valid & ~prevalent & onset.notna() & onset.gt(baseline) & onset.le(censor)
    stop = censor.where(~incident, onset)
    duration = (stop - baseline).dt.days / 365.25
    out["time"] = duration
    out["event"] = incident.astype(int)
    out["eligible_cox"] = valid & ~prevalent & ~undated & duration.gt(0)
    out["prevalent_case"] = prevalent.astype(int)
    out["ever_case"] = (prevalent | incident).astype(int)
    out["eligible_prevalent"] = valid & ~undated
    out["eligible_ever"] = valid & ~undated
    out["outcome_reason"] = np.select(
        [baseline.isna(), before_birth, censor.lt(baseline), undated, prevalent, duration.le(0)],
        ["missing_baseline", "date_before_birth", "censor_before_baseline",
         "undated_disease_evidence", "baseline_prevalent_including_same_day", "zero_followup"],
        default="eligible_or_later_censored")
    out["censor_date"] = censor
    audit.update(status="derived", prevalent=int(prevalent.sum()),
                 same_day=int((onset.notna() & onset.eq(baseline)).sum()),
                 incident=int(incident.sum()), after_censor=int((valid & onset.gt(censor)).sum()),
                 unresolved_undated=int(undated.sum()), invalid_followup=int((~valid).sum()),
                 zero_followup=int((valid & duration.eq(0)).sum()),
                 eligible_cox=int(out["eligible_cox"].sum()),
                 eligible_prevalent=int(out["eligible_prevalent"].sum()))
    return out, audit


def phenotype_main(argv=None):
    parser = argparse.ArgumentParser(description=_description_assoc_phenotypes, allow_abbrev=False)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--cohort", help="Portable tab-delimited cohort (or an explicit RDS)")
    sources.add_argument("--dir-ukb-phe", "--ukb-phe", dest="ukb_phe")
    parser.add_argument("--Y", default=",".join(DEFAULT_Y))
    parser.add_argument("--config")
    parser.add_argument("--out", required=True)
    parser.add_argument("--data-label", choices=("REAL", "SIMULATED"), default="REAL")
    args = parser.parse_args(argv)
    cfg = load_config(args.config)
    outcomes = csv_words(args.Y)
    if not outcomes or any(not re.fullmatch(r"[A-Za-z0-9_.-]+", y) for y in outcomes):
        raise ValueError("Invalid outcome names")
    source = Path(args.cohort) if args.cohort else Path(args.ukb_phe)
    if not args.cohort and source.is_dir():
        source = source / "Rdata" / "all.rds"
    data = phenotype_read_source(source)
    inputs = [file_record(source)]
    bbc = None
    absent_bbc = set(cfg["phenotypes"]["bbc_columns"]) - set(data)
    if absent_bbc:
        bbc_path = source.parent / "bbc.rds"
        if not bbc_path.is_file():
            raise ValueError(f"Requested BBC variables absent from source and no {bbc_path}")
        bbc = phenotype_read_source(bbc_path)
        inputs.append(file_record(bbc_path))
    prepared, audit = phenotype_normalize_cohort(data, outcomes, cfg, bbc)
    write_tsv(args.out, prepared)
    write_tsv(str(args.out) + ".audit.tsv", audit)
    write_json(str(args.out) + ".manifest.json", {
        "version": VERSION, "stage": "phenotypes", "data_label": args.data_label,
        "status": "completed", "outcomes": outcomes, "config": cfg,
        "inputs": inputs, "output": file_record(args.out), "n": len(prepared),
        "sample_id_hash": sample_hash(prepared.sample_id),
        "columns": list(prepared),
        "missing_values": "Preserved; no mean/zero replacement",
        "date_contract": "LE8 same-day prevalent; onset in (baseline,censor] incident; positive undated disease excluded",
        "time_scale": "Years since baseline; cause-specific Cox with death/loss/admin censoring",
        "source_code": [file_record(__file__), file_record(Path(__file__).with_name('burden.py'))],
    })
    print(f"[ASSOC] phenotype bridge: {len(prepared)} participants, {len(outcomes)} outcomes; {args.data_label}")
    return 0


def _entry_assoc_phenotypes():
    try:
        raise SystemExit(phenotype_main())
    except (ValueError, RuntimeError, OSError) as phenotype_exc:
        print(f"ERROR: {phenotype_exc}", file=sys.stderr)
        raise SystemExit(2)


# ========================================================================
# assoc models
# ========================================================================
_description_assoc_models = 'Audited individual extra-allele burden and bounded multi-allelic association.\n\nPrimary: incident, baseline disease-free Cox HR per complete-case SD of\nextra_copy_rate. Secondary: multi_het_rate, declared quantiles, logistic case\nstatus, and optional covariate enhancement. No mixed-model claim is made:\nprovide an externally verified unrelated keep list or explicitly choose\nfixed_effect_exploratory. Every planned contrast, including failures, is saved.'
#!/usr/bin/env python3
"""Audited individual extra-allele burden and bounded multi-allelic association.

Primary: incident, baseline disease-free Cox HR per complete-case SD of
extra_copy_rate. Secondary: multi_het_rate, declared quantiles, logistic case
status, and optional covariate enhancement. No mixed-model claim is made:
provide an externally verified unrelated keep list or explicitly choose
fixed_effect_exploratory. Every planned contrast, including failures, is saved.
"""
import warnings
from scipy import optimize, stats


def model_typed_header(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        try:
            header = next(csv.reader(handle, delimiter="\t"))
        except StopIteration as exc:
            raise ValueError(f"Empty table: {path}") from exc
    if len(header) != len(set(header)):
        raise ValueError(f"Duplicate columns in {path}")
    return header


def model_read_exposures(path):
    """Keep IDs as text while avoiding an object allocation for every number.

    Full-cohort RAM has not been benchmarked by the synthetic verification.
    This reader makes numeric counts/rates native arrays and repeated set/cohort
    labels categorical; it does not load a 20-million-row matrix as all strings.
    """
    header = model_typed_header(path)
    types = {"sample_id": "string", "set_id": "category", "catalog_id": "category",
             "cohort_id": "category", "build": "category", "source_scope": "category", "eligibility_reason": "category",
             "n_callable": "int64", "n_targets": "int64", "eligible": "int8"}
    for col in ("n_missing", "extra_copies", "extra_carrier_sites", "heterozygous_sites",
                "extra_copy_rate", "extra_carrier_rate", "multi_het_rate", "call_rate", "missing_fraction"):
        types[col] = "float64"
    return pd.read_csv(path, sep="\t", dtype={k: v for k, v in types.items() if k in header},
                       keep_default_na=False, na_values=["", ".", "NA", "NaN", "nan"], low_memory=False)


def model_read_unrelated_keep(path):
    """Read a one-column sample_id/eid list or a PLINK FID/IID keep file."""
    opener = __import__("gzip").open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        lines = [line.strip().split() for line in handle if line.strip()]
    if not lines:
        raise ValueError("Unrelated keep list is empty")
    header = [x.lstrip("#").lower() for x in lines[0]]
    candidates = [i for i, col in enumerate(header) if col in ("sample_id", "eid", "iid")]
    if candidates:
        column = candidates[-1]
        lines = lines[1:]
    else:
        column = 1 if len(lines[0]) == 2 else 0
        if len(lines[0]) > 2:
            raise ValueError("Unrelated keep file needs sample_id or a FID/IID header")
    if any(len(line) <= column for line in lines):
        raise ValueError("Malformed unrelated keep file")
    ids = validate_ids(pd.Series([line[column] for line in lines]), "unrelated keep")
    if not len(ids):
        raise ValueError("Unrelated keep contains no IDs")
    return set(ids)


def model_validate_exposures(data):
    # Reassignment below never mutates a shared column in place; a shallow
    # frame avoids duplicating every numeric array in a large exposure table.
    data = data.copy(deep=False)
    need = {"sample_id", "set_id", "n_callable", "n_targets"}
    if need - set(data):
        raise ValueError(f"Missing exposure columns: {sorted(need - set(data))}")
    if isinstance(data.sample_id.dtype, pd.StringDtype):
        validate_ids(data.sample_id.drop_duplicates(), "unique exposure IDs")
    else:
        data["sample_id"] = validate_ids(data.sample_id, "exposures", unique=False)
    if data.set_id.isna().any() or data.set_id.astype(str).str.strip().eq("").any():
        raise ValueError("Missing set_id in exposures")
    if data.duplicated(["sample_id", "set_id"]).any():
        raise ValueError("Duplicate sample_id/set_id exposure records")
    for col in ("cohort_id", "build", "catalog_id"):
        if col in data and data[col].dropna().nunique() > 1:
            raise ValueError(f"Exposure file mixes {col}; analyze explicitly separated cohorts/catalogues")
    for col in ("n_callable", "n_targets"):
        data[col] = numeric(data[col], col, allow_missing=False)
        if ((data[col] < 0) | (data[col] != np.floor(data[col]))).any():
            raise ValueError(f"{col} must contain nonnegative integer locus counts")
    if data.n_callable.gt(data.n_targets).any():
        raise ValueError("n_callable exceeds n_targets")
    rate = data.n_callable / data.n_targets.replace(0, np.nan)
    if "call_rate" in data:
        supplied = numeric(data.call_rate, "call_rate")
        comparable = supplied.notna() & rate.notna()
        if not np.allclose(supplied[comparable], rate[comparable], atol=1e-8, rtol=1e-6):
            raise ValueError("call_rate disagrees with callable/target locus counts")
    data["call_rate"] = rate
    for rate_col, count_col, upper in (("extra_copy_rate", "extra_copies", 2),
                                       ("extra_carrier_rate", "extra_carrier_sites", 1),
                                       ("multi_het_rate", "heterozygous_sites", 1)):
        if rate_col not in data:
            continue
        data[rate_col] = numeric(data[rate_col], rate_col)
        r = data[rate_col]
        if ((r.dropna() < 0) | (r.dropna() > upper)).any():
            raise ValueError(f"{rate_col} outside [0,{upper}]")
        if (data.n_callable.eq(0) & r.notna()).any():
            raise ValueError("Zero callable loci must have NA exposure rates, never zero")
        if count_col in data:
            count = numeric(data[count_col], count_col, allow_missing=False)
            expected = count / data.n_callable.replace(0, np.nan)
            valid = r.notna() & expected.notna()
            if ((count < 0) | (count > upper * data.n_callable)).any() or not np.allclose(
                    r[valid], expected[valid], atol=1e-8, rtol=1e-6):
                raise ValueError(f"{rate_col} disagrees with {count_col}/n_callable")
    if "eligible" in data:
        data["eligible"] = numeric(data.eligible, "eligible", allow_missing=False)
        if not data.eligible.isin([0, 1]).all():
            raise ValueError("Exposure eligible must be 0/1")
    else:
        data["eligible"] = 1
    return data


def model_encode_covariates(data, columns, categorical, intercept=False):
    arrays, names, dropped, references = [], [], [], {}
    if intercept:
        arrays.append(np.ones(len(data)))
        names.append("intercept")
    for col in columns:
        if col in categorical:
            values = data[col].astype(str)
            levels = sorted(values.unique())
            if len(levels) <= 1:
                dropped.append(col)
                continue
            references[col] = levels[0]
            for level in levels[1:]:
                arrays.append(values.eq(level).to_numpy(float))
                names.append(f"{col}={level}")
        else:
            values = numeric(data[col], col, allow_missing=False).to_numpy(float)
            sd = float(np.std(values, ddof=1))
            if sd < 1e-12:
                dropped.append(col)
                continue
            arrays.append((values - values.mean()) / sd)
            names.append(col)
    matrix = np.column_stack(arrays) if arrays else np.empty((len(data), 0))
    return matrix, names, dropped, references


def model_fit_array_model(x, event, time, endpoint, settings):
    """Fit calibrated-to-the-declared-fixed-effect-design MLE, or abstain.

    Wald inference is asymptotic. Support gates do not turn it into SAIGE, a
    kinship mixed model, Firth regression or a causal estimate.
    """
    import statsmodels.api as sm
    from statsmodels.duration.hazard_regression import PHReg
    from statsmodels.tools.sm_exceptions import PerfectSeparationError

    x, y = np.asarray(x, float), np.asarray(event, float)
    n, k = x.shape
    events = int(y.sum())
    base = {"status": "not_estimable", "reason": "", "n": n, "n_events": events,
            "n_controls": int(n - events), "n_parameters": k, "warnings": ""}
    if n < settings["min_n"]:
        return dict(base, reason="few_complete_cases")
    if events < settings["min_events"] or n - events < settings["min_controls"]:
        return dict(base, reason="few_events_or_controls")
    information = min(events, n - events) if endpoint == "logistic" else events
    if information < settings["min_events_per_parameter"] * k:
        return dict(base, reason="events_per_parameter_gate")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        return dict(base, reason="nonfinite_design")
    rank = int(np.linalg.matrix_rank(x))
    condition = float(np.linalg.cond(x))
    base.update(design_rank=rank, condition_number=condition)
    if rank < k or not np.isfinite(condition) or condition > settings["condition_limit"]:
        return dict(base, reason="rank_deficient_or_ill_conditioned_design")
    if endpoint == "logistic" and n <= 20000 and k <= 50:
        # A feasible separating hyperplane proves MLE does not exist. This
        # bounded LP is an additional check, never a substitute for convergence.
        separated = optimize.linprog(np.zeros(k), A_ub=-x * (2 * y - 1)[:, None],
                                     b_ub=-np.ones(n), bounds=[(None, None)] * k,
                                     method="highs")
        if separated.success:
            return dict(base, reason="complete_separation")
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            if endpoint == "cox":
                model = PHReg(np.asarray(time, float), x, status=y,
                              ties=settings["ties"], missing="raise")
                result = model.fit(maxiter=settings["max_iter"], disp=0)
                converged = np.max(np.abs(model.score(result.params))) / max(n, 1) < 1e-6
            else:
                model = sm.GLM(y, x, family=sm.families.Binomial(), missing="raise")
                result = model.fit(maxiter=settings["max_iter"], tol=1e-9)
                converged = bool(result.converged)
            beta = np.asarray(result.params)
            covariance = np.asarray(result.cov_params())
        messages = [str(w.message) for w in caught]
        base["warnings"] = "; ".join(dict.fromkeys(messages))
        warning_failure = any(any(word in text.lower() for word in
                                  ("separation", "failed to converge", "overflow", "invalid value", "singular"))
                              for text in messages)
        if warning_failure or not converged:
            return dict(base, reason="fit_warning_or_nonconvergence")
        if not np.isfinite(beta).all() or not np.isfinite(covariance).all():
            return dict(base, reason="nonfinite_fit")
        eig = np.linalg.eigvalsh(covariance)
        if eig.min() <= 0 or eig.max() / eig.min() > settings["condition_limit"]**2:
            return dict(base, reason="singular_information_matrix")
        se = np.sqrt(np.diag(covariance))
        if (np.abs(beta) > 30).any() or (se > 30).any():
            return dict(base, reason="unstable_large_coefficients_or_se")
        if endpoint == "logistic":
            eta = x @ beta
            if eta[y == 1].min() >= eta[y == 0].max() - 1e-10:
                return dict(base, reason="complete_or_quasi_separation")
        return dict(base, status="estimated", beta_vector=beta, covariance=covariance,
                    log_likelihood=float(result.llf), fit_method="PHReg_fixed_effect" if endpoint == "cox" else "binomial_GLM_MLE")
    except (ValueError, np.linalg.LinAlgError, FloatingPointError, ZeroDivisionError, PerfectSeparationError) as exc:
        return dict(base, reason="fit_failure:" + type(exc).__name__, warnings=str(exc))


def model_contrast_result(fit, vector=None, indices=None):
    out = {k: v for k, v in fit.items() if k not in ("beta_vector", "covariance")}
    out.update(beta=np.nan, se=np.nan, ratio=np.nan, ci_low=np.nan, ci_high=np.nan, p=np.nan, df=1)
    if fit["status"] != "estimated":
        return out
    b, v = fit["beta_vector"], fit["covariance"]
    if indices is not None:
        ii = list(indices)
        part = v[np.ix_(ii, ii)]
        rank = np.linalg.matrix_rank(part)
        if rank != len(ii):
            return dict(out, status="not_estimable", reason="joint_covariance_singular")
        value = float(b[ii] @ np.linalg.solve(part, b[ii]))
        return dict(out, p=float(stats.chi2.sf(max(value, 0), len(ii))), df=len(ii),
                    joint_statistic=value)
    delta = float(np.asarray(vector) @ b)
    variance = float(np.asarray(vector) @ v @ np.asarray(vector))
    if variance <= 0 or not np.isfinite(variance):
        return dict(out, status="not_estimable", reason="contrast_variance_invalid")
    se = math.sqrt(variance)
    low, high = delta - 1.959963984540054 * se, delta + 1.959963984540054 * se
    if max(abs(low), abs(high)) > 700:
        return dict(out, status="not_estimable", reason="contrast_interval_overflow")
    return dict(out, beta=delta, se=se, ratio=math.exp(delta), ci_low=math.exp(low), ci_high=math.exp(high),
                p=float(2 * stats.norm.sf(abs(delta / se))))


def model_base_record(outcome, set_id, exposure, endpoint, adjustment, representation, contrast,
                settings, data_label, inference_scope, locus_id=""):
    key = [outcome, str(set_id), locus_id, exposure, endpoint, adjustment, representation, contrast]
    family = "primary" if (endpoint == "cox" and adjustment == "standard" and
                            representation == "continuous" and exposure == settings["primary_exposure"] and not locus_id) else (
                                "locus" if locus_id else "sensitivity")
    return {"test_id": hashlib.sha256(json.dumps(key).encode()).hexdigest()[:24],
            "outcome": outcome, "set_id": set_id, "locus_id": locus_id,
            "exposure": exposure, "endpoint": endpoint, "adjustment": adjustment,
            "representation": representation, "contrast": contrast, "family": family,
            "data_label": data_label, "inference_scope": inference_scope,
            "effect_measure": "HR" if endpoint == "cox" else "OR",
            "case_definition": "incident" if endpoint == "cox" else settings["logistic_kind"],
            "unit": "per_complete_case_SD" if representation == "continuous" else (
                "per_allele_copy" if representation == "locus" else "declared_quantile_contrast"),
            "status": "not_estimable", "reason": "", "p": np.nan,
            "beta": np.nan, "se": np.nan, "ratio": np.nan, "ci_low": np.nan, "ci_high": np.nan,
            "n": 0, "n_events": 0, "n_controls": 0}


def model_endpoint_columns(endpoint, settings):
    if endpoint == "cox":
        return "__event", "__eligible_cox"
    kind = settings["logistic_kind"]
    if kind == "explicit":
        return "__explicit_case", "__eligible_explicit"
    return "__prevalent_case" if kind == "prevalent" else "__ever_case", "__eligible_" + kind


def model_prepare_frame(data, exposures, endpoint, adjustment, settings, keep_ids):
    """Complete-case filtering is per declared model and is fully counted."""
    covars = list(dict.fromkeys(settings["covariates"] + settings["technical_covariates"] + (
        settings["enhanced_covariates"] if adjustment == "enhanced" else [])))
    audit = {"n_input": len(data), "n_unrelated": 0, "n_outcome_eligible": 0,
             "n_exposure_qc": 0, "n_complete": 0, "covariates": ";".join(covars)}
    if keep_ids is None and settings["relatedness_policy"] == "require_unrelated":
        return None, audit, "unrelated_keep_required"
    data = data.copy()
    if keep_ids is not None:
        data = data.loc[data.sample_id.isin(keep_ids)].copy()
    audit["n_unrelated"] = len(data)
    ycol, ecol = model_endpoint_columns(endpoint, settings)
    data = data.loc[data[ecol]].copy()
    audit["n_outcome_eligible"] = len(data)
    data = data.loc[data.eligible.eq(1) & data.n_callable.gt(0) & data.call_rate.ge(settings["min_call_rate"])].copy()
    audit["n_exposure_qc"] = len(data)
    absent = [c for c in list(exposures) + covars if c not in data]
    if absent:
        return None, audit, "missing_covariates_or_exposure:" + ",".join(absent)
    need = list(dict.fromkeys(list(exposures) + covars + [ycol] + (["__time"] if endpoint == "cox" else [])))
    try:
        for col in need:
            if col not in settings["categorical"]:
                data[col] = numeric(data[col], col)
    except ValueError as exc:
        return None, audit, str(exc)
    data = data.loc[data[need].notna().all(axis=1)].copy()
    audit["n_complete"] = len(data)
    audit["sample_id_hash"] = sample_hash(data.sample_id)
    audit["n_missing_covariate_or_exposure"] = audit["n_exposure_qc"] - len(data)
    if not len(data):
        return None, audit, "no_complete_eligible_participants"
    cov, names, dropped, refs = model_encode_covariates(data, covars, settings["categorical"], endpoint == "logistic")
    audit.update(constant_covariates=";".join(dropped), category_references=json.dumps(refs, sort_keys=True),
                 n_events=int(data[ycol].sum()), n_controls=int(len(data) - data[ycol].sum()))
    return (data, cov, names, ycol), audit, ""


def model_burden_tests(data, outcome, set_id, exposure, endpoint, adjustment, settings, keep_ids,
                 data_label, inference_scope):
    prepared, audit, reason = model_prepare_frame(data, [exposure], endpoint, adjustment, settings, keep_ids)
    common = dict(outcome=outcome, set_id=set_id, exposure=exposure, endpoint=endpoint,
                  adjustment=adjustment, settings=settings, data_label=data_label, inference_scope=inference_scope)
    records, flows = [], []
    representations = ["continuous"] + (["quantile"] if settings["quantiles"] else [])
    for rep in representations:
        contrasts = ["per_SD"] if rep == "continuous" else [f"Q{k}_vs_Q1" for k in range(2, settings["quantiles"] + 1)] + ["joint_quantiles"]
        rows = [model_base_record(**common, representation=rep, contrast=c) for c in contrasts]
        for row in rows:
            row.update({k: v for k, v in audit.items() if k not in ("n_input", "n_unrelated", "n_outcome_eligible", "n_exposure_qc")})
            row["n"] = audit.get("n_complete", 0)
        flow = dict(outcome=outcome, set_id=set_id, exposure=exposure, endpoint=endpoint,
                    adjustment=adjustment, representation=rep, **audit)
        local_reason = reason
        if prepared is not None:
            frame, cov, covnames, ycol = prepared
            xraw = frame[exposure].to_numpy(float)
            mean, sd = float(np.mean(xraw)), float(np.std(xraw, ddof=1)) if len(xraw) > 1 else 0
            for row in rows:
                row.update(exposure_mean=mean, exposure_sd=sd)
            if not np.isfinite(sd) or sd < 1e-12:
                local_reason = "constant_exposure"
            elif rep == "continuous":
                matrix = ((xraw - mean) / sd)[:, None]
                labels = ["exposure_per_SD"]
            else:
                q = settings["quantiles"]
                edges = np.quantile(xraw, np.linspace(0, 1, q + 1))
                for row in rows:
                    row["quantile_edges"] = json.dumps(edges.tolist())
                if len(np.unique(edges)) != q + 1:
                    local_reason = "nonunique_quantile_edges_no_rank_jitter"
                else:
                    groups = np.searchsorted(edges[1:-1], xraw, side="left")
                    group_counts = []
                    for g in range(q):
                        use = groups == g
                        ne = int(frame.loc[use, ycol].sum())
                        group_counts.append({"group": f"Q{g+1}", "n": int(use.sum()), "events": ne})
                        if use.sum() < settings["min_group_n"] or ne < settings["min_group_events"] or (
                                endpoint == "logistic" and use.sum() - ne < settings["min_group_events"]):
                            local_reason = "quantile_group_support_gate"
                    for row in rows:
                        row["group_counts"] = json.dumps(group_counts)
                    matrix = np.column_stack([groups == g for g in range(1, q)]).astype(float)
                    labels = [f"Q{g+1}" for g in range(1, q)]
            if not local_reason:
                design = np.column_stack([matrix, cov])
                fit = model_fit_array_model(design, frame[ycol], frame["__time"] if endpoint == "cox" else None,
                                      endpoint, settings)
                for k, row in enumerate(rows):
                    if row["contrast"] == "joint_quantiles":
                        effect = model_contrast_result(fit, indices=range(matrix.shape[1]))
                    else:
                        vector = np.zeros(design.shape[1])
                        vector[k] = 1
                        effect = model_contrast_result(fit, vector=vector)
                    row.update(effect, design_columns=";".join(labels + covnames))
        if local_reason:
            for row in rows:
                row["reason"] = local_reason
        flow["status"] = rows[0]["status"]
        flow["reason"] = rows[0]["reason"]
        flows.append(flow)
        records.extend(rows)
    return records, flows


def model_load_locus_dosages(path, settings):
    required = {"sample_id", "site_id", "baseline_base", "allele", "allele_rank", "dosage"}
    header = model_typed_header(path)
    types = {"sample_id": "string", "site_id": "category", "baseline_base": "category",
             "allele": "category", "allele_rank": "float64", "dosage": "float64",
             "allele_role": "category", "measurement_status": "category", "catalog_id": "category",
             "cohort_id": "category", "build": "category", "evaluable": "float64"}
    chunks, sites, count = [], set(), 0
    for chunk in pd.read_csv(path, sep="\t", dtype={k: v for k, v in types.items() if k in header}, keep_default_na=False,
                             na_values=["", ".", "NA", "NaN", "nan"], chunksize=200000):
        if required - set(chunk):
            raise ValueError(f"Locus dosage columns absent: {sorted(required-set(chunk))}")
        sites.update(chunk.site_id.dropna())
        count += len(chunk)
        if len(sites) > settings["max_loci"] or count > settings["max_dosage_rows"]:
            raise ValueError("Locus panel exceeds declared bounds; explicitly select fewer sites or change bounds")
        chunks.append(chunk)
    if not chunks:
        return pd.DataFrame(columns=sorted(required))
    data = pd.concat(chunks, ignore_index=True)
    data["sample_id"] = validate_ids(data.sample_id, "locus dosages", unique=False)
    if data.site_id.isna().any() or data.duplicated(["sample_id", "site_id", "allele_rank"]).any():
        raise ValueError("Missing site ID or duplicate sample/site/allele locus dosage")
    data["allele_rank"] = numeric(data.allele_rank, "allele_rank", allow_missing=False)
    if not data.allele_rank.isin([2, 3, 4]).all():
        raise ValueError("Locus dosages must omit baseline and retain ranks 2/3/4")
    data["allele_rank"] = data.allele_rank.astype(int)
    data["dosage"] = numeric(data.dosage, "locus dosage")
    if ((data.dosage.dropna() < 0) | (data.dosage.dropna() > 2)).any():
        raise ValueError("Locus dosage outside [0,2]")
    return data


def model_locus_tests(dosages, phenotype, outcome, settings, keep_ids, data_label, inference_scope, all_exposures):
    records, flows = [], []
    # Global genetic-QC eligibility remains in force for optional locus analyses.
    if "ALL" in set(all_exposures.set_id):
        global_qc = all_exposures.loc[all_exposures.set_id.eq("ALL"),
                                     ["sample_id", "eligible", "n_callable", "call_rate"]].copy()
    else:
        raise ValueError("Optional locus tests require the ALL exposure set for shared sample QC")
    for site_id, group in dosages.groupby("site_id", sort=True, observed=True):
        ranks = sorted(int(k) for k in group.allele_rank.unique())
        if 2 not in ranks or not any(k >= 3 for k in ranks):
            raise ValueError(f"Locus {site_id} is not an externally defined multi-allelic panel")
        if ranks != list(range(2, max(ranks) + 1)):
            raise ValueError("Locus allele ranks must be contiguous; an omitted alternate cannot be recoded to baseline")
        baseline = group.baseline_base.dropna().unique()
        alleles = {}
        if len(baseline) != 1 or len(str(baseline[0])) != 1 or str(baseline[0]) not in "ACGT":
            raise ValueError("Locus baseline must be frozen across people")
        for k in ranks:
            names = group.loc[group.allele_rank.eq(k), "allele"].dropna().unique()
            if len(names) != 1 or names[0] not in "ACGT" or len(names[0]) != 1:
                raise ValueError("Each locus rank must have one fixed nucleotide")
            alleles[k] = names[0]
        if len(set(alleles.values()) | set(baseline)) != len(ranks) + 1:
            raise ValueError("Duplicate nucleotide labels in locus panel")
        wide = group.pivot(index="sample_id", columns="allele_rank", values="dosage")
        observed = wide.notna().sum(axis=1)
        if ((observed != 0) & (observed != len(ranks))).any():
            raise ValueError("All nonbaseline dosage columns at one locus must be jointly callable or jointly missing")
        if wide.sum(axis=1, min_count=len(ranks)).gt(2 + 1e-8).any():
            raise ValueError("Nonbaseline allele dosages exceed diploid ploidy")
        columns = [f"dosage_rank{k}" for k in ranks]
        wide.columns = columns
        base = phenotype.merge(wide.reset_index(), on="sample_id", how="inner", validate="one_to_one")
        base = base.drop(columns=[x for x in ("eligible", "n_callable", "call_rate") if x in base])
        base = base.merge(global_qc, on="sample_id", how="left", validate="one_to_one")
        contrasts = [(f"rank{k}_vs_baseline", "single", k) for k in ranks]
        contrasts += [("joint_all_nonbaseline", "joint", ranks),
                      ("joint_extra_rank3_rank4", "joint", [k for k in ranks if k >= 3])]
        contrasts += [(f"rank{k}_vs_rank2", "difference", k) for k in ranks if k >= 3]
        for endpoint in settings["endpoints"]:
            for adjustment in settings["adjustments"]:
                prepared, audit, reason = model_prepare_frame(base, columns, endpoint, adjustment, settings, keep_ids)
                rows = [model_base_record(outcome, "LOCUS", "multiallelic_dosages", endpoint, adjustment,
                                    "locus", c, settings, data_label, inference_scope, str(site_id))
                        for c, _, _ in contrasts]
                for row in rows:
                    row.update(audit, n=audit.get("n_complete", 0), baseline_base=baseline[0],
                               allele_map=json.dumps(alleles, sort_keys=True))
                if prepared is not None:
                    frame, cov, covnames, ycol = prepared
                    x = frame[columns].to_numpy(float)
                    ac = x.sum(axis=0)
                    mac = np.minimum(ac, 2 * len(frame) - ac)
                    if np.any(mac < settings["min_locus_mac"]):
                        reason = "locus_minor_allele_count_gate"
                    elif np.any(np.std(x, axis=0) < 1e-12):
                        reason = "constant_locus_dosage"
                    if not reason:
                        design = np.column_stack([x, cov])
                        fit = model_fit_array_model(design, frame[ycol], frame["__time"] if endpoint == "cox" else None,
                                              endpoint, settings)
                        for row, (_, kind, target) in zip(rows, contrasts):
                            if kind == "joint":
                                effect = model_contrast_result(fit, indices=[ranks.index(k) for k in target])
                            else:
                                vector = np.zeros(design.shape[1])
                                vector[ranks.index(target)] = 1
                                if kind == "difference":
                                    vector[ranks.index(2)] = -1
                                effect = model_contrast_result(fit, vector=vector)
                            row.update(effect, design_columns=";".join(columns + covnames),
                                       allele_counts=json.dumps(dict(zip(columns, ac.tolist()))))
                if reason:
                    for row in rows:
                        row["reason"] = reason
                records.extend(rows)
                flows.append(dict(outcome=outcome, locus_id=site_id, endpoint=endpoint, adjustment=adjustment,
                                  representation="locus", status=rows[0]["status"], reason=rows[0]["reason"], **audit))
    return records, flows


def model_run_models(phenotype, exposures, outcome, config, keep_ids=None, data_label="REAL", dosages=None):
    settings = config["models"]
    phenotype = phenotype.copy()
    phenotype["sample_id"] = validate_ids(phenotype.sample_id, "prepared phenotypes")
    exposures = model_validate_exposures(exposures)
    derived, outcome_audit = phenotype_derive_outcome(phenotype, outcome, config)
    for col in derived:
        phenotype["__" + col] = derived[col]
    phenotype["__eligible_explicit"] = phenotype["__explicit_case"].notna()
    scope = "fixed_effect_unrelated_keep_supplied" if keep_ids is not None else "fixed_effect_exploratory_relatedness_unresolved"
    overlap = set(phenotype.sample_id) & set(exposures.sample_id)
    if not overlap:
        raise ValueError("No matching sample IDs between phenotypes and exposures; row-order joins are prohibited")
    records, flows = [], []
    for set_id, sub in exposures.groupby("set_id", sort=True, observed=True):
        collision = set(sub) & set(phenotype) - {"sample_id"}
        # Exposure-derived callability has one unambiguous provenance. Other
        # collisions can obscure an accidentally reused phenotype/exposure.
        permitted = {"call_rate", "eligible", "n_callable", "n_targets"}
        if collision - permitted:
            raise ValueError(f"Ambiguous columns in phenotype/exposure join: {sorted(collision-permitted)}")
        p = phenotype.drop(columns=list(collision))
        data = p.merge(sub, on="sample_id", how="inner", validate="one_to_one", sort=False)
        for exposure in settings["exposures"]:
            for endpoint in settings["endpoints"]:
                for adjustment in settings["adjustments"]:
                    rows, flow = model_burden_tests(data, outcome, set_id, exposure, endpoint, adjustment,
                                              settings, keep_ids, data_label, scope)
                    records.extend(rows)
                    flows.extend(flow)
    if dosages is not None and not dosages.empty:
        extra, extra_flow = model_locus_tests(dosages, phenotype, outcome, settings, keep_ids,
                                        data_label, scope, exposures)
        records.extend(extra)
        flows.extend(extra_flow)
    results = pd.DataFrame(records)
    if results.test_id.duplicated().any():
        raise ValueError("Duplicate test identities")
    results, families = adjust_all_tests(results)
    exposure_metadata = {key: sorted(str(x) for x in exposures[key].dropna().unique()) if key in exposures else ["not_recorded"]
                         for key in ("source_scope", "catalog_id", "cohort_id", "build")}
    set_sizes = [{"set_id": str(name), "n_samples": int(group.sample_id.nunique()),
                  "n_targets_min": int(group.n_targets.min()), "n_targets_max": int(group.n_targets.max())}
                 for name, group in exposures.groupby("set_id", sort=True, observed=True)]
    audit = {"phenotype_n": len(phenotype), "exposure_n_unique": exposures.sample_id.nunique(),
             "matched_ids": len(overlap), "phenotype_without_exposure": len(set(phenotype.sample_id) - overlap),
             "exposure_without_phenotype": len(set(exposures.sample_id) - overlap),
             "outcome": outcome_audit, "inference_scope": scope,
             "exposure_metadata": exposure_metadata, "catalogue_set_sizes": set_sizes,
             "ALL_definition": "ALL means all supplied external catalogue loci, not necessarily the whole genome",
             "n_planned_tests": len(results), "n_estimable": int(results.status.eq("estimated").sum()),
             "limitations": [
                 "Fixed-effect models; no kinship mixed-model correction is claimed",
                 "Complete-case sample and exposure SD can differ across adjustment models; see sample hashes",
                 "Cause-specific Cox censors death; proportional-hazards and censoring assumptions are not established by this scan",
                 "Allele-rank burden is not a causal intervention or the inverse of ROH",
                 "Default call-rate adjustment alone does not establish depth/batch independence; declare additional technical covariates",
                 "Genetic ancestry PCs and an ancestry label do not eliminate all population structure",
             ]}
    if not settings["technical_covariates"]:
        audit["limitations"].append("No technical covariates were requested: sequencing-confounded exploratory association")
    if keep_ids is None:
        audit["limitations"].append("Relatedness unresolved: inference is explicitly exploratory, not confirmatory")
    return results, pd.DataFrame(flows), families, audit


def model_main(argv=None):
    parser = argparse.ArgumentParser(description=_description_assoc_models, allow_abbrev=False)
    parser.add_argument("--phenotypes", "--cohort", dest="phenotypes", required=True)
    parser.add_argument("--exposures", required=True)
    parser.add_argument("--Y", required=True)
    parser.add_argument("--config")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--unrelated-keep")
    parser.add_argument("--locus-dosages")
    parser.add_argument("--data-label", choices=("REAL", "SIMULATED"), default="REAL")
    args = parser.parse_args(argv)
    if len(csv_words(args.Y)) != 1:
        raise ValueError("Each model stage takes exactly one outcome; report stage combines all outcomes")
    cfg = load_config(args.config)
    input_manifest = Path(str(args.phenotypes) + ".manifest.json")
    if input_manifest.exists():
        pm = json.loads(input_manifest.read_text())
        if pm.get("data_label") != args.data_label:
            raise ValueError("Phenotype data label disagrees with REAL/SIMULATED model request")
        if pm.get("output", {}).get("sha256") != file_record(args.phenotypes)["sha256"]:
            raise ValueError("Phenotype contents changed since the bridge manifest")
    data = read_tsv(args.phenotypes)
    if "sample_id" not in data:
        raise ValueError('Prepared phenotypes require sample_id; run burden.py first')
    exposure = model_read_exposures(args.exposures)
    keep = model_read_unrelated_keep(args.unrelated_keep) if args.unrelated_keep else None
    dosages = model_load_locus_dosages(args.locus_dosages, cfg["models"]) if args.locus_dosages else None
    results, flow, families, audit = model_run_models(data, exposure, args.Y, cfg, keep, args.data_label, dosages)
    out = Path(args.out_dir)
    write_tsv(out / "association_results.tsv", results)
    write_tsv(out / "cohort_flow.tsv", flow)
    write_tsv(out / "analysis_families.tsv", families)
    import statsmodels
    inputs = [file_record(args.phenotypes), file_record(args.exposures)]
    inputs.extend(file_record(p) for p in (args.unrelated_keep, args.locus_dosages, args.config) if p)
    manifest = {"version": VERSION, "stage": "models", "outcome": args.Y, "data_label": args.data_label,
                "status": "completed" if audit["n_estimable"] else "completed_without_estimable_tests",
                "config": cfg, "inputs": inputs, "audit": audit,
                "statsmodels_version": statsmodels.__version__, "p_value_sidedness": "two_sided",
                "multiplicity_scope": "Per-outcome interim; report recomputes across every requested outcome and failed planned test",
                "outputs": [file_record(out / p) for p in ("association_results.tsv", "cohort_flow.tsv", "analysis_families.tsv")],
                "source_code": [file_record(__file__), file_record(Path(__file__).with_name('burden.py')),
                                file_record(Path(__file__).with_name('burden.py'))]}
    write_json(out / "model_manifest.json", manifest)
    print(f"[ASSOC] {args.Y}: {audit['n_estimable']}/{len(results)} estimable planned tests; {args.data_label}")
    return 0


def _entry_assoc_models():
    try:
        raise SystemExit(model_main())
    except (ValueError, RuntimeError, OSError) as model_exc:
        print(f"ERROR: {model_exc}", file=sys.stderr)
        raise SystemExit(2)


# ========================================================================
# assoc report
# ========================================================================
_description_assoc_report = 'Combine all requested outcomes before BH/BY; retain both effect directions.\n\nCreates a portable aggregate HTML report and publication-exportable forest\nplots. Individual sample IDs never enter the report. SIMULATED is conspicuous.'
#!/usr/bin/env python3
"""Combine all requested outcomes before BH/BY; retain both effect directions.

Creates a portable aggregate HTML report and publication-exportable forest
plots. Individual sample IDs never enter the report. SIMULATED is conspicuous.
"""
import html


def report_gather_results(root, outcomes, out_dir):
    root, out_dir = Path(root).resolve(), Path(out_dir).resolve()
    paths = sorted(p for p in root.rglob("association_results.tsv")
                   if out_dir not in p.resolve().parents and p.parent != out_dir)
    frames, manifests, records = [], [], []
    for path in paths:
        data = read_tsv(path)
        if "outcome" not in data or not set(data.outcome.dropna()) & set(outcomes):
            continue
        manifest_path = path.parent / "model_manifest.json"
        if not manifest_path.is_file():
            raise ValueError(f"Missing model manifest for {path}; do not combine unverified loose result tables")
        manifest = json.loads(manifest_path.read_text())
        expected = [x for x in manifest.get("outputs", []) if Path(x["path"]).name == path.name]
        observed = file_record(path)
        if len(expected) != 1 or expected[0].get("sha256") != observed["sha256"]:
            raise ValueError(f"Model result changed after its manifest: {path}")
        data = data.loc[data.outcome.isin(outcomes)].copy()
        frames.append(data)
        manifests.append(manifest)
        records.extend([observed, file_record(manifest_path)])
    if not frames:
        raise ValueError("No complete model outputs for the requested outcomes")
    data = pd.concat(frames, ignore_index=True)
    missing = set(outcomes) - set(data.outcome)
    if missing:
        raise ValueError(f"Missing outcome result families: {sorted(missing)}; do not reduce the correction family silently")
    if data.test_id.duplicated().any():
        raise ValueError("Duplicate test identities across result directories; point report at one run, not mixed reruns")
    labels = set(data.data_label) | {m.get("data_label") for m in manifests}
    if len(labels) != 1 or next(iter(labels)) not in ("REAL", "SIMULATED"):
        raise ValueError("REAL and SIMULATED outputs cannot be mixed")
    # Discard interim per-Y q-values and recompute the full declared family.
    data = data.drop(columns=[c for c in data if c.startswith("q_bh") or c.startswith("q_by")])
    for col in ("p", "ratio", "ci_low", "ci_high", "beta", "se", "n", "n_events", "n_controls"):
        if col in data:
            data[col] = pd.to_numeric(data[col], errors="coerce")
    data, families = adjust_all_tests(data)
    return data, families, manifests, records, next(iter(labels))


def report_forest_plot(data, out_dir, label, max_rows=60):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    # Prespecified primary family first; never select rows by protective sign
    # or statistical significance. Failed rows remain visible in the forest.
    primary = data.loc[data.family.eq("primary")].copy()
    is_primary = not primary.empty
    if primary.empty:
        primary = data.loc[data.representation.eq("continuous")].copy()
    primary = primary.sort_values(["outcome", "set_id", "exposure", "adjustment", "endpoint"])
    pages = max(1, math.ceil(len(primary) / max_rows))
    outputs = []
    for page in range(pages):
        d = primary.iloc[page * max_rows:(page + 1) * max_rows]
        height = max(4.2, 1.9 + .34 * len(d))
        fig, ax = plt.subplots(figsize=(12, height))
        ax.axvline(1, color="#526079", linestyle="--", linewidth=.9)
        labels = []
        valid_intervals = []
        for j, (_, row) in enumerate(d.iterrows()):
            labels.append(f"{row['outcome']}  |  {row['set_id']}" +
                          ("" if is_primary else f"  |  {row['exposure']} / {row['endpoint']}"))
            good = row.status == "estimated" and all(pd.notna(row.get(c)) and row[c] > 0
                                                    for c in ("ratio", "ci_low", "ci_high"))
            if good:
                color = "#1a778a" if row.ratio < 1 else "#ad4d43"
                ax.plot([row.ci_low, row.ci_high], [j, j], color=color, lw=1.35)
                ax.scatter([row.ratio], [j], color=color, s=29, zorder=3)
                valid_intervals.extend([row.ci_low, row.ci_high])
            else:
                ax.text(.02, j, str(row.reason or "not estimable"), transform=ax.get_yaxis_transform(),
                        color="#77818e", fontsize=8, va="center")
        ax.set_yticks(range(len(d)), labels)
        ax.invert_yaxis()
        ax.set_xscale("log")
        if valid_intervals:
            lo = max(min(valid_intervals) / 1.2, 1e-8)
            hi = min(max(valid_intervals) * 1.2, 1e8)
            ax.set_xlim(min(lo, .8), max(hi, 1.2))
        else:
            ax.set_xlim(.5, 2)
        measures = set(d.effect_measure)
        measure_label = "Hazard ratio" if measures == {"HR"} else ("Odds ratio" if measures == {"OR"} else "HR / OR as labelled by endpoint")
        ax.set_xlabel(measure_label + " per complete-case SD (nominal 95% CI)")
        title = "Extra-allele burden: declared primary incident associations" if is_primary else "Declared continuous associations"
        if label == "SIMULATED":
            title = "SIMULATED — " + title
        ax.set_title(title, loc="left", fontsize=13, weight="bold", pad=18)
        ax.grid(axis="x", color="#e4e9ef", linewidth=.65)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="y", length=0, labelsize=8.5)
        fig.text(.02, .012, "Both directions retained. Colors show point-estimate direction, not significance. Full BH/BY results and failed tests are in the table.",
                 fontsize=8, color="#526079")
        fig.tight_layout(rect=(0, .035, 1, 1))
        stem = "forest" if page == 0 else f"forest_{page + 1}"
        for extension in ("png", "svg"):
            path = out / f"{stem}.{extension}"
            fig.savefig(path, dpi=190, facecolor="white", bbox_inches="tight")
            outputs.append(path)
        plt.close(fig)
    return outputs


def report_render_html(data, families, manifests, plots, label, out):
    fields = ["outcome", "set_id", "locus_id", "exposure", "endpoint", "case_definition", "adjustment",
              "representation", "contrast", "family", "status", "reason", "n", "n_events", "n_controls",
              "effect_measure", "ratio", "ci_low", "ci_high", "p", "q_bh_all", "q_by_all",
              "q_bh_family", "q_by_family", "inference_scope", "constant_covariates", "exposure_sd"]
    export = data.reindex(columns=fields).astype(object).where(pd.notna(data.reindex(columns=fields)), None)
    # HTML script escaping avoids </script> in any external annotation text.
    payload = json.dumps(json_safe(export.to_dict("records")), ensure_ascii=False, allow_nan=False).replace("<", "\\u003c")
    limitations = sorted(set(x for m in manifests for x in m.get("audit", {}).get("limitations", [])))
    estimated = int(data.status.eq("estimated").sum())
    q = pd.to_numeric(data.q_by_all, errors="coerce")
    n_by = int((q < .05).sum())
    summaries = []
    for y, g in data.groupby("outcome", sort=False):
        estimates = g.status.eq("estimated")
        summaries.append({"Outcome": y, "Planned contrasts": len(g), "Estimated": int(estimates.sum()),
                          "Not estimable": int((~estimates).sum()),
                          "BY < 0.05 (all tests)": int((pd.to_numeric(g.q_by_all, errors="coerce") < .05).sum())})
    table_summary = pd.DataFrame(summaries).to_html(index=False, border=0, escape=True)
    family_table = families.to_html(index=False, border=0, escape=True)
    scopes = []
    for manifest in manifests:
        audit = manifest.get("audit", {})
        metadata = audit.get("exposure_metadata", {})
        row = {key: ";".join(metadata.get(key, ["not_recorded"]))
               for key in ("source_scope", "catalog_id", "cohort_id", "build")}
        all_sizes = [r for r in audit.get("catalogue_set_sizes", []) if r.get("set_id") == "ALL"]
        row["ALL target loci (min/max)"] = ";".join(f"{r['n_targets_min']}/{r['n_targets_max']}" for r in all_sizes) or "not_recorded"
        if row not in scopes:
            scopes.append(row)
    scopes_html = pd.DataFrame(scopes).to_html(index=False, border=0, escape=True)
    limitations_html = "".join(f"<li>{html.escape(x)}</li>" for x in limitations)
    plot_html = "".join(f'<img src="{html.escape(p.name)}" alt="All prespecified primary effects, both directions" loading="lazy">'
                        for p in plots if p.suffix == ".png")
    note = '<div class="simulation">SIMULATED DATA — software demonstration only. These are not human base_occupancy findings.</div>' if label == "SIMULATED" else ""
    page = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Allele3 individual association report</title><style>
:root{color-scheme:light;--ink:#172b40;--muted:#516579;--line:#dce4eb;--blue:#176c84;--risk:#a44b43}*{box-sizing:border-box}
body{margin:0;background:#f4f7fa;color:var(--ink);font:15px/1.55 system-ui,sans-serif}main{max-width:1450px;margin:auto;padding:34px 30px 70px}
h1{font-size:34px;line-height:1.2;letter-spacing:-.6px;margin:12px 0}h2{font-size:22px;margin:5px 0 16px}p{max-width:1100px}.eyebrow{text-transform:uppercase;color:var(--blue);font-size:12px;font-weight:700;letter-spacing:1.2px}
section{background:white;border:1px solid var(--line);border-radius:12px;padding:24px;margin:22px 0}.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:15px}.card{background:white;border:1px solid var(--line);border-radius:10px;padding:18px}.card b{font-size:30px;display:block}.muted{color:var(--muted)}.simulation{padding:16px 20px;border:2px solid #ad4d43;background:#fff1eb;font-weight:700;margin:18px 0;border-radius:8px}
table{border-collapse:collapse;white-space:nowrap;font-size:12px;width:100%;text-align:left}th,td{border-bottom:1px solid var(--line);padding:9px 11px;max-width:430px;overflow:hidden;text-overflow:ellipsis}th{background:#edf3f8;position:sticky;top:0}tr:hover{background:#f7fafc}.scroll{overflow:auto;max-height:680px}img{width:100%;height:auto}.filters{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:14px 0}input,select,button{font:inherit;border:1px solid #b8c8d5;border-radius:6px;background:white;padding:8px}input{min-width:260px}button{cursor:pointer}a{color:var(--blue)}.downloads{display:flex;gap:20px}.protective{color:var(--blue)}.higher{color:var(--risk)}code{font-size:13px}.status{font-weight:600}@media(max-width:750px){main{padding:18px 12px}.cards{grid-template-columns:1fr}section{padding:16px}h1{font-size:27px}}
</style><main><div class="eyebrow">Allele3 · Individual association module</div>
<h1>Extra alleles and recorded disease outcomes</h1>__NOTE__
<p class="muted">Population-defined extra alleles are externally ranked third and fourth supported nucleotides. A diploid person carries at most two alleles at one ordinary autosomal SNV. This report studies inherited copies and heterozygosity, not “three alleles within one person”.</p>
<div class="cards"><div class="card"><span>All planned contrasts</span><b>__PLANNED__</b></div><div class="card"><span>Estimable contrasts</span><b>__ESTIMATED__</b></div><div class="card"><span>BY &lt; 0.05 across all tests</span><b>__BY__</b></div></div>
<section><h2>Outcome coverage</h2>__SUMMARY__<p class="muted">Failed fits and support-gated tests remain in the table. No direction or significance filter was used to construct this report.</p></section>
<section><h2>Assayed catalogue and cohort</h2>__SCOPES__<p><strong>ALL means all supplied external catalogue loci.</strong> A HERV panel or a selected regional catalogue remains a panel, even when its set is named ALL. It must not be interpreted as whole-genome burden. The external source scope and the actual target-locus counts are shown above.</p></section>
<section><h2>Prespecified continuous effects</h2><p>For incident analyses, baseline and same-day cases are excluded. Death, loss to follow-up and the declared administrative end censor follow-up. Events on the censor date are retained. Positive disease evidence without a diagnosis date is unresolved and excluded from incident and date-derived case-status analyses.</p>__PLOTS__</section>
<section><h2>Every planned association</h2><p>HR and OR answer different questions. A ratio below 1 is an inverse association for that specific outcome; it is not evidence of general health, balancing selection or a causal effect. Intervals are nominal 95% intervals. Both BH and BY corrections are recalculated over every requested outcome.</p>
<div class="downloads"><a href="association_results.tsv">Full results TSV</a><a href="analysis_families.tsv">Test families TSV</a><a href="manifest.json">Provenance manifest</a><a href="forest.svg">Vector forest</a></div>
<div class="filters"><input id="query" aria-label="Search all result fields" placeholder="Search outcome, set, exposure or reason"><select id="status" aria-label="Result status"><option value="">All statuses</option><option>estimated</option><option>not_estimable</option></select><select id="endpoint" aria-label="Endpoint"><option value="">All endpoints</option><option>cox</option><option>logistic</option></select><select id="family" aria-label="Test family"><option value="">All families</option><option>primary</option><option>sensitivity</option><option>locus</option></select></div>
<p id="count" class="muted"></p><div class="scroll"><table><thead id="head"></thead><tbody id="rows"></tbody></table></div><div class="filters"><button id="previous">Previous</button><button id="next">Next</button><span id="page"></span></div></section>
<section><h2>Multiplicity and interpretation</h2>__FAMILIES__<p>Every predeclared failed hypothesis is included in the adjustment denominator as P=1, while its public P and Q entries remain missing. “ALL_ATTEMPTED” includes all outcomes, sets, exposure definitions, adjustment sets, encodings and optional locus contrasts. Family-specific values are secondary views, not replacements for the full correction.</p><ul>__LIMITATIONS__</ul></section>
<footer class="muted">Version __VERSION__. External catalogue ranks are fixed; REF is not assumed to be the major or ancestral allele. No PCA or candidate-window outputs were modified.</footer>
</main><script>
const data=__DATA__, cols=['outcome','set_id','exposure','endpoint','case_definition','adjustment','representation','contrast','status','n','n_events','effect_measure','ratio','ci_low','ci_high','p','q_bh_all','q_by_all','reason'];let page=0,filtered=[];const size=100;
const esc=x=>String(x??'NA').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function fmt(x,k){if(x===null||x===undefined||x==='')return 'NA';if(['p','q_bh_all','q_by_all'].includes(k))return Number(x).toExponential(3);if(['ratio','ci_low','ci_high'].includes(k))return Number(x).toPrecision(4);return esc(x)}
document.getElementById('head').innerHTML='<tr>'+cols.map(x=>'<th>'+esc(x)+'</th>').join('')+'</tr>';
function render(reset=false){if(reset)page=0;const q=document.getElementById('query').value.toLowerCase();filtered=data.filter(d=>(!q||Object.values(d).some(x=>String(x??'').toLowerCase().includes(q)))&&['status','endpoint','family'].every(k=>!document.getElementById(k).value||d[k]===document.getElementById(k).value));page=Math.min(page,Math.max(0,Math.ceil(filtered.length/size)-1));document.getElementById('rows').innerHTML=filtered.slice(page*size,(page+1)*size).map(d=>'<tr>'+cols.map(k=>'<td title="'+esc(d[k])+'"'+(k==='ratio'&&d.ratio!==null?' class="'+(d.ratio<1?'protective':'higher')+'"':'')+'>'+fmt(d[k],k)+'</td>').join('')+'</tr>').join('');document.getElementById('count').textContent=filtered.length+' of '+data.length+' planned results match';document.getElementById('page').textContent='Page '+(page+1)+' of '+Math.max(1,Math.ceil(filtered.length/size));document.getElementById('previous').disabled=page===0;document.getElementById('next').disabled=(page+1)*size>=filtered.length}
['query','status','endpoint','family'].forEach(id=>document.getElementById(id).addEventListener('input',()=>render(true)));document.getElementById('previous').onclick=()=>{page--;render()};document.getElementById('next').onclick=()=>{page++;render()};render();
</script></html>"""
    for key, value in {"NOTE": note, "PLANNED": str(len(data)), "ESTIMATED": str(estimated), "BY": str(n_by),
                       "SUMMARY": table_summary, "SCOPES": scopes_html, "PLOTS": plot_html, "FAMILIES": family_table,
                       "LIMITATIONS": limitations_html, "VERSION": VERSION, "DATA": payload}.items():
        page = page.replace("__" + key + "__", value)
    (Path(out) / "report.html").write_text(page, encoding="utf-8")


def report_main(argv=None):
    parser = argparse.ArgumentParser(description=_description_assoc_report, allow_abbrev=False)
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--Y", default=",".join(DEFAULT_Y))
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args(argv)
    outcomes = csv_words(args.Y)
    data, families, manifests, inputs, label = report_gather_results(args.results_root, outcomes, args.out_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "association_results.tsv", data)
    write_tsv(out / "analysis_families.tsv", families)
    plots = report_forest_plot(data, out, label)
    report_render_html(data, families, manifests, plots, label, out)
    write_json(out / "manifest.json", {
        "version": VERSION, "stage": "report", "status": "completed", "data_label": label,
        "outcomes": outcomes, "n_planned": len(data), "n_estimable": int(data.status.eq("estimated").sum()),
        "multiplicity": "BH and BY over all planned tests across all outcomes; failed tests enter denominator as p=1",
        "inputs": inputs,
        "outputs": [file_record(p) for p in [out / "association_results.tsv", out / "analysis_families.tsv",
                                             out / "report.html", *plots]],
        "source_code": [file_record(__file__), file_record(Path(__file__).with_name('burden.py'))],
    })
    print(f"[ASSOC] report: {len(outcomes)} outcomes, {len(data)} planned tests; {label}")
    return 0


def _entry_assoc_report():
    try:
        raise SystemExit(report_main())
    except (ValueError, RuntimeError, OSError) as report_exc:
        print(f"ERROR: {report_exc}", file=sys.stderr)
        raise SystemExit(2)


# ========================================================================
# assoc demo
# ========================================================================
_description_assoc_demo = 'Generate an explicitly SIMULATED, deterministic, two-chromosome association run.\n\nThis is an engineering/measurement example, never evidence about human disease.\nThe full public workflow is exercised, including genotype QC, chromosome merge,\nLE8 date semantics, unrelated filtering, three outcomes and combined correction.'
#!/usr/bin/env python3
"""Generate an explicitly SIMULATED, deterministic, two-chromosome association run.

This is an engineering/measurement example, never evidence about human disease.
The full public workflow is exercised, including genotype QC, chromosome merge,
LE8 date semantics, unrelated filtering, three outcomes and combined correction.
"""
import shutil
import subprocess


def demo_write_table(path, rows, fields=None):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def demo_make_fixture(out, participants, seed):
    import numpy as np
    import pysam

    fixture = out / "fixtures"
    marker = fixture / "SIMULATED.json"
    if (out / "run").exists() and not marker.is_file():
        raise ValueError("Choose a new demo directory; existing run has no SIMULATED fixture marker")
    if fixture.exists() and any(fixture.iterdir()):
        if not marker.is_file():
            raise ValueError("Refusing to replace a non-demo fixture directory")
        old = json.loads(marker.read_text())
        if old.get("seed") != seed or old.get("participants") != participants:
            raise ValueError("Choose a new --out for a different simulation seed or size")
        if (fixture / "assoc.demo.json").is_file():
            return fixture / "assoc.demo.json"
    fixture.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    ids = [f"SIM{i:07d}" for i in range(1, participants + 1)]
    latent = rng.normal(size=participants)
    n_sites = 96
    # Externally ranked nucleotide identities change among sites. At site zero,
    # REF=G is rank3; at some quad sites REF is rank4. ALT order is deliberately
    # scrambled. This makes accidental ALT-index exposure definitions fail.
    counts, specs = [], []
    extra = np.zeros(participants, dtype=float)
    called = np.zeros(participants, dtype=int)
    for j in range(n_sites):
        chrom = str(1 + j // 48)
        pos = 100010 + (j % 48) * 101
        ranked = list("ACGT")
        ranked = ranked[j % 4:] + ranked[:j % 4]
        quad = j % 4 == 1
        k = 4 if quad else 3
        supported = ranked[:k]
        ref = ranked[2] if j % 3 == 0 else (ranked[3] if quad else ranked[0])
        alts = list(reversed([b for b in supported if b != ref]))
        alleles = [ref] + alts
        source_ac = [1300, 500, 160, 40] if quad else [1300, 500, 200, 0]
        sid = f"SIM_chr{chrom}_{pos}"
        row = {"site_id": sid, "chrom": "chr" + chrom, "pos": pos, "ref": ref,
               "population": "EUR", "build": "GRCh38", "AN": 2000, "call_rate": 1.0}
        row.update({b + "_AC": source_ac[ranked.index(b)] for b in "ACGT"})
        counts.append(row)
        # Shared individual latent variation creates a measurable genome-set
        # burden; no disease is used to select loci or define their ranking.
        p_extra = 1 / (1 + np.exp(-(-2.0 + .7 * latent + .12 * np.sin(j))))
        probs = np.column_stack(((1-p_extra)*.72, (1-p_extra)*.28,
                                 p_extra * (.8 if quad else 1.0),
                                 p_extra * (.2 if quad else 0.0)))
        draws = rng.uniform(size=(participants, 2))
        ranks = (draws[:, :, None] > np.cumsum(probs, axis=1)[:, None, :]).sum(axis=2)
        ranks = np.minimum(ranks, k - 1)
        passing = rng.uniform(size=participants) >= .008
        extra += ((ranks >= 2).sum(axis=1) * passing)
        called += passing
        specs.append((chrom, pos, sid, ref, alts, ranked, ranks, passing))
    exposure = extra / called
    z = (exposure - exposure.mean()) / exposure.std(ddof=0)
    vcf_plain = fixture / "SIMULATED.joint.vcf"
    with vcf_plain.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n##reference=GRCh38\n")
        for chrom in (1, 2):
            handle.write(f"##contig=<ID=chr{chrom},length=1000000>\n")
        handle.write('##FILTER=<ID=PASS,Description="Simulation QC pass">\n')
        handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Diploid native genotype">\n')
        handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        handle.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        handle.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n')
        handle.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + ids) + "\n")
        for chrom, pos, sid, ref, alts, ranked, ranks, passing in specs:
            alleles = [ref] + alts
            values = []
            for pair, valid in zip(ranks, passing):
                if not valid:
                    values.append("./.:0:0:" + ",".join("0" for _ in alleles))
                    continue
                a, b = (alleles.index(ranked[int(rank)]) for rank in pair)
                ad = [0] * len(alleles)
                ad[a] += 20
                ad[b] += 20
                values.append(f"{a}/{b}:40:99:" + ",".join(str(v) for v in ad))
            handle.write("\t".join(["chr" + chrom, str(pos), sid, ref, ",".join(alts),
                                     "100", "PASS", ".", "GT:DP:GQ:AD"] + values) + "\n")
    vcf = fixture / "SIMULATED.joint.vcf.gz"
    pysam.tabix_compress(str(vcf_plain), str(vcf), force=True)
    pysam.tabix_index(str(vcf), preset="vcf", force=True)
    vcf_plain.unlink()
    demo_write_table(fixture / "external_counts.tsv", counts)
    demo_write_table(fixture / "unrelated_keep.tsv", [{"sample_id": sid} for sid in ids])
    demo_write_table(fixture / "dosage_sites.tsv", [{"site_id": row["site_id"]} for row in counts[:2]])
    demo_write_table(fixture / "expected_exposures.tsv", [
        {"sample_id": sid, "n_callable": int(c), "extra_copies": int(e), "extra_copy_rate": float(x)}
        for sid, c, e, x in zip(ids, called, extra, exposure)])
    (fixture / "sets.bed").write_text(
        "chr1\t100000\t102434\tHERV_SIMULATED\n"
        "chr1\t102434\t105000\tMATCHED_SIMULATED\n"
        "chr2\t100000\t105000\tHOTSPOT_SIMULATED\n")
    (fixture / "upstream_qc.txt").write_text(
        "SIMULATED DATA ONLY. Native, unsplit, diploid autosomal VCF with known bases.\n"
        "PASS records; each nonmissing GT has DP=40, GQ=99, balanced called-allele AD.\n"
        "No real read mapping, participant consent, ancestry or biological evidence is implied.\n")

    age = rng.uniform(40, 70, participants)
    sex = rng.integers(0, 2, participants)
    tdi = rng.normal(size=participants)
    pcs = rng.normal(size=(participants, 10))
    baseline = [dt.date(2008, 1, 1) + dt.timedelta(days=int(i)) for i in rng.integers(0, 731, participants)]
    end = dt.date(2023, 4, 1)
    rows = []
    for i, sid in enumerate(ids):
        row = {"sample_id": sid, "age": float(age[i]), "sex": int(sex[i]), "tdi": float(tdi[i]),
               "ethnic.c": "White", "date_attend": str(baseline[i]),
               "birth_date": str(baseline[i] - dt.timedelta(days=int(age[i] * 365.25))),
               "date_lost": "NA", "date_death": "NA",
               "F_ROH": float(rng.uniform(0, .03)), "heterozygosity": float(rng.normal(.29, .02)),
               "other_rare_burden": int(rng.poisson(500))}
        row.update({f"PC{k+1}": float(pcs[i, k]) for k in range(10)})
        # Loss/death are independent of exposure in this engineering fixture.
        if rng.uniform() < .025:
            row["date_lost"] = str(baseline[i] + dt.timedelta(days=int(rng.uniform(1200, 4000))))
        if rng.uniform() < .025:
            row["date_death"] = str(baseline[i] + dt.timedelta(days=int(rng.uniform(1200, 4500))))
        for outcome, beta in (("cvd_cad", -.7), ("ra", .7), ("depression", 0.0)):
            predictor = beta*z[i] + .025*(age[i]-55) + .12*sex[i] + .08*tdi[i]
            prevalence = 1 / (1 + np.exp(-(-.85 + predictor)))
            onset = None
            if rng.uniform() < prevalence:
                onset = baseline[i] - dt.timedelta(days=int(rng.uniform(0, 3000)))
            else:
                years = rng.exponential(1 / (.07 * np.exp(predictor)))
                proposed = baseline[i] + dt.timedelta(days=max(1, int(years * 365.25)))
                if proposed <= end:
                    onset = proposed
            row["fod_icd10_" + outcome] = str(onset) if onset else "NA"
            row["icd10Ct_" + outcome] = int(onset is not None)
        rows.append(row)
    demo_write_table(fixture / "cohort.tsv", rows)
    config = json.loads((Path(__file__).resolve().parents[1] / "config" / "burden.json").read_text())
    config["data_label"] = "SIMULATED"
    config["runtime"].update(python=sys.executable, cores=2, memory_limit_gb=4, hours=1)
    config["catalog"].update(counts=str(fixture / "external_counts.tsv"),
                             source_cohort="SIMULATED_EXTERNAL", source_scope="simulation_target_panel",
                             sets_bed=str(fixture / "sets.bed"))
    config["exposures"].update(vcf=str(vcf), chromosomes="1,2", cohort_id="SIMULATED_UKB_LIKE",
                               upstream_qc=str(fixture / "upstream_qc.txt"),
                               dosage_sites=str(fixture / "dosage_sites.tsv"))
    config["phenotype_source"].update(cohort=str(fixture / "cohort.tsv"), dir_ukb_phe=None,
                                     unrelated_keep=str(fixture / "unrelated_keep.tsv"))
    # Full default PC1..PC10 adjustment and the real disease-model support gates
    # remain in force. The catalogue contains known tri and quad loci.
    config_path = fixture / "assoc.demo.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    marker.write_text(json.dumps({"data_label": "SIMULATED", "seed": seed,
        "participants": participants, "loci": n_sites, "tri_loci": 72, "quad_loci": 24,
        "expected_log_effect_per_sd": {"cvd_cad": -.7, "ra": .7, "depression": 0.0},
        "purpose": "Implementation validation only; no real disease finding"}, indent=2) + "\n")
    return config_path


def demo_verify_exposures(out):
    import numpy as np
    import pandas as pd

    expected = pd.read_csv(out / "fixtures" / "expected_exposures.tsv", sep="\t")
    actual_path = out / "run" / "exposures" / "exposures.tsv.gz"
    if not actual_path.is_file():
        raise ValueError(f"Merged exposures not produced: {actual_path}")
    actual = pd.read_csv(actual_path, sep="\t")
    actual = actual[actual.set_id.eq("ALL")]
    joined = expected.merge(actual, on="sample_id", how="outer", suffixes=("_truth", "_actual"), indicator=True)
    if not joined._merge.eq("both").all():
        raise ValueError("Simulation sample IDs do not match the merged exposure output")
    for name in ("n_callable", "extra_copies", "extra_copy_rate"):
        if not np.allclose(joined[name + "_truth"], joined[name + "_actual"], atol=1e-9, rtol=1e-7):
            raise ValueError(f"Decoded simulation exposure differs from independent nucleotide truth: {name}")
    receipt = {"data_label": "SIMULATED", "status": "passed",
               "participants": len(joined), "checks": ["sample_id_one_to_one",
               "two_chromosome_denominator_sum", "native_extra_copies_including_REF_rank3_and_rank4",
               "known_rate_recovered"], "scope": "Measurement validation; model results are synthetic"}
    (out / "SIMULATED_VERIFICATION.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def demo_main(argv=None):
    parser = argparse.ArgumentParser(description=_description_assoc_demo, allow_abbrev=False)
    parser.add_argument("--out", required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--participants", type=int, default=1800)
    parser.add_argument("--seed", type=int, default=20261106)
    parser.add_argument("--jobs-N", type=int, default=2)
    parser.add_argument("--fixtures-only", action="store_true")
    args = parser.parse_args(argv)
    selected_python = shutil.which(args.python) or args.python
    if Path(selected_python).absolute() != Path(sys.executable).absolute():
        # The launcher/planner need only stdlib; scientific imports happen after
        # choosing the intended environment. No shell interpolation is involved.
        return subprocess.call([selected_python, str(Path(__file__).resolve()),'demo',] + (argv if argv is not None else sys.argv[1:]))
    if args.participants < 1200:
        raise ValueError("Use at least 1200 participants for the default full-adjustment demonstration")
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    config = demo_make_fixture(out, args.participants, args.seed)
    print(f"SIMULATED ONLY: {config}", flush=True)
    if args.fixtures_only:
        return 0
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1")
    command = ["bash", str(Path(__file__).resolve().parents[1] / "allele3.sh"), "burden",
               "--config", str(config), "--python", sys.executable,
               "--out", str(out / "run"), "--submit", "local", "--jobs-N", str(args.jobs_N),
               "--wait", "TRUE", "--data-label", "SIMULATED"]
    result = subprocess.run(command, env=env, check=False)
    if result.returncode:
        return result.returncode
    print(json.dumps(demo_verify_exposures(out), indent=2))
    return 0


def _entry_assoc_demo():
    try:
        raise SystemExit(demo_main())
    except (ValueError, OSError, ImportError) as demo_exc:
        print(f"SIMULATED demo failed: {demo_exc}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == '__main__':
    routes = {'catalog': _entry_assoc_catalog, 'exposure': _entry_assoc_exposure, 'phenotype': _entry_assoc_phenotypes, 'model': _entry_assoc_models, 'report': _entry_assoc_report, 'demo': _entry_assoc_demo}
    arguments = list(sys.argv[1:])
    action = arguments.pop(0) if arguments and arguments[0] in routes else None
    if action is None:
        raise SystemExit('Choose a stage: ' + ', '.join(routes))
    sys.argv[1:] = arguments
    routes[action]()
