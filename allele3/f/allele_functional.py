#!/usr/bin/env python3
"""Compare matched, signed molecular profiles of the four possible SNV bases.

This optional analysis consumes a deliberately explicit raw-profile export. It
does not consume the AVI score, AVI SHAP attributions, maximum-effect summaries,
or the merged splicing download. A feature_id must identify the SAME molecular
readout, gene, tissue, genomic scoring region and track for all alleles at a
site. Its scale and activity threshold must be declared independently of allele
frequencies and outcomes. Every schema feature belongs to the intended panel;
an omitted profile value remains missing, even if absent for all four bases.

Required profile columns:
  chrom pos ref allele feature_id value semantics background_id model_version
Allowed semantics are signed_delta_to_ref or allele_activity. Within one site,
semantics, background_id and model_version must be identical. For signed deltas
only, REF=0 is a mathematical baseline, not a prediction of reference activity.
For allele_activity, an actual REF profile is needed to calculate REF-relative
breadth. All rows for a site must be contiguous; their order need not be sorted.

Required schema columns:
  feature_id modality scale active_threshold provenance
Thresholds apply to abs(delta / scale), not to an AVI attribution or probability.

Optional count columns:
  chrom pos ref population AN A_AC C_AC G_AC T_AC
Counts must describe actual nucleotides in a single common called sample set.
Zero means measured zero; missing any count or AN prohibits frequency ranking.
The numerical checks cannot establish the underlying genotype provenance.

All inputs are declared GRCh38 using --build GRCh38. Exact REF consistency is
checked across inputs; this module does not verify REF against a FASTA. Nothing
here establishes fitness, antagonistic pleiotropy, mutation age or selection.

Example (tab-separated files; the feature label stands for the complete identity):
  profiles.tsv:
    chrom pos ref allele feature_id value semantics background_id model_version
    chr1 100 A C gene=G1;tissue=liver;region=chr1:90-110;track=rna1 1 signed_delta_to_ref hg38_window_1 model_v1
    chr1 100 A G gene=G1;tissue=liver;region=chr1:90-110;track=rna1 -1 signed_delta_to_ref hg38_window_1 model_v1
    chr1 100 A T gene=G1;tissue=liver;region=chr1:90-110;track=rna1 0.2 signed_delta_to_ref hg38_window_1 model_v1
  feature_schema.tsv:
    feature_id modality scale active_threshold provenance
    gene=G1;tissue=liver;region=chr1:90-110;track=rna1 RNA_SEQ 1 0.5 raw_signed_log2FC_fixed_readout
  counts.tsv:
    chrom pos ref population AN A_AC C_AC G_AC T_AC
    chr1 100 A EUR 100 5 80 15 0
Run:
  python allele_functional.py --profiles profiles.tsv --feature-schema feature_schema.tsv \
      --counts counts.tsv --build GRCh38 --min-ac 5 --out functional
The example counts make C rank1, G rank2, A rank3; REF is not the major allele.
These numbers are illustrative only, not measured AlphaGenome predictions.
"""
from __future__ import annotations

# Allele Functional
# ======================================================================
import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
import gzip
import itertools
import json
import math
from pathlib import Path
import re

from analysis_integrity import identity, transactional

BASES = "ACGT"
SEMANTICS = {"signed_delta_to_ref", "allele_activity"}
MISSING = {"", ".", "NA"}
PROFILE_COLUMNS = {
    "chrom", "pos", "ref", "allele", "feature_id", "value", "semantics",
    "background_id", "model_version",
}
SCHEMA_COLUMNS = {"feature_id", "modality", "scale", "active_threshold", "provenance"}
COUNT_COLUMNS = {"chrom", "pos", "ref", "population", "AN", *(b + "_AC" for b in BASES)}
SUMMARY_PATTERN = re.compile(
    r"(?:^|[^A-Z0-9])(?:SHAP|FEATURE_IMPORTANCE|FEATURE_ATTRIBUTIONS?|"
    r"AVI(?:_[A-Z0-9]+)*|PHRED|MAX_ABS(?:_[A-Z0-9]+)*|MAXABS|"
    r"MERGED_SPLICING|ALPHAMISSENSE|CACTUS_241_WAY|PHASTCONS_470_WAY|"
    r"PROTEIN_TERMINATION|START_LOST|STOP_LOST|IS_INSERTION|IS_DELETION)"
    r"(?:$|[^A-Z0-9])"
)


@dataclass(frozen=True)
class Feature:
    feature_id: str
    modality: str
    scale: float
    active_threshold: float
    provenance: str


@dataclass
class Site:
    chrom: str
    pos: int
    ref: str
    semantics: str
    background_id: str
    model_version: str
    values: dict[str, dict[str, float]] = field(default_factory=lambda: {b: {} for b in BASES})
    n_input_rows: int = 0

    @property
    def coordinate(self):
        return self.chrom, self.pos

    @property
    def labels(self):
        return dict(chrom=self.chrom, pos=self.pos, ref=self.ref,
                    semantics=self.semantics, background_id=self.background_id,
                    model_version=self.model_version)


def text_value(value, name):
    value = str(value).strip()
    if value in MISSING:
        raise ValueError(f"{name} must be explicitly supplied")
    return value


def finite(value, name, *, positive=False):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric; omit missing profile rows instead of imputing") from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return result


def integer(value, name, *, missing_allowed=False):
    if missing_allowed and value in MISSING:
        return None
    try:
        result = Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise ValueError(f"{name} must be a nonnegative integer") from exc
    if not result.is_finite() or result < 0 or result != result.to_integral_value():
        raise ValueError(f"{name} must be a finite nonnegative integer")
    return int(result)


def coordinate(row):
    chrom = text_value(row["chrom"], "chrom")
    # Only conventional primary-chromosome aliases are normalised; no liftover,
    # complementing, REF/ALT swapping or population-major recoding is performed.
    if chrom in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        chrom = "chr" + chrom
    elif chrom in {"MT", "M", "chrMT"}:
        chrom = "chrM"
    if not re.fullmatch(r"[A-Za-z0-9_.:-]+", chrom):
        raise ValueError(f"Invalid chromosome label: {chrom!r}")
    pos = integer(row["pos"], "pos")
    if pos < 1:
        raise ValueError("pos must be 1-based and positive")
    ref = row["ref"]
    if ref not in BASES or len(ref) != 1:
        raise ValueError("ref must be an explicit uppercase A/C/G/T nucleotide")
    return chrom, pos, ref


def table_rows(path, required, build, metadata):
    """Read tab-separated data, preserving optional leading ##key=value metadata."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        header = None
        line_number = 0
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            if line.startswith("##"):
                item = line[2:].strip()
                if "=" in item:
                    key, value = item.split("=", 1)
                    if key in metadata:
                        raise ValueError(f"Duplicate metadata key {key}: {path}")
                    metadata[key] = value
                continue
            header = next(csv.reader([line], delimiter="\t"))
            break
        if header is None:
            raise ValueError(f"Empty input: {path}")
        if len(header) != len(set(header)):
            raise ValueError(f"Duplicate column names in {path}")
        if required - set(header):
            raise ValueError(f"{path}: missing columns {sorted(required-set(header))}")
        for label in ("build", "genome_build", "reference_build"):
            if label in metadata and metadata[label] != build:
                raise ValueError(f"{path}: {label} conflicts with declared build {build}")
        reader = csv.DictReader(handle, fieldnames=header, delimiter="\t")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed TSV row near line {line_number+reader.line_num}: {path}")
            if "build" in row and row["build"] != build:
                raise ValueError(f"Row build conflicts with {build}: {path}")
            yield row


def load_schema(path, build, metadata):
    features = {}
    for row in table_rows(path, SCHEMA_COLUMNS, build, metadata):
        fid = text_value(row["feature_id"], "feature_id")
        modality = text_value(row["modality"], "modality")
        provenance = text_value(row["provenance"], "provenance")
        if any(SUMMARY_PATTERN.search(x.upper()) for x in (fid, modality, provenance)):
            raise ValueError(f"Unsupported aggregate/attribution feature {fid}: use matched raw molecular readouts")
        if fid in features:
            raise ValueError(f"Duplicate feature_id in schema: {fid}")
        features[fid] = Feature(fid, modality, finite(row["scale"], "scale", positive=True),
                                finite(row["active_threshold"], "active_threshold", positive=True), provenance)
    if not features:
        raise ValueError("Feature schema is empty")
    return features


def iter_sites(path, features, build, metadata):
    current = None
    seen = set()
    for row in table_rows(path, PROFILE_COLUMNS, build, metadata):
        chrom, pos, ref = coordinate(row)
        key = chrom, pos
        semantics = row["semantics"]
        if semantics not in SEMANTICS:
            raise ValueError("semantics must be signed_delta_to_ref or allele_activity; SHAP is not a molecular effect")
        background = text_value(row["background_id"], "background_id")
        model = text_value(row["model_version"], "model_version")
        if current is None or key != current.coordinate:
            if current is not None:
                yield current
            if key in seen:
                raise ValueError(f"Site {key} reappears; all profile rows for a site must be contiguous")
            seen.add(key)
            current = Site(chrom, pos, ref, semantics, background, model)
        if (ref, semantics, background, model) != (current.ref, current.semantics, current.background_id, current.model_version):
            raise ValueError(f"{key}: REF, semantics, background_id and model_version must match within site")
        allele = row["allele"]
        if allele not in BASES or len(allele) != 1:
            raise ValueError("allele must be an uppercase A/C/G/T nucleotide, never ALT1/ALT2/rank labels")
        fid = row["feature_id"]
        if fid not in features:
            raise ValueError(f"Unknown feature_id (not in declared schema): {fid}")
        if fid in current.values[allele]:
            raise ValueError(f"Duplicate profile for {key}, {allele}, {fid}")
        value = finite(row["value"], "profile value")
        if semantics == "signed_delta_to_ref" and allele == ref and value != 0:
            raise ValueError("A signed_delta_to_ref REF row must equal zero; it is a baseline, not an activity")
        current.values[allele][fid] = value
        current.n_input_rows += 1
    if current is not None:
        yield current


def load_counts(path, build, metadata):
    counts = defaultdict(list)
    seen = set()
    refs = {}
    if not path:
        return counts
    for row in table_rows(path, COUNT_COLUMNS, build, metadata):
        chrom, pos, ref = coordinate(row)
        pop = text_value(row["population"], "population")
        key = chrom, pos
        unique = chrom, pos, pop
        if unique in seen:
            raise ValueError(f"Duplicate site/population count row: {unique}")
        if key in refs and refs[key] != ref:
            raise ValueError(f"Inconsistent REF across population counts: {key}")
        seen.add(unique)
        refs[key] = ref
        if row.get("count_semantics", "joint_nucleotide_counts") != "joint_nucleotide_counts":
            raise ValueError("counts must be joint_nucleotide_counts on one common sample denominator")
        values = {b: integer(row[b+"_AC"], b+"_AC", missing_allowed=True) for b in BASES}
        an = integer(row["AN"], "AN", missing_allowed=True)
        missing = [b+"_AC" for b in BASES if values[b] is None] + (["AN"] if an is None else [])
        present_sum = sum(x for x in values.values() if x is not None)
        if an is not None and (present_sum > an or (not missing and present_sum != an)):
            raise ValueError(f"{key}/{pop}: AN must equal A_AC+C_AC+G_AC+T_AC")
        counts[key].append(dict(chrom=chrom, pos=pos, ref=ref, population=pop,
                               AN=an, counts=values, missing_fields=missing))
    if metadata.get("count_semantics", "joint_nucleotide_counts") != "joint_nucleotide_counts":
        raise ValueError("Count metadata does not declare joint nucleotide counts")
    return counts


def scaled_profiles(site, features):
    result = {b: {fid: finite(value/features[fid].scale, "scaled profile")
                  for fid, value in site.values[b].items()} for b in BASES}
    if site.semantics == "signed_delta_to_ref":
        result[site.ref] = {fid: 0.0 for fid in features}
    return result


def distance(left, right, fids):
    fids = sorted(fids)
    if not fids:
        return dict(distance_l2_scaled=None, distance_rms_scaled=None, mean_abs_scaled_difference=None)
    diff = [finite(left[fid]-right[fid], "pairwise molecular difference") for fid in fids]
    l2 = finite(math.hypot(*diff), "pairwise distance")
    return dict(distance_l2_scaled=l2, distance_rms_scaled=l2/math.sqrt(len(fids)),
                mean_abs_scaled_difference=math.fsum(abs(x)/len(fids) for x in diff))


def allele_summaries(site, features, scaled):
    rows = []
    for base in BASES:
        values = scaled[base]
        delta_features = set(values) & set(scaled[site.ref])
        deltas = {fid: finite(values[fid]-scaled[site.ref][fid], "REF-relative delta") for fid in delta_features}
        active = [fid for fid, value in deltas.items() if abs(value) >= features[fid].active_threshold]
        modalities = {features[fid].modality for fid in active}
        max_abs = max((abs(x) for x in deltas.values()), default=None)
        # Participation ratio of squared feature effects. Normalise first so
        # extreme but finite measurements cannot overflow the energy sums.
        pr = None
        if max_abs is not None and max_abs > 0:
            energy = [(value/max_abs)**2 for value in deltas.values()]
            pr = math.fsum(energy)**2 / math.fsum(e*e for e in energy)
        is_baseline = base == site.ref and site.semantics == "signed_delta_to_ref"
        complete = len(deltas) == len(features)
        rows.append(dict(site.labels, allele=base,
                         profile_source="reference_zero_baseline_not_activity_prediction" if is_baseline else
                         ("provided_molecular_profile" if site.values[base] else "missing"),
                         n_provided_features=len(site.values[base]), n_profile_features=len(values),
                         n_expected_features=len(features), n_delta_features=len(deltas),
                         n_active_features_relative_to_ref=len(active) if deltas else None,
                         n_active_modalities_relative_to_ref=len(modalities) if deltas else None,
                         active_modalities=";".join(sorted(modalities)) if deltas else None,
                         max_abs_scaled_delta_relative_to_ref=max_abs,
                         effect_energy_participation_ratio=pr,
                         breadth_status="complete_declared_panel" if complete else
                         ("observed_subset_lower_bound" if deltas else "unavailable"),
                         breadth_interpretation="molecular_effect_breadth_relative_to_REF; not_phenotypic_pleiotropy_or_fitness"))
    return rows


def pair_summaries(site, features, scaled):
    rows = []
    for a, b in itertools.combinations(BASES, 2):
        common = set(scaled[a]) & set(scaled[b])
        rows.append(dict(site.labels, allele_a=a, allele_b=b, n_common_features=len(common),
                         n_expected_features=len(features),
                         coverage_fraction=len(common)/len(features),
                         distance_status="complete_declared_panel" if len(common) == len(features) else
                         ("observed_common_subset" if common else "unavailable"),
                         reference_zero_baseline_involved=int(site.semantics == "signed_delta_to_ref" and site.ref in (a, b)),
                         **distance(scaled[a], scaled[b], common)))
    return rows


def rank_summary(site, features, scaled, count, min_ac):
    if count["ref"] != site.ref:
        raise ValueError(f"{site.coordinate}: count REF does not match profile REF; no allele swapping is permitted")
    result = dict(site.labels, population=count["population"], AN=count["AN"], min_ac=min_ac,
                  missing_count_fields=";".join(count["missing_fields"]),
                  rank_status="incomplete_counts", tie_groups="", n_observed_bases=None,
                  n_bases_ac_ge_min=None, n_common_rank_features=0, n_expected_features=len(features),
                  d_rank3_rank1_rms_scaled=None, d_rank3_rank2_rms_scaled=None,
                  d_rank2_rank1_rms_scaled=None, difference_d31_minus_d21=None,
                  difference_d31_minus_d32=None, nearest_of_rank1_rank2=None,
                  rank_interpretation="current_population_count_order; not_mutation_age_or_ancestry")
    for i in range(1, 4):
        result.update({f"rank{i}_base": None, f"rank{i}_ac": None, f"rank{i}_af": None})
    if count["missing_fields"]:
        return result
    values = count["counts"]
    an = count["AN"]
    result.update(n_observed_bases=sum(x > 0 for x in values.values()),
                  n_bases_ac_ge_min=sum(x >= min_ac for x in values.values()))
    if an == 0:
        result["rank_status"] = "no_called_copies"
        return result
    by_count = defaultdict(list)
    for base, ac in values.items():
        by_count[ac].append(base)
    result["tie_groups"] = ";".join(f"{ac}:{','.join(sorted(bb))}" for ac, bb in sorted(by_count.items(), reverse=True) if len(bb) > 1)
    order = sorted(BASES, key=lambda b: (-values[b], b))
    if values[order[2]] < min_ac:
        result["rank_status"] = "fewer_than_three_count_supported_bases"
        return result
    # Every equality between neighbouring ranks affects rank 1, 2 or 3 in a
    # four-state locus. Lexical ordering is never used to interpret a tie.
    if any(values[order[i]] == values[order[i+1]] for i in range(3)):
        result["rank_status"] = "ambiguous_count_rank"
        return result
    for i, base in enumerate(order[:3], 1):
        result.update({f"rank{i}_base": base, f"rank{i}_ac": values[base], f"rank{i}_af": values[base]/an})
    a, b, c = order[:3]
    common = set(scaled[a]) & set(scaled[b]) & set(scaled[c])
    result["n_common_rank_features"] = len(common)
    if not common:
        result["rank_status"] = "no_common_molecular_features"
        return result
    # Both distances MUST use the same three-way feature intersection; pairwise
    # intersections can otherwise make a missing-data artefact look like novelty.
    d31 = distance(scaled[c], scaled[a], common)["distance_rms_scaled"]
    d32 = distance(scaled[c], scaled[b], common)["distance_rms_scaled"]
    d21 = distance(scaled[b], scaled[a], common)["distance_rms_scaled"]
    result.update(d_rank3_rank1_rms_scaled=d31, d_rank3_rank2_rms_scaled=d32,
                  d_rank2_rank1_rms_scaled=d21, difference_d31_minus_d21=d31-d21,
                  difference_d31_minus_d32=d31-d32)
    if len(common) != len(features):
        result["rank_status"] = "incomplete_declared_molecular_panel"
        return result
    result["rank_status"] = "complete_descriptive"
    result["nearest_of_rank1_rank2"] = "equidistant" if math.isclose(d31, d32, rel_tol=1e-12, abs_tol=1e-12) else ("rank1" if d31 < d32 else "rank2")
    return result


def rank4_summary(site, features, scaled, count, min_ac):
    """Fourth-base geometry on one common four-way molecular feature panel.

    The comparison with rank2 or rank3 uses the SAME feature denominator. A
    frequency tie is not broken alphabetically and an unobserved fourth base is
    not interpreted as an extinct or future mutation.
    """
    rank3 = rank_summary(site, features, scaled, count, min_ac)
    result = {k: rank3[k] for k in (
        *site.labels, "population", "AN", "min_ac", "missing_count_fields",
        "tie_groups", "n_observed_bases", "n_bases_ac_ge_min", "n_expected_features",
        "rank_interpretation",
    )}
    result.update(rank4_status=rank3["rank_status"], n_common_rank4_features=0,
                  nearest_of_rank1_rank2_rank3=None)
    for i in range(1, 5):
        result.update({f"rank{i}_base": None, f"rank{i}_ac": None, f"rank{i}_af": None})
    for name in ("d_rank4_rank1_rms_scaled", "d_rank4_rank2_rms_scaled",
                 "d_rank4_rank3_rms_scaled", "d_rank3_rank1_four_panel_rms_scaled",
                 "d_rank2_rank1_four_panel_rms_scaled", "difference_d41_minus_d21",
                 "difference_d41_minus_d31"):
        result[name] = None
    if count["missing_fields"] or rank3["rank_status"] in {
        "no_called_copies", "fewer_than_three_count_supported_bases", "ambiguous_count_rank"
    }:
        return result
    values, an = count["counts"], count["AN"]
    order = sorted(BASES, key=lambda b: (-values[b], b))
    if values[order[3]] < min_ac:
        result["rank4_status"] = "fewer_than_four_count_supported_bases"
        return result
    for i, base in enumerate(order, 1):
        result.update({f"rank{i}_base": base, f"rank{i}_ac": values[base], f"rank{i}_af": values[base]/an})
    common = set.intersection(*(set(scaled[b]) for b in order))
    result["n_common_rank4_features"] = len(common)
    if not common:
        result["rank4_status"] = "no_common_molecular_features"
        return result
    a, b, c, d = order
    dist = lambda u, v: distance(scaled[u], scaled[v], common)["distance_rms_scaled"]
    d41, d42, d43, d31, d21 = dist(d, a), dist(d, b), dist(d, c), dist(c, a), dist(b, a)
    result.update(d_rank4_rank1_rms_scaled=d41, d_rank4_rank2_rms_scaled=d42,
                  d_rank4_rank3_rms_scaled=d43, d_rank3_rank1_four_panel_rms_scaled=d31,
                  d_rank2_rank1_four_panel_rms_scaled=d21,
                  difference_d41_minus_d21=d41-d21, difference_d41_minus_d31=d41-d31)
    if len(common) != len(features):
        result["rank4_status"] = "incomplete_declared_molecular_panel"
        return result
    result["rank4_status"] = "complete_descriptive"
    smallest = min(d41, d42, d43)
    nearest = [f"rank{i}" for i, value in enumerate((d41, d42, d43), 1)
               if math.isclose(value, smallest, rel_tol=1e-12, abs_tol=1e-12)]
    result["nearest_of_rank1_rank2_rank3"] = nearest[0] if len(nearest) == 1 else "equidistant:" + ",".join(nearest)
    return result


ALLELE_FIELDS = ["chrom", "pos", "ref", "semantics", "background_id", "model_version", "allele", "profile_source",
                 "n_provided_features", "n_profile_features", "n_expected_features", "n_delta_features",
                 "n_active_features_relative_to_ref", "n_active_modalities_relative_to_ref", "active_modalities",
                 "max_abs_scaled_delta_relative_to_ref", "effect_energy_participation_ratio", "breadth_status", "breadth_interpretation"]
PAIR_FIELDS = ["chrom", "pos", "ref", "semantics", "background_id", "model_version", "allele_a", "allele_b",
               "n_common_features", "n_expected_features", "coverage_fraction", "distance_status", "reference_zero_baseline_involved",
               "distance_l2_scaled", "distance_rms_scaled", "mean_abs_scaled_difference"]
RANK_FIELDS = ["chrom", "pos", "ref", "semantics", "background_id", "model_version", "population", "AN", "min_ac",
               "missing_count_fields", "rank_status", "tie_groups", "n_observed_bases", "n_bases_ac_ge_min",
               *[f"rank{i}_{suffix}" for i in range(1, 4) for suffix in ("base", "ac", "af")],
               "n_common_rank_features", "n_expected_features", "d_rank3_rank1_rms_scaled", "d_rank3_rank2_rms_scaled",
               "d_rank2_rank1_rms_scaled", "difference_d31_minus_d21",
               "difference_d31_minus_d32", "nearest_of_rank1_rank2", "rank_interpretation"]
RANK4_FIELDS = ["chrom", "pos", "ref", "semantics", "background_id", "model_version", "population", "AN", "min_ac",
                "missing_count_fields", "rank4_status", "tie_groups", "n_observed_bases", "n_bases_ac_ge_min",
                *[f"rank{i}_{suffix}" for i in range(1, 5) for suffix in ("base", "ac", "af")],
                "n_common_rank4_features", "n_expected_features", "d_rank4_rank1_rms_scaled",
                "d_rank4_rank2_rms_scaled", "d_rank4_rank3_rms_scaled", "d_rank3_rank1_four_panel_rms_scaled",
                "d_rank2_rank1_four_panel_rms_scaled", "difference_d41_minus_d21", "difference_d41_minus_d31",
                "nearest_of_rank1_rank2_rank3", "rank_interpretation"]
SITE_FIELDS = ["chrom", "pos", "ref", "semantics", "background_id", "model_version", "n_input_rows",
               "provided_bases", "reference_baseline_inferred", "n_expected_features", "n_common_four_features",
               "four_profile_status", "n_count_populations", "n_unavailable_pairs"]
COUNT_AUDIT_FIELDS = ["chrom", "pos", "ref", "population", "status", "missing_fields"]


class TableWriter:
    def __init__(self, path, columns):
        self.handle = Path(path).open("w", newline="")
        self.writer = csv.DictWriter(self.handle, columns, delimiter="\t", lineterminator="\n")
        self.writer.writeheader()

    def write(self, row):
        self.writer.writerow({k: "NA" if v is None else v for k, v in row.items()})

    def close(self):
        self.handle.close()


@transactional("manifest.json", input_files=("profiles", "feature_schema", "counts"))
def run(args):
    if args.build != "GRCh38":
        raise ValueError("Only explicitly declared GRCh38 inputs are supported")
    if args.min_ac < 1:
        raise ValueError("min-ac must be positive; 5 is the default and 10 is the prespecified sensitivity option")
    inputs = {key: identity(getattr(args, key), full=True) for key in ("profiles", "feature_schema", "counts") if getattr(args, key, None)}
    metadata = {key: {} for key in ("profiles", "feature_schema", "counts")}
    features = load_schema(args.feature_schema, args.build, metadata["feature_schema"])
    counts = load_counts(args.counts, args.build, metadata["counts"])
    output = Path(args.out)
    output.mkdir(parents=True, exist_ok=True)
    writers = {
        "allele": TableWriter(output/"molecular_allele_profiles.tsv", ALLELE_FIELDS),
        "pair": TableWriter(output/"molecular_pairwise_distances.tsv", PAIR_FIELDS),
        "rank": TableWriter(output/"molecular_rank_geometry.tsv", RANK_FIELDS),
        "rank4": TableWriter(output/"molecular_rank4_geometry.tsv", RANK4_FIELDS),
        "site": TableWriter(output/"molecular_profile_audit.tsv", SITE_FIELDS),
        "count": TableWriter(output/"molecular_count_audit.tsv", COUNT_AUDIT_FIELDS),
    }
    audit = Counter()
    seen = set()
    model_versions = set()
    try:
        for site in iter_sites(args.profiles, features, args.build, metadata["profiles"]):
            seen.add(site.coordinate)
            model_versions.add(site.model_version)
            audit["sites"] += 1
            scaled = scaled_profiles(site, features)
            for row in allele_summaries(site, features, scaled):
                writers["allele"].write(row)
                audit["allele_breadth:" + row["breadth_status"]] += 1
            pairs = pair_summaries(site, features, scaled)
            for row in pairs:
                writers["pair"].write(row)
                audit["pair_distance:" + row["distance_status"]] += 1
            common_four = set.intersection(*(set(scaled[b]) for b in BASES))
            matched = counts.get(site.coordinate, [])
            for count in matched:
                row = rank_summary(site, features, scaled, count, args.min_ac)
                writers["rank"].write(row)
                audit["rank:" + row["rank_status"]] += 1
                row4 = rank4_summary(site, features, scaled, count, args.min_ac)
                writers["rank4"].write(row4)
                audit["rank4:" + row4["rank4_status"]] += 1
                writers["count"].write({k: count[k] for k in ("chrom", "pos", "ref", "population")} |
                                       dict(status=row["rank_status"], missing_fields=";".join(count["missing_fields"])))
            writers["site"].write(dict(site.labels, n_input_rows=site.n_input_rows,
                                       provided_bases=",".join(b for b in BASES if site.values[b]),
                                       reference_baseline_inferred=int(site.semantics == "signed_delta_to_ref"),
                                       n_expected_features=len(features), n_common_four_features=len(common_four),
                                       four_profile_status="complete_declared_panel" if len(common_four) == len(features) else "incomplete",
                                       n_count_populations=len(matched),
                                       n_unavailable_pairs=sum(r["distance_status"] == "unavailable" for r in pairs)))
        if not seen:
            raise ValueError("No molecular profile rows were provided")
        for coord, group in counts.items():
            if coord not in seen:
                for count in group:
                    audit["count_rows_without_profile"] += 1
                    writers["count"].write({k: count[k] for k in ("chrom", "pos", "ref", "population")} |
                                           dict(status="no_profile_for_count_site", missing_fields=";".join(count["missing_fields"])))
    finally:
        for writer in writers.values():
            writer.close()
    for key, before in inputs.items():
        if before != identity(getattr(args, key), full=True):
            raise ValueError(f"Input changed during molecular analysis: {key}")
    incomplete = (audit["pair_distance:observed_common_subset"] + audit["pair_distance:unavailable"] +
                  sum(n for k, n in audit.items() if k.startswith("rank:") and k != "rank:complete_descriptive") +
                  audit["count_rows_without_profile"])
    manifest = dict(
        status="complete_with_exclusions" if incomplete else "complete",
        analysis="matched_raw_molecular_profiles_v2", build=args.build,
        reference_validation="exact_input_REF_consistency; not_checked_against_reference_FASTA",
        input_checksums=inputs, input_metadata=metadata, settings=vars(args),
        declared_feature_panel=[vars(feature) for feature in features.values()],
        schema_contract="feature_id must encode exact gene/tissue/region/track identity; scales and thresholds declared independently of counts/outcomes",
        distance_definition="sqrt(mean(((value_a-value_b)/declared_scale)^2)); rank3 contrasts use one three-way common panel; rank4 contrasts use one four-way common panel",
        rank_hypothesis_contrasts="d(rank3,rank1)-d(rank2,rank1); d(rank4,rank1)-d(rank2,rank1); positive means further from current major on the same molecular panel, not newer or less fit",
        breadth_definition="number of observed molecular modalities with >=1 abs(REF-relative delta/scale) >= feature threshold",
        participation_ratio_definition="(sum squared scaled deltas)^2 / sum fourth powers; correlated features are not independent traits",
        ranking_definition="descending actual nucleotide counts within population; all 4 AC values and common AN required; any tie affecting top 3 excludes nearest interpretation",
        min_ac=args.min_ac, sensitivity_min_ac=10,
        counts_status="provided" if args.counts else "not_provided; theoretical_allele_geometry_only",
        count_provenance_limit="numeric coherence validated; common called sample and genotype provenance are asserted by input provider",
        ref_zero_meaning="signed-delta reference zero is a contrast baseline, not zero biological function or a prediction of reference activity",
        model_versions=sorted(model_versions), audit=dict(audit),
        interpretation="descriptive molecular effects only; no fitness, phenotypic pleiotropy, antagonistic pleiotropy, selection or mutation-age inference",
    )
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False, default=str)+"\n")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profiles", required=True, help="Raw profiles TSV/TSV.GZ; site-contiguous rows")
    parser.add_argument("--feature-schema", required=True, help="Predeclared exact-feature panel, positive scales and thresholds")
    parser.add_argument("--counts", help="Optional coherent population-specific A/C/G/T nucleotide counts")
    parser.add_argument("--build", required=True, choices=["GRCh38"], help="Explicit declaration applying to every input")
    parser.add_argument("--min-ac", type=int, default=5, help="Minimum count of current rank3; rerun at 10 as sensitivity")
    parser.add_argument("--out", required=True, help="Dedicated output directory, disjoint from inputs")
    result = run(parser.parse_args(argv))
    print(json.dumps({"status": result["status"], "audit": result["audit"]}, sort_keys=True))


# Allele Strata
# ======================================================================
_description_allele_strata = 'AVI/region stratification of observed three- and four-base SNV sites.\n\nThis streaming analysis joins actual joint A/C/G/T counts to the four-base\noutput of atlas_scores.py. The site predictor is the maximum (or mean) AVI\nPHRED over ALL three theoretical ALT substitutions, selected without looking\nat which ALT is observed. REF has no substitution score. AVI is frequency-\ntrained, so these associations are descriptive, not independent selection\ntests. Raw matched molecular profiles belong in allele_functional.py instead.\n\nCounts columns: chrom pos ref population AN A_AC C_AC G_AC T_AC\n               callable_all_four call_rate\nAll four counts must be jointly measured on one denominator, including measured\nzeros. Four-base callability is an explicit input assertion, not inferred from\nan absent VCF ALT. No counts are reconstructed from split AC/AN summary rows.\n\nOptional annotations: chrom pos ref, plus caller-selected strata and binary\nannotation columns. Missing annotations never become zero. All input tables\nmust be sorted by canonical chromosome order 1..22,X,Y,MT then position; rows\nfor populations/alleles at the same site must be contiguous. This permits\nwhole-chromosome batches without loading site-level tables into memory.\n\nCounts are summarized within populations. Comparative odds ratios are\nMantel-Haenszel estimates over chromosome, log2 AN bins and optional prespecified\nstrata (e.g. trinucleotide context and mappability). Uncertainty uses genomic-\nblock bootstrap, with blocks resampled within chromosome. No P values are\nreported. A sample selected for multi-allelic status is refused as an enrichment\ndenominator. An independent/prespecified target subset remains target-conditional.\n'

import html
import random


STRATA_BASES = "ACGT"
STRATA_MISSING = {"", ".", "NA", "NAN", "NONE", "NULL"}
STRATA_VERSION = "allele-strata-1.0"
STRATA_LIMITS = [
    "The denominator is the declared callable sampling frame after count/callability filters; it is not automatically the whole genome.",
    "Three/four-base support uses actual A/C/G/T counts >= min_ac within one population and common AN. REF is neither assumed major nor ancestral.",
    "An absent ALT row is not an observed zero. callable_all_four is an input-provider assertion that requires joint WGS/gVCF evidence.",
    "The site AVI statistic uses all three theoretical nonreference substitutions, independently of their observed counts; the maximum's bins are not the single-variant genome-wide percentiles.",
    "AVI was trained using gnomAD frequency-derived labels and includes conservation information. AVI/frequency associations, especially in gnomAD, are not independent tests of selection or conservation.",
    "Annotation overlaps are associations. HERV overlap does not make an SNV a viral insertion, and introgressed-tract overlap does not prove an allele's archaic origin.",
    "Block bootstrap accounts for dependence within supplied/fixed genomic blocks, not arbitrary long-range LD, family structure or cohort overlap. Repeat with larger blocks, independent cohorts and technical masks.",
    "AN bins, supplied context/mappability strata and population separation address only measured confounding. Coverage, mutation supply and ascertainment can still produce enrichment.",
    "Current count rank does not date mutation origin; no fitness, phenotypic pleiotropy or balancing-selection mechanism is inferred.",
]


def strata_number(value, label, *, integer=False, missing=False):
    if str(value).strip().upper() in STRATA_MISSING:
        if missing:
            return None
        raise ValueError(f"Missing {label}")
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {label}: {value!r}") from exc
    if not math.isfinite(x) or x < 0 or integer and (x != int(x) or x > 2**53):
        raise ValueError(f"Invalid {label}: {value!r}")
    return int(x) if integer else x


def strata_location(row):
    chrom = re.sub(r"^chr", "", row["chrom"], flags=re.I).upper()
    chrom = "MT" if chrom == "M" else chrom
    if chrom not in {*(str(i) for i in range(1, 23)), "X", "Y", "MT"}:
        raise ValueError(f"Unsupported chromosome {chrom}")
    order = int(chrom) if chrom.isdigit() else {"X": 23, "Y": 24, "MT": 25}[chrom]
    pos = strata_number(row["pos"], "POS", integer=True)
    if not pos or row["ref"] not in STRATA_BASES or len(row["ref"]) != 1:
        raise ValueError("Require positive 1-based POS and explicit uppercase SNV REF")
    if row.get("build", "GRCh38").lower() not in {"grch38", "hg38", "38"}:
        raise ValueError("Every input must be GRCh38; implicit liftover is forbidden")
    return (order, pos), chrom, row["ref"]


def strata_table(path, required):
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        first = None
        for raw in stream:
            if raw.startswith("##") or not raw.strip():
                continue
            first = raw
            break
        if first is None:
            raise ValueError(f"Empty TSV: {path}")
        reader = csv.DictReader(itertools.chain([first], stream), delimiter="\t")
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or set(required) - set(fields):
            raise ValueError(f"{path}: missing/duplicate fields; required {sorted(required)}")
        previous = None
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f"Malformed TSV row in {path}")
            row = {k: v.strip() for k, v in row.items()}
            key, chrom, ref = strata_location(row)
            if previous is not None and key < previous:
                raise ValueError(f"{path}: rows must be sorted in chromosome/position order")
            previous = key
            row.update(_key=key, _chrom=chrom, _ref=ref)
            yield row


def strata_grouped(path, required):
    for key, records in itertools.groupby(strata_table(path, required), key=lambda r: r["_key"]):
        rows = list(records)
        if len({r["_ref"] for r in rows}) != 1:
            raise ValueError(f"{path}: conflicting REF at {key}")
        yield key, rows


class strata_Join:
    def __init__(self, groups):
        self.groups = iter(groups)
        self.current = next(self.groups, None)

    def get(self, key, ref):
        while self.current is not None and self.current[0] < key:
            self.current = next(self.groups, None)
        if self.current is None or self.current[0] != key:
            return []
        rows = self.current[1]
        if rows[0]["_ref"] != ref:
            raise ValueError(f"Joined inputs disagree on REF at {key}; no allele swaps")
        return rows


def strata_avi_site(rows, statistic):
    if not rows:
        return None
    bases = {}
    sources = set()
    for row in rows:
        base, ref = row["allele"], row["_ref"]
        if len(base) != 1 or base not in STRATA_BASES or base in bases:
            raise ValueError("Atlas requires one unique row per A/C/G/T base")
        if row["allele_role"] != ("REF" if base == ref else "ALT"):
            raise ValueError("Atlas allele_role disagrees with REF")
        value = strata_number(row["avi_phred"], "AVI PHRED", missing=True)
        if base == ref and value is not None:
            raise ValueError("REF is not a substitution and cannot have AVI=0 or another score")
        if base != ref and value is not None and row["status"] not in {"available", "available_partial"}:
            raise ValueError("Atlas AVI availability contradicts status")
        if base != ref and value is None and row["status"] == "available":
            raise ValueError("Atlas marks a missing PHRED score available")
        bases[base] = value
        sources.add((row["resource_id"], row["source_version"]))
    if len(sources) != 1:
        raise ValueError("Mixed AVI resources or versions at one site")
    if set(bases) != set(STRATA_BASES):
        raise ValueError("Atlas derivative must list all four bases including explicit missing/REF rows")
    values = [v for b, v in bases.items() if b != rows[0]["_ref"]]
    if any(v is None for v in values):
        return None
    return (max(values) if statistic == "max" else math.fsum(values)/3), next(iter(sources))


def strata_validate_count(row, args):
    if row.get("is_preview", "0") not in {"0", "1"}:
        raise ValueError("is_preview must be 0 or 1")
    if row.get("is_preview", "0") == "1":
        raise ValueError("A first-N preview cannot provide the prespecified enrichment sampling frame")
    if row.get("count_semantics", "joint_nucleotide_counts") != "joint_nucleotide_counts":
        raise ValueError("Require joint_nucleotide_counts; split population AC/AN is not accepted")
    raw = [strata_number(row[b+"_AC"], b+"_AC", integer=True, missing=True) for b in STRATA_BASES]
    an = strata_number(row["AN"], "AN", integer=True, missing=True)
    if an is not None and sum(v for v in raw if v is not None) > an:
        raise ValueError("Sum of nucleotide counts exceeds AN")
    if an is None or any(v is None for v in raw):
        return None, "missing_joint_counts"
    if sum(raw) != an:
        raise ValueError("A_AC+C_AC+G_AC+T_AC must equal AN")
    if row["callable_all_four"] not in {"0", "1", "NA", ".", ""}:
        raise ValueError("callable_all_four must be 0, 1 or missing")
    if row["callable_all_four"] != "1":
        return None, "four_base_callability_not_asserted"
    call_rate = strata_number(row["call_rate"], "call_rate", missing=True)
    if call_rate is not None and call_rate > 1:
        raise ValueError("call_rate must be in [0,1]")
    if call_rate is None or call_rate < args.min_call_rate:
        return None, "missing_or_low_call_rate"
    if an < args.min_an:
        return None, "low_AN"
    n = sum(v >= args.min_ac for v in raw)
    if args.denominator == "polymorphic" and n < 2:
        return None, "fewer_than_two_supported_bases"
    return (an, n), "eligible"


def strata_labels(cutpoints):
    edges = [0.0, *cutpoints]
    return [f"[{edges[i]:g},{edges[i+1]:g})" for i in range(len(edges)-1)] + [f"[{edges[-1]:g},inf)"]


def strata_bin_avi(value, cutpoints):
    return sum(value >= cut for cut in cutpoints)


def strata_mh_estimate(cells, outcome, multiplicities=None):
    """Cells: (stratum, block) -> exposed/control arrays [n,>=3,4]."""
    merged = defaultdict(lambda: [[0, 0, 0], [0, 0, 0]])
    for (stratum, block), pair in cells.items():
        weight = 1 if multiplicities is None else multiplicities.get(block, 0)
        if not weight:
            continue
        for e in (0, 1):
            for i in range(3):
                merged[stratum][e][i] += pair[e][i]*weight
    numerator, denominator, informative = 0.0, 0.0, 0
    for pair in merged.values():
        if not pair[0][0] or not pair[1][0]:
            continue
        a, b = pair[1][outcome], pair[1][0]-pair[1][outcome]
        c, d = pair[0][outcome], pair[0][0]-pair[0][outcome]
        n = a+b+c+d
        if not a+c or not b+d:
            continue
        informative += 1
        numerator += a*d/n
        denominator += b*c/n
    if not informative or numerator == 0 or denominator == 0:
        return None, informative
    return numerator/denominator, informative


def strata_quantile(values, p):
    values = sorted(values)
    index = (len(values)-1)*p
    lo, hi = math.floor(index), math.ceil(index)
    return values[lo] + (values[hi]-values[lo])*(index-lo)


def strata_contrast(cells, outcome, reps, seed, min_blocks):
    estimate, informative = strata_mh_estimate(cells, outcome)
    blocks = sorted({block for _, block in cells})
    result = dict(mh_odds_ratio=estimate, informative_strata=informative, genomic_blocks=len(blocks),
                  ci95_low=None, ci95_high=None, bootstrap_replicates=reps,
                  bootstrap_valid=0, uncertainty_status="not_estimable")
    if estimate is None:
        return result
    if len(blocks) < min_blocks:
        result["uncertainty_status"] = "too_few_blocks_descriptive_estimate_only"
        return result
    if not reps:
        result["uncertainty_status"] = "bootstrap_not_requested"
        return result
    rng = random.Random(seed)
    by_chrom = defaultdict(list)
    for block in blocks:
        by_chrom[block[0]].append(block)
    estimates = []
    for _ in range(reps):
        weights = Counter()
        for group in by_chrom.values():
            weights.update(rng.choices(group, k=len(group)))
        value, _ = strata_mh_estimate(cells, outcome, weights)
        if value is not None:
            estimates.append(value)
    result["bootstrap_valid"] = len(estimates)
    if len(estimates) >= max(20, math.ceil(0.8*reps)):
        result.update(ci95_low=strata_quantile(estimates, .025), ci95_high=strata_quantile(estimates, .975),
                      uncertainty_status="block_bootstrap_descriptive")
    else:
        result["uncertainty_status"] = "too_many_degenerate_bootstrap_replicates"
    return result


def strata_write_table(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if row.get(k) is None else row.get(k) for k in fields})


@transactional("qc/manifest.json", input_files=("counts", "atlas_variants", "annotations"))
def strata_run(args):
    if args.sampling_frame == "multi_allelic_selected":
        raise ValueError("A multi-allelic-selected target set cannot provide the enrichment denominator")
    if args.sampling_frame == "observed_variant_sites" and args.denominator != "polymorphic":
        raise ValueError("A variant-only sampling frame supports only the polymorphic denominator, not all callable bases")
    if args.min_ac < 1 or args.min_an < 1 or not 0 <= args.min_call_rate <= 1:
        raise ValueError("Require min-ac>=1, min-an>=1 and call rate in [0,1]")
    if args.block_bp < 1 or args.bootstrap < 0 or args.min_blocks < 20:
        raise ValueError("Require positive block-bp, nonnegative bootstrap and min-blocks>=20")
    cutpoints = [float(x) for x in args.avi_cutpoints.split(",")]
    if not cutpoints or any(not math.isfinite(x) or x <= 0 for x in cutpoints) or sorted(set(cutpoints)) != cutpoints:
        raise ValueError("AVI cutpoints must be distinct positive increasing finite numbers")
    columns = lambda text: [x.strip() for x in text.split(",") if x.strip()]
    annotations, strata_columns = columns(args.annotation_columns), columns(args.strata_columns)
    if len(set(annotations)) != len(annotations) or len(set(strata_columns)) != len(strata_columns):
        raise ValueError("Annotation/stratum columns must not repeat")
    if set(annotations) & set(strata_columns):
        raise ValueError("An exposure cannot also be an exact-matching stratum")
    if (annotations or strata_columns or args.block_column) and not args.annotations:
        raise ValueError("Selected annotation/stratum/block columns require --annotations")
    if args.annotations_only and args.atlas_variants:
        raise ValueError("--annotations-only cannot be combined with --atlas-variants")
    if not args.atlas_variants and not annotations:
        raise ValueError("Supply --atlas-variants or --annotations with at least one --annotation-columns exposure")
    inputs = {key: identity(getattr(args, key)) for key in ("counts", "atlas_variants", "annotations") if getattr(args, key)}
    atlas = strata_Join(strata_grouped(args.atlas_variants, {"chrom", "pos", "ref", "allele", "allele_role", "avi_phred",
                                            "status", "resource_id", "source_version"})) if args.atlas_variants else None
    annotation_join = strata_Join(strata_grouped(args.annotations, {"chrom", "pos", "ref", *annotations, *strata_columns,
                                                      *([args.block_column] if args.block_column else [])})) if args.annotations else None
    required = {"chrom", "pos", "ref", "population", "AN", *(b+"_AC" for b in STRATA_BASES), "callable_all_four", "call_rate"}
    # Aggregate by genomic block and stratum, retaining no individual site rows.
    aggregates = defaultdict(lambda: [0, 0, 0])
    audit, sources, cohorts = Counter(), set(), set()
    for loc, rows in strata_grouped(args.counts, required):
        audit["count_sites"] += 1
        pops = [r["population"] for r in rows]
        if len(set(pops)) != len(pops) or any(p.upper() in STRATA_MISSING for p in pops):
            raise ValueError("Exactly one nonmissing count row per site/population is required")
        ref, chrom = rows[0]["_ref"], rows[0]["_chrom"]
        avi = strata_avi_site(atlas.get(loc, ref), args.avi_statistic) if atlas else None
        if avi is not None:
            sources.add(avi[1])
            if len(sources) > 1:
                raise ValueError("Use one comparable AVI resource/version for a strata analysis")
        ann_rows = annotation_join.get(loc, ref) if annotation_join else []
        if len(ann_rows) > 1:
            raise ValueError("Annotations must have one row per site")
        ann = ann_rows[0] if ann_rows else {}
        if args.block_column:
            block_id = ann.get(args.block_column, "NA")
            block_ok = block_id.upper() not in STRATA_MISSING
        else:
            block_id, block_ok = str((loc[1]-1)//args.block_bp), True
        extras = tuple(ann.get(k, "NA") for k in strata_columns)
        strata_ok = all(x.upper() not in STRATA_MISSING for x in extras)
        for row in rows:
            audit["count_site_population_rows"] += 1
            if row.get("cohort", "").upper() not in STRATA_MISSING:
                cohorts.add(row["cohort"])
                if len(cohorts) > 1:
                    raise ValueError("Analyze each cohort separately; mixed cohort counts are not pooled")
            count, status = strata_validate_count(row, args)
            audit["count_filter:"+status] += 1
            if count is None:
                continue
            an, nbase = count
            population = row["population"]
            # Coarse AN strata remain explicit in the result contract.
            stratum = (chrom, math.floor(math.log2(an)), *extras)
            block = chrom, block_id
            states = [1, int(nbase >= 3), int(nbase == 4)]
            categories = []
            if avi is not None:
                categories.append(("site_avi_"+args.avi_statistic, str(strata_bin_avi(avi[0], cutpoints))))
            elif atlas:
                audit["eligible_without_complete_three_ALT_AVI"] += 1
            for column in annotations:
                value = ann.get(column, "NA")
                if value.upper() in STRATA_MISSING:
                    audit["missing_annotation:"+column] += 1
                    continue
                if value not in {"0", "1"}:
                    raise ValueError(f"Annotation {column} must be 0/1/NA")
                categories.append((column, value))
            for dimension, category in categories:
                # Summary remains available without optional stratification fields.
                key = population, dimension, category, stratum if strata_ok else None, block if block_ok else None
                target = aggregates[key]
                for i in range(3):
                    target[i] += states[i]
            if not strata_ok:
                audit["eligible_without_complete_matching_strata"] += 1
            if not block_ok:
                audit["eligible_without_block_id"] += 1
    if not audit["count_sites"]:
        raise ValueError("No count rows")
    summary_counts = defaultdict(lambda: [0, 0, 0, 0])
    dimensions = defaultdict(set)
    by_dimension = defaultdict(list)
    for (pop, dim, cat, stratum, block), values in aggregates.items():
        sums = summary_counts[pop, dim, cat]
        for i in range(3):
            sums[i] += values[i]
        if stratum is not None and block is not None:
            sums[3] += values[0]
        dimensions[pop, dim].add(cat)
        by_dimension[pop, dim].append((cat, stratum, block, values))
    summary, contrasts = [], []
    bin_labels = strata_labels(cutpoints)
    for (pop, dim, cat), (total, multi, four, matched) in sorted(summary_counts.items()):
        summary.append(dict(population=pop, dimension=dim, category=bin_labels[int(cat)] if dim.startswith("site_avi_") else cat,
                            eligible_sites=total, three_or_four_base_sites=multi, four_base_sites=four,
                            fraction_three_or_four=multi/total, fraction_four=four/total,
                            sites_with_matching_strata_and_blocks=matched, sampling_frame=args.sampling_frame,
                            denominator=args.denominator, min_ac=args.min_ac))
    for (pop, dim), categories in sorted(dimensions.items()):
        for cat in sorted(categories-{"0"}, key=int):
            cells = defaultdict(lambda: [[0, 0, 0], [0, 0, 0]])
            for c, stratum, block, values in by_dimension[pop, dim]:
                if c in {cat, "0"} and stratum is not None and block is not None:
                    target = cells[stratum, block][int(c == cat)]
                    for i in range(3):
                        target[i] += values[i]
            for outcome, label in ((1, "at_least_three"), (2, "four")):
                stats = strata_contrast(cells, outcome, args.bootstrap, args.seed, args.min_blocks)
                contrasts.append(dict(population=pop, dimension=dim,
                    category=bin_labels[int(cat)] if dim.startswith("site_avi_") else cat,
                    reference=bin_labels[0] if dim.startswith("site_avi_") else "0", outcome=label,
                    **stats, sampling_frame=args.sampling_frame, inference="descriptive_association_not_selection"))
    output = Path(args.out)
    for directory in ("results", "reports", "qc"):
        (output/directory).mkdir(parents=True, exist_ok=True)
    strata_write_table(output/"results/strata_summary.tsv", summary, ["population", "dimension", "category", "eligible_sites",
                "three_or_four_base_sites", "four_base_sites", "fraction_three_or_four", "fraction_four",
                "sites_with_matching_strata_and_blocks", "sampling_frame", "denominator", "min_ac"])
    strata_write_table(output/"results/strata_contrasts.tsv", contrasts, ["population", "dimension", "category", "reference", "outcome",
                "mh_odds_ratio", "ci95_low", "ci95_high", "informative_strata", "genomic_blocks", "bootstrap_replicates",
                "bootstrap_valid", "uncertainty_status", "sampling_frame", "inference"])
    for key, before in inputs.items():
        if before != identity(getattr(args, key)):
            raise ValueError(f"Input changed during analysis: {key}")
    rows_html = "".join("<tr>"+"".join(f"<td>{html.escape(str(r[k]))}</td>" for k in
        ("population", "dimension", "category", "eligible_sites", "three_or_four_base_sites", "four_base_sites"))+"</tr>" for r in summary)
    notes = "".join(f"<li>{html.escape(x)}</li>" for x in STRATA_LIMITS)
    (output/"reports/strata.html").write_text("<!doctype html><meta charset='utf-8'><title>Allele3 strata</title>"
        "<style>body{font:16px system-ui;max-width:1150px;margin:32px auto;line-height:1.5}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:left}table{border-collapse:collapse}li{margin:8px 0}</style>"
        "<h1>Observed multi-allelic sites across AVI and annotation strata</h1>"
        f"<p>Sampling frame: <b>{html.escape(args.sampling_frame)}</b>; denominator: {html.escape(args.denominator)}; min AC: {args.min_ac}.</p>"
        "<p><a href='../results/strata_summary.tsv'>Summary TSV</a> · <a href='../results/strata_contrasts.tsv'>Stratified contrasts TSV</a></p>"
        "<table><tr><th>Population</th><th>Dimension</th><th>Category</th><th>Eligible</th><th>≥3 bases</th><th>4 bases</th></tr>"
        +rows_html+"</table><h2>Interpretation</h2><ul>"+notes+"</ul>", encoding="utf-8")
    result = dict(status="complete" if summary else "complete_with_exclusions", version=STRATA_VERSION, build=args.build,
        input_identity=inputs, settings=vars(args), audit=dict(audit), avi_sources=sorted(sources),
        avi_analysis="requested" if atlas else "not_requested_annotations_only", declared_cohorts=sorted(cohorts),
        matching_strata=["chromosome", "floor(log2(AN))", *strata_columns],
        block_policy=("annotation_column:"+args.block_column) if args.block_column else f"fixed_nonoverlapping_{args.block_bp}_bp",
        sampling_declaration="The caller must establish the sampling frame and joint four-base callability from original data; this table cannot prove them.",
        missingness="Missing score/count/annotation/stratum is retained as unavailable, never imputed as zero.", limitations=STRATA_LIMITS)
    (output/"qc/manifest.json").write_text(json.dumps(result, indent=2, allow_nan=False, default=str)+"\n")
    return result


def strata_parser():
    p = argparse.ArgumentParser(description=_description_allele_strata, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True)
    p.add_argument("--atlas-variants", help="Optional sorted atlas_scores.py atlas_variants.tsv.gz; omit for annotation-only analysis")
    p.add_argument("--annotations-only", action="store_true", help="Explicit region-only analysis without AlphaGenome inputs")
    p.add_argument("--annotations", help="Sorted per-site strata/binary annotation TSV[.gz]")
    p.add_argument("--annotation-columns", default="")
    p.add_argument("--strata-columns", default="", help="Predeclared context/mappability etc. Exact matching; do not include outcome/rank")
    p.add_argument("--block-column", help="Optional annotation column defining genomic LD blocks")
    p.add_argument("--build", required=True, choices=["GRCh38"])
    p.add_argument("--sampling-frame", required=True, choices=["all_eligible_sites", "observed_variant_sites", "prespecified_targets", "multi_allelic_selected"])
    p.add_argument("--denominator", choices=["polymorphic", "callable"], default="polymorphic")
    p.add_argument("--min-ac", type=int, default=5)
    p.add_argument("--min-an", type=int, default=100)
    p.add_argument("--min-call-rate", type=float, default=.95)
    p.add_argument("--avi-statistic", choices=["max", "mean"], default="max")
    p.add_argument("--avi-cutpoints", default="10,20,30")
    p.add_argument("--block-bp", type=int, default=1000000)
    p.add_argument("--bootstrap", type=int, default=500)
    p.add_argument("--min-blocks", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261110)
    p.add_argument("--out", required=True)
    return p


def strata_main(argv=None):
    result = strata_run(strata_parser().parse_args(argv))
    print(json.dumps({"status": result["status"], "audit": result["audit"]}, sort_keys=True))


def cli(argv=None):
    """Dispatch auxiliary commands owned by this module."""
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    routes = {'strata': strata_main}
    entry = routes[argv.pop(0)] if argv and argv[0] in routes else main
    return entry(argv)



if __name__ == "__main__":
    cli()

