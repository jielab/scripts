#!/usr/bin/env python3
"""Query local AlphaGenome Atlas Tabix scores for bounded GRCh38 SNV targets.

The official AVI table contains one reference-relative score per substitution.
The merged per-chromosome VCF retains these values in Number=A INFO arrays.
Every target emits the structural REF (score NA) and all three theoretical ALT
bases. Listing an ALT in --targets is not a new population observation. This
module neither discovers alleles nor estimates a distance between ALT alleles.
It performs local CPU-only indexed reads; it never downloads or scans the atlas.
"""
from __future__ import annotations

# Atlas Resource Paths
# ======================================================================
_description_atlas_resource_paths = 'Metadata-only expansion of Atlas source paths for frozen batch plans.\n\nNo variant, index, genotype or matrix is opened here. Before genome targets\nexist, all declared canonical chromosome paths (including missing ones) are\nbound so that subsequent conversion/replacement cannot go unnoticed.\n'

import csv
import json
import os
from pathlib import Path

CHROMS = [str(n) for n in range(1, 23)] + ["X", "Y", "M"]


def require_atlas_inputs(dir_data):
    """Accept imported indexed resources without inventing conversion provenance.

    A present conversion receipt must still pass its original checks. Imported
    resources have no such receipt: require every autosome and its index, and
    leave header, build, Number=A and exact nucleotide checks to atlas_scores
    in the worker. File presence is never labelled a completed conversion.
    """
    root = Path(dir_data) / "annot/allele3/alphagenome/38"
    summary = root / "split2vcf.summary.json"
    if summary.is_file():
        require_complete_conversion(dir_data)
        return summary, "completed_conversion_receipt"
    for c in range(1, 23):
        vcf = root / f"chr{c}.vcf.bgz"
        indexes = [Path(str(vcf) + s) for s in (".tbi", ".csi")]
        if not vcf.is_file() or not vcf.stat().st_size or not any(p.is_file() and p.stat().st_size for p in indexes):
            raise ValueError("Missing/empty imported Atlas VCF or index: " + str(vcf))
    # Bind the absent receipt too: introducing it later invalidates stale plans.
    return summary, "imported_indexed_inputs_conversion_provenance_unverified"


def require_complete_conversion(dir_data):
    """Check only completion metadata and file sizes, without reading variants."""
    root = Path(dir_data) / "annot/allele3/alphagenome/38"
    summary = root / "split2vcf.summary.json"
    message = "AlphaGenome conversion is not complete; wait for alphagenome.split2vcf, then rerun the same genome command. "
    if not summary.is_file():
        raise ValueError(message + "Missing " + str(summary))
    result = json.loads(summary.read_text())
    autosomes = {f"chr{i}" for i in range(1, 23)}
    if result.get("status") != "complete" or not autosomes <= set(result.get("chromosomes", [])):
        raise ValueError(message + "Incomplete autosome inventory in " + str(summary))
    for chrom in sorted(autosomes):
        vcf = root / (chrom + ".vcf.bgz")
        indexes = [Path(str(vcf) + suffix) for suffix in (".tbi", ".csi")]
        if not vcf.is_file() or vcf.stat().st_size == 0 or not any(p.is_file() and p.stat().st_size for p in indexes):
            raise ValueError(message + "Missing/empty VCF or index: " + str(vcf))
    return summary


def manifest_source_paths(descriptor, dir_data):
    descriptor = Path(descriptor).expanduser().resolve()
    seen = set()
    with descriptor.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("enabled", "0") not in {"1", "TRUE", "true"}:
                continue
            for key in ("path", "index_path"):
                template = row.get(key, "").replace("{dir_data}", str(dir_data))
                if not template:
                    continue
                values = [template.replace("{chrom}", c) for c in CHROMS] if "{chrom}" in template else [template]
                for value in values:
                    path = Path(os.path.expandvars(value)).expanduser()
                    if not path.is_absolute():
                        path = descriptor.parent / path
                    paths = [path]
                    if key == "path":
                        paths += [Path(str(path) + suffix) for suffix in (".tbi", ".csi")]
                    for candidate in paths:
                        candidate = candidate.resolve()
                        if candidate not in seen:
                            seen.add(candidate)
                            yield candidate


# Atlas Scores
# ======================================================================
# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path

import argparse
from collections import Counter
from contextlib import ExitStack
import gzip
import hashlib
import math
import re
import struct
import sys

from analysis_integrity import InputPath, check_output_paths, identity, transactional, validate_manifest

VERSION = "atlas-scores-1.0"
SEMANTICS_VERSION = "atlas-metric-semantics-20261106"
STORAGE_VERSION = "atlas-indexed-tsv-or-merged-vcf-20261107"
BASES = "ACGT"
MISSING = {"", ".", "NA", "NAN", "NONE", "NULL"}
PORTAL = "https://alphagenome.google/downloads"
DEFAULT_RESOURCES = Path(__file__).resolve().parent.parent / "config" / "alphagenome_resources.tsv"
VCF_RESOURCES = DEFAULT_RESOURCES.with_name("alphagenome_resources.vcf.tsv")
AVI_COLUMNS = ["CHROM", "POS", "REF", "ALT", "raw_score", "PHRED"]
SPLICING_COLUMNS = ["CHROM", "POS", "REF", "ALT", "alphagenome_splicing"]
ATTRIBUTION_METRICS = [
    "MERGED_SPLICING", "MAX_ABS_ATAC", "MAX_ABS_CONTACT_MAPS", "MAX_ABS_DNASE",
    "MAX_ABS_CHIP_TF", "MAX_ABS_CHIP_HISTONE", "MAX_ABS_CAGE", "MAX_ABS_PROCAP",
    "MAX_ABS_RNA_SEQ", "MAX_ABS_POLYADENYLATION", "ALPHAMISSENSE", "CACTUS_241_WAY",
    "PROTEIN_TERMINATION", "START_LOST", "STOP_LOST", "PHASTCONS_470_WAY",
    "IS_INSERTION", "IS_DELETION",
]
ATTRIBUTION_COLUMNS = ["CHROM", "POS", "REF", "ALT", *ATTRIBUTION_METRICS]
VCF_METRICS = {
    "avi": {"raw_score": "AG_AVI_RAW", "PHRED": "AG_AVI_PHRED"},
    "splicing": {"alphagenome_splicing": "AG_SPLICE"},
    "features": {name: "AG_SHAP_" + name for name in ATTRIBUTION_METRICS},
}
VCF_SOURCE_BITS = {"avi": 1, "splicing": 2, "features": 4}
VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
AVI_TRAINING_CAVEAT = (
    "AVI was trained on gnomAD v4.1 genome FAF95_GRPMAX frequency-derived proxy labels. "
    "An AVI-versus-gnomAD frequency/retention association is not independent evidence of selection. "
    "External-cohort analysis still needs cohort-overlap and ascertainment checks."
)
SHAP_LIMIT = (
    "Signed contribution to the AVI prediction, not a raw molecular effect. The sign does not "
    "encode expression up/down; nonzero columns do not count biological functions, and differences "
    "between ALT SHAP vectors do not establish ALT-to-ALT functional distance."
)
VARIANT_FIELDS = ["site_id", "chrom", "pos", "build", "ref", "target_alts", "allele", "allele_role",
                  "observed_from_targets", "avi_raw", "avi_phred", "status", "resource_id", "source_version"]
METRIC_FIELDS = ["site_id", "chrom", "pos", "build", "ref", "allele", "allele_role", "observed_from_targets",
                 "resource_id", "source_version", "metric", "value", "status", "metric_semantics",
                 "training_label_caveat", "interpretation_limit"]
AVAILABILITY_FIELDS = ["site_id", "chrom", "pos", "ref", "target_alts", "avi_available_substitutions",
                       "avi_complete_substitutions", "avi_unavailable_substitutions", "reference_check"]
RESOURCE_FIELDS = ["resource_id", "kind", "status", "required", "enabled", "path", "index_path", "build",
                   "coordinate_system", "columns", "metric_columns", "source_version", "portal_url",
                   "download_url", "url_kind", "index_check", "header_check", "compression_check",
                   "metric_semantics", "training_label_caveat", "interpretation_limit", "storage_format",
                   "info_fields", "source_mask_bit", "physical_files"]
LIMITS = [
    "AVI raw and PHRED are reference-relative scalar predictions, not allele-to-allele functional distances, fitness, pleiotropy counts or tests of balancing selection. No pairwise score differences are reported.",
    "All three possible nonreference nucleotides are queried. observed_from_targets=1 means the nucleotide was listed in the supplied targets; it does not establish AC, AN, AF or a population observation.",
    "The structural REF is a sequence label and has no REF-to-REF substitution score. Its score and observed_from_targets are NA. REF is not assumed ancestral or the most frequent allele.",
    "Missing rows, missing metrics and absent indexed contigs remain NA. A numeric zero is retained as zero. Neither missing nor negative raw scores imply absence of the allele.",
    "GRCh38 declarations, recognized source headers, Tabix indexing convention and exact REF at returned target positions are checked. This is not a whole-genome FASTA audit; a target with no returned row has no sequence-level REF verification.",
    "User-supplied source versions and data/index pairing are recorded as declarations. Local file identity uses size and mtime for large atlas files; it is not an official release checksum, an exhaustive stale-index test or proof that an atlas was downloaded intact.",
    "Optional generic metrics require an explicit local SNV schema and retain their original units and semantics. They are not combined into a pleiotropy or functional-distance score.",
    "The official feature_importance SNV table stores signed SHAP contributions to AVI, including columns named MAX_ABS_*. These are not the raw absolute molecular predictions supplied to the AVI model. " + SHAP_LIMIT,
    AVI_TRAINING_CAVEAT,
    "Splicing scores summarize reference-relative predicted disruption. The sparse splicing table is not a whole-genome absence mask: missing rows remain NA rather than zero.",
    "Merged VCF ALT lists enumerate substitutions present in the prediction sources, not observed population alleles. Number=A scores are matched by the actual ALT base and AG_SOURCE_MASK retains per-source row absence. The same physical VCF batch is fetched and parsed once for its AVI, splicing and attribution resources.",
    "Requested positions are queried through bounded Tabix windows; intervening unrequested rows are discarded. No full atlas scan, API calls, network downloads, GPU inference, liftover, strand complementation, mutation-age inference or third-allele discovery are performed.",
]


def resource_semantics(row):
    """Assign explicit scientific meaning; known AVI attribution data cannot be relabeled.

    Atlas Methods: https://doi.org/10.64898/2026.09.16.26363192 . The exported
    feature_importance columns name AVI inputs, but their values are SHAP
    contributions to AVI. Do not infer their semantics from the numeric sign.
    """
    columns = [v.strip() for v in row["columns"].split(",")]
    basename = Path(row.get("path", "")).name.lower()
    if row["kind"] == "avi":
        semantics = "gnomad_frequency_trained_variant_impact"
        caveat = AVI_TRAINING_CAVEAT
        limit = "Reference-relative impact ranking; neither a fitness estimate nor a molecular direction, function count or ALT-to-ALT distance."
    elif columns == ATTRIBUTION_COLUMNS or "feature_importance" in basename:
        semantics = "signed_shap_contribution_to_avi"
        caveat, limit = AVI_TRAINING_CAVEAT, SHAP_LIMIT
    elif columns == SPLICING_COLUMNS or basename == "combined_alphagenome_splicing_snvs.tsv.gz":
        semantics = "reference_relative_splicing_prediction"
        caveat = "Not the frequency-trained AVI score. Splicing model predictions and their training provenance remain separate from population observations."
        limit = "Nonnegative splicing disruption summary; not effect direction, tissue breadth, a function count or an ALT-to-ALT distance. Missing splicing rows are NA."
    else:
        semantics = "user_declared_unverified_metric"
        caveat = "Training-label provenance not established for this user-supplied metric."
        limit = "Generic values retain user-supplied meaning and units; no biological function count or ALT-to-ALT distance is inferred."
    declared = row.get("metric_semantics", "").strip()
    if declared and declared != semantics:
        raise ValueError(f"Resource {row['resource_id']}: metric_semantics={declared!r} conflicts with the recognized source semantics {semantics!r}")
    return dict(metric_semantics=semantics, training_label_caveat=caveat, interpretation_limit=limit)


def metric_semantics(spec, metric):
    if spec["kind"] == "avi":
        return "gnomad_frequency_trained_avi_phred_rank" if metric == "PHRED" else "gnomad_frequency_trained_avi_raw_score"
    return spec["metric_semantics"]


def genome_build(value):
    normalized = str(value).strip().lower()
    aliases = {"38": "GRCh38", "grch38": "GRCh38", "hg38": "GRCh38",
               "37": "GRCh37", "grch37": "GRCh37", "hg19": "GRCh37"}
    if normalized not in aliases:
        raise ValueError(f"Unrecognized genome build {value!r}; declare GRCh38 explicitly")
    return aliases[normalized]


def chromosome(value):
    value = str(value).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    if value in {"M", "MT"}:
        return "MT"
    if value in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        return value
    raise ValueError(f"Unsupported chromosome {value!r}; expected 1-22, X, Y or MT")


def chromosome_key(value):
    return int(value) if value.isdigit() else {"X": 23, "Y": 24, "MT": 25}[value]


def chromosome_filter(value):
    if not value or value.lower() == "all":
        return None
    selected = set()
    for part in value.split(","):
        part = part.strip()
        if re.fullmatch(r"(?:chr)?\d+-(?:chr)?\d+", part, re.I):
            a, b = [int(re.sub(r"^chr", "", item, flags=re.I)) for item in part.split("-")]
            if not 1 <= a <= b <= 22:
                raise ValueError("Chromosome ranges must be within 1-22 and ascending")
            selected.update(str(i) for i in range(a, b + 1))
        else:
            selected.add(chromosome(part))
    return selected


def region_filter(value):
    if not value:
        return None
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", value)
    if not match:
        raise ValueError("--region must be chromosome:start-end, using 1-based closed coordinates")
    c, begin, end = chromosome(match[1]), int(match[2]), int(match[3])
    if not 1 <= begin <= end:
        raise ValueError("--region requires 1 <= start <= end")
    return c, begin, end


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).endswith((".gz", ".bgz")) else Path(path).open(encoding="utf-8", newline="")


def table(path, required):
    """Read a TSV without collapsing empty cells, which carry missingness."""
    names = None
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("##"):
                continue
            cells = raw.rstrip("\r\n").split("\t")
            if names is None:
                names = [v.lstrip("#").strip().lower() for v in cells]
                if len(set(names)) != len(names) or not set(required) <= set(names):
                    raise ValueError(f"{path}: need unique tab-separated columns {','.join(required)}")
                continue
            if len(cells) != len(names):
                raise ValueError(f"{path}:{line_number}: expected {len(names)} TSV cells, got {len(cells)}")
            yield line_number, dict(zip(names, cells))
    if names is None:
        raise ValueError(f"Empty table: {path}")


def load_targets(path, declared_build, maximum=5000, chromosomes=None, region=None):
    targets, ids, coordinates, counts = [], set(), set(), Counter()
    for line_number, row in table(path, ("chrom", "pos", "ref", "alts")):
        counts["input_rows"] += 1
        for column in ("build", "genome_build", "assembly", "target_build"):
            if row.get(column, "").strip().upper() not in MISSING:
                if genome_build(row[column]) != declared_build:
                    raise ValueError(f"{path}:{line_number}: target {column}={row[column]} conflicts with {declared_build}; no relabeling/liftover")
        c = chromosome(row["chrom"])
        if not re.fullmatch(r"[0-9]+", row["pos"].strip()):
            raise ValueError(f"{path}:{line_number}: pos must be a positive 1-based integer")
        pos, ref = int(row["pos"]), row["ref"].strip().upper()
        alts = [v.strip().upper() for v in row["alts"].split(",")]
        if pos < 1 or ref not in BASES or len(ref) != 1 or not 1 <= len(alts) <= 3 or len(set(alts)) != len(alts) or any(len(a) != 1 or a not in BASES or a == ref for a in alts):
            raise ValueError(f"{path}:{line_number}: targets require a 1-based SNV REF and 1-3 distinct nonreference A/C/G/T ALTs")
        if (chromosomes is not None and c not in chromosomes) or (region is not None and not (c == region[0] and region[1] <= pos <= region[2])):
            counts["filtered_rows"] += 1
            continue
        sid = row.get("site_id", "").strip() or f"{c}:{pos}:{ref}"
        if sid.upper() in MISSING or sid in ids or (c, pos) in coordinates:
            raise ValueError(f"{path}:{line_number}: selected targets require unique nonmissing site_id and coordinates")
        ids.add(sid)
        coordinates.add((c, pos))
        targets.append(dict(site_id=sid, chrom=c, pos=pos, ref=ref, alts=alts))
        if maximum and len(targets) > maximum:
            raise ValueError(f"Target cap {maximum} exceeded; shard --targets or explicitly increase --max-targets")
    if not targets:
        raise ValueError("No target sites remain after chromosome/region filters")
    counts["selected_sites"] = len(targets)
    return sorted(targets, key=lambda r: (chromosome_key(r["chrom"]), r["pos"])), dict(counts)


def boolean(value, name):
    if value.lower() in {"1", "true", "yes"}:
        return True
    if value.lower() in {"0", "false", "no"}:
        return False
    raise ValueError(f"{name} must be 0/1 or true/false")


def local_path(value, dir_data, config_dir):
    if not value:
        return ""
    value = value.replace("{dir_data}", str(Path(dir_data).expanduser().resolve()))
    remaining = value.replace("{chrom}", "")
    if "{" in remaining or "}" in remaining or re.match(r"[A-Za-z][A-Za-z0-9+.-]*://", value):
        raise ValueError(f"Atlas paths must be local files, with only {{dir_data}} and {{chrom}} placeholders: {value}")
    path = Path(value).expanduser()
    return str((path if path.is_absolute() else config_dir / path).resolve())


def default_resources(dir_data=None, all_resources=False):
    """Prefer finished chromosome VCFs, without scanning any score data.

    One existing data/index pair selects the VCF layout. Missing requested
    chromosomes then fail explicitly; a partial conversion is never silently
    mixed with the old whole-genome TSV source. Legacy descriptors remain
    available through --resources and are the fallback before conversion.
    """
    root = Path(dir_data or os.environ.get("DIR_DATA", configured_path("DIR_BIG"))).expanduser()
    directory = root / "annot" / "allele3" / "alphagenome" / "38"
    for chrom in [*[str(i) for i in range(1, 23)], "X", "Y", "M"]:
        path = directory / f"chr{chrom}.vcf.bgz"
        if path.is_file() and any(Path(str(path) + suffix).is_file() for suffix in (".tbi", ".csi")):
            return VCF_RESOURCES
    return DEFAULT_RESOURCES.with_name("alphagenome_resources.all3.tsv") if all_resources else DEFAULT_RESOURCES


def requested_chromosomes(args, targets=None):
    if targets is not None:
        return sorted({target["chrom"] for target in targets}, key=chromosome_key)
    if getattr(args, "targets", None):
        selected, _ = load_targets(Path(args.targets).expanduser().resolve(), "GRCh38",
                                   getattr(args, "max_targets", 5000),
                                   chromosome_filter(getattr(args, "chromosomes", "all")),
                                   region_filter(getattr(args, "region", None)))
        return sorted({target["chrom"] for target in selected}, key=chromosome_key)
    selected = chromosome_filter(getattr(args, "chromosomes", "all"))
    if selected:
        return sorted(selected, key=chromosome_key)
    raise ValueError("Atlas {chrom} templates require target sites or an explicit chromosome selection; no whole-genome scan is inferred")


def vcf_resource_schema(row):
    """Map original metric names to the converter's fixed lossless INFO schema."""
    if row["column_names"] == AVI_COLUMNS and row["kind"] == "avi":
        source = "avi"
    elif row["column_names"] == SPLICING_COLUMNS and row["kind"] == "metrics":
        source = "splicing"
    elif row["column_names"] == ATTRIBUTION_COLUMNS and row["kind"] == "metrics":
        source = "features"
    else:
        raise ValueError(f"Resource {row['resource_id']}: merged_vcf requires the recognized AVI, splicing or attribution schema")
    mapping = VCF_METRICS[source]
    if row["metrics"] != list(mapping):
        raise ValueError(f"Resource {row['resource_id']}: merged_vcf metric order must retain the original source columns")
    declared = row.get("info_fields", "").strip()
    expected = ",".join(f"{metric}={tag}" for metric, tag in mapping.items())
    if declared and declared != expected:
        raise ValueError(f"Resource {row['resource_id']}: INFO mapping conflicts with the merged AlphaGenome VCF schema")
    bit = VCF_SOURCE_BITS[source]
    if str(row.get("source_mask_bit", "")).strip() not in {"", str(bit)}:
        raise ValueError(f"Resource {row['resource_id']}: source_mask_bit conflicts with {source}")
    row.update(info_fields=expected, info_mapping=dict(mapping), source_mask_bit=bit, source_kind=source)


def resource_input_paths(spec):
    """Yield concrete data/index inputs, never an unresolved path template.

    Both implicit index candidates are retained so frozen plans and identities
    notice a later ambiguous .tbi/.csi sibling. Disabled rows have no inputs.
    """
    if not spec.get("enabled"):
        return
    seen = set()
    for item in spec.get("physical_files", [spec]):
        paths = [item.get("path"), item.get("index_path")]
        if item.get("path") and not item.get("index_explicit", bool(item.get("index_path"))):
            paths.extend(item["path"] + suffix for suffix in (".tbi", ".csi"))
        for path in paths:
            if path and path not in seen:
                seen.add(path)
                yield path


def choose_index(path, explicit):
    if explicit:
        if not str(explicit).endswith((".tbi", ".csi")):
            raise ValueError("--avi-index/index_path must name the matching .tbi or .csi file")
        return explicit
    existing = [str(path) + suffix for suffix in (".tbi", ".csi") if Path(str(path) + suffix).is_file()]
    if len(existing) > 1:
        raise ValueError(f"Both .tbi and .csi exist for {path}; select the matching index explicitly")
    return existing[0] if existing else ""


def load_resources(args, targets=None):
    config = Path(args.resources or default_resources(args.dir_data)).expanduser().resolve()
    required = ("resource_id", "kind", "path", "build", "coordinate_system", "columns", "chrom_column", "pos_column", "ref_column", "alt_column", "metric_columns", "required", "enabled")
    specs, ids, missing = [], set(), []
    for line_number, raw in table(config, required):
        row = {k: v.strip() for k, v in raw.items()}
        rid = row["resource_id"]
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", rid) or rid in ids:
            raise ValueError(f"{config}:{line_number}: resource_id must be unique and filesystem-safe")
        ids.add(rid)
        if row["kind"] not in {"avi", "metrics"}:
            raise ValueError(f"Resource {rid}: kind must be avi or metrics")
        row["required"] = boolean(row["required"], "required")
        row["enabled"] = boolean(row["enabled"], "enabled")
        if row["enabled"] and getattr(args, "require_all_resources", False):
            row["required"] = True
        row["path"] = local_path(row["path"], args.dir_data, config.parent)
        row["index_path"] = local_path(row.get("index_path", ""), args.dir_data, config.parent)
        row["storage_format"] = row.get("storage_format", "").strip() or "tsv"
        if row["storage_format"] not in {"tsv", "merged_vcf"}:
            raise ValueError(f"Resource {rid}: storage_format must be tsv or merged_vcf")
        if row["kind"] == "avi":
            if args.avi_file:
                row["path"] = local_path(args.avi_file, args.dir_data, config.parent)
                row["index_path"] = ""  # A configured index must never follow a different data override.
                row["storage_format"] = "merged_vcf" if str(args.avi_file).endswith((".vcf.gz", ".vcf.bgz")) else "tsv"
                row["info_fields"], row["source_mask_bit"] = "", ""
            if args.avi_index:
                row["index_path"] = local_path(args.avi_index, args.dir_data, config.parent)
            if args.source_version:
                row["source_version"] = args.source_version
            row["enabled"], row["required"] = True, True
        row.setdefault("source_version", "user-declared-unverified")
        row["source_version"] = row["source_version"] or "user-declared-unverified"
        row.update(resource_semantics(row))
        if not row["enabled"]:
            row["status"] = "disabled"
            specs.append(row)
            continue
        if genome_build(row["build"]) != "GRCh38" or row["coordinate_system"] != "1-based":
            raise ValueError(f"Resource {rid}: only explicitly declared GRCh38 1-based SNV resources are supported")
        row["build"] = "GRCh38"
        columns = [v.strip() for v in row["columns"].split(",")]
        metrics = [v.strip() for v in row["metric_columns"].split(",")]
        key_columns = [row[v + "_column"] for v in ("chrom", "pos", "ref", "alt")]
        if any(not v for v in columns) or len({v.lower() for v in columns}) != len(columns) or not metrics or any(not v for v in metrics) or len(set(metrics)) != len(metrics) or len(set(key_columns)) != 4 or not set(key_columns + metrics) <= set(columns) or set(key_columns) & set(metrics):
            raise ValueError(f"Resource {rid}: invalid explicit column/metric schema")
        if row["kind"] == "avi" and (columns != AVI_COLUMNS or key_columns != AVI_COLUMNS[:4] or metrics != AVI_COLUMNS[4:]):
            raise ValueError("AVI must use the official six-column CHROM,POS,REF,ALT,raw_score,PHRED schema; generic tables use kind=metrics")
        row["column_names"], row["metrics"] = columns, metrics
        row["column_indices"] = {v: columns.index(row[v + "_column"]) for v in ("chrom", "pos", "ref", "alt")}
        if row["storage_format"] == "merged_vcf":
            vcf_resource_schema(row)
        templated = "{chrom}" in row["path"]
        if "{chrom}" in row["index_path"] and not templated:
            raise ValueError(f"Resource {rid}: index_path uses {{chrom}} but path does not")
        if templated and row["index_path"] and "{chrom}" not in row["index_path"]:
            raise ValueError(f"Resource {rid}: a chromosome path template needs a matching chromosome index template")
        selections = requested_chromosomes(args, targets) if templated else [None]
        row["physical_files"] = []
        for chrom in selections:
            label = "M" if chrom == "MT" else chrom
            path = row["path"].replace("{chrom}", label) if templated else row["path"]
            explicit_index = row["index_path"].replace("{chrom}", label) if templated else row["index_path"]
            item = dict(chrom=chrom, path=path, index_explicit=bool(explicit_index))
            if not path or not Path(path).is_file():
                if row["required"]:
                    missing.append(f"Required Atlas file missing ({rid}): {path or 'path not configured in ' + str(config)}")
                item["status"] = "deferred_missing_file"
            selected_index = choose_index(path, explicit_index)
            if not selected_index or not Path(selected_index).is_file():
                if row["required"]:
                    expected = explicit_index or f"{path}.tbi OR {path}.csi"
                    missing.append(f"Matching Tabix .tbi/.csi missing ({rid}): {expected}")
                item.setdefault("status", "deferred_missing_index")
            else:
                item.setdefault("status", "ready")
            item["index_path"] = selected_index
            row["physical_files"].append(item)
        row["status"] = "ready" if any(item["status"] == "ready" for item in row["physical_files"]) else row["physical_files"][0]["status"]
        row["has_exclusions"] = any(item["status"] != "ready" for item in row["physical_files"])
        if not templated:
            row["index_path"] = row["physical_files"][0]["index_path"]
        elif row["storage_format"] != "merged_vcf":
            raise ValueError(f"Resource {rid}: {{chrom}} templates currently require storage_format=merged_vcf")
        specs.append(row)
    if sum(r["kind"] == "avi" for r in specs) != 1:
        raise ValueError("Resource manifest requires exactly one AVI row")
    if missing:
        raise ValueError("\n".join(missing) + "\nUse split2vcf.sh for chromosome BGZF VCFs or an indexed legacy BGZF TSV. Ordinary gzip and ZIP archives cannot be queried directly. Override the manifest with --resources if needed.")
    return specs


def compression_metadata(path):
    """Check only the first compressed block header, without decompressing a source.

    A '.gz' suffix does not establish BGZF. Ordinary gzip cannot be queried
    randomly with Tabix; no attempt is made to rewrite or re-index any input.
    This bounded check does not certify later blocks or a whole-file checksum.
    """
    with Path(path).open("rb") as handle:
        prefix = handle.read(12)
        if len(prefix) < 12 or prefix[:3] != b"\x1f\x8b\x08":
            raise ValueError(f"Atlas source is not a BGZF gzip stream: {path}")
        if not prefix[3] & 4:
            raise ValueError(f"Atlas source uses ordinary gzip, not BGZF: {path}. Use the downloaded BGZF TSV with its matching Tabix index; an ordinary .gz cannot support indexed lookup.")
        extra_length = struct.unpack("<H", prefix[10:12])[0]
        extra = handle.read(extra_length)
        if len(extra) != extra_length:
            raise ValueError(f"Truncated BGZF extra header: {path}")
    offset, block_size = 0, None
    while offset < len(extra):
        if offset + 4 > len(extra):
            raise ValueError(f"Malformed gzip extra header: {path}")
        name, length = extra[offset:offset + 2], struct.unpack("<H", extra[offset + 2:offset + 4])[0]
        offset += 4
        if offset + length > len(extra):
            raise ValueError(f"Truncated gzip extra field: {path}")
        if name == b"BC":
            if length != 2 or block_size is not None:
                raise ValueError(f"Invalid BGZF block-size field: {path}")
            block_size = struct.unpack("<H", extra[offset:offset + 2])[0] + 1
        offset += length
    if block_size is None:
        raise ValueError(f"Atlas source uses ordinary gzip without a BGZF BC field: {path}. Random-access Tabix requires the BGZF source and its matching index.")
    if block_size < 22 + extra_length or block_size > 65536:
        raise ValueError(f"Invalid first BGZF block size: {path}")
    return dict(format="BGZF", first_block_bytes=block_size, verification="first_block_header_only_no_full_scan")


def index_metadata(path, spec):
    """Read only the bounded TBI/CSI prefix, including the Tabix coordinate flag."""
    with Path(path).open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rb") as handle:
        magic = handle.read(4)
        if magic == b"TBI\x01":
            data = handle.read(32)
            if len(data) != 32:
                raise ValueError(f"Truncated Tabix index: {path}")
            n_ref, fmt, seq, begin, end, meta, skip, name_bytes = struct.unpack("<8i", data)
            kind = "TBI"
            if n_ref < 1 or name_bytes < 1:
                raise ValueError(f"Empty/invalid Tabix index: {path}")
        elif magic == b"CSI\x01":
            prefix = handle.read(12)
            if len(prefix) != 12:
                raise ValueError(f"Truncated CSI index: {path}")
            min_shift, depth, aux_len = struct.unpack("<3i", prefix)
            data = handle.read(28)
            if min_shift < 0 or depth < 0 or aux_len < 28 or len(data) != 28:
                raise ValueError(f"CSI lacks valid Tabix column metadata: {path}")
            fmt, seq, begin, end, meta, skip, name_bytes = struct.unpack("<7i", data)
            kind = "CSI"
        else:
            raise ValueError(f"Not a TBI/CSI Tabix index: {path}")
    if fmt & 0x10000:
        raise ValueError(f"Resource {spec['resource_id']}: index uses 0-based coordinates; expected a 1-based SNV index")
    if fmt not in {0, 2}:
        raise ValueError(f"Resource {spec['resource_id']}: expected generic TSV or SNV-compatible VCF Tabix format, got preset {fmt}")
    merged_vcf = spec.get("storage_format") == "merged_vcf"
    if merged_vcf:
        if fmt != 2 or (seq, begin, end) != (1, 2, 0):
            raise ValueError(f"Resource {spec['resource_id']}: merged VCF requires a VCF-preset index with columns (1,2,0)")
    elif fmt == 2:
        # Official AVI/splicing SNV downloads use TBX_VCF indexes over the
        # CHROM,POS,REF,ALT score table. The fourth field is a one-base ALT,
        # so HTSlib's VCF interval decoder still yields a one-base interval.
        # Do not generalize this to reordered schemas or non-SNV records;
        # query() also checks the exact returned coordinate and REF/ALT bases.
        columns = [spec[key + "_column"] for key in ("chrom", "pos", "ref", "alt")]
        if spec["column_names"][:4] != columns or (seq, begin, end) != (1, 2, 0):
            raise ValueError(f"Resource {spec['resource_id']}: VCF-preset SNV index requires leading CHROM,POS,REF,ALT and index columns (1,2,0)")
    expected_seq, expected_pos = spec["column_indices"]["chrom"] + 1, spec["column_indices"]["pos"] + 1
    if seq != expected_seq or begin != expected_pos or end not in {0, expected_pos}:
        raise ValueError(f"Resource {spec['resource_id']}: Tabix index columns {(seq, begin, end)} conflict with SNV schema {(expected_seq, expected_pos, expected_pos)}")
    if skip < 0 or name_bytes < 1:
        raise ValueError(f"Resource {spec['resource_id']}: invalid Tabix metadata")
    return dict(format=kind, preset="merged_vcf_snv" if merged_vcf else "vcf_snv" if fmt == 2 else "generic",
                coordinate_system="1-based", chrom_column=seq, start_column=begin, end_column=end,
                meta_character=meta, skip_lines=skip)


def normalized_column(value):
    value = re.sub(r"[^a-z0-9]", "", value.lower().lstrip("#"))
    aliases = {"chr": "chrom", "chromosome": "chrom", "position": "pos", "start": "pos",
               "reference": "ref", "alternate": "alt", "raw": "rawscore", "aviraw": "rawscore",
               "avirawscore": "rawscore", "avi": "rawscore", "aviscore": "rawscore",
               "aviphred": "phred", "aviphredscore": "phred", "phredscore": "phred"}
    return aliases.get(value, value)


def check_headers(lines, spec):
    if spec.get("storage_format") == "merged_vcf":
        return check_vcf_headers(lines, spec)
    recognized_builds, observed_header = set(), None
    for line in lines:
        text = line.decode() if isinstance(line, bytes) else line
        lowered = text.lower()
        if re.search(r"(?:reference|genome_build|assembly|build)\s*[:=]", lowered):
            recognized_builds.update(genome_build(v) for v in re.findall(r"\b(?:grch37|grch38|hg19|hg38)\b", lowered))
        if re.search(r"(?:coordinate_system|coordinates)\s*[:=]\s*['\"]?(?:0-based|zero-based|bed)\b", lowered):
            raise ValueError(f"Resource {spec['resource_id']}: header declares zero-based coordinates")
        cells = text.rstrip("\r\n").lstrip("#").split("\t")
        if cells and normalized_column(cells[0]) == normalized_column(spec["column_names"][0]):
            expected = spec["column_names"]
            if spec["kind"] == "avi":
                matches = [normalized_column(v) for v in cells] == [normalized_column(v) for v in expected]
            else:
                matches = [v.lower() for v in cells] == [v.lower() for v in expected]
            if not matches:
                raise ValueError(f"Resource {spec['resource_id']}: source column header conflicts with the explicit schema")
            observed_header = cells
    if recognized_builds - {"GRCh38"}:
        raise ValueError(f"Resource {spec['resource_id']}: source header build {sorted(recognized_builds)} conflicts with GRCh38")
    return dict(recognized_builds=sorted(recognized_builds), columns=observed_header,
                verification="GRCh38_header_and_declared" if recognized_builds else "declared_only_no_recognized_build_header")


def check_vcf_headers(lines, spec):
    """Require the fixed site-only annotation schema before decoding INFO arrays."""
    recognized_builds, columns, fileformat, declarations = set(), None, None, {}
    for line in lines:
        text = line.decode() if isinstance(line, bytes) else line
        lowered = text.lower()
        if re.search(r"(?:reference|genome_build|assembly|build)\s*[:=]", lowered):
            recognized_builds.update(genome_build(v) for v in re.findall(r"\b(?:grch37|grch38|hg19|hg38)\b", lowered))
        if re.search(r"(?:coordinate_system|coordinates)\s*[:=]\s*['\"]?(?:0-based|zero-based|bed)\b", lowered):
            raise ValueError(f"Resource {spec['resource_id']}: header declares zero-based coordinates")
        if text.startswith("##fileformat="):
            if fileformat is not None:
                raise ValueError(f"Resource {spec['resource_id']}: duplicate VCF fileformat header")
            fileformat = text.strip().split("=", 1)[1]
        if text.startswith("##INFO=<"):
            # Read only the machine fields; quoted Description may contain commas.
            fields = dict(re.findall(r"(?:^##INFO=<|,)(ID|Number|Type)=([^,>]+)", text))
            if "ID" not in fields or not {"Number", "Type"} <= fields.keys():
                raise ValueError(f"Resource {spec['resource_id']}: malformed VCF INFO declaration")
            tag = fields["ID"]
            if tag in declarations:
                raise ValueError(f"Resource {spec['resource_id']}: duplicate VCF INFO declaration {tag}")
            declarations[tag] = fields
        if text.startswith("#CHROM\t"):
            if columns is not None:
                raise ValueError(f"Resource {spec['resource_id']}: duplicate VCF column header")
            columns = text.rstrip("\r\n").lstrip("#").split("\t")
    if fileformat not in {"VCFv4.2", "VCFv4.3", "VCFv4.4", "VCFv4.5"} or columns != VCF_COLUMNS:
        raise ValueError(f"Resource {spec['resource_id']}: merged VCF requires a VCFv4 site-only CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO header")
    if recognized_builds != {"GRCh38"}:
        raise ValueError(f"Resource {spec['resource_id']}: merged VCF source header build must explicitly declare GRCh38, got {sorted(recognized_builds)}")
    expected = {tag: "Float" for mapping in VCF_METRICS.values() for tag in mapping.values()}
    expected["AG_SOURCE_MASK"] = "Integer"
    for tag, dtype in expected.items():
        if tag not in declarations or declarations[tag]["Number"] != "A" or declarations[tag]["Type"] != dtype:
            raise ValueError(f"Resource {spec['resource_id']}: merged VCF INFO/{tag} requires Number=A,Type={dtype}")
    return dict(recognized_builds=["GRCh38"], columns=columns, fileformat=fileformat,
                verification="GRCh38_header_and_declared", info_number="A", source_mask_bits=VCF_SOURCE_BITS,
                expected_info_fields=sorted(expected))


def parse_metric(value, metric, spec):
    value = value.strip()
    if value.upper() in MISSING:
        return None
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"Resource {spec['resource_id']}: {metric} is not numeric: {value!r}") from exc
    if not math.isfinite(number) or spec["kind"] == "avi" and metric == "PHRED" and number < 0:
        raise ValueError(f"Resource {spec['resource_id']}: invalid {metric} value {value!r}")
    if spec.get("metric_semantics") == "reference_relative_splicing_prediction" and number < 0:
        raise ValueError(f"Resource {spec['resource_id']}: splicing disruption summary must be nonnegative, got {value!r}")
    return value


class TSVAtlasReader:
    def __init__(self, spec):
        try:
            import pysam
        except ImportError as exc:
            raise RuntimeError("atlas requires pysam for bounded local Tabix reads; install the declared CPU environment") from exc
        self.spec, self.audit = spec, Counter()
        spec["compression_check"] = compression_metadata(spec["path"])
        spec["index_check"] = index_metadata(spec["index_path"], spec)
        self.handle = pysam.TabixFile(spec["path"], index=spec["index_path"])
        try:
            spec["header_check"] = check_headers(self.handle.header, spec)
            self.contigs = {}
            for actual in self.handle.contigs:
                try:
                    canonical = chromosome(actual)
                except ValueError:
                    continue  # Noncanonical contigs are outside the target contract.
                if canonical in self.contigs:
                    raise ValueError(f"Resource {spec['resource_id']}: ambiguous indexed contig aliases for {canonical}")
                self.contigs[canonical] = actual
        except BaseException:
            self.handle.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.handle.close()

    def query(self, target):
        return self.query_batch([target])[target["pos"]]

    def query_batch(self, targets):
        """Fetch one bounded interval and retain only exact requested substitutions."""
        if not targets:
            return {}
        spec, c = self.spec, targets[0]["chrom"]
        requested = {target["pos"]: target for target in targets}
        if len(requested) != len(targets) or any(target["chrom"] != c for target in targets):
            raise ValueError("Atlas query_batch requires distinct target positions on one chromosome")
        begin, end = min(requested), max(requested)
        found_by_pos = {pos: {} for pos in requested}
        self.audit["target_queries"] += len(targets)
        if c not in self.contigs:
            self.audit["absent_contig_queries"] += len(targets)
            return {pos: ({}, "not_indexed_contig", "unverified_no_returned_rows") for pos in requested}
        self.audit["tabix_fetches"] += 1
        self.audit["fetched_interval_bases"] += end - begin + 1
        self.audit["largest_fetch_interval_bases"] = max(self.audit["largest_fetch_interval_bases"], end - begin + 1)
        for line in self.handle.fetch(self.contigs[c], begin - 1, end):
            self.audit["returned_rows"] += 1
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) != len(spec["column_names"]):
                raise ValueError(f"Resource {spec['resource_id']}: expected {len(spec['column_names'])} columns in {c}:{begin}-{end}")
            values = dict(zip(spec["column_names"], cells))
            chrom = chromosome(values[spec["chrom_column"]])
            pos_text = values[spec["pos_column"]]
            if not re.fullmatch(r"[0-9]+", pos_text) or int(pos_text) < 1:
                raise ValueError(f"Resource {spec['resource_id']}: invalid 1-based position {pos_text!r}")
            pos = int(pos_text)
            if chrom != c or not begin <= pos <= end:
                raise ValueError(f"Resource {spec['resource_id']}: fetched row coordinate differs from the bounded query; index/data mismatch")
            ref, alt = values[spec["ref_column"]], values[spec["alt_column"]]
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES or ref == alt:
                raise ValueError(f"Resource {spec['resource_id']}: source is not a canonical nonreference SNV at {chrom}:{pos}")
            if pos not in requested:
                self.audit["off_target_rows_discarded"] += 1
                continue
            target, found = requested[pos], found_by_pos[pos]
            self.audit["exact_target_rows"] += 1
            if ref != target["ref"]:
                raise ValueError(f"Resource {spec['resource_id']}: REF mismatch at {chrom}:{pos}: target={target['ref']} atlas={ref}; verify assembly and reference")
            metrics = {m: parse_metric(values[m], m, spec) for m in spec["metrics"]}
            if alt in found:
                old = {k: None if v is None else float(v) for k, v in found[alt].items()}
                new = {k: None if v is None else float(v) for k, v in metrics.items()}
                if old != new:
                    raise ValueError(f"Resource {spec['resource_id']}: conflicting duplicate substitution {chrom}:{pos}:{ref}:{alt}")
                self.audit["identical_duplicate_rows"] += 1
            else:
                found[alt] = metrics
        self.audit["queried_substitutions_found"] += sum(len(found) for found in found_by_pos.values())
        return {pos: (found, "not_scored", "exact_REF_checked" if found else "unverified_no_returned_rows")
                for pos, found in found_by_pos.items()}


class PhysicalVCFReader:
    """One open chromosome source and one decoded target batch, shared by views.

    Only requested positions are retained. Score strings are decoded from Tabix
    text rather than VariantFile float32 values, preserving the TSV precision.
    """
    def __init__(self, item, spec):
        try:
            import pysam
        except ImportError as exc:
            raise RuntimeError("atlas requires pysam for bounded local Tabix reads; install the declared CPU environment") from exc
        self.path, self.audit = item["path"], Counter()
        self.spec = spec
        self.compression_check = compression_metadata(self.path)
        self.index_check = index_metadata(item["index_path"], spec)
        self.handle = pysam.TabixFile(self.path, index=item["index_path"])
        self.cache_key, self.cache = None, None
        try:
            self.header_check = check_vcf_headers(self.handle.header, spec)
            self.contigs = {}
            for actual in self.handle.contigs:
                try:
                    canonical = chromosome(actual)
                except ValueError:
                    continue
                if canonical in self.contigs:
                    raise ValueError(f"Resource {spec['resource_id']}: ambiguous indexed contig aliases for {canonical}")
                self.contigs[canonical] = actual
            if item.get("chrom") and set(self.contigs) != {item["chrom"]}:
                raise ValueError(f"Resource {spec['resource_id']}: chromosome shard {self.path} indexes {sorted(self.contigs)} instead of {item['chrom']}")
        except BaseException:
            self.handle.close()
            raise

    def close(self):
        self.cache_key, self.cache = None, None
        self.handle.close()

    def _decode(self, cells, target):
        c, pos, ref = chromosome(cells[0]), int(cells[1]), cells[3]
        if ref != target["ref"]:
            raise ValueError(f"Resource {self.spec['resource_id']}: REF mismatch at {c}:{pos}: target={target['ref']} atlas={ref}; verify assembly and reference")
        alts = cells[4].split(",")
        info = {}
        for field in cells[7].split(";"):
            if not field or field == ".":
                continue
            key, separator, value = field.partition("=")
            if key in info:
                raise ValueError(f"Resource {self.spec['resource_id']}: duplicate INFO/{key} at {c}:{pos}")
            info[key] = value if separator else None
        mask = info.get("AG_SOURCE_MASK")
        if mask is None:
            raise ValueError(f"Resource {self.spec['resource_id']}: missing INFO/AG_SOURCE_MASK at {c}:{pos}")
        masks = mask.split(",")
        if len(masks) != len(alts) or any(not re.fullmatch(r"[1-7]", value) for value in masks):
            raise ValueError(f"Resource {self.spec['resource_id']}: AG_SOURCE_MASK must have one integer in 1..7 per ALT at {c}:{pos}")
        masks = [int(value) for value in masks]
        decoded = {source: {} for source in VCF_METRICS}
        for source, mapping in VCF_METRICS.items():
            semantic_spec = dict(resource_id=self.spec["resource_id"], kind="avi" if source == "avi" else "metrics",
                                 metric_semantics="reference_relative_splicing_prediction" if source == "splicing" else "signed_shap_contribution_to_avi")
            by_allele = {alt: {} for i, alt in enumerate(alts) if masks[i] & VCF_SOURCE_BITS[source]}
            for metric, tag in mapping.items():
                raw = info.get(tag, ".")
                if raw is None:
                    raise ValueError(f"Resource {self.spec['resource_id']}: INFO/{tag} cannot be a flag at {c}:{pos}")
                values = ["."] * len(alts) if raw == "." else raw.split(",")
                if len(values) != len(alts):
                    raise ValueError(f"Resource {self.spec['resource_id']}: INFO/{tag} Number=A has {len(values)} values for {len(alts)} ALTs at {c}:{pos}")
                for i, (alt, value) in enumerate(zip(alts, values)):
                    parsed = parse_metric(value, metric, semantic_spec)
                    if masks[i] & VCF_SOURCE_BITS[source]:
                        by_allele[alt][metric] = parsed
                    elif parsed is not None:
                        raise ValueError(f"Resource {self.spec['resource_id']}: INFO/{tag} has a score but AG_SOURCE_MASK excludes {source} at {c}:{pos}:{ref}:{alt}")
            decoded[source] = by_allele
        return decoded

    def query_batch(self, targets):
        c = targets[0]["chrom"]
        requested = {target["pos"]: target for target in targets}
        if len(requested) != len(targets) or any(target["chrom"] != c for target in targets):
            raise ValueError("Atlas query_batch requires distinct target positions on one chromosome")
        key = (c, tuple((target["pos"], target["ref"]) for target in targets))
        if key == self.cache_key:
            self.audit["shared_batch_cache_hits"] += 1
            return self.cache, True
        self.audit["target_queries"] += len(targets)
        found = {pos: None for pos in requested}
        if c not in self.contigs:
            self.audit["absent_contig_queries"] += len(targets)
            result = (found, "not_indexed_contig")
        else:
            begin, end = min(requested), max(requested)
            self.audit["tabix_fetches"] += 1
            self.audit["fetched_interval_bases"] += end - begin + 1
            self.audit["largest_fetch_interval_bases"] = max(self.audit["largest_fetch_interval_bases"], end - begin + 1)
            last_position = None
            for line in self.handle.fetch(self.contigs[c], begin - 1, end):
                self.audit["returned_site_rows"] += 1
                cells = line.rstrip("\r\n").split("\t")
                if len(cells) != len(VCF_COLUMNS):
                    raise ValueError(f"Resource {self.spec['resource_id']}: expected 8 site-only VCF columns in {c}:{begin}-{end}")
                if not re.fullmatch(r"[0-9]+", cells[1]) or int(cells[1]) < 1:
                    raise ValueError(f"Resource {self.spec['resource_id']}: invalid 1-based VCF position {cells[1]!r}")
                pos, ref, alts = int(cells[1]), cells[3], cells[4].split(",")
                if chromosome(cells[0]) != c or not begin <= pos <= end:
                    raise ValueError(f"Resource {self.spec['resource_id']}: fetched VCF coordinate differs from bounded query; index/data mismatch")
                if last_position is not None and pos <= last_position:
                    raise ValueError(f"Resource {self.spec['resource_id']}: merged VCF requires one sorted site record per position; duplicate/out-of-order {c}:{pos}")
                last_position = pos
                if len(ref) != 1 or ref not in BASES or not 1 <= len(alts) <= 3 or len(set(alts)) != len(alts) or any(len(alt) != 1 or alt not in BASES or alt == ref for alt in alts):
                    raise ValueError(f"Resource {self.spec['resource_id']}: merged VCF is not a canonical nonreference SNV at {c}:{pos}")
                if pos not in requested:
                    self.audit["off_target_site_rows_discarded"] += 1
                    continue
                found[pos] = self._decode(cells, requested[pos])
                self.audit["exact_target_site_rows"] += 1
            result = (found, "not_scored")
        self.cache_key, self.cache = key, result
        return result, False


class VCFReaderPool:
    def __init__(self):
        self.readers = {}
        self.active_reader = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        for reader in self.readers.values():
            reader.close()

    def reader(self, item, spec):
        key = (item["path"], item["index_path"])
        if key not in self.readers:
            self.readers[key] = PhysicalVCFReader(item, spec)
        reader = self.readers[key]
        if item.get("chrom") and set(reader.contigs) != {item["chrom"]}:
            raise ValueError(f"Resource {spec['resource_id']}: a shared VCF shard has incompatible chromosome aliases")
        return reader

    def audits(self):
        return {path: dict(reader.audit) for (path, _), reader in sorted(self.readers.items())}

    def activate(self, reader):
        if self.active_reader is not None and self.active_reader is not reader:
            self.active_reader.cache_key, self.active_reader.cache = None, None
        self.active_reader = reader


class MergedVCFAtlasReader:
    def __init__(self, spec, pool=None):
        self.spec, self.audit, self.owns_pool = spec, Counter(), pool is None
        self.pool = pool or VCFReaderPool()
        self.sources = {}
        try:
            for item in spec["physical_files"]:
                if item["status"] == "ready":
                    self.sources[item["chrom"]] = self.pool.reader(item, spec)
            formats = sorted({reader.index_check["format"] for reader in self.sources.values()})
            spec["compression_check"] = dict(format="BGZF", verification="first_block_header_only_no_full_scan",
                                              files={reader.path: reader.compression_check for reader in self.sources.values()})
            spec["index_check"] = dict(format=formats[0] if len(formats) == 1 else "mixed", preset="merged_vcf_snv", coordinate_system="1-based",
                                        files={reader.path: reader.index_check for reader in self.sources.values()})
            spec["header_check"] = dict(verification="GRCh38_header_and_declared", recognized_builds=["GRCh38"], columns=VCF_COLUMNS,
                                         files={reader.path: reader.header_check for reader in self.sources.values()})
        except BaseException:
            if self.owns_pool:
                self.pool.__exit__()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *_):
        if self.owns_pool:
            self.pool.__exit__()

    def query(self, target):
        return self.query_batch([target])[target["pos"]]

    def query_batch(self, targets):
        if not targets:
            return {}
        self.audit["target_queries"] += len(targets)
        source = self.sources.get(targets[0]["chrom"]) or self.sources.get(None)
        if source is None:
            self.audit["absent_chromosome_file_queries"] += len(targets)
            return {target["pos"]: ({}, "not_indexed_contig", "unverified_no_returned_rows") for target in targets}
        before = source.audit.copy()
        self.pool.activate(source)
        (found, missing_status), cached = source.query_batch(targets)
        if cached:
            self.audit["shared_batch_cache_hits"] += 1
        else:
            for key, value in source.audit.items():
                if key == "largest_fetch_interval_bases":
                    self.audit[key] = max(self.audit[key], value)
                elif key != "target_queries":
                    self.audit[key] += value - before[key]
        result = {}
        for pos, entry in found.items():
            values = {} if entry is None else entry[self.spec["source_kind"]]
            self.audit["queried_substitutions_found"] += len(values)
            result[pos] = (values, missing_status, "exact_REF_checked" if entry is not None else "unverified_no_returned_rows")
        return result


def AtlasReader(spec, pool=None):
    """Open one logical resource while sharing physical merged VCF reads."""
    return MergedVCFAtlasReader(spec, pool) if spec.get("storage_format") == "merged_vcf" else TSVAtlasReader(spec)


def query_batches(targets, gap_bp=64, span_bp=16384):
    """Coalesce nearby sorted targets without ever fetching a chromosome-wide range."""
    if gap_bp < 0 or span_bp < 1:
        raise ValueError("Atlas query gap must be >= 0 and maximum query span must be >= 1")
    batch = []
    for target in targets:
        if batch and (target["chrom"] != batch[-1]["chrom"] or target["pos"] - batch[-1]["pos"] > gap_bp
                      or target["pos"] - batch[0]["pos"] + 1 > span_bp):
            yield batch
            batch = []
        if batch and target["pos"] <= batch[-1]["pos"]:
            raise ValueError("Atlas targets must be sorted and have unique positions within each chromosome")
        batch.append(target)
    if batch:
        yield batch


def score_status(metrics, absent_status):
    if metrics is None:
        return absent_status
    available = sum(v is not None for v in metrics.values())
    return "available" if available == len(metrics) else "available_partial" if available else "score_missing"


def writer(handle, fields):
    result = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    result.writeheader()
    return result


def protected_inputs(request):
    return [InputPath(path) for path in sorted({path for spec in request["resources"] for path in resource_input_paths(spec)})]


def query_targets(request, readers, targets=None):
    """Yield exact results using a bounded interval cache shared by nearby targets."""
    parameters = request["parameters"]
    for batch in query_batches(request["targets"] if targets is None else targets,
                               parameters.get("query_gap_bp", 64), parameters.get("query_span_bp", 16384)):
        results = {rid: reader.query_batch(batch) for rid, reader in readers.items()}
        for target in batch:
            yield target, {rid: entries[target["pos"]] for rid, entries in results.items()}


def preflight(request, maximum=16):
    """Read-only format/index/header/REF checks on a deterministic target subset."""
    targets = request["targets"]
    n = min(len(targets), maximum)
    indices = [0] if n == 1 else [i * (len(targets) - 1) // (n - 1) for i in range(n)]
    selected = [targets[i] for i in indices]
    counts = Counter(target_sites=len(targets), sampled_target_sites=n)
    readers = {}
    with ExitStack() as stack:
        pool = stack.enter_context(VCFReaderPool())
        for spec in request["resources"]:
            if spec["status"] == "ready":
                readers[spec["resource_id"]] = stack.enter_context(AtlasReader(spec, pool))
                counts["resources_checked"] += 1
                counts["resources_with_missing_shards"] += bool(spec.get("has_exclusions"))
                print(f"[preflight] atlas={spec['resource_id']} compression=BGZF index={spec['index_check']['format']}:1-based semantics={spec['metric_semantics']}", flush=True)
            elif spec["enabled"]:
                counts["resources_deferred"] += 1
                print(f"[preflight] atlas={spec['resource_id']} status={spec['status']}", flush=True)
        for _, queries in query_targets(request, readers, selected):
            for rid, (entries, _, ref_check) in queries.items():
                counts[rid + ":substitutions_found"] += len(entries)
                counts[rid + ":sites_REF_verified"] += ref_check == "exact_REF_checked"
    return dict(status="preflight_with_exclusions" if counts["resources_deferred"] or counts["resources_with_missing_shards"] else "preflight_passed",
                schema=VERSION, semantics_version=SEMANTICS_VERSION, counts=dict(counts), reused=False,
                resources=request["resources"], query_audits={rid: dict(reader.audit) for rid, reader in readers.items()},
                physical_query_audits=pool.audits(), storage_version=STORAGE_VERSION,
                limitations=LIMITS + ["Preflight sampled only the recorded target subset and compressed headers. It did not verify every source row, exhaustively check index pairing or publish analysis files."],
                read_only=True, source_identity=request["source_identity"])


@transactional("manifest.json", input_files=("targets_path", "resources_path"), input_prefixes=(),
               resolve_inputs=lambda raw: protected_inputs(raw["request"]))
def publish(out, targets_path, resources_path, request):
    out = Path(out)
    specs = request["resources"]
    counts = Counter(target_sites=len(request["targets"]))
    readers = {}
    with ExitStack() as stack:
        pool = stack.enter_context(VCFReaderPool())
        for spec in specs:
            if spec["status"] == "ready":
                reader = stack.enter_context(AtlasReader(spec, pool))
                readers[spec["resource_id"]] = reader
                print(f"[build-check] atlas={spec['resource_id']} declared=GRCh38 source={spec['header_check']['verification']} index={spec['index_check']['format']}:1-based exact_REF=checked_per_returned_target", flush=True)
            print(f"[dataset] atlas={spec['resource_id']} status={spec['status']} source_version={spec['source_version']} file={spec['path'] or 'not_configured'}", flush=True)
        avi_spec = next(s for s in specs if s["kind"] == "avi")
        vh = stack.enter_context(gzip.open(out / "atlas_variants.tsv.gz", "wt", encoding="utf-8", newline=""))
        mh = stack.enter_context(gzip.open(out / "atlas_metrics.tsv.gz", "wt", encoding="utf-8", newline=""))
        ah = stack.enter_context((out / "target_availability.tsv").open("w", encoding="utf-8", newline=""))
        vw, mw, aw = writer(vh, VARIANT_FIELDS), writer(mh, METRIC_FIELDS), writer(ah, AVAILABILITY_FIELDS)
        for target, queries in query_targets(request, readers):
            common = dict(site_id=target["site_id"], chrom=target["chrom"], pos=target["pos"], build="GRCh38", ref=target["ref"], target_alts=",".join(target["alts"]))
            found, absent_status, ref_check = queries[avi_spec["resource_id"]]
            available, complete = 0, 0
            for allele in [target["ref"], *[base for base in BASES if base != target["ref"]]]:
                reference = allele == target["ref"]
                record = dict(common, allele=allele, allele_role="REF" if reference else "ALT",
                              observed_from_targets="NA" if reference else int(allele in target["alts"]),
                              resource_id=avi_spec["resource_id"], source_version=avi_spec["source_version"])
                metrics = None if reference else found.get(allele)
                status = "reference_label_no_substitution" if reference else score_status(metrics, absent_status)
                vw.writerow(dict(record, avi_raw="NA" if metrics is None or metrics["raw_score"] is None else metrics["raw_score"],
                                 avi_phred="NA" if metrics is None or metrics["PHRED"] is None else metrics["PHRED"], status=status))
                counts["variant_rows"] += 1
                if not reference:
                    counts["avi_status_" + status] += 1
                    available += status in {"available", "available_partial"}
                    complete += status == "available"
                    counts["theoretical_alt_rows"] += 1
                    counts["alts_listed_in_targets"] += allele in target["alts"]
                for rid, reader in readers.items():
                    spec = reader.spec
                    entries, missing_status, _ = queries[rid]
                    entry = None if reference else entries.get(allele)
                    for metric in spec["metrics"]:
                        value = None if entry is None else entry[metric]
                        metric_status = "reference_label_no_substitution" if reference else missing_status if entry is None else "score_missing" if value is None else "available"
                        mw.writerow(dict(record, resource_id=rid, source_version=spec["source_version"], metric=metric,
                                         value="NA" if value is None else value, status=metric_status,
                                         metric_semantics=metric_semantics(spec, metric),
                                         training_label_caveat=spec["training_label_caveat"], interpretation_limit=spec["interpretation_limit"]))
                        counts["metric_rows"] += 1
            counts["avi_sites_with_any_score"] += available > 0
            counts["avi_sites_with_all_three_complete"] += complete == 3
            aw.writerow(dict(common, avi_available_substitutions=available, avi_complete_substitutions=complete,
                             avi_unavailable_substitutions=3 - available, reference_check=ref_check))
    with (out / "resources.tsv").open("w", encoding="utf-8", newline="") as handle:
        rw = writer(handle, RESOURCE_FIELDS)
        for spec in specs:
            rendered = dict(spec)
            for key in ("index_check", "header_check", "compression_check", "physical_files"):
                rendered[key] = json.dumps(spec.get(key, {}), sort_keys=True, separators=(",", ":"))
            rw.writerow(rendered)
    usable = counts["avi_status_available"] + counts["avi_status_available_partial"]
    deferred = any(s["enabled"] and (s["status"] != "ready" or s.get("has_exclusions")) for s in specs)
    status = "not_available" if not usable else "complete" if counts["avi_status_available"] == 3 * len(request["targets"]) and not deferred else "complete_with_exclusions"
    return dict(status=status, schema=VERSION, semantics_version=SEMANTICS_VERSION, signature=request["signature"], counts=dict(counts),
                target_selection=request["selection"], parameters=request["parameters"], source_identity=request["source_identity"],
                resources=specs, query_audits={rid: dict(reader.audit) for rid, reader in readers.items()},
                physical_query_audits=pool.audits(), storage_version=STORAGE_VERSION,
                limitations=LIMITS, genome_build="GRCh38", coordinate_system="1-based",
                observed_from_targets_semantics="ALT membership in input alts; not an inferred population observation; REF is NA",
                source_version_semantics="user declaration, not independently verified official release",
                memory_policy="bounded target list plus at most three substitutions per resource and requested position in one bounded Tabix window; a shared VCF batch is fetched/decoded once for its logical resources; unrequested intervening rows discarded; no full atlas scan")


def run(args):
    declared_build = genome_build(args.build)
    if declared_build != "GRCh38":
        raise ValueError("AlphaGenome Atlas SNV resources require GRCh38; harmonize targets with verified liftover/REF checks first")
    if args.max_targets < 0:
        raise ValueError("--max-targets must be >= 0; 0 explicitly removes the candidate cap")
    if args.query_gap_bp < 0 or not 1 <= args.query_span_bp <= 1048576:
        raise ValueError("--query-gap-bp must be >= 0 and --query-span-bp must be within 1..1048576")
    if args.preflight_sites < 1 or args.preflight_sites > 1000:
        raise ValueError("--preflight-sites must be within 1..1000")
    args.resources = str(Path(args.resources or default_resources(args.dir_data)).expanduser().resolve())
    targets_path, resources_path = Path(args.targets).expanduser().resolve(), Path(args.resources)
    critical_sources = {"targets": identity(targets_path, full=True), "resources_manifest": identity(resources_path, full=True)}
    selected_chromosomes, selected_region = chromosome_filter(args.chromosomes), region_filter(args.region)
    targets, selection = load_targets(targets_path, declared_build, args.max_targets, selected_chromosomes, selected_region)
    specs = load_resources(args, targets)
    parameters = dict(build=declared_build, chromosomes=sorted(selected_chromosomes, key=chromosome_key) if selected_chromosomes is not None else "all",
                      region=list(selected_region) if selected_region is not None else None, max_targets=args.max_targets,
                      theoretical_alts="all_three", reference_score="NA", query_gap_bp=args.query_gap_bp,
                      query_span_bp=args.query_span_bp, semantics_version=SEMANTICS_VERSION, storage_version=STORAGE_VERSION)
    if critical_sources["targets"] != identity(targets_path, full=True) or critical_sources["resources_manifest"] != identity(resources_path, full=True):
        raise ValueError("Target/schema inputs changed while preparing the Atlas request")
    sources = dict(critical_sources, module=identity(__file__, full=True),
                   analysis_integrity=identity(Path(__file__).with_name("analysis_integrity.py"), full=True))
    for path in sorted({path for spec in specs for path in resource_input_paths(spec)}):
        sources[f"atlas_file:{path}"] = identity(path)
    # Include resolved schemas and paths, not only config identity: CLI overrides matter.
    signature_payload = dict(schema=VERSION, parameters=parameters, sources=sources, resources=specs)
    signature = hashlib.sha256(json.dumps(signature_payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    request = dict(targets=targets, resources=specs, selection=selection, parameters=parameters, source_identity=sources, signature=signature)
    check_output_paths(args.out, [InputPath(targets_path), InputPath(resources_path), *protected_inputs(request)])
    print(f"[build-check] atlas_targets={declared_build} coordinates=1-based selected_sites={len(targets)} theoretical_ALT_per_site=3 REF_score=NA", flush=True)
    if args.preflight:
        return preflight(request, args.preflight_sites)
    if not args.no_reuse and (Path(args.out) / "manifest.json").is_file():
        try:
            previous = validate_manifest(args.out)
        except (ValueError, OSError, KeyError) as exc:
            print(f"[reuse] atlas=refresh reason={exc}", flush=True)
        else:
            if previous.get("schema") == VERSION and previous.get("signature") == signature:
                print(f"[reuse] atlas=verified source_signature={signature} output={Path(args.out).resolve()}", flush=True)
                return dict(previous, reused=True)
            print("[reuse] atlas=refresh reason=source_or_parameters_changed", flush=True)
    result = publish(args.out, targets_path, resources_path, request)
    return dict(result, reused=False)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--targets", required=True, help="Bounded TSV: chrom,pos,ref,alts; optional site_id and build/genome_build/assembly/target_build; alts comma separated")
    p.add_argument("--build", "--target-build", dest="build", required=True, choices=("38", "GRCh38", "hg38", "37", "GRCh37", "hg19"), help="Explicit target assembly; only GRCh38 proceeds, no implicit liftover")
    p.add_argument("--avi-file", help="Local official six-column BGZF TSV or split2vcf .vcf.bgz; override AVI only, supports {chrom} for merged VCFs")
    p.add_argument("--avi-index", help="Explicit index paired with --avi-file (.tbi or .csi); required if both siblings exist")
    p.add_argument("--resources", help="Local resource schema TSV; auto prefers indexed chr*.vcf.bgz under alphagenome/38, then the legacy TSV descriptor")
    p.add_argument("--dir-data", default=os.environ.get("DIR_DATA", configured_path('DIR_BIG')), help="Expand {dir_data}; defaults to the shared local/HPC data root")
    p.add_argument("--source-version", help="User-declared AVI release identifier from the download portal or receipt; never guessed from filename")
    p.add_argument("--chromosomes", default="all", help="Filter targets, e.g. 1-22,X,Y or chr1; does not scan chromosomes")
    p.add_argument("--region", help="Filter targets with 1-based closed chr:start-end; does not scan the region")
    p.add_argument("--max-targets", type=int, default=5000, help="Selected target cap; 0 explicitly removes it. Use independent target batches on HPC")
    p.add_argument("--query-gap-bp", type=int, default=64, help="Coalesce nearby target positions only up to this inter-target gap; 0 uses exact point queries")
    p.add_argument("--query-span-bp", type=int, default=16384, help="Maximum span of one Tabix fetch (1..1048576); unrequested intervening rows are discarded")
    p.add_argument("--require-all-resources", action="store_true", help="Treat every enabled optional resource as required, failing on a missing file/index")
    p.add_argument("--preflight", action="store_true", help="Read-only BGZF/header/TBI-or-CSI checks and sampled exact REF lookups; creates no analysis files")
    p.add_argument("--preflight-sites", type=int, default=16, help="Maximum deterministic, spread-out target sites checked by --preflight (1..1000)")
    p.add_argument("--no-reuse", action="store_true", help="Force a fresh atomic publication instead of reusing matching source/output identities")
    p.add_argument("--out", required=True)
    return p


def main(argv=None):
    try:
        result = run(parser().parse_args(argv))
        print(json.dumps(dict(status=result["status"], counts=result["counts"], reused=result["reused"]), sort_keys=True))
        return 0 if result["status"] != "not_available" else 2
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        raise SystemExit(f"atlas: {exc}")


# Atlas Allofus
# ======================================================================
_description_atlas_allofus = 'Join completed All of Us AF-only queries to exact local Atlas substitutions.\n\nInputs must use the same target-file content and GRCh38 coordinates. Preserve\npopulation frequencies and availability separately from AVI scores; do not\nreconstruct REF frequencies, counts, denominators or functional distances.\n'


from allofus_shards import VERSION as ALLOFUS_VERSION, AF_COLUMNS, TARGET_FIELDS
ATLAS_VERSION = VERSION

ALLOFUS_JOIN_VERSION = "atlas-allofus-join-1"
ALLOFUS_FIELDS = ["site_id", "chrom", "pos", "build", "ref", "alt", "allele_role", "target_alts", "observed_from_targets",
          "avi_raw", "avi_phred", "atlas_status", "atlas_resource_id", "atlas_source_version",
          "allofus_availability", "allofus_match_status", "allofus_requested_alt", "allofus_observed_in_target",
          "allofus_additional_reported_alt", "allofus_source_ref", "allofus_evidence", *AF_COLUMNS, "gvs_max_af"]
ALLOFUS_LIMITS = [
    "This is an exact chrom/pos/REF/ALT join of completed derivatives from the same target-file content. It is not an allele discovery scan or a source-cohort independence test.",
    "Eight All of Us frequency columns and gvs_max_af retain their source meanings. gvs_max_af is a cross-population maximum indicator, not a ninth population. No denominators, counts, REF frequencies or joint allele probabilities are reconstructed.",
    "Unreported theoretical ALTs have NA frequencies, not zero. A reported AF=0 is preserved and does not demonstrate AC=0. Missing metrics, missing rows and uncallable alleles are not extinction.",
    "observed_from_targets, allofus_requested_alt and allofus_observed_in_target describe target-list membership. allofus_additional_reported_alt describes an additional summary row, not independently established genotype-level tri-allelism.",
    "The structural REF row has NA AF and NA AVI scores. REF is not assumed ancestral or most frequent. AVI raw/PHRED are reference-relative predictions; no ALT-to-ALT distance, pleiotropy count or selection statistic is inferred.",
    "Frequency provenance and prediction provenance remain separate. This join does not combine WGS counts with gnomAD denominators, ancient observations or imputed genotypes.",
]


def allofus_is_missing(value):
    return str(value).strip().upper() in MISSING


def allofus_number(value, label, frequency=False, nonnegative=False):
    if allofus_is_missing(value):
        return "NA"
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a numeric value or NA") from exc
    if not math.isfinite(parsed) or (frequency and not 0 <= parsed <= 1) or (nonnegative and parsed < 0):
        raise ValueError(f"Invalid {label}: {value}")
    return value


def allofus_bit(value, label):
    if value not in {"0", "1"}:
        raise ValueError(f"{label} must be an explicit 0/1")
    return int(value)


def allofus_row_key(row):
    if genome_build(row["build"]) != "GRCh38":
        raise ValueError("All join rows must declare GRCh38")
    if not re.fullmatch(r"[0-9]+", row["pos"]) or int(row["pos"]) < 1:
        raise ValueError("All join rows need positive 1-based positions")
    if len(row["ref"]) != 1 or row["ref"] not in BASES or not row["site_id"] or allofus_is_missing(row["site_id"]):
        raise ValueError("Join rows require a canonical REF and a site_id")
    return chromosome(row["chrom"]), int(row["pos"])


def allofus_require_inventory(manifest, filenames):
    if not set(filenames) <= set(manifest["output_inventory"]):
        raise ValueError("Required join tables are absent from the validated source output inventory")


def allofus_validate_sources(allofus_dir, atlas_dir):
    aou = validate_manifest(allofus_dir, require_evidence=True)
    atlas = validate_manifest(atlas_dir)
    if (aou.get("schema") != ALLOFUS_VERSION or aou.get("status") != "complete_descriptive_unverified"
            or aou.get("evidence") != "frequency_only" or "targets" not in aou or "target_count" not in aou):
        raise ValueError("--allofus-dir must be a completed allofus_shards.query AF-only result")
    if atlas.get("schema") != ATLAS_VERSION or atlas.get("status") not in {"complete", "complete_with_exclusions", "not_available"}:
        raise ValueError("--atlas-dir must be a completed atlas_scores result, including an explicit no-score audit if unavailable")
    if genome_build(aou.get("build")) != "GRCh38" or genome_build(atlas.get("genome_build")) != "GRCh38" or atlas.get("coordinate_system") != "1-based":
        raise ValueError("Both completed sources must use GRCh38, with 1-based Atlas coordinates")
    allofus_require_inventory(aou, ["target_alleles.tsv.gz", "target_sites.tsv"])
    allofus_require_inventory(atlas, ["atlas_variants.tsv.gz", "target_availability.tsv"])
    target_a, target_b = aou["targets"], atlas.get("source_identity", {}).get("targets", {})
    if not re.fullmatch(r"[a-fA-F0-9]{64}", str(target_a.get("sha256", ""))) or not re.fullmatch(r"[a-fA-F0-9]{64}", str(target_b.get("sha256", ""))):
        raise ValueError("Both source manifests must retain full target-file SHA256 identities")
    if any(target_a.get(k) != target_b.get(k) for k in ("sha256", "size")) or target_a.get("size", 0) <= 0:
        raise ValueError("Target-file content identities differ between All of Us and Atlas; do not join unrelated target requests")
    return aou, atlas


def allofus_read_atlas(directory, manifest):
    sites, site_ids = {}, {}
    avi_resources = [r for r in manifest.get("resources", []) if r.get("kind") == "avi" and r.get("status") == "ready"]
    if len(avi_resources) != 1:
        raise ValueError("Atlas manifest must identify its one active AVI resource")
    source = avi_resources[0]
    for _, row in table(Path(directory) / "atlas_variants.tsv.gz", VARIANT_FIELDS):
        key = allofus_row_key(row)
        sid, ref, allele = row["site_id"], row["ref"], row["allele"]
        if row["resource_id"] != source["resource_id"] or row["source_version"] != source["source_version"]:
            raise ValueError("Atlas row provenance conflicts with its resource manifest")
        alts = row["target_alts"].split(",")
        if not alts or len(alts) != len(set(alts)) or any(len(v) != 1 or v not in BASES or v == ref for v in alts):
            raise ValueError("Invalid target ALT membership in Atlas rows")
        if sid in site_ids and site_ids[sid] != key:
            raise ValueError("Atlas site_id refers to multiple coordinates")
        site_ids[sid] = key
        site = sites.setdefault(key, dict(site_id=sid, ref=ref, alts=set(alts), rows={}))
        if site["site_id"] != sid or site["ref"] != ref or site["alts"] != set(alts):
            raise ValueError("Atlas rows have inconsistent site identity, REF or target ALTs")
        if len(allele) != 1 or allele not in BASES or allele in site["rows"]:
            raise ValueError("Atlas requires exactly one row for each of the four bases")
        reference = allele == ref
        if row["allele_role"] != ("REF" if reference else "ALT"):
            raise ValueError("Atlas REF/ALT role does not match its exact nucleotide")
        if reference:
            if not all(allofus_is_missing(row[k]) for k in ("avi_raw", "avi_phred", "observed_from_targets")) or row["status"] != "reference_label_no_substitution":
                raise ValueError("Atlas REF must retain NA metrics and target-observation flag")
        elif allofus_bit(row["observed_from_targets"], "observed_from_targets") != int(allele in site["alts"]):
            raise ValueError("Atlas target membership flags conflict with target_alts")
        row["avi_raw"] = allofus_number(row["avi_raw"], "avi_raw")
        row["avi_phred"] = allofus_number(row["avi_phred"], "avi_phred", nonnegative=True)
        if not reference:
            numeric = sum(row[k] != "NA" for k in ("avi_raw", "avi_phred"))
            expected = {"available"} if numeric == 2 else {"available_partial"} if numeric else {"score_missing", "not_scored", "not_indexed_contig"}
            if row["status"] not in expected:
                raise ValueError("Atlas availability status conflicts with stored scores")
        site["rows"][allele] = row
    if not sites or len(sites) != manifest.get("counts", {}).get("target_sites") or any(set(s["rows"]) != set(BASES) for s in sites.values()):
        raise ValueError("Atlas is missing target sites or required theoretical-base rows")
    return sites


def allofus_read_allofus(directory, manifest, sites):
    rows, coordinates, requested_seen = {}, set(), {}
    for _, row in table(Path(directory) / "target_alleles.tsv.gz", TARGET_FIELDS):
        key = allofus_row_key(row)
        if key not in sites:
            raise ValueError("All of Us and Atlas target coordinate sets differ; region-filtered requests need matching target subsets")
        site = sites[key]
        if row["site_id"] != site["site_id"] or row["ref"] != site["ref"]:
            raise ValueError("All of Us / Atlas site_id or REF mismatch at the exact coordinate")
        alt = row["alt"]
        if len(alt) != 1 or alt not in BASES or alt == row["ref"] or (key, alt) in rows:
            raise ValueError("All of Us query requires unique canonical nonreference substitutions")
        if (not allofus_is_missing(row["source_ref"]) and row["source_ref"] != row["ref"]) or row["match_status"] == "reference_mismatch":
            raise ValueError("All of Us source REF conflicts with the harmonized target REF")
        if row["match_status"] not in {"matched_frequency_only", "matched_metrics_missing", "outside_prepared_windows", "no_reported_snv_rows_in_window", "not_in_summary"}:
            raise ValueError("Unrecognized All of Us query availability status")
        requested = allofus_bit(row["requested_alt"], "requested_alt")
        observed = allofus_bit(row["observed_in_target"], "observed_in_target")
        additional = allofus_bit(row["additional_reported_alt"], "additional_reported_alt")
        if requested != observed or requested != int(alt in site["alts"]) or requested and additional:
            raise ValueError("All of Us / Atlas target-list membership flags disagree")
        matched = row["match_status"] in {"matched_frequency_only", "matched_metrics_missing"}
        if additional != int(not requested and matched):
            raise ValueError("All of Us additional ALT flag conflicts with its availability status")
        for col in [*AF_COLUMNS, "gvs_max_af"]:
            row[col] = allofus_number(row[col], col, frequency=True)
        known_af = any(row[col] != "NA" for col in AF_COLUMNS)
        if (matched and allofus_is_missing(row["source_ref"])) or (not matched and any(row[col] != "NA" for col in [*AF_COLUMNS, "gvs_max_af"])):
            raise ValueError("All of Us AF values or source REF conflict with row availability")
        if (row["match_status"] == "matched_frequency_only" and not known_af) or (row["match_status"] == "matched_metrics_missing" and known_af):
            raise ValueError("All of Us metric availability is inconsistent")
        expected_evidence = "frequency_only" if known_af else "metrics_missing" if matched else "unavailable"
        if row["evidence"] != expected_evidence:
            raise ValueError("All of Us row evidence differs from its AF-only contract")
        rows[(key, alt)] = row
        coordinates.add(key)
        if requested:
            requested_seen.setdefault(key, set()).add(alt)
    if coordinates != set(sites) or len(coordinates) != manifest["target_count"]:
        raise ValueError("All of Us and Atlas completed target-site sets differ")
    if any(requested_seen.get(key, set()) != site["alts"] for key, site in sites.items()):
        raise ValueError("All of Us query omitted a target ALT instead of preserving an unavailable row")
    return rows


def allofus_protected(raw):
    return [InputPath(Path(raw[arg]) / filename) for arg, filenames in (
        ("allofus_dir", ["manifest.json", "target_alleles.tsv.gz", "target_sites.tsv"]),
        ("atlas_dir", ["manifest.json", "atlas_variants.tsv.gz", "target_availability.tsv"])) for filename in filenames]


@transactional("manifest.json", input_dirs=("allofus_dir", "atlas_dir"), input_prefixes=(), resolve_inputs=allofus_protected)
def allofus_run(allofus_dir, atlas_dir, out):
    aou, atlas = allofus_validate_sources(allofus_dir, atlas_dir)
    sites = allofus_read_atlas(atlas_dir, atlas)
    aou_rows = allofus_read_allofus(allofus_dir, aou, sites)
    counts = Counter(target_sites=len(sites))
    with gzip.open(Path(out) / "allofus_atlas.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        out_rows = writer(handle, ALLOFUS_FIELDS)
        for key in sorted(sites, key=lambda v: (chromosome_key(v[0]), v[1])):
            site = sites[key]
            for allele in [site["ref"], *[v for v in BASES if v != site["ref"]]]:
                ar = site["rows"][allele]
                aour = aou_rows.get((key, allele))
                reference = allele == site["ref"]
                record = dict(site_id=site["site_id"], chrom=key[0], pos=key[1], build="GRCh38", ref=site["ref"], alt=allele,
                              allele_role=ar["allele_role"], target_alts=ar["target_alts"], observed_from_targets=ar["observed_from_targets"],
                              avi_raw=ar["avi_raw"], avi_phred=ar["avi_phred"], atlas_status=ar["status"],
                              atlas_resource_id=ar["resource_id"], atlas_source_version=ar["source_version"])
                for source in ("match_status", "requested_alt", "observed_in_target", "additional_reported_alt", "source_ref", "evidence"):
                    record["allofus_" + source] = aour[source] if aour else "NA"
                for col in [*AF_COLUMNS, "gvs_max_af"]:
                    record[col] = aour[col] if aour else "NA"
                availability = "reference_AF_unavailable" if reference else "no_exact_alt_row_in_query" if aour is None else aour["match_status"]
                record["allofus_availability"] = availability
                out_rows.writerow(record)
                counts["rows"] += 1
                counts["allofus_status_" + availability] += 1
                if not reference:
                    af_available = any(record[col] != "NA" for col in AF_COLUMNS)
                    avi_available = ar["status"] in {"available", "available_partial"}
                    counts["substitutions_with_AF"] += af_available
                    counts["substitutions_with_AVI"] += avi_available
                    counts["substitutions_with_AF_and_AVI"] += af_available and avi_available
                    counts["substitutions_with_complete_metrics"] += ar["status"] == "available" and all(record[col] != "NA" for col in [*AF_COLUMNS, "gvs_max_af"])
    # Revalidate complete inventories, not only the subset of files read in this join.
    allofus_validate_sources(allofus_dir, atlas_dir)
    print(f"[build-check] atlas_allofus=GRCh38 same_target_SHA256=checked sites={len(sites)} exact_site_REF_ALT=checked", flush=True)
    print(f"[dataset] allofus=frequency_only populations=8 max_AF=separate AVI=reference_relative jointly_available_substitutions={counts['substitutions_with_AF_and_AVI']}", flush=True)
    return dict(status="complete" if counts["substitutions_with_complete_metrics"] == 3 * len(sites) else "complete_with_exclusions",
                schema=ALLOFUS_JOIN_VERSION, build="GRCh38", coordinate_system="1-based", counts=dict(counts), limits=ALLOFUS_LIMITS,
                target_identity=aou["targets"], target_identity_check="identical full SHA256 and size; exact site/REF/ALT and membership cross-check",
                source_manifests=dict(allofus=identity(Path(allofus_dir) / "manifest.json", full=True),
                                      atlas=identity(Path(atlas_dir) / "manifest.json", full=True)),
                allofus_provenance=dict(source_identity=aou.get("source_identity"), plan_id=aou.get("plan_id"), evidence="frequency_only"),
                atlas_provenance=dict(source_identity=atlas.get("source_identity"), resources=atlas.get("resources")),
                denominator_policy="not available in AF-only source; never reconstructed, pooled or transferred")


def allofus_parser():
    p = argparse.ArgumentParser(description=_description_atlas_allofus)
    p.add_argument("--allofus-dir", required=True, help="Completed allofus_shards query directory")
    p.add_argument("--atlas-dir", required=True, help="Completed atlas_scores directory for exactly the same target-file content")
    p.add_argument("--out", required=True)
    return p


def allofus_main(argv=None):
    args = allofus_parser().parse_args(argv)
    try:
        result = allofus_run(args.allofus_dir, args.atlas_dir, args.out)
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        raise SystemExit(f"atlas-allofus: {exc}")
    print(json.dumps(dict(status=result["status"], counts=result["counts"]), sort_keys=True))
    return 0


def cli(argv=None):
    """Dispatch auxiliary commands owned by this module."""
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    routes = {'allofus': allofus_main}
    entry = routes[argv.pop(0)] if argv and argv[0] in routes else main
    return entry(argv)



if __name__ == "__main__":
    raise SystemExit(cli())

