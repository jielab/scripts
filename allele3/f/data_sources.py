#!/usr/bin/env python3
"""Register existing sources and audit bounded coordinate/build evidence.

Reads headers and a bounded SNV sample, never scans an entire cohort, creates an
index, lifts coordinates, or rewrites inputs. Dataset/path names are declarations,
not evidence of an assembly or of individual genotypes. Matching a small sample
does not validate every locus. The supplied ANNOVAR All of Us header identifies
AF-only summaries; unknown All of Us schemas still require explicit mappings.

Public API: load_config(config) -> list[dict]; audit_entries(entries, out=None,
fasta37=None, fasta38=None, max_records=5000, max_sites=200, chromosomes=None).
Missing sources are deferred individually. The CLI exits 2 on build/reference
conflicts; unrelated missing datasets do not make the audit fail.
"""
from __future__ import annotations

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path

import argparse
from collections import Counter
from contextlib import contextmanager
import csv
import gzip
import io
import json
import os
from pathlib import Path
import re
import string
import sys

from analysis_integrity import InputPath, check_output_paths, transactional

SCHEMA = "allele3-data-sources-1"
TYPES = {"vcf", "vcf_directory", "pgen", "pvar", "scan", "gu", "allofus_text",
         "ukb_mfi", "annotation", "metadata"}
REQUIRED = ("id", "type", "path", "build", "evidence", "role")
AUTOSOMES = tuple(str(c) for c in range(1, 23))


def default_fasta(build):
    """Honor shared shell defaults, while keeping explicit CLI paths untouched."""
    key = "ALLELE_FASTA37" if build == "37" else "ALLELE_FASTA"
    return os.environ.get(key) or str(Path(os.environ.get("DIR_DATA") or configured_path('DIR_BIG'))/
                                      "gen"/"1hgp"/build/f"GRCH{build}.fasta.gz")


DEFAULT_FASTA37 = default_fasta("37")
DEFAULT_FASTA38 = default_fasta("38")
DEFAULT_CONFIG = Path(__file__).resolve().parents[1]/"config/data_sources.tsv"
MAX_LINE_BYTES = 16 * 1024 * 1024
ROOT_ENV = {"dir_data": ("DIR_DATA", "DIR_BIG"),
            "dir_work": ("DIR_WORK", "DIR_DAT"),
            "analysis_root": ("ANALYSIS_ROOT", "ANALYSIS_ROOT")}
ALLOFUS_COLUMNS = ["Chr", "Start", "End", "Ref", "Alt", "gvs_all_af", "gvs_max_af",
                   "gvs_afr_af", "gvs_amr_af", "gvs_eas_af", "gvs_eur_af",
                   "gvs_mid_af", "gvs_oth_af", "gvs_sas_af"]
# Primary-assembly nuclear chromosome lengths. Mitochondria are deliberately
# excluded: hg19 and later GRCh37 analysis references differ there.
LENGTH_SOURCES = {
    "37": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",
    "38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
}
LENGTHS = {
    "37": dict(zip((*AUTOSOMES, "X", "Y"),
        (249250621,243199373,198022430,191154276,180915260,171115067,159138663,
         146364022,141213431,135534747,135006516,133851895,115169878,107349540,
         102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566))),
    "38": dict(zip((*AUTOSOMES, "X", "Y"),
        (248956422,242193529,198295559,190214555,181538259,170805979,159345973,
         145138636,138394717,133797422,135086622,133275309,114364328,107043718,
         101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415))),
}


class BuildEvidenceConflict(ValueError):
    """Contradictory build metadata is a blocking result, not a missing input."""


def canonical_chrom(value):
    value = str(value).removeprefix("chr")
    return value if value in {*AUTOSOMES, "X", "Y"} else None


def parse_chromosomes(value=None):
    if value is None or value == "autosomes":
        return AUTOSOMES
    items = value.replace(",", " ").split() if isinstance(value, str) else list(value)
    result = tuple(canonical_chrom(c) for c in items)
    if not result or None in result or len(result) != len(set(result)):
        raise ValueError("Chromosomes must be distinct canonical 1..22, X or Y; default is autosomes")
    return result


def normalize_build(value):
    value = str(value).strip()
    aliases = {"37":"37", "grch37":"37", "hg19":"37", "b37":"37",
               "38":"38", "grch38":"38", "hg38":"38", "b38":"38",
               "auto":"metadata", "metadata":"metadata", "unknown":"unknown",
               "na":"not_applicable", "n/a":"not_applicable", "not_applicable":"not_applicable"}
    if value.lower() not in aliases:
        raise ValueError(f"Unsupported build declaration: {value}")
    return aliases[value.lower()]


def check_entries(entries):
    result, ids = [], set()
    for raw in entries:
        row = {k: str(v).strip() if v is not None else "" for k, v in raw.items()}
        missing = [key for key in REQUIRED if not row.get(key)]
        if missing:
            raise ValueError("Missing source fields: " + ", ".join(missing))
        if row["id"] in ids:
            raise ValueError("Duplicate source id: " + row["id"])
        ids.add(row["id"])
        if row["type"] not in TYPES:
            raise ValueError("Unsupported source type: " + row["type"])
        row["build"] = normalize_build(row["build"])
        for _, field, spec, conversion in string.Formatter().parse(row["path"]):
            if field is not None and (field != "chrom" or spec or conversion):
                raise ValueError("Only a literal {chrom} path placeholder is supported")
        result.append(row)
    if not result:
        raise ValueError("No data sources were selected")
    return result


def expand_registry_path(value):
    """Expand only explicit registry root tokens, leaving {chrom} for units().

    No environment-variable interpolation, format expressions, path-prefix
    rewriting or manifest rewriting occurs. A concrete path stays concrete.
    DIR_* normally comes from 0f/0os.sh; direct Python calls use local defaults.
    """
    path = str(value)
    parts = list(string.Formatter().parse(path))
    for _, field, spec, conversion in parts:
        if field is not None and (field not in {*ROOT_ENV, "chrom"} or spec or conversion):
            raise ValueError("Registry path permits only {dir_data}, {dir_work}, {analysis_root} and literal {chrom}")
    # Preserve old concrete paths exactly; only root-bearing templates expand.
    if not any(field in ROOT_ENV for _, field, _, _ in parts):
        return path
    roots = {key: os.environ.get(env) or configured_path(default) for key, (env, default) in ROOT_ENV.items()}
    if any("{" in value or "}" in value for value in roots.values()):
        raise ValueError("Registry root environment values must be concrete paths, without brace placeholders")
    return "".join(literal + ("{chrom}" if field == "chrom" else roots[field] if field else "")
                   for literal, field, _, _ in parts)


def load_config(config):
    """Read registry root templates; callers may override actual paths afterward."""
    with Path(config).open(encoding="utf-8-sig") as handle:
        reader = csv.DictReader((line for line in handle if line.strip() and not line.startswith("#")), delimiter="\t")
        if not set(REQUIRED) <= set(reader.fieldnames or []):
            raise ValueError("Registry needs columns: " + ", ".join(REQUIRED))
        entries = list(reader)
        for row in entries:
            row["path"] = expand_registry_path(row["path"] or "")
        return check_entries(entries)


@contextmanager
def text_reader(path):
    path = Path(path)
    if str(path).endswith(".zst"):
        try:
            import zstandard
        except ImportError as exc:
            raise RuntimeError("requires_zstandard_reader") from exc
        with path.open("rb") as raw, zstandard.ZstdDecompressor().stream_reader(raw) as stream, io.TextIOWrapper(stream) as handle:
            yield handle
    elif str(path).endswith((".gz", ".bgz")):
        with gzip.open(path, "rt") as handle:
            yield handle
    else:
        with path.open() as handle:
            yield handle


def bounded_lines(handle, maximum):
    for _ in range(maximum):
        line = handle.readline(MAX_LINE_BYTES + 1)
        if not line:
            break
        if len(line) > MAX_LINE_BYTES:
            raise ValueError("Input line exceeds bounded audit size")
        yield line.rstrip("\r\n")


def build_mentions(value):
    value = str(value)
    if value.strip() in {"37", "38"}:
        return {value.strip()}
    found = set()
    for token in re.findall(r"(?i)(?<![a-z0-9])(?:GRCh(?:37|38)|hg(?:19|38)|b(?:37|38)|build[ _-]?(?:37|38))(?![a-z0-9])", value):
        found.add("37" if token.lower().endswith(("37", "19")) else "38")
    return found


def length_evidence(lengths):
    matches = {build: sorted(c for c, n in lengths.items() if LENGTHS[build].get(c) == n) for build in ("37", "38")}
    unknown = sorted(c for c, n in lengths.items() if all(LENGTHS[b].get(c) != n for b in ("37", "38")))
    return dict(matches=matches, unrecognized=unknown)


def sample_variant_metadata(path, kind, chromosomes, max_records, max_sites):
    """One bounded sequential read; no sample IDs or individual GTs are retained."""
    result = dict(header_builds=[], reference_headers=[], contig_lengths={}, samples=[],
                  sample_count=None, genotype_scope="not_assessed", info_fields=[], lines_read=0,
                  variants_examined=0, header_complete=False, provisional_ref_sites_skipped=0)
    builds, infos, columns, samples = set(), set(), None, []
    with text_reader(path) as handle:
        for line in bounded_lines(handle, max_records):
            result["lines_read"] += 1
            if line.startswith("##"):
                key, _, value = line[2:].partition("=")
                if key.lower() in {"reference", "assembly", "genome_build", "build"}:
                    builds.update(build_mentions(value))
                    result["reference_headers"].append(line[:1000])
                if key.lower() == "contig":
                    fields = dict(re.findall(r'(\w+)=([^,>]+)', value))
                    chrom = canonical_chrom(fields.get("ID", "").strip('"'))
                    length = fields.get("length", fields.get("Length", "")).strip('"')
                    if chrom in chromosomes and length:
                        n = int(length)
                        if n < 1:
                            raise ValueError("Invalid contig length")
                        if chrom in result["contig_lengths"] and result["contig_lengths"][chrom] != n:
                            raise BuildEvidenceConflict("Contradictory lengths for the same canonical contig")
                        result["contig_lengths"][chrom] = n
                    if "assembly" in fields:
                        builds.update(build_mentions(fields["assembly"]))
                if key == "INFO":
                    match = re.search(r"(?:^<|,)ID=([^,>]+)", value)
                    if match:
                        infos.add(match.group(1))
                continue
            if line.startswith("#CHROM"):
                columns = line.lstrip("#").split()
                if not {"CHROM", "POS", "REF", "ALT"} <= set(columns):
                    raise ValueError("Variant table lacks coordinate/REF/ALT fields")
                result["header_complete"] = True
                if kind == "vcf":
                    count = max(0, len(columns)-9)
                    result.update(sample_count=count, genotype_scope="individuals" if count else "sites_only")
                else:
                    result["genotype_scope"] = "variant_metadata_only"
                continue
            if not line or line.startswith("#"):
                continue
            if columns is None:
                raise ValueError("No #CHROM header before variant records")
            result["variants_examined"] += 1
            # Split only the fields before genotypes; 500k sample strings need not
            # become 500k Python objects for a coordinate audit.
            fields = line.split("\t", 9) if "\t" in line else line.split(maxsplit=9)
            values = {name: fields[index] for index, name in enumerate(columns[:9]) if index < len(fields)}
            chrom = canonical_chrom(values.get("CHROM", ""))
            if chrom not in chromosomes:
                continue
            pos = int(values["POS"])
            if pos < 1:
                raise ValueError("Variant POS must be positive and 1-based")
            ref = values.get("REF", "").upper()
            alts = values.get("ALT", "").upper().split(",")
            if ref not in {"A", "C", "G", "T"} or not alts or any(a not in {"A", "C", "G", "T"} for a in alts):
                continue
            # PLINK can emit INFO/PR in both PVAR and VCF. These are provisional
            # allele labels, not trustworthy REF tests (official format docs:
            # https://www.cog-genomics.org/plink/2.0/formats#pvar).
            if "PR" in values.get("INFO", "").split(";"):
                result["provisional_ref_sites_skipped"] += 1
                continue
            samples.append((chrom, pos, ref))
            if len(samples) >= max_sites:
                break
    result.update(header_builds=sorted(builds), info_fields=sorted(infos), samples=samples,
                  sampled_snvs=len(samples), bounded_read_limit_reached=result["lines_read"] == max_records)
    return result


class IndexedReference:
    """Only open pre-indexed FASTA; this class never calls faidx or builds indexes."""
    def __init__(self, path, label, max_records):
        self.reader, self.names, self.lengths = None, {}, {}
        self.info = dict(label=label, path=str(path) if path else None, status="not_supplied")
        if not path:
            return
        p = Path(path).expanduser()
        self.info["path"] = str(p)
        required = [p, Path(str(p)+".fai")]
        if str(p).endswith((".gz", ".bgz")):
            required.append(Path(str(p)+".gzi"))
        missing = [str(x) for x in required if not x.is_file()]
        if missing:
            self.info.update(status="deferred_missing_fasta_or_index", missing=missing)
            return
        try:
            count = 0
            with required[1].open() as handle:
                for line in bounded_lines(handle, max_records):
                    count += 1
                    fields = line.split("\t")
                    chrom = canonical_chrom(fields[0])
                    if chrom:
                        if chrom in self.names:
                            raise ValueError("Ambiguous canonical aliases in FASTA index")
                        self.names[chrom], self.lengths[chrom] = fields[0], int(fields[1])
            if count == max_records:
                self.info.update(status="deferred_fai_read_limit", index_lines_read=count)
                return
            evidence = length_evidence(self.lengths)
            other = "38" if label == "37" else "37"
            if evidence["matches"][other]:
                self.info.update(status="reference_label_conflict", length_evidence=evidence)
                return
            self.info["assembly_identity"] = "canonical_lengths_support_label" if evidence["matches"][label] else "unverified_custom_or_partial_reference"
            self.info["length_evidence"] = evidence
            import pysam
            self.reader = pysam.FastaFile(str(p))
            self.info["status"] = "ready"
        except ImportError:
            self.info["status"] = "deferred_requires_pysam"
        except (ValueError, OSError, IndexError) as exc:
            self.info.update(status="deferred_unreadable_reference", error=str(exc))

    def base(self, chrom, pos):
        if self.reader is None or chrom not in self.names:
            return None
        if not 1 <= pos <= self.lengths[chrom]:
            return "out_of_range"
        value = self.reader.fetch(self.names[chrom], pos-1, pos).upper()
        return value if value in {"A", "C", "G", "T"} else None

    def close(self):
        if self.reader is not None:
            self.reader.close()


def compare_refs(samples, references):
    counts = {b: Counter() for b in ("37", "38")}
    joint = Counter()
    examples = []
    for chrom, pos, ref in samples:
        bases = {b: references[b].base(chrom, pos) for b in references}
        for b, base in bases.items():
            counts[b]["unavailable" if base is None else "out_of_range" if base == "out_of_range" else "match" if base == ref else "mismatch"] += 1
        if all(bases[b] in {"A", "C", "G", "T"} for b in references):
            joint["comparable"] += 1
            m37, m38 = bases["37"] == ref, bases["38"] == ref
            joint["both_match" if m37 and m38 else "only_37" if m37 else "only_38" if m38 else "neither_match"] += 1
            if m37 != m38 and len(examples) < 5:
                examples.append(dict(chrom=chrom, pos=pos, REF=ref, fasta37=bases["37"], fasta38=bases["38"]))
    support = None
    informative = joint["only_37"] + joint["only_38"]
    if joint["comparable"] >= 20 and informative >= 5:
        for b in ("37", "38"):
            called = counts[b]["match"] + counts[b]["mismatch"]
            if called and counts[b]["match"]/called >= .95 and joint["only_"+b]/informative >= .95:
                support = b
    return dict(per_reference={k: dict(v) for k, v in counts.items()}, joint=dict(joint),
                discriminating_support=support, diagnostic_examples=examples,
                rule="At least 20 dual-reference comparable SNVs, 5 discriminating SNVs, >=95% REF match and >=95% directional support; limited prefix sample only")


def assess_build(declared, metadata, references):
    headers = set(metadata.get("header_builds", []))
    lengths = length_evidence(metadata.get("contig_lengths", {}))
    length_builds = {b for b, cs in lengths["matches"].items() if cs}
    seq = compare_refs(metadata.get("samples", []), references)
    reasons = []
    observed = headers | length_builds
    if len(headers) > 1:
        reasons.append("Header/provenance declares both builds")
    if len(length_builds) > 1:
        reasons.append("Canonical contig lengths mix GRCh37 and GRCh38")
    if len(observed) > 1:
        reasons.append("Header build and canonical contig lengths disagree")
    if declared in {"37", "38"} and observed - {declared}:
        reasons.append("Observed header/contig evidence conflicts with declared build " + declared)
    effective = declared if declared in {"37", "38"} else next(iter(observed)) if len(observed) == 1 else None
    if effective:
        out_of_range = sum(pos > LENGTHS[effective][chrom] for chrom, pos, _ in metadata.get("samples", []))
        if out_of_range:
            reasons.append(f"{out_of_range} sampled positions exceed canonical build {effective} chromosome lengths")
        if seq["discriminating_support"] and seq["discriminating_support"] != effective:
            reasons.append("REF/coordinate sample supports the opposite supplied FASTA; verify assembly and REF normalization")
        if references[effective].info["status"] == "reference_label_conflict":
            reasons.append("Supplied FASTA canonical lengths conflict with its requested build label")
    if declared == "not_applicable":
        status = "not_applicable"
    elif reasons:
        status = "conflict"
    elif seq["discriminating_support"]:
        status = "sample_supports_supplied_reference_" + seq["discriminating_support"]
    elif len(observed) == 1:
        status = "metadata_supports_" + next(iter(observed))
    elif seq["joint"].get("comparable", 0):
        status = "ambiguous_two_builds"
    elif declared in {"37", "38"}:
        status = "declared_only"
    else:
        status = "unverified_requires_metadata"
    # Even metadata-supported rows can have an ambiguous SNV sequence check.
    sequence_status = "discriminating_support" if seq["discriminating_support"] else "ambiguous_or_insufficient" if seq["joint"].get("comparable") else "unavailable"
    return dict(build_status=status, conflict=bool(reasons), conflict_reasons=reasons,
                observed_builds=sorted(observed), contig_evidence=lengths, ref_check=seq,
                sequence_check_status=sequence_status, evidence_scope="bounded metadata and first eligible SNVs; never whole-file validation")


def existing_pvar(prefix):
    variants = [Path(str(prefix)+s) for s in (".pvar", ".pvar.zst", ".pvar.gz")]
    existing = [p for p in variants if p.is_file()]
    if len(existing) > 1:
        raise ValueError("Multiple PVAR versions exist; supply type=pvar with an explicit file for metadata-only audit")
    return existing[0] if existing else variants[0]


def units(entry, chromosomes):
    raw, kind = str(Path(entry["path"]).expanduser()), entry["type"]
    if "{chrom}" in raw:
        return [(chrom, Path(raw.format(chrom=chrom))) for chrom in chromosomes]
    path = Path(raw)
    if kind == "vcf_directory":
        files = sorted(p for p in path.iterdir() if p.is_file() and p.name.endswith((".vcf", ".vcf.gz", ".vcf.bgz"))) if path.is_dir() else []
        return [(p.name, p) for p in files] if files else [("directory", path)]
    if kind == "gu" and path.is_dir():
        if (path/"run.meta.tsv").is_file():
            return [(path.name, path/"run.meta.tsv")]
        module = path if path.name in {"ibdmix", "ibdmix-GPU"} else path/entry.get("module", "ibdmix")
        return [(chrom, module/f"chr{chrom}"/"run.meta.tsv") for chrom in chromosomes]
    if kind == "gu" and path.name != "run.meta.tsv":
        module = path if path.name in {"ibdmix", "ibdmix-GPU"} else path/entry.get("module", "ibdmix")
        return [(chrom, module/f"chr{chrom}"/"run.meta.tsv") for chrom in chromosomes]
    return [("combined", path)]


def bounded_json(path, max_records):
    limit = min(MAX_LINE_BYTES, max_records*8192)
    with Path(path).open() as handle:
        text = handle.read(limit+1)
    if len(text) > limit or text.count("\n") >= max_records:
        raise RuntimeError("requires_larger_metadata_read_limit")
    return json.loads(text)


def inspect_metadata(path, kind, max_records):
    result = dict(header_builds=[], samples=[], contig_lengths={})
    if kind == "scan":
        manifest = path/"run.json" if path.is_dir() else path
        data = bounded_json(manifest, max_records)
        options = data.get("options", data.get("identity", {}).get("options", {}))
        builds = set()
        for record in (data, options):
            for field in ("build", "genome_build", "assembly"):
                if field in record:
                    builds.update(build_mentions(record[field]))
        inputs = data.get("inputs", data.get("identity", {}).get("inputs", []))
        result.update(header_builds=sorted(builds), run_status=data.get("status", "unknown"),
                      recorded_input_paths=[x.get("path") if isinstance(x, dict) else x for x in inputs][:max_records],
                      scan_identity_checked=False, note="Reads existing run metadata; original scan/store validator remains authoritative")
    elif kind == "gu":
        builds, fields = set(), {}
        with path.open() as handle:
            for line in bounded_lines(handle, max_records):
                pair = line.replace("\\t", "\t").split("\t", 1)
                if len(pair) != 2:
                    continue
                if pair[0] in {"genome_build", "build", "assembly"}:
                    builds.update(build_mentions(pair[1]))
                if pair[0] in {"genome_build", "coordinate_system", "pipeline_version", "refs", "lod_cut", "len_cut"}:
                    fields[pair[0]] = pair[1]
        result.update(header_builds=sorted(builds), gu_metadata=fields,
                      note="Build comes from actual GU metadata; this audit does not validate tracts, rosters, calibration or phase")
    else:
        result["note"] = "No nucleotide REF schema is assumed for this source"
        first, sampled_rows = None, []
        if path.suffix.lower() not in {".xlsx", ".xls", ".pdf"}:
            counts, first, lines_read = [], None, 0
            with text_reader(path) as handle:
                for line in bounded_lines(handle, min(max_records, 10)):
                    lines_read += 1
                    tokens = line.split("\t") if "\t" in line else line.split()
                    counts.append(len(tokens))
                    if first is None:
                        first = tokens
                    elif kind == "allofus_text":
                        sampled_rows.append(tokens)
            # Do not put individual records/IDs from an unknown table into logs.
            result["sniff"] = dict(lines_read=lines_read, column_counts=counts,
                first_row_field_types=["integer" if re.fullmatch(r"[+-]?\d+", x) else "numeric" if re.fullmatch(r"[+-]?\d*\.\d+(?:[eE][+-]?\d+)?", x) else "text" for x in (first or [])[:30]])
        if kind == "allofus_text":
            result.update(schema_status="requires_schema", index_status="unclassified_idx_not_tabix" if Path(str(path)+".idx").is_file() else "no_idx_seen",
                          note="Unknown All of Us schema: no AC/AN/AF, sample count, coordinate or callable denominator is inferred; .idx is not assumed to be tabix")
            header = [first[0].lstrip("#"), *first[1:]] if first else []
            if header == ALLOFUS_COLUMNS:
                result.update(schema_status="annovar-allofus-hg38", coordinate_system="1-based-inclusive",
                    populations=["ALL", "AFR", "AMR", "EAS", "EUR", "MID", "OTH", "SAS"],
                    additional_metrics=["gvs_max_af"], counts_available=False, genotypes_available=False,
                    note="Known AF-only header: preserve empty fields as NA; gvs_max_af is not a population. No AC/AN, reference frequency or callable denominator is inferred. The header does not independently establish the declared build; .idx is not assumed to be tabix.")
                for values in sampled_rows:
                    if len(values) != len(ALLOFUS_COLUMNS):
                        raise ValueError("All of Us row field count differs from supplied tab-delimited schema")
                    chrom = canonical_chrom(values[0])
                    start, end, ref, alt = values[1:5]
                    if chrom is not None and ref in {"A", "C", "G", "T"} and alt in {"A", "C", "G", "T"} and ref != alt:
                        pos = int(start)
                        if pos > 0 and pos == int(end):
                            result["samples"].append((chrom, pos, ref))
        elif kind == "ukb_mfi":
            result.update(schema_status="metadata_only", note="Imputation metadata is not individual GT, an observed third-base call, or a WGS callability denominator")
    return result


def audit_unit(entry, unit, path, references, chromosomes, max_records, max_sites):
    kind = entry["type"]
    row = dict(id=entry["id"], unit=unit, type=kind, path=str(path), declared_build=entry["build"],
               evidence=entry["evidence"], role=entry["role"], availability="available", status="audited",
               build_status="unverified", conflict=False, conflict_reasons=[], genotype_scope="not_applicable",
               sample_count=None, genotype_files_present=False, missing_files=[], details={})
    read_path, meta_kind = path, kind
    if kind == "scan":
        read_path = path/"run.json" if path.is_dir() or path.suffix != ".json" else path
    elif kind == "pgen":
        pvar = existing_pvar(path)
        required = [Path(str(path)+".pgen"), pvar, Path(str(path)+".psam")]
        row["missing_files"] = [str(p) for p in required if not p.is_file()]
        row["genotype_files_present"] = not row["missing_files"]
        row["details"]["fileset"] = {"pgen":str(required[0]), "pvar":str(pvar), "psam":str(required[2]),
                                    "acount_present":Path(str(path)+".acount").is_file(),
                                    "genotype_content_verified":False}
        if row["missing_files"]:
            row.update(availability="partial" if any(p.is_file() for p in required) else "missing", status="deferred_missing_components")
        read_path, meta_kind = pvar, "pvar"
    elif kind == "vcf_directory":
        meta_kind = "vcf"
    if not read_path.is_file():
        row.update(availability="missing" if row["availability"] == "available" else row["availability"],
                   status="deferred_missing", build_status="not_checked_missing")
        row["missing_files"] = sorted(set(row["missing_files"]+[str(read_path)]))
        return row
    if meta_kind in {"vcf", "pvar"}:
        metadata = sample_variant_metadata(read_path, meta_kind, chromosomes, max_records, max_sites)
        row["genotype_scope"] = metadata["genotype_scope"]
        row["sample_count"] = metadata["sample_count"]
        if kind == "pgen":
            row["genotype_scope"] = "individual_fileset_present_not_decoded" if row["genotype_files_present"] else "variant_metadata_only"
        if not metadata["header_complete"]:
            row["status"] = "deferred_header_limit"
        row["details"]["index_files_present"] = [str(p) for p in (Path(str(read_path)+".tbi"), Path(str(read_path)+".csi")) if p.is_file()]
        if meta_kind == "vcf":
            row["genotype_files_present"] = bool(metadata["sample_count"])
    else:
        metadata = inspect_metadata(read_path if kind != "scan" else path, kind, max_records)
        if kind == "allofus_text":
            row["status"] = "audited_af_only" if metadata.get("schema_status") == "annovar-allofus-hg38" else "requires_schema"
    row.update(assess_build(entry["build"], metadata, references))
    metadata.pop("samples", None)
    row["details"]["metadata"] = metadata
    if row["conflict"]:
        row["status"] = "conflict"
    if entry["build"] == "not_applicable":
        row["build_status"] = "not_applicable"
    return row


def source_plan(entry, rows):
    evidence, kind = entry["evidence"], entry["type"]
    if kind == "allofus_text":
        if evidence == "aggregate_frequencies" or any(r.get("details", {}).get("metadata", {}).get("schema_status") == "annovar-allofus-hg38" for r in rows):
            use = "AF-only cohort corroboration: reported positive ALT identities, population frequency contrasts and indexed candidate lookup after wide-cache preparation"
            limit = "No AC/AN or callability: never pool with count-supported AC>=5 gnomAD enrichment; missing entries remain unknown, gvs_max_af is not a population, cohort overlap is unverified"
        else:
            use = "Confirm the actual summary schema before selecting frequency or count analyses"
            limit = "Unknown current schema and denominator; do not infer GT, zero counts or index type"
    elif kind == "ukb_mfi":
        use = "Imputation quality and allele-representation covariates"
        limit = "No individual GT; no WGS or third-base observation inferred"
    elif kind in {"annotation"}:
        use = "Reference-coordinate overlaps and documented annotation partitions"
        limit = "Annotation coverage is not individual insertion burden or a mutation rate"
    elif kind in {"metadata"}:
        use = entry["role"]
        limit = entry.get("notes", "Metadata only; no allele measurements")
    elif kind == "gu":
        use = "Read actual caller/build metadata before per-person tract cooccurrence"
        limit = "Unphased IBDmix intervals cannot assign the third base to an archaic homolog"
    elif "imputed" in evidence:
        use = "Descriptive frequencies for explicitly encoded imputed alleles; fixed-panel sensitivity analysis"
        limit = "Neither .acount nor imputed PGEN is WGS; unencoded third bases remain unavailable"
    elif any(r.get("genotype_scope") == "sites_only" for r in rows):
        use = "Aggregate allele-count analysis if actual INFO fields and denominators support it"
        limit = "Observed VCF has no individual columns; individual analyses need a separate genotype source"
    elif evidence in {"aggregate_counts", "scan_counters"}:
        use = "Population aggregate multi-ALT occupancy and matched enrichment"
        limit = "No individual genotype, haplotype, mutation age or per-base callability is supplied by counts alone"
    else:
        use = "Individual allele measurement after genotype representation and callability checks"
        limit = "Fileset/header audit does not validate GT content, phase, imputation status or absent alleles"
    return dict(id=entry["id"], type=kind, evidence=evidence, role=entry["role"], proposed_use=use,
                constraint=limit, build_states=",".join(sorted({r["build_status"] for r in rows})),
                availability_counts=dict(Counter(r["availability"] for r in rows)),
                next_step="Resolve conflicting build/reference evidence before using this source" if any(r["conflict"] for r in rows)
                else "Provide missing native files or schema; other available sources remain independent" if any(r["status"].startswith("deferred") or r["status"] == "requires_schema" for r in rows)
                else "Use only this source's supported role; retain the original downstream strict checks")


def _audit(entries, fasta37, fasta38, max_records, max_sites, chromosomes):
    references = {"37":IndexedReference(fasta37, "37", max_records), "38":IndexedReference(fasta38, "38", max_records)}
    datasets, plans = [], []
    try:
        for entry in entries:
            entry_rows = []
            for unit, path in units(entry, chromosomes):
                try:
                    row = audit_unit(entry, unit, path, references, chromosomes, max_records, max_sites)
                except BuildEvidenceConflict as exc:
                    row = dict(id=entry["id"], unit=unit, type=entry["type"], path=str(path), declared_build=entry["build"],
                               evidence=entry["evidence"], role=entry["role"], availability="available", status="conflict",
                               build_status="conflict", conflict=True, conflict_reasons=[str(exc)])
                except (OSError, ValueError, RuntimeError, UnicodeError) as exc:
                    row = dict(id=entry["id"], unit=unit, type=entry["type"], path=str(path), declared_build=entry["build"],
                               evidence=entry["evidence"], role=entry["role"], availability="unreadable", status="deferred_read_error",
                               build_status="unverified_read_error", conflict=False, conflict_reasons=[], error=str(exc))
                entry_rows.append(row)
                print(f"[dataset] {row['id']} {unit}: {row['availability']}; status={row['status']}; {row.get('genotype_scope','unassessed')}; "
                      f"sample_count={row.get('sample_count') if row.get('sample_count') is not None else 'NA'}; {path}", flush=True)
                reasons = "; ".join(row.get("conflict_reasons", [])) or row.get("error", "limited audit; no whole-file validation")
                print(f"[build-check] {row['id']} {unit}: declared={entry['build']} status={row['build_status']} "
                      f"dual_ref={row.get('sequence_check_status','unavailable')}; {reasons}", flush=True)
            datasets.extend(entry_rows)
            plans.append(source_plan(entry, entry_rows))
        n_conflicts = sum(bool(row["conflict"]) for row in datasets)
        return dict(status="complete", schema=SCHEMA, has_conflicts=bool(n_conflicts), n_conflicts=n_conflicts,
                    chromosomes=list(chromosomes), datasets=datasets, source_plan=plans,
                    references={b:r.info for b,r in references.items()}, settings=dict(max_records=max_records,max_sites=max_sites),
                    length_sources=LENGTH_SOURCES,
                    limitations=["No path/name-based build inference; headers are declarations and contig lengths are metadata evidence.",
                                 "First eligible SNVs within a bounded prefix are not a random genome-wide validation sample.",
                                 "Ambiguous or unavailable dual-FASTA checks never establish a validated build.",
                                 "Only existing FASTA indexes are used; no liftover, downloads, genotype scan or input edits occur.",
                                 "Presence of a PGEN fileset or VCF sample columns does not validate genotype content or phase."])
    finally:
        for reference in references.values():
            reference.close()


def protected_inputs(raw):
    inputs = []
    for entry in raw["entries"]:
        for _, p in units(entry, raw["chromosomes"]):
            inputs.append(InputPath(p, "source"))
            if entry["type"] == "pgen":
                inputs.extend(InputPath(str(p)+s) for s in (".pgen", ".pvar", ".pvar.zst", ".pvar.gz", ".psam", ".acount"))
    for key in ("fasta37", "fasta38"):
        if raw.get(key):
            inputs.extend(InputPath(str(raw[key])+s) for s in ("", ".fai", ".gzi"))
    return inputs


def write_tsv(path, rows):
    columns = list(dict.fromkeys(k for row in rows for k in row)) or ["status"]
    with Path(path).open("w") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if v is None else json.dumps(v, sort_keys=True) if isinstance(v, (dict,list)) else v for k,v in row.items()})


@transactional("audit.json", input_prefixes=(), resolve_inputs=protected_inputs)
def _audit_to_disk(entries, out, fasta37, fasta38, max_records, max_sites, chromosomes):
    result = _audit(entries, fasta37, fasta38, max_records, max_sites, chromosomes)
    write_tsv(Path(out)/"audit.tsv", result["datasets"])
    write_tsv(Path(out)/"source_plan.tsv", result["source_plan"])
    return result


def audit_entries(entries, out=None, fasta37=None, fasta38=None, max_records=5000, max_sites=200, chromosomes=None):
    """Audit only selected sources; callers decide which workflow a conflict blocks."""
    entries = check_entries(entries)
    if not isinstance(max_records, int) or not isinstance(max_sites, int) or max_records < 1 or max_sites < 1:
        raise ValueError("Positive integer record and sample limits are required")
    chromosomes = parse_chromosomes(chromosomes)
    if out is not None:
        return _audit_to_disk(entries, out, fasta37, fasta38, max_records, max_sites, chromosomes)
    return _audit(entries, fasta37, fasta38, max_records, max_sites, chromosomes)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("audit",))
    p.add_argument("--config", default=str(DEFAULT_CONFIG))
    p.add_argument("--datasets", help="Comma-separated registered IDs; omitted audits every configured source")
    p.add_argument("--out", required=True)
    p.add_argument("--fasta37", default=default_fasta("37"))
    p.add_argument("--fasta38", default=default_fasta("38"))
    p.add_argument("--chromosomes", default="autosomes")
    p.add_argument("--max-records", type=int, default=5000, help="Maximum text lines per file, including headers")
    p.add_argument("--max-sites", type=int, default=200, help="Maximum first eligible SNVs per file for REF comparisons")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    entries = load_config(args.config)
    check_output_paths(args.out, [InputPath(args.config, "registry")])
    if args.datasets:
        selected = args.datasets.replace(",", " ").split()
        missing = set(selected)-{r["id"] for r in entries}
        if missing:
            raise ValueError("Unknown source IDs: " + ",".join(sorted(missing)))
        entries = [r for r in entries if r["id"] in selected]
    result = audit_entries(entries, args.out, args.fasta37, args.fasta38, args.max_records, args.max_sites, args.chromosomes)
    print(f"[data-audit] {len(result['datasets'])} file units; conflicts={result['n_conflicts']}; out={args.out}", flush=True)
    return 2 if result["has_conflicts"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
