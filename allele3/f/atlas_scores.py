#!/usr/bin/env python3
"""Query local AlphaGenome Atlas Tabix scores for bounded GRCh38 SNV targets.

The official AVI table contains one reference-relative score per substitution.
Every target emits the structural REF (score NA) and all three theoretical ALT
bases. Listing an ALT in --targets is not a new population observation. This
module neither discovers alleles nor estimates a distance between ALT alleles.
It performs local CPU-only indexed reads; it never downloads or scans the atlas.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import sys

from analysis_integrity import InputPath, check_output_paths, identity, transactional, validate_manifest

VERSION = "atlas-scores-1.0"
BASES = "ACGT"
MISSING = {"", ".", "NA", "NAN", "NONE", "NULL"}
PORTAL = "https://alphagenome.google/downloads"
DEFAULT_RESOURCES = Path(__file__).resolve().parent.parent / "config" / "alphagenome_resources.tsv"
AVI_COLUMNS = ["CHROM", "POS", "REF", "ALT", "raw_score", "PHRED"]
VARIANT_FIELDS = ["site_id", "chrom", "pos", "build", "ref", "target_alts", "allele", "allele_role",
                  "observed_from_targets", "avi_raw", "avi_phred", "status", "resource_id", "source_version"]
METRIC_FIELDS = ["site_id", "chrom", "pos", "build", "ref", "allele", "allele_role", "observed_from_targets",
                 "resource_id", "source_version", "metric", "value", "status"]
AVAILABILITY_FIELDS = ["site_id", "chrom", "pos", "ref", "target_alts", "avi_available_substitutions",
                       "avi_complete_substitutions", "avi_unavailable_substitutions", "reference_check"]
RESOURCE_FIELDS = ["resource_id", "kind", "status", "required", "enabled", "path", "index_path", "build",
                   "coordinate_system", "columns", "metric_columns", "source_version", "portal_url",
                   "download_url", "url_kind", "index_check", "header_check"]
LIMITS = [
    "AVI raw and PHRED are reference-relative scalar predictions, not allele-to-allele functional distances, fitness, pleiotropy counts or tests of balancing selection. No pairwise score differences are reported.",
    "All three possible nonreference nucleotides are queried. observed_from_targets=1 means the nucleotide was listed in the supplied targets; it does not establish AC, AN, AF or a population observation.",
    "The structural REF is a sequence label and has no REF-to-REF substitution score. Its score and observed_from_targets are NA. REF is not assumed ancestral or the most frequent allele.",
    "Missing rows, missing metrics and absent indexed contigs remain NA. A numeric zero is retained as zero. Neither missing nor negative raw scores imply absence of the allele.",
    "GRCh38 declarations, recognized source headers, Tabix indexing convention and exact REF at returned target positions are checked. This is not a whole-genome FASTA audit; a target with no returned row has no sequence-level REF verification.",
    "User-supplied source versions and data/index pairing are recorded as declarations. Local file identity uses size and mtime for large atlas files; it is not an official release checksum, an exhaustive stale-index test or proof that an atlas was downloaded intact.",
    "Optional generic metrics require an explicit local SNV schema and retain their original units and semantics. They are not combined into a pleiotropy or functional-distance score.",
    "Only requested target positions are read through Tabix. No API calls, network downloads, GPU inference, liftover, strand complementation, mutation-age inference or third-allele discovery are performed.",
]


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
    if "{" in value or "}" in value or re.match(r"[A-Za-z][A-Za-z0-9+.-]*://", value):
        raise ValueError(f"Atlas paths must be local files, with only {{dir_data}} as a placeholder: {value}")
    path = Path(value).expanduser()
    return str((path if path.is_absolute() else config_dir / path).resolve())


def choose_index(path, explicit):
    if explicit:
        if not str(explicit).endswith((".tbi", ".csi")):
            raise ValueError("--avi-index/index_path must name the matching .tbi or .csi file")
        return explicit
    existing = [str(path) + suffix for suffix in (".tbi", ".csi") if Path(str(path) + suffix).is_file()]
    if len(existing) > 1:
        raise ValueError(f"Both .tbi and .csi exist for {path}; select the matching index explicitly")
    return existing[0] if existing else ""


def load_resources(args):
    config = Path(args.resources).expanduser().resolve()
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
        row["path"] = local_path(row["path"], args.dir_data, config.parent)
        row["index_path"] = local_path(row.get("index_path", ""), args.dir_data, config.parent)
        if row["kind"] == "avi":
            if args.avi_file:
                row["path"] = str(Path(args.avi_file).expanduser().resolve())
                row["index_path"] = ""  # A configured index must never follow a different data override.
            if args.avi_index:
                row["index_path"] = str(Path(args.avi_index).expanduser().resolve())
            if args.source_version:
                row["source_version"] = args.source_version
            row["enabled"], row["required"] = True, True
        row.setdefault("source_version", "user-declared-unverified")
        row["source_version"] = row["source_version"] or "user-declared-unverified"
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
        if not row["path"] or not Path(row["path"]).is_file():
            if row["required"]:
                missing.append(f"Required Atlas file missing ({rid}): {row['path'] or 'path not configured in ' + str(config)}")
            row["status"] = "deferred_missing_file"
        selected_index = choose_index(row["path"], row["index_path"])
        if not selected_index or not Path(selected_index).is_file():
            if row["required"]:
                expected = row["index_path"] or f"{row['path']}.tbi OR {row['path']}.csi"
                missing.append(f"Matching Tabix .tbi/.csi missing ({rid}): {expected}")
            row.setdefault("status", "deferred_missing_index")
        elif "status" not in row:
            row["status"] = "ready"
        row["index_path"] = selected_index
        specs.append(row)
    if sum(r["kind"] == "avi" for r in specs) != 1:
        raise ValueError("Resource manifest requires exactly one AVI row")
    if missing:
        raise ValueError("\n".join(missing) + "\nUse the extracted BGZF TSV and its matching index; ZIP archives cannot be queried directly. Override paths with --avi-file/--avi-index if needed.")
    return specs


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
    if fmt != 0:
        raise ValueError(f"Resource {spec['resource_id']}: expected generic TSV Tabix format, got preset {fmt}")
    expected_seq, expected_pos = spec["column_indices"]["chrom"] + 1, spec["column_indices"]["pos"] + 1
    if seq != expected_seq or begin != expected_pos or end not in {0, expected_pos}:
        raise ValueError(f"Resource {spec['resource_id']}: Tabix index columns {(seq, begin, end)} conflict with SNV schema {(expected_seq, expected_pos, expected_pos)}")
    if skip < 0 or name_bytes < 1:
        raise ValueError(f"Resource {spec['resource_id']}: invalid Tabix metadata")
    return dict(format=kind, coordinate_system="1-based", chrom_column=seq, start_column=begin, end_column=end, meta_character=meta, skip_lines=skip)


def normalized_column(value):
    value = re.sub(r"[^a-z0-9]", "", value.lower().lstrip("#"))
    aliases = {"chr": "chrom", "chromosome": "chrom", "position": "pos", "start": "pos",
               "reference": "ref", "alternate": "alt", "raw": "rawscore", "aviraw": "rawscore",
               "avirawscore": "rawscore", "avi": "rawscore", "aviscore": "rawscore",
               "aviphred": "phred", "aviphredscore": "phred", "phredscore": "phred"}
    return aliases.get(value, value)


def check_headers(lines, spec):
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
    return value


class AtlasReader:
    def __init__(self, spec):
        try:
            import pysam
        except ImportError as exc:
            raise RuntimeError("atlas requires pysam for bounded local Tabix reads; install the declared CPU environment") from exc
        self.spec, self.audit = spec, Counter()
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
        spec, found = self.spec, {}
        self.audit["target_queries"] += 1
        if target["chrom"] not in self.contigs:
            self.audit["absent_contig_queries"] += 1
            return found, "not_indexed_contig", "unverified_no_returned_rows"
        self.audit["tabix_fetches"] += 1
        for line in self.handle.fetch(self.contigs[target["chrom"]], target["pos"] - 1, target["pos"]):
            self.audit["returned_rows"] += 1
            cells = line.rstrip("\r\n").split("\t")
            if len(cells) != len(spec["column_names"]):
                raise ValueError(f"Resource {spec['resource_id']}: expected {len(spec['column_names'])} columns at {target['chrom']}:{target['pos']}")
            values = dict(zip(spec["column_names"], cells))
            chrom = chromosome(values[spec["chrom_column"]])
            pos_text = values[spec["pos_column"]]
            if not re.fullmatch(r"[0-9]+", pos_text) or int(pos_text) < 1:
                raise ValueError(f"Resource {spec['resource_id']}: invalid 1-based position {pos_text!r}")
            pos = int(pos_text)
            if chrom != target["chrom"] or pos != target["pos"]:
                raise ValueError(f"Resource {spec['resource_id']}: fetched row coordinate differs from the exact point query; index/data mismatch")
            ref, alt = values[spec["ref_column"]], values[spec["alt_column"]]
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES or ref == alt:
                raise ValueError(f"Resource {spec['resource_id']}: source is not a canonical nonreference SNV at {chrom}:{pos}")
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
        self.audit["queried_substitutions_found"] += len(found)
        return found, "not_scored", "exact_REF_checked" if found else "unverified_no_returned_rows"


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
    result = []
    for spec in request["resources"]:
        # Disabled configurations are documented but never opened or protected as inputs.
        if spec["enabled"]:
            result.extend(InputPath(spec[k]) for k in ("path", "index_path") if spec.get(k))
    return result


@transactional("manifest.json", input_files=("targets_path", "resources_path"), input_prefixes=(),
               resolve_inputs=lambda raw: protected_inputs(raw["request"]))
def publish(out, targets_path, resources_path, request):
    out = Path(out)
    specs = request["resources"]
    counts = Counter(target_sites=len(request["targets"]))
    readers = {}
    with ExitStack() as stack:
        for spec in specs:
            if spec["status"] == "ready":
                reader = stack.enter_context(AtlasReader(spec))
                readers[spec["resource_id"]] = reader
                print(f"[build-check] atlas={spec['resource_id']} declared=GRCh38 source={spec['header_check']['verification']} index={spec['index_check']['format']}:1-based exact_REF=checked_per_returned_target", flush=True)
            print(f"[dataset] atlas={spec['resource_id']} status={spec['status']} source_version={spec['source_version']} file={spec['path'] or 'not_configured'}", flush=True)
        avi_spec = next(s for s in specs if s["kind"] == "avi")
        vh = stack.enter_context(gzip.open(out / "atlas_variants.tsv.gz", "wt", encoding="utf-8", newline=""))
        mh = stack.enter_context(gzip.open(out / "atlas_metrics.tsv.gz", "wt", encoding="utf-8", newline=""))
        ah = stack.enter_context((out / "target_availability.tsv").open("w", encoding="utf-8", newline=""))
        vw, mw, aw = writer(vh, VARIANT_FIELDS), writer(mh, METRIC_FIELDS), writer(ah, AVAILABILITY_FIELDS)
        for target in request["targets"]:
            common = dict(site_id=target["site_id"], chrom=target["chrom"], pos=target["pos"], build="GRCh38", ref=target["ref"], target_alts=",".join(target["alts"]))
            queries = {rid: reader.query(target) for rid, reader in readers.items()}
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
                                         value="NA" if value is None else value, status=metric_status))
                        counts["metric_rows"] += 1
            counts["avi_sites_with_any_score"] += available > 0
            counts["avi_sites_with_all_three_complete"] += complete == 3
            aw.writerow(dict(common, avi_available_substitutions=available, avi_complete_substitutions=complete,
                             avi_unavailable_substitutions=3 - available, reference_check=ref_check))
    with (out / "resources.tsv").open("w", encoding="utf-8", newline="") as handle:
        rw = writer(handle, RESOURCE_FIELDS)
        for spec in specs:
            rendered = dict(spec)
            for key in ("index_check", "header_check"):
                rendered[key] = json.dumps(spec.get(key, {}), sort_keys=True, separators=(",", ":"))
            rw.writerow(rendered)
    usable = counts["avi_status_available"] + counts["avi_status_available_partial"]
    deferred = any(s["enabled"] and s["status"] != "ready" for s in specs)
    status = "not_available" if not usable else "complete" if counts["avi_status_available"] == 3 * len(request["targets"]) and not deferred else "complete_with_exclusions"
    return dict(status=status, schema=VERSION, signature=request["signature"], counts=dict(counts),
                target_selection=request["selection"], parameters=request["parameters"], source_identity=request["source_identity"],
                resources=specs, query_audits={rid: dict(reader.audit) for rid, reader in readers.items()},
                limitations=LIMITS, genome_build="GRCh38", coordinate_system="1-based",
                observed_from_targets_semantics="ALT membership in input alts; not an inferred population observation; REF is NA",
                source_version_semantics="user declaration, not independently verified official release",
                memory_policy="bounded target list plus three substitutions per resource and target; local Tabix point queries only; no full atlas scan")


def run(args):
    declared_build = genome_build(args.build)
    if declared_build != "GRCh38":
        raise ValueError("AlphaGenome Atlas SNV resources require GRCh38; harmonize targets with verified liftover/REF checks first")
    if args.max_targets < 0:
        raise ValueError("--max-targets must be >= 0; 0 explicitly removes the candidate cap")
    targets_path, resources_path = Path(args.targets).expanduser().resolve(), Path(args.resources).expanduser().resolve()
    critical_sources = {"targets": identity(targets_path, full=True), "resources_manifest": identity(resources_path, full=True)}
    specs = load_resources(args)
    selected_chromosomes, selected_region = chromosome_filter(args.chromosomes), region_filter(args.region)
    targets, selection = load_targets(targets_path, declared_build, args.max_targets, selected_chromosomes, selected_region)
    parameters = dict(build=declared_build, chromosomes=sorted(selected_chromosomes, key=chromosome_key) if selected_chromosomes is not None else "all",
                      region=list(selected_region) if selected_region is not None else None, max_targets=args.max_targets,
                      theoretical_alts="all_three", reference_score="NA")
    if critical_sources["targets"] != identity(targets_path, full=True) or critical_sources["resources_manifest"] != identity(resources_path, full=True):
        raise ValueError("Target/schema inputs changed while preparing the Atlas request")
    sources = dict(critical_sources, module=identity(__file__, full=True),
                   analysis_integrity=identity(Path(__file__).with_name("analysis_integrity.py"), full=True))
    for spec in specs:
        if spec["enabled"]:
            for kind in ("path", "index_path"):
                if spec.get(kind):
                    sources[f"{spec['resource_id']}:{kind}"] = identity(spec[kind])
    # Include resolved schemas and paths, not only config identity: CLI overrides matter.
    signature_payload = dict(schema=VERSION, parameters=parameters, sources=sources, resources=specs)
    signature = hashlib.sha256(json.dumps(signature_payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    request = dict(targets=targets, resources=specs, selection=selection, parameters=parameters, source_identity=sources, signature=signature)
    check_output_paths(args.out, [InputPath(targets_path), InputPath(resources_path), *protected_inputs(request)])
    print(f"[build-check] atlas_targets={declared_build} coordinates=1-based selected_sites={len(targets)} theoretical_ALT_per_site=3 REF_score=NA", flush=True)
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
    p.add_argument("--avi-file", help="Official local six-column BGZF TSV; override the configured AVI path")
    p.add_argument("--avi-index", help="Explicit index paired with --avi-file (.tbi or .csi); required if both siblings exist")
    p.add_argument("--resources", default=str(DEFAULT_RESOURCES), help="Explicit local resource schema TSV; optional metrics can remain disabled/missing")
    p.add_argument("--dir-data", default=os.environ.get("DIR_DATA", "/mnt/f"), help="Expand {dir_data}; defaults to the shared local/HPC data root")
    p.add_argument("--source-version", help="User-declared AVI release identifier from the download portal or receipt; never guessed from filename")
    p.add_argument("--chromosomes", default="all", help="Filter targets, e.g. 1-22,X,Y or chr1; does not scan chromosomes")
    p.add_argument("--region", help="Filter targets with 1-based closed chr:start-end; does not scan the region")
    p.add_argument("--max-targets", type=int, default=5000, help="Selected target cap; 0 explicitly removes it. Use independent target batches on HPC")
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


if __name__ == "__main__":
    raise SystemExit(main())
