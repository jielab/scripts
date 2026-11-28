#!/usr/bin/env python3
"""Shared I/O, interval masks, input caching and genome-coordinate mapping."""
from __future__ import annotations


# ========================================================================
# Core
# ========================================================================

# Small, dependency-free utilities for the local allele3 workflow.


import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

VERSION = "4.0.0-focused-20261024"


def open_text(path, mode="rt"):
    path = Path(path)
    if path.name.endswith((".gz", ".bgz", ".bgzip")):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode.replace("t", ""), encoding="utf-8", newline="")


def dump_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(value, f, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")
    os.replace(tmp, path)


def fingerprint(path, sample_bytes=65536):
    """Detect ordinary input changes without hashing hundreds of GB twice.

    This is deliberately labelled a metadata/edge fingerprint, NOT a full-file
    cryptographic integrity check. Download checksums remain the user's source
    integrity check. Small files are hashed completely.
    """
    p = Path(path).resolve()
    s = p.stat()
    h = hashlib.sha256()
    with p.open("rb") as f:
        h.update(f.read(sample_bytes))
        if s.st_size > sample_bytes:
            f.seek(max(sample_bytes, s.st_size - sample_bytes))
            h.update(f.read(sample_bytes))
    return {"path": str(p), "size": s.st_size, "mtime_ns": s.st_mtime_ns,
            "edge_sha256": h.hexdigest()}


def stable_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def code_hash():
    h = hashlib.sha256()
    for path in sorted(Path(__file__).parent.glob("*.py")):
        h.update(path.name.encode())
        h.update(path.read_bytes())
    return h.hexdigest()


def tsv_writer(path, fields):
    """Return a handle/writer; caller closes handle before committing a shard."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.name.endswith((".gz", ".bgz")):
        f = gzip.open(p, "wt", encoding="utf-8", newline="", compresslevel=2)
    else:
        f = p.open("w", encoding="utf-8", newline="")
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n",
                            extrasaction="raise", restval=".")
    writer.writeheader()
    return f, writer


def fmt(value):
    if value is None:
        return "."
    if isinstance(value, float):
        return format(value, ".10g")
    return value


# ========================================================================
# Gnomad io
# ========================================================================

# Bounded-memory local I/O for gnomAD sites VCFs.
# 
# VCF positions are one-based.  No record filtering, allele merging, inference of
# reference counts, or recalculation of INFO counts occurs in this module.
# ``bcftools query`` is used when available; a standard-library VCF/gzip reader is
# provided for portability.  Neither backend writes a decompressed intermediate.


import gzip
import math
import os
from collections import OrderedDict
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import warnings


AUTOSOMES = tuple(str(i) for i in range(1, 23))
_FIELD = re.compile(r"^[A-Za-z_][A-Za-z0-9_.]*$")
_SUFFIXES = (".vcf", ".vcf.gz", ".vcf.bgz", ".bcf")


def canonical_chrom(chrom):
    """Normalize only conventional chromosome aliases; retain other contig IDs."""
    chrom = str(chrom)
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom.isdigit():
        return str(int(chrom))
    if chrom in ("M", "MT"):
        return "MT"
    if not chrom:
        raise ValueError("Empty chromosome identifier")
    return chrom


def chrom_key(chrom):
    chrom = canonical_chrom(chrom)
    if chrom.isdigit():
        return (0, int(chrom), "")
    if chrom in ("X", "Y", "MT"):
        return (1, {"X": 23, "Y": 24, "MT": 25}[chrom], "")
    return (2, 0, chrom)


def parse_region(region):
    """Return (canonical chromosome, start, end), inclusive and one-based."""
    if region is None or region == "":
        return None
    match = re.fullmatch(r"([^:\s,]+):([1-9][0-9]*)-([1-9][0-9]*)", str(region))
    if not match or int(match[2]) > int(match[3]):
        raise ValueError("Region must be CHR:START-END with 1 <= START <= END")
    return canonical_chrom(match[1]), int(match[2]), int(match[3])


def _chromosomes(value):
    if value is None or value == "":
        return None
    if isinstance(value, str):
        value = value.split(",")
    result = []
    for item in value:
        item = canonical_chrom(str(item).strip())
        if item == "autosomes":
            result.extend(AUTOSOMES)
        elif item == "all":
            result.extend((*AUTOSOMES, "X", "Y"))
        elif item in (*AUTOSOMES, "X", "Y", "MT"):
            result.append(item)
        else:
            raise ValueError("Chromosomes must be autosomes, all, or a comma list of 1..22,X,Y,MT")
    return list(dict.fromkeys(result))


def _file_chrom(path):
    name = Path(path).name
    match = re.search(r"(?:^|[._-])chr([0-9]+|X|Y|MT|M)(?=[._-]|$)", name)
    if not match:
        match = re.match(r"([0-9]+|X|Y|MT|M)\.vcf(?:\.(?:gz|bgz))?$", name)
    return canonical_chrom(match[1]) if match else None


def discover_inputs(vcf_file=None, vcf_dir=None, chromosomes=None):
    """Resolve a single file or nonoverlapping per-chromosome files from one cohort.

    Explicit chromosome selections must be complete.  Multiple files for the
    same chromosome, or multiple files whose chromosome cannot be inferred from
    the filename, are rejected instead of silently combining overlapping data.
    """
    if vcf_file and vcf_dir:
        raise ValueError("Use vcf_file or vcf_dir, not both")
    if vcf_file:
        if chromosomes:
            raise ValueError("Chromosome selection is for directory input; use region with one VCF")
        paths = [Path(vcf_file).expanduser()]
    else:
        if not vcf_dir:
            raise ValueError("A VCF file or directory is required")
        directory = Path(vcf_dir).expanduser()
        if not directory.is_dir():
            raise ValueError(f"VCF directory does not exist: {directory}")
        paths = [p for p in directory.iterdir() if p.is_file() and p.name.endswith(_SUFFIXES)]
        chosen = _chromosomes(chromosomes)
        if chosen is not None:
            paths = [p for p in paths if _file_chrom(p) in chosen]
            missing = set(chosen) - {_file_chrom(p) for p in paths}
            if missing:
                raise ValueError("Missing requested chromosomes: " + ",".join(sorted(missing, key=chrom_key)))
        if len(paths) > 1 and any(_file_chrom(p) is None for p in paths):
            raise ValueError("Cannot establish disjoint chromosomes for directory inputs; use chr*.vcf names or one --vcf")
        seen = {}
        for path in paths:
            chrom = _file_chrom(path)
            if chrom in seen:
                raise ValueError(f"Multiple VCFs for chromosome {chrom}: {seen[chrom]}, {path}; do not merge cohorts/releases")
            seen[chrom] = path
        paths.sort(key=lambda p: (chrom_key(_file_chrom(p) or p.name), p.name))
    if not paths:
        raise ValueError(f"No VCF files found: {vcf_dir}")
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty VCF: {path}")
        if Path(str(path) + ".aria2").exists():
            raise ValueError(f"Download is marked incomplete (.aria2): {path}")
    return paths


def _backend(path, backend):
    if backend not in ("auto", "bcftools", "python"):
        raise ValueError("Backend must be auto, bcftools, or python")
    executable = shutil.which("bcftools")
    if backend == "bcftools" and executable is None:
        raise ValueError("bcftools was requested but is not on PATH")
    selected = "bcftools" if executable and backend != "python" else "python"
    if selected == "python" and str(path).endswith(".bcf"):
        raise ValueError("BCF input requires bcftools; the Python fallback reads VCF text only")
    return selected


def _open_text(path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def _structured_meta(text):
    """Split VCF angle-bracket metadata, respecting quoted commas and escapes."""
    fields, current = [], []
    quoted = escaped = False
    for char in text:
        if escaped:
            current.append(char)
            escaped = False
        elif char == "\\" and quoted:
            current.append(char)
            escaped = True
        elif char == '"':
            quoted = not quoted
            current.append(char)
        elif char == "," and not quoted:
            fields.append("".join(current))
            current = []
        else:
            current.append(char)
    if quoted:
        raise ValueError("Unclosed quote in VCF metadata")
    fields.append("".join(current))
    result = {}
    for field in fields:
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        if value.startswith('"') and value.endswith('"'):
            value = value[1:-1].replace('\\"', '"').replace("\\\\", "\\")
        result[key.strip()] = value
    return result


def _parse_header(raw):
    info, contigs, samples = {}, [], None
    for line in raw.splitlines():
        if line.startswith("##INFO=<") and line.endswith(">"):
            values = _structured_meta(line[len("##INFO=<"):-1])
            if "ID" not in values or "Number" not in values or "Type" not in values:
                raise ValueError("Malformed INFO definition in VCF header")
            entry = {"number": values["Number"], "type": values["Type"], "description": values.get("Description", "")}
            if values["ID"] in info and info[values["ID"]] != entry:
                raise ValueError(f"Conflicting INFO header definitions: {values['ID']}")
            info[values["ID"]] = entry
        elif line.startswith("##contig=<") and line.endswith(">"):
            values = _structured_meta(line[len("##contig=<"):-1])
            if values.get("ID") and values["ID"] not in contigs:
                contigs.append(values["ID"])
        elif line.startswith("#CHROM\t"):
            columns = line.split("\t")
            if columns[:8] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]:
                raise ValueError("Malformed VCF column header")
            samples = max(0, len(columns) - 9)
    if samples is None:
        raise ValueError("VCF #CHROM header is missing")
    return {"info": info, "contigs": contigs, "samples": samples, "raw": raw}


def _error_tail(handle):
    handle.flush()
    handle.seek(0, os.SEEK_END)
    handle.seek(max(0, handle.tell() - 12000))
    return handle.read().decode("utf-8", errors="replace").strip()


def read_header(path, backend="auto"):
    path = Path(path)
    if _backend(path, backend) == "bcftools":
        with tempfile.TemporaryFile() as errors:
            result = subprocess.run(["bcftools", "view", "--no-version", "-h", str(path)],
                                    stdout=subprocess.PIPE, stderr=errors, text=True, check=False)
            if result.returncode:
                raise ValueError(f"Cannot read VCF header ({path}): {_error_tail(errors)}")
            raw = result.stdout
    else:
        lines = []
        with _open_text(path) as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                lines.append(line)
                if line.startswith("#CHROM\t"):
                    break
        raw = "".join(lines)
    return _parse_header(raw)


def validate_count_header(header, ac="AC", an="AN", af="AF"):
    """Validate tags needed to interpret allele counts without guessing shapes."""
    info = header["info"]
    for key in (ac, an):
        if key not in info:
            raise ValueError(f"Required INFO/{key} is absent from the VCF header")
    if info[ac]["number"] != "A" or info[ac]["type"] != "Integer":
        raise ValueError(f"INFO/{ac} must have Number=A, Type=Integer")
    if info[an]["number"] not in ("1", "A") or info[an]["type"] != "Integer":
        raise ValueError(f"INFO/{an} must have Number=1 or A, Type=Integer")
    if af in info and (info[af]["number"] != "A" or info[af]["type"] != "Float"):
        raise ValueError(f"INFO/{af} must have Number=A, Type=Float")


def allele_info_values(record, key, header, *, kind="float", required=False, allow_number_one=False):
    """Return ALT-aligned numeric values or None for missing entries.

    Number=A is checked against the original ALT count, including nonsingle-
    nucleotide ALT entries.  A scalar is broadcast only for a declared Number=1
    field when the caller explicitly sets allow_number_one=True (e.g. AN).
    ``required`` requires the header definition; per-record missing values remain
    None so the caller can count/exclude them rather than invent zeros.
    """
    count = len(record["alts"])
    if key not in header["info"]:
        if required or record["info"].get(key) not in (None, "", "."):
            raise ValueError(f"INFO/{key} has no header definition")
        return [None] * count
    number = header["info"][key]["number"]
    if number != "A" and not (number == "1" and allow_number_one):
        raise ValueError(f"INFO/{key} Number={number} cannot be aligned to ALTs")
    raw = record["info"].get(key, ".")
    if raw is None or raw in ("", "."):
        return [None] * count
    values = str(raw).split(",")
    expected = count if number == "A" else 1
    if len(values) != expected:
        raise ValueError(f"INFO/{key}/ALT cardinality mismatch at {record['chrom']}:{record['pos']}: expected {expected}, got {len(values)}")
    if kind not in ("int", "float"):
        raise ValueError("Numeric kind must be int or float")
    parsed = []
    for value in values:
        if value == ".":
            parsed.append(None)
            continue
        if kind == "int":
            if not re.fullmatch(r"[+-]?[0-9]+", value):
                raise ValueError(f"INFO/{key} is not an integer at {record['chrom']}:{record['pos']}: {value}")
            parsed.append(int(value))
        else:
            try:
                numeric = float(value)
            except ValueError as error:
                raise ValueError(f"Invalid INFO/{key} at {record['chrom']}:{record['pos']}: {value}") from error
            if not math.isfinite(numeric):
                raise ValueError(f"Nonfinite INFO/{key} at {record['chrom']}:{record['pos']}")
            parsed.append(numeric)
    return parsed * count if number == "1" else parsed


def _record(columns, fields):
    if len(columns) != 5 + len(fields):
        raise ValueError("Unexpected number of columns from VCF/query reader")
    try:
        pos = int(columns[1])
    except ValueError as error:
        raise ValueError(f"Invalid VCF position: {columns[1]}") from error
    if pos < 1:
        raise ValueError("VCF positions must be positive")
    return {"chrom": columns[0], "pos": pos, "ref": columns[2],
            "alts": [] if columns[3] == "." else columns[3].split(","),
            "filter": columns[4], "info": dict(zip(fields, columns[5:]))}


def _python_records(path, fields):
    selected = set(fields)
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            values = line.rstrip("\r\n").split("\t", 8)
            if len(values) < 8:
                raise ValueError(f"VCF record has fewer than eight columns: {path}")
            info = {}
            for item in values[7].split(";"):
                key, sep, value = item.partition("=")
                if key in selected:
                    if key in info:
                        raise ValueError(f"Duplicate INFO/{key} at {values[0]}:{values[1]}")
                    info[key] = value if sep else "1"
            columns = [values[0], values[1], values[3], values[4], values[6]]
            yield _record(columns + [info.get(key, ".") for key in fields], fields)


def _bcftools_records(path, fields, region):
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER" + "".join("\t%INFO/" + key for key in fields) + "\n"
    command = ["bcftools", "query", "--allow-undef-tags", "-f", fmt]
    if region:
        command.extend(["-r", region])
    command.append(str(path))
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors, text=True)
        exhausted = False
        try:
            for line in process.stdout:
                yield _record(line.rstrip("\r\n").split("\t"), fields)
            exhausted = True
            status = process.wait()
            if status:
                raise ValueError(f"bcftools query failed ({path}, exit {status}): {_error_tail(errors)}")
        finally:
            process.stdout.close()
            if not exhausted and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()


def iterate_records(path, fields, region=None, backend="auto", threads=1):
    """Yield raw selected INFO strings with no whole-file materialization.

    Sorting follows ##contig order if supplied, otherwise natural chromosome
    order.  Equal positions are expected for split multiallelic records.  A
    requested region is inclusive; records overlapping it but starting outside
    it are not returned.  ``threads`` is accepted for a uniform caller API;
    bcftools query remains streaming and parallelism is across input files.
    """
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    fields = list(fields)
    if len(set(fields)) != len(fields) or any(not _FIELD.fullmatch(x) for x in fields):
        raise ValueError("Selected INFO names must be unique literal VCF tags")
    path = Path(path)
    selected_backend = _backend(path, backend)
    bounds = parse_region(region)
    if bounds and selected_backend == "bcftools" and not any(Path(str(path) + x).is_file() for x in (".tbi", ".csi")):
        if backend == "bcftools" or str(path).endswith(".bcf"):
            raise ValueError(f"Indexed region query requires a .tbi/.csi index: {path}")
        selected_backend = "python"
    if bounds and selected_backend == "python":
        warnings.warn("Region selected with Python fallback: scanning the VCF sequentially (no indexed seek)", RuntimeWarning, stacklevel=2)
    header = read_header(path, selected_backend)
    ranks = {canonical_chrom(c): i for i, c in enumerate(header["contigs"])}
    if len(ranks) != len(header["contigs"]):
        raise ValueError("VCF header has ambiguous chromosome aliases (e.g. 1 and chr1)")
    query_region = None
    if bounds:
        chrom, start, end = bounds
        matches = [c for c in header["contigs"] if canonical_chrom(c) == chrom]
        if header["contigs"] and not matches:
            raise ValueError(f"Region chromosome is absent from VCF header: {chrom}")
        raw_chrom = matches[0] if matches else str(region).split(":", 1)[0]
        query_region = f"{raw_chrom}:{start}-{end}"
    records = (_bcftools_records(path, fields, query_region) if selected_backend == "bcftools"
               else _python_records(path, fields))
    previous = None
    try:
        for record in records:
            chrom = canonical_chrom(record["chrom"])
            order = (0, ranks[chrom], "", record["pos"]) if chrom in ranks else (1, *chrom_key(chrom), record["pos"])
            if previous is not None and order < previous:
                raise ValueError(f"VCF must be coordinate sorted: {path} at {record['chrom']}:{record['pos']}")
            previous = order
            if bounds and (chrom != bounds[0] or not bounds[1] <= record["pos"] <= bounds[2]):
                continue
            yield record
    finally:
        records.close()


class FastaContext:
    """Random-access trinucleotide contexts with a bounded 1 MiB block cache.

    An existing .fai is required.  Compressed FASTA additionally requires pysam,
    a bgzip-compressed file and its .gzi; no implicit indexing/decompression is
    performed.  The standard-library backend supports ordinary FASTA with LF
    or CRLF line endings.  Two blocks are retained, including boundary lookups.
    """
    def __init__(self, fasta):
        self.path = Path(fasta)
        self._cache = OrderedDict()
        self._block_size = 1 << 20
        self._max_blocks = 2
        self._handle = self._pysam = None
        fai = Path(str(self.path) + ".fai")
        if not self.path.is_file() or not fai.is_file():
            raise ValueError(f"FASTA and existing .fai are required: {self.path}; run samtools faidx first")
        self._entries = {}
        with fai.open() as handle:
            for line in handle:
                values = line.rstrip("\r\n").split("\t")
                if len(values) < 5:
                    raise ValueError(f"Invalid FASTA index: {fai}")
                name = values[0]
                length, offset, bases, width = map(int, values[1:5])
                canonical = canonical_chrom(name)
                if canonical in self._entries:
                    raise ValueError(f"Ambiguous FASTA chromosome aliases: {name}")
                if length < 1 or offset < 0 or bases < 1 or width < bases:
                    raise ValueError(f"Invalid FASTA index entry: {name}")
                self._entries[canonical] = (name, length, offset, bases, width)
        if not self._entries:
            raise ValueError("FASTA index contains no sequences")
        self.lengths = {key: entry[1] for key, entry in self._entries.items()}
        compressed = str(self.path).endswith((".gz", ".bgz"))
        try:
            import pysam
        except ImportError:
            pysam = None
        if compressed:
            if pysam is None or not Path(str(self.path) + ".gzi").is_file():
                raise ValueError("Compressed FASTA requires pysam and pre-existing .fai/.gzi indexes; alternatively provide uncompressed FASTA")
            try:
                self._pysam = pysam.FastaFile(str(self.path))
            except (OSError, ValueError) as error:
                raise ValueError(f"Cannot seek compressed FASTA (must be bgzip, not ordinary gzip): {error}") from error
        elif pysam is not None:
            self._pysam = pysam.FastaFile(str(self.path))
        else:
            self._handle = self.path.open("rb")

    def _read(self, chrom, start, end):
        name, length, offset, bases, width = self._entries[chrom]
        if not 0 <= start <= end <= length:
            raise ValueError("Internal FASTA interval is out of range")
        if start == end:
            return ""
        if self._pysam is not None:
            return self._pysam.fetch(name, start, end).upper()
        first = offset + (start // bases) * width + start % bases
        last = offset + ((end - 1) // bases) * width + (end - 1) % bases
        self._handle.seek(first)
        sequence = self._handle.read(last - first + 1).replace(b"\r", b"").replace(b"\n", b"").decode("ascii").upper()
        if len(sequence) != end - start:
            raise ValueError(f"FASTA/index disagreement at {name}:{start+1}-{end}")
        return sequence

    def _block(self, chrom, number):
        key = (chrom, number)
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        start = number * self._block_size
        sequence = self._read(chrom, start, min(start + self._block_size, self.lengths[chrom]))
        self._cache[key] = sequence
        while len(self._cache) > self._max_blocks:
            self._cache.popitem(last=False)
        return sequence

    def context(self, chrom, pos):
        chrom = canonical_chrom(chrom)
        if chrom not in self.lengths:
            raise KeyError(f"Chromosome absent from FASTA: {chrom}")
        if not isinstance(pos, int) or not 1 <= pos <= self.lengths[chrom]:
            raise ValueError(f"Position outside FASTA: {chrom}:{pos}")
        result = []
        for index in (pos - 2, pos - 1, pos):
            if index < 0 or index >= self.lengths[chrom]:
                result.append("N")
            else:
                number, offset = divmod(index, self._block_size)
                result.append(self._block(chrom, number)[offset])
        return "".join(result)

    def close(self):
        if self._pysam is not None:
            self._pysam.close()
            self._pysam = None
        if self._handle is not None:
            self._handle.close()
            self._handle = None
        self._cache.clear()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# ========================================================================
# Common
# ========================================================================

# Adapters reusing the original scanner's chromosome and I/O conventions.

import math
chrom_sort_key = chrom_key

def af_bin(x):
    try:
        v = float(x)
    except (ValueError, TypeError):
        return 'af_missing'
    if not math.isfinite(v) or not 0 <= v <= 1:
        return 'af_missing'
    for edge, name in [(1e-5,'lt_1e-5'),(1e-4,'1e-5_to_1e-4'),(1e-3,'1e-4_to_1e-3'),(.01,'0.001_to_0.01'),(.05,'0.01_to_0.05'),(.5,'0.05_to_0.5')]:
        if v < edge: return name
    return 'ge_0.5'

def an_bin(x):
    try:
        v = float(x)
        if not math.isfinite(v) or v < 0: return 'an_missing'
    except (ValueError, TypeError):
        return 'an_missing'
    for edge in (100,1000,10000,50000,100000,120000,140000,150000,200000):
        if v < edge: return 'lt_' + str(edge)
    return 'ge_200000'

def ac_tier(x):
    try: v=float(x)
    except (ValueError,TypeError): return 'ac_missing'
    if not math.isfinite(v) or v < 0: return 'ac_missing'
    if v == 0: return '0'
    if v == 1: return '1'
    if v < 5: return '2-4'
    if v < 10: return '5-9'
    if v < 100: return '10-99'
    return '100+'


# ========================================================================
# Intervals
# ========================================================================

# Callable BED handling and streaming conversion of gnomAD all-sites AN.


from array import array
from bisect import bisect_right
import csv
import math
import os
from pathlib import Path


def canonical(chrom):
    value = str(chrom).removeprefix("chr")
    if value.isdigit():
        return str(int(value))
    return "MT" if value in ("M", "MT") else value


def order(chrom):
    c = canonical(chrom)
    return (0, int(c)) if c.isdigit() else (1, {"X": 23, "Y": 24, "MT": 25}.get(c, 100), c)


class IntervalMask:
    """Merged 0-based half-open intervals, compact arrays for point lookup."""

    def __init__(self, path, allowed_chroms=None):
        allowed = None if allowed_chroms is None else {canonical(c) for c in allowed_chroms}
        self.starts, self.ends = {}, {}
        intervals = {}
        with open_text(path) as f:
            for line_no, line in enumerate(f, 1):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                cells = line.rstrip("\n").split("\t")
                if len(cells) < 3:
                    raise ValueError(f"Invalid BED at {path}:{line_no}")
                chrom = canonical(cells[0])
                try:
                    start, end = int(cells[1]), int(cells[2])
                except ValueError as e:
                    raise ValueError(f"BED must use integer 0-based half-open coordinates: {path}:{line_no}") from e
                if start < 0 or end <= start:
                    raise ValueError(f"Invalid BED interval at {path}:{line_no}")
                if allowed is None or chrom in allowed:
                    intervals.setdefault(chrom, []).append((start, end))
        for chrom, rows in intervals.items():
            rows.sort()
            merged = []
            for start, end in rows:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.starts[chrom] = array("Q", (s for s, _ in merged))
            self.ends[chrom] = array("Q", (e for _, e in merged))

    def contains(self, chrom, pos1):
        c, point = canonical(chrom), int(pos1) - 1
        starts = self.starts.get(c)
        if starts is None:
            return False
        i = bisect_right(starts, point) - 1
        return i >= 0 and point < self.ends[c][i]


def make_callable(all_sites_an, output_bed, min_an=100000, chromosomes=None):
    """Stream a scalar AN table into merged BED; never builds a per-base map.

    Supported headers: locus,AN (Hail export: chr1:123) or chrom,pos,AN.
    Input coordinates are always 1-based; output BED is 0-based half-open.
    AN is a coverage criterion, not a reconstructed joint allele denominator.
    """
    if min_an < 1:
        raise ValueError("min_an must be positive")
    source, output = Path(all_sites_an), Path(output_bed)
    if source.resolve() == output.resolve():
        raise ValueError("Input and output must be different files")
    if output.exists():
        raise FileExistsError(f"Output exists: {output}; choose a new path")
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(output.name + ".partial")
    allowed = None if not chromosomes else {canonical(c) for c in chromosomes}
    n_rows = n_kept = n_intervals = 0
    previous = None
    current = None
    try:
        with open_text(source) as inp, tmp.open("w", encoding="utf-8") as out:
            lines = (line for line in inp if line.strip() and not line.startswith("##"))
            reader = csv.DictReader(lines, delimiter="\t")
            names = {n.lstrip("#").lower(): n for n in (reader.fieldnames or [])}
            an_key = names.get("an")
            locus_key = names.get("locus")
            chrom_key = next((names[x] for x in ("chrom", "chr", "contig", "chromosome") if x in names), None)
            pos_key = next((names[x] for x in ("pos", "position") if x in names), None)
            if not an_key or not (locus_key or (chrom_key and pos_key)):
                raise ValueError("all-sites AN table requires header locus+AN or chrom+pos+AN; AN must be scalar")
            for row in reader:
                n_rows += 1
                if locus_key:
                    try:
                        raw_chrom, raw_pos = row[locus_key].rsplit(":", 1)
                    except ValueError as e:
                        raise ValueError(f"Invalid Hail locus on data row {n_rows}") from e
                else:
                    raw_chrom, raw_pos = row[chrom_key], row[pos_key]
                chrom, pos = canonical(raw_chrom), int(raw_pos)
                if pos < 1:
                    raise ValueError("all-sites positions must be 1-based positive integers")
                key = (order(chrom), pos)
                if previous is not None and key <= previous:
                    raise ValueError("all-sites AN must be naturally chromosome/position sorted, without duplicates")
                previous = key
                try:
                    an_value = float(row[an_key])
                except ValueError as e:
                    raise ValueError("AN must be a scalar number, not a Hail array; use the official TSV export") from e
                if not math.isfinite(an_value) or an_value < 0 or not an_value.is_integer():
                    raise ValueError(f"Invalid AN on row {n_rows}")
                if allowed is not None and chrom not in allowed:
                    continue
                if an_value < min_an:
                    continue
                n_kept += 1
                if current is not None and chrom == current[0] and pos - 1 == current[2]:
                    current[2] = pos
                else:
                    if current is not None:
                        out.write("\t".join(map(str, current)) + "\n")
                        n_intervals += 1
                    current = [chrom, pos - 1, pos]
            if current is not None:
                out.write("\t".join(map(str, current)) + "\n")
                n_intervals += 1
        os.replace(tmp, output)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    result = {"input": fingerprint(source), "output": str(output.resolve()),
              "min_an": min_an, "input_rows": n_rows, "callable_bases": n_kept,
              "merged_intervals": n_intervals, "coordinates": "BED 0-based half-open",
              "warning": "Overall pre-split AN defines the mask only; it does not recover joint REF/ALT counts."}
    dump_json(str(output) + ".json", result)
    return result


# ========================================================================
# Input cache
# ========================================================================

# Verified optional PGEN execution copies on a fast local filesystem.
# 
# Scientific outputs stay in the analysis directory. Original inputs are never
# modified. A complete cache requires source fingerprints and full copy hashes.

import fcntl
import hashlib
import os
from pathlib import Path
import shutil


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        released = consumed = 0
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
            consumed += len(block)
            if consumed - released >= 256 * 1024 * 1024 and hasattr(os, 'posix_fadvise'):
                os.posix_fadvise(handle.fileno(), released, consumed - released, os.POSIX_FADV_DONTNEED)
                released = consumed
        if hasattr(os, 'posix_fadvise'):
            os.posix_fadvise(handle.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
    return h.hexdigest()


def cache_pgen_prefix(prefix, cache_dir):
    prefix, root = Path(prefix).resolve(), Path(cache_dir).expanduser().resolve()
    sources = {suffix: Path(str(prefix) + suffix) for suffix in ('.pgen', '.pvar', '.psam')}
    identity = {suffix: fingerprint(path) for suffix, path in sources.items()}
    root.mkdir(parents=True, exist_ok=True)
    destination = root / (prefix.name + '.' + stable_hash(identity)[:20])
    with (root / '.input_cache.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        manifest = destination / 'cache.json'
        if destination.exists():
            if not manifest.is_file():
                raise ValueError(f'Incomplete/unrecognized input cache: {destination}')
            saved = __import__('json').loads(manifest.read_text())
            if saved['source_fingerprints'] != identity:
                raise ValueError('PGEN cache source identity changed')
            for suffix, spec in saved['files'].items():
                path = destination / (prefix.name + suffix)
                if fingerprint(path) != spec['fingerprint'] or _sha(path) != spec['sha256']:
                    raise ValueError(f'PGEN cache checksum differs: {path}')
            return str(destination / prefix.name), saved
        staging = destination.with_name(destination.name + '.partial')
        if staging.exists():
            shutil.rmtree(staging)  # only this helper's exact incomplete cache
        required = sum(x['size'] for x in identity.values())
        if shutil.disk_usage(root).free < required + 1024**3:
            raise ValueError(f'PGEN cache needs {required / 1024**3:.1f} GiB plus 1 GiB free at {root}')
        staging.mkdir()
        copied = {}
        try:
            for suffix, source in sources.items():
                target = staging / (prefix.name + suffix)
                h = hashlib.sha256()
                count = dirty = 0
                with source.open('rb') as src, target.open('wb') as dst:
                    for block in iter(lambda: src.read(8 * 1024 * 1024), b''):
                        dst.write(block); h.update(block)
                        count += len(block); dirty += len(block)
                        if dirty >= 256 * 1024 * 1024:
                            dst.flush(); os.fdatasync(dst.fileno())
                            if hasattr(os, 'posix_fadvise'):
                                os.posix_fadvise(dst.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
                            dirty = 0
                        if count % (1024**3) < len(block):
                            print(f'[PGEN cache] {suffix}: {count / 1024**3:.1f} GiB copied', flush=True)
                    dst.flush(); os.fsync(dst.fileno())
                if fingerprint(source) != identity[suffix]:
                    raise ValueError(f'Source changed while copying: {source}')
                if _sha(target) != h.hexdigest():
                    raise ValueError(f'PGEN execution copy failed checksum: {target}')
                copied[suffix] = {'sha256': h.hexdigest()}
            os.replace(staging, destination)
            for suffix in copied:
                copied[suffix]['fingerprint'] = fingerprint(destination / (prefix.name + suffix))
            record = dict(schema=1, source_prefix=str(prefix), cached_prefix=str(destination / prefix.name),
                          source_fingerprints=identity, files=copied,
                          verification='Full SHA-256 of streamed source and independent reread of copied file; original inputs retained')
            dump_json(manifest, record)
            return str(destination / prefix.name), record
        except BaseException:
            if staging.exists():
                shutil.rmtree(staging)
            raise


# ========================================================================
# Harmonize
# ========================================================================

# Local, reciprocal, single-nucleotide GRCh37/38 harmonization.
# 
# No rsID-only joins, network downloads, indel lift-over, or strand guessing.
# Reference changes between assemblies are retained as nucleotide identities.


BASES = frozenset("ACGT")
MAPPING_AUTOSOMES = frozenset(str(i) for i in range(1, 23))
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def orient(base, strand):
    if base not in BASES or strand not in ("+", "-"):
        raise ValueError("orient requires one canonical base and + or - strand")
    return base if strand == "+" else base.translate(COMPLEMENT)


class MappingError(ValueError):
    def __init__(self, reason, detail=""):
        self.reason = reason
        super().__init__(reason + (": " + detail if detail else ""))


class GenomeMapper:
    """Map a validated source base to both builds; all reported positions are 1-based.

    Both chains are deliberately local files. The caller supplies GRCh38 and
    GRCh37/hs37d5 FASTAs, including indexes. Reject all ambiguous/non-reciprocal
    mappings rather than accepting a highest-scoring chain by default.
    """
    def __init__(self, fasta38, fasta37, chain38to37, chain37to38):
        try:
            from pyliftover import LiftOver
        except ImportError as exc:
            raise ValueError("Ancient integration requires pyliftover; pip install pyliftover==0.4.1") from exc
        if not chain38to37 or not chain37to38:
            raise ValueError("Both local reciprocal chain files are required")
        self.f38 = FastaContext(fasta38)
        try:
            self.f37 = FastaContext(fasta37)
            self.to37 = LiftOver(str(chain38to37))
            self.to38 = LiftOver(str(chain37to38))
        except BaseException:
            self.f38.close()
            if hasattr(self, "f37"):
                self.f37.close()
            raise

    @staticmethod
    def _unique(chain, chrom, pos):
        hits = chain.convert_coordinate("chr" + canonical_chrom(chrom), pos - 1, "+") or []
        # Two chains with the same coordinate still constitute an ambiguous map.
        if not hits:
            raise MappingError("unmapped")
        if len(hits) != 1:
            raise MappingError("ambiguous_mapping")
        c, p, strand, _ = hits[0]
        c = canonical_chrom(c)
        if c not in MAPPING_AUTOSOMES:
            raise MappingError("non_autosomal_target")
        if strand not in ("+", "-"):
            raise MappingError("invalid_mapping_strand")
        return c, p + 1, strand

    def map_coordinates(self, chrom, pos, source_build="GRCh38"):
        """Return reciprocal coordinates and strand without reading either FASTA.

        This is an exact candidate prefilter, not allele harmonization. A result
        must subsequently pass ``map_from_coordinates`` with the modern source
        REF before it can contribute to an intersection. Reversing the direction
        of this reciprocal, unique mapping preserves the accepted coordinate set.
        """
        chrom = canonical_chrom(chrom)
        if chrom not in MAPPING_AUTOSOMES:
            raise MappingError("non_autosomal_source")
        if source_build not in ("GRCh37", "GRCh38"):
            raise ValueError("source_build must be GRCh37 or GRCh38")
        if not isinstance(pos, int) or pos < 1:
            raise MappingError("invalid_source_position")
        forward, backward = ((self.to37, self.to38) if source_build == "GRCh38"
                             else (self.to38, self.to37))
        target_c, target_p, strand = self._unique(forward, chrom, pos)
        back_c, back_p, back_strand = self._unique(backward, target_c, target_p)
        if (back_c, back_p) != (chrom, pos) or back_strand != strand:
            raise MappingError("non_reciprocal")
        if source_build == "GRCh38":
            c38, p38, c37, p37 = chrom, pos, target_c, target_p
        else:
            c38, p38, c37, p37 = target_c, target_p, chrom, pos
        return dict(chrom38=c38, pos38=p38, chrom37=c37, pos37=p37,
                    strand37to38=strand)

    def _source_context(self, chrom, pos, ref, source_build):
        if chrom not in MAPPING_AUTOSOMES:
            raise MappingError("non_autosomal_source")
        if ref not in BASES:
            raise MappingError("non_snv_source")
        if source_build not in ("GRCh37", "GRCh38"):
            raise ValueError("source_build must be GRCh37 or GRCh38")
        source_fasta = self.f38 if source_build == "GRCh38" else self.f37
        try:
            context = source_fasta.context(chrom, pos)
        except (KeyError, ValueError) as exc:
            raise MappingError("source_fasta_unavailable", str(exc)) from exc
        if context[1] != ref:
            raise MappingError("source_ref_mismatch", f"{chrom}:{pos} {ref}!={context[1]}")
        return context

    def _with_reference(self, coordinates, ref, source_build, source_context):
        c38, p38 = coordinates["chrom38"], coordinates["pos38"]
        c37, p37 = coordinates["chrom37"], coordinates["pos37"]
        strand = coordinates["strand37to38"]
        to38 = "+" if source_build == "GRCh38" else strand
        try:
            context38 = (source_context if source_build == "GRCh38" else self.f38.context(c38, p38))
            context37 = (source_context if source_build == "GRCh37" else self.f37.context(c37, p37))
        except (KeyError, ValueError) as exc:
            raise MappingError("target_fasta_unavailable", str(exc)) from exc
        ref38, ref37 = context38[1], context37[1]
        if ref38 not in BASES or ref37 not in BASES:
            raise MappingError("noncanonical_reference")
        canonical_context = (context38 if ref38 in "CT" else
                             context38.translate(COMPLEMENT)[::-1])
        return dict(chrom38=c38, pos38=p38, ref38=ref38, chrom37=c37,
                    pos37=p37, ref37=ref37, strand37to38=strand,
                    source_to38=to38, source_ref38=orient(ref, to38),
                    context38=context38, canonical_context=canonical_context,
                    reference_changed=int(orient(ref37, strand) != ref38),
                    site_id=f"{c38}:{p38}:{ref38}")

    def map_from_coordinates(self, coordinates, ref, source_build="GRCh38"):
        """Validate both FASTAs using coordinates returned by this mapper.

        The caller must match the modern coordinate to the appropriate build in
        ``coordinates`` exactly. This method deliberately makes no chain queries;
        it is only for reusing the reciprocal prefilter from ``map_coordinates``.
        """
        if source_build not in ("GRCh37", "GRCh38"):
            raise ValueError("source_build must be GRCh37 or GRCh38")
        suffix = "38" if source_build == "GRCh38" else "37"
        source_context = self._source_context(coordinates["chrom" + suffix],
                                              coordinates["pos" + suffix], ref, source_build)
        return self._with_reference(coordinates, ref, source_build, source_context)

    def map(self, chrom, pos, ref, source_build="GRCh38"):
        # Keep source-REF validation first for the reference modern-first path.
        chrom = canonical_chrom(chrom)
        source_context = self._source_context(chrom, pos, ref, source_build)
        coordinates = self.map_coordinates(chrom, pos, source_build)
        return self._with_reference(coordinates, ref, source_build, source_context)

    def close(self):
        self.f38.close()
        self.f37.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()



def load_module(filename, module_name):
    """Load a numbered entry filename under a stable Python module name."""
    import importlib.util
    import sys
    if module_name in sys.modules:
        return sys.modules[module_name]
    spec = importlib.util.spec_from_file_location(module_name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    return module
