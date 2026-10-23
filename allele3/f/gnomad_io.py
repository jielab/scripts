"""Bounded-memory local I/O for gnomAD sites VCFs.

VCF positions are one-based.  No record filtering, allele merging, inference of
reference counts, or recalculation of INFO counts occurs in this module.
``bcftools query`` is used when available; a standard-library VCF/gzip reader is
provided for portability.  Neither backend writes a decompressed intermediate.
"""
from __future__ import annotations

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
