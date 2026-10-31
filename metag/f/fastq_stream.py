"""Bounded-memory FASTQ validation and optional query-only sequence staging.

Original FASTQ files are read-only.  A prepared directory is published only after
all selected records, gzip trailers (for full scans), and paired IDs validate.
"""
from __future__ import annotations

from contextlib import ExitStack
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import struct
import tempfile
import time
import zlib

from common import atomic_json

SCHEMA = "metag-prepared-reads-v1"
FNV_OFFSET = 14695981039346656037
FNV_PRIME = 1099511628211
UINT64_MASK = (1 << 64) - 1
_U64 = struct.Struct("<Q")
_SRA_MATE = re.compile(rb"^([SED]RR\d+\.\d+)\.([12])$")


def input_state(path):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise RuntimeError(f"FASTQ input is not a regular file: {path}")
    st = path.stat()
    return {"path": str(path), "bytes": st.st_size, "mtime_ns": st.st_mtime_ns,
            "ctime_ns": st.st_ctime_ns, "device": st.st_dev, "inode": st.st_ino}


def native_query_header(header):
    """Bytes printed by pinned MetaCache: R1 header, sans @, before first space.

Mate suffixes are retained.  The compact merge validates this exact header via a
per-ordinal FNV-1a fingerprint, rather than relying on globally unique read names.
    """
    value = header.rstrip(b"\r\n")
    if not value.startswith(b"@"):
        raise RuntimeError("FASTQ header must start with @")
    token = value[1:].split(b" ", 1)[0]
    if not token or any(c in token for c in (b"\t", b"|", b"\r", b"\n")):
        raise RuntimeError("FASTQ read ID is empty or contains a tab/pipe unsupported by MetaCache mapping output")
    return token


def header_fingerprint(header):
    value = FNV_OFFSET
    for char in native_query_header(header):
        value = ((value ^ char) * FNV_PRIME) & UINT64_MASK
    return value


def new_header_digest():
    return hashlib.sha256()


def update_header_digest(digest, rec1, rec2=None):
    """Hash both exact raw headers, each prefixed by its uint64-LE byte length."""
    for record in (rec1, rec2) if rec2 is not None else (rec1,):
        digest.update(_U64.pack(len(record[0])))
        digest.update(record[0])


def validate_pair_ids(header1, header2):
    first, second = native_query_header(header1), native_query_header(header2)
    # Explicit mate labels cannot both name mate 1, nor be swapped.  Plain
    # identical identifiers remain valid (e.g. CASAVA mate labels after a space).
    for token, mate in ((first, b"1"), (second, b"2")):
        if token.endswith((b"/1", b"/2")) and token[-1:] != mate:
            raise RuntimeError(f"Paired FASTQ mate label is reversed or duplicated: {first!r} vs {second!r}")
        sra = _SRA_MATE.fullmatch(token)
        if sra and sra[2] != mate:
            raise RuntimeError(f"Paired FASTQ SRA mate label is reversed or duplicated: {first!r} vs {second!r}")
    a = first[:-2] if first.endswith((b"/1", b"/2")) else first
    b = second[:-2] if second.endswith((b"/1", b"/2")) else second
    if a == b:
        return
    ma, mb = _SRA_MATE.fullmatch(a), _SRA_MATE.fullmatch(b)
    if ma and mb and ma[1] == mb[1] and ma[2] == b"1" and mb[2] == b"2":
        return
    raise RuntimeError(f"Paired FASTQ IDs differ: {first!r} vs {second!r}")


def fastq_records(path):
    """Read conventional four-line FASTQ with large, bounded I/O buffers.

CRLF and a missing final quality newline are preserved.  A full iteration
consumes gzip members/trailers, so damaged or incomplete gzip input fails.
    """
    path = Path(path)
    with ExitStack() as stack:
        raw = stack.enter_context(path.open("rb", buffering=8 << 20))
        if str(path).lower().endswith(".gz"):
            compressed = stack.enter_context(gzip.GzipFile(fileobj=raw, mode="rb"))
            handle = stack.enter_context(io.BufferedReader(compressed, buffer_size=1 << 20))
        else:
            handle = raw
        number = 0
        try:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
                number += 1
                if not sequence or not plus or not quality or not header.startswith(b"@") or not plus.startswith(b"+"):
                    raise RuntimeError(f"Invalid/truncated FASTQ {path} record {number}")
                bases = sequence.rstrip(b"\r\n")
                if not bases:
                    raise RuntimeError(f"Empty FASTQ sequence {path} record {number}")
                if len(bases) != len(quality.rstrip(b"\r\n")):
                    raise RuntimeError(f"Sequence/quality length mismatch {path} record {number}")
                native_query_header(header)
                yield header, sequence, plus, quality
        except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
            raise RuntimeError(f"Invalid/truncated gzip FASTQ {path} near record {number}: {exc}") from exc


def _disk_guard(path, reserve, required=0):
    free = shutil.disk_usage(path).free
    if free < reserve + required:
        raise RuntimeError(f"Insufficient free work space at {path}: {free / (1 << 30):.2f} GiB free, "
                           f"need {(reserve + required) / (1 << 30):.2f} GiB including reserve")


def _output_state(path):
    stat = path.stat()
    return {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _paths(receipt, directory):
    result = dict(receipt)
    result["header_fingerprints"] = str(directory / "candidate_headers.bin")
    result["candidate_headers"] = result["header_fingerprints"]
    for key in ("query_r1", "query_r2", "split_r1", "split_r2"):
        value = receipt[key]
        result[key] = str(directory / value) if value and not Path(value).is_absolute() else value
    result["receipt"] = str(directory / "prepared.json")
    return result


def prepare_reads(r1, r2, out, *, limit=0, cache_mode="auto", min_free_bytes=10 << 30,
                  max_cache_bytes=100 << 30):
    """Validate/count reads once, write ordinal header hashes, and optionally stage.

``limit=0`` validates the entire input.  A positive limit creates a deterministic
FASTQ prefix for a pilot; the remainder of the source is deliberately not read or
claimed valid.  ``cache_mode='fasta'`` writes uncompressed query-only FASTA files
to avoid repeated quality parsing/gzip work; species bins still use FASTQ via
``split_r1``/``split_r2``.  ``auto`` attempts that cache in the same pass, dropping
only its own partial cache if it exceeds the budget or disk reserve; a pilot uses
its bounded FASTQ prefix directly.  Explicit ``fasta`` fails instead of falling
back.  Output is a receipt dict with absolute usable paths.

Existing completed preparation is reused only when its inputs, options, and all
derived-file states match.  An incompatible directory is left intact and refused.
    """
    if not isinstance(limit, int) or isinstance(limit, bool) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")
    if cache_mode not in ("none", "fasta", "auto"):
        raise ValueError("cache_mode must be 'none', 'fasta', or 'auto'")
    if min_free_bytes < 0:
        raise ValueError("min_free_bytes must be nonnegative")
    if max_cache_bytes < 0:
        raise ValueError("max_cache_bytes must be nonnegative")
    r1 = Path(r1).resolve(strict=True)
    r2 = Path(r2).resolve(strict=True) if r2 is not None else None
    if r2 is not None and r1.samefile(r2):
        raise RuntimeError("Paired FASTQ inputs refer to the same file")
    inputs = [input_state(r1)] + ([input_state(r2)] if r2 else [])
    identity = {"schema": SCHEMA, "inputs": inputs, "limit": limit, "cache_mode": cache_mode,
                "max_cache_bytes": max_cache_bytes}
    out = Path(out).absolute()
    if out.exists():
        try:
            receipt = json.loads((out / "prepared.json").read_text())
            if receipt.get("identity") != identity or receipt.get("status") != "complete":
                raise RuntimeError("input/options mismatch")
            for name, state in receipt["outputs"].items():
                if Path(name).name != name or _output_state(out / name) != state:
                    raise RuntimeError(f"derived-file state mismatch: {name}")
            if receipt["outputs"]["candidate_headers.bin"]["bytes"] != receipt["fragments"] * 8:
                raise RuntimeError("header fingerprint count mismatch")
            print(f"Reuse validated reads: {receipt['fragments']:,} fragments", flush=True)
            return _paths(receipt, out)
        except (OSError, ValueError, KeyError, TypeError, RuntimeError) as exc:
            raise RuntimeError(f"Prepared reads are incomplete or incompatible: {out}; use a new work/name. {exc}") from exc
    out.parent.mkdir(parents=True, exist_ok=True)
    _disk_guard(out.parent, min_free_bytes)
    stage = Path(tempfile.mkdtemp(prefix=out.name + ".partial-", dir=out.parent))
    count = total_bases = total_bytes = 0
    digest = new_header_digest()
    start = last_report = time.monotonic()
    outputs = ["candidate_headers.bin"]
    split_names = ["prefix.R1.fastq", "prefix.R2.fastq"] if limit else [str(r1), str(r2) if r2 else None]
    if r2 is None:
        split_names[1] = None
    actual_cache_mode = "fasta" if cache_mode == "fasta" or (cache_mode == "auto" and not limit and max_cache_bytes) else "none"
    fallback_reason = "pilot_prefix" if cache_mode == "auto" and limit else ("cache_disabled_by_budget" if cache_mode == "auto" and not max_cache_bytes else "")
    cache_bytes = 0
    query_names = ["query.R1.fasta", "query.R2.fasta" if r2 else None] if actual_cache_mode == "fasta" else list(split_names)
    try:
        with ExitStack() as stack:
            fingerprints = stack.enter_context((stage / outputs[0]).open("wb", buffering=1 << 20))
            prefix_handles = []
            fasta_handles = []
            for name in split_names:
                if limit and name:
                    prefix_handles.append(stack.enter_context((stage / name).open("wb", buffering=1 << 20)))
                    outputs.append(name)
            for name in query_names:
                if actual_cache_mode == "fasta" and name:
                    fasta_handles.append(stack.enter_context((stage / name).open("wb", buffering=1 << 20)))
                    outputs.append(name)

            def drop_auto_cache(reason):
                nonlocal query_names, actual_cache_mode, fallback_reason
                if cache_mode != "auto":
                    raise RuntimeError(f"Requested FASTA cache cannot be completed: {reason}; increase its budget or use --read-cache none")
                for handle in fasta_handles:
                    handle.close()
                fasta_handles.clear()
                for name in query_names:
                    if name:
                        (stage / name).unlink()
                        outputs.remove(name)
                query_names = list(split_names)
                actual_cache_mode = "none"
                fallback_reason = reason
                print(f"Query cache disabled: {reason}; continuing input validation", flush=True)
            reads1, reads2 = fastq_records(r1), fastq_records(r2) if r2 else None
            stack.callback(reads1.close)
            if reads2:
                stack.callback(reads2.close)
            while not limit or count < limit:
                rec1 = next(reads1, None)
                rec2 = next(reads2, None) if reads2 else None
                if rec1 is None:
                    if rec2 is not None:
                        raise RuntimeError("FASTQ mates contain different numbers of records (R2 is longer)")
                    break
                if reads2 is not None and rec2 is None:
                    raise RuntimeError("FASTQ mates contain different numbers of records (R1 is longer)")
                if rec2 is not None:
                    validate_pair_ids(rec1[0], rec2[0])
                update_header_digest(digest, rec1, rec2)
                fingerprints.write(_U64.pack(header_fingerprint(rec1[0])))
                records = (rec1, rec2) if rec2 is not None else (rec1,)
                if fasta_handles:
                    extra = sum(len(record[0].rstrip(b"\r\n")) + len(record[1].rstrip(b"\r\n")) + 2 for record in records)
                    if cache_bytes + extra > max_cache_bytes:
                        drop_auto_cache("cache_budget_exceeded")
                for mate, record in enumerate(records):
                    total_bases += len(record[1].rstrip(b"\r\n"))
                    total_bytes += sum(map(len, record))
                    if prefix_handles:
                        prefix_handles[mate].write(b"".join(record))
                    if fasta_handles:
                        fasta_handles[mate].write(b">" + record[0][1:].rstrip(b"\r\n") + b"\n")
                        fasta_handles[mate].write(record[1].rstrip(b"\r\n") + b"\n")
                if fasta_handles:
                    cache_bytes += extra
                count += 1
                if count % 10000 == 0:
                    # Keep reserve beyond buffered bytes; check on actual output FS.
                    required = ((len(outputs) + 1) << 20) + count * 64
                    if fasta_handles and shutil.disk_usage(stage).free < min_free_bytes + required:
                        drop_auto_cache("low_disk_space")
                    _disk_guard(stage, min_free_bytes, required=required)
                    now = time.monotonic()
                    if now - last_report >= 30:
                        print(f"Validate reads: {count:,} fragments; {count / max(1, now-start):,.0f}/s", flush=True)
                        last_report = now
        if not count:
            raise RuntimeError("FASTQ input contains no fragments")
        if inputs != [input_state(r1)] + ([input_state(r2)] if r2 else []):
            raise RuntimeError("FASTQ input changed while it was being validated")
        if fasta_handles and shutil.disk_usage(stage).free < min_free_bytes + count * 64:
            drop_auto_cache("low_disk_space")
        _disk_guard(stage, min_free_bytes, required=count * 64)
        receipt = {"schema": SCHEMA, "status": "complete", "identity": identity,
                   "fragments": count, "total_reads": count * (2 if r2 else 1), "paired": bool(r2),
                   "total_bases": total_bases, "sequence_bases": total_bases,
                   "fastq_uncompressed_bytes": total_bytes, "input_bytes": sum(item["bytes"] for item in inputs),
                   "header_digest": digest.hexdigest(),
                   "header_digest_format": "SHA256: uint64-LE byte-length + raw header, R1 then R2 per fragment",
                   "header_fingerprint_format": "uint64-LE FNV1a64 of exact native first-mate query_header",
                   "validated_scope": "prefix" if limit else "entire_input", "limit": limit,
                   "query_limit": 0, "cache_mode": actual_cache_mode,
                   "requested_cache_mode": cache_mode, "cache_fallback_reason": fallback_reason,
                   "cache_bytes": cache_bytes if actual_cache_mode == "fasta" else 0,
                   "query_r1": query_names[0], "query_r2": query_names[1],
                   "split_r1": split_names[0], "split_r2": split_names[1],
                   "inputs": inputs, "elapsed_seconds": time.monotonic() - start,
                   "outputs": {name: _output_state(stage / name) for name in outputs}}
        atomic_json(stage / "prepared.json", receipt)
        stage.rename(out)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(f"Validated {count:,} fragments ({receipt['total_reads']:,} reads), query cache={actual_cache_mode}", flush=True)
    return _paths(receipt, out)
