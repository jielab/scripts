"""Strict, bounded FASTQ streaming; records retain all original bytes."""
from contextlib import ExitStack
import gzip
import hashlib
import io
from pathlib import Path
import re
import struct
import zlib

_U64 = struct.Struct("<Q")
_SRA_MATE = re.compile(rb"^([SED]RR\d+\.\d+)\.([12])$")
_CASAVA_MATE = re.compile(rb"^([12]):[YN]:")


def read_id(header):
    if not header.startswith(b"@"):
        raise RuntimeError("FASTQ header must begin with @")
    fields = header[1:].rstrip(b"\r\n").split(None, 1)
    if not fields or not fields[0]:
        raise RuntimeError("FASTQ read identifier is empty")
    return fields[0]


def validate_pair_ids(first, second):
    ids = (read_id(first), read_id(second))
    for header, name, mate in zip((first, second), ids, (b"1", b"2")):
        if name.endswith((b"/1", b"/2")) and name[-1:] != mate:
            raise RuntimeError(f"Reversed/duplicated FASTQ mate label: {ids!r}")
        sra = _SRA_MATE.fullmatch(name)
        if sra and sra[2] != mate:
            raise RuntimeError(f"Reversed/duplicated SRA mate label: {ids!r}")
        fields = header[1:].rstrip(b"\r\n").split(None, 1)
        casava = _CASAVA_MATE.match(fields[1]) if len(fields) > 1 else None
        if casava and casava[1] != mate:
            raise RuntimeError(f"Reversed/duplicated CASAVA mate label: {ids!r}")
    a, b = (name[:-2] if name.endswith((b"/1", b"/2")) else name for name in ids)
    if a == b:
        return
    ma, mb = _SRA_MATE.fullmatch(a), _SRA_MATE.fullmatch(b)
    if ma and mb and ma[1] == mb[1] and ma[2] == b"1" and mb[2] == b"2":
        return
    raise RuntimeError(f"FASTQ paired identifiers differ: {ids!r}")


def fastq_records(path):
    """Four-line FASTQ; a complete scan also validates every gzip trailer.

    CRLF, plus-line comments and a missing final quality newline are preserved.
    Wrapped FASTQ is rejected rather than silently changing record boundaries.
    """
    path = Path(path)
    with ExitStack() as stack:
        raw = stack.enter_context(path.open("rb", buffering=4 << 20))
        if path.name.lower().endswith(".gz"):
            compressed = stack.enter_context(gzip.GzipFile(fileobj=raw, mode="rb"))
            handle = stack.enter_context(io.BufferedReader(compressed, 1 << 20))
        else:
            handle = raw
        number = 0
        try:
            while True:
                header = handle.readline()
                if not header:
                    return
                sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
                number += 1
                if not sequence or not plus or not quality or not plus.startswith(b"+"):
                    raise RuntimeError(f"Invalid/truncated FASTQ {path}, record {number}")
                read_id(header)
                bases, scores = sequence.rstrip(b"\r\n"), quality.rstrip(b"\r\n")
                if not bases or len(bases) != len(scores):
                    raise RuntimeError(f"Sequence/quality length mismatch in {path}, record {number}")
                yield header, sequence, plus, quality
        except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
            raise RuntimeError(f"Damaged gzip FASTQ {path}, near record {number}: {exc}") from exc


def fragments(r1, r2=None, limit=0):
    if isinstance(limit, bool) or not isinstance(limit, int) or limit < 0:
        raise ValueError("limit must be a nonnegative integer")
    with ExitStack() as stack:
        left = fastq_records(r1)
        right = fastq_records(r2) if r2 is not None else None
        stack.callback(left.close)
        if right is not None:
            stack.callback(right.close)
        count = 0
        while not limit or count < limit:
            a = next(left, None)
            b = next(right, None) if right is not None else None
            if a is None:
                if b is not None:
                    raise RuntimeError("FASTQ mate counts differ: R2 is longer")
                return
            if right is not None:
                if b is None:
                    raise RuntimeError("FASTQ mate counts differ: R1 is longer")
                validate_pair_ids(a[0], b[0])
            count += 1
            yield a, b


def new_digest():
    return hashlib.sha256()


def update_digest(digest, a, b=None):
    """Bind ordinal order to every original header, sequence, plus and quality byte."""
    for record in (a, b) if b is not None else (a,):
        for line in record:
            digest.update(_U64.pack(len(line)))
            digest.update(line)
