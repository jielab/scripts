#!/usr/bin/env python3
"""One-pass, ID-checked Kraken2/FASTQ demultiplexer. No all-read ID dictionary."""
from __future__ import annotations
from collections import Counter, OrderedDict
from pathlib import Path
import gzip
import re

ASSIGNMENT_FIELDS = ["ordinal", "read_id", "assigned_taxid", "group", "species_taxid", "bucket",
    "canonical_taxid", "base_group", "virus_genome_class", "virus_genome_type",
    "virus_resolution_status", "evidence_taxid"]


def fastq(path: Path):
    """ENA generated four-line FASTQ. Wrapped/malformed/truncated records fail closed."""
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rb') as f:
        n = 0
        while True:
            header = f.readline()
            if not header:
                break
            n += 1
            seq, plus, qual = f.readline(), f.readline(), f.readline()
            if (not header.startswith(b'@') or not header[1:].strip() or not plus.startswith(b'+') or not qual
                    or len(seq.rstrip(b'\r\n')) != len(qual.rstrip(b'\r\n'))
                    or not seq.rstrip(b'\r\n')):
                raise RuntimeError(f'Invalid/truncated/non-four-line FASTQ: {path}, record {n}')
            yield (header, seq, plus, qual)


def rid(header: bytes) -> bytes:
    token = header.strip().split()[0]
    if token.startswith(b'@'):
        token = token[1:]
    if token.endswith((b'/1', b'/2')):
        token = token[:-2]
    return token


def pair_ids(a: bytes, b: bytes) -> bytes:
    for header, expected in ((a, b'1'), (b, b'2')):
        token = header.strip().split()[0]
        if token.endswith((b'/1', b'/2')) and token[-1:] != expected:
            raise RuntimeError(f'Mate suffix contradicts file role: {token!r}')
        parts = header.strip().split()
        if len(parts)>1 and parts[1].startswith((b'1:', b'2:')) and parts[1][:1] != expected:
            raise RuntimeError('CASAVA mate suffix contradicts file role')
    x, y = rid(a), rid(b)
    if x == y:
        return x
    # Some SRA-derived names explicitly suffix the mate after the spot number.
    ma = re.fullmatch(rb'([SED]RR[0-9]+\.[0-9]+)\.1', x)
    mb = re.fullmatch(rb'([SED]RR[0-9]+\.[0-9]+)\.2', y)
    if ma and mb and ma[1] == mb[1]:
        return ma[1]
    raise RuntimeError(f'Mate IDs do not match: {x!r}, {y!r}')


class GzipBins:
    """Bounded RAM and file handles; standard concatenated gzip members.

    Flush >=256 KiB batches instead of spawning pigz once per species or opening
    thousands of persistent gzip writers. gzip/pigz/seqkit can read all members.
    """
    def __init__(self, root: Path, memory_mb: int = 64, max_open: int = 32):
        self.root, self.limit, self.max_open = root, memory_mb << 20, max_open
        self.buffers = OrderedDict()
        self.handles = OrderedDict()
        self.buffered = 0
        self.created = set()

    def _flush(self, key):
        buf = self.buffers.pop(key)
        self.buffered -= len(buf)
        if key in self.handles:
            handle = self.handles.pop(key)
        else:
            path = self.root / key
            path.parent.mkdir(parents=True, exist_ok=True)
            # The whole output root is a NEW staging directory; final files never append.
            mode = 'ab' if key in self.created else 'xb'
            handle = path.open(mode)
            self.created.add(key)
        self.handles[key] = handle
        handle.write(gzip.compress(buf, compresslevel=1, mtime=0))
        while len(self.handles) > self.max_open:
            _, h = self.handles.popitem(last=False)
            h.close()

    def write(self, key: str, record):
        data = b''.join(record)
        buf = self.buffers.pop(key, bytearray())
        buf.extend(data)
        self.buffers[key] = buf
        self.buffered += len(data)
        if len(buf) >= 262144:
            self._flush(key)
        while self.buffered > self.limit:
            self._flush(next(iter(self.buffers)))

    def close(self):
        for key in list(self.buffers):
            self._flush(key)
        for h in self.handles.values():
            h.close()
        self.handles.clear()


def split_stream(lines, r1: Path, r2: Path | None, stream: str, tax, bins: GzipBins,
                 assignments, full_kraken=None, virus_counter=None) -> tuple[Counter, dict]:
    """Kraken2 ordered output is checked against EVERY ID and length; never blind zip.

    PE assignments apply to the fragment, not independent certification of each mate.
    """
    a, b = fastq(r1), (fastq(r2) if r2 else None)
    counter = Counter()
    totals = {'fragments': 0, 'reads': 0, 'bases': 0}
    for number, line in enumerate(lines, 1):
        if isinstance(line, str):
            line = line.encode()
        p = line.rstrip(b'\r\n').split(b'\t')
        if len(p) != 5 or p[0] not in (b'C', b'U'):
            raise RuntimeError(f'Unexpected Kraken2 output line {number}: {line[:100]!r}')
        rec1, rec2 = next(a, None), (next(b, None) if b is not None else None)
        if rec1 is None or (b is not None and rec2 is None):
            raise RuntimeError('More Kraken2 assignments than FASTQ fragments')
        key = pair_ids(rec1[0], rec2[0]) if rec2 else rid(rec1[0])
        if rid(p[1]) not in (key, rid(rec1[0])):
            raise RuntimeError(f'Kraken2/FASTQ ID/order mismatch at {number}; STOP, never silently join by position.')
        lengths = [len(rec1[1].rstrip(b'\r\n'))]
        if rec2:
            lengths.append(len(rec2[1].rstrip(b'\r\n')))
        observed = [int(x) for x in p[3].split(b'|')]
        if observed != lengths:
            raise RuntimeError(f'Kraken2/FASTQ length mismatch at {number}: {observed} != {lengths}')
        tid = int(p[2])
        if (p[0] == b'U') != (tid == 0):
            raise RuntimeError('Inconsistent classified flag / taxid')
        group, bucket, species, name = tax.route(tid)
        ann = tax.annotation(tid)
        if virus_counter is not None and ann['base_group'] == 'virus':
            virus_counter[(tid, stream, len(lengths))] += 1
        if rec2:
            bins.write(f'{group}/{bucket}/{stream}.R1.fastq.gz', rec1)
            bins.write(f'{group}/{bucket}/{stream}.R2.fastq.gz', rec2)
        else:
            bins.write(f'{group}/{bucket}/{stream}.fastq.gz', rec1)
        counter[(group, bucket, species, name, stream, len(lengths))] += 1
        assignments.write(b'\t'.join([str(number).encode(), key, p[2], group.encode(),
                                    str(species).encode(), bucket.encode(),
                                    *[str(ann[k]).encode() for k in ('canonical_taxid', 'base_group',
                                      'genome_class', 'genome_type_detail', 'resolution_status', 'evidence_taxid')]]) + b'\n')
        if full_kraken is not None:
            full_kraken.write(line)
        totals['fragments'] += 1
        totals['reads'] += len(lengths)
        totals['bases'] += sum(lengths)
    if next(a, None) is not None or (b is not None and next(b, None) is not None):
        raise RuntimeError('Fewer Kraken2 assignments than FASTQ fragments')
    return counter, totals
