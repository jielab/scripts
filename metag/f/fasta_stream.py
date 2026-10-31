"""Bounded FASTA reader: never buffer an entire chromosome/record."""
import gzip
from contextlib import ExitStack

SKETCH_BATCH_BASES = 4 << 20


def chunks(path, block_size=4 << 20):
    """Yield None at each header, then DNA chunks; headers may span blocks."""
    if block_size < 2:
        raise ValueError("FASTA block size must be >= 2")
    with ExitStack() as stack:
        # WSL-mounted drives make gzip's default 8 KiB source reads expensive.
        raw = stack.enter_context(open(path, "rb", buffering=8 << 20))
        handle = stack.enter_context(gzip.GzipFile(fileobj=raw)) if str(path).lower().endswith(".gz") else raw
        if handle.read(1) != b">":
            raise ValueError(f"Expected FASTA header: {path}")
        yield None
        in_header = True
        tail = b""
        while True:
            block = handle.read(block_size)
            data = tail + block
            if block:
                data, tail = data[:-1], data[-1:]
            else:
                tail = b""
            if block and data.endswith(b"\n"):
                data, tail = data[:-1], b"\n" + tail
            for i, section in enumerate(data.split(b"\n>")):
                if i:
                    in_header = True
                    yield None
                if in_header:
                    newline = section.find(b"\n")
                    if newline < 0:
                        continue
                    section = section[newline+1:]
                    in_header = False
                dna = section.translate(None, b"\r\n \t")
                if b">" in dna:
                    raise ValueError(f"Malformed FASTA header inside a sequence: {path}")
                if dna:
                    yield dna
            if not block:
                if in_header:
                    raise ValueError(f"FASTA ends in an unfinished/empty header: {path}")
                break


def sketch(path, ksize, scaled):
    from sourmash import MinHash
    mh = MinHash(n=0, ksize=ksize, scaled=scaled, track_abundance=False, seed=42)
    pending = mh.copy_and_clear()
    pending_bases = 0
    overlap = b""
    records = bases = 0
    for sequence in chunks(path):
        if sequence is None:
            overlap = b""
            records += 1
            continue
        bases += len(sequence)
        joined = overlap + sequence
        # Insert into a small sketch, then merge sorted hashes in bulk. Adding
        # each new hash directly to a genome-sized sorted vector becomes costly
        # for large genomes. Accumulate across short records to bound merge work.
        pending.add_sequence(joined, force=True)
        pending_bases += len(sequence)
        if pending_bases >= SKETCH_BATCH_BASES:
            mh.merge(pending)
            pending = mh.copy_and_clear()
            pending_bases = 0
        overlap = joined[-(ksize-1):] if ksize > 1 else b""
    if not bases:
        raise ValueError(f"No sequence bases in {path}")
    if pending_bases:
        mh.merge(pending)
    return mh, records, bases
