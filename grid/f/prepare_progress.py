"""Crash-safe dosage progress. The caller must hold the cache preparation lock."""
from pathlib import Path
import json
import os
import numpy as np
from core2 import atomic_json


def write_matrix(cache, signature, shape, blocks, release_pages, *, commit_columns=None):
    """Resume only durably committed columns with identical inputs and ordering.

    blocks(start) must yield variant-by-sample arrays beginning at start.
    A failure may leave an uncommitted slice; the next call overwrites it.
    """
    cache = Path(cache)
    partial = cache / "genotypes.inprogress.npy"
    target = cache / "genotypes.npy"
    checkpoint = cache / "dosage-progress.json"
    shape = tuple(map(int, shape))
    identity = {"schema": "grid-dosage-progress-1", "signature": signature,
                "shape": list(shape), "dtype": "float16"}
    state = None
    if checkpoint.is_file():
        state = json.loads(checkpoint.read_text())
        if any(state.get(k) != v for k, v in identity.items()):
            raise ValueError("Dosage progress belongs to different inputs; use a new cache")
        start = state.get("columns_completed")
        if type(start) is not int or not 0 <= start <= shape[1]:
            raise ValueError("Invalid dosage progress offset")
        source = partial if partial.exists() else target
        if not source.is_file() or source.is_symlink():
            raise ValueError("Dosage progress file is missing or unsafe")
        stat = source.stat()
        if [stat.st_dev, stat.st_ino] != state.get("file_identity") or stat.st_size != state.get("file_size"):
            raise ValueError("Dosage progress file was replaced or truncated")
        matrix = np.lib.format.open_memmap(source, mode="r+")
        if matrix.shape != shape or matrix.dtype != np.float16 or not matrix.flags.c_contiguous:
            matrix._mmap.close()
            raise ValueError("Dosage progress matrix has an incompatible header")
        if source == target:
            matrix._mmap.close()
            if start != shape[1]:
                raise ValueError("Published dosage matrix has incomplete progress")
            return target
        print(f"PRSformer: resume dosage export at {start:,}/{shape[1]:,} SNPs", flush=True)
    else:
        if partial.exists() or target.exists():
            raise ValueError("Unrecorded dosage file exists; use a new cache")
        matrix = np.lib.format.open_memmap(partial, mode="w+", dtype=np.float16, shape=shape, fortran_order=False)
        start = 0
    capacity = min(shape[1], int(commit_columns)) if commit_columns is not None else None
    if capacity is not None and capacity < 1:
        raise ValueError("commit_columns must be positive")
    offset = start
    try:
        # Publish a durable zero-column record before any writes, too.
        if state is None:
            release_pages(matrix)
            with partial.open("rb") as handle: os.fsync(handle.fileno())
            stat = partial.stat()
            state = {**identity, "columns_completed": 0,
                     "file_identity": [stat.st_dev, stat.st_ino], "file_size": stat.st_size}
            atomic_json(state, checkpoint)
        def buffered_blocks():
            if capacity is None:
                yield from blocks(start)
                return
            # Sample-major RAM tile: one cast/copy, then wide contiguous row writes.
            buffer = np.empty((shape[0], capacity), dtype=np.float16)
            filled = 0
            seen = start
            for values in blocks(start):
                values = np.asarray(values)
                if values.ndim != 2 or values.shape[1] != shape[0] or not len(values) or seen + len(values) > shape[1]:
                    raise ValueError("Invalid dosage export block")
                seen += len(values)
                pos = 0
                while pos < len(values):
                    take = min(capacity - filled, len(values) - pos)
                    buffer[:, filled:filled + take] = values[pos:pos + take].T
                    filled += take
                    pos += take
                    if filled == capacity:
                        yield buffer.T
                        filled = 0
            if filled:
                yield buffer[:, :filled].T

        for values in buffered_blocks():
            values = np.asarray(values)
            if values.ndim != 2 or values.shape[1] != shape[0] or not len(values) or offset + len(values) > shape[1]:
                raise ValueError("Invalid dosage export block")
            matrix[:, offset:offset + len(values)] = values.T
            release_pages(matrix)
            with partial.open("rb") as handle: os.fsync(handle.fileno())
            offset += len(values)
            atomic_json({**state, "columns_completed": offset}, checkpoint)
            print(f"PRSformer dosage export: {offset:,}/{shape[1]:,} SNPs committed", flush=True)
        if offset != shape[1]:
            raise RuntimeError("Incomplete dosage export; committed progress retained")
    finally:
        matrix._mmap.close()
    partial.replace(target)
    # Retain the checkpoint until all companion tables and the receipt commit.
    return target
