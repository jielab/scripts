"""Ordered process parallelism with one shared CPU and host-memory budget."""
from __future__ import annotations

import math
import os


def worker_count(threads, tasks, *, memory_per_worker_gib=.5, reserve_gib=0.):
    if int(threads) != threads or threads < 1:
        raise ValueError("threads must be a positive integer")
    count = min(int(threads), max(1, int(tasks)))
    cap = os.environ.get("GRID_MEMORY_CAP_GB")
    if cap:
        available = max(0., float(cap)*.8-float(reserve_gib))
        count = min(count, max(1, int(available/max(.01, memory_per_worker_gib))))
    return max(1, count)


def ordered_map(function, arguments, workers, *, inner_threads=1):
    """Yield in input order; process exceptions propagate before publication.

    Arrays above 1 MiB are shared through read-only joblib memmaps. Dispatch is
    bounded; inner BLAS/OpenMP pools cannot each claim the whole allocation.
    The caller retains all random-number generation and deterministic ordering.
    """
    if workers == 1:
        for args in arguments:
            yield function(*args)
        return
    from joblib import Parallel, delayed, parallel_config
    with parallel_config(backend="loky", inner_max_num_threads=int(inner_threads)):
        with Parallel(n_jobs=int(workers), return_as="generator", pre_dispatch=workers,
                      batch_size=1, max_nbytes="1M", mmap_mode="r") as pool:
            yield from pool(delayed(function)(*args) for args in arguments)


def row_slices(length, workers):
    width = max(1, math.ceil(length/max(1, workers)))
    return [(start, min(length, start+width)) for start in range(0, length, width)]
