"""Ordered process parallelism with one shared CPU and host-memory budget."""
from __future__ import annotations

import math
import os
from pathlib import Path

def available_memory():
    """Conservative host/cgroup headroom; zero disables RAM staging on unknown hosts."""
    try:
        fields=dict(line.split(":",1) for line in Path("/proc/meminfo").read_text().splitlines())
        available=int(fields["MemAvailable"].split()[0])*1024
        # Resolve this process's cgroup, including nested container limits.
        cgroup=next((line.split(":",2)[2] for line in Path("/proc/self/cgroup").read_text().splitlines()
                     if line.startswith("0::")), "")
        root=Path("/sys/fs/cgroup")
        current=root/cgroup.lstrip("/")
        while current.is_relative_to(root):
            limit=current/"memory.max"
            if limit.is_file() and limit.read_text().strip()!="max":
                available=min(available,max(0,int(limit.read_text())-int((current/"memory.current").read_text())))
            if current==root:break
            current=current.parent
        return available
    except (OSError,ValueError,KeyError):
        return 0



def worker_count(threads, tasks, *, memory_per_worker_gib=.5, reserve_gib=0.):
    if int(threads) != threads or threads < 1:
        raise ValueError("threads must be a positive integer")
    count = min(int(threads), max(1, int(tasks)))
    cap = os.environ.get("GRID_MEMORY_CAP_GB")
    headroom=available_memory()/2**30
    if cap or headroom:
        budget=min([v for v in (float(cap) if cap else 0,headroom) if v>0])
        available = max(0., budget*.8-float(reserve_gib))
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
