#!/usr/bin/env python3
"""Bounded genotype staging and sequential publication without shared mmap writes."""
from __future__ import annotations

from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import threading
import time

import numpy as np

GIB = 1024 ** 3
LOCAL_FILESYSTEMS = {"xfs", "ext2", "ext3", "ext4", "btrfs", "zfs", "overlay"}


def memory_headroom():
    """Available RAM bounded by the job budget and Linux cgroup, never host RAM alone."""
    info = {k: int(v.split()[0]) * 1024 for k, v in
            (line.split(":", 1) for line in Path('/proc/meminfo').read_text().splitlines())}
    status = Path('/proc/self/status').read_text().splitlines()
    rss = next((int(s.split()[1]) * 1024 for s in status if s.startswith('VmRSS:')), 0)
    limits = [info.get('MemAvailable', info['MemFree'])]
    budget = os.environ.get('GRID_MEMORY_CAP_GB')
    if budget:
        value = float(budget)
        if not np.isfinite(value) or value <= 0:
            raise ValueError('GRID_MEMORY_CAP_GB must be positive and finite')
        limits.append(max(0, int(value * GIB) - rss))
    elif os.environ.get('LSB_JOBID'):
        # An HPC process without a known allocation must not allocate host-wide RAM.
        limits.append(0)
    # Resolve this process's nested cgroup and every ancestor limit.
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        _, controllers, group = line.split(':', 2)
        if controllers == '':
            base, names = Path('/sys/fs/cgroup'), ('memory.max', 'memory.current')
        elif 'memory' in controllers.split(','):
            base, names = Path('/sys/fs/cgroup/memory'), ('memory.limit_in_bytes', 'memory.usage_in_bytes')
        else:
            continue
        current = base / group.lstrip('/')
        if '..' in current.parts:
            current = base
        while current == base or base in current.parents:
            try:
                ceiling = (current / names[0]).read_text().strip()
                if ceiling != 'max':
                    used = int((current / names[1]).read_text())
                    limits.append(max(0, int(ceiling) - used))
            except FileNotFoundError:
                pass
            if current == base:
                break
            current = current.parent
    return min(limits)


def local_staging_directory(required):
    """Never silently put the transpose staging file back onto GPFS/NFS/tmpfs."""
    explicit = os.environ.get('GRID_PREPARE_TMPDIR')
    candidates = [explicit] if explicit else [os.environ.get('TMPDIR'), '/tmp']
    reasons = []
    for name in dict.fromkeys(x for x in candidates if x):
        folder = Path(name).expanduser().resolve()
        if not folder.is_dir():
            reasons.append(f'{folder}: missing directory')
            continue
        result = subprocess.run(['findmnt', '-n', '-T', str(folder), '-o', 'FSTYPE'],
                                text=True, capture_output=True, check=True)
        kind = result.stdout.strip()
        if kind not in LOCAL_FILESYSTEMS:
            reasons.append(f'{folder}: {kind} is not node-local disk')
            continue
        if shutil.disk_usage(folder).free < required + max(GIB, required // 20):
            reasons.append(f'{folder}: insufficient free disk')
            continue
        return folder
    raise ValueError('No safe local genotype staging space: ' + '; '.join(reasons)
                     + '. Set GRID_PREPARE_TMPDIR to node-local disk or request more RAM.')


@contextmanager
def variant_buffer(shape, log=print):
    size = int(shape[0]) * int(shape[1]) * np.dtype('float16').itemsize
    mode = os.environ.get('GRID_PREPARE_STORAGE', 'auto').lower()
    if mode not in ('auto', 'ram', 'local'):
        raise ValueError('GRID_PREPARE_STORAGE must be auto, ram or local')
    available = memory_headroom()
    # Keep at least 25% headroom (and 1 GiB) for decoded blocks, cohort and libraries.
    fits = size + max(GIB, available // 4) <= available
    if mode == 'ram' and not fits:
        raise ValueError(f'RAM staging needs {size/GIB:.2f} GiB plus headroom; only {available/GIB:.2f} GiB available')
    if mode != 'local' and fits:
        log(f'STAGING: RAM {size/GIB:.2f} GiB; job headroom {available/GIB:.2f} GiB')
        values = np.empty(shape, dtype=np.float16)
        try:
            yield values
        finally:
            del values
        return
    folder = local_staging_directory(size)
    with tempfile.TemporaryDirectory(prefix='grid-genotype-', dir=folder) as name:
        work = Path(name)
        (work / 'OWNER.json').write_text(json.dumps(dict(pid=os.getpid(), job=os.getenv('LSB_JOBID'))))
        log(f'STAGING: node-local disk {size/GIB:.2f} GiB at {work}')
        values = np.memmap(work / 'variants.bin', mode='w+', dtype=np.float16, shape=shape)
        try:
            yield values
        finally:
            values._mmap.close()


def _sync(stream, done, total, log):
    # Give a visible heartbeat even if the filesystem stalls inside fsync.
    started = time.monotonic()
    stopped = threading.Event()
    def report():
        while not stopped.wait(60):
            log(f'WRITE: waiting for filesystem sync at {done:,}/{total:,} rows ({time.monotonic()-started:.0f}s)')
    worker = threading.Thread(target=report, daemon=True)
    worker.start()
    try:
        stream.flush()
        os.fsync(stream.fileno())
    finally:
        stopped.set()
        worker.join()


def write_sample_major(path, raw, variants, *, log=print, buffer_bytes=64*1024**2,
                       sync_bytes=GIB):
    """Write a C-order float16 .npy in large sequential writes, with bounded dirty data.

    The genotype/QC values and ordering are unchanged. No writable mapping of
    the shared destination is created, and publication waits for durable data.
    """
    samples = raw.shape[1]
    if not 0 < variants <= raw.shape[0] or samples <= 0:
        raise ValueError('Invalid genotype matrix dimensions')
    batch = max(1, buffer_bytes // (2 * variants))
    started = last = time.monotonic()
    since_sync = 0
    log(f'WRITE: sequential sample-major genotypes {samples:,} x {variants:,}; {(samples*variants*2)/GIB:.2f} GiB')
    with open(path, 'wb', buffering=0) as stream:
        np.lib.format.write_array_header_1_0(stream, dict(descr=np.dtype('float16').str,
                                                        fortran_order=False, shape=(samples, variants)))
        for start in range(0, samples, batch):
            end = min(start + batch, samples)
            block = np.ascontiguousarray(raw[:variants, start:end].T, dtype=np.float16)
            data = memoryview(block).cast('B')
            while data:
                written = stream.write(data)
                if not written:
                    raise OSError('No progress writing genotype cache')
                data = data[written:]
            since_sync += block.nbytes
            if since_sync >= sync_bytes:
                _sync(stream, end, samples, log)
                since_sync = 0
            now = time.monotonic()
            if now - last >= 60 or end == samples:
                rate = end * variants * 2 / max(now-started, .001) / 1024**2
                log(f'WRITE: {end:,}/{samples:,} rows ({end/samples:.1%}); {rate:.1f} MiB/s')
                last = now
        _sync(stream, samples, samples, log)
        expected = stream.tell()
    if Path(path).stat().st_size != expected:
        raise ValueError('Truncated genotype output')
    log(f'WRITE OK: durable genotype matrix in {time.monotonic()-started:.1f}s')
