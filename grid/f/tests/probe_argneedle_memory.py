"""Bounded real-input probe; intentionally stops after a few threaded samples.

Uses the full HAPS and the requested scaffold size, so it exercises the original
large allocation. Never writes an inference checkpoint or claims a full build.
Run inside the public 32 GiB / 4 GiB resource guard.
"""
import argparse
import logging
from pathlib import Path
import random
import sys

import numpy as np
import psutil

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from argneedle_memory import install

p = argparse.ArgumentParser()
p.add_argument('--root', required=True)
p.add_argument('--map', required=True)
p.add_argument('--samples', type=int, default=4000)
p.add_argument('--probe-samples', type=int, default=3)
p.add_argument('--probe-windows', type=int, default=0,
               help='Stress representative windows across the full chromosome; no sample threading completion')
p.add_argument('--mode', default='sequence')
a = p.parse_args()
install()
from arg_needle import inference
from arg_needle.decoders import ASMCDecoder

class ProbeComplete(Exception):
    pass

original = ASMCDecoder.compute_with_hashing
class RepresentativeHasher:
    def __init__(self, native):
        self.native = native
        self.windows = []
    def __getattr__(self, name):
        return getattr(self.native, name)
    def get_closest_cousins(self, *args):
        windows = self.native.get_closest_cousins(*args)
        selected = np.unique(np.linspace(0, len(windows) - 1,
                                        min(a.probe_windows, len(windows)), dtype=int))
        self.windows = [windows[j] for j in selected]
        return self.windows

def measured(self, i, *args, **kwargs):
    if a.probe_windows:
        self.hasher = RepresentativeHasher(self.hasher)
        if self.backup_hasher is not None:
            self.backup_hasher = RepresentativeHasher(self.backup_hasher)
        previous = None
        for repeat in range(a.probe_samples):
            result = original(self, i, *args, **kwargs)
            selected = self.backup_hasher.windows if self.backup_hasher is not None else self.hasher.windows
            for begin, end, _ in selected:
                assert np.isfinite(args[2][begin:end + 1]).all()
            if previous is not None:
                np.testing.assert_array_equal(previous[0], result[0])
                np.testing.assert_array_equal(previous[1], result[1])
            previous = tuple(x.copy() for x in result)
            logging.info('PROBE window pass=%s windows=%s rss_gib=%.3f',
                         repeat + 1, len(selected), psutil.Process().memory_info().rss / 2**30)
        raise ProbeComplete()
    if i > a.probe_samples:
        raise ProbeComplete()
    result = original(self, i, *args, **kwargs)
    logging.info('PROBE decoded sample=%s rss_gib=%.3f', i,
                 psutil.Process().memory_info().rss / 2**30)
    return result
ASMCDecoder.compute_with_hashing = measured
settings = argparse.ArgumentParser()
inference.add_default_arg_building_arguments(settings)
args = settings.parse_args([])
args.num_sequence_samples = a.samples
args.num_snp_samples = a.samples
random.seed(20260904)
np.random.seed(20260904)
try:
    inference.build_arg(args, a.root, a.map, a.mode, verbose=True)
except ProbeComplete:
    group = Path('/proc/self/cgroup').read_text().strip().split('::', 1)[1]
    cg = Path('/sys/fs/cgroup' + group)
    for name in ('memory.peak', 'memory.swap.peak', 'memory.events'):
        if (cg / name).exists():
            print(name + ': ' + (cg / name).read_text().strip(), flush=True)
    scope = (f'{a.probe_windows} representative windows x {a.probe_samples} passes; no complete sample'
             if a.probe_windows else f'{a.probe_samples} samples threaded')
    print(f'PASS: full input, scaffold={a.samples}, {scope}; '
          'probe only, no full-build checkpoint', flush=True)
else:
    raise RuntimeError('Probe did not reach its stopping point')
