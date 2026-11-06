"""Fixed, trait-weighted LD-block geometry with no PCA or clinical coordinates.

Distances average standardized squared differences within each LD map, then
average across maps.  Means, variance floors and retained blocks use build only.
This is a declared genetic-profile distance, not IBD or genome identity.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import hashlib
import json
from pathlib import Path

import numpy as np


@dataclass
class BlockGeometry:
    input_width: int
    keep: list
    means: list
    scales: list
    divisors: list
    retained_maps: list
    floor_quantile: float
    representation_id: str = "csx_ldblocks_fixed"
    sd_floors: dict | None = None
    input_maps: list | None = None

    @classmethod
    def fit(cls, build_blocks, block_maps, *, floor_quantile=0.1,
            row_batch=1024, variance_epsilon=1e-12):
        x = build_blocks if hasattr(build_blocks, "shape") else np.asarray(build_blocks)
        maps = np.asarray(block_maps, dtype=str)
        if len(x.shape) != 2 or len(x) < 2 or maps.shape != (x.shape[1],):
            raise ValueError("build_blocks must be N x blocks with one map label per block")
        if not 0 <= floor_quantile <= 1 or row_batch < 1:
            raise ValueError("Invalid variance-floor quantile or row batch")
        # Merge batch moments: stable and bounded memory for a memmapped matrix.
        n = 0
        mean = np.zeros(x.shape[1], dtype=np.float64)
        m2 = np.zeros(x.shape[1], dtype=np.float64)
        for start in range(0, len(x), row_batch):
            b = np.asarray(x[start:start + row_batch], dtype=np.float64)
            if not np.isfinite(b).all():
                raise ValueError("Block cache contains missing/nonfinite values; repair upstream QC")
            bn = len(b)
            bm = b.mean(axis=0)
            delta = bm - mean
            m2 += ((b - bm) ** 2).sum(axis=0) + delta ** 2 * n * bn / (n + bn)
            mean += delta * bn / (n + bn)
            n += bn
        sd = np.sqrt(m2 / n)
        keep = np.flatnonzero(sd > variance_epsilon)
        if not len(keep):
            raise ValueError("No variable blocks remain in build")
        remaining_maps = np.unique(maps[keep])
        scale = sd[keep].copy()
        divisors = np.empty(len(keep), dtype=float)
        floors = {}
        for label in remaining_maps:
            where = maps[keep] == label
            floor = np.quantile(scale[where], floor_quantile)
            floors[str(label)] = float(max(floor, variance_epsilon))
            scale[where] = np.maximum(scale[where], max(floor, variance_epsilon))
            divisors[where] = np.sqrt(where.sum() * len(remaining_maps))
        return cls(x.shape[1], keep.tolist(), mean[keep].tolist(),
                   scale.tolist(), divisors.tolist(), maps[keep].tolist(),
                   float(floor_quantile), sd_floors=floors, input_maps=maps.tolist())

    def transform(self, blocks):
        x = np.asarray(blocks)
        if x.ndim != 2 or x.shape[1] != self.input_width:
            raise ValueError("Block schema differs from the frozen geometry")
        a = np.asarray(x[:, self.keep], dtype=np.float64)
        if not np.isfinite(a).all():
            raise ValueError("Invalid genotype block values cannot be assigned a support score")
        return (a - np.asarray(self.means)) / (np.asarray(self.scales) * self.divisors)

    def transform_to_npy(self, blocks, output, *, row_batch=512, dtype="float32"):
        # Transform/query preprocessing is label-free and uses the frozen scale.
        if dtype not in {"float32", "float64"} or row_batch < 1:
            raise ValueError("Invalid storage dtype or batch")
        x = blocks if hasattr(blocks, "shape") else np.asarray(blocks)
        if len(x.shape) != 2 or x.shape[1] != self.input_width:
            raise ValueError("Block schema differs from the frozen geometry")
        out = np.lib.format.open_memmap(output, mode="w+", dtype=dtype,
                                        shape=(len(x), len(self.keep)))
        try:
            for start in range(0, len(x), row_batch):
                out[start:start + row_batch] = self.transform(x[start:start + row_batch])
            out.flush()
        finally:
            out._mmap.close()
        return Path(output)

    def fingerprint(self):
        data = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(data.encode()).hexdigest()
