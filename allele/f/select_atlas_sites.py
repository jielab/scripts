#!/usr/bin/env python3
"""Keep all multi-allelic sites and sample controls without loading all candidates."""
import argparse
from collections import Counter
import os
import numpy as np
import pandas as pd
from parquet_io import batches

COLUMNS = ["site_id", "chrom", "pos", "ref", "context", "occupancy", "primary_alt"]


def site_batches(path):
    previous = None
    for part in batches(path, COLUMNS):
        part = part.drop_duplicates("site_id")
        part = part[part.site_id != previous]
        if len(part):
            previous = part.site_id.iloc[-1]
            yield part


def keep_smallest(pools, key, group, limit):
    if limit < 1:
        return
    old = pools.get(key)
    pool = pd.concat([old, group], ignore_index=True) if old is not None else group
    pools[key] = pool.nsmallest(limit, "_draw")


def main():
    p = argparse.ArgumentParser()
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--sites")
    source.add_argument("--candidates")
    p.add_argument("--ratio", type=float, default=2)
    p.add_argument("--seed", type=int, default=20261015)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    if not np.isfinite(a.ratio) or a.ratio < 0:
        raise SystemExit("--ratio must be nonnegative and finite")
    path = a.sites or a.candidates
    multi_parts = []
    counts = Counter()
    for s in site_batches(path):
        multi = s[s.occupancy >= 2].copy()
        if len(multi):
            multi_parts.append(multi)
            counts.update(multi.groupby(["ref", "context"]).size().to_dict())
    if not multi_parts:
        raise SystemExit("No multiallelic sites found")
    multi = pd.concat(multi_parts, ignore_index=True)
    del multi_parts
    wanted = {key: int(np.ceil(n * a.ratio)) for key, n in counts.items()}
    per_ref = {
        ref: int(np.ceil(n * a.ratio)) for ref, n in multi.groupby("ref").size().items()
    }
    rng = np.random.default_rng(a.seed)
    exact, relaxed = {}, {}
    # Random priorities give a uniform sample within each stratum, independent of batching.
    for s in site_batches(path):
        bi = s[s.occupancy == 1].copy()
        bi["_draw"] = rng.random(len(bi))
        for key, group in bi.groupby(["ref", "context"]):
            keep_smallest(exact, key, group, wanted.get(key, 0))
        for ref, group in bi.groupby("ref"):
            keep_smallest(relaxed, ref, group, per_ref.get(ref, 0))
    ctrl = (
        pd.concat(exact.values(), ignore_index=True) if exact else multi.iloc[:0].copy()
    )
    parts = [
        multi.assign(atlas_group="multiallelic"),
        ctrl.assign(atlas_group="biallelic_control"),
    ]
    chosen = set(ctrl.site_id)
    for ref, limit in per_ref.items():
        need = max(0, limit - int((ctrl.ref == ref).sum()))
        if need and ref in relaxed:
            pool = relaxed[ref]
            fill = pool[~pool.site_id.isin(chosen)].nsmallest(need, "_draw")
            parts.append(fill.assign(atlas_group="biallelic_control_relaxed"))
    target = pd.concat(parts, ignore_index=True).drop(columns="_draw", errors="ignore")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    target.to_parquet(a.out, index=False)
    print(target.atlas_group.value_counts().to_string())
    print(f"selected loci={len(target):,}; Atlas candidate alleles={3*len(target):,}")


if __name__ == "__main__":
    main()
