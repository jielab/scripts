#!/usr/bin/env python3
"""Compute AlphaGenome functional geometry for all 3 candidate alleles per locus."""
from __future__ import annotations
import argparse
import json
import os
import numpy as np
import pandas as pd
import zarr
from parquet_io import selected_rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--atlas-sites", required=True)
    p.add_argument("--atlas-dir", required=True)
    p.add_argument("--out-dir", required=True)
    a = p.parse_args()
    sel = pd.read_parquet(a.atlas_sites)
    cand = selected_rows(a.candidates, sel.site_id)
    rows = pd.read_parquet(os.path.join(a.atlas_dir, "rows.parquet"))
    manifest_path = os.path.join(a.atlas_dir, "atlas_manifest.json")
    if not os.path.isfile(manifest_path) or not json.load(open(manifest_path)).get(
        "complete"
    ):
        raise SystemExit(
            "Atlas cache is incomplete or unversioned; complete the atlas stage first"
        )
    z = zarr.open(os.path.join(a.atlas_dir, "vectors.zarr"), mode="r")
    if len(rows) != z.shape[0]:
        raise SystemExit("Atlas rows and vectors differ")
    if (
        rows.variant_id.duplicated().any()
        or cand.variant_id.duplicated().any()
        or sel.site_id.duplicated().any()
    ):
        raise SystemExit("Duplicate candidate, Atlas row or selected site identifiers")
    features = (
        pd.read_csv(os.path.join(a.atlas_dir, "feature_names.tsv"), sep="\t")
        .feature.astype(str)
        .tolist()
    )
    if len(features) != z.shape[1] or len(set(features)) != len(features):
        raise SystemExit("Atlas feature metadata and vectors differ")
    scorers = np.array([x.split("|", 1)[0] for x in features])
    # Feature-wise RMS standardization, then equal-weight modalities so a scorer with more tracks does not dominate.
    sums = np.zeros(z.shape[1], dtype=np.float64)
    block_rows = max(1, min(4096, 64 * 1024**2 // (8 * z.shape[1])))
    for start in range(0, z.shape[0], block_rows):
        block = np.asarray(z[start : start + block_rows], dtype=np.float64)
        if not np.isfinite(block).all():
            raise SystemExit("Atlas vectors contain nonfinite scores")
        sums += np.einsum("ij,ij->j", block, block)
    rms = np.sqrt(sums / z.shape[0])
    rms[~np.isfinite(rms) | (rms < 1e-8)] = 1.0
    weights = 1 / rms
    for s in np.unique(scorers):
        idx = np.where(scorers == s)[0]
        weights[idx] /= np.sqrt(len(idx))
    idx_map = {v: i for i, v in enumerate(rows.variant_id)}
    base = cand[cand.site_id.isin(sel.site_id)].merge(
        sel[["site_id", "atlas_group"]], on="site_id", how="left"
    )
    if set(base.variant_id) != set(rows.variant_id) or set(base.site_id) != set(
        sel.site_id
    ):
        raise SystemExit(
            "Selected sites, candidates and Atlas cache differ; re-run select/atlas"
        )
    out = []
    siteout = []
    eps = 1e-12
    for site_id, g in base.groupby("site_id", sort=False):
        g = g.copy().sort_values("alt")
        inds = [idx_map[v] for v in g.variant_id]
        M = np.asarray(z.oindex[inds, :], dtype=np.float64) * weights
        if (
            len(g) != 3
            or g.is_primary.sum() != 1
            or g.observed.sum() != g.occupancy.iloc[0]
        ):
            raise SystemExit(f"Invalid candidate/occupancy/primary layout at {site_id}")
        norms = np.linalg.norm(M, axis=1)
        obs = g.observed.to_numpy().astype(bool)
        pidx = np.flatnonzero(g.is_primary.to_numpy())[0]
        if not obs[pidx]:
            raise SystemExit(f"Primary ALT is not observed at {site_id}")
        b = M[pidx]
        nb = max(norms[pidx], eps)
        pair = []
        for i in range(3):
            for j in range(i + 1, 3):
                pair.append(float(np.linalg.norm(M[i] - M[j])))
        obsinds = np.where(obs)[0]
        opair = []
        for ii, i in enumerate(obsinds):
            for j in obsinds[ii + 1 :]:
                opair.append(float(np.linalg.norm(M[i] - M[j])))
        siteout.append(
            dict(
                site_id=site_id,
                chrom=g.chrom.iloc[0],
                pos=int(g.pos.iloc[0]),
                ref=g.ref.iloc[0],
                occupancy=int(g.occupancy.iloc[0]),
                atlas_group=g.atlas_group.iloc[0],
                potential_mean_magnitude=float(norms.mean()),
                potential_max_magnitude=float(norms.max()),
                potential_min_magnitude=float(norms.min()),
                potential_pairwise_diversity=float(np.mean(pair)),
                observed_pairwise_diversity=float(np.mean(opair)) if opair else np.nan,
                primary_alt=g.alt.iloc[pidx],
                primary_magnitude=float(norms[pidx]),
            )
        )
        for i, r in enumerate(g.itertuples(index=False)):
            dref = float(norms[i])
            dprim = float(np.linalg.norm(M[i] - b))
            cos = float(np.clip(np.dot(M[i], b) / (max(norms[i], eps) * nb), -1, 1))
            proj = (np.dot(M[i], b) / (nb * nb)) * b
            novelty = float(np.linalg.norm(M[i] - proj) / max(norms[i], eps))
            affinity = (dref - dprim) / (dref + dprim + eps)
            d_obs = []
            for j in obsinds:
                if j != i:
                    d_obs.append(float(np.linalg.norm(M[i] - M[j])))
            out.append(
                {
                    **r._asdict(),
                    "effect_magnitude": dref,
                    "distance_to_primary": dprim,
                    "cosine_to_primary": cos,
                    "novelty_perpendicular": novelty,
                    "affinity_to_primary": float(affinity),
                    "distance_to_nearest_other_observed": (
                        min(d_obs) if d_obs else np.nan
                    ),
                }
            )
    os.makedirs(a.out_dir, exist_ok=True)
    pd.DataFrame(out).to_parquet(
        os.path.join(a.out_dir, "geometry_candidates.parquet"), index=False
    )
    pd.DataFrame(siteout).to_parquet(
        os.path.join(a.out_dir, "geometry_sites.parquet"), index=False
    )
    pd.DataFrame({"feature": features, "rms": rms}).to_csv(
        os.path.join(a.out_dir, "geometry_feature_scaling.tsv"), sep="\t", index=False
    )
    print(f"geometry loci={len(siteout):,} candidate alleles={len(out):,}")


if __name__ == "__main__":
    main()
