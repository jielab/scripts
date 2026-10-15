#!/usr/bin/env python3
"""Query signed fixed-dimensional Atlas vectors with resumable chunk checkpoints."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
from scipy import sparse
import zarr
from tqdm import tqdm

SCORERS = ("ATAC", "DNASE", "CHIP_TF", "CHIP_HISTONE", "CAGE", "PROCAP")


def canonical_chrom(chrom):
    chrom = str(chrom).removeprefix("chr")
    return "chr" + ("M" if chrom == "MT" else chrom)


def variant_id(variant):
    if all(
        hasattr(variant, x)
        for x in ["chromosome", "position", "reference_bases", "alternate_bases"]
    ):
        return f"{canonical_chrom(variant.chromosome)}:{variant.position}:{variant.reference_bases}>{variant.alternate_bases}"
    match = re.fullmatch(r"([^:]+):(\d+):([ACGT])>([ACGT])", str(variant))
    if not match:
        raise ValueError(f"Unrecognized Atlas variant identifier: {variant}")
    chrom, pos, ref, alt = match.groups()
    return f"{canonical_chrom(chrom)}:{int(pos)}:{ref}>{alt}"


def response_matrix(answer, scorers, expected):
    """Align strictly on chromosome, position, REF, ALT; reject duplicate rows."""
    names, blocks = [], []
    for scorer in scorers:
        if scorer not in answer:
            raise ValueError(f"Missing Atlas scorer: {scorer}")
        ad = answer[scorer]
        values = ad.X.toarray() if sparse.issparse(ad.X) else np.asarray(ad.X)
        if (
            "variant" not in ad.obs
            or values.ndim != 2
            or len(ad.obs) != values.shape[0]
        ):
            raise ValueError(f"{scorer}: expected fixed-dimensional per-variant rows")
        ids = [variant_id(v) for v in ad.obs["variant"]]
        if len(set(ids)) != len(ids):
            raise ValueError(
                f"{scorer}: duplicate variant rows; gene-centric scorers are unsupported"
            )
        if set(ids) != set(expected):
            raise ValueError(
                f"{scorer}: response variants differ from the requested chromosome/position/alleles"
            )
        if not np.isfinite(values).all():
            raise ValueError(f"{scorer}: nonfinite scores")
        features = [f"{scorer}|{v}" for v in ad.var.index.astype(str)]
        if not features or len(set(features)) != len(features):
            raise ValueError(f"{scorer}: missing or duplicate track IDs")
        names.extend(features)
        mapping = {v: i for i, v in enumerate(ids)}
        blocks.append(values[[mapping[v] for v in expected]])
    return np.concatenate(blocks, axis=1).astype("f4"), names


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def query(sites, out, scorers, chunk_size, max_workers, key):
    from alphagenome.atlas import atlas
    from alphagenome.data import genome

    if sites.empty or sites.site_id.duplicated().any():
        raise ValueError("Atlas sites must be nonempty and unique")
    rows = []
    for row in sites.itertuples(index=False):
        if row.ref not in "ACGT" or len(row.ref) != 1:
            raise ValueError(f"Invalid reference base: {row.ref}")
        for alt in "ACGT":
            if alt != row.ref:
                rows.append(
                    dict(
                        site_id=row.site_id,
                        chrom=row.chrom,
                        pos=int(row.pos),
                        ref=row.ref,
                        alt=alt,
                        variant_id=f"{row.chrom}:{int(row.pos)}:{row.ref}>{alt}",
                        occupancy=int(row.occupancy),
                        atlas_group=row.atlas_group,
                    )
                )
    request = pd.DataFrame(rows)
    canonical_ids = [variant_id(v) for v in request.variant_id]
    if len(set(canonical_ids)) != len(canonical_ids):
        raise ValueError("Duplicate requested variants after chromosome normalization")
    identity = dict(
        schema=1,
        request_sha256=hashlib.sha256(
            request.to_json(orient="records").encode()
        ).hexdigest(),
        scorers=scorers,
        sdk_version=importlib.metadata.version("alphagenome"),
        total_rows=len(request),
    )
    manifest_path = out / "atlas_manifest.json"
    row_path = out / "rows.parquet"
    feature_path = out / "feature_names.tsv"
    zpath = out / "vectors.zarr"
    completed, z, features = 0, None, None
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if any(manifest.get(k) != v for k, v in identity.items()):
            raise ValueError(
                "Atlas cache belongs to different sites/scorers/SDK; use a new output directory"
            )
        if row_path.exists():
            saved = pd.read_parquet(row_path)
            completed = len(saved)
            if not saved.equals(request.iloc[:completed].reset_index(drop=True)):
                raise ValueError("Atlas cached rows do not match the request prefix")
            features = pd.read_csv(feature_path, sep="\t").feature.tolist()
            z = zarr.open(str(zpath), mode="a")
            if z.shape[0] < completed or z.shape[1] != len(features):
                raise ValueError("Atlas checkpoint is inconsistent with vectors")
            # Discard a chunk written before a crash but not committed in rows.parquet.
            z.resize((completed, len(features)))
    elif any(p.exists() for p in [row_path, feature_path, zpath]):
        raise ValueError(
            "Unversioned Atlas cache found; preserve it and use a new output directory"
        )

    manifest = {
        **identity,
        "completed_rows": completed,
        "complete": completed == len(request),
    }
    write_json(manifest_path, manifest)
    if completed == len(request):
        print(f"Atlas cache complete: {completed:,} rows; no API calls needed")
        return
    client = atlas.create(key)
    available = client.scorer_metadata()
    for scorer in scorers:
        if scorer not in available or not available[scorer].is_signed:
            raise ValueError(f"{scorer}: scorer unavailable or not signed")
    meta = {
        s: {
            "is_signed": bool(available[s].is_signed),
            "n_tracks": len(available[s].track_metadata),
        }
        for s in scorers
    }
    metadata_path = out / "atlas_metadata.json"
    if completed and json.loads(metadata_path.read_text()) != meta:
        raise ValueError("Atlas scorer metadata changed since the checkpoint")
    write_json(metadata_path, meta)

    for start in tqdm(range(completed, len(request), chunk_size), desc="Atlas chunks"):
        chunk = request.iloc[start : start + chunk_size]
        variants = [
            genome.Variant(canonical_chrom(r.chrom), int(r.pos), r.ref, r.alt)
            for r in chunk.itertuples(index=False)
        ]
        for attempt in range(4):
            try:
                answer = client.query_variants(
                    variants,
                    requested_scorers=scorers,
                    progress_bar=False,
                    max_workers=max_workers,
                )
                break
            except (PermissionError, ValueError, NotImplementedError):
                raise
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2**attempt)
        matrix, names = response_matrix(
            answer, scorers, canonical_ids[start : start + len(chunk)]
        )
        if features is None:
            features = names
            z = zarr.open(
                str(zpath),
                mode="w",
                shape=(0, len(features)),
                chunks=(min(chunk_size, 128), len(features)),
                dtype="f4",
            )
            pd.DataFrame({"feature": features}).to_csv(
                feature_path, sep="\t", index=False
            )
        elif features != names:
            raise ValueError("Atlas feature metadata changed between chunks")
        end = start + len(chunk)
        z.resize((end, len(features)))
        z[start:end, :] = matrix
        temporary = row_path.with_suffix(".parquet.tmp")
        request.iloc[:end].to_parquet(temporary, index=False)
        temporary.replace(row_path)
        manifest.update(completed_rows=end, complete=end == len(request))
        write_json(manifest_path, manifest)
    print(f"Atlas rows={len(request):,}; features={len(features):,}; zarr={zpath}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--scorers", default=",".join(SCORERS))
    parser.add_argument("--chunk", type=int, default=128)
    parser.add_argument("--max-workers", type=int, default=10)
    args = parser.parse_args()
    key = os.getenv("ALPHA_GENOME_API_KEY")
    if not key:
        raise SystemExit("ALPHA_GENOME_API_KEY is required")
    scorers = [s.strip() for s in args.scorers.split(",") if s.strip()]
    if not scorers or len(set(scorers)) != len(scorers) or set(scorers) - set(SCORERS):
        raise SystemExit(
            "Choose unique signed CenterMask scorers from: " + ",".join(SCORERS)
        )
    if args.chunk < 1 or args.max_workers < 1:
        raise SystemExit("--chunk and --max-workers must be positive")
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".query.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit("Another Atlas query is writing this output directory")
        query(
            pd.read_parquet(args.sites), out, scorers, args.chunk, args.max_workers, key
        )


if __name__ == "__main__":
    main()
