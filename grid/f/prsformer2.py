#!/usr/bin/env python3
"""Opt-in official PRSformer adapter for the immutable GRID four-role cohort.

This module prepares ALT-dosage inputs and delegates neural training/inference
to ``prsformer.py``, which imports the actual 23andMe implementation:
https://github.com/23andMe/PRSformer (src/model.py, modules.py, utils.py).
No fitted weights are distributed by that repository. This is a new training
experiment using its architecture, not a replication of its very large 23andMe
training cohort or a replacement algorithm with the same name.

TRAIN fits genotype QC, covariate transforms/baselines and the neural model;
TUNE chooses the checkpoint. CALIBRATION and TEST outcomes never reach the
worker. Both roles receive inference, and their evaluation stays with GRID.
The original participants and family groups are preserved exactly.

Genotypes use a reference/HM3-to-GRCh38 target map built independently from PVARs. This
universe is shared across multitask traits; each GWAS's available posterior
SNP set can differ. PGEN is decoded twice in bounded variant chunks: TRAIN-only
QC, then direct writes to a sample-major float16 disk memmap. There is no dense
N-by-M RAM array or variant-major scratch copy. Disk needs still scale as 2NM
bytes and can approach a terabyte; resource estimates and configurable limits
are checked before allocation. Direct strided writes trade disk traffic for
avoiding a second huge temporary file. Decode chunks are coalesced in a bounded
RAM tile (write_buffer_variants) before flushing wide row slices. Increase buffers only within the
reported RAM budget. The memory estimate includes pages touched in the mmap;
OS file cache and GPU activations are not a guaranteed upper bound.

API: run(cfg, root, traits, stage='all'), stages plan/prepare/train/predict/all.
Dependencies at preparation: numpy, pandas, pgenlib (zstandard for .pvar.zst).
Training additionally needs the dependencies of the official model, CUDA and
NATTEN for genome-scale neighborhood attention. There is no CPU/global-attention
fallback at genome scale. A tiny fake-reader test verifies orchestration only.
"""
from __future__ import annotations
from output_layout import project_path, trait_path, module_root
from prsformer_inputs import cohort_dir, target_manifest, ensure_inputs

import argparse
from contextlib import contextmanager
from copy import deepcopy
import csv
import hashlib
import json
import mmap
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

from core2 import aligned, atomic_json, atomic_table, check_roles, identified, read_table, stamp, text_open

TRAITS = {"height": "continuous", "ldl": "continuous", "t2dm": "binary"}
GIB = 1024 ** 3
ROLE_SPLIT = {"train": "train", "tune": "validation", "calibration": "test", "test": "test"}
MAP_COLUMNS = ["SNP", "CHR", "BP38", "TARGET_ID", "REF38", "ALT38"]
PREPARED_FORMAT = "grid-prsformer-20261127"


WORKER_VALUE_OPTIONS = {
    "device": "cuda", "attention": "neighborhood", "embed_dim": 64,
    "heads": 4, "layers": 2, "ff_dim": 128, "kernel_size": 385,
    "dilation": "1", "batch_size": 1, "accumulation_steps": 64,
    "epochs": 30, "patience": 3, "min_delta": 1e-4, "lr": 5e-4,
    "weight_decay": .05, "clip_grad": 1., "amp": "fp16", "workers": 0,
    "threads": 4,
}


def _section(options, name):
    value = options.get(name, {})
    if isinstance(value, bool):
        return {"enabled": value}
    if not isinstance(value, dict):
        raise ValueError(f"prsformer.{name} must be a mapping or boolean")
    return value


def _options(cfg, name):
    base = cfg.get("prsformer", {})
    return {**{k: v for k, v in base.items() if k not in {"prepare", "train", "predict"}},
            **_section(base, name)}


def _traits(traits):
    traits = list(traits)
    if not traits or len(set(traits)) != len(traits) or not set(traits) <= set(TRAITS):
        raise ValueError("PRSformer traits must be a unique nonempty subset of height, ldl, t2dm")
    return traits


def worker_cohort(cfg, root, traits):
    """Read development outcomes only, retaining labels solely for TRAIN/TUNE.

    Calibration values are masked before numeric conversion and validation.
    The private test-outcome path is never opened by this module.
    """
    root = Path(root)
    traits = _traits(traits)
    path = cohort_dir(root, traits)/"participants.tsv.gz"
    roster = identified(read_table(path), "PRSformer participants")
    check_roles(roster)
    if roster.ancestry.isna().any() or roster.group_id.astype(str).str.strip().eq("").any():
        raise ValueError("PRSformer participants need complete ancestry and family groups")
    covariates = list(cfg.get("benchmark", {}).get("covariates", []))
    if len(set(covariates)) != len(covariates) or not set(covariates) <= set(roster):
        raise ValueError("Missing or repeated common benchmark covariates")
    data = roster[["id", "ancestry", "group_id", "role", *covariates]].copy()
    data = data.rename(columns={"id": "eid", "ancestry": "target", "group_id": "group", "role": "grid_role"})
    data["split"] = data.grid_role.map(ROLE_SPLIT)
    allowed = roster.role.isin(["train", "tune"]).to_numpy()
    records = [stamp(path)]
    for trait in traits:
        source = cohort_dir(root, traits)/f"{trait}.development.tsv.gz"
        outcomes = identified(read_table(source), f"{trait} development outcomes")
        if "y" not in outcomes:
            raise ValueError(f"{source} must contain y")
        if set(outcomes.id) - set(roster.id):
            raise ValueError(f"{trait} development outcomes contain IDs outside the fixed roster")
        test_ids = set(roster.loc[roster.role.eq("test"), "id"])
        if outcomes.loc[outcomes.id.isin(test_ids), "y"].notna().any():
            raise ValueError(f"{trait} development outcome file contains forbidden TEST labels")
        # Only TRAIN/TUNE labels are handed to the learner; CALIBRATION is held out.
        selected = outcomes.loc[outcomes.id.isin(roster.loc[allowed, "id"]), ["id", "y"]]
        y = aligned(selected, roster.id, trait, allow_missing=True).y
        values = pd.to_numeric(y, errors="raise").to_numpy(dtype=float, copy=True)
        if np.isinf(values).any():
            raise ValueError(f"{trait}: infinite development outcomes")
        values[~allowed] = np.nan
        for role in ("train", "tune"):
            observed = values[roster.role.eq(role).to_numpy() & np.isfinite(values)]
            if len(observed) < 2:
                raise ValueError(f"{trait} needs at least two observed labels in {role}")
            if TRAITS[trait] == "binary" and (not np.isin(observed, [0, 1]).all() or len(np.unique(observed)) != 2):
                raise ValueError(f"{trait} needs both 0/1 classes in {role}")
        data[trait] = values
        records.append(stamp(source))
    for column in covariates:
        data[column] = pd.to_numeric(data[column], errors="raise")
        if np.isinf(data[column].to_numpy(float)).any():
            raise ValueError(f"Infinite covariate: {column}")
        if not np.isfinite(data.loc[data.split.eq("train"), column].to_numpy(float)).any():
            raise ValueError(f"No observed TRAIN value for covariate {column}")
    return roster, data, covariates, records


def _pvar_indices(path, wanted, chromosome):
    """Stream the full PVAR; only selected reference variants enter memory."""
    found = {}
    with text_open(path) as stream:
        reader = csv.reader((line for line in stream if not line.startswith("##") and line.strip()), delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"Empty PVAR: {path}") from exc
        header = [x.lstrip("#") for x in header]
        if not {"CHROM", "POS", "ID", "REF", "ALT"} <= set(header):
            raise ValueError(f"Invalid PVAR header: {path}")
        index = {x: header.index(x) for x in ("CHROM", "POS", "ID", "REF", "ALT")}
        for row_number, row in enumerate(reader):
            target_id = row[index["ID"]]
            if target_id not in wanted:
                continue
            if target_id in found:
                raise ValueError(f"Duplicate selected TARGET_ID in PVAR: {target_id}")
            current = (int(row[index["CHROM"]].lower().removeprefix("chr")),
                       int(row[index["POS"]]), row[index["REF"]].upper(), row[index["ALT"]].upper())
            expected = wanted[target_id]
            if current != expected or current[0] != chromosome:
                raise ValueError(f"Target mapping does not match actual GRCh38 PVAR: {target_id}")
            found[target_id] = row_number
    if len(found) != len(wanted):
        raise ValueError(f"{len(wanted) - len(found)} selected TARGET_IDs missing from {path}")
    return found


def variant_sources(cfg, root, traits, roster):
    root = Path(root)
    options = _options(cfg, "prepare")
    manifest_path = target_manifest(cfg, root, traits)
    manifest = read_table(manifest_path)
    if not {"chromosome", "path", "target_build", "n_variants"} <= set(manifest):
        raise ValueError("PRSformer needs its prepared target/manifest.tsv")
    if not manifest.target_build.isin(["GRCh38", "38", "hg38"]).all():
        raise ValueError("PRSformer target map must describe GRCh38 coordinates")
    manifest["chromosome"] = pd.to_numeric(manifest.chromosome, errors="raise").astype(int)
    chromosomes = [int(c) for c in cfg.get("chromosomes", range(1, 23))]
    if len(set(chromosomes)) != len(chromosomes) or not set(chromosomes) <= set(range(1, 23)):
        raise ValueError("Use unique autosomes 1..22")
    if manifest.chromosome.duplicated().any() or set(manifest.chromosome) != set(chromosomes):
        raise ValueError("Target manifest chromosomes must equal the configured chromosome set")
    directory = Path(cfg["inputs"]["genotype_dir"])
    tables, sources, records = [], {}, [stamp(manifest_path)]
    for row in manifest.sort_values("chromosome").itertuples(index=False):
        chromosome = int(row.chromosome)
        mapping_path = Path(row.path)
        if not mapping_path.is_absolute():
            mapping_path = manifest_path.parent / mapping_path
        mapping = read_table(mapping_path)
        if not set(MAP_COLUMNS) <= set(mapping) or mapping[MAP_COLUMNS].isna().any().any():
            raise ValueError(f"Incomplete target map: {mapping_path}")
        if len(mapping) != int(row.n_variants) or mapping.TARGET_ID.duplicated().any() or mapping.SNP.duplicated().any():
            raise ValueError(f"Target map counts or identifiers are inconsistent: {mapping_path}")
        for column in ("CHR", "BP38"):
            numeric = pd.to_numeric(mapping[column], errors="raise")
            if not np.equal(numeric, np.floor(numeric)).all():
                raise ValueError(f"Noninteger target coordinate: {column}")
            mapping[column] = numeric.astype(int)
        if not mapping.CHR.eq(chromosome).all() or not mapping.BP38.gt(0).all():
            raise ValueError("Invalid target coordinates")
        for column in ("REF38", "ALT38"):
            mapping[column] = mapping[column].str.upper()
            if not mapping[column].isin(list("ACGT")).all():
                raise ValueError("PRSformer target universe must contain only biallelic SNPs")
        if mapping.REF38.eq(mapping.ALT38).any():
            raise ValueError("REF38 and ALT38 cannot match")
        prefix = directory / f"chr{chromosome}"
        pgen, psam = Path(str(prefix) + ".pgen"), Path(str(prefix) + ".psam")
        pvar = Path(str(prefix) + ".pvar")
        if not pvar.is_file():
            pvar = Path(str(prefix) + ".pvar.zst")
        records.extend(stamp(path) for path in (mapping_path, pgen, pvar, psam))
        samples = identified(read_table(psam), f"chr{chromosome} PSAM")
        lookup = {eid: i for i, eid in enumerate(samples.id)}
        if not set(roster.id) <= set(lookup):
            raise ValueError(f"chr{chromosome}: fixed roster contains absent genotype samples")
        original = np.asarray([lookup[eid] for eid in roster.id], dtype=np.uint32)
        subset = np.sort(original)
        reorder = np.searchsorted(subset, original)
        wanted = {r.TARGET_ID: (int(r.CHR), int(r.BP38), r.REF38, r.ALT38) for r in mapping.itertuples(index=False)}
        indices = _pvar_indices(pvar, wanted, chromosome)
        mapping["source_index"] = mapping.TARGET_ID.map(indices).astype(np.uint32)
        mapping = mapping.rename(columns={"BP38": "BP", "REF38": "REF", "ALT38": "ALT"})
        tables.append(mapping[["CHR", "BP", "SNP", "REF", "ALT", "TARGET_ID", "source_index"]])
        sources[chromosome] = {"chromosome": chromosome, "pgen": pgen, "pvar": pvar, "psam": psam,
                               "raw_sample_ct": len(samples), "subset": subset, "reorder": reorder}
    variants = pd.concat(tables, ignore_index=True).sort_values(["CHR", "BP", "SNP"], kind="stable").reset_index(drop=True)
    if variants.SNP.duplicated().any() or not len(variants):
        raise ValueError("The shared PRSformer universe must have unique, nonempty rsIDs")
    limit = int(options.get("max_variants", 0))
    if limit < 0:
        raise ValueError("max_variants must be nonnegative")
    if limit:
        variants = variants.iloc[:limit].copy()
        print("PRSformer: explicit max_variants is a small-data experiment, not full-genome benchmarking", flush=True)
    return variants, sources, records


@contextmanager
def pgen_reader(source, subset):
    try:
        import pgenlib
    except ImportError as exc:
        raise RuntimeError("PRSformer preparation needs pgenlib; no BED/hard-call fallback is used") from exc
    with pgenlib.PvarReader(os.fsencode(source["pvar"])) as pvar:
        with pgenlib.PgenReader(os.fsencode(source["pgen"]), raw_sample_ct=source["raw_sample_ct"],
                               sample_subset=subset, pvar=pvar) as reader:
            yield reader


def dosage_blocks(variants, sources, chunk_variants, reader_factory=pgen_reader):
    """Yield variant-by-sample float32 blocks in the original frozen row order."""
    for chromosome, block_variants in variants.groupby("CHR", sort=True):
        source = sources[int(chromosome)]
        with reader_factory(source, source["subset"]) as reader:
            for start in range(0, len(block_variants), chunk_variants):
                block = block_variants.iloc[start:start + chunk_variants]
                values = np.empty((len(block), len(source["subset"])), dtype=np.float32, order="C")
                indices = block.source_index.to_numpy(dtype=np.uint32)
                order = np.argsort(indices)
                reader.read_dosages_list(np.ascontiguousarray(indices[order]), values, allele_idx=1)
                # PVAR files need not be ordered by BP. Decode in native order,
                # then restore both the model's positional order and frozen IDs.
                values = np.ascontiguousarray(values[np.argsort(order)][:, source["reorder"]])
                missing = (values == -9) | ~np.isfinite(values)
                if np.any((~missing) & ((values < 0) | (values > 2))):
                    raise ValueError(f"chr{chromosome}: unrecognized dosage outside [0, 2]")
                values[missing] = -1
                yield block, values


def resource_plan(n, m, options, n_traits=3, metadata_bytes=0):
    """Conservative preparation working-set estimate plus model lower bounds.

    GPU states and activation estimates are *lower bounds*, never a promise
    that a particular device can fit the official network.
    """
    n, m = int(n), int(m)
    chunk = int(options.get("chunk_variants", 128))
    if min(n, m, chunk) < 1:
        raise ValueError("Nonempty cohort/variants and positive chunk_variants required")
    block = min(chunk, m)
    page = int(os.sysconf("SC_PAGE_SIZE"))
    # Decoding/reorder/QC and output casting may coexist. mmap spans are bounded
    # by the number of per-row pages touched, then explicitly released.
    decode = 20 * n * block
    commit = min(m, int(options.get("write_buffer_variants", chunk)))
    if commit < 1:
        raise ValueError("write_buffer_variants must be positive")
    write_buffer = 2 * n * commit
    dirty_pages = min(2 * n * m, n * (2 * commit + 2 * page))
    reserve = float(options.get("disk_reserve_gib", 2.0)) * GIB
    prep_ram = decode + write_buffer + dirty_pages + int(metadata_bytes) + 256 * 1024 ** 2
    dimension, layers = int(options.get("embed_dim", 64)), int(options.get("layers", 2))
    feedforward = int(options.get("ff_dim", 128))
    params = int(m * dimension * (2 + n_traits) + layers *
                 (4 * dimension ** 2 + 2 * dimension * feedforward + 9 * dimension + feedforward)
                 + 2 * dimension + n_traits)
    batch, workers = int(options.get("batch_size", 1)), int(options.get("workers", 0))
    return {"n_samples": n, "n_variants": m, "dtype": "float16", "dosage": "ALT", "missing": -1,
            "chunk_variants": chunk, "write_buffer_variants": commit, "estimated_write_buffer_bytes": write_buffer, "genotype_bytes": 2 * n * m,
            "required_disk_bytes": int(2 * n * m + reserve + 65536),
            "estimated_preparation_ram_bytes": int(prep_ram), "estimated_decode_buffer_bytes": int(decode),
            "estimated_mmap_touched_bytes": int(dirty_pages),
            "estimated_training_cpu_batch_bytes": int(4 * m * batch * (2 * max(workers, 1) + 2)),
            "estimated_parameters": params, "gpu_weights_gradients_adam_lower_bound_bytes": 16 * params,
            "gpu_activations_lower_bound_bytes": int(2 * batch * m * dimension),
            "memory_note": "Preparation estimate excludes OS cache/allocator overhead; GPU values exclude most activations and kernels"}


def check_resources(plan, cache, options, *, check_disk=True, allocated=0):
    ram_limit = float(options.get("max_ram_gib", 8.0))
    disk_limit = options.get("max_disk_gib")
    if not np.isfinite(ram_limit) or ram_limit <= 0:
        raise ValueError("max_ram_gib must be positive and finite")
    if plan["estimated_preparation_ram_bytes"] > ram_limit * GIB:
        raise ValueError(f"PRSformer preparation estimate {plan['estimated_preparation_ram_bytes']/GIB:.2f} GiB exceeds max_ram_gib={ram_limit}; reduce chunk_variants or request more RAM")
    if check_disk:
        if disk_limit is not None and (not np.isfinite(float(disk_limit)) or float(disk_limit) <= 0):
            raise ValueError("max_disk_gib must be positive and finite")
        if disk_limit is not None and plan["required_disk_bytes"] > float(disk_limit) * GIB:
            raise ValueError("PRSformer disk estimate exceeds configured max_disk_gib")
        free = shutil.disk_usage(cache).free
        if plan["required_disk_bytes"] - allocated > free:
            raise ValueError(f"PRSformer needs {(plan['required_disk_bytes']-allocated)/GIB:.2f} GiB additional disk, only {free/GIB:.2f} GiB available; choose a larger cache filesystem")


def _print_plan(plan, label):
    print(f"PRSformer {label}: {plan['n_samples']:,} participants x {plan['n_variants']:,} SNPs; "
          f"float16 matrix {plan['genotype_bytes']/GIB:.2f} GiB, disk including reserve "
          f"{plan['required_disk_bytes']/GIB:.2f} GiB, preparation RAM estimate "
          f"{plan['estimated_preparation_ram_bytes']/GIB:.2f} GiB; "
          f"model state lower bound {plan['gpu_weights_gradients_adam_lower_bound_bytes']/GIB:.2f} GiB GPU (activations additional)", flush=True)


def _signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _stable(records):
    if any(stamp(record["path"]) != record for record in records):
        raise RuntimeError("A PRSformer input changed during preparation; no cache was marked complete")


def _release_pages(matrix):
    matrix.flush()
    # Linux supports MADV_DONTNEED; it avoids retaining every written row page in
    # this process across successive strided writes. The file remains intact.
    if hasattr(matrix._mmap, "madvise") and hasattr(mmap, "MADV_DONTNEED"):
        matrix._mmap.madvise(mmap.MADV_DONTNEED)


def _prepared_paths(root, options):
    cache = Path(options.get("cache_dir", project_path(root, "tmp/prsformer/cache"))).expanduser().resolve()
    return cache, {"data": cache / "data.tsv.gz", "variants": cache / "variants.tsv.gz",
                   "genotypes": cache / "genotypes.npy", "receipt": cache / "preparation.json"}


def prepare(cfg, root, traits, *, plan_only=False, reader_factory=pgen_reader):
    """Prepare bounded disk-backed dosages; reader_factory is a test seam only."""
    traits = _traits(traits)
    options = _options(cfg, "prepare")
    ensure_inputs(cfg, root, traits)
    roster, data, covariates, records = worker_cohort(cfg, root, traits)
    variants, sources, source_records = variant_sources(cfg, root, traits, roster)
    records += source_records
    cache, paths = _prepared_paths(root, options)
    cache.mkdir(parents=True, exist_ok=True)
    train = data.split.eq("train").to_numpy()
    maf_cut, call_cut = float(options.get("maf", .01)), float(options.get("call_rate", .98))
    if not 0 <= maf_cut <= .5 or not 0 < call_cut <= 1:
        raise ValueError("Use 0 <= maf <= .5 and 0 < call_rate <= 1")
    metadata_bytes = sum(int(d.memory_usage(index=True, deep=True).sum()) for d in (roster, data, variants))
    combined = {**_options(cfg, "train"), **options}
    upper_plan = resource_plan(len(data), len(variants), combined, len(traits), metadata_bytes)
    _print_plan(upper_plan, "candidate plan")
    check_resources(upper_plan, cache, options, check_disk=plan_only)
    if plan_only:
        return {"status": "planned", "resources": upper_plan, "qc_not_yet_applied": True}
    contract = {"format": PREPARED_FORMAT, "records": records, "traits": traits,
                "covariates": covariates, "maf": maf_cut, "call_rate": call_cut,
                "max_variants": int(options.get("max_variants", 0)), "schema": 1,
                "genotype_dir": str(Path(cfg["inputs"]["genotype_dir"]).expanduser().resolve()),
                "target_manifest": str(target_manifest(cfg, root, traits).expanduser().resolve()),
                "chromosomes": sorted(int(c) for c in cfg.get("chromosomes", range(1, 23)))}
    signature = _signature(contract)
    import fcntl
    with (cache / ".prepare.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if paths["receipt"].is_file() and not options.get("replace", False):
            receipt = json.loads(paths["receipt"].read_text())
            if receipt.get("signature") != signature:
                raise ValueError("PRSformer prepared inputs/configuration changed; use a new cache or explicit prepare.replace")
            if any(not paths[key].is_file() or stamp(paths[key]) != receipt["outputs"].get(key)
                   for key in ("data", "variants", "genotypes")):
                raise ValueError("Prepared PRSformer cache is incomplete or changed")
            print(f"PRSformer: reuse verified prepared cache {cache}", flush=True)
            return receipt
        if not options.get("replace", False) and not (cache / "dosage-progress.json").is_file() and any(paths[key].exists() for key in ("data", "variants", "genotypes")):
            raise ValueError("Unrecorded PRSformer cache exists; use a new cache or explicit prepare.replace")
        qc_table, qc_state = cache / "qc_variants.tsv.gz", cache / "qc-progress.json"
        qc = json.loads(qc_state.read_text()) if qc_state.is_file() else None
        if qc and not options.get("replace", False):
            if qc.get("signature") != signature or not qc_table.is_file() or stamp(qc_table) != qc.get("output"):
                raise ValueError("PRSformer QC progress changed; use a new cache or explicit prepare.replace")
            final = pd.read_csv(qc_table, sep="\t", dtype={"CHR": int, "BP": int, "SNP": str,
                "REF": str, "ALT": str, "TARGET_ID": str, "source_index": np.uint32}, float_precision="round_trip")
            print(f"PRSformer: reuse verified TRAIN QC for {len(final):,} SNPs", flush=True)
        else:
            kept = []
            last = time.monotonic()
            for block, values in dosage_blocks(variants, sources, upper_plan["chunk_variants"], reader_factory):
                observed = values[:, train]
                valid = observed >= 0
                calls = valid.sum(axis=1)
                freq = np.divide(np.where(valid, observed, 0).sum(axis=1, dtype=np.float64), 2 * calls,
                                 out=np.zeros(len(block)), where=calls > 0)
                maf = np.minimum(freq, 1 - freq)
                call_rate = calls / int(train.sum())
                use = (calls > 0) & (maf > 0) & (maf >= maf_cut) & (call_rate >= call_cut)
                if use.any():
                    chosen = block.loc[use].copy()
                    chosen["MAF_train"], chosen["call_rate_train"] = maf[use], call_rate[use]
                    kept.append(chosen)
                if time.monotonic() - last > 30:
                    print(f"PRSformer TRAIN genotype QC: chr{int(block.CHR.iloc[0])}, {sum(len(x) for x in kept):,} SNPs retained so far", flush=True)
                    last = time.monotonic()
            if not kept:
                raise ValueError("No SNP passed TRAIN-only MAF/call-rate QC")
            final = pd.concat(kept, ignore_index=True)
            _stable(records)
            atomic_table(final, qc_table)
            atomic_json({"signature": signature, "output": stamp(qc_table)}, qc_state)
        actual_plan = resource_plan(len(data), len(final), combined, len(traits), metadata_bytes)
        _print_plan(actual_plan, "QC-passed plan")
        # Credit disk blocks already allocated to this exact resumable matrix.
        progress=cache/"dosage-progress.json"
        allocated=0
        if progress.is_file() and not options.get("replace",False):
            saved=json.loads(progress.read_text())
            partial=cache/"genotypes.inprogress.npy"
            matrix=partial if partial.is_file() else paths["genotypes"]
            if matrix.is_file() and saved.get("signature")==signature:
                stat=matrix.stat()
                if saved.get("file_identity")==[stat.st_dev,stat.st_ino]:
                    allocated=min(stat.st_blocks*512, actual_plan["genotype_bytes"])
        check_resources(actual_plan, cache, options, allocated=allocated)
        _stable(records)
        from prepare_progress import write_matrix
        # prepare.replace is explicit permission to rebuild this cache.
        if options.get("replace", False):
            for name in ("genotypes.inprogress.npy", "dosage-progress.json", "genotypes.npy", "preparation.json"):
                (cache / name).unlink(missing_ok=True)
        def remaining_blocks(start):
            for _, values in dosage_blocks(final.iloc[start:], sources, actual_plan["chunk_variants"], reader_factory):
                yield values
        write_matrix(cache, signature, (len(data), len(final)), remaining_blocks, _release_pages,
                     commit_columns=actual_plan["write_buffer_variants"])
        _stable(records)
        atomic_table(data, paths["data"])
        atomic_table(final.drop(columns="source_index"), paths["variants"])
        receipt = {**contract, "signature": signature, "status": "prepared", "resources": actual_plan,
                   "shape": [len(data), len(final)], "genome_build": "GRCh38", "dosage": "ALT",
                   "missing": -1, "role_counts": roster.role.value_counts().to_dict(),
                   "phenotype_use": {"train": "fit", "tune": "checkpoint", "calibration": "masked", "test": "masked"},
                   "variant_universe": "common LD/HM3 reference intersect target GRCh38 PVAR; TRAIN-only QC; trait GWAS availability can differ",
                   "small_data_experiment": bool(options.get("max_variants", 0)),
                   "outputs": {key: stamp(paths[key]) for key in ("data", "variants", "genotypes")}}
        atomic_json(receipt, paths["receipt"])
        (cache / "dosage-progress.json").unlink(missing_ok=True)
        return receipt


def validate_prepared(cfg, root, traits):
    options = _options(cfg, "prepare")
    cache, paths = _prepared_paths(root, options)
    if not paths["receipt"].is_file():
        raise ValueError("Run PRSformer prepare before train/predict")
    receipt = json.loads(paths["receipt"].read_text())
    if receipt.get("format") != PREPARED_FORMAT:
        raise ValueError("This training entry requires dosages prepared by this release")
    if receipt.get("traits") != list(traits):
        raise ValueError("Prepared PRSformer task order differs from requested traits")
    _stable(receipt["records"])
    for key in ("data", "variants", "genotypes"):
        if stamp(paths[key]) != receipt["outputs"][key]:
            raise ValueError(f"Prepared PRSformer {key} changed")
    if receipt.get("covariates") != list(cfg.get("benchmark", {}).get("covariates", [])):
        raise ValueError("PRSformer covariates differ from the prepared common benchmark covariates")
    if float(options.get("maf", .01)) != receipt["maf"] or float(options.get("call_rate", .98)) != receipt["call_rate"]:
        raise ValueError("PRSformer QC thresholds differ from the prepared cache")
    expected_source = {
        "genotype_dir": str(Path(cfg["inputs"]["genotype_dir"]).expanduser().resolve()),
        "target_manifest": str(target_manifest(cfg, root, traits).expanduser().resolve()),
        "chromosomes": sorted(int(c) for c in cfg.get("chromosomes", range(1, 23))),
        "max_variants": int(options.get("max_variants", 0)),
    }
    if any(receipt.get(key) != value for key, value in expected_source.items()):
        raise ValueError("PRSformer genotype source, target universe or SNP limit differs from the prepared cache")
    return paths, receipt


def worker_command(cfg, root, traits, stage):
    if stage not in {"train", "predict"}:
        raise ValueError("Worker stage must be train or predict")
    root = Path(root).resolve()
    paths, receipt = validate_prepared(cfg, root, traits)
    options = _options(cfg, stage)
    upstream = Path(options.get("upstream_dir", "")).expanduser().resolve()
    for filename in ("model.py", "modules.py", "utils.py"):
        if not (upstream / "src" / filename).is_file():
            raise FileNotFoundError(f"Set prsformer.upstream_dir to the official 23andMe/PRSformer checkout: missing src/{filename}")
    run_dir = Path(options.get("run_dir", project_path(root, "tmp/prsformer/run"))).expanduser().resolve()
    worker = Path(__file__).with_name("prsformer.py")
    if not worker.is_file():
        raise FileNotFoundError(worker)
    command = [str(options.get("python", sys.executable)), str(worker), stage,
               "--analysis-root", str(root), "--genotypes", str(paths["genotypes"]),
               "--data", str(paths["data"]), "--variants", str(paths["variants"]),
               "--upstream-dir", str(upstream), "--out-dir", str(run_dir),
               "--fixed-roster", str(cohort_dir(root, traits)/"participants.tsv.gz"),
               "--traits", ",".join(traits), "--covariates", ",".join(receipt["covariates"])]
    budget=max(1,int(cfg["runtime"]["threads"]))
    workers=min(int(options.get("workers",2)),max(0,budget-1))
    for name, default in WORKER_VALUE_OPTIONS.items():
        value=options.get(name,default)
        if name=="workers": value=workers
        if name=="threads": value=min(budget-workers,int(value))
        command += ["--" + name.replace("_", "-"), str(value)]
    if options.get("autotune", True) and stage == "train":
        command.append("--autotune")
    command += ["--seed", str(options.get("seed", cfg.get("cohort", {}).get("seed", 20261115)))]
    command += ["--gradient-checkpointing" if options.get("gradient_checkpointing", True) else "--no-gradient-checkpointing"]
    if options.get("check", False):
        command.append("--check")
    if stage == "train" and options.get("resume", False):
        command.append("--resume")
    if stage == "predict":
        checkpoint = Path(options.get("checkpoint", run_dir / "model.pt"))
        if not checkpoint.is_file():
            raise FileNotFoundError(checkpoint)
        command += ["--checkpoint", str(checkpoint.resolve())]
    return command, run_dir


def export_predictions(cfg, root, traits, run_dir):
    root, run_dir = Path(root), Path(run_dir)
    roster = identified(read_table(cohort_dir(root, traits)/"participants.tsv.gz"), "PRSformer participants")
    check_roles(roster)
    source = run_dir / "test_predictions.tsv.gz"
    predictions = read_table(source)
    if not {"eid", "trait", "prediction", "y", "split"} <= set(predictions):
        raise ValueError("Official PRSformer worker output is incomplete")
    if predictions.y.notna().any() or not predictions.split.eq("test").all():
        raise ValueError("PRSformer worker predictions must contain only masked inference-role outcomes")
    if set(predictions.trait) != set(traits):
        raise ValueError("PRSformer output task set differs")
    held = roster.role.isin(["calibration", "test"])
    expected = set(roster.loc[held, "id"])
    outputs = {}
    for trait in traits:
        current = identified(predictions.loc[predictions.trait.eq(trait)], f"PRSformer {trait} predictions")
        if set(current.id) != expected:
            raise ValueError(f"PRSformer {trait} predictions must cover exactly CALIBRATION and TEST IDs")
        values = pd.to_numeric(current.prediction, errors="raise").to_numpy(float)
        if not np.isfinite(values).all() or (TRAITS[trait] == "binary" and not ((values >= 0) & (values <= 1)).all()):
            raise ValueError(f"Invalid PRSformer {trait} predictions")
        current["prsformer.prediction"] = values
        table = aligned(current[["id", "prsformer.prediction"]], roster.id, trait, allow_missing=True)
        path = project_path(root, f"tmp/{trait}/prsformer.predictions.tsv.gz")
        atomic_table(table, path)
        outputs[trait] = str(path)
    atomic_json({"status": "predicted", "method": "official PRSformer architecture newly trained on GRID",
                 "source": stamp(source), "roles": "TRAIN/TUNE predictions absent; CALIBRATION/TEST inference only",
                 "predictions": outputs}, module_root(root,"prsformer")/"qc"/("export-"+"-".join(traits)+".json"))
    return outputs


@contextmanager
def adapter_lock(root, trait=None):
    """Share the repository's cross-node OFD locking convention."""
    import fcntl
    from file_lock import flock
    folder = project_path(root, "tmp/prsformer")
    if trait is not None:
        folder = folder / trait
    path = folder / ".adapter.lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as handle:
        try:
            flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another PRSformer adapter is using this analysis root") from exc
        yield


def archive_run(root, run_dir, cache, label):
    """Preserve an explicitly replaced model run, without moving its data cache."""
    root, run_dir, cache = Path(root).resolve(), Path(run_dir).resolve(), Path(cache).resolve()
    if run_dir == root or run_dir in root.parents or run_dir == cache or run_dir in cache.parents:
        raise ValueError("PRSformer run_dir must be separate from the analysis root and genotype cache before force_retrain")
    if not run_dir.is_dir() or not any(run_dir.iterdir()):
        return None
    import fcntl
    from file_lock import flock
    # Contend with a directly launched model worker as well as this adapter.
    # Keep this inode intact when the directory is renamed.
    with (run_dir / ".output.lock").open("a+") as handle:
        try:
            flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Cannot archive a PRSformer run while its model worker is active") from exc
        from workflow_layout import new_run_dir
        archive=new_run_dir(module_root(root,"prsformer")/"tmp"/"previous"/label)/"model"
        run_dir.rename(archive)
    atomic_json({"requested_action": "force_retrain", "previous_run": str(archive),
                 "new_run": str(run_dir), "genotype_cache_preserved": str(cache)},
                module_root(root,"prsformer")/"qc"/(label+".retrain.json"))
    print(f"PRSformer: preserved previous model run at {archive}", flush=True)
    return archive


def training_or_recovery_command(cfg, root, traits):
    """Recover completed inference or resume an interrupted original worker.

    The worker always performs checkpoint/data/upstream validation itself. A
    completed-fit marker leads to predict, never another optimizer iteration.
    An interrupted optimizer checkpoint uses its original signature guard.
    """
    _, preparation = validate_prepared(cfg, root, traits)
    command, run_dir = worker_command(cfg, root, traits, "train")
    model, training = run_dir / "model.pt", run_dir / "training.pt"
    evaluation, trained = run_dir / "evaluation.json", run_dir / "TRAINED.json"
    if evaluation.is_file() or trained.is_file():
        if not model.is_file():
            raise ValueError("PRSformer completion marker exists without model.pt; restore the model or set force_retrain")
        if evaluation.is_file():
            record = json.loads(evaluation.read_text())
            if record.get("preparation_signature") != preparation["signature"]:
                raise ValueError("Completed PRSformer model uses a different prepared cache; use a new run_dir or force_retrain")
        if trained.is_file():
            from file_state import file_record
            record = json.loads(trained.read_text())
            if record.get("model_file_record") != file_record(model):
                raise ValueError("PRSformer TRAINED.json does not match model.pt; restore the checkpoint or set force_retrain")
        context = run_dir / "adapter.json"
        if context.is_file() and json.loads(context.read_text()).get("preparation_signature") != preparation["signature"]:
            raise ValueError("PRSformer adapter and prepared cache identities differ")
        # In all/train recovery, infer from the checkpoint just completed in the
        # training directory. Explicit predict output overrides are for the
        # standalone predict stage; they must not redirect recovery elsewhere.
        recovered = deepcopy(cfg)
        recovered["prsformer"]["predict"] = {
            **_options(cfg, "train"), **_section(cfg["prsformer"], "predict"),
            "run_dir": str(run_dir), "checkpoint": str(model),
        }
        command, run_dir = worker_command(recovered, root, traits, "predict")
        print("PRSformer: completed fit found; validating the saved checkpoint and recovering inference/export", flush=True)
        return command, run_dir, "recovered_completed_fit"
    if model.is_file() != training.is_file():
        raise ValueError("Incomplete PRSformer checkpoint pair; restore model.pt/training.pt or set force_retrain")
    if model.is_file():
        if not _options(cfg, "train").get("auto_resume", True):
            raise ValueError("Interrupted PRSformer run exists; enable auto_resume or set force_retrain")
        if "--resume" not in command:
            command.append("--resume")
        print("PRSformer: resuming interrupted training with the original worker's signature checks", flush=True)
        return command, run_dir, "resumed_interrupted_fit"
    run_dir.mkdir(parents=True, exist_ok=True)
    atomic_json({"preparation_signature": preparation["signature"],
                 "traits": list(traits), "worker_command": command}, run_dir / "adapter.json")
    return command, run_dir, "new_fit"


def trait_config(cfg, root, traits):
    """Isolate single-trait caches and model runs, including custom base paths."""
    from copy import deepcopy
    cfg = deepcopy(cfg)
    if len(traits) != 1:
        return cfg
    trait = traits[0]
    base = cfg.setdefault("prsformer", {})
    folder=module_root(root,"prsformer")
    cache=folder/"tmp"/trait
    model=folder/"main"/trait
    for key, default in (("cache_dir",cache),("run_dir",model)):
        value = base.get(key)
        base[key] = str(Path(value).expanduser() / trait) if value else str(default)
        for stage in ("prepare", "train", "predict"):
            section = base.get(stage)
            if isinstance(section, dict) and section.get(key):
                section[key] = str(Path(section[key]).expanduser() / trait)
    return cfg


def run(cfg, root, traits, stage="all"):
    """Entrypoint used by grid.py; completed fits recover through real inference.

    ``force_retrain: true`` explicitly archives an existing model run before
    training afresh. It never deletes a genotype cache. ``prepare.replace`` is
    a separate decision to rebuild that potentially very large data cache.
    """
    if stage not in {"plan", "prepare", "train", "predict", "all"}:
        raise ValueError("PRSformer stage must be plan, prepare, train, predict, or all")
    if not cfg.get("prsformer", {}).get("enabled", False):
        print("PRSformer disabled; set prsformer.enabled: true to opt in", flush=True)
        return {"status": "disabled"}
    traits = _traits(traits)
    root = Path(root)
    cfg = trait_config(cfg, root, traits)
    if stage == "plan":
        return prepare(cfg, root, traits, plan_only=True)
    with adapter_lock(root, traits[0] if len(traits) == 1 else None):
        if stage == "prepare":
            return prepare(cfg, root, traits)
        return _run_locked(cfg, root, traits, stage)


def _run_locked(cfg, root, traits, stage):
    base = cfg.get("prsformer", {})
    fitting = stage == "train" or (stage == "all" and _section(base, "train").get("enabled", True))
    archive = None
    if base.get("force_retrain", False):
        if not fitting:
            raise ValueError("force_retrain requires a train/all stage with training enabled")
        train_options = _options(cfg, "train")
        run_dir = Path(train_options.get("run_dir", project_path(root, "tmp/prsformer/run")))
        cache, _ = _prepared_paths(root, _options(cfg, "prepare"))
        archive = archive_run(root, run_dir, cache, "-".join(traits))
        cfg = deepcopy(cfg)
        cfg["prsformer"]["train"] = {**_section(base, "train"), "resume": False}
    if stage == "all" and _section(base, "prepare").get("enabled", True):
        prepare(cfg, root, traits)
    ran_fit_worker = False
    run_dir = None
    disposition = "prediction"
    if fitting:
        command, run_dir, disposition = training_or_recovery_command(cfg, root, traits)
        subprocess.run(command, check=True)
        if "--check" in command:
            return {"status": "official_model_preflight_passed", "full_training_performed": False}
        ran_fit_worker = True
    if stage == "predict" or (stage == "all" and not ran_fit_worker and _section(base, "predict").get("enabled", True)):
        command, run_dir = worker_command(cfg, root, traits, "predict")
        subprocess.run(command, check=True)
        if "--check" in command:
            return {"status": "official_model_preflight_passed", "full_training_performed": False}
    if run_dir is not None:
        return {"status": "predicted", "run_action": disposition,
                "archived_run": str(archive) if archive else None,
                "outputs": export_predictions(cfg, root, traits, run_dir)}
    return {"status": "prepared_without_training"}


def main(argv=None):
    from core2 import load_config
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--root", required=True)
    parser.add_argument("--traits", help="Comma-separated subset; default config.traits")
    parser.add_argument("--stage", choices=["plan", "prepare", "train", "predict", "all"], default="all")
    args = parser.parse_args(argv)
    cfg = load_config(args.config)
    traits = args.traits.split(",") if args.traits else cfg.get("traits", list(TRAITS))
    result = run(cfg, args.root, traits, args.stage)
    print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
