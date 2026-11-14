#!/usr/bin/env python3
"""Independent-cohort validation against a frozen primary background.

No external latent model is fitted. Primary, observed external and matched-null
profiles are projected onto the same primary-only background and genomic-distance
nuisance features on exactly the same target subset. Empirical contrasts remain
conditional and exploratory.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

from graph_model import load_matrix, atom_json, checksum, log
from profile_utils import normalise_cohort, project_profile


SCHEMA = "phole-validation-v3"
OUTPUTS = ("replication.tsv", "validated-candidates.tsv", "summary.json", "null-contrasts.tsv.gz")


def output_receipt(out):
    receipt = {}
    for name in OUTPUTS:
        path = Path(out) / name
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing/empty validation output: {path}")
        if name.endswith(".json"):
            read_json(path)
        else:
            opener = gzip.open if name.endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8", newline="") as stream:
                reader = csv.reader(stream, delimiter="\t")
                header = next(reader)
                if not header or len(set(header)) != len(header):
                    raise ValueError(f"Invalid validation table header: {path}")
                rows = 0
                for row in reader:
                    if len(row) != len(header):
                        raise ValueError(f"Invalid validation table row: {path}")
                    rows += 1
            receipt[name] = dict(columns=header, rows=rows)
        receipt.setdefault(name, {}).update(size=path.stat().st_size, sha256=checksum(path))
    return receipt

DATASETS = ("ukb", "ckb", "decode", "scot")
COLUMNS = [
    "CANDIDATE", "DATASET", "N_PROTEINS", "N_COMMON_TARGETS", "PROFILE_CORR",
    "SIGN_AGREEMENT", "STANDARDIZED_BETA_CORR", "EMPIRICAL_P", "EMPIRICAL_Q",
    "N_MATCHED_NULL", "MIN_ATTAINABLE_P", "MAF_MATCH", "DEGREE_MATCH",
    "MEDIAN_NULL_DEGREE_DISTANCE", "MIN_NULL_COMMON_TARGETS",
    "PAIRED_OBS_CORR_MIN", "PAIRED_OBS_CORR_MAX", "PANEL_SCOPE", "STATUS", "REASON",
]
NULL_COLUMNS = [
    "CANDIDATE", "DATASET", "NULL_SOURCE_INDEX", "NULL_LOCUS", "NULL_BLOCK",
    "N_COMMON_TARGETS", "N_CHILDREN", "PAIRED_OBS_CORR", "NULL_CORR",
    "NULL_GE_OBS", "DEGREE_DISTANCE",
]
PANEL_SCOPE = (
    "frozen primary-only panel, child set and background; paired common-target "
    "and genomic-distance projections; conditional exploratory empirical contrast"
)


def bh(values):
    p = np.asarray(values, float)
    q = np.full(len(p), np.nan)
    ids = np.flatnonzero(np.isfinite(p))
    if not len(ids):
        return q
    ids = ids[np.argsort(p[ids])]
    adjusted = p[ids] * len(ids) / np.arange(1, len(ids) + 1)
    q[ids] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return q


def corr(x, y):
    good = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[good], np.asarray(y)[good]
    if len(x) < 3 or np.std(x) <= 1e-12 or np.std(y) <= 1e-12:
        return None
    value = float(np.corrcoef(x, y)[0, 1])
    return float(np.clip(value, -1, 1)) if np.isfinite(value) else None


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def strings(array):
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in array.tolist()]


def atomic_frame(frame, path):
    temporary = Path(str(path) + ".tmp")
    frame.to_csv(temporary, sep="\t", index=False, na_rep=".")
    os.replace(temporary, path)


def record(candidate, dataset):
    return dict(CANDIDATE=str(candidate), DATASET=dataset, N_PROTEINS=0,
                N_COMMON_TARGETS=0, PROFILE_CORR=None, SIGN_AGREEMENT=None,
                STANDARDIZED_BETA_CORR=None, EMPIRICAL_P=None, EMPIRICAL_Q=None,
                N_MATCHED_NULL=0, MIN_ATTAINABLE_P=None, MAF_MATCH="not_evaluated",
                DEGREE_MATCH="not_evaluated", MEDIAN_NULL_DEGREE_DISTANCE=None,
                MIN_NULL_COMMON_TARGETS=None, PAIRED_OBS_CORR_MIN=None,
                PAIRED_OBS_CORR_MAX=None, PANEL_SCOPE=PANEL_SCOPE,
                STATUS="untestable", REASON="")


def frozen_programs(root, primary, model_dir, discovery_dir):
    candidate_path = discovery_dir / "candidates.tsv"
    program_path = discovery_dir / "candidate_programs.npz"
    candidates = pd.read_csv(candidate_path, sep="\t", keep_default_na=False,
                             dtype={"CANDIDATE": str, "LOCUS": str, "CHR": str, "LD_BLOCK": str})
    required_table = {"CANDIDATE", "SOURCE_INDEX", "LOCUS", "LD_BLOCK"}
    if not required_table.issubset(candidates):
        raise ValueError("candidates.tsv is missing required identity columns")
    if candidates.CANDIDATE.duplicated().any():
        raise ValueError("Duplicate candidate IDs in candidates.tsv")
    dm = read_json(discovery_dir / "discovery.manifest.json")
    fm = read_json(Path(model_dir) / "fit.manifest.json")
    mm = read_json(root / primary / "manifest.json")
    if dm.get("schema") != "phole-discovery-v2" or fm.get("schema") != "phole-signed-fit-v2":
        raise ValueError("Signed v2 discovery and fit manifests are required")
    if mm.get("schema", mm.get("format")) != "phole-signed-matrix-v2" or not mm.get("directional") or mm.get("input_scope") != "full":
        raise ValueError("Primary input must be a signed full-summary v2 matrix")
    if dm.get("source_hash") != mm.get("source_hash"):
        raise ValueError("Primary matrix differs from the frozen discovery source panel")
    if checksum(program_path) != dm.get("program_signature"):
        raise ValueError("Discovery programs changed after freezing")
    if checksum(candidate_path) != dm.get("candidate_table_signature"):
        raise ValueError("Candidate table changed after freezing, or discovery needs regeneration")
    source_path, target_path = root / primary / "source.tsv", root / primary / "target.tsv"
    source_hash, target_hash = checksum(source_path), checksum(target_path)
    if source_hash != dm.get("source_hash") or source_hash != dm.get("source_signature", source_hash):
        raise ValueError("Primary source.tsv differs from frozen discovery")
    if target_hash != dm.get("target_signature", target_hash):
        raise ValueError("Primary target.tsv differs from frozen discovery")
    source = pd.read_csv(source_path, sep="\t", keep_default_na=False, dtype={"CHR": str, "LD_BLOCK": str, "LOCUS": str})
    target = pd.read_csv(target_path, sep="\t", keep_default_na=False, dtype={"PROT": str, "ASSAY": str})
    required = {"candidates", "source_index", "proteins", "primary_values", "primary_se",
                "background_basis", "geometry_profiles", "children", "excluded_blocks", "min_abs_z"}
    with np.load(program_path, allow_pickle=False) as stored:
        if not required.issubset(stored.files):
            raise ValueError("Frozen-background program fields are absent; regenerate discover before validate")
        p = {name: stored[name].copy() for name in required}
    ids, names = strings(p["candidates"]), strings(p["proteins"])
    if len(ids) != len(set(ids)) or len(names) != len(set(names)):
        raise ValueError("Frozen candidate or protein IDs are duplicated")
    wanted = candidates.CANDIDATE.tolist()
    if set(ids) != set(wanted):
        raise ValueError("NPZ candidate IDs and candidates.tsv differ")
    if set(names) != set(target.PROT.astype(str)):
        raise ValueError("Frozen protein names differ from primary target metadata")
    lookup = {name: i for i, name in enumerate(ids)}
    order = np.array([lookup[name] for name in wanted], dtype=int)
    shape = (len(ids), len(names))
    for name in ("primary_values", "primary_se", "children"):
        if p[name].shape != shape:
            raise ValueError(f"Frozen {name} has inconsistent dimensions")
        p[name] = p[name][order]
    if p["source_index"].shape != (len(ids),) or not np.issubdtype(p["source_index"].dtype, np.integer):
        raise ValueError("NPZ source_index must be one integer per candidate")
    npz_index = p["source_index"][order].astype(int)
    table_index = pd.to_numeric(candidates.SOURCE_INDEX, errors="raise").to_numpy()
    if np.any(~np.isfinite(table_index)) or np.any(table_index != table_index.astype(int)):
        raise ValueError("Candidate SOURCE_INDEX is not integral")
    table_index = table_index.astype(int)
    if not np.array_equal(npz_index, table_index):
        raise ValueError("NPZ source_index differs from candidate table")
    if np.any(table_index < 0) or np.any(table_index >= len(source)):
        raise ValueError("Candidate source index is outside the primary matrix")
    for k, row in candidates.iterrows():
        actual = source.iloc[table_index[k]]
        if str(actual.LOCUS) != str(row.LOCUS) or str(actual.LD_BLOCK) != str(row.LD_BLOCK):
            raise ValueError("Candidate locus or block identity differs from source metadata")
    basis = np.asarray(p["background_basis"], float)
    if basis.ndim != 2 or basis.shape[0] != len(names) or basis.shape[1] < 1:
        raise ValueError("Frozen background basis has inconsistent dimensions")
    if not np.all(np.isfinite(basis)) or not np.allclose(basis[:, 0], 1):
        raise ValueError("Frozen background basis must be finite and start with an intercept")
    geometry = p["geometry_profiles"]
    if geometry.shape != (len(source), len(names), 4) or not np.issubdtype(geometry.dtype, np.number):
        raise ValueError("Frozen genomic geometry has inconsistent dimensions or dtype")
    for first in range(0, len(source), 256):
        if not np.all(np.isfinite(geometry[first:first + 256])):
            raise ValueError("Frozen genomic geometry must be finite")
    if p["children"].dtype.kind != "b":
        raise ValueError("Frozen children must be a Boolean matrix")
    if np.any(p["children"] & (~np.isfinite(p["primary_values"]) |
                              ~np.isfinite(p["primary_se"]) | (p["primary_se"] <= 0))):
        raise ValueError("Frozen child has no valid primary effect/SE")
    excluded = set(strings(p["excluded_blocks"])) | set(candidates.LD_BLOCK.astype(str))
    if not excluded.issubset(set(source.LD_BLOCK.astype(str))):
        raise ValueError("An excluded discovery block is absent from the source panel")
    min_abs_z = float(p["min_abs_z"].item())
    if not np.isfinite(min_abs_z) or min_abs_z <= 0:
        raise ValueError("Invalid frozen discovery screening threshold")
    return candidates, p, names, basis, excluded, min_abs_z, dm, fm, mm


def cohort_maf(folder, shape):
    path = folder / "eaf.npy"
    result = np.full(shape[0], np.nan)
    if not path.is_file():
        return result
    a = np.load(path, mmap_mode="r", allow_pickle=False)
    if a.shape != shape:
        raise ValueError("eaf.npy and beta matrix shapes differ")
    for i in range(shape[0]):
        v = np.asarray(a[i])
        v = v[np.isfinite(v) & (v >= 0) & (v <= 1)]
        if len(v):
            f = float(np.median(v))
            result[i] = min(f, 1 - f)
    return result


def input_signature(root, model_dir, discovery_dir, primary, datasets, parameters):
    files = [Path(__file__), Path(__file__).with_name("profile_utils.py"),
             Path(__file__).with_name("graph_model.py"),
             Path(model_dir) / "fit.manifest.json",
             discovery_dir / "discovery.manifest.json", discovery_dir / "candidates.tsv",
             discovery_dir / "candidate_programs.npz"]
    for d in (primary, *datasets):
        folder = root / d
        files.extend(folder / f for f in ("manifest.json", "source.tsv", "target.tsv",
                     "beta.npy", "se.npy", "observed.npy", "cis.npy", "target_cis_known.npy", "eaf.npy"))
    identities = []
    for path in files:
        if path.is_file():
            stat = path.stat()
            identity = dict(path=str(path.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns)
            if path.suffix != ".npy":
                identity["sha256"] = checksum(path)
            identities.append(identity)
        else:
            identities.append(dict(path=str(path.resolve()), missing=True))
    signature = hashlib.sha256(json.dumps([parameters, identities], sort_keys=True).encode()).hexdigest()
    return signature, identities


def validate(matrix_root, model_dir, discovery_dir, out, primary="ukb",
             datasets=("ckb", "decode", "scot"), min_targets=10, min_null=30,
             max_null=1000, seed=20261105, replace=False):
    if primary not in DATASETS or min_targets < 3 or min_null < 10 or max_null < min_null:
        raise ValueError("Invalid primary or replication thresholds")
    datasets = tuple(dict.fromkeys(d for d in datasets if d and d != primary))
    if any(d not in DATASETS for d in datasets):
        raise ValueError("Unknown validation dataset")
    root, out, discovery_dir = Path(matrix_root), Path(out), Path(discovery_dir)
    out.mkdir(parents=True, exist_ok=True)
    candidates, programs, ppnames, full_basis, excluded, min_abs_z, dm, fm, mm = frozen_programs(
        root, primary, model_dir, discovery_dir)
    parameters = dict(primary=primary, datasets=datasets, min_targets=min_targets, min_null=min_null,
                      max_null=max_null, seed=seed, maf_tolerance=.10, common_target_fraction=.80,
                      minimum_profile_corr=.30, minimum_sign_agreement=.60, fdr=.05)
    signature, identities = input_signature(root, model_dir, discovery_dir, primary, datasets, parameters)
    old_path = out / "validation.manifest.json"
    if not replace and old_path.is_file():
        old = read_json(old_path)
        if old.get("signature") == signature:
            if old.get("outputs") != output_receipt(out):
                raise ValueError("Validation output receipt missing or changed; inspect outputs and use --replace to recompute")
            log("Reuse matching frozen-background validation " + str(out))
            return read_json(out / "summary.json")
    (out / "validate.done").unlink(missing_ok=True)
    old_path.unlink(missing_ok=True)
    rows, normalisation_summary = [], {}
    null_path = out / "null-contrasts.tsv.gz"
    null_tmp = Path(str(null_path) + ".tmp")
    try:
        with gzip.open(null_tmp, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, NULL_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for d in datasets:
                if not len(candidates):
                    continue
                folder = root / d
                if not (folder / "manifest.json").is_file():
                    for candidate in candidates.CANDIDATE:
                        rec = record(candidate, d)
                        rec["REASON"] = "full summary matrix unavailable; missing is not failed replication"
                        rows.append(rec)
                    continue
                beta, se, observed, cis, source, target, meta = load_matrix(folder)
                if meta.get("input_scope") != "full" or meta.get("source_hash") != dm.get("source_hash"):
                    raise ValueError(f"{d}: full-summary source panel differs from frozen discovery")
                if checksum(folder / "source.tsv") != dm["source_hash"]:
                    raise ValueError(f"{d}: source.tsv content differs from frozen discovery")
                normal = normalise_cohort(beta, se, observed, cis, source, target, excluded)
                scales = np.asarray(normal["scales"])
                if np.any(np.isfinite(scales) & (scales <= 0)):
                    raise ValueError("Cohort normalisation returned a nonpositive scale")
                values, errors, names = normal["values"], normal["errors"], list(normal["proteins"])
                if len(names) != len(set(names)) or values.shape != errors.shape or values.shape != (len(source), len(names)):
                    raise ValueError("Normalised protein profile dimensions or identities are invalid")
                external_name_set = set(names)
                shared = [name for name in ppnames if name in external_name_set]
                pmap, emap = {name: i for i, name in enumerate(ppnames)}, {name: i for i, name in enumerate(names)}
                pi, ei = np.array([pmap[n] for n in shared], int), np.array([emap[n] for n in shared], int)
                basis = full_basis[pi]
                ev = np.asarray(values[:, ei], float)
                normalisation_summary[d] = dict(
                    shared_proteins=len(shared), usable_assays=int(np.sum(normal.get("usable", np.isfinite(scales)))),
                    excluded_blocks=len(excluded), method="positive scale and background median from noncandidate blocks only",
                    latent_refit=False, negative_sign_flips=False)
                maf = cohort_maf(folder, beta.shape)
                ratio = np.divide(values, errors, out=np.full(values.shape, np.nan),
                                  where=np.isfinite(errors) & (errors > 0))
                degree = np.sum(np.isfinite(ratio) & (np.abs(ratio) >= min_abs_z), axis=1)
                del ratio
                blocks = source.LD_BLOCK.astype(str).to_numpy()
                chromosomes = source.CHR.astype(str).to_numpy()
                chrclass = np.array([c if c in ("X", "Y", "MT") else "autosome" for c in chromosomes])
                noncandidate = ~np.isin(blocks, list(excluded))
                for k, candidate in candidates.iterrows():
                    idx = int(candidate.SOURCE_INDEX)
                    if idx >= len(source) or str(source.iloc[idx].LOCUS) != str(candidate.LOCUS):
                        raise ValueError("External candidate row identity differs from frozen source panel")
                    rec = record(candidate.CANDIDATE, d)
                    x, y = programs["primary_values"][k, pi], ev[idx]
                    children = programs["children"][k, pi]
                    # Remove conservative genomic-distance nuisance structure even
                    # when its coefficient was not selected by the fitted model.
                    candidate_basis = np.column_stack((basis, programs["geometry_profiles"][idx, pi, :]))
                    common = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(candidate_basis), axis=1)
                    n_common, n_children = int(common.sum()), int(np.sum(children & common))
                    rec["N_COMMON_TARGETS"], rec["N_PROTEINS"] = n_common, n_children
                    if n_children < min_targets:
                        rec["REASON"] = "too few observed frozen child proteins; missing is not failed replication"
                        rows.append(rec); continue
                    rx = project_profile(x, candidate_basis, common)
                    ry = project_profile(y, candidate_basis, common)
                    statistic = corr(rx[common], ry[common])
                    rec["PROFILE_CORR"] = statistic
                    rec["STANDARDIZED_BETA_CORR"] = corr(x[common], y[common])
                    child_common = children & common & np.isfinite(rx) & np.isfinite(ry)
                    if statistic is None or int(child_common.sum()) < min_targets:
                        rec["REASON"] = "common-target projection is constant or has too few residual degrees of freedom"
                        rows.append(rec); continue
                    rec["SIGN_AGREEMENT"] = float(np.mean(np.sign(rx[child_common]) == np.sign(ry[child_common])))
                    pool = np.flatnonzero(noncandidate & (chrclass == chrclass[idx]))
                    if np.isfinite(maf[idx]):
                        matched = pool[np.isfinite(maf[pool]) & (np.abs(maf[pool] - maf[idx]) <= .10)]
                        if len(set(blocks[matched])) >= min_null:
                            pool = matched; rec["MAF_MATCH"] = "within_0.10"
                        else:
                            rec["MAF_MATCH"] = "fallback_insufficient_MAF_matched_blocks"
                    else:
                        rec["MAF_MATCH"] = "fallback_candidate_MAF_unavailable"
                    rec["DEGREE_MATCH"] = "nearest_degree_first; no fixed caliper; one accepted row per block"
                    salt = int.from_bytes(hashlib.sha256((d + "|" + str(candidate.CANDIDATE)).encode()).digest()[:8], "little")
                    rng = np.random.default_rng(int(seed) ^ salt)
                    distance = np.abs(degree[pool] - degree[idx])
                    order = pool[np.lexsort((rng.random(len(pool)), distance))]
                    seen, comparisons, paired_obs, common_sizes, degree_distances = set(), [], [], [], []
                    required_common = max(3, int(math.ceil(.8 * n_common)))
                    for ri in order:
                        if blocks[ri] in seen:
                            continue
                        yy = ev[ri]
                        subset = common & np.isfinite(yy)
                        n_shared, n_child = int(subset.sum()), int(np.sum(children & subset))
                        if n_shared < required_common or n_child < min_targets:
                            continue
                        # The candidate and null have different source positions.
                        # Remove both positions' nuisance spaces from all three
                        # vectors, preserving one target set and one projector.
                        paired_basis = np.column_stack((candidate_basis, programs["geometry_profiles"][ri, pi, :]))
                        pxx = project_profile(x, paired_basis, subset)
                        pyy = project_profile(y, paired_basis, subset)
                        pnn = project_profile(yy, paired_basis, subset)
                        observed_pair = corr(pxx[subset], pyy[subset])
                        null_value = corr(pxx[subset], pnn[subset])
                        if observed_pair is None or null_value is None:
                            continue
                        exceed = bool(null_value >= observed_pair)
                        seen.add(blocks[ri]); comparisons.append(exceed); paired_obs.append(observed_pair)
                        common_sizes.append(n_shared); degree_distances.append(int(abs(degree[ri] - degree[idx])))
                        writer.writerow(dict(
                            CANDIDATE=candidate.CANDIDATE, DATASET=d, NULL_SOURCE_INDEX=int(ri),
                            NULL_LOCUS=str(source.iloc[ri].LOCUS), NULL_BLOCK=blocks[ri],
                            N_COMMON_TARGETS=n_shared, N_CHILDREN=n_child,
                            PAIRED_OBS_CORR=observed_pair, NULL_CORR=null_value,
                            NULL_GE_OBS=int(exceed), DEGREE_DISTANCE=degree_distances[-1]))
                        if len(comparisons) >= max_null:
                            break
                    count = len(comparisons)
                    rec["N_MATCHED_NULL"] = count
                    rec["MIN_ATTAINABLE_P"] = 1 / (count + 1) if count else None
                    if count:
                        rec["MEDIAN_NULL_DEGREE_DISTANCE"] = float(np.median(degree_distances))
                        rec["MIN_NULL_COMMON_TARGETS"] = min(common_sizes)
                        rec["PAIRED_OBS_CORR_MIN"], rec["PAIRED_OBS_CORR_MAX"] = min(paired_obs), max(paired_obs)
                    if count < min_null:
                        rec["REASON"] = "insufficient noncandidate block matches after paired common-target requirements"
                        rows.append(rec); continue
                    rec["EMPIRICAL_P"] = (1 + int(np.sum(comparisons))) / (count + 1)
                    rec["STATUS"] = "tested"
                    rec["REASON"] = (
                        "frozen-background and distance-adjusted paired-target whole-locus contrast; conditional exploratory inference; "
                        + rec["MAF_MATCH"] + "; genomic blocks are not LD proof")
                    rows.append(rec)
                log(f"{d}: frozen-background profile validation complete; no external latent fit")
        os.replace(null_tmp, null_path)
    finally:
        null_tmp.unlink(missing_ok=True)
    result = pd.DataFrame(rows, columns=COLUMNS)
    if len(result):
        result["EMPIRICAL_Q"] = bh(pd.to_numeric(result.EMPIRICAL_P, errors="coerce"))
        for i, row in result.iterrows():
            if row.STATUS == "tested":
                supported = row.EMPIRICAL_Q <= .05 and row.PROFILE_CORR >= .3 and row.SIGN_AGREEMENT >= .6
                result.loc[i, "STATUS"] = "profile_replicated" if supported else "not_confirmed"
    atomic_frame(result, out / "replication.tsv")
    supports = (result[result.STATUS == "profile_replicated"].groupby("CANDIDATE").DATASET.agg(list).to_dict()
                if len(result) else {})
    combined = candidates.copy()
    combined["REPLICATED_COHORTS"] = combined.CANDIDATE.map(lambda x: ";".join(supports.get(x, [])) or ".")
    combined["N_REPLICATED_COHORTS"] = combined.CANDIDATE.map(lambda x: len(supports.get(x, []))).astype(int)
    combined["INTERPRETATION"] = combined.CANDIDATE.map(
        lambda x: "replicated frozen-background unexplained profile; causal node unresolved"
        if len(supports.get(x, [])) >= 2 else "requires additional validation")
    atomic_frame(combined, out / "validated-candidates.tsv")
    tested = int(result.STATUS.isin(["profile_replicated", "not_confirmed"]).sum()) if len(result) else 0
    summary = dict(
        schema=SCHEMA, primary=primary, validation_datasets=list(datasets), tested=tested,
        replicated=int((result.STATUS == "profile_replicated").sum()) if len(result) else 0,
        untestable=int((result.STATUS == "untestable").sum()) if len(result) else 0,
        candidate_count=len(candidates), at_least_two_external=int((combined.N_REPLICATED_COHORTS >= 2).sum()),
        method="shared frozen primary background plus source-distance nuisance projection; no external latent refitting",
        profile_statistic="correlation over all shared projected targets; sign agreement and min_targets use frozen primary children",
        null_method="one whole profile per noncandidate block; paired observed/null common targets; deterministic degree order with seeded ties",
        null_exclusions="all initially nominated and retained candidate blocks",
        n_excluded_blocks=len(excluded), fdr_family="BH across tested candidate x external contrasts; conditional exploratory inference",
        minimum_attainable_p="1 / (N_MATCHED_NULL + 1); no zero P; resolution may prevent an isolated BH rejection",
        isolated_BH_first_rejection_threshold=.05 / tested if tested else None,
        maf_matching="within-cohort MAF difference <=0.10 when enough distinct blocks; fallback recorded per contrast",
        geometry="main profile uses [background,G_candidate]; paired null uses [background,G_candidate,G_null] for all three vectors on identical targets; distance features always removed as conservative nuisance; residual degrees of freedom use actual matrix rank",
        normalisation=normalisation_summary, causal_claim=False,
        synthetic=bool(mm.get("synthetic", False) or fm.get("synthetic", False) or dm.get("synthetic", False)))
    atom_json(out / "summary.json", summary)
    atom_json(out / "validation.manifest.json", dict(
        schema=SCHEMA, signature=signature, parameters=parameters, inputs=identities, outputs=output_receipt(out),
        program_signature=dm["program_signature"], candidate_table_signature=dm["candidate_table_signature"],
        source_hash=dm["source_hash"], method=summary["method"], synthetic=summary["synthetic"]))
    (out / "validate.done").write_text(signature + "\n")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix-root", required=True)
    p.add_argument("--model-dir", required=True)
    p.add_argument("--discovery-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--primary", choices=DATASETS, default="ukb")
    p.add_argument("--datasets", default="ckb,decode,scot")
    p.add_argument("--min-targets", type=int, default=10)
    p.add_argument("--min-null", type=int, default=30)
    p.add_argument("--max-null", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261105)
    p.add_argument("--replace", choices=["TRUE", "FALSE"], default="TRUE")
    a = p.parse_args()
    if a.min_targets < 3 or a.min_null < 10 or a.max_null < a.min_null:
        p.error("min-targets>=3, min-null>=10 and max-null>=min-null are required")
    ds = tuple(d.strip() for d in a.datasets.split(",") if d.strip())
    print(json.dumps(validate(a.matrix_root, a.model_dir, a.discovery_dir, a.out, a.primary, ds,
                             a.min_targets, a.min_null, a.max_null, a.seed, a.replace == "TRUE"), indent=2))


if __name__ == "__main__":
    main()
