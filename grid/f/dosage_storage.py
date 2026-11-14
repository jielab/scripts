#!/usr/bin/env python3
"""Keep the 32 real dosage rows used by ABM audits, release full staging matrices.

Frozen block receipts are never rewritten. Training can rebuild released dosage
from PGEN; unfinished benchmark checkpoints keep their original input files.
"""
from __future__ import annotations


# Shared readable path registry; fingerprints remain in QC metadata.
import sys as _layout_sys
from pathlib import Path as _LayoutPath
_layout_root = next(p for p in _LayoutPath(__file__).resolve().parents if (p / "0f" / "workflow_layout.py").is_file())
_layout_sys.path.append(str(_layout_root / "0f"))
from workflow_layout import named_directory

from file_lock import flock as shared_flock
from contextlib import contextmanager, ExitStack
import fcntl
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from resume import atomic_json, file_record

PROBE = "dosage_probe.npy"
MARKER = "DOSAGE_PROBE.json"


@contextmanager
def dosage_lock(config, *, exclusive=False, wait=False):
    root = Path(config["cache_root"])
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".dosage.lock").open("a+") as handle:
        shared_flock(handle, (fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH) | (0 if wait else fcntl.LOCK_NB))
        yield


@contextmanager
def prsformer_cache_lock(cache, *, exclusive=False, wait=True):
    """Use the shell pipeline's cache key for readers, producers and release."""
    root = Path(cache)
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".pipeline.lock").open("a+") as handle:
        shared_flock(handle, (fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH) | (0 if wait else fcntl.LOCK_NB))
        yield


@contextmanager
def dosage_supplement_lock(config, snps, *, wait=True):
    """Protect equivalent supplement writers while distinct SNP sets coexist.

    A supplement target binds its complete variant records, so equal targets
    necessarily have equal SNP sets. Ignore order, duplicates and trait names:
    traits can share the same target. The outer dosage read lock must be held
    throughout preparation and use to exclude cleanup and input replacement.
    """
    from file_lock import file_lock
    identity = json.dumps(sorted(set(map(str, snps))), separators=(",", ":"))
    key = hashlib.sha256(identity.encode()).hexdigest()
    root = Path(config["cache_root"]) / ".supplement_locks"
    with file_lock(root / (key + ".lock"), wait=wait):
        yield


def probe_rows(n):
    return np.sort(np.random.default_rng(20261020).choice(n, min(32, n), replace=False))


def probe_identity(receipt, rows):
    return dict(schema="grid.dosage_probe.1", rows=rows.tolist(),
                signature_payload=receipt["signature_payload"],
                array_binding=receipt["array_binding"],
                dosage_supplement=receipt.get("dosage_supplement"))


def save_probe(base, receipt, rows, genotype):
    """Caller has verified the full dosage hashes and independent block scoring."""
    base = Path(base)
    expected = probe_identity(receipt, rows)
    if (base / MARKER).is_file():
        load_probe(base, receipt, rows)
        return
    temporary = base / (PROBE + ".tmp")
    try:
        values = genotype[np.ix_(rows, np.arange(genotype.shape[1]))]
        with temporary.open("wb") as handle:
            np.save(handle, values, allow_pickle=False)
        os.replace(temporary, base / PROBE)
        atomic_json(base / MARKER, dict(identity=expected, file_record=file_record(base / PROBE),
                                       shape=list(values.shape)))
    finally:
        temporary.unlink(missing_ok=True)


def load_probe(base, receipt, rows):
    base = Path(base)
    if any((base / name).is_symlink() for name in (PROBE, MARKER)):
        raise ValueError("Symlinked dosage probe")
    saved = json.loads((base / MARKER).read_text())
    if saved["identity"] != probe_identity(receipt, rows) or file_record(base / PROBE) != saved["file_record"]:
        raise ValueError("Dosage probe changed or belongs to different block inputs")
    values = np.load(base / PROBE, allow_pickle=False)
    if list(values.shape) != saved["shape"] or len(values) != len(rows):
        raise ValueError("Dosage probe dimensions changed")
    return values


def _safe_file(path, root):
    path, root = Path(path), Path(root)
    if path.is_symlink() or root.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Refusing cleanup outside the declared cache: " + str(path))
    for parent in path.parents:
        if parent == root:
            break
        if parent.is_symlink():
            raise ValueError("Refusing symlinked cache ancestor: " + str(parent))


def _release(target, expected_hash, marker, details):
    from block_inputs import file_record_cached
    from abm.io import busy_reason
    if not target.is_file():
        return 0
    _safe_file(target, target.parent)
    reason = busy_reason([target])
    if reason:
        print(f"[GRID cleanup] RETAIN {target}: {reason}", flush=True)
        return 0
    if file_record_cached(target) != expected_hash:
        raise ValueError("Dosage changed before release: " + str(target))
    stat = target.stat()
    atomic_json(marker, dict(details, genotype_file_record=expected_hash,
                            released_bytes=stat.st_size, genotype_mtime_ns=stat.st_mtime_ns))
    target.unlink()
    print(f"[GRID cleanup] Released {stat.st_size / 2**30:.2f} GiB: {target}", flush=True)
    return stat.st_size


def _benchmark_evaluations(root, cache):
    """Return completed evaluations, or None if a checkpoint still needs dosage."""
    from prsformer import validate_evaluation
    evaluations = {}
    from prsformer_cache import benchmark_runs
    for run in benchmark_runs(root):
        if not any((run / name).exists() for name in ("model.pt", "training.pt", "TRAINED.json", "evaluation.json")):
            continue
        if not (run / "evaluation.json").is_file():
            print(f"[GRID cleanup] RETAIN shared dosage: unfinished benchmark {run}", flush=True)
            return None
        validate_evaluation(run, cache)
        evaluations[str((run / "evaluation.json").resolve())] = file_record(run / "evaluation.json")
    return evaluations


def release_dosages(config_path, *, traits, analysis_root, apply=False, best_effort=False):
    """Release only inputs with verified compact probes and intact block outputs.

    All configured traits must have probes before the shared PRSformer matrix is
    released, so a multi-trait build never decodes the same PGEN twice.
    """
    import block_inputs as b
    from abm.io import verify_array_binding, binding_metadata
    config = b.read_config(config_path)
    root = Path(config["cache_root"])
    cache = Path(config["prsformer_cache"])
    analysis = Path(analysis_root).resolve()
    if not root.resolve().is_relative_to(analysis) or not cache.resolve().is_relative_to(analysis):
        print("[GRID cleanup] RETAIN dosage outside --analysis-root; use matching block config", flush=True)
        return 0
    if not root.is_dir() or not cache.is_dir():
        return 0
    if not (cache / "genotypes.npy").exists() and not any((root / "dosage_supplements").glob("*/genotypes.npy")):
        return 0
    total = 0
    try:
        with ExitStack() as stack:
            stack.enter_context(dosage_lock(config, exclusive=True))
            for module in ("prsformer", "prsformer2"):
                folder = Path(analysis_root) / module
                if folder.is_dir():
                    handle = stack.enter_context((folder / ".resume.lock").open("a+"))
                    shared_flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            for lock in (cache / "grid.wrapper.lock", cache / ".pipeline.lock"):
                handle = stack.enter_context(lock.open("a+"))
                shared_flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            entries = []
            for marker in sorted((root / b.BINDING_NAMESPACE).glob("*/*/_LOCAL_ONLY/BLOCKS_COMPLETE.json")):
                _safe_file(marker, root)
                receipt = json.loads(marker.read_text())
                base = marker.parent
                if receipt.get("status") != "complete" or not (base / MARKER).is_file():
                    # A missing probe is not proof that the dosage is unused.
                    entries.append((base, receipt, False))
                    continue
                rows = probe_rows(receipt["shape"][0])
                load_probe(base, receipt, rows)
                payload = receipt["signature_payload"]
                for name, key in (("prepare.json", "prsformer_prepare"),
                                  ("data.tsv.gz", "cohort_file_record"), ("variants.tsv.gz", "variants")):
                    if file_record(cache / name) != payload[key]:
                        raise ValueError("Probe source metadata changed: " + str(cache / name))
                for pop, path in b.locate_weights(config, receipt["trait"])[0].items():
                    if file_record(path) != payload["weights"][pop]:
                        raise ValueError("Block weights changed before dosage cleanup")
                if file_record(Path(config["score_root"]) / receipt["trait"] / "1csx.scores.rds") != payload["whole_score_file_record"]:
                    raise ValueError("Full CSx scores changed before dosage cleanup")
                binding = receipt["array_binding"]
                score = base / "block_scores.npy"
                _safe_file(score, root)
                from abm.grid_inputs import table
                cohort = table(base / "cohort.tsv.gz")
                tokens = table(base / "tokens.tsv.gz")
                verify_array_binding(score, binding_metadata(cohort), binding,
                                     tokens=tokens, channels=receipt["channels"])
                entries.append((base, receipt, True))
            supplements = {}
            for base, receipt, valid in entries:
                supplement = receipt.get("dosage_supplement")
                if supplement:
                    supplements.setdefault(supplement["signature"], []).append((base, supplement, valid))
            for key, users in supplements.items():
                if not all(valid for _, _, valid in users):
                    continue
                folder = named_directory(root / "dosage_supplements", key, "supplement")
                target = folder / "genotypes.npy"
                _safe_file(target, root)
                if not target.exists():
                    continue
                saved = json.loads((folder / "COMPLETE.json").read_text())
                if any(saved != supplement for _, supplement, _ in users):
                    raise ValueError("Supplement producer receipt changed")
                if apply:
                    total += _release(target, saved["genotype_file_record"], folder / "genotypes.released.json",
                                      dict(probes=[str(base / MARKER) for base, _, _ in users],
                                           complete_file_record=file_record(folder / "COMPLETE.json")))
                else:
                    total += target.stat().st_size
            layout = b.resolve_block_files(config["ld_block_dir"], config["ld_block_files"])
            ready_traits = {r["trait"] for _, r, valid in entries if valid and r["layout"]["signature"] == layout["signature"]}
            target = cache / "genotypes.npy"
            if (entries and all(valid for _, _, valid in entries) and set(traits) <= ready_traits
                    and target.is_file()):
                expected = {r["signature_payload"]["genotype_file_record"] for _, r, _ in entries}
                if len(expected) != 1:
                    raise ValueError("Block probes refer to different shared genotype matrices")
                evaluations = _benchmark_evaluations(Path(analysis_root), cache)
                if evaluations is not None:
                    if apply:
                        total += _release(target, expected.pop(), cache / "genotypes.released.json",
                                          dict(schema="grid.dosage_release.1", prepare_file_record=file_record(cache / "prepare.json"),
                                               evaluations=evaluations,
                                               probes=[str(base / MARKER) for base, _, _ in entries]))
                    else:
                        total += target.stat().st_size
    except (OSError, ValueError, KeyError, TypeError) as exc:
        if not best_effort:
            raise
        print(f"[GRID cleanup] RETAIN remaining dosage: {exc}", flush=True)
    return total
