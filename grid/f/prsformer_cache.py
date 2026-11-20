#!/usr/bin/env python3
"""Resolve only independent benchmark outputs; joint layouts are never reused."""
from __future__ import annotations

# benchmark_paths
import json
from pathlib import Path

TRAITS = ("height", "ldl", "t2dm")


def prsformer_directory(root, trait):
    return Path(root) / "tmp/prsformer" / trait


def prsformer_publication(root, trait):
    base = prsformer_directory(root, trait)
    return base / "split.rds", base / "predictions.rds"


def prsformer2_current(root, trait):
    base = Path(root) / "tmp/prsformer2" / trait
    pointer = base / "CURRENT.json"
    if not pointer.is_file():
        return None
    name = json.loads(pointer.read_text())["run"]
    if not isinstance(name, str) or name in ("", ".", "..") or Path(name).name != name:
        raise ValueError("PRSformer2 current run must stay inside its trait directory")
    run = base / name
    if run.is_symlink():
        raise ValueError("PRSformer2 current run cannot link to another model directory")
    return run


def benchmark_runs(root):
    """Include legacy and all independent runs when protecting unfinished work."""
    root = Path(root)
    runs = [root / "tmp/prsformer/run"]
    runs.extend(root / "tmp/prsformer" / trait / "run" for trait in TRAITS)
    runs.extend((root / "tmp/prsformer2").glob("*/_LOCAL_ONLY"))
    runs.extend((root / "tmp/prsformer2").glob("*/*/_LOCAL_ONLY"))
    return runs

# prsformer_cleanup
"""Release evaluated dosage caches while preserving verifiable report inputs."""
from file_lock import flock as shared_flock
import argparse
import contextlib
import fcntl
import hashlib
import os


def restore_released(cache, work):
    """Install an exact rebuilt matrix while keeping frozen metadata bytes."""
    import pandas as pd
    from resume import file_record
    cache, work = Path(cache), Path(work)
    marker = cache / 'genotypes.released.json'
    released = json.loads(marker.read_text()) if marker.is_file() else {}
    if released.get('schema') != 'grid.dosage_release.1':
        return False
    if released['prepare_file_record'] != file_record(cache / 'prepare.json'):
        raise ValueError('Released preparation metadata changed')
    for name in ('data.tsv.gz', 'variants.tsv.gz'):
        pd.testing.assert_frame_equal(pd.read_csv(work / name, sep='\t', dtype=str),
                                      pd.read_csv(cache / name, sep='\t', dtype=str))
    if file_record(work / 'genotypes.npy') != released['genotype_file_record']:
        raise ValueError('Rebuilt dosage differs from released original; use a new cache')
    stamp = released['genotype_mtime_ns']
    os.utime(work / 'genotypes.npy', ns=(stamp, stamp))
    os.replace(work / 'genotypes.npy', cache / 'genotypes.npy')
    return True


def validate_release(cache, evaluation):
    from resume import file_record
    receipt = json.loads((Path(cache) / 'genotypes.released.json').read_text())
    evaluation = Path(evaluation).resolve()
    if receipt.get('schema') == 'grid.dosage_release.1':
        if (receipt.get('prepare_file_record') != file_record(Path(cache) / 'prepare.json') or
                receipt.get('evaluations', {}).get(str(evaluation)) != file_record(evaluation)):
            raise ValueError('Released genotype receipt does not match this evaluation/cache')
        return
    if (receipt.get('evaluation') != str(evaluation)
            or receipt.get('evaluation_file_record') != file_record(evaluation)
            or receipt.get('prepare_file_record') != file_record(Path(cache) / 'prepare.json')):
        raise ValueError('Released genotype receipt does not match this evaluation/cache')


@contextlib.contextmanager
def cache_lock(cache):
    with (Path(cache).resolve() / '.pipeline.lock').open('a+') as handle:
        shared_flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def release(run, cache):
    """Caller holds cache and output locks; validate everything before unlinking."""
    from prsformer import validate_evaluation
    from resume import atomic_json, file_record
    run, cache = Path(run).resolve(), Path(cache).resolve()
    validate_evaluation(run, cache)
    target = cache / 'genotypes.npy'
    if not target.exists():
        return 0
    if target.is_symlink():
        raise ValueError('Refusing to release a symlinked genotype cache')
    # ABM also consumes this shared matrix, outside the PRSformer lock.
    # Retain it whenever another process has an open/mapped reference.
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit() or int(proc.name) == os.getpid():
            continue
        try:
            if proc.stat().st_uid != os.getuid():
                continue
            if (proc/'stat').read_text().rsplit(')', 1)[1].split()[0] == 'Z':
                continue
            if (proc/'comm').read_text().strip() == '(sd-pam)':
                continue
            busy = str(target) in (proc/'maps').read_text()
            for fd in (proc/'fd').iterdir():
                try:
                    busy = busy or os.readlink(fd) == str(target)
                except FileNotFoundError:
                    pass
            if busy:
                print(f'CLEANUP: genotype cache still used by PID {proc.name}; retained', flush=True)
                return 0
        except (FileNotFoundError, ProcessLookupError):
            continue
        except PermissionError:
            print(f'CLEANUP: cannot inspect PID {proc.name}; genotype cache retained', flush=True)
            return 0
    size = target.stat().st_size
    atomic_json(cache / 'genotypes.released.json', {
        'evaluation': str(run / 'evaluation.json'),
        'evaluation_file_record': file_record(run / 'evaluation.json'),
        'prepare_file_record': file_record(cache / 'prepare.json'),
        'released_bytes': size,
    })
    target.unlink()
    print(f'CLEANUP: released {size / 2**30:.2f} GiB: {target}', flush=True)
    return size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--cache-dir', required=True)
    parser.add_argument('--cache-lock-held', action='store_true', help='Internal shell wrapper only')
    args = parser.parse_args()
    from prsformer import output_lock
    with contextlib.ExitStack() as stack:
        if not args.cache_lock_held:
            stack.enter_context(cache_lock(args.cache_dir))
        stack.enter_context(output_lock(args.run_dir))
        release(args.run_dir, args.cache_dir)



if __name__ == '__main__':
    main()
