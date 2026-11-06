"""Release evaluated dosage caches while preserving verifiable report inputs."""
import argparse
import contextlib
import fcntl
import hashlib
import json
import os
from pathlib import Path


def restore_released(cache, work):
    """Install an exact rebuilt matrix while keeping frozen metadata bytes."""
    import pandas as pd
    from resume import sha
    cache, work = Path(cache), Path(work)
    marker = cache / 'genotypes.released.json'
    released = json.loads(marker.read_text()) if marker.is_file() else {}
    if released.get('schema') != 'grid.dosage_release.1':
        return False
    if released['prepare_sha256'] != sha(cache / 'prepare.json'):
        raise ValueError('Released preparation metadata changed')
    for name in ('data.tsv.gz', 'variants.tsv.gz'):
        pd.testing.assert_frame_equal(pd.read_csv(work / name, sep='\t', dtype=str),
                                      pd.read_csv(cache / name, sep='\t', dtype=str))
    if sha(work / 'genotypes.npy') != released['genotype_sha256']:
        raise ValueError('Rebuilt dosage differs from released original; use a new cache')
    stamp = released['genotype_mtime_ns']
    os.utime(work / 'genotypes.npy', ns=(stamp, stamp))
    os.replace(work / 'genotypes.npy', cache / 'genotypes.npy')
    return True


def validate_release(cache, evaluation):
    from resume import sha
    receipt = json.loads((Path(cache) / 'genotypes.released.json').read_text())
    evaluation = Path(evaluation).resolve()
    if receipt.get('schema') == 'grid.dosage_release.1':
        if (receipt.get('prepare_sha256') != sha(Path(cache) / 'prepare.json') or
                receipt.get('evaluations', {}).get(str(evaluation)) != sha(evaluation)):
            raise ValueError('Released genotype receipt does not match this evaluation/cache')
        return
    if (receipt.get('evaluation') != str(evaluation)
            or receipt.get('evaluation_sha256') != sha(evaluation)
            or receipt.get('prepare_sha256') != sha(Path(cache) / 'prepare.json')):
        raise ValueError('Released genotype receipt does not match this evaluation/cache')


@contextlib.contextmanager
def cache_lock(cache):
    with (Path(cache).resolve() / '.pipeline.lock').open('a+') as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def release(run, cache):
    """Caller holds cache and output locks; validate everything before unlinking."""
    from prsformer import validate_evaluation
    from resume import atomic_json, sha
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
        'evaluation_sha256': sha(run / 'evaluation.json'),
        'prepare_sha256': sha(cache / 'prepare.json'),
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
