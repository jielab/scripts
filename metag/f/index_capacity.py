"""Cache the unchanged Kraken2 capacity scan, optionally ahead of index building.

Precomputation waits for each frozen library in its configured order. Its output
is reusable only with the same executable, arguments, file order and identities.
No estimates from subsets are extrapolated or added together.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from reference_batches import identity
from taxonomy import atomic_json

SCHEMA = 'metag-capacity-cache-v1'


def configuration(db):
    raw = (db / 'metag_build_state/config.txt').read_bytes()
    rows = [line.partition('=')[2].split() for line in raw.decode().splitlines()
            if line.startswith('libraries=')]
    if len(rows) != 1 or not rows[0] or len(set(rows[0])) != len(rows[0]):
        raise RuntimeError('Missing/ambiguous library configuration')
    return rows[0], hashlib.sha256(raw).hexdigest()


def executable_signature(command):
    executable = Path(shutil.which(command[0]) or command[0]).resolve(strict=True)
    return dict(executable=str(executable),
                sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                arguments=command[1:])


def snapshot(paths):
    return [dict(path=str(p.resolve()), identity=identity(p)) for p in paths]


def valid_cache(saved, parameters, rows):
    return (saved.get('schema') == SCHEMA and saved.get('parameters') == parameters
            and saved.get('inputs') == rows and type(saved.get('estimate')) is int
            and saved['estimate'] > 0)


def capacity(db, command, *, precompute=False, interval=10):
    from index_inputs import one_library, merged_inputs
    db = Path(db).resolve()
    libraries, config_sha = configuration(db)
    parameters = dict(config_sha256=config_sha, **executable_signature(command))
    folder = db / 'metag_build_state'
    receipt = folder / 'capacity.cache.json'
    with (folder / 'capacity.cache.lock').open('a') as lock:
        print('[metag] START reference/capacity: validated full-input scan/cache', file=sys.stderr, flush=True)
        fcntl.flock(lock, fcntl.LOCK_EX)
        # A foreground build must already have every validated library.
        if not precompute:
            paths, _ = merged_inputs(db)
            rows = snapshot(paths)
            try:
                saved = json.loads(receipt.read_text())
            except (OSError, ValueError):
                saved = {}
            if valid_cache(saved, parameters, rows):
                print('[metag] SKIP reference/capacity: reuse identical complete scan', file=sys.stderr, flush=True)
                return saved['estimate']
        paths, rows = [], []
        total = 0
        started = reported = time.monotonic()
        output = folder / 'capacity.stdout.partial'
        process = None
        try:
            with output.open('wb') as result:
                process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=result)
                for library in libraries:
                    marker = folder / f'{library}.complete'
                    while not marker.exists():
                        if not precompute:
                            raise RuntimeError(f'Library is incomplete: {library}')
                        if process.poll() is not None:
                            raise RuntimeError('Capacity estimator exited before all libraries were ready')
                        if time.monotonic() - reported >= 60:
                            print(f'[capacity] waiting for {library}; scanned={total / 2**30:.1f} GiB', file=sys.stderr, flush=True)
                            reported = time.monotonic()
                        time.sleep(interval)
                    inputs = one_library(db, library)
                    for path in inputs:
                        before = identity(path)
                        with path.open('rb') as source:
                            while block := source.read(8 << 20):
                                process.stdin.write(block)
                                total += len(block)
                                now = time.monotonic()
                                if now - reported >= 60:
                                    print(f'[capacity] library={library}; scanned={total / 2**30:.1f} GiB; '
                                          f'elapsed={(now-started)/60:.1f} min', file=sys.stderr, flush=True)
                                    reported = now
                        if identity(path) != before:
                            raise RuntimeError(f'Capacity input changed while reading: {path}')
                        paths.append(path)
                        rows.append(dict(path=str(path.resolve()), identity=before))
                process.stdin.close()
                if process.wait() != 0:
                    raise RuntimeError('Capacity estimator failed; no cache published')
            actual, _ = merged_inputs(db)
            if (configuration(db)[1] != config_sha or snapshot(actual) != rows
                    or executable_signature(command) != {k: v for k, v in parameters.items() if k != 'config_sha256'}):
                raise RuntimeError('Capacity inputs/configuration/executable changed; cache rejected')
            value = int(output.read_text().strip())
            if value <= 0:
                raise RuntimeError('Invalid capacity estimate')
            atomic_json(receipt, dict(schema=SCHEMA, parameters=parameters, inputs=rows,
                                     estimate=value, bytes=total, seconds=time.monotonic()-started))
            print(f'[metag] DONE reference/capacity: {total / 2**30:.2f} GiB; raw estimate={value}', file=sys.stderr, flush=True)
            return value
        finally:
            if process is not None and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            if process is not None and process.stdin is not None and not process.stdin.closed:
                try:
                    process.stdin.close()
                except (BrokenPipeError, OSError):
                    pass
            output.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, required=True)
    parser.add_argument('--precompute', action='store_true')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('An unchanged estimate_capacity command is required')
    try:
        print(capacity(args.db, command, precompute=args.precompute), flush=True)
    except (OSError, ValueError, RuntimeError) as error:
        print(f'[capacity] failed: {error}', file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
