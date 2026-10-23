"""Verified optional PGEN execution copies on a fast local filesystem.

Scientific outputs stay in the analysis directory. Original inputs are never
modified. A complete cache requires source fingerprints and full copy hashes.
"""
from __future__ import annotations
import fcntl
import hashlib
import os
from pathlib import Path
import shutil
from core import dump_json, fingerprint, stable_hash


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        released = consumed = 0
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
            consumed += len(block)
            if consumed - released >= 256 * 1024 * 1024 and hasattr(os, 'posix_fadvise'):
                os.posix_fadvise(handle.fileno(), released, consumed - released, os.POSIX_FADV_DONTNEED)
                released = consumed
        if hasattr(os, 'posix_fadvise'):
            os.posix_fadvise(handle.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
    return h.hexdigest()


def cache_pgen_prefix(prefix, cache_dir):
    prefix, root = Path(prefix).resolve(), Path(cache_dir).expanduser().resolve()
    sources = {suffix: Path(str(prefix) + suffix) for suffix in ('.pgen', '.pvar', '.psam')}
    identity = {suffix: fingerprint(path) for suffix, path in sources.items()}
    root.mkdir(parents=True, exist_ok=True)
    destination = root / (prefix.name + '.' + stable_hash(identity)[:20])
    with (root / '.input_cache.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        manifest = destination / 'cache.json'
        if destination.exists():
            if not manifest.is_file():
                raise ValueError(f'Incomplete/unrecognized input cache: {destination}')
            saved = __import__('json').loads(manifest.read_text())
            if saved['source_fingerprints'] != identity:
                raise ValueError('PGEN cache source identity changed')
            for suffix, spec in saved['files'].items():
                path = destination / (prefix.name + suffix)
                if fingerprint(path) != spec['fingerprint'] or _sha(path) != spec['sha256']:
                    raise ValueError(f'PGEN cache checksum differs: {path}')
            return str(destination / prefix.name), saved
        staging = destination.with_name(destination.name + '.partial')
        if staging.exists():
            shutil.rmtree(staging)  # only this helper's exact incomplete cache
        required = sum(x['size'] for x in identity.values())
        if shutil.disk_usage(root).free < required + 1024**3:
            raise ValueError(f'PGEN cache needs {required / 1024**3:.1f} GiB plus 1 GiB free at {root}')
        staging.mkdir()
        copied = {}
        try:
            for suffix, source in sources.items():
                target = staging / (prefix.name + suffix)
                h = hashlib.sha256()
                count = dirty = 0
                with source.open('rb') as src, target.open('wb') as dst:
                    for block in iter(lambda: src.read(8 * 1024 * 1024), b''):
                        dst.write(block); h.update(block)
                        count += len(block); dirty += len(block)
                        if dirty >= 256 * 1024 * 1024:
                            dst.flush(); os.fdatasync(dst.fileno())
                            if hasattr(os, 'posix_fadvise'):
                                os.posix_fadvise(dst.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
                            dirty = 0
                        if count % (1024**3) < len(block):
                            print(f'[PGEN cache] {suffix}: {count / 1024**3:.1f} GiB copied', flush=True)
                    dst.flush(); os.fsync(dst.fileno())
                if fingerprint(source) != identity[suffix]:
                    raise ValueError(f'Source changed while copying: {source}')
                if _sha(target) != h.hexdigest():
                    raise ValueError(f'PGEN execution copy failed checksum: {target}')
                copied[suffix] = {'sha256': h.hexdigest()}
            os.replace(staging, destination)
            for suffix in copied:
                copied[suffix]['fingerprint'] = fingerprint(destination / (prefix.name + suffix))
            record = dict(schema=1, source_prefix=str(prefix), cached_prefix=str(destination / prefix.name),
                          source_fingerprints=identity, files=copied,
                          verification='Full SHA-256 of streamed source and independent reread of copied file; original inputs retained')
            dump_json(manifest, record)
            return str(destination / prefix.name), record
        except BaseException:
            if staging.exists():
                shutil.rmtree(staging)
            raise
