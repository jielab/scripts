"""Bounded native-filesystem cache for read-only LD HDF5 references on /mnt drives."""
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile


def fingerprint(path):
    s = path.stat()
    return {'source': str(path.resolve()), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}


@contextmanager
def cached_ld_file(source):
    source = Path(source)
    directory = os.getenv('GRID_CSX_LD_CACHE', '')
    if not directory:
        yield source
        return
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    limit = int(float(os.getenv('GRID_CSX_LD_CACHE_GIB', '16')) * 2**30)
    if limit < 0:
        raise ValueError('GRID_CSX_LD_CACHE_GIB must be nonnegative')
    identity = fingerprint(source)
    if identity['bytes'] > limit:
        yield source
        return
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    dest = root / f'ld-{key}.hdf5'
    metadata = dest.with_suffix('.json')
    # Pin the file through HDF5 close. Eviction never removes another reader's
    # reference, and two workers never race to publish the same cache entry.
    with dest.with_suffix('.lock').open('a') as pin:
        fcntl.flock(pin, fcntl.LOCK_EX)
        valid = False
        if dest.is_file() and metadata.is_file():
            try:
                valid = json.loads(metadata.read_text()) == identity and dest.stat().st_size == identity['bytes']
            except (OSError, ValueError):
                pass
        if not valid:
            with (root / '.cache.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                entries = sorted((p for p in root.glob('ld-*.hdf5') if p.with_suffix('.json').is_file()),
                                 key=lambda p: p.stat().st_mtime_ns)
                used = sum(p.stat().st_size for p in entries if p != dest)
                for old in entries:
                    if used + identity['bytes'] <= limit:
                        break
                    if old == dest:
                        continue
                    with old.with_suffix('.lock').open('a') as old_lock:
                        try:
                            fcntl.flock(old_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        except BlockingIOError:
                            continue
                        size = old.stat().st_size
                        old.unlink()
                        old.with_suffix('.json').unlink(missing_ok=True)
                        used -= size
                if used + identity['bytes'] > limit or shutil.disk_usage(root).free < identity['bytes'] + 512*2**20:
                    fcntl.flock(lock, fcntl.LOCK_UN)
                    yield source
                    return
                fd, temporary = tempfile.mkstemp(prefix='.ld-', suffix='.part', dir=root)
                try:
                    print('... LD cache COPY: %s ...' % source, flush=True)
                    with os.fdopen(fd, 'wb') as output, source.open('rb') as handle:
                        shutil.copyfileobj(handle, output, length=4*2**20)
                    if fingerprint(source) != identity or Path(temporary).stat().st_size != identity['bytes']:
                        raise ValueError('LD reference changed during cache copy')
                    os.replace(temporary, dest)
                    # Metadata is published last, under the same per-file lock.
                    metadata.write_text(json.dumps(identity, sort_keys=True)+'\n')
                finally:
                    Path(temporary).unlink(missing_ok=True)
        else:
            print('... LD cache HIT: %s ...' % source, flush=True)
        os.utime(dest, None)
        yield dest
