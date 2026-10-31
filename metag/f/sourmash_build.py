#!/usr/bin/env python3
"""Build exact FracMinHash sketches with bounded FASTA reads and checkpoints."""
from __future__ import annotations
import csv
import fcntl
import io
import json
import multiprocessing
import shutil
import time
import zipfile
from pathlib import Path

from common import DIR_LIB, SOURMASH_K, SOURMASH_SCALED, SOURMASH_BUILD_WORKERS, atomic_json
from library import active_references, file_state, digest
from fasta_stream import sketch

VERSION = 1


def read_manifest(path):
    with zipfile.ZipFile(path) as archive:
        with archive.open('SOURMASH-MANIFEST.csv') as raw:
            lines = io.TextIOWrapper(raw, encoding='utf-8')
            rows = list(csv.DictReader(line for line in lines if not line.startswith('#')))
        members = set(archive.namelist())
        if any(row['internal_location'] not in members for row in rows):
            raise ValueError(f'Missing signatures in {path}')
        return rows


def validate(path, expected, ksize, scaled, crc=False):
    rows = read_manifest(path)
    names = [row['name'] for row in rows]
    if len(names) != len(expected) or set(names) != set(expected):
        raise ValueError(f'Incomplete/duplicate reference signatures in {path}')
    if any(int(r['ksize']) != ksize or int(r['scaled']) != scaled or r['moltype'] != 'DNA'
           or r['with_abundance'].lower() not in ('0', 'false') for r in rows):
        raise ValueError(f'Unexpected sketch parameters in {path}')
    if crc:
        with zipfile.ZipFile(path) as archive:
            bad = archive.testzip()
            if bad:
                raise ValueError(f'Corrupt signature member {bad} in {path}')
    return rows


def batches(rows, folder, params, max_bytes=256 << 20, max_files=250):
    group, size = [], 0
    for row in rows:
        state = file_state(row['path'])
        if state[0] != int(row['bytes']):
            raise RuntimeError(f"Reference changed since check-ref: {row['path']}")
        if group and (size + state[0] > max_bytes or len(group) >= max_files):
            yield batch_task(group, folder, params)
            group, size = [], 0
        group.append({'path': row['path'], 'name': row['filename'], 'state': state})
        size += state[0]
    if group:
        yield batch_task(group, folder, params)


def batch_task(group, folder, params):
    key = digest({'inputs': group, 'params': params, 'version': VERSION})
    return {'inputs': group, 'params': params, 'key': key, 'folder': str(folder)}


def cached_batch(task):
    path = Path(task['folder'])/(task['key'] + '.zip')
    try:
        receipt = json.loads(path.with_suffix('.json').read_text())
        if receipt.get('key') != task['key'] or receipt.get('output') != file_state(path):
            return False
        validate(path, [r['name'] for r in task['inputs']], **task['params'])
        return True
    except (OSError, ValueError, KeyError, zipfile.BadZipFile):
        return False


def build_batch(task):
    from sourmash import SourmashSignature
    from sourmash.sourmash_args import SaveSignaturesToLocation
    path = Path(task['folder'])/(task['key'] + '.zip')
    temp = path.with_name(task['key'] + '.partial.zip')
    temp.unlink(missing_ok=True)  # Only this derived batch's unfinished output.
    stats = []
    current = '(opening output)'
    try:
        with SaveSignaturesToLocation(str(temp)) as output:
            for item in task['inputs']:
                current = item['path']
                if file_state(current) != item['state']:
                    raise RuntimeError(f'Reference changed: {current}')
                mh, records, bases = sketch(current, **task['params'])
                output.add(SourmashSignature(mh, name=item['name'], filename=current))
                stats.append({'name': item['name'], 'records': records, 'bases': bases})
                if file_state(current) != item['state']:
                    raise RuntimeError(f'Reference changed during sketching: {current}')
        validate(temp, [r['name'] for r in task['inputs']], **task['params'], crc=True)
        temp.replace(path)
        atomic_json(path.with_suffix('.json'), {'key': task['key'], 'output': file_state(path), 'references': stats})
    except Exception as exc:
        raise RuntimeError(f'Sourmash batch failed at {current}: {exc}') from exc
    return task['key'], len(task['inputs'])


def merge_batches(tasks, target, params):
    temporary = target.with_name('references.partial.zip')
    manifest = []
    with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_STORED, allowZip64=True) as dest:
        for task in tasks:
            source = Path(task['folder'])/(task['key'] + '.zip')
            rows = validate(source, [r['name'] for r in task['inputs']], **params)
            with zipfile.ZipFile(source) as archive:
                copied = set()
                for row in rows:
                    original = row['internal_location']
                    location = f"signatures/{task['key']}/{Path(original).name}"
                    if original not in copied:
                        with archive.open(original) as src, dest.open(location, 'w', force_zip64=True) as out:
                            shutil.copyfileobj(src, out, length=4 << 20)
                        copied.add(original)
                    manifest.append({**row, 'internal_location': location})
        data = io.StringIO()
        data.write('# SOURMASH-MANIFEST-VERSION: 1.0\n')
        writer = csv.DictWriter(data, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
        dest.writestr('SOURMASH-MANIFEST.csv', data.getvalue())
    validate(temporary, [r['name'] for t in tasks for r in t['inputs']], **params)
    temporary.replace(target)


def main(workers=None):
    out = DIR_LIB/'sourmash'
    out.mkdir(parents=True, exist_ok=True)
    with (out/'build.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Another sourmash build is already running') from exc
        return build(out, workers)


def run_pending(tasks, workers, heartbeat=None):
    if not tasks:
        return
    with multiprocessing.Pool(processes=workers) as pool:
        # Keep original worker handles: Pool replaces an OOM-killed worker but
        # would otherwise wait forever for the lost task's result (Python 3.11).
        processes = tuple(pool._pool)
        results = pool.imap_unordered(build_batch, tasks, chunksize=1)
        remaining = len(tasks)
        while remaining:
            try:
                result = results.next(timeout=1)
            except multiprocessing.TimeoutError:
                if any(p.exitcode not in (None, 0) for p in processes):
                    raise RuntimeError("Sourmash worker exited unexpectedly, possibly due to memory pressure; "
                                       "completed batches are retained. Reduce SOURMASH_BUILD_WORKERS and rerun.")
                if heartbeat:
                    heartbeat()
                continue
            remaining -= 1
            yield result


def progress_line(done, total, completed_bytes, total_bytes, elapsed, new_bytes):
    rate = new_bytes / elapsed if elapsed > 0 else 0
    eta = f'; rough ETA {(total_bytes-completed_bytes)/rate/3600:.1f} h' if rate else ''
    return (f'[sourmash] completed {done:,}/{total:,} references; '
            f'completed input {completed_bytes/2**30:.1f}/{total_bytes/2**30:.1f} GiB '
            f'({completed_bytes/total_bytes if total_bytes else 0:.1%} of compressed bytes); '
            f'this run {elapsed/60:.1f} min; {rate/2**20:.1f} MiB/s{eta}')


def build(out, workers=None):
    workers = SOURMASH_BUILD_WORKERS if workers is None else workers
    if not 1 <= workers <= 32:
        raise ValueError('Sourmash workers must be between 1 and 32')
    rows = active_references()
    params = {'ksize': SOURMASH_K, 'scaled': SOURMASH_SCALED}
    cache = out/'batches'
    cache.mkdir(exist_ok=True)
    tasks = list(batches(rows, cache, params))
    identity = digest([t['key'] for t in tasks])
    target, receipt = out/'references.zip', out/'build.json'
    try:
        saved = json.loads(receipt.read_text())
        if (saved.get('status') == 'complete' and saved.get('identity') == identity
                and saved.get('output') == file_state(target)):
            validate(target, [r['filename'] for r in rows], **params)
            print(f'REUSE completed sourmash library: {target}', flush=True)
            return 0
    except (OSError, ValueError, KeyError, zipfile.BadZipFile):
        pass
    atomic_json(receipt, {'status': 'building', 'identity': identity, 'references': len(rows), 'params': params})
    pending = [task for task in tasks if not cached_batch(task)]
    done = len(rows) - sum(len(t['inputs']) for t in pending)
    sizes = {t['key']:sum(r['state'][0] for r in t['inputs']) for t in tasks}
    total_bytes = sum(sizes.values())
    completed_bytes = total_bytes - sum(sizes[t['key']] for t in pending)
    reused_bytes = completed_bytes
    # Only scheduling changes; checkpoint keys and final library identity stay stable.
    pending.sort(key=lambda task:sizes[task['key']], reverse=True)
    print(f'sourmash: {len(rows):,} references; {workers} CPU worker processes; 4 MiB FASTA blocks; '
          f'{done:,} references reused from completed batches\nDatabase: {target}', flush=True)
    print('Progress/ETA use completed compressed input bytes; genome expansion and speed vary. '
          'This sourmash stage uses CPU; the following MetaCache stage can use GPU.', flush=True)
    started = time.monotonic()
    next_report = (done // 1000 + 1) * 1000

    def report(force=False):
        nonlocal next_report
        if not force and done < next_report:
            return
        now=time.monotonic()
        print(progress_line(done,len(rows),completed_bytes,total_bytes,now-started,completed_bytes-reused_bytes),flush=True)
        next_report = (done // 1000 + 1) * 1000

    report(True)
    try:
        # Terminate workers on failure/Ctrl-C; only completed atomic batches are reused.
        for key, count in run_pending(pending, workers, heartbeat=report):
            done += count
            completed_bytes += sizes[key]
            report(done == len(rows))
        if any(file_state(r['path']) != r['state'] for t in tasks for r in t['inputs']):
            raise RuntimeError('References changed during building; rerun check-ref')
        print('[sourmash] All batches complete; assemble and validate permanent library', flush=True)
        merge_batches(tasks, target, params)
        atomic_json(receipt, {'status': 'complete', 'identity': identity, 'references': len(rows),
                              'params': params, 'output': file_state(target), 'batches': len(tasks)})
    except BaseException as exc:
        atomic_json(receipt, {'status': 'incomplete', 'identity': identity, 'error': str(exc) or type(exc).__name__,
                              'references': len(rows), 'params': params})
        raise
    print(f'PASS: complete sourmash library: {target}', flush=True)
    return 0


if __name__ == '__main__':
    main()
