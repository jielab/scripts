"""Shared immutable reference/index caches and small per-panel checkpoints.

Cache only exact inputs. Never merge separately built indexes: competitive
mapping and within-panel MAPQ retain the original reference order and scope.
"""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from common import atomic_json, fingerprint, locked, sha256
from plant_validate import build_panel, file_state

SCHEMA = 2
GIB = 2**30


def read_json(path):
    try:
        value = json.loads(Path(path).read_text())
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def observation(path):
    p = Path(path)
    s = p.stat()
    # st_dev is a node-local mount number, even on the same shared filesystem.
    return dict(size=s.st_size, mtime_ns=s.st_mtime_ns, ctime_ns=s.st_ctime_ns,
                inode=s.st_ino)


def valid_entry(folder, key):
    rec = read_json(folder / 'complete.json')
    if rec.get('key') != key or rec.get('schema') != SCHEMA or not rec.get('files'):
        return None
    try:
        if any(Path(name).name != name or observation(folder / name) != state
               for name, state in rec['files'].items()):
            return None
    except OSError:
        return None
    return rec


def publish(folder, key, names, **extra):
    atomic_json(folder / 'complete.json', dict(schema=SCHEMA, key=key,
        files={name: observation(folder / name) for name in names}, **extra))


def cache_root(cfg):
    root = Path(cfg['work_dir']).absolute() / 'cache/alignment'
    root.mkdir(parents=True, exist_ok=True)
    if root.resolve() != root:
        raise RuntimeError('Alignment cache must not use a symlink: ' + str(root))
    return root


def cache_slot(root, kind, key, label):
    """Keep fingerprints in metadata and use human-readable paths on disk."""
    pool = root / kind
    registry = pool / 'entries.json'
    with locked(pool / 'entries.lock'):
        slots = json.loads(registry.read_text()) if registry.exists() else {}
        if key not in slots:
            slots[key] = f'{label}-{len(slots) + 1}'
            atomic_json(registry, slots)
        return pool / slots[key]


@contextmanager
def admission(root, kind, key, needed, cap_gib):
    """Reserve disk budget across nodes; interrupted reservations are conservative.

Same-key retries reclaim their abandoned reservation while holding the key lock.
Capacity exhaustion falls back to the original uncached operation, never changes
reference selection. No completed cache is evicted while another process reads it.
"""
    pool = root / kind
    pool.mkdir(exist_ok=True)
    ticket = pool / (key + '.reservation.json')
    allowed = False
    with locked(root / 'budget.lock'):
        ticket.unlink(missing_ok=True)
        used = 0
        for p in pool.glob('*/complete.json'):
            used += sum(int(s['size']) for s in read_json(p).get('files', {}).values())
        used += sum(int(read_json(p).get('bytes', 0)) for p in pool.glob('*.reservation.json'))
        if used + needed <= int(float(cap_gib) * GIB) and shutil.disk_usage(root).free >= needed + 2 * GIB:
            atomic_json(ticket, dict(bytes=needed, created_at=time.time()))
            allowed = True
    try:
        yield allowed
    finally:
        if allowed:
            with locked(root / 'budget.lock'):
                ticket.unlink(missing_ok=True)


@contextmanager
def normalized_reference(cfg, row, max_bases):
    """Keep the key lock until the caller has finished reading both files."""
    root = cache_root(cfg)
    source = file_state(row['resolved_path'])
    key = fingerprint(dict(schema=SCHEMA, source=source, species_taxid=int(row['species_taxid']),
                           implementation=sha256(Path(__file__).with_name('plant_validate.py'))))
    folder = cache_slot(root, 'references', key, 'taxon-' + str(int(row['species_taxid'])))
    with locked(root / 'locks' / (folder.name + '.lock')):
        rec = valid_entry(folder, key)
        if rec:
            yield folder, rec, True
            return
        # Existing incomplete files are private to this key and never published.
        if folder.exists():
            shutil.rmtree(folder)
        expected = max(int(row.get('bases', 0)), 1)
        needed = int(expected * 1.25) + 16 * 2**20
        cap = (cfg.get('stage_alignment') or {}).get('reference_cache_gib', 512)
        with admission(root, 'references', folder.name, needed, cap) as allowed:
            if not allowed:
                yield None, None, False
                return
            folder.mkdir(parents=True)
            try:
                contigs, actual = build_panel([row], folder / 'reference.fa', max_bases)
                atomic_json(folder / 'contigs.json', contigs)
                if file_state(source['path']) != source:
                    raise RuntimeError('Reference changed during cache preparation')
                if sum(p.stat().st_size for p in folder.iterdir()) > needed:
                    # Bad declared sizes must not exceed the reserved cache budget.
                    shutil.rmtree(folder)
                    yield None, None, False
                    return
                publish(folder, key, ['reference.fa', 'contigs.json'], bases=actual, source=source)
            except BaseException:
                shutil.rmtree(folder, ignore_errors=True)
                raise
            yield folder, valid_entry(folder, key), False


def build_cached_panel(cfg, selected, path, max_bases):
    if not (cfg.get('stage_alignment') or {}).get('cache_enabled', True):
        contigs, total = build_panel(selected, path, max_bases)
        return contigs, total, dict(reference_hits=0, reference_builds=0, reference_fallbacks=len(selected))
    contigs, total = {}, 0
    metrics = dict(reference_hits=0, reference_builds=0, reference_fallbacks=0)
    with Path(path).open('wb') as out:
        for row in selected:
            with normalized_reference(cfg, row, max_bases) as (folder, rec, hit):
                if folder is None:
                    temporary = Path(str(path) + '.reference.tmp')
                    try:
                        part, count = build_panel([row], temporary, max_bases - total)
                        with temporary.open('rb') as src:
                            shutil.copyfileobj(src, out, 8 << 20)
                        metrics['reference_fallbacks'] += 1
                    finally:
                        temporary.unlink(missing_ok=True)
                else:
                    count = rec['bases']
                    part = read_json(folder / 'contigs.json')
                    if total + count > max_bases:
                        raise RuntimeError('Actual FASTA exceeds declared panel budget')
                    with (folder / 'reference.fa').open('rb') as src:
                        shutil.copyfileobj(src, out, 8 << 20)
                    metrics['reference_hits' if hit else 'reference_builds'] += 1
            total += count
            if total > max_bases:
                raise RuntimeError('Actual FASTA exceeds declared panel budget')
            contigs.update(part)
    return contigs, total, metrics


@contextmanager
def panel_index(cfg, selected, executable, version, threads, budget, reference, *, progress=None):
    """Lock identical panels through mapping so repairs cannot race active readers."""
    options = cfg.get('stage_alignment') or {}
    if not options.get('cache_enabled', True):
        if progress:
            progress.stage('panel_fasta', cache='disabled')
        contigs, actual, metrics = build_cached_panel(cfg, selected, reference, budget)
        yield Path(reference), contigs, actual, dict(metrics, index_hit=False, index_cached=False)
        return
    root = cache_root(cfg)
    sources = [dict(state=file_state(row['resolved_path']), species_taxid=int(row['species_taxid'])) for row in selected]
    context = dict(schema=SCHEMA, sources=sources, mapper=file_state(executable), version=version,
                   preset='sr', index_batch='actual_bases+1000000',
                   implementation=sha256(Path(__file__).with_name('plant_validate.py')))
    key = fingerprint(context)
    folder = cache_slot(root, 'indexes', key, 'panel')
    with locked(root / 'locks' / (folder.name + '.lock')):
        rec = valid_entry(folder, key)
        if rec:
            if rec['bases'] > budget:
                raise RuntimeError('Cached reference exceeds panel budget')
            yield folder / 'reference.mmi', read_json(folder / 'contigs.json'), rec['bases'], dict(
                index_hit=True, index_cached=True, index_key=key,
                reference_hits=0, reference_builds=0, reference_fallbacks=0)
            return
        if folder.exists():
            shutil.rmtree(folder)
        if progress:
            progress.stage('panel_fasta', cache='miss')
        contigs, actual, metrics = build_cached_panel(cfg, selected, reference, budget)
        # Conservative index reservation. If full, map the original FASTA normally.
        needed = actual * 8 + 64 * 2**20
        with admission(root, 'indexes', folder.name, needed, options.get('index_cache_gib', 512)) as allowed:
            if not allowed:
                yield Path(reference), contigs, actual, dict(metrics, index_hit=False, index_cached=False,
                    index_fallback_reason='disk_or_cache_budget', index_required_bytes=needed)
                return
            folder.mkdir(parents=True)
            try:
                if progress:
                    progress.stage('build_mmi', bases=actual, threads=threads)
                command = [executable, '-x', 'sr', '-v', '2', '-I', str(actual + 1000000), '-t', str(threads),
                           '-d', str(folder / 'reference.mmi'), str(reference)]
                with (folder / 'build.log').open('w') as log:
                    subprocess.run(command, stdout=log, stderr=log, check=True)
                atomic_json(folder / 'contigs.json', contigs)
                if any(file_state(s['state']['path']) != s['state'] for s in sources) or file_state(executable) != context['mapper']:
                    raise RuntimeError('Reference or mapper changed during index preparation')
                if sum(p.stat().st_size for p in folder.iterdir()) > needed:
                    shutil.rmtree(folder)
                    yield Path(reference), contigs, actual, dict(metrics, index_hit=False, index_cached=False, index_fallback_reason='index_exceeded_reservation')
                    return
                publish(folder, key, ['reference.mmi', 'contigs.json', 'build.log'], bases=actual)
            except BaseException:
                shutil.rmtree(folder, ignore_errors=True)
                raise
            yield folder / 'reference.mmi', contigs, actual, dict(metrics, index_hit=False, index_cached=True, index_key=key)


def checkpoint_read(path, key):
    rec = read_json(path)
    if rec.get('key') == key and isinstance(rec.get('evidence'), list) and rec.get('digest') == fingerprint(rec['evidence']):
        return rec
    return None


def checkpoint_write(path, key, evidence, command):
    atomic_json(path, dict(key=key, evidence=evidence, digest=fingerprint(evidence), command=command))
