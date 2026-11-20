"""LE8 output locations, frozen cohort paths and result publication."""
from __future__ import annotations

# Shared output locations
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).resolve().parents[2] / "0f"))
from workflow_layout import named_directory, run_name

# Frozen cohort and report locations


def cohort_directory(root):
    return Path(root) / "logs/freeze"


def resolve_cohort_path(path):
    """Use the requested current path; no old-output relocation is performed."""
    return Path(path).expanduser()


def cohort_lock_path(path):
    path = resolve_cohort_path(path)
    if path.parent.name == "freeze" and path.parent.parent.name == "logs":
        return path.parent.parent / "locks" / (path.name + ".lock")
    return Path(str(path) + ".lock")


def report_directory(root):
    return Path(root) / "reports/c1"


def summary_lock_directory(root):
    return Path(root) / "logs/locks"

# Result publication
import errno
import json
import os
import shutil
import tempfile

from file_lock import file_lock

PUBLIC_SUFFIXES = {'.png', '.xlsx', '.pdf', '.svg'}
PRIVATE_PARTS = {'tmp', 'inputs', 'donor_evidence', 'checkpoints', 'qc', '_history'}


def flat_result_names(scope):
    path = Path(scope) / 'qc/published_results.json'
    return set(json.loads(path.read_text()).get('files', {})) if path.is_file() else set()


def _sources(scope):
    for label in ('evaluation', 'reliability'):
        directory = scope / 'qc' / label
        if not (directory / 'COMPLETE.json').is_file():
            continue
        books = sorted(directory.glob('*.xlsx'))
        # Packing supersedes the original workbook; publish one current book.
        packed = [p for p in books if p.name.endswith('.results.xlsx')]
        candidates = packed or books
        preferred = max(candidates, key=lambda p: (p.stat().st_mtime_ns, p.name)) if candidates else None
        if preferred:
            yield preferred, label + '.xlsx'
        for path in sorted(directory.glob('*')):
            if path.suffix.lower() in PUBLIC_SUFFIXES - {'.xlsx'}:
                yield path, label + '.' + path.name
    models = scope / 'tmp/models'
    index = scope / 'qc/C1_FINAL_INDEX.json'
    if not index.is_file():
        index = scope / 'qc/C1_MODELS_FROZEN.json'
    recorded = json.loads(index.read_text()).get('methods', {}) if index.is_file() else None
    selected = []
    if recorded is not None:
        for method, entry in recorded.items():
            if not isinstance(entry, dict) or not entry.get('path'):
                continue
            version = Path(entry['path'])
            if not version.is_absolute():
                version = scope / version
            if version.is_relative_to(models) and (version / 'COMPLETE.json').is_file():
                selected.append((method, version))
    else:
        for method in sorted(models.glob('*')):
            versions = [v for v in method.glob('*') if (v / 'COMPLETE.json').is_file()]
            if versions:
                selected.append((method.name, max(versions, key=lambda p: p.stat().st_mtime_ns)))
    for method, version in selected:
        for path in sorted(version.rglob('*')):
            relative = path.relative_to(version)
            if (path.is_file() and path.suffix.lower() in PUBLIC_SUFFIXES
                    and not PRIVATE_PARTS.intersection(relative.parts[:-1])
                    and 'LOCAL_ONLY' not in str(relative)
                    and not any(part.startswith(('.', '_')) for part in relative.parts)):
                yield path, '.'.join((method, *relative.parts))
    for path in sorted((scope / 'performance').glob('*')):
        if path.is_file() and path.suffix.lower() in PUBLIC_SUFFIXES:
            yield path, path.name


def staging_directory(scope):
    path = Path(scope) / 'tmp/staging'
    path.mkdir(parents=True, exist_ok=True)
    return path


def publish_flat_results(scope):
    scope = Path(scope)
    qc = scope / 'qc'
    manifest = qc / 'published_results.json'
    sources = list(_sources(scope))
    if not sources and not manifest.is_file():
        return 0
    qc.mkdir(parents=True, exist_ok=True)
    with file_lock(qc / 'publish_results.lock'):
        old = json.loads(manifest.read_text()).get('files', {}) if manifest.is_file() else {}
        files = {}
        for source, name in sources:
            if Path(name).name != name or name in files:
                raise ValueError('Ambiguous flat result name: ' + name)
            target = scope / name
            if not target.exists() or not os.path.samefile(source, target):
                fd, staging = tempfile.mkstemp(prefix='.publish-', dir=qc)
                os.close(fd)
                os.unlink(staging)
                try:
                    try:
                        os.link(source, staging)
                    except OSError as error:
                        if error.errno != errno.EXDEV:
                            raise
                        shutil.copy2(source, staging)
                    os.replace(staging, target)
                finally:
                    if os.path.exists(staging):
                        os.unlink(staging)
            stat = target.stat()
            files[name] = dict(source=str(source.relative_to(scope)), device=stat.st_dev,
                               inode=stat.st_ino, size=stat.st_size, mtime_ns=stat.st_mtime_ns)
        # Remove only obsolete publisher-owned links, never a user's replacement.
        for name, saved in old.items():
            if name in files:
                continue
            if Path(name).name != name:
                raise ValueError('Unsafe published result manifest')
            target = scope / name
            if target.is_file():
                stat = target.stat()
                if (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns) == tuple(
                        saved.get(k) for k in ('device', 'inode', 'size', 'mtime_ns')):
                    target.unlink()
        with tempfile.NamedTemporaryFile(mode='w', dir=qc, delete=False) as stream:
            json.dump(dict(files=files), stream, indent=2, sort_keys=True)
            temporary = stream.name
        os.replace(temporary, manifest)
    return len(files)
