"""Validated stage receipts for local GRID runs; reports never invalidate a fit."""
from contextlib import contextmanager
from pathlib import Path
import fcntl
import hashlib
import json
import os
import uuid
import shutil
import joblib


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp.' + uuid.uuid4().hex)
    try:
        tmp.write_text(json.dumps(value, sort_keys=True, indent=2))
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


@contextmanager
def run_lock(root):
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    with (root / '.resume.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield


class Checkpoints:
    def __init__(self, root, contract):
        self.root = Path(root)
        self.path = self.root / 'RESUME.json'
        self.data = json.loads(self.path.read_text()) if self.path.is_file() else dict(contract=contract, stages={})
        if self.data['contract'] != contract:
            raise ValueError('Run inputs/configuration differ; use --new-run or a new --run-id: ' + str(root))
        self.root.mkdir(parents=True, exist_ok=True)
        atomic_json(self.path, self.data)

    def valid(self, name):
        files = self.data['stages'].get(name, {})
        return bool(files) and all((self.root / p).is_file() and sha(self.root / p) == digest
                                   for p, digest in files.items())

    def mark(self, name, files):
        files = list(files)
        if not files:
            raise ValueError('Cannot mark an empty stage: ' + name)
        self.data['stages'][name] = {str(Path(p).relative_to(self.root)): sha(p) for p in files}
        atomic_json(self.path, self.data)

    def adopt(self, name, path):
        target = self.root / '_LOCAL_ONLY' / ('resume_' + name + '.joblib')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        self.mark(name, [target])

    def object(self, name, compute):
        path = self.root / '_LOCAL_ONLY' / ('resume_' + name + '.joblib')
        if self.valid(name):
            saved = next(iter(self.data['stages'][name]))
            print('[GRID] SKIP ' + name + ': verified checkpoint', flush=True)
            return joblib.load(self.root / saved)
        # A recorded but corrupt scientific checkpoint must not silently retrain.
        if name in self.data['stages']:
            raise ValueError('Checkpoint missing or changed: ' + name + '; restore it or use --new-run')
        value = compute()
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + '.tmp.' + uuid.uuid4().hex)
        try:
            joblib.dump(value, tmp, compress=3)
            os.replace(tmp, path)
        finally:
            tmp.unlink(missing_ok=True)
        self.mark(name, [path])
        return value

    def report(self, name, directory, write):
        if self.valid(name):
            print('[GRID] SKIP ' + name + ': verified outputs', flush=True)
            return
        write()
        self.mark(name, [p for p in Path(directory).rglob('*') if p.is_file()])


def choose_run(root, run_id=None, new_run=False):
    root = Path(root)
    if run_id:
        if not __import__('re').fullmatch(r'[A-Za-z0-9_-]+', run_id):
            raise ValueError('Unsafe run name')
        path = root / run_id
        if new_run and path.exists():
            raise FileExistsError(path)
        return path
    if not new_run:
        # Resume a newer interrupted attempt before returning to the published run.
        candidates = [p for p in root.iterdir() if p.is_dir() and
                      ((p / 'RESUME.json').is_file() or (p / 'protocol.json').is_file())]
        if candidates:
            return max(candidates, key=lambda p: (p/'protocol.json').stat().st_mtime_ns if (p/'protocol.json').exists() else (p/'RESUME.json').stat().st_mtime_ns)
    from datetime import datetime, timezone
    return root / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:8])


def validate_manifests(root):
    root = Path(root)
    manifests = list(root.rglob('table_manifest.json'))
    if not manifests:
        return False
    for path in manifests:
        for item in json.loads(path.read_text()):
            target = path.parent / item['path']
            if not target.is_file() or sha(target) != item['sha256']:
                return False
    return True
