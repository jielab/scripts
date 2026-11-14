"""Create immutable source snapshots for scheduled workers and their shared helpers."""
from pathlib import Path
import os
import shutil
import tempfile

from common import CODE_ROOT, SHARED_RUNTIME, atomic_json, code_identity, locked, sha256


def freeze_source(cfg):
    identity = code_identity()
    root = Path(cfg['work_dir']) / 'releases' / identity
    root.parent.mkdir(parents=True, exist_ok=True)
    with locked(root.parent / 'freeze.lock'):
        if root.exists():
            import json
            receipt = json.loads((root / 'release.json').read_text())
            if receipt['code_identity'] != identity:
                raise RuntimeError('Frozen source identity mismatch')
            for relative, digest in receipt['files'].items():
                if sha256(root / relative) != digest:
                    raise RuntimeError('Frozen source changed: ' + relative)
            return root / 'metag'
        stage = Path(tempfile.mkdtemp(prefix='.source-', dir=root.parent))
        try:
            ignore = shutil.ignore_patterns('__pycache__', '*.pyc', '.git', '.pytest_cache')
            shutil.copytree(CODE_ROOT, stage / 'metag', ignore=ignore)
            shutil.copytree(SHARED_RUNTIME, stage / '0f', ignore=ignore)
            if code_identity() != identity:
                raise RuntimeError('Source changed while freezing a worker release')
            files = {str(p.relative_to(stage)): sha256(p) for p in sorted(stage.rglob('*')) if p.is_file()}
            atomic_json(stage / 'release.json', dict(schema='metag-source-release-v1', code_identity=identity, files=files))
            for path in stage.rglob('*'):
                path.chmod(0o555 if path.is_dir() or os.access(path,os.X_OK) else 0o444)
            stage.chmod(0o555)
            stage.rename(root)
        except BaseException:
            # Retain an incomplete snapshot for diagnosis; never reuse it.
            raise
    return root / 'metag'
