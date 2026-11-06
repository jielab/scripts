"""Fail-closed scientific checks, cached by source/config and runtime identity."""

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]


def fingerprint(root=ROOT):
    h=hashlib.sha256()
    files=sorted(p for folder in ('f','config','2test') for p in (root/folder).rglob('*')
                 if p.is_file() and p.suffix in ('.py','.sh','.json','.R'))
    files += [root/'grid.sh',root/'4.prsformer.sh',root/'environment.yml']
    for p in files:
        h.update(str(p.relative_to(root)).encode());h.update(p.read_bytes())
    for name in ('abm.py', 'phenotype.R'):
        shared=root.parent/'0f'/name
        if shared.is_file():h.update(shared.read_bytes())
    h.update(sys.executable.encode());h.update(sys.version.encode())
    from importlib.metadata import version, PackageNotFoundError
    for name in ('numpy','pandas','scipy','scikit-learn','torch','natten'):
        try: value=version(name)
        except PackageNotFoundError: value='missing'
        h.update((name+value).encode())
    return h.hexdigest()


def ensure_self_check(root=ROOT,cache=None):
    digest=fingerprint(root)
    cache=Path(cache or Path(os.environ.get('GRID_CACHE_ROOT',configured_path('ANALYSIS_ROOT', 'grid/cache')))/'self_check')
    cache.mkdir(parents=True,exist_ok=True,mode=0o700)
    receipt=cache/'self-check.json'
    if receipt.is_file():
        try: saved=json.loads(receipt.read_text())
        except (ValueError,OSError): saved={}
        if saved.get('status')=='passed' and saved.get('fingerprint')==digest:
            print('[GRID] REUSE scientific self-check '+digest[:12],flush=True);return
    print('[GRID] Running scientific self-check '+digest[:12],flush=True)
    subprocess.run([sys.executable,str(root/'f/self_check.py')],check=True)
    if fingerprint(root)!=digest:raise RuntimeError('Source changed during self-check; rerun')
    fd,tmp=tempfile.mkstemp(dir=cache,suffix='.tmp')
    with os.fdopen(fd,'w') as f:json.dump(dict(status='passed',fingerprint=digest),f)
    os.replace(tmp,receipt)
