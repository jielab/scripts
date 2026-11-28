"""Small shared runtime utilities. No data or GPU access on import."""
from __future__ import annotations
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import random
import time

OUTPUT_FORMAT = 2
ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.append(str(ROOT.parent / "0f"))
from ygb_followup import resolve_followup, validate_record
from workflow_cpu import core_budget, worker_cores
from workflow_paths import resolved_table


def merge(left, right):
    result = copy.deepcopy(left)
    for key, value in right.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def load_config(path=None, *, overrides=None):
    import yaml
    fu = resolve_followup()
    with open(ROOT / 'config/flux.yml') as handle:
        cfg = yaml.safe_load(handle)
    validate_record(cfg, "data", "admin_end", fu)
    if path and Path(path).resolve() != (ROOT / 'config/flux.yml').resolve():
        with open(path) as handle:
            supplied = yaml.safe_load(handle) or {}
        validate_record(supplied, "data", "admin_end", fu)
        cfg = merge(cfg, supplied)
    if overrides is not None:
        if not isinstance(overrides, dict):
            raise ValueError('Configuration overrides must be a YAML mapping')
        validate_record(overrides, "data", "admin_end", fu)
        cfg = merge(cfg, overrides)
    for key, value in resolved_table().items():
        os.environ.setdefault(key, value)
    os.environ.setdefault('FLUX_GPU_QUEUE', '2v100-32-e5')

    def expand(value):
        if isinstance(value, dict): return {k: expand(v) for k, v in value.items()}
        if isinstance(value, list): return [expand(v) for v in value]
        if isinstance(value, str): return os.path.expandvars(os.path.expanduser(value))
        return value
    cfg = expand(cfg)
    workers = cfg.setdefault('train', {}).setdefault('workers', 0)
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 0:
        raise ValueError('train.workers must be a nonnegative integer')
    cfg['resources'].setdefault('cores', core_budget())
    cfg["data"]["admin_end"] = fu["date"]
    cfg["followup_provenance"] = fu
    for key in ('ygb', 'out'):
        cfg['paths'][key] = str(Path(cfg['paths'][key]).resolve())
    return cfg


def validate_resources(resources):
    for key in ('cores', 'gpus'):
        try:
            value = float(resources[key])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f'resources.{key} must be a positive finite number') from exc
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'resources.{key} must be a positive finite number')
        if key in ('cores', 'gpus') and not value.is_integer():
            raise ValueError(f'resources.{key} must be a positive integer')
        resources[key] = int(value) if key in ('cores', 'gpus') else value


def save_config(cfg, path):
    import yaml
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile('w', dir=path.parent, delete=False) as stream:
        yaml.safe_dump(cfg,stream,sort_keys=False)
    Path(stream.name).replace(path)


def log(message):
    print(time.strftime('[%Y-%m-%d %H:%M:%S]'), message, flush=True)


def outdir(cfg, sub=''):
    p = Path(cfg['paths']['out']) / sub
    p.mkdir(parents=True, exist_ok=True)
    return p


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


# Shared helper; keep the existing project API and checksum format.
import sys as _digest_sys
from pathlib import Path as _DigestPath
_digest_sys.path.append(str(_DigestPath(__file__).resolve().parents[2] / '0f'))
from file_digest import sha256_file as _shared_sha256_file

def file_digest(path):
    return _shared_sha256_file(path)


def seed_all(seed):
    import numpy as np
    import torch
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(seed)


def choose_device(requested):
    import torch
    if requested == 'gpu' and not torch.cuda.is_available():
        raise RuntimeError('GPU explicitly requested but CUDA is unavailable in this allocation. Check the wheel, driver and LSF GPU request.')
    return torch.device('cuda' if requested != 'cpu' and torch.cuda.is_available() else 'cpu')


def to_device(value, device):
    import torch
    if isinstance(value, torch.Tensor): return value.to(device, non_blocking=True)
    if isinstance(value, dict): return {k: to_device(v, device) for k, v in value.items()}
    return value


def torch_load(path, device='cpu'):
    # Only load checkpoints created in this trusted local workflow.
    import torch
    return torch.load(path, map_location=device, weights_only=False)


def atomic_checkpoint(value, path):
    import torch
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    torch.save(value, tmp)
    tmp.replace(path)
