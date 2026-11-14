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

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.append(str(ROOT.parent / "0f"))
from ukb_followup import resolve_followup, validate_record


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
    os.environ.setdefault('DIR_UKB', '/mnt/d/phe/ukb')
    os.environ.setdefault('ANALYSIS_ROOT', '/mnt/d/analysis')
    os.environ.setdefault('HPC_GPU_QUEUE', '4a100-40')
    os.environ.setdefault('FLUX_GPU_QUEUE', '2v100-32-e5')

    def expand(value):
        if isinstance(value, dict): return {k: expand(v) for k, v in value.items()}
        if isinstance(value, list): return [expand(v) for v in value]
        if isinstance(value, str): return os.path.expandvars(os.path.expanduser(value))
        return value
    cfg = expand(cfg)
    cfg["data"]["admin_end"] = fu["date"]
    cfg["followup_provenance"] = fu
    for key in ('ukb', 'out'):
        cfg['paths'][key] = str(Path(cfg['paths'][key]).resolve())
    return cfg


def validate_resources(resources):
    for key in ('cores', 'gpus', 'cpu_memory_gib', 'gpu_memory_gib', 'hours'):
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
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))


def log(message):
    print(time.strftime('[%Y-%m-%d %H:%M:%S]'), message, flush=True)


def outdir(cfg, sub=''):
    p = Path(cfg['paths']['out']) / sub
    p.mkdir(parents=True, exist_ok=True)
    return p


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def file_digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


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
