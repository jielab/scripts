"""Content identities and atomic, locked publication for derived analyses.

Permanent inputs are read only. A failed refresh invalidates current evidence,
while preserving the previous files and receipt in an attempt history.
"""
from __future__ import annotations
import copy
from dataclasses import dataclass
import datetime as dt
import fcntl
import functools
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import tempfile
import uuid

TERMINAL = {
    'complete', 'complete_descriptive_imputed_panel', 'complete_descriptive_all_allele_panel',
    'complete_descriptive_unverified', 'complete_adjustment_unavailable', 'planned',
    'no_usable_exact_high_quality_contrast', 'no_usable_contrast', 'not_estimable',
    'not_available', 'not_run',
}
EVIDENCE = TERMINAL - {'planned', 'not_available', 'not_run'}
SMALL_LIMIT = 64 * 1024 * 1024


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def identity(path, full=False):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        return {'path': str(path), 'status': 'missing'}
    s = path.stat()
    result = {'path': str(path), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns}
    if full:
        result.update(sha256=sha256(path), verification='full_sha256_critical_input_only')
    else:
        result.update(verification='size_mtime_only')
    if path.stat().st_mtime_ns != s.st_mtime_ns or path.stat().st_size != s.st_size:
        raise ValueError(f'Input changed while hashing: {path}')
    return result


def atomic_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.partial')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False, default=str) + '\n')
    os.replace(temporary, path)


def code_identity():
    root = Path(__file__).resolve().parent.parent
    files = {str(p.relative_to(root)): {'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns}
             for p in sorted(root.rglob('*')) if p.is_file() and p.suffix in {'.py', '.sh', '.env', '.yml'} and '__pycache__' not in p.parts}
    return {'files': files, 'identity_policy': 'file_size_mtime; no routine source hashing',
            'python': sys.version, 'executable': sys.executable, 'git_commit': 'unavailable_local_checkout_without_git'}


def inventory(root, manifest):
    root = Path(root)
    return {str(p.relative_to(root)): {k: v for k, v in identity(p).items() if k != 'path'}
            for p in sorted(root.rglob('*')) if p.is_file() and p != root / manifest}


def validate_manifest(root, name='manifest.json', require_evidence=False):
    root = Path(root)
    result = json.loads((root / name).read_text())
    allowed = EVIDENCE if require_evidence else TERMINAL
    if result.get('status') not in allowed:
        raise ValueError(f'Nonterminal/unusable result: {root}: {result.get("status", "missing")}')
    outputs = result.get('output_inventory')
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError(f'Missing verified output inventory: {root / name}')
    for rel, expected in outputs.items():
        path = root / rel
        if Path(rel).is_absolute() or '..' in Path(rel).parts or path.is_symlink():
            raise ValueError(f'Unsafe inventory entry: {rel}')
        actual = identity(path, full='sha256' in expected)
        keys = ('size', 'sha256') if 'sha256' in expected else ('size', 'mtime_ns')
        if any(actual.get(k) != expected.get(k) for k in keys):
            raise ValueError(f'Output identity changed: {path}')
    return result


def _relocate(value, source, destination):
    if isinstance(value, dict):
        return {_relocate(k, source, destination): _relocate(v, source, destination) for k, v in value.items()}
    if isinstance(value, list): return [_relocate(v, source, destination) for v in value]
    if isinstance(value, str) and (value == source or value.startswith(source + '/')):
        return destination + value[len(source):]
    return value


@dataclass(frozen=True)
class InputPath:
    path: str | Path
    kind: str = "file"


def check_output_paths(output, inputs):
    """Call before mkdir, lock files, running manifests, staging or rotation."""
    out = Path(output).expanduser().resolve()
    if out == Path(out.anchor) or out.is_file():
        raise ValueError(f"Invalid output directory: {out}")
    for item in inputs:
        p = Path(item.path).expanduser().resolve()
        if out == p or out in p.parents or p in out.parents:
            raise ValueError(f"Output overlaps protected {item.kind}: {out} <-> {p}")
    return out


def declared_inputs(raw, *, input_dirs=(), input_files=(), input_prefixes=(), resolver=None):
    result = []
    for kind, keys in (("directory", input_dirs), ("file", input_files)):
        for key in keys:
            if raw.get(key): result.append(InputPath(raw[key], kind))
    for key in input_prefixes:
        if raw.get(key):
            prefix = Path(raw[key]).expanduser()
            result += [InputPath(str(prefix) + suffix, "file") for suffix in (".pgen", ".pvar", ".pvar.zst", ".psam")]
    # Preserve file tracking and protect directory-valued inputs before writes.
    excluded = {'out', 'out_dir', 'dir_out', 'analysis_root', 'dir_analysis_root'}
    for key, value in raw.items():
        if key in excluded or not isinstance(value, (str, Path)) or not value:
            continue
        path = Path(value).expanduser()
        if path.is_file(): result.append(InputPath(path, 'file'))
        elif path.is_dir(): result.append(InputPath(path, 'directory'))
    # Shards and annotation roots can live outside the scan root after relocation.
    modern = raw.get('modern_dir') or raw.get('modern')
    base = Path(modern).expanduser().resolve() if modern else None
    if base is not None and (base / 'run.json').is_file():
        from scan_store import load_run_metadata
        run = load_run_metadata(base)
        result.append(InputPath(base / 'run.json', 'file'))
        for entry in run.get('shards', []):
            shard = Path(entry).expanduser()
            if not shard.is_absolute(): shard = base / shard
            result.append(InputPath(shard, 'directory'))
            for filename in ('COMPLETE.json', 'scan.json', 'multiallelic_sites.tsv.gz'):
                if (shard / filename).is_file(): result.append(InputPath(shard / filename, 'file'))
        options = run.get('options', run.get('identity', {}).get('options', {}))
        for key in ('architecture_annotation_cache', 'ld_block_dir'):
            if options.get(key): result.append(InputPath(options[key], 'directory'))
        for key in ('architecture_lcr_bed', 'architecture_segdup_bed', 'architecture_low_mappability_bed'):
            if options.get(key): result.append(InputPath(options[key], 'file'))
    if resolver is not None:
        result.extend(resolver(raw))
    dedup = {}
    for item in result:
        p = Path(item.path).expanduser().resolve()
        dedup[(p, item.kind)] = InputPath(p, item.kind)
    return tuple(dedup.values())


def transactional(manifest, out_arg="out", *, input_dirs=(), input_files=(),
                  input_prefixes=("pgen_prefix",), resolve_inputs=None):
    """Publish derived outputs after read-only input/output path validation."""
    def decorate(function):
        signature = inspect.signature(function)
        @functools.wraps(function)
        def wrapped(*args, **kwargs):
            bound = signature.bind(*args, **kwargs); bound.apply_defaults()
            namespace = out_arg not in bound.arguments
            original = bound.arguments["args"] if namespace else None
            raw = dict(vars(original)) if namespace else dict(bound.arguments)
            output_value = raw[out_arg]
            protected = declared_inputs(raw, input_dirs=input_dirs, input_files=input_files,
                                        input_prefixes=input_prefixes, resolver=resolve_inputs)
            output = check_output_paths(output_value, protected)
            # No filesystem mutation has occurred above this line.
            output.parent.mkdir(parents=True, exist_ok=True)
            with (output.parent / ("." + output.name + ".lock")).open("a") as lock:
                try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as exc: raise RuntimeError(f"Output in use: {output}") from exc
                # Re-evaluate aliases/dynamic source paths while holding the lock.
                protected = declared_inputs(raw, input_dirs=input_dirs, input_files=input_files,
                                            input_prefixes=input_prefixes, resolver=resolve_inputs)
                check_output_paths(output, protected)
                input_files_now = {str(Path(p.path)): identity(p.path) for p in protected if Path(p.path).is_file()}
                run_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S") + "-" + uuid.uuid4().hex[:8]
                history = output.parent / ("." + output.name + ".history")
                history.mkdir(exist_ok=True)
                staging = Path(tempfile.mkdtemp(prefix=run_id + "-", dir=history))
                old_manifest = output / manifest
                if old_manifest.is_file():
                    (staging / "previous_manifest.receipt").write_bytes(old_manifest.read_bytes())
                request = {k: str(v) for k, v in raw.items()}
                started = {"status": "running", "run_id": run_id, "request": request,
                           "inputs": input_files_now, "protected_inputs": [dict(path=str(p.path), kind=p.kind) for p in protected],
                           "code": code_identity(), "validation_policy": "declared input paths before writes; per-run identities"}
                output.mkdir(parents=True, exist_ok=True)
                atomic_json(old_manifest, started)
                try:
                    if namespace:
                        cloned = copy.copy(original); setattr(cloned, out_arg, str(staging)); bound.arguments["args"] = cloned
                    else:
                        bound.arguments[out_arg] = staging
                    result = function(*bound.args, **bound.kwargs)
                    if not isinstance(result, dict): result = json.loads((staging / manifest).read_text())
                    if result.get("status") not in TERMINAL:
                        raise ValueError(f"Stage did not finish: {result.get('status', 'missing')}")
                    if any(identity(path) != saved for path, saved in input_files_now.items()):
                        raise ValueError("Scientific inputs changed during analysis")
                    for path in staging.rglob("*.json"):
                        try: value = json.loads(path.read_text())
                        except (ValueError, UnicodeError): continue
                        atomic_json(path, _relocate(value, str(staging), str(output)))
                    result = _relocate(result, str(staging), str(output))
                    result.update(run_id=run_id, request=request, inputs=input_files_now,
                                  protected_inputs=started["protected_inputs"], code=started["code"],
                                  output_inventory=inventory(staging, manifest))
                    atomic_json(staging / manifest, result)
                    archive = history / (run_id + "-previous")
                    os.replace(output, archive)
                    receipt = staging / "previous_manifest.receipt"
                    if receipt.is_file(): (archive / manifest).write_bytes(receipt.read_bytes())
                    os.replace(staging, output)
                    return result
                except BaseException as exc:
                    failure = dict(started, status="failed", error=str(exc), diagnostics=str(staging))
                    atomic_json(output / manifest, failure)
                    atomic_json(staging / "failure.json", failure)
                    raise
        return wrapped
    return decorate
