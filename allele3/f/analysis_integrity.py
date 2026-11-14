#!/usr/bin/env python3
"""Content identities and atomic, locked publication for derived analyses.

Permanent inputs are read only. A failed refresh invalidates current evidence,
while preserving the previous files and receipt in an attempt history.
"""
from __future__ import annotations

# Analysis Integrity
# ======================================================================
import copy
from dataclasses import dataclass
import fcntl
import functools
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import uuid


@dataclass(frozen=True)
class ModuleLayout:
    """Public results and sibling support directories for one analysis module."""
    root: Path
    key: Path

    @property
    def output(self): return self.root / self.key
    @property
    def runtime(self): return self.root / 'tmp' / self.key / 'work'
    @property
    def logs(self): return self.root / 'logs' / self.key
    @property
    def jobs(self): return self.root / 'jobs' / self.key
    @property
    def qc(self): return self.root / 'qc' / self.key


PUBLIC_MODULES = {'pca3d', 'herv_census', 'gnomad_scan', 'hotspots', 'density',
                  'herv_enrichment', 'base_occupancy', 'base_simulation'}


def module_layout(output, analysis_root=None):
    """Resolve paths without creating directories or loading genomic data."""
    output = Path(output).expanduser().absolute()
    configured = analysis_root or os.environ.get('ALLELE_ANALYSIS_ROOT')
    root = Path(configured).expanduser().absolute() if configured else None
    if root is None or not output.is_relative_to(root):
        module = next((p for p in reversed([output, *output.parents]) if p.name in PUBLIC_MODULES), None)
        root = module.parent if module else output.parent
    key = output.relative_to(root)
    if not key.parts or key.parts[0] in {'logs', 'jobs', 'qc', 'tmp'} or '..' in key.parts:
        raise ValueError('Choose a module result directory, not an analysis/support root: ' + str(output))
    return ModuleLayout(root, key)


def resolve_result_input(path):
    """Resolve published results/bases to their receipt-bound native artifacts.

    A layout receipt maps public report paths to their native artifacts.
    It never invents scientific completion.
    """
    if not path or not isinstance(path, (str, Path)):
        return path
    candidate = Path(path).expanduser().absolute()
    for public in [candidate, *candidate.parents]:
        if not any(p.name in PUBLIC_MODULES for p in [public, *public.parents]):
            continue
        try:
            layout = module_layout(public)
        except ValueError:
            # Inputs may already be inside tmp/MODULE/work. Only public
            # result directories have layout receipts to resolve here.
            continue
        receipt = layout.qc / 'layout.json'
        if not receipt.is_file():
            continue
        saved = json.loads(receipt.read_text())
        if saved.get('schema') != 'allele3-module-layout-1' or saved.get('public_root') != str(public):
            raise ValueError('Invalid module layout receipt: ' + str(receipt))
        native = Path(saved['runtime_root'])
        if native != layout.runtime:
            raise ValueError('Unexpected runtime location in ' + str(receipt))
        relative = str(candidate.relative_to(public))
        mapped = saved.get('public_files', {}).get(relative, {}).get('native')
        if mapped:
            target = native / mapped
        else:
            aliases = saved.get('directory_map', {})
            target = native / aliases.get(relative, relative)
        if target.is_relative_to(native) and target.exists():
            return target
    return path


def resolve_cli_inputs(argv):
    """Translate input artifacts; output destinations remain public API paths."""
    result = list(argv)
    inputs = {'--scan-dir', '--store', '--modern-dir', '--dir-modern', '--basis-dir', '--reuse-basis-dir',
              '--metadata', '--coordinates', '--eigenvec', '--source-html', '--time-dir',
              '--alleles-dir', '--allofus-dir', '--atlas-dir', '--summary', '--tiles',
              '--catalog', '--targets', '--sites', '--observations', '--genotypes',
              '--coordinate-receipt', '--provenance-manifest', '--reference-metadata',
              '--target-info', '--source', '--reuse-from', '--pca-basis'}
    for i, arg in enumerate(result):
        if arg in inputs and i+1 < len(result):
            result[i+1] = str(resolve_result_input(result[i+1]))
        elif '=' in arg and arg.split('=', 1)[0] in inputs:
            key, value = arg.split('=', 1)
            result[i] = key + '=' + str(resolve_result_input(value))
    return result

TERMINAL = {
    'complete', 'complete_descriptive_imputed_panel', 'complete_descriptive_all_allele_panel',
    'complete_descriptive_unverified', 'complete_adjustment_unavailable', 'planned',
    'no_usable_exact_high_quality_contrast', 'no_usable_contrast', 'not_estimable',
    'not_available', 'not_run', 'complete_with_exclusions',
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
    index = 1
    while True:
        temporary = path.with_name(path.name + ('.tmp' if index == 1 else f'.tmp{index}'))
        try:
            stream = temporary.open('x')
            break
        except FileExistsError:
            index += 1
    try:
        with stream:
            stream.write(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False, default=str) + '\n')
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)



def code_identity():
    root = Path(__file__).resolve().parent.parent
    files = {str(p.relative_to(root)): {'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns}
             for p in sorted(root.rglob('*')) if p.is_file() and p.suffix in {'.py', '.sh', '.env', '.yml', '.html', '.js'} and '__pycache__' not in p.parts}
    return {'files': files, 'identity_policy': 'file_size_mtime; no routine source hashing',
            'python': sys.version, 'executable': sys.executable, 'git_commit': 'unavailable_local_checkout_without_git'}


def inventory(root, manifest):
    root = Path(root)
    return {str(p.relative_to(root)): {k: v for k, v in identity(p).items() if k != 'path'}
            for p in sorted(root.rglob('*')) if p.is_file() and p != root / manifest
            and not set(p.relative_to(root).parts[:-1]) & {"jobs", "logs"}}


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
        actual = identity(path)
        keys = ('size',)
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
            lock_dir = output.parent / "qc/locks"
            lock_dir.mkdir(parents=True, exist_ok=True)
            with (lock_dir / (output.name + ".lock")).open("a") as lock:
                try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError as exc: raise RuntimeError(f"Output in use: {output}") from exc
                # Re-evaluate aliases/dynamic source paths while holding the lock.
                protected = declared_inputs(raw, input_dirs=input_dirs, input_files=input_files,
                                            input_prefixes=input_prefixes, resolver=resolve_inputs)
                check_output_paths(output, protected)
                input_files_now = {str(Path(p.path)): identity(p.path) for p in protected if Path(p.path).is_file()}
                run_id = uuid.uuid4().hex
                history = output.parent / "qc" / "history" / output.name
                history.mkdir(parents=True, exist_ok=True)
                sequence = 1
                while (history / str(sequence)).exists():
                    sequence += 1
                attempt = history / str(sequence)
                attempt.mkdir()
                staging = attempt / "staging"
                staging.mkdir()
                old_manifest = output / manifest
                if old_manifest.is_file():
                    (attempt / "previous_manifest.receipt").write_bytes(old_manifest.read_bytes())
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
                    if Path(manifest).parent == Path("qc"):
                        result["inventory_base"] = "module"
                    atomic_json(staging / manifest, result)
                    archive = attempt / "previous"
                    # Module job ledgers and open logs survive atomic scientific publication.
                    for support in ('jobs', 'logs'):
                        source = output / support
                        if source.exists():
                            if (staging / support).exists():
                                raise ValueError('Scientific stage wrote reserved support directory: ' + support)
                            os.replace(source, staging / support)
                    os.replace(output, archive)
                    receipt = attempt / "previous_manifest.receipt"
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


# Output Layout
# ======================================================================
_description_output_layout = 'Publish short, explicit result paths without altering scientific receipts.\n\nOnly aggregate tables selected below and scan summary figures are published.\nIndividual-level exports and engine intermediates are never auto-discovered.\n'

import html
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / '0f'))
from workflow_layout import publish_result

GENOME_RESULTS = {
    'windows': ('nominated_windows.tsv', 'all_windows.tsv.gz', 'loci.tsv', 'loci.bed', 'calibration.tsv'),
    'characterize': ('selected_windows.tsv', 'window_matches.tsv', 'window_match_balance.tsv',
                    'window_annotations.tsv', 'annotation_contrasts.tsv', 'site_sampling_audit.tsv'),
    'impact': ('window_atlas_summary.tsv', 'annotation_contrasts.tsv', 'window_match_balance.tsv'),
    'selection': ('window_selection_summary.tsv', 'window_selection_contrasts.tsv'),
}
SCAN_RESULTS = {
    'chromosome_summary.tsv', 'population_summary.tsv', 'population_classes.tsv',
    'window_summary.tsv', 'context_summary.tsv', 'expected_context.tsv', 'occupancy_thresholds.tsv',
    'annotation_counts.tsv', 'frequency_counts.tsv', 'population_counts.tsv', 'relationship_counts.tsv',
    'singleton_sensitivity.tsv', 'site_ac_spectrum.tsv', 'enrichment_summary.tsv',
    'ld_block_ranked.tsv', 'ld_block_enrichment.tsv', 'erv_enrichment.tsv',
    'ld_hotspot_consensus.tsv', 'ld_hotspot_consensus_regions.tsv',
}
FIGURES = {'.png', '.pdf', '.svg', '.xlsx'}


def write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.layout-', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            handle.write(text)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def publish_hotspots(out):
    """Expose selected aggregate outputs; scientific stage files stay in runtime."""
    out = Path(out)
    for stage, names in GENOME_RESULTS.items():
        for name in names:
            source = out / 'runtime' / stage / name
            if source.is_file():
                publish_result(source, out / 'results' / (stage + '.' + name))


def publish_scan(out):
    out = Path(out)
    summary = out / 'summary'
    if not summary.is_dir():
        raise FileNotFoundError(summary)
    for name in sorted(SCAN_RESULTS):
        source = summary / name
        if source.is_file():
            publish_result(source, out / 'results' / name)
    figures = []
    receipt = out / 'qc/published_scan.json'
    previous = json.loads(receipt.read_text()).get('files', []) if receipt.is_file() else []
    for source in sorted(summary.iterdir()):
        if source.is_file() and source.suffix.lower() in FIGURES:
            publish_result(source, out / 'reports' / source.name)
            figures.append(source.name)
    current = ['reports/' + name for name in figures]
    current += ['results/' + name for name in sorted(SCAN_RESULTS) if (summary / name).is_file()]
    if (summary / 'REPORT.md').is_file():
        publish_result(summary / 'REPORT.md', out / 'reports/README.md')
        current.append('reports/README.md')
    for relative in set(previous) - set(current):
        path = Path(relative)
        if len(path.parts) != 2 or path.parts[0] not in {'results', 'reports'}:
            raise ValueError('Unsafe publication receipt entry: ' + relative)
        (out / path).unlink(missing_ok=True)
    write_text(receipt, json.dumps({'files': current}, indent=2) + '\n')
    links = [(name, name) for name in figures]
    links += [('../results/' + p.name, p.name) for p in sorted((out / 'results').glob('*'))]
    write_text(out / 'reports/index.html', '<!doctype html><meta charset="utf-8"><title>gnomAD scan</title>'
        '<h1>gnomAD scan results</h1><ul>' + ''.join(
        '<li><a href="' + html.escape(url, quote=True) + '">' + html.escape(label) + '</a></li>'
        for url, label in links) + '</ul>')


def publish_linked_report(module, source_root, document, name):
    """Publish the aggregate files explicitly linked by a module's report."""
    import re
    from urllib.parse import urlsplit, unquote
    module, source_root = Path(module), Path(source_root).resolve()
    used = {}
    def replace(match):
        attribute, quote, url = match.groups()
        parsed = urlsplit(html.unescape(url))
        if parsed.scheme or parsed.netloc or not parsed.path:
            return match.group(0)
        source = (source_root / unquote(parsed.path)).resolve()
        if source.is_relative_to(source_root) and source.is_file() and source.suffix.lower() in FIGURES | {'.tsv', '.csv', '.bed', '.gz'}:
            rel = source.relative_to(source_root)
            folder = 'reports' if source.suffix.lower() in FIGURES else 'results'
            destination = module / folder / '.'.join(rel.parts)
            if destination in used and used[destination] != source:
                raise ValueError('Colliding report export: ' + str(destination))
            used[destination] = source
            publish_result(source, destination)
            target = os.path.relpath(destination, module / 'reports')
        else:
            target = os.path.relpath(source, module / 'reports')
        if parsed.fragment:
            target += '#' + parsed.fragment
        return attribute + '=' + quote + html.escape(target, quote=True) + quote
    document = re.sub(r'''(href|src)=(["'])(.*?)\2''', replace, document)
    destination = module / 'reports' / name
    write_text(destination, document)
    return destination


# Fresh Run
# ======================================================================
_description_fresh_run = 'Fresh output preparation without deleting previous calculations or inputs.'

# Shared root defaults come from 0f/0os.sh.
import sys as _paths_sys
from pathlib import Path as _PathsPath
_paths_root = next(p for p in _PathsPath(__file__).resolve().parents if (p / "0f" / "workflow_paths.py").is_file())
_paths_sys.path.insert(0, str(_paths_root / "0f"))
from workflow_paths import configured_path
from contextlib import contextmanager


def boolean(value):
    value=str(value).upper()
    if value not in {'TRUE','FALSE'}:raise ValueError('--replace must be TRUE or FALSE')
    return value=='TRUE'


@contextmanager
def output_run(output, replace, protected=()):
    output=Path(output).expanduser().resolve()
    analysis=Path(os.environ.get('ALLELE_ANALYSIS_ROOT',configured_path('ANALYSIS_ROOT', 'allele3'))).resolve()
    if output==analysis or output in analysis.parents:raise ValueError('Choose a separate task output, not the analysis root or its parent')
    for original in protected:
        if not original:continue
        original=Path(original).expanduser().resolve()
        if output==original or original in output.parents or output in original.parents:
            raise ValueError(f'Output overlaps a preserved input/cache: {original}; choose a new --out')
    output.parent.mkdir(parents=True,exist_ok=True)
    with (output.parent/('.'+output.name+'.recalculation.lock')).open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc:raise RuntimeError('Another calculation owns this output directory') from exc
        held=[]
        try:
            for name in ('job.lock','run.lock','.allele3.lock','.pca3d.lock'):
                p=output/name
                if p.exists():
                    h=p.open('a');held.append(h)
                    try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
                    except BlockingIOError as exc:raise RuntimeError('Another job is active: '+str(p)) from exc
            old=[p for p in output.iterdir() if p.name not in {'jobs','logs','job.lock','run.lock','.allele3.lock','.pca3d.lock'}] if output.is_dir() else []
            if replace and old:
                history = output.parent / 'tmp' / output.name
                history.mkdir(parents=True, exist_ok=True)
                sequence = 1
                while (history / str(sequence)).exists():
                    sequence += 1
                archive = history / str(sequence)
                archive.mkdir()
                for p in old:p.rename(archive/p.name)
                print('[preserved previous output] '+str(archive),flush=True)
            output.mkdir(parents=True,exist_ok=True)
        finally:
            for handle in held:handle.close()
        print('[mode] '+('FRESH RECALCULATION (--replace TRUE)' if replace else 'CHECKPOINT RESUME (--replace FALSE)'),flush=True)
        print('[new output] '+str(output),flush=True)
        for original in protected:
            if original:print('[permanent input/cache; retained] '+str(original),flush=True)
        yield output


