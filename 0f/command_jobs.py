"""Common execution of standalone .cmd files, including local LSF arrays."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
from submission_io import atomic_json, utcnow, receipt_lock
from lsf_submit import submit_attempt, require_submitted, submission_timeout
from submit_jobs import run_local
from workflow_mode import resolve_mode


class SubmissionError(RuntimeError):
    def __init__(self, message, result):
        super().__init__(message)
        self.result = result


def lsf_arguments(name, resources, stdout, stderr):
    return ['bsub', '-J', name, '-q', str(resources['queue']), '-n', str(resources['cores']),
            '-R', 'span[hosts=1]' +
            (f" rusage[mem={resources['memory_gib']}GB/host]" if resources.get('memory_gib') is not None else ''),
            '-oo', str(stdout), '-eo', str(stderr)]


def local_submit(cfg, path, expected_sha256=None):
    path = Path(path)
    raw = path.read_bytes(); digest = hashlib.sha256(raw).hexdigest()
    if expected_sha256 is not None and digest != expected_sha256:
        raise SubmissionError('Generated command changed before local execution',
                              dict(status='failed', job_id=None, command_sha256=digest))
    text = raw.decode('utf-8')
    match = re.search(r'^#BSUB -J .*\[([0-9,:-]+)\](?:%([1-9][0-9]*))?', text, re.M)
    indices, array_limit = None, None
    if match:
        indices = []; array_limit = int(match[2]) if match[2] else None
        for token in match[1].split(','):
            span, _, step = token.partition(':'); bounds = list(map(int, span.split('-')))
            indices.extend(range(bounds[0], bounds[-1]+1, int(step or 1)))
        if not indices: raise ValueError('Empty job array')
    required = re.search(r'rusage\[mem=([0-9.]+)([KMGT]?B)(/host)?\]', text)
    if required is None: raise RuntimeError('Cannot establish local memory request from .cmd')
    units = {'B': 1, 'KB': 1024, 'MB': 1024**2, 'GB': 1024**3, 'TB': 1024**4}
    job_bytes = float(required[1]) * units[required[2]]
    if not required[3] and cfg.get('lsf_memory_scope') == 'site_task':
        cores = re.search(r'^#BSUB -n ([1-9][0-9]*)$', text, re.M)
        if cores is None: raise RuntimeError('Per-task memory needs an explicit CPU request')
        job_bytes *= int(cores[1])
    from workflow_cpu import local_jobs
    jobs = local_jobs()
    if jobs < 1: raise ValueError('jobs_N must be positive')
    concurrent = min(jobs, len(indices), array_limit or jobs) if indices is not None else 1
    available = os.sysconf('SC_PHYS_PAGES') * os.sysconf('SC_PAGE_SIZE')
    if job_bytes * concurrent > available: raise RuntimeError('Requested memory exceeds local RAM')
    results = run_local([path], concurrent, indices)
    complete = len(results) == len(indices or [None]) and all(r['returncode'] == 0 for r in results)
    result = dict(schema=1, job_id='local', status='submitted' if complete else 'failed',
                  command_path=str(path), command_sha256=digest, results=results, submitted_at=utcnow())
    receipt = path.with_suffix('.local.json'); atomic_json(receipt, result)
    if not complete: raise SubmissionError('Local stage failed: ' + str(receipt), result)
    return result


def _submit_cmd(cfg, path, dependency=None, expected_sha256=None):
    """Compatibility API for generated #BSUB scripts; receipt logic is shared."""
    path = Path(path).expanduser().resolve()
    if resolve_mode() == 'local': return local_submit(cfg, path, expected_sha256)
    argv = ['bsub']
    if dependency is not None:
        if isinstance(dependency, int) and not isinstance(dependency, bool) and dependency > 0:
            dependency = f'done({dependency})'
        dependency = str(dependency)
        if not dependency.strip() or any(ord(c) < 32 for c in dependency): raise ValueError('Invalid dependency')
        argv += ['-w', dependency]
    timeout = min(60.0, submission_timeout(cfg.get('submit_timeout_seconds', 30)))
    receipt = path.with_suffix('.submission.json')
    history = json.loads(receipt.read_text()) if receipt.exists() else {'attempts': []}
    if history['attempts'] and history['attempts'][-1]['status'] in {'submitting', 'unknown', 'submitted'}:
        raise SubmissionError('Previous submission already accepted or uncertain; inspect ' + str(receipt), history['attempts'][-1])
    history['attempts'].append({})
    def persist(row):
        row = dict(row, schema=1, command_path=str(path), dependency=dependency, receipt_path=str(receipt))
        if row['job_id'] is not None: row['job_id'] = int(row['job_id'])
        history['attempts'][-1] = row
        atomic_json(receipt, history)
    result = submit_attempt(argv, path=path, expected_sha256=expected_sha256, timeout=timeout, persist=persist)
    try: require_submitted(result, receipt)
    except RuntimeError as exc: raise SubmissionError(str(exc), history['attempts'][-1]) from exc
    return history['attempts'][-1]


def _execute_script(path, mode, *, argv=None, receipt, metadata=None, plan_only=False):
    """Run an already defined script and store one durable submission receipt."""
    path, receipt = Path(path), Path(receipt)
    metadata = dict(metadata or {}, command=argv, script=str(path))
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if old.get('status') != 'planned': raise RuntimeError('Receipt already exists; inspect ' + str(receipt))
    atomic_json(receipt, dict(metadata, status='planned'))
    if plan_only: return dict(status='PLANNED', plan=str(path), jobs=1)
    if mode == 'local':
        atomic_json(receipt, dict(metadata, status='running', mode='local'))
        result = subprocess.run(['bash', str(path)], env=dict(os.environ, submit='local', SUBMIT='local', WORKFLOW_WORKER='1'))
        atomic_json(receipt, dict(metadata, status='complete' if result.returncode == 0 else 'failed', exit_code=result.returncode, mode='local'))
        return dict(status='DONE' if result.returncode == 0 else 'EXIT', exit_code=result.returncode)
    result = submit_attempt(argv or ['bsub'], path=path, persist=lambda row: atomic_json(receipt, dict(metadata, **row)))
    require_submitted(result, receipt)
    return result


def submit_cmd(cfg, path, dependency=None, expected_sha256=None):
    with receipt_lock(Path(path).with_suffix('.submission.json')):
        return _submit_cmd(cfg, path, dependency, expected_sha256)


def execute_script(path, mode, *, argv=None, receipt, metadata=None, plan_only=False):
    with receipt_lock(receipt):
        return _execute_script(path, mode, argv=argv, receipt=receipt, metadata=metadata, plan_only=plan_only)
