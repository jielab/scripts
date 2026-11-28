#!/usr/bin/env python3
"""Shared LSF client, acknowledgement parsing and durable submission attempts."""
from __future__ import annotations
import hashlib
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from submission_io import atomic_json, utcnow, receipt_lock

CAPACITY_RETRY = re.compile(r'[^\n]*(?:Pending (?:job|job slots?)|Active job) threshold reached[^\n]*Retrying in (\d+) seconds[^\n]*', re.I)


def decoded(value):
    return value.decode('utf-8', errors='replace') if isinstance(value, bytes) else value or ''


def submission_timeout(value=None):
    value = os.getenv('WORKFLOW_BSUB_TIMEOUT_SECONDS', '300') if value is None else value
    try:
        timeout = float(value)
    except (ValueError, TypeError):
        raise ValueError('Submission timeout must be a finite positive number') from None
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('Submission timeout must be a finite positive number')
    return timeout


def job_ids(stdout, stderr):
    return set(re.findall(r'\bJob\s*<([1-9][0-9]*)>', decoded(stdout) + '\n' + decoded(stderr)))


def run_bsub(argv, *, timeout=None, wait_capacity=True, **kwargs):
    """One client; capacity progress extends its deadline, never duplicates it."""
    timeout = submission_timeout(timeout)
    if not kwargs.get('capture_output'):
        return subprocess.run(argv, timeout=timeout, **kwargs)
    kwargs.pop('capture_output')
    payload = kwargs.pop('input', None)
    check = kwargs.pop('check', False)
    if payload is not None:
        if kwargs.get('stdin') is not None:
            raise ValueError('Use stdin or input, not both')
        kwargs['stdin'] = subprocess.PIPE
    kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with subprocess.Popen(argv, **kwargs) as process:
        deadline, seen, first = time.monotonic() + timeout, 0, True
        try:
            while True:
                try:
                    supplied = payload if first else None
                    first = False
                    stdout, stderr = process.communicate(input=supplied, timeout=max(.001, min(10, deadline-time.monotonic())))
                    result = subprocess.CompletedProcess(argv, process.returncode, stdout, stderr)
                    if check:
                        result.check_returncode()
                    return result
                except subprocess.TimeoutExpired as exc:
                    retries = list(CAPACITY_RETRY.finditer(decoded(exc.stderr))) if wait_capacity else []
                    if len(retries) > seen:
                        for retry in retries[seen:]:
                            print('[LSF WAIT] ' + retry[0].strip(), flush=True)
                        seen = len(retries)
                        deadline = time.monotonic() + max(timeout, int(retries[-1][1]) + 30)
                    if time.monotonic() >= deadline:
                        process.kill()
                        stdout, stderr = process.communicate()
                        raise subprocess.TimeoutExpired(argv, timeout, output=stdout, stderr=stderr) from None
        except BaseException:
            if process.poll() is None:
                process.kill()
                process.communicate()
            raise


def bsub_process(argv, handle, timeout=None):
    """Compatibility for .cmd callers; retry only an explicit NOT-submitted reply."""
    position = handle.tell()
    while True:
        result = run_bsub(argv, stdin=handle, timeout=timeout, capture_output=True,
                          text=True, encoding='utf-8', errors='replace', env={**os.environ, 'LC_ALL': 'C'})
        response = result.stdout + result.stderr
        if (result.returncode and CAPACITY_RETRY.search(response)
                and 'User permission denied. Job not submitted.' in response
                and not job_ids(result.stdout, result.stderr)):
            handle.seek(position)
            continue
        return result


def named_argv(argv, raw):
    name = os.getenv('WORKFLOW_JOB_NAME')
    argv = list(argv)
    if not name: return argv
    if not re.fullmatch(r'[A-Za-z0-9_.-]{1,120}', name): raise ValueError('Invalid workflow job name')
    old = ''
    if '-J' in argv:
        index = argv.index('-J'); old = argv[index+1]; del argv[index:index+2]
    else:
        match = re.search(r'^#BSUB\s+-J\s+(.+)$', decoded(raw), re.M)
        if match: old = match[1]
    array = re.search(r'(\[[0-9,:-]+\](?:%[0-9]+)?)', old)
    return argv + ['-J', name + (array[1] if array else '')]


def submit_attempt(argv, *, persist, path=None, script=None, timeout=None,
                   expected_sha256=None, wait_capacity=True, refresh_rejected=False):
    """Persist intent before bsub and its result afterwards; ambiguous IDs survive.

    Receipt adapters preserve existing schemas. A failed intent write prevents
    submission; a failed result write leaves the durable intent blocking retries.
    """
    if (path is None) == (script is None):
        raise ValueError('Supply one command path or script')
    timeout = submission_timeout(timeout)
    record = dict(status='submitting', argv=list(argv), submitted_at=utcnow(),
                  job_id=None, returncode=None, stdout='', stderr='', timeout_seconds=timeout)
    handle = Path(path).open('rb') if path is not None else None
    try:
        raw = handle.read() if handle else (script.encode('utf-8') if isinstance(script, str) else script)
        argv = named_argv(argv, raw)
        record['argv'] = argv
        record['command_sha256'] = hashlib.sha256(raw).hexdigest()
        if expected_sha256 is not None and record['command_sha256'] != expected_sha256:
            record.update(status='failed', error='Generated command checksum changed; not submitted')
            persist(dict(record))
            return record
        persist(dict(record))
        if handle:
            handle.seek(0)
        try:
            if refresh_rejected and handle:
                result = bsub_process(argv, handle, timeout)
            else:
                source = {'stdin': handle} if handle else {'input': decoded(raw)}
                result = run_bsub(argv, **source, timeout=timeout, wait_capacity=wait_capacity,
                                  capture_output=True, text=True, encoding='utf-8', errors='replace',
                                  env={**os.environ, 'LC_ALL': 'C'})
            record.update(returncode=result.returncode, stdout=decoded(result.stdout), stderr=decoded(result.stderr))
        except subprocess.TimeoutExpired as exc:
            record.update(stdout=decoded(exc.stdout), stderr=decoded(exc.stderr), error='LSF acknowledgement timed out; acceptance unknown')
        except OSError as exc:
            record.update(status='failed', error=str(exc))
        ids = job_ids(record['stdout'], record['stderr'])
        if len(ids) == 1:
            record['job_id'] = next(iter(ids))
        if record['returncode'] == 0 and len(ids) == 1:
            record['status'] = 'submitted'
        elif record['returncode'] not in (None, 0) and not ids:
            record['status'] = 'failed'
        elif record['status'] != 'failed':
            record['status'] = 'unknown'
        record['finished_at'] = utcnow()
        persist(dict(record))
        return record
    finally:
        if handle:
            handle.close()


def require_submitted(result, receipt):
    if result['status'] != 'submitted':
        raise RuntimeError(f"LSF status={result['status']}, job_id={result.get('job_id')}; inspect {receipt} before retrying. "
                           + result.get('error', '') + '\n' + result['stdout'] + result['stderr'])
    return result['job_id']


def submit_graph(tasks, saved, *, persist, receipt, timeout=None):
    """Submit a topologically ordered graph using confirmed IDs and one transaction.

    Tasks define key, path, optional argv/checksum and (parent, done|ended) edges.
    persist(key, row) adapts the canonical receipt to existing on-disk schemas.
    """
    if any(row.get('status') in {'submitting', 'unknown'} for row in saved.values()):
        raise RuntimeError('An earlier submission is uncertain; inspect ' + str(receipt))
    results = dict(saved)
    for task in tasks:
        key = task['key']
        old = results.get(key, {})
        if old.get('status') == 'submitted':
            if task.get('sha256') and old.get('command_sha256') != task['sha256']:
                raise ValueError('Recorded command identity changed: ' + key)
            continue
        terms = []
        for parent, condition in task.get('depends', []):
            if condition not in {'done', 'ended'}:
                raise ValueError('Unsupported dependency condition')
            row = results.get(parent, {})
            if row.get('status') != 'submitted' or not row.get('job_id'):
                raise RuntimeError('Dependency has no confirmed scheduler ID: ' + parent)
            terms.append(f"{condition}({row['job_id']})")
        dependency = ' && '.join(dict.fromkeys(terms))
        argv = list(task.get('argv', ['bsub']))
        if dependency:
            argv += ['-w', dependency]
        def save(row):
            row = dict(row, dependency=dependency)
            persist(key, row)
            results[key] = row
        result = submit_attempt(argv, path=task['path'], expected_sha256=task.get('sha256'),
                                timeout=timeout, persist=save)
        job_id = require_submitted(result, receipt)
        print(f'[SUBMITTED] {key} {job_id}', flush=True)
    return results


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--script', required=True, type=Path)
    p.add_argument('options', nargs=argparse.REMAINDER)
    a = p.parse_args()
    options = a.options[1:] if a.options[:1] == ['--'] else a.options
    receipt = a.script.with_suffix('.submission.json')
    with receipt_lock(receipt):
        if receipt.exists():
            p.error(f'Submission receipt already exists: {receipt}; inspect it before retrying')
        result = submit_attempt(['bsub', *options], path=a.script, persist=lambda row: atomic_json(receipt, row))
    print(result['stdout'], end='')
    print(result['stderr'], end='', file=sys.stderr)
    require_submitted(result, receipt)


if __name__ == '__main__':
    main()
