#!/usr/bin/env python3
"""Shared .cmd runner: a bounded local process pool or independent LSF jobs.

Local runs wait, record exit codes and stop scheduling new work on failure.
LSF runs record each accepted ID; --wait verifies DONE, never just disappearance.
"""
from __future__ import annotations
import argparse
import concurrent.futures
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
import uuid


def submission_timeout():
    value = os.environ.get('WORKFLOW_BSUB_TIMEOUT_SECONDS', '300')
    try:
        timeout = float(value)
    except ValueError:
        raise ValueError('WORKFLOW_BSUB_TIMEOUT_SECONDS must be a finite positive number') from None
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('WORKFLOW_BSUB_TIMEOUT_SECONDS must be a finite positive number')
    return timeout


PENDING_RETRY = re.compile(r'[^\n]*Pending (?:job|job slots?) threshold reached[^\n]*Retrying in (\d+) seconds[^\n]*', re.I)


def bsub_process(argv, handle, timeout=None):
    timeout = submission_timeout() if timeout is None else timeout
    position = handle.tell()
    while True:
        result = _bsub_client(argv, handle, timeout)
        response = result.stdout + result.stderr
        # Some LSF installations reject the old credentials after a long
        # capacity wait. Refresh only after an explicit NOT-submitted reply;
        # a lost reply or a reported job ID must never trigger a new client.
        if (result.returncode and PENDING_RETRY.search(response)
                and 'User permission denied. Job not submitted.' in response
                and not re.search(r'Job <\d+>', response)):
            print('[LSF WAIT] LSF explicitly rejected the old request after waiting for capacity; '
                  'refreshing the submission client (Job not submitted).', flush=True)
            handle.seek(position)
            continue
        return result


def _bsub_client(argv, handle, timeout):
    # LSF keeps the same client alive when MAX_PEND_JOBS is reached. Its
    # explicit rejection/retry messages are progress, not a lost acknowledgement.
    # Keep that client (never launch a duplicate) while it waits for capacity.
    with subprocess.Popen(argv, stdin=handle, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True, encoding='utf-8', errors='replace',
                          env={**os.environ, 'LC_ALL': 'C'}) as process:
        deadline = time.monotonic() + timeout
        seen_retries = 0
        while True:
            try:
                stdout, stderr = process.communicate(timeout=max(.001, min(10, deadline - time.monotonic())))
                return subprocess.CompletedProcess(argv, process.returncode, stdout, stderr)
            except subprocess.TimeoutExpired as exc:
                partial = exc.stderr or b''
                if isinstance(partial, bytes):
                    partial = partial.decode('utf-8', errors='replace')
                retries = list(PENDING_RETRY.finditer(partial))
                if len(retries) > seen_retries:
                    for retry in retries[seen_retries:]:
                        print('[LSF WAIT] ' + retry[0].strip(), flush=True)
                    seen_retries = len(retries)
                    deadline = time.monotonic() + max(timeout, int(retries[-1][1]) + 30)
                if time.monotonic() >= deadline:
                    process.kill()
                    stdout, stderr = process.communicate()
                    raise subprocess.TimeoutExpired(argv, timeout, output=stdout, stderr=stderr) from None
            except BaseException:
                process.kill()
                process.communicate()
                raise


def local_one(path, index=None):
    path = Path(path).resolve()
    suffix = f'.{index}' if index is not None else ''
    env = {**os.environ, 'submit': 'local', 'SUBMIT': 'local', 'WORKFLOW_WORKER': '1'}
    # This is only an array index, never a fake LSB_JOBID or allocation.
    if index is not None:
        env['LSB_JOBINDEX'] = str(index)
    with Path(str(path) + suffix + '.out').open('a') as out, Path(str(path) + suffix + '.err').open('a') as err:
        process = subprocess.Popen(['bash', str(path)], stdout=out, stderr=err,
                                   env=env, start_new_session=True)
        try:
            code = process.wait()
        except BaseException:
            os.killpg(process.pid, signal.SIGTERM)
            process.wait()
            raise
    return {'cmd': str(path), 'index': index, 'returncode': code, 'pid': process.pid,
            'status': 'DONE' if code == 0 else 'EXIT'}


def run_local(paths, jobs=1, indices=None):
    if not 1 <= jobs <= 8:
        raise ValueError('jobs_N must be 1..8')
    tasks = iter((p, i) for p in paths for i in (indices or [None]))
    results, pending = [], {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        for _ in range(jobs):
            task = next(tasks, None)
            if task is not None:
                pending[pool.submit(local_one, *task)] = task
        failed = False
        while pending:
            done, _ = concurrent.futures.wait(pending, return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                pending.pop(future)
                result = future.result()
                results.append(result)
                failed |= result['returncode'] != 0
                print(f"[{result['status']}] {result['cmd']} index={result['index']}", flush=True)
            if not failed:
                for _ in done:
                    task = next(tasks, None)
                    if task is not None:
                        pending[pool.submit(local_one, *task)] = task
    return results


def wait_lsf(ids):
    states = {}
    while len(states) < len(ids):
        remaining = [i for i in ids if i not in states]
        probe = subprocess.run(['bjobs', '-a', '-noheader', '-o', 'jobid stat', *remaining],
                               capture_output=True, text=True, timeout=30)
        seen = set()
        for line in probe.stdout.splitlines():
            fields = line.split()
            if len(fields) == 2 and fields[0] in remaining:
                seen.add(fields[0])
                if fields[1] in {'DONE', 'EXIT'}:
                    states[fields[0]] = fields[1]
        if set(remaining) - seen:
            raise RuntimeError('Cannot verify some job states; inspect bjobs/bacct. Missing jobs are NOT considered successful.')
        if len(states) < len(ids):
            time.sleep(5)
    return states


def run(paths, mode=None, jobs=None, wait=False, merge=None):
    mode = mode or os.environ.get('submit', 'local')
    jobs = int(jobs or os.environ.get('jobs_N', '8'))
    if mode not in {'local', 'hpc'}:
        raise ValueError('submit must be local or hpc')
    limit = 10000 if mode == 'hpc' else 8
    if not 1 <= jobs <= limit:
        raise ValueError(f'jobs_N must be 1..{limit} for {mode}')
    paths = [Path(p).resolve() for p in paths]
    if not paths or any(not p.is_file() for p in paths):
        raise ValueError('Every command file must exist')
    timeout = submission_timeout() if mode == 'hpc' else None
    journal = paths[0].parent / ('submission-' + uuid.uuid4().hex[:12] + '.jsonl')
    print(f'Submission journal: {journal}', flush=True)
    def record(value):
        with journal.open('a') as out:
            out.write(json.dumps(value) + '\n'); out.flush(); os.fsync(out.fileno())
    if mode == 'local':
        results = run_local(paths, jobs)
        for result in results:
            record(result)
        if any(r['returncode'] for r in results):
            raise RuntimeError('Local command failed; later unscheduled work was stopped. See .cmd.err and journal.')
        if merge:
            result = local_one(merge); record(result)
            if result['returncode']:
                raise RuntimeError('Merge failed')
        return results
    ids = []
    for n, path in enumerate(paths + ([Path(merge).resolve()] if merge else [])):
        argv = ['bsub']
        if n == len(paths):
            argv += ['-w', ' && '.join(f'done({i})' for i in ids)]
        elif n >= jobs:
            argv += ['-w', f'ended({ids[n-jobs]})']
        record({'cmd': str(path), 'status': 'submitting', 'argv': argv,
                'timeout_seconds': timeout})
        try:
            with path.open() as handle:
                process = bsub_process(argv, handle, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            def decoded(value):
                return value.decode('utf-8', errors='replace') if isinstance(value, bytes) else (value or '')
            stdout, stderr = decoded(exc.stdout), decoded(exc.stderr)
            match = re.search(r'Job <(\d+)>', stdout + stderr)
            record({'cmd': str(path), 'status': 'uncertain', 'argv': argv,
                    'job_id': match[1] if match else None, 'timeout_seconds': timeout,
                    'stdout': stdout, 'stderr': stderr})
            raise RuntimeError(
                f'LSF submission acknowledgement timed out after {timeout:g}s for {path.name}; '
                f'the job may already be accepted. Inspect bjobs and {journal} before retrying. '
                'No automatic resubmission was attempted. '
                'WORKFLOW_BSUB_TIMEOUT_SECONDS controls this wait.'
            ) from exc
        match = re.search(r'Job <(\d+)>', process.stdout + process.stderr)
        result = {'cmd': str(path), 'job_id': match[1] if match else None,
                  'returncode': process.returncode, 'stdout': process.stdout, 'stderr': process.stderr,
                  'status': 'submitted' if process.returncode == 0 and match else 'uncertain'}
        record(result)
        print((process.stdout + process.stderr).strip(), flush=True)
        if process.returncode or not match:
            raise RuntimeError(f'Submission stopped; inspect {journal} before retrying')
        ids.append(match[1])
    if wait:
        # Wait parents first: a failed parent leaves a done()-dependent merge pending.
        states = wait_lsf(ids[:len(paths)])
        record({'states': states})
        if any(v != 'DONE' for v in states.values()):
            raise RuntimeError('An LSF worker failed; dependent merge was not run')
        if merge:
            states = wait_lsf(ids[-1:]); record({'merge_states': states})
            if any(v != 'DONE' for v in states.values()):
                raise RuntimeError('LSF merge failed')
    return ids


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--submit', choices=['local', 'hpc'], default=os.environ.get('submit', 'local'))
    p.add_argument('--jobs-N', type=int, default=int(os.environ.get('jobs_N', '8')))
    p.add_argument('--wait', action='store_true')
    p.add_argument('--merge')
    p.add_argument('commands', nargs='+')
    a = p.parse_args()
    try:
        run(a.commands, a.submit, a.jobs_N, a.wait, a.merge)
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as e:
        p.exit(2, f'ERROR: {e}\n')
