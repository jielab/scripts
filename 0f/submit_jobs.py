#!/usr/bin/env python3
"""Shared .cmd runner: a bounded local process pool or independent LSF jobs.

Local runs wait, record exit codes and stop scheduling new work on failure.
LSF runs record each accepted ID; --wait verifies DONE, never just disappearance.
"""
from __future__ import annotations
import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
import uuid


def bsub_process(argv, handle, timeout=30):
    return subprocess.run(argv, stdin=handle, capture_output=True, text=True,
                          encoding='utf-8', errors='replace', timeout=timeout,
                          env={**os.environ, 'LC_ALL': 'C'})


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
        record({'cmd': str(path), 'status': 'submitting', 'argv': argv})
        with path.open() as handle:
            process = bsub_process(argv, handle)
        match = re.search(r'Job <(\d+)>', process.stdout + process.stderr)
        result = {'cmd': str(path), 'job_id': match[1] if match else None,
                  'returncode': process.returncode, 'stdout': process.stdout, 'stderr': process.stderr}
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
