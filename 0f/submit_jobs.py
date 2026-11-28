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
from workflow_mode import resolve_mode
from submission_io import atomic_json
from lsf_submit import submission_timeout, bsub_process, submit_attempt, require_submitted


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



def local_foreground_one(path, index=None):
    """GU local workers belong to this terminal and clean up their descendants."""
    path = Path(path).resolve()
    env = {**os.environ, 'submit': 'local', 'SUBMIT': 'local',
           'WORKFLOW_WORKER': '1', 'SCRIPT_BACKGROUND': '0'}
    if index is not None:
        env['LSB_JOBINDEX'] = str(index)
    guard = [sys.executable, str(Path(__file__).with_name('process_lifetime.py')),
             '--parent-pid', str(os.getpid()), '--', 'bash', str(path)]
    process = subprocess.Popen(guard, env=env, start_new_session=False)
    try:
        code = process.wait()
    except BaseException:
        if process.poll() is None:
            process.terminate()
        process.wait()
        raise
    return {'cmd': str(path), 'index': index, 'returncode': code,
            'pid': process.pid, 'status': 'DONE' if code == 0 else 'EXIT'}

def run_local(paths, jobs=1, indices=None):
    from workflow_cpu import local_jobs
    jobs = local_jobs(jobs)
    if jobs < 1:
        raise ValueError('jobs_N must be positive')
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


def run(paths, mode=None, jobs=None, wait=False, merge=None, journal=None, *,
        merge_condition="done", state=None, state_schema="command.submission.1",
        plan=None, before_run=None, replace_pending_merge=False):
    """One command-batch executor; projects supply dependency and reuse policies."""
    mode = resolve_mode(mode)
    jobs = int(jobs or os.environ.get('jobs_N', '1'))
    if jobs < 1:
        raise ValueError('Invalid jobs_N')
    if merge_condition not in {'done', 'ended'}:
        raise ValueError('merge_condition must be done or ended')
    paths = [Path(p).resolve() for p in paths]
    merge = Path(merge).resolve() if merge else None
    if not paths or any(not p.is_file() for p in paths + ([merge] if merge else [])):
        raise ValueError('Every command file must exist')
    if plan:
        frozen = json.loads(Path(plan).read_text())
        if [str(p) for p in paths] != [j['command'] for j in frozen['jobs']] or str(merge) != frozen['merge']:
            raise ValueError('Commands differ from the exact plan')
    journal = Path(journal) if journal else paths[0].parent / ('submission-' + uuid.uuid4().hex[:12] + '.jsonl')
    if journal.exists():
        latest = {}
        for line in journal.read_text().splitlines():
            row = json.loads(line)
            if 'cmd' in row: latest[row['cmd']] = row
        if any(row.get('status') in {'submitting', 'uncertain', 'unknown'} for row in latest.values()):
            raise RuntimeError('Uncertain previous submission; inspect ' + str(journal))
    print('Submission journal: ' + str(journal), flush=True)
    def record(value):
        with journal.open('a') as out:
            out.write(json.dumps(value) + '\n'); out.flush(); os.fsync(out.fileno())
    state_path = Path(state) if state else None
    saved = json.loads(state_path.read_text()) if state_path and state_path.is_file() else {'schema': state_schema, 'workers': {}}
    if saved.get('schema') != state_schema:
        raise RuntimeError('Unrecognized submission state; use a clean output root')
    if saved.get('uncertain'):
        raise RuntimeError('Previous submission acknowledgement is uncertain; inspect the scheduler before retrying')
    def save():
        if state_path: atomic_json(state_path, saved)
    def status(job):
        result = subprocess.run(['bjobs', '-a', '-noheader', '-o', 'jobid stat', str(job)], capture_output=True, text=True, timeout=30)
        for line in result.stdout.splitlines():
            fields = line.split()
            if len(fields) == 2 and fields[0] == str(job): return fields[1]
        if 'not found' in (result.stdout + result.stderr).lower(): return None
        raise RuntimeError(f'Cannot query prior job {job}: {result.stdout} {result.stderr}')
    active = {'PEND', 'RUN', 'PSUSP', 'USUSP', 'SSUSP', 'WAIT', 'UNKWN', 'PROV'}
    for label, job in saved['workers'].items():
        if status(job) in active: raise RuntimeError(f'{label} already has active job {job}')
    if saved.get('merge'):
        previous = saved['merge']; observed = status(previous)
        if observed == 'PEND' and replace_pending_merge:
            result = subprocess.run(['bkill', str(previous)], capture_output=True, text=True, timeout=30)
            record({'job_id': previous, 'status': 'cancel_superseded_pending_merge', 'returncode': result.returncode})
            if result.returncode: raise RuntimeError('Could not cancel the previous pending merge')
        elif observed in active:
            raise RuntimeError('The prior merge is active')
    if before_run: before_run()
    saved = {'schema': state_schema, 'workers': {}, 'plan': str(plan) if plan else None}; save()
    if mode == 'local':
        results = run_local(paths, jobs)
        for result in results: record(result)
        if len(results) != len(paths) or any(r['returncode'] for r in results):
            raise RuntimeError('Local worker failed; remaining work and merge were not run')
        if merge:
            result = local_one(merge); record(result)
            if result['returncode']: raise RuntimeError('Merge failed')
        return results
    ids = []
    for n, path in enumerate(paths + ([merge] if merge else [])):
        is_merge = n == len(paths)
        argv = ['bsub']
        if is_merge: argv += ['-w', ' && '.join(f'{merge_condition}({i})' for i in ids)]
        def persist(row):
            record(dict(row, cmd=str(path), status='uncertain' if row['status'] == 'unknown' else row['status']))
            if row['status'] == 'submitted':
                saved.pop('uncertain', None)
                if is_merge: saved['merge'] = row['job_id']
                else: saved['workers'][path.stem] = row['job_id']
            else:
                saved['uncertain'] = {'cmd': str(path), 'journal': str(journal)}
            save()
        result = submit_attempt(argv, path=path, persist=persist, refresh_rejected=True)
        print((result['stdout'] + result['stderr']).strip(), flush=True)
        ids.append(require_submitted(result, journal))
    if wait:
        states = wait_lsf(ids[:len(paths)]); record({'states': states})
        if any(v != 'DONE' for v in states.values()): raise RuntimeError('An LSF worker failed')
        if merge:
            states = wait_lsf(ids[-1:]); record({'merge_states': states})
            if any(v != 'DONE' for v in states.values()): raise RuntimeError('LSF merge failed')
    return ids


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--submit', choices=['auto', 'local', 'hpc'], default='auto')
    p.add_argument('--jobs-N', type=int, default=int(os.environ.get('jobs_N', '1')))
    p.add_argument('--wait', action='store_true')
    p.add_argument('--merge')
    p.add_argument('--journal', type=Path)
    p.add_argument('commands', nargs='+')
    a = p.parse_args()
    try:
        run(a.commands, a.submit, a.jobs_N, a.wait, a.merge, a.journal)
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as e:
        p.exit(2, f'ERROR: {e}\n')
