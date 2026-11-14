#!/usr/bin/env python3
"""Shared .cmd runner: a bounded local process pool or independent LSF jobs.

Local runs wait, record exit codes and stop scheduling new work on failure.
LSF runs record each accepted ID; --wait verifies DONE, never just disappearance.
"""
from __future__ import annotations
import argparse
import contextlib
import fcntl
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
def resolve_mode(mode=None):
    import shutil
    mode = mode or os.environ.get('SUBMIT', 'auto')
    if mode == 'auto':
        return 'hpc' if shutil.which('bsub') else 'local'
    return mode


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


def run(paths, mode=None, jobs=None, wait=False, merge=None, journal=None, state=None, plan=None):
    mode=resolve_mode(mode);jobs=int(jobs or 1)
    if mode not in {'local','hpc'} or not 1 <= jobs <= (10000 if mode=='hpc' else 8):
        raise ValueError('Invalid submit mode or jobs_N (local 1..8, HPC 1..10000)')
    paths=[Path(p).resolve() for p in paths]
    if not paths or any(not p.is_file() for p in paths):raise ValueError('Every command file must exist')
    if merge and not Path(merge).is_file():raise ValueError('Merge command does not exist')
    if plan:
        frozen=json.loads(Path(plan).read_text())
        if [str(p) for p in paths]!=[j['command'] for j in frozen['jobs']] or str(merge)!=frozen['merge']:
            raise ValueError('Submission command set differs from the exact plan')
    journal=Path(journal) if journal else paths[0].parent/('submission-'+uuid.uuid4().hex[:12]+'.jsonl')
    print('Submission journal: '+str(journal),flush=True)
    def record(value):
        with journal.open('a') as stream:
            stream.write(json.dumps(value)+'\n');stream.flush();os.fsync(stream.fileno())
    state_path=Path(state) if state else None
    saved=json.loads(state_path.read_text()) if state_path and state_path.is_file() else {'schema':'ld.submission.1','workers':{}}
    if saved.get('schema')!='ld.submission.1':raise RuntimeError('Unrecognized submission state; use a clean output root for this release')
    def save():
        if state_path:
            state_path.parent.mkdir(parents=True,exist_ok=True)
            temporary=state_path.with_suffix('.tmp')
            with temporary.open('w') as stream:
                json.dump(saved,stream,indent=2);stream.flush();os.fsync(stream.fileno())
            temporary.replace(state_path)
    if saved.get('uncertain'):
        raise RuntimeError('Previous submission acknowledgement is uncertain; inspect bjobs and the submission journal before retrying')
    def status(job):
        result=subprocess.run(['bjobs','-a','-noheader','-o','jobid stat',str(job)],capture_output=True,text=True,timeout=30)
        for line in result.stdout.splitlines():
            fields=line.split()
            if len(fields)==2 and fields[0]==str(job):return fields[1]
        if 'not found' in (result.stdout+result.stderr).lower():return None
        raise RuntimeError(f'Cannot query prior job {job}: {result.stdout} {result.stderr}')
    active={'PEND','RUN','PSUSP','USUSP','SSUSP','WAIT','UNKWN','PROV'}
    if mode=='hpc' or saved['workers'] or saved.get('merge'):
        for label,job in saved['workers'].items():
            if status(job) in active:raise RuntimeError(f'{label} already has active job {job}; finish the current plan before replacing the same output root')
        previous=saved.get('merge')
        if previous:
            observed=status(previous)
            if observed=='PEND':
                result=subprocess.run(['bkill',str(previous)],capture_output=True,text=True,timeout=30)
                record({'job_id':previous,'status':'cancel_superseded_pending_merge','returncode':result.returncode})
                if result.returncode:raise RuntimeError("Could not cancel this workflow's prior pending merge")
            elif observed in active:raise RuntimeError('The prior merge is active; wait for it before replacing these outputs')
    if plan:
        from merge import activate_plan
        activate_plan(plan)
    saved={'schema':'ld.submission.1','workers':{},'plan':str(plan) if plan else None};save()
    if mode=='local':
        results=run_local(paths,jobs)
        for result in results:record(result)
        if len(results)!=len(paths) or any(r['returncode'] for r in results):
            raise RuntimeError('Local worker failed; later unscheduled work and merge were not run')
        if merge:
            result=local_one(merge);record(result)
            if result['returncode']:raise RuntimeError('Merge failed; no current complete result was published')
        return results
    timeout=submission_timeout();ids=[]
    for n,path in enumerate(paths+([Path(merge).resolve()] if merge else [])):
        argv=['bsub']
        if n==len(paths):argv+=['-w',' && '.join(f'ended({i})' for i in ids)]
        record({'cmd':str(path),'status':'submitting','argv':argv})
        saved['uncertain']={'cmd':str(path),'journal':str(journal)};save()
        try:
            with path.open() as handle:process=bsub_process(argv,handle,timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            record({'cmd':str(path),'status':'uncertain','timeout_seconds':timeout})
            raise RuntimeError(f'LSF submission acknowledgement timed out for {path.name}; no automatic duplicate submission; inspect {journal}') from exc
        match=re.search(r'Job <(\d+)>',process.stdout+process.stderr)
        result={'cmd':str(path),'job_id':match[1] if match else None,'returncode':process.returncode,
                'stdout':process.stdout,'stderr':process.stderr,'status':'submitted' if process.returncode==0 and match else 'uncertain'}
        record(result);print((process.stdout+process.stderr).strip(),flush=True)
        if process.returncode or not match:raise RuntimeError('Submission stopped; inspect '+str(journal))
        ids.append(match[1]);saved.pop('uncertain',None)
        if n==len(paths):saved['merge']=match[1]
        else:saved['workers'][path.stem]=match[1]
        save()
    if wait:
        states=wait_lsf(ids[:len(paths)]);record({'states':states})
        if any(value!='DONE' for value in states.values()):raise RuntimeError('An LSF population worker failed; the final merge validation cannot publish incomplete results')
        if merge:
            states=wait_lsf(ids[-1:]);record({'merge_states':states})
            if any(value!='DONE' for value in states.values()):raise RuntimeError('LSF merge failed')
    return ids


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--submit', choices=['auto', 'local', 'hpc'], default='auto')
    p.add_argument('--jobs-N', type=int, default=int(os.environ.get('jobs_N', '8')))
    p.add_argument('--wait', action='store_true')
    p.add_argument('--merge')
    p.add_argument('--journal', type=Path)
    p.add_argument('--state', type=Path)
    p.add_argument('--plan',type=Path)
    p.add_argument('commands', nargs='+')
    a = p.parse_args()
    try:
        with contextlib.ExitStack() as stack:
            if a.state:
                a.state.parent.mkdir(parents=True, exist_ok=True)
                lock=stack.enter_context(a.state.with_suffix('.lock').open('a'))
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            run(a.commands, a.submit, a.jobs_N, a.wait, a.merge, a.journal, a.state, a.plan)
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as e:
        p.exit(2, f'ERROR: {e}\n')
