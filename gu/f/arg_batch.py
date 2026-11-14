#!/usr/bin/env python3
"""Independent ARG method x chromosome LSF jobs; no inference on the login node."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import time

from hpc_batch import get, put, snapshot_code, worker_environment
from processor import validate_queue
from qiming_resources import cpu_slots
import submit_jobs

ROOT = Path(__file__).resolve().parents[1]
METHODS = ('threads', 'needle', 'tsinfer', 'singer')


def atomic(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def chromosomes(value):
    result = []
    for token in re.split(r'[,;\s]+', value.strip()):
        token = re.sub(r'^chr', '', token, flags=re.I).upper()
        if re.fullmatch(r'\d+-\d+', token):
            a, b = map(int, token.split('-'))
            if not 1 <= a <= b <= 22:
                raise ValueError('Invalid chromosome range')
            result.extend(map(str, range(a, b + 1)))
        else:
            result.append('X' if token == '23' else token)
    if not result or not set(result) <= set(map(str, range(1, 23))) | {'X'}:
        raise ValueError('--chr must select chromosomes 1..22 or X')
    if len(result) != len(set(result)):
        raise ValueError('Duplicate chromosome selection')
    return result


def directory(value):
    path = Path(value).expanduser().resolve()
    if any(c in str(path) for c in ('\n', '\r', '"')):
        raise ValueError('Directory names cannot contain newlines or double quotes')
    return path


def record(task_dir, phase, code=0):
    job = os.environ['LSB_JOBID']
    path = Path(task_dir) / f'attempt.{job}.json'
    now = time.time()
    data = json.loads(path.read_text()) if path.exists() else {'job_id': job, 'started': now}
    data.update(state='RUN' if phase == 'start' else ('DONE' if code == 0 else 'EXIT'))
    if phase != 'start':
        data.update(finished=now, elapsed_seconds=round(now - data['started'], 3), exit_code=code)
    atomic(path, data)


def plan(raw):
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument('--queue', default=os.environ.get('GU_QUEUE', '33'))
    parser.add_argument('--slots', type=int)
    parser.add_argument('--walltime', default=os.environ.get('HPC_WALLTIME', '120:00'))
    parser.add_argument('--run-dir')
    parser.add_argument('--plan-only', action='store_true')
    opts, args = parser.parse_known_args(raw)
    if not args or args[0] not in ('build', 'prep_gen'):
        raise ValueError('Expected build or prep_gen')
    action = args[0]
    methods = get(args, '--method', 'needle').lower().split(',')
    if not methods or not set(methods) <= set(METHODS) or len(methods) != len(set(methods)):
        raise ValueError('--method must contain distinct threads,needle,tsinfer,singer names')
    chrs = chromosomes(get(args, '--chr', '1-22,X'))
    target = get(args, '--target', '1kg')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', target):
        raise ValueError('Invalid --target')
    build = get(args, '--grch', '37').lower()
    build = {'b37': '37', 'grch37': '37', 'hg19': '37', 'b38': '38', 'grch38': '38', 'hg38': '38'}.get(build, build)
    if build not in ('37', '38'):
        raise ValueError('--grch must be 37 or 38')
    if '--dir-target' in args and '--dir-gen' in args:
        raise ValueError('Use only --dir-target (or its legacy alias --dir-gen)')
    args = ['--dir-target' if x == '--dir-gen' else x for x in args]
    source = directory(get(args, '--dir-target', str(Path(os.environ['GU_REF_ROOT']) / target / build)))
    if not source.is_dir():
        raise ValueError(f'Missing target directory: {source}')
    threads, jobs = int(get(args, '--threads', '8')), int(get(args, '--jobs', '1'))
    memory = int(get(args, '--memory-limit-gb', os.environ.get('ARG_MEMORY_LIMIT_GB', '32')))
    if min(threads, jobs, memory) < 1:
        raise ValueError('--threads, --jobs, --memory-limit-gb must be positive')
    slots = opts.slots if opts.slots is not None else threads
    if slots < threads:
        raise ValueError('--slots must cover --threads')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', opts.queue) or not re.fullmatch(r'\d+:[0-5]\d', opts.walltime):
        raise ValueError('Invalid --queue or --walltime')
    hours, minutes = map(int, opts.walltime.split(':'))
    queue, slots = cpu_slots(opts.queue, slots, memory, hours + minutes / 60)
    validate_queue(queue, False)
    if get(args, '--format', 'native') not in ('native', 'trace'):
        raise ValueError('--format must be native or trace')
    # Fail typos on the login node, using the same public option names as workers.
    shell = (ROOT / 'f/arg.sh').read_text()
    py = (ROOT / 'f/arg.py').read_text()
    known = set(re.findall(r'--[a-z][a-z0-9-]*', shell + py))
    for arg in args[1:]:
        if arg.startswith('--') and arg not in known:
            raise ValueError(f'Unknown ARG option: {arg}')
    # Resolve user paths now; compute nodes need not start in the submission CWD.
    for flag in ('--dir-target', '--dir-pfile', '--dir-vcf', '--dir-map', '--sample-panel',
                 '--ancestry-file', '--keep', '--keep-snv', '--demography', '--map-pattern'):
        value = get(args, flag)
        if value:
            args = put(args, flag, str(directory(value)))
    base = directory(os.environ['GU_ANALYSIS_BASE']) / target / build / 'arg-batches'
    if opts.run_dir:
        run = directory(opts.run_dir)
        if run.exists():
            raise ValueError(f'Run directory already exists; use submit/retry: {run}')
        run.mkdir(parents=True)
    else:
        base.mkdir(parents=True, exist_ok=True)
        run = Path(tempfile.mkdtemp(prefix=time.strftime('%Y%m%d-%H%M%S-'), dir=base))
    frozen = snapshot_code(run)
    environment = worker_environment()
    original_root = str(ROOT.parent)
    for key, value in environment.items():
        if value == original_root or value.startswith(original_root + '/'):
            candidate = frozen.parent / Path(value).relative_to(ROOT.parent)
            if candidate.exists():
                environment[key] = str(candidate)
    environment.update(GU_CODE_ROOT=str(frozen.parent), GU_ARG_BATCH_WORKER='1',
                       GU_WORKER_CODE_ROOT=str(frozen), submit='hpc', SUBMIT='hpc',
                       WORKFLOW_WORKER='1', ARG_RESOURCE_GUARD='scheduler')
    output_base = directory(get(args, '--dir-arg', str(run / 'outputs')))
    prep_base = directory(get(args, '--dir-gen4arg', str(run / 'prepared')))
    if output_base == prep_base:
        raise ValueError('--dir-arg and --dir-gen4arg must differ')
    scratch_value = get(args, '--scratch')
    scratch_base = directory(scratch_value) if scratch_value else None
    if scratch_base in (output_base, prep_base):
        raise ValueError('--scratch must differ from output/preparation directories')
    q = shlex.quote
    python = str(Path(sys.executable).resolve())
    tasks = []
    for method in methods:
        for chrom in chrs:
            key = f'{method}.chr{chrom}'
            task = run / 'jobs' / key
            task.mkdir(parents=True)
            out, prep = output_base / method / f'chr{chrom}', prep_base / method / f'chr{chrom}'
            cmd = list(args)
            values = {'--method': method, '--chr': chrom, '--target': target, '--grch': build,
                      '--dir-target': str(source), '--dir-arg': str(out), '--dir-gen4arg': str(prep),
                      '--threads': str(threads), '--jobs': str(jobs if method == 'singer' else 1),
                      '--memory-limit-gb': str(memory), '--resource-guard': 'scheduler'}
            if scratch_base:
                if method == 'singer':
                    # SINGER's window scratch lives under its isolated output directory.
                    raise ValueError('SINGER does not accept --scratch; omit it for mixed-method runs')
                values['--scratch'] = str(scratch_base / method / f'chr{chrom}')
            for flag, value in values.items():
                cmd = put(cmd, flag, value)
            name = f'gu.arg.{target}.{key}'
            script = task / 'job.cmd'
            recorder = shlex.join([python, str(frozen / 'f/arg_batch.py'), '_record', str(task)])
            lines = ['#!/bin/bash -l', f'#BSUB -q {queue}', f'#BSUB -n {slots}',
                     f'#BSUB -W {opts.walltime}', f'#BSUB -J {name}',
                     f'#BSUB -R "span[hosts=1] rusage[mem={memory}GB/host]"',
                     f'#BSUB -o "{task}/%J.out"', f'#BSUB -e "{task}/%J.err"',
                     'set -euo pipefail', ': "${LSB_JOBID:?Submit this file through LSF}"',
                     f'exec 7>{q(str(task / ".worker.lock"))}',
                     'flock -n 7 || { echo "Another worker is using this task" >&2; exit 75; }',
                     f'{recorder} start',
                     'finish() { code=$?; trap - EXIT; ' + recorder + ' finish "$code" || true; exit "$code"; }',
                     'trap finish EXIT', 'trap "exit 143" TERM', 'trap "exit 130" INT',
                     *[f'export {key}={q(value)}' for key, value in sorted(environment.items())],
                     f'export SCRIPT_LOG_DIR={q(str(task / "console"))}',
                     f'source {q(str(frozen / "f/hpc.runtime.env"))}',
                     f'cd {q(str(frozen))}',
                     shlex.join(['bash', str(frozen / 'arg.sh')] + cmd), '']
            script.write_text('\n'.join(lines))
            script.chmod(0o755)
            tasks.append(dict(key=key, method=method, chromosome=chrom, cmd=str(script),
                              directory=str(task), output=str(out), prepared=str(prep), argv=cmd))
    manifest = dict(version=1, run_dir=str(run), created=time.time(), target=target, build=build,
                    queue=queue, slots=slots, threads=threads, memory_gb=memory, tasks=tasks)
    atomic(run / 'plan.json', manifest)
    print(f'RUN_DIR={run}\nGenerated {len(tasks)} independent LSF jobs ({len(methods)} methods x {len(chrs)} chromosomes).', flush=True)
    print(f'Each job: queue={queue}, slots={slots}, threads={threads}, RAM={memory} GiB', flush=True)
    if opts.plan_only:
        print('Plan only: no jobs submitted. Use: ./arg.sh submit ' + shlex.quote(str(run)))
    else:
        manage('submit', run)
    return run


def task_states(plan):
    run = Path(plan['run_dir'])
    journal = run / 'submission.jsonl'
    latest = {}
    if journal.exists():
        for line in journal.read_text().splitlines():
            row = json.loads(line)
            if 'cmd' in row:
                latest[row['cmd']] = row
    ids = [r['job_id'] for r in latest.values() if r.get('job_id')]
    live = {}
    if ids:
        try:
            probe = subprocess.run(['bjobs', '-a', '-noheader', '-o', 'jobid stat', *ids],
                                   capture_output=True, text=True, timeout=30)
            for line in probe.stdout.splitlines():
                fields = line.split()
                if len(fields) == 2 and fields[0] in ids:
                    live[fields[0]] = fields[1]
        except (OSError, subprocess.SubprocessError):
            pass  # An absent job is UNKNOWN, never an implicit success/failure.
    cache_path = run / 'terminal.json'
    terminal = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    terminal.update({job: state for job, state in live.items() if state in ('DONE', 'EXIT')})
    atomic(cache_path, terminal)
    result = []
    for task in plan['tasks']:
        entry = latest.get(task['cmd'])
        job = entry.get('job_id') if entry else None
        receipt = Path(task['directory']) / f'attempt.{job}.json'
        data = json.loads(receipt.read_text()) if job and receipt.exists() else {}
        state = 'NOT_SUBMITTED' if entry is None else 'UNKNOWN'
        if job:
            state = live.get(job, terminal.get(job, 'UNKNOWN'))
            if state == 'UNKNOWN' and data.get('state') in ('DONE', 'EXIT'):
                state = data['state']
        elapsed = data.get('elapsed_seconds')
        if elapsed is None and data.get('started') and state == 'RUN':
            elapsed = round(time.time() - data['started'], 1)
        result.append(dict(task=task, job_id=job, state=state, elapsed_seconds=elapsed))
    return result


def manage(action, run, plan_only=False):
    run = directory(run)
    plan_data = json.loads((run / 'plan.json').read_text())
    if Path(plan_data['run_dir']) != run:
        raise ValueError('Saved batches contain absolute paths; do not move the run directory')
    with (run / '.submit.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        rows = task_states(plan_data)
        print('METHOD\tCHR\tJOB_ID\tSTATUS\tELAPSED_SECONDS')
        for row in rows:
            t = row['task']
            elapsed = row['elapsed_seconds']
            print(f"{t['method']}\t{t['chromosome']}\t{row['job_id'] or '-'}\t{row['state']}\t{elapsed if elapsed is not None else '-'}")
        if action == 'status':
            return rows
        wanted = 'EXIT' if action == 'retry' else 'NOT_SUBMITTED'
        paths = [Path(r['task']['cmd']) for r in rows if r['state'] == wanted]
        print(f'{action}: {len(paths)} {wanted} jobs; active, completed and UNKNOWN jobs are excluded.', flush=True)
        if plan_only:
            for p in paths:
                print(p)
        elif paths:
            submit_jobs.run(paths, mode='hpc', wait=os.environ.get('SUBMIT_WAIT') == 'TRUE',
                            journal=run / 'submission.jsonl')
        return rows


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if args and args[0] == '_record':
        record(args[1], args[2], int(args[3]) if len(args) > 3 else 0)
    elif args and args[0] in ('status', 'submit', 'retry'):
        p = argparse.ArgumentParser()
        p.add_argument('action', choices=('status', 'submit', 'retry'))
        p.add_argument('run_dir')
        p.add_argument('--plan-only', action='store_true')
        a = p.parse_args(args)
        manage(a.action, a.run_dir, a.plan_only)
    else:
        args = [part for arg in args for part in (arg.split('=', 1) if arg.startswith('--') and '=' in arg else [arg])]
        plan(args)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(2)
