"""Fail-closed cgroup isolation for external MetaCache builders."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import uuid
import fcntl


def scope_path():
    group = next((line[3:] for line in Path('/proc/self/cgroup').read_text().splitlines()
                  if line.startswith('0::')), '')
    return group, Path('/sys/fs/cgroup') / group.lstrip('/')


def verify_scope(size, prefix):
    group, base = scope_path()
    expected = {'memory.max': str(size), 'memory.swap.max': '0', 'memory.oom.group': '1'}
    if not group or prefix not in group or any(
            (base / key).read_text().strip() != value for key, value in expected.items()):
        raise RuntimeError('Memory isolation verification failed; refusing unrestricted work')
    return group


def in_workflow_scope(memory_gb=16):
    group, _ = scope_path()
    if 'metag-workflow-' not in group:
        return False
    verify_scope(memory_gib(memory_gb) * 2**30, 'metag-workflow-')
    return True


def require_workflow_scope(memory_gb=16):
    if not in_workflow_scope(memory_gb):
        raise RuntimeError('The complete workflow must run inside its shared memory scope; use ./1.library.sh')


def run_workflow_limited(command, memory_gb=16, stdout=None, stderr=None):
    return _run_limited(command, stdout=stdout, stderr=stderr,
                        memory_gb=memory_gb, scope_kind='workflow')


def memory_gib(value):
    number = int(value)
    if not 1 <= number <= 16:
        raise ValueError('MetaCache memory limit must be between 1 and 16 GiB on this shared host')
    return number


def run_limited(command, *, stdout, stderr, memory_gb=16):
    # Only one builder per user, even if two terminals launch library builds.
    runtime = Path(f'/run/user/{os.getuid()}')
    try:
        lock = (runtime / 'metacache-build.lock').open('a')
    except OSError as exc:
        raise RuntimeError('MetaCache was not started: user runtime directory is unavailable') from exc
    with lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Another MetaCache build is running; concurrent builds are refused') from exc
        if in_workflow_scope(memory_gb):
            # Stay in the inherited outer scope. A new systemd scope here would
            # be a sibling, escaping the combined preparer + GPU memory cap.
            result = subprocess.run(command, stdout=stdout, stderr=stderr)
            check_result(result, memory_gb)
            return result
        return _run_limited(command, stdout=stdout, stderr=stderr, memory_gb=memory_gb)


def _run_limited(command, *, stdout, stderr, memory_gb, scope_kind='build'):
    size = memory_gib(memory_gb) * 2**30
    env = dict(os.environ)
    env.setdefault('XDG_RUNTIME_DIR', f'/run/user/{os.getuid()}')
    env.setdefault('DBUS_SESSION_BUS_ADDRESS', f"unix:path={env['XDG_RUNTIME_DIR']}/bus")
    try:
        subprocess.run(['systemctl', '--user', 'show-environment'], env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                       check=True, timeout=5)
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError('MetaCache was not started: user systemd is required for memory isolation') from exc
    unit = ('metag-workflow-' if scope_kind == 'workflow' else 'metacache-build-') + uuid.uuid4().hex
    args = ['systemd-run', '--user', '--scope', '--quiet', '--unit=' + unit,
            '-p', f'MemoryMax={size}', '-p', 'MemorySwapMax=0',
            '-p', 'OOMPolicy=kill', '-p', 'TasksMax=256',
            sys.executable, str(Path(__file__).resolve()),
            *(['--workflow'] if scope_kind == 'workflow' else []), str(size), *map(str, command)]
    result = subprocess.run(args, stdout=stdout, stderr=stderr, env=env)
    check_result(result, memory_gb)
    return result


def check_result(result, memory_gb):
    if result.returncode == 125:
        raise RuntimeError('MetaCache was not started: cgroup verification failed; see build log')
    if result.returncode in (-9, 137):
        raise RuntimeError(f'MetaCache was killed inside its {memory_gb} GiB limit; '
                           'build stopped without CPU retry. Completed shards are retained; '
                           'inspect the bounded shard and build log before retrying.')


def guarded_exec(size, command, workflow=False):
    group = verify_scope(size, 'metag-workflow-' if workflow else 'metacache-build-')
    print(f'{"Workflow" if workflow else "MetaCache"} verified cgroup={group} RAM={size} bytes swap=0 oom.group=1', flush=True)
    os.execvp(command[0], command)


if __name__ == '__main__':
    try:
        workflow = sys.argv[1:2] == ['--workflow']
        args = sys.argv[2:] if workflow else sys.argv[1:]
        guarded_exec(int(args[0]), args[1:], workflow=workflow)
    except (OSError, ValueError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(125)
