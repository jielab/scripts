#!/usr/bin/env python3
"""Keep a command and every descendant owned by its foreground launcher (Linux/WSL).

This supervisor stays in the caller's terminal session/process group. It adopts
orphaned grandchildren and exits only after its workload has stopped. A kernel
parent-death signal also covers a launcher killed with SIGKILL, where Python
finally blocks cannot run. No daemon or detached monitoring service is created.
"""
from __future__ import annotations

import argparse
import ctypes
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def prctl(option, value):
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(option, value, 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))


def arm_parent_death(parent_pid=None):
    expected = os.getppid() if parent_pid is None else parent_pid
    prctl(1, signal.SIGTERM)  # PR_SET_PDEATHSIG, survives exec of this Python process.
    if os.getppid() != expected:
        os.kill(os.getpid(), signal.SIGTERM)


def process_table():
    entries = {}
    for path in Path('/proc').glob('[0-9]*/stat'):
        try:
            fields = path.read_text().rsplit(')', 1)[1].split()
            entries[int(path.parent.name)] = (int(fields[1]), fields[19], fields[0])
        except (OSError, ValueError, IndexError):
            pass
    return entries


def descendants(parent):
    entries = process_table()
    selected = {parent}
    while True:
        children = {pid for pid, (ppid, _, _) in entries.items() if ppid in selected}
        if children <= selected:
            break
        selected |= children
    return {pid: entries[pid][1] for pid in selected - {parent} if pid in entries}


def send(pids, signum):
    entries = process_table()
    for pid, birth in pids.items():
        if pid in entries and entries[pid][1] == birth and entries[pid][2] != 'Z':
            try:
                os.kill(pid, signum)
            except ProcessLookupError:
                pass


def cleanup(proc, grace=3.):
    """Repeated snapshots catch workers that fork or are reparented during exit."""
    deadline = time.monotonic() + grace
    sent = set()
    killed = False
    while True:
        proc.poll()  # Let Popen retain the direct child's real exit status.
        children = descendants(os.getpid())
        for pid in children:
            if pid != proc.pid:
                try:
                    os.waitpid(pid, os.WNOHANG)
                except ChildProcessError:
                    pass
        children = descendants(os.getpid())
        if not children:
            break
        if time.monotonic() >= deadline:
            send(children, signal.SIGKILL)
            killed = True
        else:
            fresh = {pid: birth for pid, birth in children.items() if (pid, birth) not in sent}
            send(fresh, signal.SIGTERM)
            send(fresh, signal.SIGCONT)  # A suspended worker must be able to terminate.
            sent.update(fresh.items())
        if killed and time.monotonic() >= deadline + 5:
            raise RuntimeError('Worker remains in the kernel after SIGKILL; check uninterruptible I/O')
        time.sleep(.05)
    proc.wait()


def supervise(command, parent_pid):
    stopped = []
    suspended = {}

    def cancel(signum, _frame):
        if not stopped:
            stopped.append(signum)

    def suspend(_signum, _frame):
        suspended.update(descendants(os.getpid()))
        send(suspended, signal.SIGSTOP)
        os.kill(os.getpid(), signal.SIGSTOP)

    def resume(_signum, _frame):
        send(suspended, signal.SIGCONT)
        suspended.clear()

    background = os.getenv('SCRIPT_BACKGROUND') == '1'
    signal.signal(signal.SIGHUP, signal.SIG_IGN if background else cancel)
    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, cancel)
    signal.signal(signal.SIGTSTP, suspend)
    signal.signal(signal.SIGCONT, resume)
    prctl(36, 1)  # PR_SET_CHILD_SUBREAPER: adopt descendants when intermediates exit.
    arm_parent_death(parent_pid)
    if stopped:
        return 128 + stopped[0]
    proc = subprocess.Popen(command, start_new_session=False)
    try:
        while proc.poll() is None and not stopped:
            time.sleep(.05)
        rc = 128 + stopped[0] if stopped else proc.returncode
    finally:
        cleanup(proc)
    return rc if rc >= 0 else 128 - rc


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent-pid', type=int, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('A command is required')
    raise SystemExit(supervise(command, args.parent_pid))
